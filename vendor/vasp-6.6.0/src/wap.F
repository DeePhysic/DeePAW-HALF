#include "symbol.inc"

#define use_wap_supercell

MODULE wap
  use prec
  use tutor, ONLY: vtutor
  IMPLICIT NONE

  integer, parameter :: FILENAME_MAX_LEN = 100

CONTAINS

  SUBROUTINE getwave(EFERMI, W, CDIJ, CQIJ, GRIDC, SV, CVTOT, T_INFO, P, NONL_S, LATT_CUR, IO)
    USE prec
    USE wave
    USE base
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl
    USE constant
    USE pseudo
    USE POSCAR
    USE vhdf5
    REAL(q), intent(in) :: EFERMI
    TYPE (wavespin), intent(in) :: W
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    OVERLAP, intent(in) :: CDIJ(P(1)%LMDIM,P(1)%LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    OVERLAP, intent(in) :: CQIJ(P(1)%LMDIM,P(1)%LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    TYPE (grid_3d), intent(in) :: GRIDC
    RGRID, intent(in) :: SV(DIMREAL(W%WDES%GRID%MPLWV),W%WDES%NCDIJ)
    COMPLEX(q), intent(in) :: CVTOT(GRIDC%MPLWV,W%WDES%NCDIJ)
    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(latt), intent(in) :: LATT_CUR
    TYPE(in_struct), intent(in) :: IO

    !
    ! Local variables
    !
    TYPE (wavedes), POINTER :: WDES
    TYPE (grid_3d), POINTER :: GRID
    COMPLEX(8), ALLOCATABLE :: POT(:,:,:,:), POT_SV(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: WAVES(:,:,:,:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: DPROJECTORS(:,:,:,:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCDIJ(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCQIJ(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCRIJ(:,:,:,:,:)
    REAL(8), ALLOCATABLE :: EIGENVALUE(:,:,:)
    REAL(8), ALLOCATABLE :: KPOINTS(:,:)
    REAL(8), ALLOCATABLE :: QTOT(:,:,:)
    INTEGER :: LMDIM, LDIM
    INTEGER :: MSIZE
    INTEGER NUNIT
    LOGICAL INUSE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                      ALIAS                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    GRID => W%WDES%GRID
    WDES => W%WDES

    LMDIM = P(1)%LMDIM
    LDIM = P(1)%LDIM

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                    PAW q_ij                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(QTOT(LDIM, LDIM, T_INFO%NIONS))
    QTOT = 0.0d0
    !CALL getwave_QTOT(P, T_INFO, QTOT)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                    PAW R_ij                       !
    !       RIJ=-(<J|nabla|I>-<J_ps|nabla|I_ps>)        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE( CCRIJ( LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ, 2))
    CCRIJ = dcmplx(0)
    CALL getwave_RIJ(WDES%NCDIJ,LMDIM,P, T_INFO, CCRIJ)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 PAW PARAMETERS                    !
    !              QIJ= <J|I>-<J_ps|I_ps>               !
    !          DIJ= <J|h|I>-<J_ps|h_ps|I_ps>            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(CCDIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ))
    CCDIJ = dcmplx(0)
    ALLOCATE(CCQIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ))
    CCQIJ = dcmplx(0)
    CALL getwave_QIJ_DIJ(WDES, T_INFO, LMDIM, CQIJ, CDIJ, CCDIJ, CCQIJ)

    MSIZE = T_INFO%NIONS * LMDIM * LMDIM * WDES%NCDIJ * 2
    CALLMPI(M_sum_d(WDES%COMM_INB, CCQIJ, MSIZE))
    CALLMPI(M_sum_d(WDES%COMM_INB, CCDIJ, MSIZE))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                  LOCAL POTENTIAL                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ALLOCATE(POT_SV(GRID%NGX, GRID%NGY, GRID%NGZ, WDES%NCDIJ))
    POT_SV = dcmplx(0)
    ALLOCATE(POT(GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ, WDES%NCDIJ))
    POT = dcmplx(0)
    CALL getwave_local_potential(WDES, GRID, SV, POT_SV)
    CALL getwave_local_potential_CVTOT(WDES, GRIDC, CVTOT, POT)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                    KPOINTS                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(KPOINTS(1:3, 1:WDES%NKPTS))
    KPOINTS = 0.0d0

#ifdef use_wap_supercell
    ! KPOINTS are overwritten in getwave_set_WAVES_DPROJECTORS unless
    ! use_wap_supercell defined. This is to make KPOINTS order be surely
    ! the same as WAVES, DPROJECTORS, and EIGENVALUE.
    call getwave_set_KPOINTS(WDES, KPOINTS)
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       DERIVATIVES OF PROJECTORS and WAVES       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef use_wap_supercell
    ALLOCATE(DPROJECTORS(0:3,1:LMDIM,1:T_INFO%NIONS,1:WDES%NRSPINORS,&
         1:WDES%ISPIN,1:WDES%NB_TOT, 1:WDES%NKPTS))
    DPROJECTORS = dcmplx(0)

    ALLOCATE(WAVES(GRID%NGX, GRID%NGY, GRID%NGZ,&
         WDES%NRSPINORS, WDES%ISPIN, WDES%NB_TOT,  WDES%NKPTS))
    WAVES = dcmplx(0)

    ALLOCATE(EIGENVALUE(WDES%ISPIN,WDES%NB_TOT,WDES%NKPTS))   ! change the args order
    EIGENVALUE = dble(0)


    CALL getwave_set_WAVES_DPROJECTORS(T_INFO, GRID, NONL_S, W,&
         WDES, LATT_CUR, LMDIM, KPOINTS, WAVES, EIGENVALUE, DPROJECTORS)

    MSIZE = WDES%ISPIN * WDES%NKPTS * WDES%NB_TOT
    CALLMPI(M_sum_d(WDES%COMM_INTER, EIGENVALUE, MSIZE))

    MSIZE = T_INFO%NIONS * LMDIM * WDES%NKPTS * WDES%NB_TOT * WDES%ISPIN * WDES%NRSPINORS * 4 * 2
    CALLMPI(M_sum_d(WDES%COMM, DPROJECTORS, MSIZE))

    MSIZE = GRID%NGX * GRID%NGY * GRID%NGZ * WDES%NRSPINORS * WDES%ISPIN * WDES%NKPTS * WDES%NB_TOT * 2
    CALLMPI(M_sum_d(WDES%COMM, WAVES, MSIZE))

#endif

    IF (WDES%COMM%NODE_ME == WDES%COMM%IONODE) THEN

       NUNIT = 71

       INQUIRE(UNIT=NUNIT, OPENED=INUSE)
       IF (INUSE) THEN
          CALL vtutor%error('Error to open file in wap.F')
       END IF

#ifndef use_wap_supercell
       !! Uncommenting below, data are written in text format which may be
       !! used for the test.
       ! CALL getwave_write_data_in_text(NUNIT, T_INFO, WDES, &
       !      LMDIM, GRID, GRIDC, CCDIJ, CCQIJ, POT_SV, POT, &
       !      DPROJECTORS, WAVES)
#endif

       CALL getwave_write_inwap_yaml(NUNIT, T_INFO, WDES, 3, LATT_CUR, &
            P, GRID, GRIDC, KPOINTS, EFERMI)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !                 Write in binary                 !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef use_wap_supercell
       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='WAVES.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) WAVES
       CLOSE(NUNIT)

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='DPROJECTORS.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) DPROJECTORS
       CLOSE(NUNIT)

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='EIGENVALUE.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) EIGENVALUE
       CLOSE(NUNIT)

       DEALLOCATE(WAVES)
       DEALLOCATE(DPROJECTORS)
       DEALLOCATE(EIGENVALUE)
#endif

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='LOCAL-POTENTIAL.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) POT
       CLOSE(NUNIT)

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='PAW-OVERLAP.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) CCQIJ
       CLOSE(NUNIT)

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='PAW-STRENGTH.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) CCDIJ
       CLOSE(NUNIT)

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='QTOT.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) QTOT
       CLOSE(NUNIT)

       OPEN(UNIT=NUNIT, STATUS='REPLACE', FILE='PAW-Rnij.bin',&
            ACCESS='STREAM', FORM='UNFORMATTED')
       WRITE(NUNIT) CCRIJ
       CLOSE(NUNIT)

    ENDIF

    DEALLOCATE(CCRIJ)
    DEALLOCATE(POT)
    DEALLOCATE(CCQIJ)
    DEALLOCATE(CCDIJ)
    DEALLOCATE(KPOINTS)
  END SUBROUTINE getwave

#ifdef VASP_HDF5
  ! TODO: For the moment this routine is a copy from about but later on it should be cleaned up and
  ! moved to the main VASP code and completely independent from wap.F
  SUBROUTINE vh5_write_vtot(W, CDIJ, CQIJ, GRIDC, SV, CVTOT, T_INFO, P, NONL_S, LATT_CUR, IO)
    USE prec
    USE wave
    USE base
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl
    USE constant
    USE pseudo
    USE POSCAR
    USE vhdf5
    TYPE (wavespin), intent(in) :: W
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    OVERLAP, intent(in) :: CDIJ(P(1)%LMDIM,P(1)%LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    OVERLAP, intent(in) :: CQIJ(P(1)%LMDIM,P(1)%LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    TYPE (grid_3d), intent(in) :: GRIDC
    RGRID, intent(in) :: SV(DIMREAL(W%WDES%GRID%MPLWV),W%WDES%NCDIJ)
    COMPLEX(q), intent(in) :: CVTOT(GRIDC%MPLWV,W%WDES%NCDIJ)
    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(latt), intent(in) :: LATT_CUR
    TYPE(in_struct), intent(in) :: IO

    !
    ! Local variables
    !
    TYPE (wavedes), POINTER :: WDES
    TYPE (grid_3d), POINTER :: GRID
    COMPLEX(8), ALLOCATABLE :: POT(:,:,:,:), POT_SV(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCDIJ(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCQIJ(:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCRIJ(:,:,:,:,:)
    REAL(8), ALLOCATABLE :: QTOT(:,:,:)
    INTEGER :: LMDIM, LDIM
    INTEGER :: MSIZE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                      ALIAS                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    GRID => W%WDES%GRID
    WDES => W%WDES

    LMDIM = P(1)%LMDIM
    LDIM = P(1)%LDIM

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                    PAW q_ij                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(QTOT(LDIM, LDIM, T_INFO%NIONS))
    QTOT = 0.0d0
    !CALL getwave_QTOT(P, T_INFO, QTOT)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                    PAW R_ij                       !
    !       RIJ=-(<J|nabla|I>-<J_ps|nabla|I_ps>)        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE( CCRIJ( LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ, 2))
    CCRIJ = dcmplx(0)
    CALL getwave_RIJ(WDES%NCDIJ,LMDIM,P, T_INFO, CCRIJ)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 PAW PARAMETERS                    !
    !              QIJ= <J|I>-<J_ps|I_ps>               !
    !          DIJ= <J|h|I>-<J_ps|h_ps|I_ps>            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(CCDIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ))
    CCDIJ = dcmplx(0)
    ALLOCATE(CCQIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ))
    CCQIJ = dcmplx(0)
    CALL getwave_QIJ_DIJ(WDES, T_INFO, LMDIM, CQIJ, CDIJ, CCDIJ, CCQIJ)

    MSIZE = T_INFO%NIONS * LMDIM * LMDIM * WDES%NCDIJ * 2
    CALLMPI(M_sum_d(WDES%COMM_INB, CCQIJ, MSIZE))
    CALLMPI(M_sum_d(WDES%COMM_INB, CCDIJ, MSIZE))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                  LOCAL POTENTIAL                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ALLOCATE(POT_SV(GRID%NGX, GRID%NGY, GRID%NGZ, WDES%NCDIJ))
    POT_SV = dcmplx(0)
    ALLOCATE(POT(GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ, WDES%NCDIJ))
    POT = dcmplx(0)
    CALL getwave_local_potential(WDES, GRID, SV, POT_SV)
    CALL getwave_local_potential_CVTOT(WDES, GRIDC, CVTOT, POT)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 WRITE TO HDF5                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL WRITE_VTOT_HDF5(IH5WAVEFILEID,POT,CCQIJ,CCDIJ,QTOT,CCRIJ)

    DEALLOCATE(CCRIJ)
    DEALLOCATE(POT)
    DEALLOCATE(CCQIJ)
    DEALLOCATE(CCDIJ)
  END SUBROUTINE vh5_write_vtot
#endif


  SUBROUTINE getwave_QIJ_DIJ(WDES, T_INFO, LMDIM, CQIJ, CDIJ, CCDIJ, CCQIJ)
    USE prec
    USE wave
    USE POSCAR

    TYPE(wavedes), intent(in):: WDES
    TYPE(type_info), intent(in) :: T_INFO
    INTEGER, intent(in):: LMDIM
    OVERLAP, intent(in):: CDIJ(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ)
    OVERLAP, intent(in):: CQIJ(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ)
    COMPLEX(8), intent(inout):: CCDIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ)
    COMPLEX(8), intent(inout):: CCQIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ)

    INTEGER:: ISP, AT, I, J, NI

    !!! See WDES_SET_NPRO, NI_GLOBAL in wave.F for index of ions
    DO ISP = 1, WDES%NCDIJ
       DO AT = 1, WDES%NIONS
          NI = NI_GLOBAL(AT, WDES%COMM_INB)
          DO I = 1, LMDIM
             DO J = 1, LMDIM
                CCQIJ(I, J, NI, ISP) = dcmplx(CQIJ(I, J, AT, ISP))
                CCDIJ(I, J, NI, ISP) = dcmplx(CDIJ(I, J, AT, ISP))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE getwave_QIJ_DIJ


  SUBROUTINE getwave_QTOT(P, T_INFO, QTOT)
    USE prec
    USE PSEUDO
    USE POSCAR

    TYPE(type_info), intent(in) :: T_INFO
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    REAL(8), intent(inout):: QTOT(P(1)%LDIM, P(1)%LDIM, T_INFO%NIONS)

    INTEGER:: NIS, NT, NI, L, LP, LDIM
    TYPE (potcar), POINTER :: PP

    LDIM = P(1)%LDIM
    NIS = 1
    DO NT = 1, T_INFO%NTYP
       DO NI = NIS, T_INFO%NITYP(NT) + NIS - 1
          ! QTOT is not distributed, so the following two lines are not needed.
          ! NIP = NI_LOCAL(NI, WDES%COMM_INB)
          ! IF (NIP == 0) CYCLE
          PP=>PP_POINTER(P, NI, NT)
          DO L = 1, LDIM
             DO LP = 1, LDIM
                QTOT(L, LP, NI) = PP%QTOT(L, LP)
             END DO
          END DO
       END DO
       NIS = NIS + T_INFO%NITYP(NT)
    END DO

  END SUBROUTINE getwave_QTOT


  SUBROUTINE getwave_local_potential(WDES, GRID, SV, POT)
    USE mgrid
    USE wave

    TYPE(wavedes), intent(in) :: WDES
    TYPE(grid_3d), intent(in) :: GRID
    RGRID, intent(in) :: SV(DIMREAL(GRID%MPLWV), WDES%NCDIJ)
    COMPLEX(8), intent(inout) :: POT(GRID%NGX, GRID%NGY, GRID%NGZ, WDES%NCDIJ)

    integer :: ISP, NC, NX, NY, NZ

    DO ISP = 1, WDES%NCDIJ
       DO NZ = 1, GRID%NGZ
          IF (GRID%RL%NFAST == 1) THEN
             DO NC = 1, GRID%RL%NCOL
                IF (GRID%RL%I3(NC) /= NZ) CYCLE
                NY = GRID%RL%I2(NC)
                DO NX = 1, GRID%NGX
                   POT(NX, NY, NZ, ISP) = dcmplx(SV(NX + (NC - 1) * GRID%NGX, ISP))
                ENDDO
             ENDDO
          ELSEIF (GRID%RL%NFAST == 2) THEN
             CALL vtutor%error('error RL%NFAST=2 not implemented')
          ELSEIF (GRID%RL%NFAST == 3) THEN
             DO NC = 1, GRID%RL%NCOL
                NX = GRID%RL%I2(NC)
                NY = GRID%RL%I3(NC)
                POT(NX, NY, NZ, ISP) = dcmplx(SV(NZ + (NC - 1) * GRID%NGZ, ISP))
             ENDDO
          ENDIF
       ENDDO

    ENDDO
    CALLMPI(M_sum_d(GRID%COMM, POT, 2 * GRID%NGX * GRID%NGY * GRID%NGZ * WDES%NCDIJ))

  END SUBROUTINE getwave_local_potential


  SUBROUTINE getwave_set_KPOINTS(WDES, KPOINTS)

    USE wave

    TYPE(wavedes), intent(in) :: WDES
    REAL(q), intent(inout) :: KPOINTS(1:3,1:WDES%NKPTS)

    integer :: I, NK
    DO NK = 1, WDES%NKPTS
       IF (MOD(NK - 1, WDES%COMM_KINTER%NCPU) .NE. WDES%COMM_KINTER%NODE_ME - 1) CYCLE
       DO I = 1, 3
          KPOINTS(I, NK) = WDES%VKPT(I, NK)
       END DO
    ENDDO
  END SUBROUTINE getwave_set_KPOINTS


  SUBROUTINE getwave_set_WAVES_DPROJECTORS(T_INFO, GRID, NONL_S, W,&
       WDES, LATT_CUR, LMDIM, KPOINTS, WAVES, EIGENVALUE, DPROJECTORS)

    USE LATTICE
    USE mgrid
    USE POSCAR
    USE PSEUDO
    USE nonl_struct_def
    USE nonl
    USE wave
    USE wave_high

    TYPE(wavespin), intent(in) :: W
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(grid_3d), intent(in) :: GRID
    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(wavedes), intent(in) :: WDES
    TYPE(latt), intent(in) :: LATT_CUR
    integer, intent(in) :: LMDIM
    REAL(q), intent(inout) :: KPOINTS(1:3,1:WDES%NKPTS)
    COMPLEX(8), intent(inout) :: WAVES(GRID%NGX, GRID%NGY, GRID%NGZ,&
         WDES%NRSPINORS,WDES%ISPIN,WDES%NB_TOT,WDES%NKPTS)
    REAL(8), intent(inout) :: EIGENVALUE(WDES%ISPIN,  WDES%NB_TOT, WDES%NKPTS)
    COMPLEX(8), intent(inout) :: DPROJECTORS(0:3,1:LMDIM,1:T_INFO%NIONS,1:WDES%NRSPINORS, &
         1:WDES%ISPIN, 1:WDES%NB_TOT, 1:WDES%NKPTS)

    integer :: I, ISP, N, NB, NK
    TYPE(wavefuna) WA
    TYPE(wavedes1) :: WDES1
    TYPE(wavefun1) :: W1(0:3)

    CALL SETWDES(WDES, WDES1, 0)
    DO I = 0, 3
       CALL NEWWAV(W1(I), WDES1, .TRUE.)
    ENDDO

    DO NK = 1, WDES%NKPTS
       IF (MOD(NK - 1, WDES%COMM_KINTER%NCPU) .NE. WDES%COMM_KINTER%NODE_ME - 1) CYCLE
       CALL SETWDES(WDES, WDES1, NK)
       CALL PHASE(WDES, NONL_S, NK)

       DO I = 1, 3
          KPOINTS(I, NK) = WDES1%VKPT(I)
       END DO

       DO ISP = 1, WDES%ISPIN
          WA = ELEMENTS(W, WDES1, ISP)

          DO NB = 1, WDES%NB_TOT
             ! NB_LOCAL returns the local storage index of a band (wave.F)
             N = NB_LOCAL(NB, WDES1)
             IF (N == 0) CYCLE

             DO I = 0, 3
                CALL W1_COPY(ELEMENT(WA, N), W1(I))
                CALL FFTWAV_W1(W1(I))
             ENDDO

             EIGENVALUE(ISP,NB, NK) = dble(W%CELEN(N, NK, ISP))
             CALL getwave_set_WAVES(GRID, W1(0), ISP, NK, NB, WDES, WAVES)
             CALL getwave_set_DPROJECTORS(T_INFO%NIONS, NONL_S, W1(:),&
                  NB, NK, ISP, LMDIM, WDES, WDES1, LATT_CUR, &
                  DPROJECTORS)
          ENDDO
       ENDDO
    ENDDO

    DO I = 0, 3
       CALL DELWAV(W1(I), .TRUE.)
    ENDDO

  END SUBROUTINE getwave_set_WAVES_DPROJECTORS


  SUBROUTINE getwave_set_WAVES(GRID, WD0, ISP, NK, NB, WDES, WAVES)
    USE mgrid
    USE wave

    integer, intent(in) :: ISP, NK, NB
    TYPE(grid_3d), intent(in) :: GRID
    TYPE(wavedes), intent(in) :: WDES
    TYPE(wavefun1), intent(in) :: WD0
    COMPLEX(8), intent(inout) :: WAVES(GRID%NGX, GRID%NGY, GRID%NGZ,&
         WDES%NRSPINORS, WDES%ISPIN, WDES%NB_TOT,  WDES%NKPTS)

    integer :: ISPINOR, NC, NX, NY, NZ

    DO ISPINOR = 0, WDES%NRSPINORS - 1
       DO NZ = 1, GRID%NGZ
          IF (GRID%RL%NFAST == 1) THEN
             DO NC = 1, GRID%RL%NCOL
                IF (GRID%RL%I3(NC) /= NZ) CYCLE
                NY = GRID%RL%I2(NC)
                DO NX = 1, GRID%NGX
                   WAVES(NX, NY, NZ, ISPINOR + 1, ISP, NB, NK) =&
                        dcmplx(WD0%CR(NX + (NC - 1) * GRID%NGX + ISPINOR * GRID%MPLWV))
                ENDDO
             ENDDO
          ELSEIF (GRID%RL%NFAST == 3) THEN
             DO NC = 1, GRID%RL%NCOL
                NX = GRID%RL%I2(NC)
                NY = GRID%RL%I3(NC)
                WAVES(NX, NY, NZ, ISPINOR + 1, ISP, NB, NK) =&
                     dcmplx(WD0%CR(NZ +(NC - 1) * GRID%NGZ + ISPINOR * GRID%MPLWV))
             ENDDO
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE getwave_set_WAVES


  SUBROUTINE getwave_set_DPROJECTORS(NIONS, NONL_S, W1,&
       NB, NK, ISP, LMDIM, WDES, WDES1, LATT_CUR, DPROJECTORS)
    USE LATTICE
    USE pseudo_struct_def
    USE nonl_struct_def
    USE PSEUDO
    USE nonl
    USE wave
    USE string, only: str

    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(wavefun1), intent(in) :: W1(0:3)
    integer, intent(in) :: LMDIM, NB, NK, ISP, NIONS
    TYPE(wavedes), intent(in) :: WDES
    TYPE(wavedes1), intent(in) :: WDES1
    TYPE(latt), intent(in) :: LATT_CUR
    COMPLEX(8), intent(inout) :: DPROJECTORS(0:3,1:LMDIM,1:NIONS, 1:WDES%NRSPINORS, &
         1:WDES%ISPIN,1:WDES%NB_TOT,1:WDES%NKPTS)

    integer :: I, II, ISPINOR, LM, LMBASE, NI, NT

    DO I = 0, 3
       LMBASE = 0
       DO ISPINOR = 1, WDES%NRSPINORS
          NI = 0
          DO NT = 1, NONL_S%NTYP
             DO II = 1, NONL_S%NITYP(NT)
                NI = NI + 1
                CALL PROJ_DER(NONL_S, WDES1, W1(I), LATT_CUR, NI, I)
                IF (NI_LOCAL(NI, WDES%COMM_INB) == 0) CYCLE
                DO LM = 1, NONL_S%LMMAX(NT)
                   DPROJECTORS(I,LM,NI,ISPINOR,ISP,NB,NK) = &
                        dcmplx(W1(I)%CPROJ(LMBASE + LM))
                ENDDO
                LMBASE = LMBASE + NONL_S%LMMAX(NT)
             ENDDO
          ENDDO
       ENDDO

       IF  (LMBASE .NE. WDES%NPRO) THEN
          CALL vtutor%bug('error in N-PRO. expected: '//str(LMBASE)//' found: '//str(WDES%NPRO),__FILE__,__LINE__)
       ENDIF
    ENDDO

  END SUBROUTINE getwave_set_DPROJECTORS


  SUBROUTINE get_DPROJ(LATT_CUR, NONL_S, W1, WDES1, DPROJECTOR)
    USE LATTICE
    USE pseudo_struct_def
    USE nonl_struct_def
    USE PSEUDO
    USE nonl
    USE wave
    USE string, only: str

    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(wavefun1), intent(in) :: W1
    TYPE(wavedes1), intent(in) :: WDES1
    TYPE(latt), intent(in) :: LATT_CUR
    GDEF, intent(inout) :: DPROJECTOR(:,0:,:)
    ! Local variables
    TYPE(wavefun1) :: W1_LOCAL
    integer :: I, II, ISPINOR, LM, LMBASE, NI, NT
    GDEF,target :: CPROJ_LOCAL(WDES1%NPRO)

    ! I don't want to mess with CPROJ from W1 as it is defined as intent(in)
    W1_LOCAL = W1
    W1_LOCAL%CPROJ => CPROJ_LOCAL

    DO I = 0, 3
       LMBASE = 0
       NI = 0
       DO NT = 1, NONL_S%NTYP
          DO II = 1, NONL_S%NITYP(NT)
             NI = NI + 1
             ! this routine replaces the CPROJ of W1 with the derivatives of the projector
             CALL PROJ_DER(NONL_S, WDES1, W1_LOCAL, LATT_CUR, NI, I)
             IF (NI_LOCAL(NI, WDES1%COMM_INB) == 0) CYCLE
             DO ISPINOR = 1, WDES1%NRSPINORS
                DO LM = 1, NONL_S%LMMAX(NT)
                   DPROJECTOR(LMBASE+LM,I,ISPINOR) = W1_LOCAL%CPROJ(LMBASE+LM+(ISPINOR-1)*WDES1%NPRO/2)
                ENDDO
             ENDDO
             LMBASE = LMBASE + NONL_S%LMMAX(NT)
          ENDDO
       ENDDO

       IF  (LMBASE .NE. WDES1%NPRO/WDES1%NRSPINORS) THEN
          CALL vtutor%bug('error in N-PRO. expected: '//str(LMBASE)//' found: '//str(WDES1%NPRO/WDES1%NRSPINORS),__FILE__,__LINE__)
       ENDIF
    ENDDO

  END SUBROUTINE get_DPROJ


  SUBROUTINE getwave_write_inwap_yaml(NUNIT, T_INFO, WDES, DER, LATT_CUR, &
       P, GRID, GRIDC, KPOINTS, EFERMI)

    USE POSCAR
    USE wave
    USE LATTICE
    USE pseudo_struct_def
    USE PSEUDO
    USE mgrid

    INTEGER, intent(in) :: NUNIT, DER
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(wavedes), intent(in) :: WDES
    TYPE(latt), intent(in) :: LATT_CUR
    TYPE(grid_3d), intent(in) :: GRID, GRIDC
    REAL(q), intent(in) :: KPOINTS(1:3, 1:WDES%NKPTS)
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    REAL(q), intent(in) ::  EFERMI

    integer :: I, II, L, LL, M, NI, NT, NIS
    TYPE (potcar), POINTER :: PP

    OPEN(UNIT=NUNIT, FILE='inwap.yaml', STATUS='REPLACE')

    write(NUNIT, '(a)') 'cell:'
    write(NUNIT, '(a, i8)') '  NIONS: ', T_INFO%NIONS
    write(NUNIT, '(a, i8)') '  NKPTS: ', WDES%NKPTS
    write(NUNIT, '(a, i8)') '  NBTOT: ', WDES%NB_TOT
    write(NUNIT, '(a, i8)') '  ISPIN: ', WDES%ISPIN
    write(NUNIT, '(a, i4)') '  NRSPINORS: ', WDES%NRSPINORS
    write(NUNIT, '(a, i8)') '  NCDIJ: ', WDES%NCDIJ
    write(NUNIT, '(a, i10)') '  DER: ', DER
    write(NUNIT, '(a, F22.14)') '  EFERMI: ',EFERMI
    ! Lattice vectors in real space correspond to column vectors of A. We output
    ! the transpose.
    write(NUNIT, '(a)') '  lattice:'
    write(NUNIT,'("  - [", F20.10, ", ", F20.10, ",", F20.10, " ] # a")') LATT_CUR%A(:, 1)
    write(NUNIT,'("  - [", F20.10, ", ", F20.10, ",", F20.10, " ] # b")') LATT_CUR%A(:, 2)
    write(NUNIT,'("  - [", F20.10, ", ", F20.10, ",", F20.10, " ] # c")') LATT_CUR%A(:, 3)
    ! Lattice vectors in reciprocal space correspond to column vectors of B. We
    ! output the transpose
    write(NUNIT,'(a)') '  rlattice:'
    write(NUNIT,'("  - [", F20.10, ", ", F20.10, ",", F20.10, " ] # a*")') LATT_CUR%B(:, 1)
    write(NUNIT,'("  - [", F20.10, ", ", F20.10, ",", F20.10, " ] # b*")') LATT_CUR%B(:, 2)
    write(NUNIT,'("  - [", F20.10, ", ", F20.10, ",", F20.10, " ] # c*")') LATT_CUR%B(:, 3)
    write(NUNIT, '(a)') '  positions:'
    DO I = 1, T_INFO%NIONS
       write(NUNIT,'("  - [", F20.14, ", ", F20.14, ",", F20.14, " ]")') &
            T_INFO%POSION(:, I)
    END DO
    write(NUNIT, '(a)') 'PAW:'
    write(NUNIT, '(a, i3)') '  LMDIM: ', P(1)%LMDIM
    write(NUNIT, '(a)', advance='no') "  LMMAX: [ "
    NIS = 1
    DO NT = 1, T_INFO%NTYP
       DO NI = NIS, T_INFO%NITYP(NT) + NIS - 1
          PP=>PP_POINTER(P, NI, NT)
          IF (NI == T_INFO%NIONS) THEN
             write(NUNIT, '(i3, " ]")') PP%LMMAX
          ELSE
             write(NUNIT, '(i3, ", ")', advance='no') PP%LMMAX
          END IF
       END DO
       NIS = NIS + T_INFO%NITYP(NT)
    END DO
    write(NUNIT, '(a, i3)') '  LDIM:  ', P(1)%LDIM
    write(NUNIT, '(a)', advance='no') "  LMAX:  [ "
    NIS = 1
    DO NT = 1, T_INFO%NTYP
       DO NI = NIS, T_INFO%NITYP(NT) + NIS - 1
          PP=>PP_POINTER(P, NI, NT)
          IF (NI == T_INFO%NIONS) THEN
             write(NUNIT, '(i3, " ]")') PP%LMAX
          ELSE
             write(NUNIT, '(i3, ", ")', advance='no') PP%LMAX
          END IF
       END DO
       NIS = NIS + T_INFO%NITYP(NT)
    END DO
    ! write(NUNIT, '(a)') '  LPS:'
    ! DO NT = 1, T_INFO%NTYP
    !    DO NI = NIS, T_INFO%NITYP(NT) + NIS - 1
    !       write(NUNIT, '(a)', advance='no') "  - [ "
    !       PP=>PP_POINTER(P, NI, NT)
    !       DO L = 1, PP%LMAX
    !          IF (L == PP%LMAX) THEN
    !             write(NUNIT, '(i3, " ]")') PP%LPS(L)
    !          ELSE
    !             write(NUNIT, '(i3, ", ")', advance='no') PP%LPS(L)
    !          END IF
    !       END DO
    !    END DO
    !    NIS = NIS + T_INFO%NITYP(NT)
    ! END DO
    write(NUNIT, '(a)') '  lm_orbitals:'
    NI = 0
    DO NT = 1, T_INFO%NTYP
       DO II = 1, T_INFO%NITYP(NT)
          NI = NI + 1
          PP=>PP_POINTER(P, NI, NT)
          write(NUNIT,'(a, i7)') '  - atom_index: ', NI
          write(NUNIT,'(a, i3)') '    num_l_channels: ', PP%LMAX
          write(NUNIT,'(a)') '    channels:'
          DO L = 1, PP%LMAX
             LL = PP%LPS(L)
             write(NUNIT,'(a, i3)') '    - l : ', LL
             write(NUNIT,'(a)', advance='no') '      m : [ '
             DO M = -LL, LL
                IF (M == LL) THEN
                   write(NUNIT, '(i3, " ]")') M
                ELSE
                   write(NUNIT, '(i3, ", ")', advance='no') M
                END IF
             END DO
          ENDDO
       ENDDO
    ENDDO
    write(NUNIT, '(a)') 'fft_grid:'
    write(NUNIT, '(a, "[ ", i8, ", ", i8, ", ", i8, " ]")') &
         '  coarse: ', GRID%NGX, GRID%NGY, GRID%NGZ
    write(NUNIT, '(a, "[ ", i8, ", ", i8, ", ", i8, " ]")') &
         '  fine:   ', GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ
    write(NUNIT, '(a)') 'KPOINTS:'
    DO I = 1, WDES%NKPTS
       write(NUNIT,'("- [", F20.14, ", ", F20.14, ",", F20.14, " ]")') &
            KPOINTS(:, I)
    END DO
    CLOSE(NUNIT)

  END SUBROUTINE getwave_write_inwap_yaml


  SUBROUTINE getwave_write_data_in_text(NUNIT, T_INFO, WDES, &
       LMDIM, GRID, GRIDC, CCDIJ, CCQIJ, POT_SV, POT, &
       DPROJECTORS, WAVES)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                Write in text                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE base
    USE POSCAR
    USE wave
    USE mgrid

    INTEGER, intent(in) :: NUNIT, LMDIM
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(wavedes), intent(in) :: WDES
    TYPE(grid_3d), intent(in) :: GRID, GRIDC
    COMPLEX(8), intent(in):: CCDIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ)
    COMPLEX(8), intent(in):: CCQIJ(LMDIM, LMDIM, T_INFO%NIONS, WDES%NCDIJ)
    COMPLEX(8), intent(in) :: POT_SV(GRID%NGX, GRID%NGY, GRID%NGZ, WDES%NCDIJ)
    COMPLEX(8), intent(in) :: POT(GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ, WDES%NCDIJ)
    COMPLEX(8), intent(in) :: DPROJECTORS(0:3,1:LMDIM,1:T_INFO%NIONS,&
         1:WDES%NRSPINORS,1:WDES%ISPIN,1:WDES%NB_TOT,  1:WDES%NKPTS)
    COMPLEX(8), intent(in) :: WAVES(GRID%NGX, GRID%NGY, GRID%NGZ,&
         WDES%NRSPINORS, WDES%ISPIN, WDES%NB_TOT, WDES%NKPTS)
    integer :: I, II, J, NI, ISP, NX, NY, NZ, NK, NB, LM
    character(len=80):: string1, formatw

    write (string1,'(I0)') LMDIM
    formatw='('//trim(string1)//'(E15.7,1X))'

    OPEN(UNIT=NUNIT, FILE='PAW-OVERLAP', STATUS='REPLACE')
    DO ISP = 1, WDES%NCDIJ
       DO NI = 1, T_INFO%NIONS
          write(NUNIT, *) 'SPIN= ', ISP, ' ION= ', NI
          DO I = 1, LMDIM
             write(NUNIT, formatw) (REAL(CCQIJ(I, J, NI, ISP)), J = 1, LMDIM)
             write(NUNIT, formatw) (AIMAG(CCQIJ(I, J, NI, ISP)), J = 1, LMDIM)
          ENDDO
          write(NUNIT, *)
       ENDDO
    ENDDO
    CLOSE(NUNIT)

    OPEN(UNIT=NUNIT, FILE='PAW-STRENGTH', STATUS='REPLACE')
    DO ISP = 1, WDES%NCDIJ
       DO NI = 1, T_INFO%NIONS
          write(NUNIT, *) 'SPIN= ', ISP, ' ION= ', NI
          DO I = 1, LMDIM
             write(NUNIT, formatw) (REAL(CCDIJ(I, J, NI, ISP)), J = 1, LMDIM)
             write(NUNIT, formatw) (AIMAG(CCDIJ(I, J, NI, ISP)), J = 1, LMDIM)
          ENDDO
          write(NUNIT, *)
       ENDDO
    ENDDO
    CLOSE(NUNIT)

    OPEN(UNIT=NUNIT, FILE='LOCAL-POTENTIAL-SV', STATUS='REPLACE')
    write(NUNIT,'(3I5)') GRID%NGX, GRID%NGY, GRID%NGZ
    DO ISP = 1,WDES%NCDIJ
       write(NUNIT,'(5(1X,E17.11))') (((REAL(POT_SV(NX, NY, NZ, ISP)),&
            NX = 1, GRID%NGX), NY = 1, GRID%NGY), NZ=1, GRID%NGZ)
       write(NUNIT, *)
    ENDDO
    CLOSE(NUNIT)

    OPEN(UNIT=NUNIT, FILE='LOCAL-POTENTIAL', STATUS='REPLACE')
    write(NUNIT,'(3I5)') GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ
    DO ISP = 1,WDES%NCDIJ
       write(NUNIT,'(5(1X,E17.11))') (((REAL(POT(NX, NY, NZ, ISP)),&
            NX = 1, GRIDC%NGX), NY = 1, GRIDC%NGY), NZ=1, GRIDC%NGZ)
       write(NUNIT, *)
    ENDDO
    CLOSE(NUNIT)

    OPEN(UNIT=NUNIT, FILE='DPROJECTORS', STATUS='REPLACE')
    DO ISP = 1, WDES%ISPIN
       DO NK = 1, WDES%NKPTS
          DO NB = 1, WDES%NB_TOT
             DO NI = 1, T_INFO%NIONS
                DO II = 0, WDES%NRSPINORS - 1
                   DO LM = 1, LMDIM
                      write(NUNIT, '(8F12.6)')&
                           REAL(DPROJECTORS(0, LM,NI,II + 1,ISP,NB, NK)),&
                           AIMAG(DPROJECTORS(0,LM,NI,II + 1,ISP,NB, NK)),&
                           REAL(DPROJECTORS(1, LM,NI,II + 1,ISP,NB, NK)),&
                           AIMAG(DPROJECTORS(1,LM,NI,II + 1,ISP,NB, NK)),&
                           REAL(DPROJECTORS(2, LM,NI,II + 1,ISP,NB, NK)),&
                           AIMAG(DPROJECTORS(2,LM,NI,II + 1,ISP,NB, NK)),&
                           REAL(DPROJECTORS(3, LM,NI,II + 1,ISP,NB, NK)),&
                           AIMAG(DPROJECTORS(3,LM,NI,II + 1,ISP,NB, NK))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CLOSE(NUNIT)

    OPEN(UNIT=NUNIT, FILE='WAVES', STATUS='REPLACE')
    DO ISP = 1, WDES%ISPIN
       DO NK = 1, WDES%NKPTS
          DO NB = 1, WDES%NB_TOT
             DO II = 0, WDES%NRSPINORS - 1
                write(NUNIT,*) 'WAVE= ', ISP, NK, NB, II
                DO NZ = 1, GRID%NGZ
                   DO NY = 1, GRID%NGY
                      DO NX = 1, GRID%NGX
                         write(NUNIT,'(E17.11,1X,E17.11)') &
                              WAVES(NX, NY, NZ, II + 1, ISP, NB, NK)
                      ENDDO
                   ENDDO
                ENDDO
                write(NUNIT, *)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CLOSE(NUNIT)

  END SUBROUTINE getwave_write_data_in_text


  SUBROUTINE getwave_RIJ(NCDIJ,LMDIM,P, T_INFO, CCRIJ)
    ! Compute R0 and R1 in Eqs D3 and 4 of PRB 100, 174304 (2019)
    !> CCRIJ(LMP,LM,NI,ISP,1) = R^1(n l, np lp)-\tilde{R}^1(n l, np lp)
    !>                        = \int dr r^2 R_{nl}(r) \frac{ R_{np lp}(r)}{r}
    !>                       = \int dr r^2 \tilde{R}_{nl}(r) \frac{ \tilde{R}_{np lp}(r)}{r}
    !>                        = \int dr \phi_L(r;NI) \phi_{LP}(r;NI) / r
    !>                        - \int dr \tilde{\phi}_L(r;NI) \tilde{\phi}_{LP}(r;NI) / r
    !>
    !> CCRIJ(LMP,LM,NI,ISP,2) = R^0(n l, np lp)-\tilde{R}^0(n l, np lp)
    !>                        = \int dr r^2 R_{nl}(r) \frac{d R_{np lp}(r)}{dr}
    !>                        - \int dr r^2 \tilde{R}_{nl}(r) \frac{d \tilde{R}_{np lp}(r)}{dr}
    !>                        = \int dr \phi_L(r;NI) \Big( \frac{d \phi_{LP}(r;NI)}{dr} - \phi_{LP}(r;NI) / r \Big)
    !>                        - \int dr \tilde{\phi}_L(r;NI) \Big( \frac{d \tilde{\phi}_{LP}(r;NI)}{dr}
    !>                                                             - \tilde{\phi}_{LP}(r;NI) / r \Big)
    !>
    !> with LM = (n l m) and LMP=(np lp mp), and R_{nl}(r) the radial part of atomic WF.
    !> then R_{nl}(r)= \phi_{nl}(r)/r, and \phi = PP%WAE(I,L) or PP%WPS
    USE prec
    USE PSEUDO
    USE POSCAR
    USE radial
    USE wave

    INTEGER, INTENT(IN) :: NCDIJ
    INTEGER, INTENT(IN) :: LMDIM
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    COMPLEX(8), intent(out):: CCRIJ(LMDIM, LMDIM, T_INFO%NIONS, NCDIJ, 2)

    INTEGER:: NIS, NT, NI, L, LP, LDIM, NMAX, I, M, MP, LM, LMP, ISP, LL, LLP
    TYPE (potcar), POINTER :: PP
    real(8), allocatable :: R(:)
    real(8), allocatable :: TMP(:)
    real(8) H
    REAL(8) val

    LDIM = P(1)%LDIM

    DO ISP = 1, NCDIJ

       NIS = 1
       DO NT = 1, T_INFO%NTYP
          DO NI = NIS, T_INFO%NITYP(NT) + NIS - 1

             PP=>PP_POINTER(P, NI, NT)
             NMAX=PP%R%NMAX
             H=PP%R%H

             allocate(R(1:NMAX))
             allocate(TMP(1:NMAX))

             LM=0
             DO L = 1, PP%LMAX
                LL=PP%LPS(L)
                DO M=-LL,LL
                   LM=LM+1

                   LMP=0
                   DO LP = 1, PP%LMAX
                      LLP=PP%LPS(LP)
                      DO MP=-LLP,LLP
                         LMP=LMP+1

                         ! first radial integral for wave function
                         DO I=1,NMAX
                            R(I)=dble(PP%R%R(I)) ! radial grid, R(I)=RSTART*EXP((I-1)*H), in Angstroms
                            IF (PP%R%R(I) .le. PP%RDEP) THEN
                               TMP(I)=dble((PP%WAE(I,L)*PP%WAE(I,LP))/(PP%R%R(I)))
                               !!TMP(I)=dble((PP%WAE(I,L)*PP%WAE(I,LP)))     !  using this instead we obtain QIJ when LM=LMP
                            ELSE
                               TMP(I)=0.0d0
                            ENDIF
                         ENDDO
                         val=int_l(R,H,TMP)
                         CCRIJ(LMP,LM,NI,ISP,1)=dcmplx(val)

                         ! second radial integral for wave function
                         DO I=1,NMAX
                            TMP(I)=dble(PP%WAE(I,LP))
                         ENDDO
                         call get_grad_l(R,H,TMP)


                         DO I=1,NMAX
                            IF (PP%R%R(I) .le. PP%RDEP) THEN
                               TMP(I)=dble(PP%WAE(I,L))*(TMP(I)-dble(PP%WAE(I,LP)/R(I)))
                            ELSE
                               TMP(I)=0.0d0
                            ENDIF
                         ENDDO
                         val=int_l(R,H,TMP)
                         CCRIJ(LMP,LM,NI,ISP,2)=dcmplx(val)


                         !                         ! first radial integral for pseudo wave function
                         DO I=1,NMAX
                            R(I)=dble(PP%R%R(I))
                            IF (PP%R%R(I) .le. PP%RDEP) THEN
                               TMP(I)=dble((PP%WPS(I,L)*PP%WPS(I,LP))/(PP%R%R(I)))
                               !TMP(I)=dble((PP%WPS(I,L)*PP%WPS(I,LP)))     !  using this instead we obtain QIJ when LM=LMP
                            ELSE
                               TMP(I)=0.0d0
                            ENDIF
                         ENDDO
                         val=int_l(R,H,TMP)
                         CCRIJ(LMP,LM,NI,ISP,1)=CCRIJ(LMP,LM,NI,ISP,1)-dcmplx(val)

                         ! second radial integral for pseudo wave function
                         DO I=1,NMAX
                            TMP(I)=dble(PP%WPS(I,LP))
                         ENDDO
                         call get_grad_l(R,h,TMP)

                         DO I=1,NMAX
                            IF (PP%R%R(I) .le. PP%RDEP) THEN
                               TMP(I)=dble(PP%WPS(I,L))*(TMP(I)-dble(PP%WPS(I,LP)/R(I)))
                            ELSE
                               TMP(I)=0.0d0
                            ENDIF
                         ENDDO
                         val=int_l(R,H,TMP)
                         CCRIJ(LMP,LM,NI,ISP,2)=CCRIJ(LMP,LM,NI,ISP,2)-dcmplx(val)

                      END DO
                   END DO
                END DO
             END DO

             deallocate(TMP)
             deallocate(R)

          END DO
          NIS = NIS + T_INFO%NITYP(NT)
       END DO

    END DO

  END SUBROUTINE getwave_RIJ


  SUBROUTINE get_grad_l(R,h,W)
    ! compute the gradient, to order 4, of the function in W, on the logarithmic grid in R with exponent h
    real(8), allocatable, intent(in) :: R(:)
    real(8), intent(in) :: h
    real(8), allocatable, intent(inout) :: W(:)
    real(8), allocatable :: DW(:)
    integer i,N

    N=size(R(:))
    IF (N==0) return
    allocate(DW(1:N))

    ! forward difference for the first 2 points
    DW(1)=-(25.0d0)/(12.0d0)*W(1)+4.0d0*W(2)-3.0d0*W(3)+(4.0d0)/(3.0d0)*W(4)-(1.0d0)/(4.0d0)*W(5)
    DW(2)=-(25.0d0)/(12.0d0)*W(2)+4.0d0*W(3)-3.0d0*W(4)+(4.0d0)/(3.0d0)*W(5)-(1.0d0)/(4.0d0)*W(6)

    ! backward difference for the last 2 points
    DW(N-1)=(25.0d0)/(12.0d0)*W(N-1)-4.0d0*W(N-2)+3.0d0*W(N-3)-(4.0d0)/(3.0d0)*W(N-4)+(1.0d0)/(4.0d0)*W(N-5)
    DW(N)=(25.0d0)/(12.0d0)*W(N)-4.0d0*W(N-1)-3.0d0*W(N-2)-(4.0d0)/(3.0d0)*W(N-3)+(1.0d0)/(4.0d0)*W(N-4)

    ! central difference for other points
    do i=3,N-2
       DW(i)=(1.0d0)/(12.0d0)*W(i-2)-(2.0d0)/(3.0d0)*W(i-1)+(2.0d0)/(3.0d0)*W(i+1)-(1.0d0)/(12.0d0)*W(i+2)
    enddo

    ! log grid
    do i=1,N
       W(i)=DW(i)/(R(i)*h)
    enddo

    deallocate(DW)
  END SUBROUTINE get_grad_l

  FUNCTION int_l(R,h,W)
    ! return the integral of the function W on the logarithmic grid R with exponent h
    real(8) :: int_l
    real(8), allocatable, intent(in) :: R(:)
    real(8), intent(in) :: h
    real(8), allocatable, intent(in) :: W(:)
    integer N,i
    real(8) v

    N=size(R(:))

    ! Simpson require N to be odd.
    ! Therefore we start the integration from the upper limit,
    ! and forgot the last innner point if N is even.

    v=0.0d0
    do i=N,3,-2
       v=v+R(i)*W(i)+4.0d0*R(i-1)*W(i-1)+R(i-2)*W(i-2)
    enddo

    int_l=v*h/(3.0d0)

  END FUNCTION int_l


  SUBROUTINE get_CCRIJ_A(NCDIJ,LMDIM,T_INFO,CCRIJ,P,CCRIJ_A)
    !> Compute RIJ matrix
    !>  CCRIJ_A(LMP,LM,NI,ISP,DIR) = R_{LM LMP}(NI)
    !>                             = - \Big(  \int d^3 \phi_{n l m}^*(r;NI) \nabla \phi_{np lp mp}(r;NI)
    !>                                       -\int d^3 \tilde{\phi}_{n l m}^*(r;NI) \nabla \tilde{\phi}_{np lp mp}(r;NI) \Big)
    !> with LM=(nlm) and LMP =(np lp mp)
    !> It is computed using formula D2 of PRB 100, 174304 (2019)
    USE prec
    USE PSEUDO
    USE POSCAR
    USE radial
    USE wave

    INTEGER, INTENT(IN) :: NCDIJ
    INTEGER, INTENT(IN) :: LMDIM
    TYPE(type_info), intent(in) :: T_INFO
    COMPLEX(8), intent(in):: CCRIJ(LMDIM, LMDIM, T_INFO%NIONS, NCDIJ, 2)
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    COMPLEX(8), intent(out):: CCRIJ_A(LMDIM, LMDIM, T_INFO%NIONS, 3, NCDIJ) ! RIJ coefficients with angular part
    REAL(8) Nlm,A(1:3,1:6),B(1:3,1:6),C(1:3,1:6),T0(1:3,1:3),T1(1:3,1:3)
    REAL(8) v0(3,1),v1(3,1),v2(3,1),v3(3,1)
    REAL(8) Norm,NormP,Plm(6,1)

    INTEGER:: NIS, NT, NI, L, LP, LDIM, NMAX, I, M, MP, LM, LMP, ISP, LL, LLP, DIR
    TYPE (potcar), POINTER :: PP
    integer i1,i2,i3

    LDIM = P(1)%LDIM

    DO ISP = 1, NCDIJ

       NIS = 1
       DO NT = 1, T_INFO%NTYP
          DO NI = NIS, T_INFO%NITYP(NT) + NIS - 1

             PP=>PP_POINTER(P, NI, NT)

             LM=0
             DO L = 1, PP%LMAX
                LL=PP%LPS(L)
                DO M=-LL,LL
                   LM=LM+1
                   call get_NYlm(LL,abs(M),Norm)

                   LMP=0
                   DO LP = 1, PP%LMAX
                      LLP=PP%LPS(LP)
                      DO MP=-LLP,LLP
                         LMP=LMP+1

                         call get_NYlm(LLP,abs(MP),NormP)


                         call get_ABC(LLP,abs(MP),A,B,C)

                         call get_T01(M,MP,T0,T1)


                         call get_vec_P(LL,abs(M),LLP,abs(MP),Plm)

                         v0=matmul(A,Plm)
                         v1=matmul(T0,v0)

                         v0=matmul(B,Plm)
                         v2=matmul(T0,v0)

                         v0=matmul(C,Plm)
                         v3=matmul(T1,v0)


                         DO DIR=1,3
                            CCRIJ_A(LMP,LM,NI,DIR,ISP)=CCRIJ(LMP,LM,NI,ISP,2)*v1(DIR,1)*((-Norm*NormP)/((2.0d0*LLP+1.0d0))) &
                                 +CCRIJ(LMP,LM,NI,ISP,1)*v2(DIR,1)*((Norm*NormP)/((4.0d0*LLP+2.0d0))) &
                                 +CCRIJ(LMP,LM,NI,ISP,1)*v3(DIR,1)*((Norm*NormP)/(2.0d0))
                         ENDDO

                      ENDDO
                   ENDDO

                ENDDO
             ENDDO

          ENDDO
          NIS = NIS + T_INFO%NITYP(NT)
       ENDDO

    ENDDO
  END SUBROUTINE get_CCRIJ_A


  SUBROUTINE check_RIJ_A(LMAX)
    integer, intent(in) :: LMAX
    integer l1,m1,l2,m2
    REAL(8) N1,N2
    REAL(8) A(1:3,1:6),B(1:3,1:6),C(1:3,1:6),T0(1:3,1:3),T1(1:3,1:3),Plm(6,1)
    REAL(8) v0(3,1),v1(3,1),v2(3,1),v3(3,1)

    OPEN(UNIT=48, FILE='check_RIJ', STATUS='REPLACE')

    do l1=0,LMAX
       do m1=-l1,l1
          call get_NYlm(l1,abs(m1),N1)
          do l2=0,LMAX
             do m2=-l2,l2
                call get_NYlm(l2,abs(m2),N2)

                call get_ABC(l2,abs(m2),A,B,C)
                call get_T01(m1,m2,T0,T1)
                call get_vec_P(l1,abs(m1),l2,abs(m2),Plm)

                v0=matmul(A,Plm)
                v1=matmul(T0,v0)
                v1=v1*N1*N2/(2.0d0*l2+1.0d0)

                v0=matmul(B,Plm)
                v2=matmul(T0,v0)
                v2=v2*N1*N2/(-4.0d0*l2-2.0d0)

                v0=matmul(C,Plm)
                v3=matmul(T1,v0)
                v3=v3*N1*N2*(-0.5d0)

                write(48,'(4I3,1X,9F8.4)') l1,m1,l2,m2,v1,v2,v3

             enddo
          enddo
       enddo
    enddo
    close(48)
  END SUBROUTINE check_RIJ_A


  SUBROUTINE get_NYlm(l,m,Nlm)
    ! define normaliation of complex spherical harmonics: Ylm=(-1)^m N_lm P_lm(cos theta) ex^(i m phi)
    ! xxx Nlm= \sqrt{\frac{2l+1}{4 \pi} \frac{(l-m)!}{(l+m)!}}
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(8), intent(out) :: Nlm
    real(8) dl,dm,PI

    PI=acos(-1.0d0)

    dl=dble(l)
    dm=dble(m)
    Nlm=sqrt((2.0d0*dl+1.0d0)/(4.0d0*PI))*exp(0.5d0*LOG_GAMMA(dl-dm+1.0d0)-0.5d0*LOG_GAMMA(dl+dm+1.0d0))

  END SUBROUTINE get_NYlm


  SUBROUTINE get_ABC(l,m,A,B,C)
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(8), intent(out) :: A(1:3,1:6)
    real(8), intent(out) :: B(1:3,1:6)
    real(8), intent(out) :: C(1:3,1:6)
    real(8) dl,dm

    dl=dble(l); dm=dble(m)

    ! Matrix A
    A(1:3,1:6)=0.0d0
    A(1,3)=-1.0d0
    A(1,6)=1.0d0
    A(2,3)=-1.0d0
    A(2,6)=1.0d0
    A(3,2)=dl+dm
    A(3,5)=dl-dm+1.0d0

    ! Matrix B
    B(1:3,1:6)=0.0d0
    B(1,1)=-(dl+dm)*(dl-dm+1.0d0)*(dl+dm-1.0d0)
    B(1,3)=dl+dm+1.0d0
    B(1,4)=-(dl+dm)*(dl-dm+1.0d0)*(dl-dm+2.0d0)
    B(1,6)=dl-dm
    B(2,1)=-(dl+dm)*(dl-dm+1.0d0)*(dl+dm-1.0d0)
    B(2,3)=dl+dm+1.0d0
    B(2,4)=-(dl+dm)*(dl-dm+1.0d0)*(dl-dm+2.0d0)
    B(2,6)=dl-dm
    B(3,2)=-(2.0d0*dl+2.0d0)*(dl+dm)
    B(3,5)=2.0d0*dl*(dl-dm+1.0d0)

    ! Matrix C
    C(1:3,1:6)=0.0d0
    C(1,4)=(dl-dm+1.0d0)*(dl-dm+2.0d0)
    C(1,6)=1.0d0
    C(2,4)=-(dl-dm+1.0d0)*(dl-dm+2.0d0)
    C(2,6)=-1.0d0
  END SUBROUTINE get_ABC


  SUBROUTINE get_T01(m1,m2,T0,T1)
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    real(8), intent(out) :: T0(1:3,1:3)
    real(8), intent(out) :: T1(1:3,1:3)
    real(8) PI

    PI=acos(-1.0d0)

    ! T0 matrix
    T0(1:3,1:3)=0.0d0
    if ((m1 .eq. 0) .and. (m2 .eq. 1)  ) then
       T0(1,1)=sqrt(2.0d0)*PI
       goto 42
    elseif((m1 .eq. 1) .and. (m2 .eq. 0)  ) then
       T0(1,1)=sqrt(2.0d0)*PI
       goto 42
    elseif((m1 .eq. 0) .and. (m2 .eq. -1)  ) then
       T0(2,2)=sqrt(2.0d0)*PI
       goto 42
    elseif((m1 .eq. -1) .and. (m2 .eq. 0)  ) then
       T0(2,2)=sqrt(2.0d0)*PI
       goto 42
    elseif((m1*m2 .gt. 0) .and. (abs(m1-m2) .eq. 1)  ) then
       T0(1,1)=PI
       goto 42
    elseif((m1*m2 .lt. 0) .and. (abs(m1+m2) .eq. 1)  ) then
       T0(2,2)=PI*sign(1.0d0,-1.0d0*(m1+m2))
       goto 42
    elseif (m1 .eq. m2) then
       T0(3,3)=2.0d0*PI
       goto 42
    endif
42  continue

    ! T1 matrix
    T1(1:3,1:3)=0.0d0
    if ((m1 .eq. 0) .and. (m2 .eq. 1)  ) then
       T1(1,1)=-sqrt(2.0d0)*PI
       goto 43
    elseif((m1 .eq. 0) .and. (m2 .eq. -1)  ) then
       T1(2,2)=sqrt(2.0d0)*PI
       goto 43
    elseif((m1*m2 .gt. 0) .and. (abs(m1-m2) .eq. 1)  ) then
       T1(1,1)=PI*sign(1.0d0,1.0d0*m1)*sign(1.0d0,1.0d0*(m1-m2))
       goto 43
    elseif((m1*m2 .lt. 0) .and. (abs(m1+m2) .eq. 1)  ) then
       T1(2,2)=PI*sign(1.0d0,1.0d0*m1)
       goto 43
    endif
43  continue
  END SUBROUTINE get_T01


  FUNCTION  Wigner3j(l1,l2,l3,m1,m2,m3)
    ! returm Wigner coefficient
    real(8) :: Wigner3j
    integer, intent(in) :: l1
    integer, intent(in) :: l2
    integer, intent(in) :: l3
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real(8) dl1,dl2,dl3,dm1,dm2,dm3,dk
    real(8) w
    integer sup_s(3),inf_s(3),minS,maxS
    integer k

    Wigner3j=0.0d0


    if (l1 .lt. 0) then
       CALL vtutor%error('Error-1a in Wigner3j')
    endif
    if (l2 .lt. 0) then
       CALL vtutor%error('Error-2 in Wigner3j')
    endif
    if (l3 .lt. 0) then
       CALL vtutor%error('Error-3 in Wigner3j')
    endif

    if ( (l1+l2-l3) .lt. 0)  goto 22
    if ( (l1-l2+l3) .lt. 0)  goto 22
    if ( (-l1+l2+l3) .lt. 0) goto 22

    if (abs(m1) .gt. l1) goto 22
    if (abs(m2) .gt. l2) goto 22
    if (abs(m3) .gt. l3) goto 22

    if (m1+m2+m3 .ne. 0) goto 22
    if ((abs(m1)+abs(m2)+abs(m3) .eq. 0) .and. (mod(abs(l1+l2+l3),2) .ne. 0)) goto 22

    ! sum lower limit
    inf_s(1)=-l3+l2-m1
    inf_s(2)=-l3+l1+m2
    inf_s(3)=0
    minS=maxval(inf_s)
    if (minS .lt. 0) then
       CALL vtutor%error('Error-4 in Wigner3j')
    endif
    ! sum upper limit
    sup_s(1)=l1+l2-l3
    sup_s(2)=l1-m1
    sup_s(3)=l2+m2
    maxS=minval(sup_s)
    if (maxS .lt. 0) then
       CALL vtutor%error('Error-5 in Wigner3j')
    endif

    dl1=dble(l1);dl2=dble(l2);dl3=dble(l3)
    dm1=dble(m1);dm2=dble(m2);dm3=dble(m3)

    w=0.0d0
    do k=minS,maxS
       dk=dble(k)
       w=w+((-1.0d0)**(k+l1-l2-m3))*exp( &
            -1.0d0*(LOG_GAMMA(dk+1.0d0) &
            +LOG_GAMMA(dl1+dl2-dl3-dk+1.0d0) &
            +LOG_GAMMA(dl1-dm1-dk+1.0d0) &
            +LOG_GAMMA(dl2+dm2-dk+1.0d0) &
            +LOG_GAMMA(dl3-dl2+dm1+dk+1.0d0)&
            +LOG_GAMMA(dl3-dl1-dm2+dk+1.0d0)) &
            + 0.5d0*(LOG_GAMMA(dl1+dl2-dl3+1.0d0) &
            +LOG_GAMMA(dl1-dl2+dl3+1.0d0)&
            +LOG_GAMMA(-dl1+dl2+dl3+1.0d0) &
            +LOG_GAMMA(dl1+dm1+1.0d0)&
            +LOG_GAMMA(dl1-dm1+1.0d0)&
            +LOG_GAMMA(dl2+dm2+1.0d0)&
            +LOG_GAMMA(dl2-dm2+1.0d0)&
            +LOG_GAMMA(dl3+dm3+1.0d0)&
            +LOG_GAMMA(dl3-dm3+1.0d0)&
            -LOG_GAMMA(dl1+dl2+dl3+2.0d0) ))
    enddo

    Wigner3j=w

22  continue
    return
  END FUNCTION  Wigner3j


  SUBROUTINE get_Gaunt(l1,m1,l2,m2,G1,G2)
    ! return the first and second Gaunt coeffients in G1 and G2. They are vectors with l-components= |l1-l2|..l1+l2
    integer, intent(in) :: l1
    integer, intent(in) :: m1
    integer, intent(in) :: l2
    integer, intent(in) :: m2
    real(8), allocatable, intent(out) :: G1(:)
    real(8), allocatable, intent(out) :: G2(:)
    real(8) w1,w2
    integer i,j,ind
    real(8) dl1,dl2,dm1,dm2,dj
    real(8) pi

    pi=acos(-1.0d0)

    allocate(G1(1:l1+l2-abs(l1-l2)+1))
    allocate(G2(1:l1+l2-abs(l1-l2)+1))

    dl1=dble(l1);dl2=dble(l2);dm1=dble(m1);dm2=dble(m2)

    ind=0
    do j=abs(l1-l2),l1+l2
       ind=ind+1
       dj=dble(j)

       if (abs(m1+m2) .gt. j) then
          G1(ind)=0.0d0
          G2(ind)=0.0d0
       else
          w1=Wigner3j(l1,l2,j,0,0,0)
          w2=Wigner3j(l1,l2,j,m1,m2,-m1-m2)
          if (w1*w2 .eq. 0.0d0) then
             G1(ind)=0; G2(ind)=0
             goto 22
          endif
          G1(ind)=((-1.0d0)**(m1+m2))* &
               sqrt((2.0d0*dl1+1)*(2.0d0*dl2+1)*(2.0d0*dj+1)/(4.0d0*pi))* &
               w1*w2
          G2(ind)=((-1.0d0)**(m1+m2))*(2.0d0*dj+1.0d0)* &
               exp(0.5d0*(LOG_GAMMA(dj-dm1-dm2+1.0d0)+LOG_GAMMA(dl1+dm1+1.0d0)+ &
               LOG_GAMMA(dl2+dm2+1.0d0)-LOG_GAMMA(dj+dm1+dm2+1.0d0)- &
               LOG_GAMMA(dl1-dm1+1.0d0)-LOG_GAMMA(dl2-dm2+1.0d0)))* &
               w1*w2
       endif
22     continue
    enddo
  END SUBROUTINE get_Gaunt


  FUNCTION overlap_Plm(l1,m1,l2,m2)
    ! return the overlap of 2 associated Legendre functions
    real(8) :: overlap_Plm
    integer, intent(in) :: l1
    integer, intent(in) :: m1
    integer, intent(in) :: l2
    integer, intent(in) :: m2
    real(8), allocatable :: G1(:)
    real(8), allocatable :: G2(:)
    real(8), allocatable :: Alm(:)
    integer l,ind,m
    real(8) dl,dm,pi
    pi=acos(-1.0d0)

    call get_Gaunt(l1,m1,l2,m2,G1,G2)
    m=abs(m1+m2)
    dm=dble(m)

    allocate(Alm(1:l1+l2-abs(l1-l2)+1))
    Alm(1:l1+l2-abs(l1-l2)+1)=0.0d0

    ind=0
    do l=abs(l1-l2),l1+l2
       dl=dble(l)
       ind=ind+1
       if ( m .gt. l) then
          Alm(ind)=0.0d0
          goto 33
       elseif ( (l .eq. 0) .and. (m .eq. 0)) then
          Alm(ind)=2.0d0
          goto 32
       elseif((mod(l,2) .eq. 0) .and. (mod(m,2) .eq. 0 )) then
          Alm(ind)=((2.0d0*dm)/dl)* &
               exp(2.0d0*LOG_GAMMA(0.5d0*dl+1.0d0)+LOG_GAMMA(dl+dm+1.0d0)- &
               LOG_GAMMA(0.5d0*(dl-dm)+1.0d0)- &
               LOG_GAMMA(0.5d0*(dl+dm)+1.0d0)-LOG_GAMMA(dl+2.0d0))
          goto 32
       elseif((mod(l,2) .ne. 0) .and. (mod(m,2) .ne. 0 )) then
          Alm(ind)=((pi*dm)/(((2.0d0)**(2*l+1))*dl))* &
               exp(      LOG_GAMMA(dl+2.0d0)+LOG_GAMMA(dl+dm+1.0d0) - &
               2.0d0*LOG_GAMMA(0.5d0*(dl+1.0d0)+1.0d0) - &
               LOG_GAMMA(0.5d0*(dl-dm)+1.0d0) - &
               LOG_GAMMA(0.5d0*(dl+dm)+1.0d0))
          goto 32
       endif
32     continue
       if (m1+m2 .lt. 0) then
          Alm(ind)=Alm(ind)*((-1.0d0)**(m))* exp( LOG_GAMMA(dl-dm+1.0d0)-LOG_GAMMA(dl+dm+1.0d0))
       endif
33     continue
    enddo

    overlap_Plm=dot_product(G2,Alm)

    deallocate(Alm)
    deallocate(G2)
    deallocate(G1)
  END FUNCTION overlap_Plm


  SUBROUTINE get_vec_P(l1,m1,l2,m2,Plm)
    integer, intent(in) :: l1
    integer, intent(in) :: m1
    integer, intent(in) :: l2
    integer, intent(in) :: m2
    real(8), intent(out) :: Plm(1:6,1:1)

    Plm(1,1)=overlap_Plm(l1,m1,l2-1,m2-1)
    Plm(2,1)=overlap_Plm(l1,m1,l2-1,m2)
    Plm(3,1)=overlap_Plm(l1,m1,l2-1,m2+1)

    Plm(4,1)=overlap_Plm(l1,m1,l2+1,m2-1)
    Plm(5,1)=overlap_Plm(l1,m1,l2+1,m2)
    Plm(6,1)=overlap_Plm(l1,m1,l2+1,m2+1)

  END SUBROUTINE get_vec_P


  SUBROUTINE get_momentum(LATT_CUR,WA,LMDIM,CCRIJ_A,p_k)
    ! compute <WA|P|WA>/m in atomic units (velocity in first Bohr orbit, v_0).
    ! the chage in energy in eV is then E(k+dk)-E(k)=27.211 * 0.529 *  <WA|P|WA>/m * dk
    ! v_0 = alpha c = 2.18769 10^6 m/s 
    ! hbar / m  = (1.0545718 10^-34 m2 kg / s) / (9.10938356 10^-31 kg) = 0.11576 10^-3 m2/s 
    ! (hbar / m) / (alpha c) = a_0
    ! to get the output velocity in m/s, just multiply by v_0.
    USE lattice
    USE constant
    USE wave
    USE wave_high
    USE mpimy
    TYPE(latt), intent(in) :: LATT_CUR
    TYPE(wavefuna), INTENT(IN)  :: WA
    INTEGER, INTENT(IN) :: LMDIM
    COMPLEX(8), INTENT(IN) :: CCRIJ_A(:,:,:,:,:)
    COMPLEX(8), INTENT(OUT) :: p_k(:,:,:)

    TYPE (wavedes1), TARGET :: WDES1
    TYPE(wavefuna)   :: PWA          !  P|WA> wavefunction
    INTEGER I,J,DIR,SIZEOFM,NPL,ISPINOR
    REAL(8) hbar_over_m,p_c
    LOGICAL LGAMMA
    COMPLEX(8) Icmplx
    REAL(8), ALLOCATABLE :: p_G(:,:)

    Icmplx=cmplx(0.0d0,1.0d0)
    !hbar_over_m=0.115767636160320e+07/(0.218769126364e+07)    ! remove  denominator to get velocity in m/s
    hbar_over_m = AUTOA

    ! check for gamma version, somehow
    SIZEOFM=sizeof( WA%CPROJ(1,1))
    IF (SIZEOFM .EQ. 8) THEN
       LGAMMA = .TRUE.
    ELSEIF (SIZEOFM .EQ. 16) THEN
       LGAMMA = .FALSE.
    ELSE
       CALL vtutor%error('Type error in get_momentum')
    ENDIF

    ! Initialization
    p_k=dcmplx(0)

    ! The descriptor for this k-point
    WDES1 = WA%WDES1

    ! local G-basis*(spinor) size
    NPL=WDES1%NPL

    ! Create the P|WA> wavefunctions; it also point to the redistributed arry shape
    CALL NEWWAVA(PWA, WDES1, WDES1%NBANDS)

    if ( size(WA%CW(1,:)) .NE. WDES1%NBANDS) then
       CALL vtutor%error('Error allocating PWA: NBANDS')
    endif
    if ( size(WA%CW(:,1)) .NE. WDES1%NRPLWV) then
       CALL vtutor%error('Error allocating PWA: G*spinors')
    endif

    IF (LGAMMA) THEN
       ALLOCATE(P_G(WDES1%NB_TOT,WDES1%NB_TOT))
    ENDIF


    DO DIR=1,3  ! loop over cartesian components of P

       ! compute hbar/m (k+G)*a_{kn}(G) in atomic units (v_0=alpha c)
       DO ISPINOR=1,WDES1%NRSPINORS
       DO I=1,WDES1%NGVECTOR
          p_c=hbar_over_m*2.0d0*pi* &
               ((WDES1%IGX(I)+WDES1%VKPT(1))*LATT_CUR%B(DIR,1)+ &
                (WDES1%IGY(I)+WDES1%VKPT(2))*LATT_CUR%B(DIR,2)+ &
                (WDES1%IGZ(I)+WDES1%VKPT(3))*LATT_CUR%B(DIR,3))
          PWA%CW(I+(ISPINOR-1)*WDES1%NGVECTOR,:)=WA%CW(I+(ISPINOR-1)*WDES1%NGVECTOR,:)*p_c
       ENDDO
       ENDDO

       ! compute sum_J R_IJ <p_J|WA>.
       ! We could use OVERL to do that, but because of COMPLEX type I use my own routine sum_over_LM
       CALL sum_over_LM(WDES1,CCRIJ_A(:,:,:,DIR,:),WA,PWA)

       ! redistribute WA and PWA wavefunctions
       CALL REDISTRIBUTE_PW(WA)
       CALL REDISTRIBUTE_PROJ(WA)

       CALL REDISTRIBUTE_PW(PWA)
       CALL REDISTRIBUTE_PROJ(PWA)

       ! compute overlap and MPI_sum.
       ! Plane wave coefficients are always complex
       CALL ZGEMM('C','N',                                &
            WDES1%NB_TOT,                          &
            WDES1%NB_TOT,                          &
            WDES1%NPL_RED,                         &
            dcmplx(1.0d0),                         &
            WA%CW_RED(1,1), WDES1%NRPLWV_RED,      &
            PWA%CW_RED(1,1), WDES1%NRPLWV_RED,     &
            dcmplx(1.0d0),p_k(1,1,DIR),WDES1%NB_TOT)

       ! projections contribution
       IF (LGAMMA) THEN
          P_G=0.0d0
          CALL DGEMM('T','N',                            &
               WDES1%NB_TOT,                        &
               WDES1%NB_TOT,                        &
               WDES1%NPRO_O_RED,                    &    ! is zero WDES1%LOVERL= .FALSE.
               hbar_over_m,                         &
               WA%CPROJ_RED(1,1), WDES1%NPROD_RED,  &
               PWA%CPROJ_RED(1,1), WDES1%NPROD_RED, &
               0.0d0,p_G(1,1),WDES1%NB_TOT)

          DO I=1,WDES1%NB_TOT
             DO J=1,WDES1%NB_TOT
                p_k(I,J,DIR)=Icmplx*(aimag(p_k(I,J,DIR))+P_G(I,J))
                ! p_k should be multiplied by 2, but in vasp normalization is
                ! over the G>0 in the Gamma case, therefore no 2.
             ENDDO
          ENDDO
       ELSE
          CALL ZGEMM('C','N',                             &
               WDES1%NB_TOT,                        &
               WDES1%NB_TOT,                        &
               WDES1%NPRO_O_RED,                    &    ! is zero WDES1%LOVERL= .FALSE.
               Icmplx*hbar_over_m,                  &
               WA%CPROJ_RED(1,1), WDES1%NPROD_RED,  &
               PWA%CPROJ_RED(1,1), WDES1%NPROD_RED, &
               dcmplx(1.0d0),p_k(1,1,DIR),WDES1%NB_TOT)
       ENDIF

       ! WA and PWA to their original distribution
       CALL REDISTRIBUTE_PW(WA)
       CALL REDISTRIBUTE_PROJ(WA)

       CALL REDISTRIBUTE_PW(PWA)
       CALL REDISTRIBUTE_PROJ(PWA)

    ENDDO

    CALLMPI( M_sum_z(WDES1%COMM,p_k, SIZE(p_k)))
    ! sum over all procs working on that k-point

    ! delete P|WA> wavefunction
    CALL DELWAVA(PWA)

    IF (LGAMMA) THEN
       DEALLOCATE(P_G)
    ENDIF
  END SUBROUTINE get_momentum


  SUBROUTINE sum_over_LM(WDES1,CCMIJ,WA1,WA2)
    ! compute <p_I|WA2> = sum_J MIJ <p_J|WA1>
    USE wave
    USE wave_high
    TYPE(wavedes1) :: WDES1
    COMPLEX(8), INTENT(IN) :: CCMIJ(:,:,:,:)
    TYPE(wavefuna), INTENT(IN)  :: WA1
    TYPE(wavefuna), INTENT(IN)  :: WA2

    INTEGER NB,ISPINOR,ISPINORP
    INTEGER NIS,NIG,NT,NI,LM,LMP
    INTEGER NPRO,NPROP,COUNTER
    INTEGER SIZEOFM
    LOGICAL M_type

    SIZEOFM=sizeof( WA2%CPROJ(1,1))
    IF (SIZEOFM .EQ. 8) THEN
       M_type = .FALSE.
       WA2%CPROJ=0.0d0
    ELSEIF (SIZEOFM .EQ. 16) THEN
       M_type = .TRUE.
       WA2%CPROJ=dcmplx(0.0d0)
    ELSE
       CALL vtutor%error('Type error in sum_over_LM')
    ENDIF


    DO NB=1,WDES1%NBANDS
       DO ISPINOR=0,WDES1%NRSPINORS-1
          DO ISPINORP=0,WDES1%NRSPINORS-1

             NPRO=ISPINOR*(WDES1%NPRO/2)
             NPROP=ISPINORP*(WDES1%NPRO/2)

             COUNTER=0
             NIS =1
             DO NT=1,WDES1%NTYP
                DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                   NIG = NI_GLOBAL(NI,WDES1%COMM_INB)
                   DO LM=1,WDES1%LMMAX(NT)
                      DO LMP=1,WDES1%LMMAX(NT)

                         IF (M_type) THEN
                            WA2%CPROJ(LM+NPRO,NB) = WA2%CPROJ(LM+NPRO,NB) + &
                                 CCMIJ(LMP,LM,NIG,1+ISPINORP+2*ISPINOR)* &
                                 WA1%CPROJ(LMP+NPROP,NB)
                         ELSE
                            WA2%CPROJ(LM+NPRO,NB) = WA2%CPROJ(LM+NPRO,NB) + &
                                 real(CCMIJ(LMP,LM,NIG,1+ISPINORP+2*ISPINOR))* &
                                 WA1%CPROJ(LMP+NPROP,NB)

                         ENDIF

                      ENDDO
                   ENDDO
                   NPRO=NPRO+WDES1%LMMAX(NT)
                   NPROP=NPROP+WDES1%LMMAX(NT)
                   COUNTER=COUNTER+WDES1%LMMAX(NT)
                ENDDO
                NIS = NIS+WDES1%NITYP(NT)
             ENDDO
             IF (COUNTER .NE. WDES1%NPRO/WDES1%NRSPINORS) THEN
                CALL vtutor%bug('Error in spinor sector',__FILE__,__LINE__)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE sum_over_LM


  SUBROUTINE getwave_RIJ_LOC(WDES,LMDIM,P,NONL_S, CCRIJ)
    ! RIJ is proportional to identity in the spinor channel
    ! LM indexes are switched to follow vaps convention
    USE prec
    USE PSEUDO
    USE POSCAR
    USE radial
    USE wave
    USE nonl
    USE nonl_struct_def

    TYPE(wavedes), intent(in):: WDES
    INTEGER, INTENT(IN) :: LMDIM
    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(potcar), intent(in) :: P(NONL_S%NTYP)
    COMPLEX(8), intent(out):: CCRIJ(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ, 2)

    INTEGER :: ISP,NI,NT,I,II,N_AT,NMAX
    INTEGER :: LM,L,LL,M,LMP,LP,LLP,MP
    INTEGER NSPINOR,NSPINORP
    TYPE (potcar), POINTER :: PP
    real(8), allocatable :: R(:)
    real(8), allocatable :: TMP(:)
    real(8) H
    REAL(8) val

    CCRIJ=0.0

    ISP=0
    DO NSPINOR=0,WDES%NRSPINORS-1
       DO NSPINORP=0,WDES%NRSPINORS-1
          ISP=ISP+1
          IF (NSPINOR .EQ. NSPINORP) THEN


             N_AT=0
             NI = 0
             DO NT = 1, NONL_S%NTYP
                DO II = 1, NONL_S%NITYP(NT)
                   NI = NI + 1

                   IF (NI_LOCAL(NI, WDES%COMM_INB) == 0) CYCLE

                   N_AT=N_AT+1
                   PP=>PP_POINTER(P, NI, NT)
                   NMAX=PP%R%NMAX
                   H=PP%R%H

                   allocate(R(1:NMAX))
                   allocate(TMP(1:NMAX))

                   LM=0
                   DO L = 1, PP%LMAX
                      LL=PP%LPS(L)
                      DO M=-LL,LL
                         LM=LM+1

                         LMP=0
                         DO LP = 1, PP%LMAX
                            LLP=PP%LPS(LP)
                            DO MP=-LLP,LLP
                               LMP=LMP+1

                               ! first radial integral for wave function
                               DO I=1,NMAX
                                  R(I)=dble(PP%R%R(I)) ! radial grid, R(I)=RSTART*EXP((I-1)*H), in Angstroms
                                  IF (PP%R%R(I) .le. PP%RDEP) THEN
                                     TMP(I)=dble((PP%WAE(I,L)*PP%WAE(I,LP))/(PP%R%R(I)))
                                     !!TMP(I)=dble((PP%WAE(I,L)*PP%WAE(I,LP)))     !  using this instead we obtain QIJ when LM=LMP
                                  ELSE
                                     TMP(I)=0.0d0
                                  ENDIF
                               ENDDO
                               val=int_l(R,H,TMP)
                               CCRIJ(LMP,LM,N_AT,ISP,1)=dcmplx(val)

                               ! second radial integral for wave function
                               DO I=1,NMAX
                                  TMP(I)=dble(PP%WAE(I,LP))
                               ENDDO
                               call get_grad_l(R,H,TMP)

                               DO I=1,NMAX
                                  IF (PP%R%R(I) .le. PP%RDEP) THEN
                                     TMP(I)=dble(PP%WAE(I,L))*(TMP(I)-dble(PP%WAE(I,LP)/R(I)))
                                  ELSE
                                     TMP(I)=0.0d0
                                  ENDIF
                               ENDDO
                               val=int_l(R,H,TMP)
                               CCRIJ(LMP,LM,N_AT,ISP,2)=dcmplx(val)


                               !                         ! first radial integral for pseudo wave function
                               DO I=1,NMAX
                                  R(I)=dble(PP%R%R(I))
                                  IF (PP%R%R(I) .le. PP%RDEP) THEN
                                     TMP(I)=dble((PP%WPS(I,L)*PP%WPS(I,LP))/(PP%R%R(I)))
                                     !TMP(I)=dble((PP%WPS(I,L)*PP%WPS(I,LP)))     !  using this instead we obtain QIJ when LM=LMP
                                  ELSE
                                     TMP(I)=0.0d0
                                  ENDIF
                               ENDDO
                               val=int_l(R,H,TMP)
                               CCRIJ(LMP,LM,N_AT,ISP,1)=CCRIJ(LMP,LM,N_AT,ISP,1)-dcmplx(val)

                               ! second radial integral for pseudo wave function
                               DO I=1,NMAX
                                  TMP(I)=dble(PP%WPS(I,LP))
                               ENDDO
                               call get_grad_l(R,h,TMP)

                               DO I=1,NMAX
                                  IF (PP%R%R(I) .le. PP%RDEP) THEN
                                     TMP(I)=dble(PP%WPS(I,L))*(TMP(I)-dble(PP%WPS(I,LP)/R(I)))
                                  ELSE
                                     TMP(I)=0.0d0
                                  ENDIF
                               ENDDO
                               val=int_l(R,H,TMP)
                               CCRIJ(LMP,LM,N_AT,ISP,2)=CCRIJ(LMP,LM,N_AT,ISP,2)-dcmplx(val)


                            ENDDO
                         ENDDO

                      ENDDO
                   ENDDO

                   deallocate(TMP)
                   deallocate(R)

                ENDDO
             ENDDO

             IF (WDES%NIONS .NE. N_AT) THEN
                CALL vtutor%error('wrong local number of ions')
             ENDIF

          ENDIF
       ENDDO
    END DO

  END SUBROUTINE getwave_RIJ_LOC


  SUBROUTINE get_CCRIJ_A_LOC(WDES,LMDIM,NONL_S,CCRIJ,P,CCRIJ_A)

    USE prec
    USE PSEUDO
    USE POSCAR
    USE radial
    USE wave
    USE nonl_struct_def
    USE nonl
    TYPE(wavedes), intent(in):: WDES
    INTEGER, INTENT(IN) :: LMDIM
    TYPE(nonl_struct), intent(in) :: NONL_S
    COMPLEX(8), intent(in):: CCRIJ(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ, 2)
    TYPE(potcar), intent(in) :: P(NONL_S%NTYP)
    COMPLEX(8), intent(out):: CCRIJ_A(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ, 3)     ! RIJ coefficients with angular part
    REAL(8) Nlm,A(1:3,1:6),B(1:3,1:6),C(1:3,1:6),T0(1:3,1:3),T1(1:3,1:3)
    REAL(8) v0(3,1),v1(3,1),v2(3,1),v3(3,1)
    REAL(8) Norm,NormP,Plm(6,1)

    INTEGER:: NIS, NT, NI, N_AT, L, LP, LDIM, NMAX, I, II, M, MP, LM, LMP, ISP, LL, LLP, DIR
    TYPE (potcar), POINTER :: PP
    integer i1,i2,i3

    DO ISP = 1, WDES%NCDIJ

       N_AT=0
       NI = 0
       DO NT = 1, NONL_S%NTYP
          DO II = 1, NONL_S%NITYP(NT)
             NI = NI + 1

             IF (NI_LOCAL(NI, WDES%COMM_INB) == 0) CYCLE

             N_AT=N_AT+1
             PP=>PP_POINTER(P, NI, NT)

             LM=0
             DO L = 1, PP%LMAX
                LL=PP%LPS(L)
                DO M=-LL,LL
                   LM=LM+1
                   call get_NYlm(LL,abs(M),Norm)

                   LMP=0
                   DO LP = 1, PP%LMAX
                      LLP=PP%LPS(LP)
                      DO MP=-LLP,LLP
                         LMP=LMP+1

                         call get_NYlm(LLP,abs(MP),NormP)
                         call get_ABC(LLP,abs(MP),A,B,C)

                         call get_T01(M,MP,T0,T1)


                         call get_vec_P(LL,abs(M),LLP,abs(MP),Plm)

                         v0=matmul(A,Plm)
                         v1=matmul(T0,v0)

                         v0=matmul(B,Plm)
                         v2=matmul(T0,v0)

                         v0=matmul(C,Plm)
                         v3=matmul(T1,v0)


                         DO DIR=1,3
                            CCRIJ_A(LMP,LM,N_AT,ISP,DIR)=CCRIJ(LMP,LM,N_AT,ISP,2)*v1(DIR,1)*((-Norm*NormP)/((2.0d0*LLP+1.0d0))) &
                                 +CCRIJ(LMP,LM,N_AT,ISP,1)*v2(DIR,1)*((Norm*NormP)/((4.0d0*LLP+2.0d0))) &
                                 +CCRIJ(LMP,LM,N_AT,ISP,1)*v3(DIR,1)*((Norm*NormP)/(2.0d0))
                         ENDDO

                      ENDDO
                   ENDDO

                ENDDO
             ENDDO

          ENDDO
       ENDDO

       IF (WDES%NIONS .NE. N_AT) THEN
          CALL vtutor%error('wrong local number of ions')
       ENDIF

    ENDDO

  END SUBROUTINE get_CCRIJ_A_LOC


  SUBROUTINE get_momentum_matrix_test(WDES,NONL_S,LATT_CUR,T_INFO,P,W)
    USE LATTICE
    USE mgrid
    USE POSCAR
    USE PSEUDO
    USE nonl_struct_def
    USE nonl
    USE wave
    USE wave_high

    TYPE(wavedes), intent(in) :: WDES
    TYPE(nonl_struct), intent(in) :: NONL_S
    TYPE(latt), intent(in) :: LATT_CUR
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(potcar), intent(in) :: P(T_INFO%NTYP)
    TYPE(wavespin), intent(in) :: W
    TYPE(wavedes1) :: WDES1
    TYPE(wavefuna) WA
    INTEGER NK,ISP,LMDIM,DIR,I,J,K
    COMPLEX(8), ALLOCATABLE :: CCRIJ(:,:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: CCRIJ_A(:,:,:,:,:)
    COMPLEX(8), ALLOCATABLE :: p_k(:,:,:,:)
    REAL(8) pi
    character*80 string
    pi=acos(-1.0d0)

    ! This is what we want to compute
    ALLOCATE(p_k(WDES%NB_TOT,WDES%NB_TOT,3,WDES%ISPIN))

    LMDIM = P(1)%LMDIM

    ALLOCATE( CCRIJ( LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ, 2))
    CCRIJ = dcmplx(0)
    CALL  getwave_RIJ_LOC(WDES,LMDIM,P,NONL_S, CCRIJ)              ! get radial part of R matrix

    ALLOCATE(CCRIJ_A(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ, 3))
    CCRIJ_A=dcmplx(0)
    call get_CCRIJ_A_LOC(WDES,LMDIM,NONL_S,CCRIJ,P,CCRIJ_A)        ! multiply R matrix with angular part
    ! and switch LM indexes to
    ! follow vasp convention
    deallocate(CCRIJ)


    DO NK = 1, WDES%NKPTS
       IF (MOD(NK - 1, WDES%COMM_KINTER%NCPU) .NE. WDES%COMM_KINTER%NODE_ME - 1) CYCLE
       ! if the local node does not belong to the NK-point, we skipp the computation.

       CALL SETWDES(WDES, WDES1, NK)    ! set descriptor for actual  WDES1%NK = NK

       DO ISP = 1, WDES%ISPIN

          WA = ELEMENTS(W, WDES1, ISP)  ! 1-k, 1-spin, but all bands

          call get_momentum(LATT_CUR,WA,LMDIM,CCRIJ_A,p_k(:,:,:,ISP))

       ENDDO

       IF (WDES1%COMM_INTER%NODE_ME==WDES1%COMM_INTER%IONODE) THEN
          write (string,'(I0)') NK

          ! output velocities intext form
          OPEN(UNIT=71, STATUS='REPLACE', FILE='VEL_'//trim(string))
          DO ISP = 1, WDES%ISPIN
             DO I=1,WDES%NB_TOT
                DO J=1,WDES%NB_TOT
                   write(71,'(6F20.10)') p_k(I,J,1:3,ISP)
                ENDDO
             ENDDO
          ENDDO
          CLOSE(71)

          ! output velocities in binary format
          OPEN(UNIT=71, STATUS='REPLACE', FILE='VEL_'//trim(string)//'.bin',&
               ACCESS='STREAM', FORM='UNFORMATTED')
          WRITE(71) p_k
          CLOSE(71)

       ENDIF

    ENDDO

    IF (WDES%COMM%NODE_ME == WDES%COMM%IONODE) THEN
       OPEN(UNIT=71, STATUS='REPLACE', FILE='EIG_VEL')
       WRITE(71,'(A1,1X,3I9)')'%', WDES%NB_TOT,WDES%ISPIN,WDES%NKPTS
       DO NK = 1, WDES%NKPTS
          DO ISP = 1, WDES%ISPIN
             DO I=1, WDES%NB_TOT
                write(71,'(4F20.10)') WDES%VKPT(1:3,NK),real(W%CELTOT(I,NK,ISP))
             ENDDO
          ENDDO
       ENDDO
       CLOSE(71)
    ENDIF

    DEALLOCATE(CCRIJ_A)
  END SUBROUTINE get_momentum_matrix_test


  SUBROUTINE grid_in_supercell(NG,M,Nsp,g_pr,prim_ind)
    integer, intent(in) :: NG(3,3) ! discretization matrix of the current(primitive) cell
    integer, intent(in) :: M(3,3)  ! primitive to conventionnal cell passage matrix
    integer, intent(in) :: Nsp(3,3) ! conventionnal to supercell passage matrix
    real(8), allocatable, intent(out) :: g_pr(:,:)
    integer(8), allocatable, intent(out) :: prim_ind(:,:)

    integer, allocatable :: A(:,:)
    integer, allocatable :: D(:,:)
    integer, allocatable :: S(:,:)
    integer, allocatable :: T(:,:)
    integer grid_size
    integer i1,i2,i3,ind
    real(8) v1(3,1),v2(3,1)

    ! passage matrix from microcell to supercell
    allocate(A(1:3,1:3))
    A=matmul(NG,matmul(M,Nsp))

    ! Smith normal form of A: D=S*A*T
    call SNF(A,D,S,T)

    ! grid size
    grid_size=D(1,1)*D(2,2)*D(3,3)
    allocate(g_pr(1:3,1:grid_size)) ! contains grid points in reduced coordinates of the primitive cell
    allocate(prim_ind(1:3,1:grid_size)) ! contains the indexes of wavefunction arrays (redundant)

    ! here we should put some checks to be sure of SNF, and that the grid size is what we want

    ! the grid
    !write(*,*) 'D=',D(1,1),D(2,2),D(3,3)
    !write(*,*) 'T=',T
    !write(*,*) 'M=',M
    !write(*,*) 'Nsp',Nsp
    ind=0
    do i3=0,D(3,3)-1
       do i2=0,D(2,2)-1
          do i1=0,D(1,1)-1
             ind=ind+1

             v1(1,1)=(1.0d0*i1)/(1.0d0*D(1,1))
             v1(2,1)=(1.0d0*i2)/(1.0d0*D(2,2))
             v1(3,1)=(1.0d0*i3)/(1.0d0*D(3,3))

             v2(1:3,1)=matmul(1.0d0*T(1:3,1:3),v1(1:3,1))

             v2(1,1)=modulo(v2(1,1),1.0d0)
             v2(2,1)=modulo(v2(2,1),1.0d0)
             v2(3,1)=modulo(v2(3,1),1.0d0)

             g_pr(1:3,ind)=matmul(1.0d0*M(1:3,1:3),matmul(1.0d0*Nsp(1:3,1:3),v2(1:3,1)))
             !write(*,*) i1,i2,i3,ind,g_pr(1:3,ind)
             prim_ind(1,ind)=modulo(nint(g_pr(1,ind)*NG(1,1)),NG(1,1))+1
             prim_ind(2,ind)=modulo(nint(g_pr(2,ind)*NG(2,2)),NG(2,2))+1
             prim_ind(3,ind)=modulo(nint(g_pr(3,ind)*NG(3,3)),NG(3,3))+1


          enddo
       enddo
    enddo
    deallocate(A)
    deallocate(D)
    deallocate(S)
    deallocate(T)
  END SUBROUTINE grid_in_supercell


  SUBROUTINE SNF(B,DD,SS,TT)
    ! Smith normal form DD=SS*A*TT
    integer, allocatable, intent(in) :: B(:,:)
    integer, allocatable, intent(out) :: DD(:,:)
    integer, allocatable, intent(out) :: SS(:,:)
    integer, allocatable, intent(out) :: TT(:,:)

    integer I,J,K,L,M,N,U,V,MPLUSN
    real(8), allocatable :: A(:,:)
    real(8) T,S,X,Y

    !
    !   Transfer data to temporaries
    !
    M=size(B(:,1))
    N=size(B(1,:))
    allocate(A(1:M+N,1:M+N))

    MPLUSN=M+N

    A=0.0d0
    A(1:M,1:N)=B(1:M,1:N)

    !
    !   Set P = I
    !
    do I=1,M
       A(I,N+I)=1.0d0
    enddo

    !
    !   Set Q = I
    !
    do J=1,N
       A(M+J,J)=1.0d0
    enddo

    !
    !   Set l = min(m,n)
    !
    L=(M+N-ABS(M-N))/2

    !
    !   Initialize k
    !
    K=1

    !
    !   Call subroutine NZM
    !
5   CONTINUE
    call NZM(N,M,A,K,I,J,T)
    IF(T.EQ. 0.0d0) then
       GOTO 70
    ENDIF
    !
    !   Interchange rows i and k,if i and k are different
    !
    if (I .NE. K) THEN
       do V=K,MPLUSN
          S=A(K,V)
          A(K,V)=A(I,V)
          A(I,V)=S
       enddo
    endif

    !
    !   Interchange columns j and k, if j and k are different
    !
    if (J .NE. K) THEN
       do U=K,MPLUSN
          S=A(U,K)
          A(U,K)=A(U,J)
          A(U,J)=S
       enddo
    endif

    !
    !   Replace row k by its negative, if A(k,k) < 0
    !
    if ( A(K,K) .LT. 0.0d0) THEN
       do V=K,MPLUSN
          A(K,V)=-A(K,V)
       enddo
    endif

    !
    !   Replace row i by row i - x row k,
    !   x = [A(i,k)/A(k,k)], i = k+1, ... ,m
    !
    do I=K,M
       if (I .EQ. K) then
          GOTO 35
       endif
       X=DINT(A(I,K)/T)
       do J=K,MPLUSN
          A(I,J)=A(I,J)-X*A(K,J)
       enddo
35     CONTINUE
    enddo

    !
    !   Replace col j by col j - y col k,
    !   where y =[A(k,j)/A(k,k)], j = k+1, ... ,n
    !

40  CONTINUE
    do J=K,N
       IF(J .EQ. K) THEN
          GOTO 45
       ENDIF
       Y=DINT(A(K,J)/T)
       do I=K,MPLUSN
          A(I,J)=A(I,J)-Y*A(I,K)
       enddo
45     CONTINUE
    enddo

    !
    !   Check to see whether or not all A(i,k) = 0,
    !   i = k+1, ... ,m
    !
    do I=K,M
       IF (I .EQ. K) THEN
          GOTO 50
       ENDIF
       IF( A(I,K) .NE. 0.0d0) THEN
          GOTO 5
       ENDIF
50     CONTINUE
    enddo

    !
    !   Check to see whether or not all A(k,j) = 0,
    !   j = k+1, ... ,n
    !
    do J=K,N
       IF( J .EQ. K) THEN
          GOTO 55
       ENDIF
       IF( A(K,J) .NE. 0.0d0) THEN
          GOTO 5
       ENDIF
55     CONTINUE
    enddo

    !
    !   Call subroutine TFD
    !
    CALL TFD(N,M,A,K,U)
    IF(U .EQ. 0) THEN
       GOTO 65
    ENDIF
    !
    !   Replace row k by row k + row p
    !
    do J=K,MPLUSN
       A(K,J)=A(K,J)+A(U,J)
    enddo
    GOTO 40

    !
    !   Check for completion
    !
65  CONTINUE
    IF (K .NE. L)THEN
       K=K+1
       GOTO 5
    ENDIF

    !
    !   Store invariant factors in DD()
    !
70  CONTINUE
    allocate(DD(1:M,1:N))
    DD=0
    do I=1,L
       DD(I,I)=NINT(A(I,I))
    enddo

    !
    !   Store left multiplier in SS()
    !
    allocate(SS(1:M,1:M))
    do I=1,M
       do J=1,M
          SS(I,J)=NINT(A(I,N+J))
       enddo
    enddo

    !
    !   Store right multiplier in TT()
    !
    allocate(TT(1:N,1:N))
    do I=1,N
       do J=1,N
          TT(I,J)=NINT(A(M+I,J))
       enddo
    enddo

    deallocate(A)
    RETURN
  END SUBROUTINE SNF

  !
  !   t is set equal to 0 if Ak = 0; otherwise
  !   t is set equal to abs(A(i,j)), where A(i,j)
  !   is a nonzero element which is least
  !   in absolute value
  !
  SUBROUTINE NZM(N,M,A,K,I,J,T)
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(8), intent(in) :: A(1:M+N,1:N+M)
    integer, intent(in) :: K
    integer, intent(out) :: I
    integer, intent(out) :: J
    real(8), intent(out) :: T

    integer U,V
    real S

    do I=K,M
       do J=K,N
          T=ABS(A(I,J))
          if (T .NE. 0.0d0) THEN
             GOTO 20
          endif
       enddo
    enddo
20  CONTINUE
    IF(T.EQ. 0.0d0) THEN
       RETURN
    ENDIF

    do U=K,M
       do V=K,N
          S=ABS(A(U,V))
          IF(S .EQ. 0.0d0) THEN
             GOTO 30
          ENDIF
          IF(S .GE. T) THEN
             GOTO 30
          ENDIF
          T=S
          I=U
          J=V
30        CONTINUE
       enddo
    enddo
    RETURN
  END SUBROUTINE NZM
  !
  !   t is set equal to 0 if Ak is congruent to 0
  !   modulo A(k,k); otherwise, t is set equal to p,
  !   where A(p,q) is an element of Ak which is
  !   not divisible by A(k,k)
  !
  SUBROUTINE TFD(N,M,A,K,T)
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(8), intent(in) :: A(1:M+N,1:M+N)
    integer, intent(in) :: K
    integer, intent(out) :: T
    integer I,J

    T=0
    do I=K,M
       do J=K,N
          IF(DMOD(A(I,J),A(K,K)) .NE. 0) THEN
             GOTO 20
          ENDIF
       enddo
    enddo
    RETURN
20  T=I
    RETURN
  END SUBROUTINE TFD

  SUBROUTINE read_dV(g_pr,NCDIJ,N_atom_conv,dV)
    !  N_atom_conv : number of atoms in conventionnal cell
    real(8), allocatable, intent(in) :: g_pr(:,:)
    integer, intent(in) :: NCDIJ
    integer, intent(in) :: N_atom_conv
    real(8), allocatable :: g_pr_read(:,:)
    complex(8), allocatable, intent(out) :: dV(:,:,:,:)
    integer(8) B_size, i, j

    INQUIRE(FILE='dV.bin', SIZE=B_size)     ! in bytes

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! only the master node should do this :
    if (B_size .NE. 16*size(g_pr(1,:))*NCDIJ*N_atom_conv*3) then
       write(*,*)'B_size',B_size,16*size(g_pr(1,:))*NCDIJ*N_atom_conv*3,16,size(g_pr(1,:)),NCDIJ,N_atom_conv,3
       write(*,*) 'The 3 potential derivatives should be given as a ',NCDIJ, &
            'entries matrix, for ',N_atom_conv,' atoms, on the grid in  GRID.bin'
       OPEN(UNIT=10, STATUS='REPLACE', FILE='GRID.bin',ACCESS='STREAM',FORM='UNFORMATTED')
       WRITE(10) g_pr
       CLOSE(10)
       call vtutor%error('error reading dV')

    else
       allocate(dV(1:size(g_pr(1,:)),1:NCDIJ,1:N_atom_conv,1:3))

       open(unit = 14,file='dV.bin',form='UNFORMATTED',ACCESS='STREAM')
       read(14) dV
       close(14)

       allocate(g_pr_read(1:3,1:size(g_pr(1,:))))

       open(unit = 14,file='GRID.bin', form='UNFORMATTED',ACCESS='STREAM')
       read(14) g_pr_read
       close(14)
       do i=1,3
          do j=1,size(g_pr(1,:))
             if (abs(g_pr(i,j)-g_pr_read(i,j)) .gt. 0.00001) then
                write(*,*) i,j,g_pr(i,j),g_pr_read(i,j)
                call vtutor%error('Error in read_dV')
             endif
          enddo
       enddo
       write(*,*) 'dV has been read'
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE read_dV



  SUBROUTINE read_and_set_dV(T_INFO,LATT_CUR,NGX,NGY,NGZ,NCDIJ,dV,g_pr,prim_ind,grids,centers_prim,R_shift,cut_off_r)
    USE POSCAR
    USE lattice
    TYPE(type_info), intent(in) :: T_INFO
    TYPE(latt), intent(in) :: LATT_CUR
    integer, intent(in) :: NGX ! discretization matrix of the current(primitive) cell
    integer, intent(in) :: NGY ! discretization matrix of the current(primitive) cell
    integer, intent(in) :: NGZ ! discretization matrix of the current(primitive) cell
    integer, intent(in) :: NCDIJ

    complex(8), allocatable, intent(out) :: dV(:,:,:,:)
    real(8), allocatable, intent(out) :: g_pr(:,:)
    integer(8), allocatable, intent(out) :: prim_ind(:,:)
    real(8), allocatable, intent(out) :: grids(:,:,:,:)
    real(8), allocatable, intent(out) :: centers_prim(:,:)
    real(8), allocatable, intent(out) :: R_shift(:,:)
    real(8), intent(out) :: cut_off_r
    ! read wap_INCAR
    integer :: M(3,3)  ! primitive to conventionnal cell passage matrix
    integer :: Nsp(3,3) ! conventionnal to supercell passage matrix
    integer :: NdV
    character*4 :: mode_dV      ! 'prim','conv','supc'
    character*4 :: mode_gkkp      ! 'prim','conv','supc'
    real(8) :: epsilon_r

    ! read dV.bin
    complex(8), allocatable :: dV_full(:,:,:,:)
    real(8), allocatable :: centers(:,:)
    real(8), allocatable :: grids_full(:,:,:,:)
    integer, allocatable :: prim_mapping(:)
    integer NG(3,3)
    integer i_prim,i

    ! grid point in primitive/current cell
    NG= 0; NG(1,1)= NGX ; NG(2,2)= NGY; NG(3,3)= NGZ

    ! read parameters which define supercell and potential derivatives
    ! cut_off_r added
    call read_INCAR_wap(M,Nsp,NdV,mode_dV,epsilon_r,cut_off_r,mode_gkkp,centers)
    write(*,*) 'M=',M
    write(*,*) 'Nsp=',Nsp
    write(*,*) 'NdV=',NdV
    write(*,*) 'mode_dV= ',mode_dV
    write(*,*) 'epsilon_r= ',epsilon_r
    write(*,*) 'cut_off_r= ',cut_off_r
    write(*,*) 'mode_gkkp= ',mode_gkkp

    ! Generate a grid in the supercell which match the primitive/current grid point mesh.
    ! The supercell has same origin than the primitive/current cell.
    ! g_pr contains grid points in reduced coordinates of the primitive cell.
    ! prim_ind are the modular integral coordinates of the grid point folded back to prim cell
    call grid_in_supercell(NG,M,Nsp,g_pr,prim_ind)


    ! read local potential and grid where it is know.
    ! If the grid read in GRID.bin is different from g_pr, then we have an error.
    call read_dV(g_pr,NCDIJ,NdV,dV_full)


    !if not read in INCAR_wap generate the positions of the atoms for which the potential derivatives in dV.bin are given can be
    !generated like this:
    ! call define_atom_centers(mode_dV,M,Nsp,T_INFO,centers)
    ! but the order has to be check
    ! -> using hdf5 format and include this info in the potential file would be nice, but i dont know how to do it.

    do i=1,size(centers(1,:))
       write(*,*) 'c=',centers(1:3,i)
    enddo


    if (NdV .ne. size(centers(1,:))) then
       call vtutor%error('Error in the number of atoms')
    endif

    ! symmetrized grid around atoms positions in "centers"
    ! grids is initialized to cut_off_r
    call sym_grid_around_atoms(epsilon_r,cut_off_r,M,Nsp,LATT_CUR%A(1:3,1:3),g_pr,centers,grids_full)

    ! generate atomic position where we want <kn|dV/dR|k'n'>
    call define_atom_centers(mode_gkkp,M,Nsp,T_INFO,centers_prim)

    ! compute the mapping between centers and centers_prim
    call map_to_prim(epsilon_r,centers,centers_prim,prim_mapping,R_shift)


    ! get grid and potentias for the atoms in centers_prim
    allocate(dV(1:size(g_pr(1,:)),1:NCDIJ,1:size(centers_prim(1,:)),1:3))
    allocate(grids(1:size(centers_prim(1,:)),1:size(grids_full(1,:,1,1)),1:12,1:3))

    do i_prim=1,size(centers_prim(1,:))
       dV(1:size(g_pr(1,:)),1:NCDIJ,i_prim,1:3)=dV_full(1:size(g_pr(1,:)),1:NCDIJ,prim_mapping(i_prim),1:3)
       grids(i_prim,1:size(grids_full(1,:,1,1)),1:12,1:3)=grids_full(prim_mapping(i_prim),1:size(grids_full(1,:,1,1)),1:12,1:3)
    enddo

    deallocate(dV_full)
    deallocate(centers)
    deallocate(grids_full)
    deallocate(prim_mapping)

  END SUBROUTINE read_and_set_dV


  SUBROUTINE map_to_prim(epsilon_r,centers,centers_prim,prim_mapping,R_shift)
    ! for each atom in centers_prim we find an equivalent one in centers
    real(8), intent(in) :: epsilon_r
    real(8), allocatable, intent(in) :: centers(:,:)
    real(8), allocatable, intent(in) :: centers_prim(:,:)
    real(8), allocatable, intent(out) :: R_shift(:,:)
    integer, allocatable, intent(out) :: prim_mapping(:)
    integer i_at,i_prim
    real(8) epsilon_l

    allocate(prim_mapping(size(centers_prim(1,:))))
    prim_mapping(:)=-1
    allocate(R_shift(1:3,size(centers_prim(1,:))))

    do i_prim=1,size(centers_prim(1,:))
       do i_at=1,size(centers(1,:))
          epsilon_l=abs(modulo(centers_prim(1,i_prim)-centers(1,i_at),1.0d0)) + &
               abs(modulo(centers_prim(2,i_prim)-centers(2,i_at),1.0d0)) + &
               abs(modulo(centers_prim(3,i_prim)-centers(3,i_at),1.0d0))
          if (epsilon_l .lt. epsilon_r) then
             prim_mapping(i_prim)=i_at
             R_shift(1:3,i_prim)=centers_prim(1:3,i_prim)-centers(1:3,i_at)
          endif
       enddo
       if (prim_mapping(i_prim) .le. 0) then
          call vtutor%error('Error in map_to_prim')
       endif
    enddo
  END SUBROUTINE map_to_prim


  SUBROUTINE sym_grid_around_atoms(tol,cut_off,M,Nsp,lattice,g_pr,centers,grids)
    real(8), intent(in) :: tol
    real(8), intent(in) :: cut_off
    integer, intent(in) :: M(3,3)  ! primitive to conventionnal cell passage matrix
    integer, intent(in) :: Nsp(3,3) ! conventionnal to supercell passage matrix
    real(8), intent(in) :: lattice(3,3) ! cartesian components of (current/primitive) lattice vectors, as columns
    real(8), allocatable, intent(in) :: g_pr(:,:)
    real(8), allocatable, intent(in) :: centers(:,:)
    real(8), allocatable, intent(out) :: grids(:,:,:,:) ! grids(atoms, grid size, add. grid point, 1:3)

    integer M0(3,3)
    real(8) lattice_sp(3,3),lattice_sp_perp(3,3)
    real(8) X(3,1),rho
    integer atoms,i,j,ind,indM
    integer, allocatable :: Lpoints(:,:)

    M0=matmul(M,Nsp)
    lattice_sp=matmul(lattice,1.0d0*M0)
    call  decomp_A_per(lattice_sp,lattice_sp_perp)

    allocate(grids(1:size(centers(1,:)),1:size(g_pr(1,:)),1:12,1:3))
    ! a clever initial guess: better idea ?
    grids=cut_off

    do atoms=1,size(centers(1,:))

       do i=1,size(g_pr(1,:))
          X(1:3,1)=matmul(lattice(1:3,1:3),centers(1:3,atoms)-g_pr(1:3,i))
          rho=sqrt(dot_product(X(1:3,1),X(1:3,1)))
          call NN_lattice_points(lattice_sp,lattice_sp_perp,X,rho,tol,Lpoints)

          if (size(Lpoints(1,:)) .gt. 12 ) then
             call vtutor%error('Error1 in sym_grid_around_atoms')
          elseif (size(Lpoints(1,:)) .le. 0 ) then
             call vtutor%error('Error2 in sym_grid_around_atoms')
          else

             do j=1,size(Lpoints(1,:))
                X(1:3,1)=matmul(1.0d0*M0,1.0d0*Lpoints(1:3,j))
                grids(atoms,i,j,1)=g_pr(1,i)+X(1,1)
                grids(atoms,i,j,2)=g_pr(2,i)+X(2,1)
                grids(atoms,i,j,3)=g_pr(3,i)+X(3,1)
             enddo

          endif

          deallocate(Lpoints)
       enddo

    enddo
  END SUBROUTINE sym_grid_around_atoms


  SUBROUTINE decomp_par_perp(A,X,X_par,X_per)
    ! decompose the vector X into a part that is parallel to the range of A, and a part that is perpendicular to it: X = X_par + X_per
    ! the vectors and matrices must be allocated before the call.
    real(8), intent(in) :: A(3,2)
    real(8), intent(in) :: X(3,1)
    real(8),  intent(out) :: X_par(3,1)
    real(8), intent(inout) :: X_per(3,1)
    integer n,i,j
    real(8)  g(2,2)
    real(8)  y(2,1)
    real(8)  z(2,1)
    real(8), allocatable :: work(:)
    integer, allocatable ::  ipiv(:)
    real(8) Xn
    integer INFO

    ! sub-metric tensor
    do i=1,2
       do j=1,2
          g(i,j)=dot_product(A(:,i),A(:,j))
       enddo
    enddo

    ! compute inverse
    allocate(ipiv(1:2))
    allocate(work(1:2))
    call DGETRF( 2, 2,g, 2, ipiv, INFO )
    call DGETRI( 2, g, 2, ipiv, work, 2, INFO )
    deallocate(work)
    deallocate(ipiv)

    ! projection  of X on the range of A
    do i=1,2
       y(i,1)=dot_product(X(:,1),A(:,i))
    enddo


    ! components  of X in the range of A
    z=matmul(g,y)

    ! X_par

    X_par=0.0d0
    do i=1,2
       do j=1,3
          X_par(j,1)=X_par(j,1)+z(i,1)*A(j,i)
       enddo
    enddo

    ! X_per
    X_per=X-X_par

  END SUBROUTINE decomp_par_perp


  SUBROUTINE decomp_A_per(A,A_per)
    ! A_per contain as column the vectors A1_per,A2_per,A3_per, which are  components of A_i perpendicular to the other 2 vectors
    real(8), intent(in) :: A(1:3,1:3)
    real(8), intent(out) :: A_per(1:3,1:3)
    real(8)  A1_par(3,1)
    real(8)  A1_per(3,1)
    real(8)  A2_par(3,1)
    real(8)  A2_per(3,1)
    real(8)  A3_par(3,1)
    real(8)  A3_per(3,1)
    real(8)  A12(3,2)
    real(8)  A13(3,2)
    real(8)  A23(3,2)
    real(8) vec(3,1)

    A12(:,1)=A(:,1)
    A12(:,2)=A(:,2)

    A13(:,1)=A(:,1)
    A13(:,2)=A(:,3)

    A23(:,1)=A(:,2)
    A23(:,2)=A(:,3)

    vec(:,1)=A(:,1)
    call decomp_par_perp(A23,vec,A1_par,A1_per)

    vec(:,1)=A(:,2)
    call decomp_par_perp(A13,vec,A2_par,A2_per)

    vec(:,1)=A(:,3)
    call decomp_par_perp(A12,vec,A3_par,A3_per)

    A_per(:,1)=A1_per(:,1)
    A_per(:,2)=A2_per(:,1)
    A_per(:,3)=A3_per(:,1)

  END SUBROUTINE decomp_A_per


  SUBROUTINE NN_lattice_points(A,A_per,X,rho,tol,Lpoints)
    ! return the nearest lattice points to vector X. If X is a lattice point, it is alqo included in the list.
    real(8), intent(in) :: A(1:3,1:3)
    real(8), intent(in) :: A_per(1:3,1:3)
    real(8), intent(in) :: X(1:3,1:1)
    real(8), intent(in) :: rho
    real(8), intent(in) :: tol
    integer, allocatable, intent(out) :: Lpoints(:,:)
    integer i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
    integer i1,i2,i3,ind
    real(8) P(3,1),norm,norm0
    integer, allocatable :: temp(:,:)

    i1_min=floor(dot_product(X(1:3,1),A_per(1:3,1))/dot_product(A_per(1:3,1),A_per(1:3,1))-rho/sqrt(dot_product(A_per(1:3,1),A_per(1:3,1))))
    i1_max=ceiling(dot_product(X(1:3,1),A_per(1:3,1))/dot_product(A_per(1:3,1),A_per(1:3,1))+rho/sqrt(dot_product(A_per(1:3,1),A_per(1:3,1))))

    i2_min=floor(dot_product(X(1:3,1),A_per(1:3,2))/dot_product(A_per(1:3,2),A_per(1:3,2))-rho/sqrt(dot_product(A_per(1:3,2),A_per(1:3,2))))
    i2_max=ceiling(dot_product(X(1:3,1),A_per(1:3,2))/dot_product(A_per(1:3,2),A_per(1:3,2))+rho/sqrt(dot_product(A_per(1:3,2),A_per(1:3,2))))

    i3_min=floor(dot_product(X(1:3,1),A_per(1:3,3))/dot_product(A_per(1:3,3),A_per(1:3,3))-rho/sqrt(dot_product(A_per(1:3,3),A_per(1:3,3))))
    i3_max=ceiling(dot_product(X(1:3,1),A_per(1:3,3))/dot_product(A_per(1:3,3),A_per(1:3,3))+rho/sqrt(dot_product(A_per(1:3,3),A_per(1:3,3))))

    allocate(temp(1:3,1:(i1_max-i1_min+1)*(i2_max-i2_min+1)*(i3_max-i3_min+1)))

    ind=0
    norm0=rho
    do i1=i1_min,i1_max
       do i2=i2_min,i2_max
          do i3=i3_min,i3_max
             P(:,1)=1.0d0*i1*A(:,1)+1.0d0*i2*A(:,2)+1.0d0*i3*A(:,3)
             P=P-X
             norm=sqrt(dot_product(P(1:3,1),P(1:3,1)))
             if (norm .le. 1.0e-10) then
                ind=ind+1
                temp(1,ind)=i1
                temp(2,ind)=i2
                temp(3,ind)=i3
             elseif (norm .lt. norm0-tol) then
                norm0=norm
                ind=1
                temp(1,ind)=i1
                temp(2,ind)=i2
                temp(3,ind)=i3
             elseif (abs(norm-norm0) .le. tol) then
                ind=ind+1
                temp(1,ind)=i1
                temp(2,ind)=i2
                temp(3,ind)=i3
             endif
          enddo
       enddo
    enddo

    allocate(Lpoints(1:3,1:ind))
    Lpoints(1:3,1:ind)=temp(1:3,1:ind)
    deallocate(temp)
  END SUBROUTINE NN_lattice_points


  SUBROUTINE define_atom_centers(mode,M,Nsp,T_INFO,centers)
    USE POSCAR
    character*4 mode       ! 'prim','conv','supc'
    TYPE(type_info), intent(in) :: T_INFO
    integer, intent(in) :: M(3,3)  ! primitive to conventionnal cell passage matrix
    integer, intent(in) :: Nsp(3,3) ! conventionnal to supercell passage matrix
    real(8), allocatable, intent(out) :: centers(:,:)
    real(8), allocatable :: pos(:,:)
    integer i,Msp(3,3)

    ! not to work with pointers; may be removed if needed
    allocate(pos(1:3,1:T_INFO%NIONS))
    pos(1:3,1:T_INFO%NIONS)=T_INFO%POSION(1:3,1:T_INFO%NIONS)


    if (mode .eq. 'prim') then
       allocate(centers(1:3,1:T_INFO%NIONS))
       do i = 1, T_INFO%NIONS
          centers(1:3,i)=pos(1:3, i)
       enddo
    elseif (mode .eq. 'conv') then
       call find_conv_cell_atoms(pos,M,centers)
    elseif (mode .eq. 'supc') then
       Msp=matmul(M,Nsp)
       call find_conv_cell_atoms(pos,Msp,centers)
    else
       call vtutor%error('Error in define_atom_centers')
    endif

  END SUBROUTINE define_atom_centers


  SUBROUTINE find_conv_cell_atoms(pos,M,pos_conv)
    ! return position in reduced coordinates of the primitive cell for : primitive, conv, or supercell
    real(8), allocatable, intent(in) :: pos(:,:)     ! not 8, probably dp
    integer, intent(in) :: M(3,3)  ! primitive to conventionnal cell passage matrix
    real(8), allocatable, intent(out) :: pos_conv(:,:)

    integer det_M, bound, n1,n2,n3,i,ind
    real(8) X0(3,1),X(3,1),Minv(3,3)
    real(8) work(3)
    integer ipiv(3),INFO

    det_M=M(1,1)*(M(2,2)*M(3,3)-M(3,2)*M(2,3)) &
         -M(1,2)*(M(2,1)*M(3,3)-M(3,1)*M(2,3)) &
         +M(1,3)*(M(2,1)*M(3,2)-M(3,1)*M(2,2))

    bound=3

    Minv(1:3,1:3)=1.0d0*M(1:3,1:3)
    call DGETRF( 3, 3, Minv, 3, ipiv, INFO )
    call DGETRI( 3, Minv, 3, ipiv, work, 3, INFO )

    allocate(pos_conv(1:3,1:det_M*size(pos(1,:))))

    do while (ind .lt. det_M*size(pos(1,:)))

       ind=0
       do n3=-bound,bound
          do n2=-bound,bound
             do n1=-bound,bound
                do i = 1, size(pos(1,:))
                   X0(1,1)=pos(1, i)+1.0d0*n1
                   X0(2,1)=pos(2, i)+1.0d0*n2
                   X0(3,1)=pos(3, i)+1.0d0*n3
                   X(1:3,1)=matmul(Minv(1:3,1:3),X0(1:3,1))
                   if ((X(1,1) .ge. 0.0d0) .and. (X(1,1) .lt. 1.0d0)) then
                      if ((X(2,1) .ge. 0.0d0) .and. (X(2,1) .lt. 1.0d0)) then
                         if ((X(3,1) .ge. 0.0d0) .and. (X(3,1) .lt. 1.0d0)) then
                            ind=ind+1
                            if (ind .gt. det_M*size(pos(1,:))) then
                               call vtutor%error('Error1 in find_conv_cell_atoms')
                            else
                               pos_conv(1:3,ind)=X0(1:3,1)
                            endif
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       if (ind .lt. det_M*size(pos(1,:))) then
          bound=bound*2
       elseif (ind .eq. det_M*size(pos(1,:))) then
          goto 22
       elseif (ind .gt. det_M*size(pos(1,:))) then
          call vtutor%error('Error2 in find_conv_cell_atoms')
       endif

       if (bound .ge. 1000) then
          call vtutor%error('Error3 in find_conv_cell_atoms: probably wrong structure')
       endif

    enddo
22  continue
  END SUBROUTINE find_conv_cell_atoms


  SUBROUTINE wap_sand_loc(cut_off_r,i_prim,alpha,isp,R_shift,NGX,NGY,NGZ,g_pr,grids,dV,k,kp,u_k,u_kp,gkkp)
    ! compute N * < psi_k| dV/dR | psi_kp >_sp, with WF normalized over V->infty, and  N number of prim cells in V.
    real(8), intent(in) :: cut_off_r
    integer, intent(in) :: i_prim
    integer, intent(in) :: alpha
    integer, intent(in) :: isp
    integer, intent(in) :: NGX
    integer, intent(in) :: NGY
    integer, intent(in) :: NGZ
    real(8),allocatable, intent(in) :: R_shift(:,:)
    real(8), allocatable, intent(in) :: g_pr(:,:)   ! just the size of that array is used in the present version
    real(8), allocatable, intent(in) :: grids(:,:,:,:)   ! last 2 fields of grids(1:size(centers(1,:)),1:indM,1:3)
    complex(8), allocatable, intent(in) :: dV(:,:,:,:)    ! a potential (derivative) on the grid in g_pr
    real(8), intent(in) :: k(1:3)
    real(8), intent(in) :: kp(1:3)
    complex(8), allocatable, intent(in) :: u_k(:)   ! periodic part of wave function: psi_k(r) =exp(i k r) u_k(r)
    complex(8), allocatable, intent(in) :: u_kp(:)
    complex(8), intent(out) :: gkkp
    integer NG(3,3)
    integer i,ir,n_grid
    complex(8) phase_at, phase_r,phase
    real(8) dot,d3r
    real(8) pi

    pi=acos(-1.0d0)
    n_grid=size(g_pr(1,:))

    NG= 0; NG(1,1)= NGX ; NG(2,2)= NGY; NG(3,3)= NGZ

    d3r=1.0d0/(NG(1,1)*NG(2,2)*NG(3,3))  ! should be V_sp/size(g_pr(1,:))= V_prim/(NGx*NGy*NGz) but because of normalization of WF, it is just 1/(NGx*NGy*NGz)

    gkkp=cmplx(0.0d0)

    do i=1,n_grid
       ! phase average over equal distance r points around center_in
       ! this is one way to average. may be not the best
       phase_r=cmplx(0.0d0)
       ir=1
       do while (grids(i_prim,i,ir,1) .lt. cut_off_r)
          dot=dot_product(kp(1:3)-k(1:3), grids(i_prim,i,ir,1:3))
          phase=cmplx(0.0d0,2.0d0*pi*dot)
          phase=exp(phase)
          phase_r=phase_r+phase
          ir=ir+1
       enddo
       if (ir .gt. 1) then
          phase_r=phase_r/(1.0d0*(ir-1))
       else
          call vtutor%error('Error in wap_sand_loc')
       endif

       gkkp=gkkp+conjg(u_k(i))*dV(i,isp,i_prim,alpha)*u_kp(i)*phase_r
    enddo

    gkkp=gkkp*d3r

    ! atomic phase
    dot=dot_product(kp(1:3)-k(1:3),R_shift(1:3,i_prim))
    phase_at=cmplx(0.0d0,2.0d0*pi*dot)
    phase_at=exp(phase_at)

    gkkp=phase_at*gkkp

  END SUBROUTINE wap_sand_loc


  SUBROUTINE wap_sand(NG,center_in,center_out,g_pr,grid,first,dV,k,kp,u_k,u_kp,gkkp)
    ! compute N * < psi_k | dV/dR | psi_kp >_sp, with WF normalized over V->infty, and  N number of prim cells in V.
    integer, intent(in) :: NG(3,3)
    real(8), intent(in) :: center_in(1:3)          ! reduce coord. wrt primitive/current lattice of the atom for which dV/dR_center is given
    real(8), intent(in) :: center_out(1:3)         ! reduce coord. wrt primitive/current lattice of the atom for which < psi_k | dV/dR |psi_kp >_sp is wanted.
    ! centers_in-centers_out must be a primitive lattice vector.
    real(8), allocatable, intent(in) :: g_pr(:,:)   ! just the size of that array is used in the present version
    real(8), allocatable, intent(in) :: grid(:,:)   ! last 2 fields of grids(1:size(centers(1,:)),1:indM,1:3)
    integer, allocatable, intent(in) :: first(:)    ! last field of  allocate(firsts(1:size(centers(1,:)),1:size(g_pr(1,:))+1))
    complex(8), allocatable, intent(in) :: dV(:)    ! a potential (derivative) on the grid in g_pr
    real(8), intent(in) :: k(1:3)
    real(8), intent(in) :: kp(1:3)
    complex(8), allocatable, intent(in) :: u_k(:)   ! periodic part of wave function: psi_k(r) =exp(i k r) u_k(r)
    complex(8), allocatable, intent(in) :: u_kp(:)
    complex(8), intent(out) :: gkkp

    integer i,ir
    complex(8) phase_at, phase_r,phase
    real(8) dot,d3r
    real(8) pi

    pi=acos(-1.0d0)

    d3r=1.0d0/(NG(1,1)*NG(2,2)*NG(3,3))  ! should be V_sp/size(g_pr(1,:))= V_prim/(NGx*NGy*NGz) but because of normalization of WF, it is just 1/(NGx*NGy*NGz)

    gkkp=cmplx(0.0d0)
    do i=1,size(g_pr(1,:))

       ! phase average over equal distance r points around center_in
       ! this is one way to average. may be not the best
       phase_r=cmplx(0.0d0)
       do ir=first(i),first(i+1)-1
          dot=dot_product(kp(1:3)-k(1:3), grid(ir,1:3))
          phase=cmplx(0.0d0,2.0d0*pi*dot)
          phase=exp(phase)
          phase_r=phase_r+phase/(1.0d0*(first(i+1)-first(i)))
       enddo

       gkkp=gkkp+conjg(u_k(i))*dV(i)*u_kp(i)*phase_r

    enddo

    gkkp=gkkp*d3r

    ! atomic phase
    dot=dot_product(kp(1:3)-k(1:3),center_out(1:3)-center_in(1:3))
    phase_at=cmplx(0.0d0,2.0d0*pi*dot)
    phase_at=exp(phase_at)

    gkkp=phase_at*gkkp

  END SUBROUTINE wap_sand


  SUBROUTINE get_u(u_k,prim_ind,WAVES,NK,spinor_k,spin_k,b_k)
    COMPLEX(8), ALLOCATABLE, intent(inout) :: u_k(:)
    integer(8), allocatable, intent(in) :: prim_ind(:,:)
    complex(8), allocatable, intent(in) :: WAVES(:,:,:,:,:,:,:)
    integer, intent(in) :: NK
    integer, intent(in) :: spinor_k
    integer, intent(in) :: spin_k
    integer, intent(in) :: b_k
    integer i,i1,i2,i3,ind

    do i=1,size(prim_ind(1,:))
       i1=prim_ind(1,i)  ! to by pass that step would be nice.
       ! But if the ordering of points is implicit it is dangerous
       i2=prim_ind(2,i)
       i3=prim_ind(3,i)
       u_k(i)=WAVES(i1,i2,i3,spinor_k, spin_k, b_k, NK)
    enddo

  END SUBROUTINE get_u


  SUBROUTINE getwave_sand_dV(NG,M,Nsp)
    integer, intent(in) :: NG(3,3) ! discretization matrix of the current(primitive) cell
    integer, intent(in) :: M(3,3)  ! primitive to conventionnal cell passage matrix
    integer, intent(in) :: Nsp(3,3) ! conventionnal to supercell passage matrix

  END SUBROUTINE getwave_sand_dV


  SUBROUTINE read_INCAR_wap(M,Nsp,NdV,mode_dV,epsilon_r,cut_off_r,mode_gkkp,centers)
    ! to be replace by keywords in INCAR
    integer, intent(out) :: M(3,3)  ! primitive to conventionnal cell passage matrix
    integer, intent(out) :: Nsp(3,3) ! conventionnal to supercell passage matrix
    integer, intent(out) :: NdV      ! number of atoms for which dV/dR is known
    character*4, intent(out) :: mode_dV
    real(8), intent(out) :: epsilon_r
    real(8), intent(out) :: cut_off_r
    character*4, intent(out) :: mode_gkkp
    real(8),allocatable, intent(out) :: centers(:,:) ! atoms where the potential derivatives are known
    NAMELIST /wap/ M,Nsp,NdV,mode_dV,epsilon_r,cut_off_r,mode_gkkp
    NAMELIST /atoms/ centers

    OPEN(11,FILE='INCAR_wap',STATUS='OLD')
    read(11,NML=wap)

    rewind(11)
    allocate(centers(1:3,1:NdV))
    read(11,NML=atoms)

    close(11)

  END SUBROUTINE read_INCAR_wap


#ifdef VASP_HDF5

  !> @brief write the vtot data to a hdf5 file
  !> This should include all the information that will be needed
  !> in the next steps of the electron-phonon calculation, including metadata
  SUBROUTINE WRITE_VTOT_HDF5(FILEID,POT,CCQIJ,CCDIJ,QTOT,CCRIJ)
    use vhdf5
    integer(HID_T), intent(in) :: FILEID
    COMPLEX(q),intent(in) :: POT(:,:,:,:)
    COMPLEX(q),intent(in) :: CCQIJ(:,:,:,:)
    COMPLEX(q),intent(in) :: CCDIJ(:,:,:,:)
    REAL(q),intent(in) :: QTOT(:,:,:)
    COMPLEX(q),intent(in) :: CCRIJ(:,:,:,:,:)
    ! local variables
    INTEGER(HID_T):: GROUPID

    ! open vtot group
    VH5_CHECK(VH5_GROUP_OPEN_OR_CREATE(FILEID, 'vtot', GROUPID))

    ! some dimensions of the arrays
    VH5_CHECK(VH5_WRITE(GROUPID,'ncdij', SIZE(POT,4)))
    VH5_CHECK(VH5_WRITE(GROUPID,'fft_mesh', [SIZE(POT,1),SIZE(POT,2),SIZE(POT,3)]))
    VH5_CHECK(VH5_WRITE(GROUPID,'lmdim', SIZE(CCQIJ,1)))
    VH5_CHECK(VH5_WRITE(GROUPID,'ldim', SIZE(QTOT,1)))
    VH5_CHECK(VH5_WRITE(GROUPID,'natoms', SIZE(CCQIJ,3)))

    ! the arrays
    VH5_CHECK(VH5_WRITE(GROUPID,'pot', POT))
    VH5_CHECK(VH5_WRITE(GROUPID,'ccqij', CCQIJ))
    VH5_CHECK(VH5_WRITE(GROUPID,'ccdij', CCDIJ))
    VH5_CHECK(VH5_WRITE(GROUPID,'ccrij', CCRIJ))
    VH5_CHECK(VH5_WRITE(GROUPID,'qtot', QTOT))

    ! close vtot group
    VH5_CHECK(VH5_GROUP_CLOSE_WRITING(GROUPID))
  END SUBROUTINE WRITE_VTOT_HDF5


  !> @brief read the contents of vtot from an hdf5 file
  SUBROUTINE READ_VTOT_HDF5(FILEID,POT,CCQIJ,CCDIJ,QTOT,CCRIJ)
    use vhdf5
    INTEGER(HID_T) :: FILEID
    COMPLEX(q) :: POT(:,:,:,:)
    COMPLEX(q) :: CCQIJ(:,:,:,:)
    COMPLEX(q) :: CCDIJ(:,:,:,:)
    REAL(q) :: QTOT(:,:,:)
    COMPLEX(q) :: CCRIJ(:,:,:,:,:)
    ! local variables
    INTEGER(HID_T) :: GROUPID

    ! open vtot group
    VH5_CHECK(VH5_GROUP_OPEN_OR_CREATE(FILEID, 'vtot', GROUPID))

    ! read the arrays
    VH5_CHECK(VH5_READ(GROUPID,'pot', POT))
    VH5_CHECK(VH5_READ(GROUPID,'ccqij', CCQIJ))
    VH5_CHECK(VH5_READ(GROUPID,'ccdij', CCDIJ))
    VH5_CHECK(VH5_READ(GROUPID,'ccrij', CCRIJ))
    VH5_CHECK(VH5_READ(GROUPID,'qtot', QTOT))

    ! close vtot group
    VH5_CHECK(VH5_GROUP_CLOSE(GROUPID))

  END SUBROUTINE READ_VTOT_HDF5


  !> @brief read the DVDU file from hdf5 produced by phelel
  SUBROUTINE READ_PHELEL_HDF5(NGRIDPOINTS,DET_P2S,NCDIJ,LMDIM,&
       NATOMS_PC,NATOMS_SC,PHONON_NATOMS_SC,&
       GRID_POINT,LATTICE_POINT,PRIMITIVE_POSITIONS,&
       SUPERCELL_MATRIX,SUPERCELL_LATTICE,SUPERCELL_POSITIONS,ATOM_INDICES_IN_DERIVATIVES,&
       PRIMITIVE_MATRIX,PRIMITIVE_LATTICE,DVDU,CCDDIJ,CCDQIJ,CCDIJ,CCQIJ,&
       BORN_EFF_CHARGES,DIELECTRIC_TENSOR,&
       PHONON_MASSES,FORCE_CONSTANTS,FILENAME)
    USE vhdf5
    INTEGER :: NGRIDPOINTS, DET_P2S, NCDIJ, LMDIM, NATOMS_PC, NATOMS_SC
    INTEGER :: PHONON_NATOMS_SC
    REAL(q) :: GRID_POINT(3,NGRIDPOINTS)
    REAL(q) :: LATTICE_POINT(3,DET_P2S)
    REAL(q) :: PRIMITIVE_POSITIONS(3,NATOMS_PC)
    INTEGER :: ATOM_INDICES_IN_DERIVATIVES(NATOMS_PC)
    REAL(q) :: SUPERCELL_POSITIONS(3,NATOMS_SC)
    REAL(q) :: SUPERCELL_MATRIX(3,3)
    REAL(q) :: PRIMITIVE_MATRIX(3,3)
    REAL(q) :: SUPERCELL_LATTICE(3,3)
    REAL(q) :: PRIMITIVE_LATTICE(3,3)
    REAL(q) :: PHONON_SUPERCELL_MATRIX(3,3)
    REAL(q) :: PHONON_MASSES(NATOMS_PC)
    REAL(q) :: FORCE_CONSTANTS(3,3,PHONON_NATOMS_SC,PHONON_NATOMS_SC)
    REAL(q) :: BORN_EFF_CHARGES(3,3,NATOMS_PC)
    REAL(q) :: DIELECTRIC_TENSOR(3,3)
    COMPLEX(q) :: CCDIJ(LMDIM, LMDIM, NATOMS_PC, NCDIJ)
    COMPLEX(q) :: CCQIJ(LMDIM, LMDIM, NATOMS_PC, NCDIJ)
    COMPLEX(q) :: CCDDIJ(LMDIM, LMDIM, NATOMS_SC, 3, NATOMS_PC, NCDIJ)
    COMPLEX(q) :: CCDQIJ(LMDIM, LMDIM, NATOMS_SC, 3, NATOMS_PC, NCDIJ)
    COMPLEX(q) :: DVDU(NGRIDPOINTS,3,NATOMS_PC,NCDIJ)
    CHARACTER(LEN=*), OPTIONAL :: FILENAME
    ! local variables
    INTEGER :: I,IERR
    REAL(q) :: LATTICE_TRANSPOSE(3,3)
    INTEGER(HID_T):: FILEID
    CHARACTER(LEN=FILENAME_MAX_LEN) :: MY_FILENAME

    MY_FILENAME = 'phelel_params.hdf5'; IF (PRESENT(FILENAME)) MY_FILENAME = FILENAME
    VH5_CHECK(VH5_FILE_OPEN_READ(MY_FILENAME, FILEID))

    ! information regarding the grids
    VH5_CHECK(VH5_READ(FILEID,'grid_point',GRID_POINT))
    VH5_CHECK(VH5_READ(FILEID,'lattice_point',LATTICE_POINT))

    ! information regarding the local potential
    VH5_CHECK(VH5_READ(FILEID,'dVdu',DVDU))

    ! information regarding the projectors
    VH5_CHECK(VH5_READ(FILEID,'Dij',CCDIJ))
    VH5_CHECK(VH5_READ(FILEID,'qij',CCQIJ))
    VH5_CHECK(VH5_READ(FILEID,'dDijdu',CCDDIJ))
    VH5_CHECK(VH5_READ(FILEID,'dqijdu',CCDQIJ))

    ! information regarding the lattice parameters
    VH5_CHECK(VH5_READ(FILEID,'supercell_matrix',LATTICE_TRANSPOSE))
    SUPERCELL_MATRIX = TRANSPOSE(LATTICE_TRANSPOSE)
    VH5_CHECK(VH5_READ(FILEID,'supercell_lattice',LATTICE_TRANSPOSE))
    SUPERCELL_LATTICE = TRANSPOSE(LATTICE_TRANSPOSE)
    VH5_CHECK(VH5_READ(FILEID,'supercell_positions',SUPERCELL_POSITIONS))
    VH5_CHECK(VH5_READ(FILEID,'atom_indices_in_derivatives',ATOM_INDICES_IN_DERIVATIVES))

    IERR = VH5_READ(FILEID,'primitive_matrix',LATTICE_TRANSPOSE)
    IF (IERR/=0) THEN
       PRIMITIVE_MATRIX = 0
       PRIMITIVE_MATRIX(1,1) = 1;
       PRIMITIVE_MATRIX(2,2) = 1;
       PRIMITIVE_MATRIX(3,3) = 1;
    ELSE
       PRIMITIVE_MATRIX = TRANSPOSE(LATTICE_TRANSPOSE)
    ENDIF
    VH5_CHECK(VH5_READ(FILEID,'primitive_lattice',LATTICE_TRANSPOSE))
    PRIMITIVE_LATTICE = TRANSPOSE(LATTICE_TRANSPOSE)
    VH5_CHECK(VH5_READ(FILEID,'primitive_positions',PRIMITIVE_POSITIONS))

    ! information regarding the force constants
    VH5_CHECK(VH5_READ(FILEID,'phonon_supercell_matrix',PHONON_SUPERCELL_MATRIX))
    VH5_CHECK(VH5_READ(FILEID,'force_constants',FORCE_CONSTANTS))
    VH5_CHECK(VH5_READ(FILEID,'primitive_masses',PHONON_MASSES))

    ! read the born effective charges
    IERR = VH5_READ(FILEID,'born_effective_charges',BORN_EFF_CHARGES)
    IF (IERR/=0) BORN_EFF_CHARGES = 0
    ! read the dielectric tensor
    IERR = VH5_READ(FILEID,'dielectric_constant',DIELECTRIC_TENSOR)
    IF (IERR/=0) THEN
       DIELECTRIC_TENSOR = 0
       DIELECTRIC_TENSOR(1,1) = 1
       DIELECTRIC_TENSOR(2,2) = 1
       DIELECTRIC_TENSOR(3,3) = 1
    ENDIF
    VH5_CHECK(VH5_FILE_CLOSE(FILEID))
  END SUBROUTINE READ_PHELEL_HDF5
#endif VASP_HDF5

  !> @brief read the dimensions of the DVDU file from hdf5 produced by phelel
  SUBROUTINE READ_PHELEL_DIMENSIONS_HDF5(NGRIDPOINTS,DET_P2S,NCDIJ,LMDIM,&
       NATOMS_PC,NATOMS_SC,FFT_MESH,&
       PRIMITIVE_LATTICE, &
       PHONON_NATOMS_SC,PHONON_DET_P2S,FILENAME)
    USE vhdf5
    USE tutor, only: vtutor, isError, ElectronPhononRequiresHDF5
    INTEGER :: NGRIDPOINTS !< PROD(FFT_MESH)*DET_P2S
    INTEGER :: DET_P2S !< Determinant of the supercell matrix WRONG = vol supercell/vol primitive cell
    INTEGER :: NCDIJ
    INTEGER :: LMDIM
    INTEGER :: FFT_MESH(3)
    REAL(q) :: PRIMITIVE_LATTICE(3,3)
    INTEGER :: NATOMS_PC !< Number of atoms in the primitive cell
    INTEGER :: NATOMS_SC !< Number of atoms in the supercell
    INTEGER :: PHONON_NATOMS_SC !< Number of atoms in the phonon supercell
    INTEGER :: PHONON_DET_P2S
    CHARACTER(LEN=*), OPTIONAL :: FILENAME
#ifdef VASP_HDF5
    ! local variables
    REAL(q) :: LATTICE_TRANSPOSE(3,3)
    INTEGER(HID_T) :: FILEID
    INTEGER :: RANK, IERR
    INTEGER :: DIMENSIONS(7) ! the variables I will read have rank smaller than 5
    CHARACTER(LEN=FILENAME_MAX_LEN) :: MY_FILENAME

    MY_FILENAME = 'phelel_params.hdf5'; IF (PRESENT(FILENAME)) MY_FILENAME = FILENAME
    IERR = VH5_FILE_OPEN_READ(MY_FILENAME, FILEID)
    IF (IERR/=0) THEN
       call vtutor%error('Failed to open '//TRIM(MY_FILENAME)//'.  Please check if the file is present in the current path.')
    ENDIF

    NGRIDPOINTS = 0; NATOMS_PC = 0; NCDIJ = 0; DET_P2S = 0
    ! first we read the dimensions
    VH5_CHECK(VH5_GET_DIMENSIONS(FILEID,'dVdu',DIMENSIONS,RANK))
    NGRIDPOINTS = DIMENSIONS(2)
    NATOMS_PC = DIMENSIONS(4)
    NCDIJ = DIMENSIONS(5)
    VH5_CHECK(VH5_GET_DIMENSIONS(FILEID,'lattice_point',DIMENSIONS,RANK))
    DET_P2S = DIMENSIONS(2)  ! it is the number of lattice points in the supercell= vol supercell/vol primitive cell
    VH5_CHECK(VH5_GET_DIMENSIONS(FILEID,'dDijdu',DIMENSIONS,RANK))
    LMDIM = DIMENSIONS(2)
    NATOMS_SC = DIMENSIONS(4)
    ! for the phonons
    VH5_CHECK(VH5_GET_DIMENSIONS(FILEID,'force_constants',DIMENSIONS,RANK))
    PHONON_NATOMS_SC = DIMENSIONS(3)
    PHONON_DET_P2S=PHONON_NATOMS_SC/NATOMS_PC

    ! read the size of the fft_grid
    VH5_CHECK(VH5_READ(FILEID,'FFT_mesh',FFT_MESH))

    ! read the primitive lattice
    VH5_CHECK(VH5_READ(FILEID,'primitive_lattice',LATTICE_TRANSPOSE))
    PRIMITIVE_LATTICE = TRANSPOSE(LATTICE_TRANSPOSE)

    VH5_CHECK(VH5_FILE_CLOSE(FILEID))
#else
    call vtutor%write(isError,ElectronPhononRequiresHDF5)
#endif
  END SUBROUTINE READ_PHELEL_DIMENSIONS_HDF5

  SUBROUTINE READ_PHELEL_PHONON_DIMENSIONS_HDF5(NATOMS_PC,PHONON_NATOMS_SC,PHONON_DET_P2S,FILENAME)
    USE vhdf5
    USE tutor, only: vtutor, isError, ElectronPhononRequiresHDF5
    INTEGER :: NATOMS_PC !< Number of atoms in the primitive cell
    INTEGER :: PHONON_NATOMS_SC !< Number of atoms in the phonon supercell
    INTEGER :: PHONON_DET_P2S !< Determinant of the supercell matrix WRONG = vol supercell/vol primitive cell
    CHARACTER(LEN=*), OPTIONAL :: FILENAME

#ifdef VASP_HDF5
    ! local variables
    INTEGER(HID_T) :: FILEID
    INTEGER :: RANK, IERR
    INTEGER :: DIMENSIONS(7) ! the variables I will read have rank smaller than 5
    CHARACTER(LEN=FILENAME_MAX_LEN) :: MY_FILENAME

    MY_FILENAME = 'phelel_params.hdf5'; IF (PRESENT(FILENAME)) MY_FILENAME = FILENAME
    IERR = VH5_FILE_OPEN_READ(MY_FILENAME, FILEID)
    IF (IERR/=0) THEN
       call vtutor%error('Failed to open '//TRIM(MY_FILENAME)//'.  Please check if the file is present in the current path.')
    ENDIF

    ! read NATOMS_PC
    VH5_CHECK(VH5_GET_DIMENSIONS(FILEID,'dVdu',DIMENSIONS,RANK))
    NATOMS_PC = DIMENSIONS(4)

    ! first we read the dimensions
    VH5_CHECK(VH5_GET_DIMENSIONS(FILEID,'force_constants',DIMENSIONS,RANK))
    PHONON_NATOMS_SC = DIMENSIONS(3)
    PHONON_DET_P2S=PHONON_NATOMS_SC/NATOMS_PC

    VH5_CHECK(VH5_FILE_CLOSE(FILEID))
#endif
    END SUBROUTINE READ_PHELEL_PHONON_DIMENSIONS_HDF5

  !@brief Transform the potential to a grid compatible with the WFs
  SUBROUTINE DVDU_TRANSFER_GRID(WDES,MAP_TO_WF)
    USE wave
    TYPE (wavedes) :: WDES
    INTEGER :: NGRIDPOINTS
    INTEGER :: NATOMS_SC
    INTEGER :: NCDIJ
    INTEGER,ALLOCATABLE :: MAP_TO_WF(:)
    ! local variables
    INTEGER :: NX,NY,NZ
    INTEGER :: NGX,NGY,NGZ
    INTEGER :: NC, NPWS, NPW, IND
    INTEGER :: ISPINOR
    NGX=WDES%GRID%NGX
    NGY=WDES%GRID%NGY
    NGZ=WDES%GRID%NGZ
    NPW=NGX*NGY*NGZ
    ALLOCATE(MAP_TO_WF(NPW*WDES%NRSPINORS))
    IF (WDES%GRID%RL%NFAST==1) THEN
       DO ISPINOR=1,WDES%NRSPINORS
          NPWS = NPW*(ISPINOR-1)
          DO NC=1,WDES%GRID%RL%NCOL
             NY=WDES%GRID%RL%I2(NC)
             NZ=WDES%GRID%RL%I3(NC)
             DO NX=1,WDES%GRID%RL%NROW
                IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
                MAP_TO_WF(IND+NPWS)=NX+(NC-1)*WDES%GRID%RL%NROW+NPWS
             ENDDO
          ENDDO
       ENDDO
    ELSEIF (WDES%GRID%RL%NFAST==3) THEN
       DO ISPINOR=1,WDES%NRSPINORS
          NPWS = NPW*(ISPINOR-1)
          DO NC=1,WDES%GRID%RL%NCOL
             NX=WDES%GRID%RL%I2(NC)
             NY=WDES%GRID%RL%I3(NC)
             DO NZ=1,WDES%GRID%RL%NROW
                IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
                MAP_TO_WF(IND+NPWS)=NX+(NC-1)*WDES%GRID%RL%NROW+NPWS
             ENDDO
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE DVDU_TRANSFER_GRID


  SUBROUTINE getwave_local_potential_CVTOT(WDES, GRIDC, CVTOT, POT)
    USE prec
    USE mgrid
    USE wave
    USE iso_c_binding

    TYPE(wavedes), intent(in) :: WDES
    TYPE(grid_3d), intent(in) :: GRIDC
    COMPLEX(q),TARGET, intent(in) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
    COMPLEX(8), intent(inout) :: POT(GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ, WDES%NCDIJ)
    ! Local variables
    RGRID, POINTER :: CVTOT_LOCAL(:)
    integer :: IDX, ISP, NC, NX, NY, NZ

    DO ISP = 1, WDES%NCDIJ
       CALL C_F_POINTER(C_LOC(CVTOT(:,ISP)),CVTOT_LOCAL,[GRIDC%RL%NP])
       DO NZ = 1, GRIDC%NGZ
          IF (GRIDC%RL%NFAST == 1) THEN
             DO NC = 1, GRIDC%RL%NCOL
                IF (GRIDC%RL%I3(NC) /= NZ) CYCLE
                NY = GRIDC%RL%I2(NC)
                DO NX = 1, GRIDC%NGX
                   IDX = NX + (NC - 1) * GRIDC%NGX
                   POT(NX, NY, NZ, ISP) = dcmplx(CVTOT_LOCAL(IDX))
                ENDDO
             ENDDO
          ELSEIF (GRIDC%RL%NFAST == 2) THEN
             CALL vtutor%error('error RL%NFAST=2 not implemented')
          ELSEIF (GRIDC%RL%NFAST == 3) THEN
             DO NC = 1, GRIDC%RL%NCOL
                NX = GRIDC%RL%I2(NC)
                NY = GRIDC%RL%I3(NC)
                IDX = NZ + (NC - 1) * GRIDC%NGZ
                POT(NX, NY, NZ, ISP) = dcmplx(CVTOT_LOCAL(IDX))
             ENDDO
          ENDIF
       ENDDO

    ENDDO
    CALLMPI(M_sum_d(GRIDC%COMM, POT, 2 * GRIDC%NGX * GRIDC%NGY * GRIDC%NGZ * WDES%NCDIJ))

  END SUBROUTINE getwave_local_potential_CVTOT

END MODULE wap
