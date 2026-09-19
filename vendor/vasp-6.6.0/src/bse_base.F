#include "symbol.inc"
#include "bse.inc"
MODULE bse_base
  USE constant
  USE prec
  USE chi_glb
  USE local_field
  USE scala
  USE wave_high
  USE bse_struct
  USE msupercell
#ifdef VASP_HDF5
    USE vhdf5
#endif

  IMPLICIT NONE

#ifdef single_prec_bse
  GDEFS, ALLOCATABLE, TARGET :: BSEMATRIX(:,:)
#else
  GDEF, ALLOCATABLE, TARGET :: BSEMATRIX(:,:)
#endif

!> the BSE matrix is preset with this value (used as an indicator for failures)
#ifdef single_prec_bse
  REAL(qs),PARAMETER:: MAGIC=-123456_qs
  INTEGER, PARAMETER :: mq =SELECTED_REAL_KIND(5)
#endif
#ifdef double_prec_bse
  REAL(q),PARAMETER::  MAGIC=-123456_q
  INTEGER, PARAMETER :: mq =SELECTED_REAL_KIND(10)
#endif
  CONTAINS

#ifdef scaLAPACK
!**********************************************************************
!
!> Writes out exciton wavefunctions
!
!> Only scaLAPACK-aware version is implemented.
!**********************************************************************

  SUBROUTINE EXCITON_WF(T_INFO,INFO,P,WHF,BSEVECTORS,LATT_CUR,LMDIM,IO,NCV,DESC,BSE_INDEX,BDS,ISP_LOW,ISP_HIGH,OMEGAMAX)

    USE pseudo, only: potcar
    USE fileio, ONLY: OUTCHG_RGRID
    USE iso_c_binding, ONLY: c_loc, c_f_pointer

    IMPLICIT NONE

    TYPE (type_info) :: T_INFO
    TYPE (info_struct) :: INFO
    TYPE (potcar) :: P(T_INFO%NTYP)
    TYPE (wavespin) :: WHF
#ifdef single_prec_bse
    GDEFS :: BSEVECTORS(:) !< BSE eigenvectors
#else
    GDEF :: BSEVECTORS(:) !< BSE eigenvectors
#endif
    TYPE (latt) :: LATT_CUR
    INTEGER :: IU6, LMDIM
    TYPE (in_struct) :: IO
    INTEGER :: NCV
    INTEGER :: DESC(DLEN_)
    TYPE(bse_matrix_index) :: BSE_INDEX
    TYPE(banddesc) :: BDS(WHF%WDES%ISPIN)
    INTEGER :: ISP_LOW, ISP_HIGH
    REAL(q) :: OMEGAMAX

    ! local variables
    TYPE(supercell),POINTER :: S=>NULL()
    TYPE (wavedes1), TARGET :: WDESK1, WDESK3
    TYPE(banddesc) :: BD
    TYPE (wavefun1), ALLOCATABLE :: W1(:), W3(:)
    TYPE (wavefun1) :: EWF
    RGRID, ALLOCATABLE :: CHG(:)
    COMPLEX(q), ALLOCATABLE, TARGET :: CWORK(:)
    GDEF, POINTER :: GWORK(:)
    COMPLEX(q), ALLOCATABLE  :: CW(:)
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: BSEVECTOR(:)
#else
    GDEF, ALLOCATABLE :: BSEVECTOR(:)
#endif
    REAL(q) :: R(3)
    INTEGER :: NG(3)
    INTEGER :: NK, N, NS, ISP
    INTEGER :: NODE
    INTEGER :: NR, NC, NV
    INTEGER :: NPROW,NPCOL,MYROW,MYCOL,LAMBDA
    INTEGER :: VBANDS, CBANDS
    INTEGER :: K1,K3, NV_INTO_TOT, NC_INTO_TOT
    INTEGER :: N1,N3,NCV13,NP,NQ, NK1_IN_KPOINTS_FULL_ORIG
    INTEGER :: DESC_V(DLEN_), ICTXC, IFAIL
    COMPLEX(q) :: PHI
    GDEF :: EPHI
    INTEGER :: NKL, NKG_FIRST
    LOGICAL :: LFIXHOLE, LFIXELECTRON
    GDEF, ALLOCATABLE :: CRHOLM(:)

    PROFILING_START('exciton_wf')

    ! setup supercell and fft
    CALL INIT_SUPERCELL(S, T_INFO, LATT_CUR, KPOINTS_FULL%NKPX, KPOINTS_FULL%NKPY,  KPOINTS_FULL%NKPZ, KPOINTS_FULL%B , IO%IU6)
    CALL INIT_SUPERCELL_FFT_WF(S, P, LMDIM,  WHF%WDES, LATT_CUR, IO%IU6)

#ifdef VASP_HDF5
    IF (IO%WRT_DENSITY%LCHARGH5) THEN
       CALL VH5_WRITE_EXCITON_START(IH5OUTFILEID, "/results/linear_response", S%GRID%NGPTAR(1),S%GRID%NGPTAR(2),S%GRID%NGPTAR(3), NBSEEIG)
    END IF
#endif

    CALL NEWWAV(EWF, S%WDES1, .TRUE.)

    CBANDS = SIZE(BSE_INDEX%INDEX,2)
    VBANDS = SIZE(BSE_INDEX%INDEX,1)

    ALLOCATE(CRHOLM(S%AUG_DES%NPROD*WHF%WDES%NRSPINORS))
    ALLOCATE(CWORK(S%GRID%MPLWV))
    ALLOCATE(CHG(S%GRID%MPLWV))
    ALLOCATE(CW(S%WDES1%NRPLWV))
    ALLOCATE(BSEVECTOR(NCV))
    ALLOCATE(W3(CBANDS))
    ALLOCATE(W1(VBANDS))

    CALL c_f_pointer(c_loc(CWORK),GWORK,[m_ SIZE(CWORK)])

    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE calculating exciton wavefunctions'

    CALL SETWDES(WHF%WDES,WDESK1,0)
    CALL SETWDES(WHF%WDES,WDESK3,0)

    DO N=1,VBANDS
       CALL NEWWAV(W1(N), WDESK1, .TRUE.)
    ENDDO

    DO N=1,CBANDS
       CALL NEWWAV(W3(N), WDESK3, .TRUE.)
    ENDDO

    ! init context for gathering eigenvectors
    CALL BLACS_GET(0, 0, ICTXC)
    CALL BLACS_GRIDINIT(ICTXC, 'R', 1, 1)
    CALL BLACS_GRIDINFO(ICTXC, NPROW, NPCOL, MYROW, MYCOL)


    ! init descriptor for gathering eigenvectors
    IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
      CALL DESCINIT(DESC_V, NCV, 1, NCV, 1, 0, 0, ICTXC, NCV, IFAIL)
    ELSE
      DESC_V(2)=-1
    ENDIF

    IF (ANY(BSEHOLE /= -123456.0_q) .AND. ANY(BSEELECTRON /= -123456.0_q)) THEN
       CALL vtutor%alert("Both electron and hole positions cannot be fixed: ignoring BSEELECTRON.")
    ENDIF

    IF (ANY(BSEHOLE /= -123456.0_q)) THEN
      R=BSEHOLE
      LFIXHOLE = .TRUE.
      LFIXELECTRON = .FALSE.
    ELSE
      R=BSEELECTRON
      LFIXHOLE= .FALSE.
      LFIXELECTRON = .TRUE.
    ENDIF

    ! find grid point of fixed position
    NG=FLOOR(R*WHF%WDES%GRID%NGPTAR)+1
    ! consider periodicity
    NG=MOD(NG,WHF%WDES%GRID%NGPTAR)
    ! find position corresponding to the grid point
    R=(1._q/WHF%WDES%GRID%NGPTAR)*(NG-1)

    IF(WHF%WDES%GRID%RL%NFAST /= 1)  THEN
      CALL vtutor%bug("The order of PW indices is incompatible with this routine", __FILE__, __LINE__)
    ENDIF

    ! index of the grid point corresponding to vector R
    NR = NG(1) + (NG(2)-1)*WHF%WDES%GRID%NGX + (NG(3)-1)*WHF%WDES%GRID%NGX*WHF%WDES%GRID%NGY

    IF (LFIXHOLE) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,'(/1X,A,3F10.6/)') 'hole position is fixed at:', R
    ELSE
       IF (IO%IU6>=0) WRITE(IO%IU6,'(/1X,A,3F10.6/)') 'electron position is fixed at:', R
    ENDIF

    DO LAMBDA=1,NBSEEIG
       EWF%CW = zero
       EWF%CR = zero
       EWF%CPROJ = zero
       BSEVECTOR = zerobse

       ! the first local k-point (k-point on this rank)
       NKL = WHF%WDES%NB_LOW
       ! the first global k-point (i.e. k-point on any rank)
       NKG_FIRST = 1

       ! gather distributed eigenvector LAMBDA onto master rank
       CALL BLACS_GRIDINFO(DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)
       IF(MYROW /= -1 .OR. MYCOL /= -1) THEN
          CALL PGGEMR2Dbse(NCV, 1, BSEVECTORS, 1, LAMBDA, DESC, BSEVECTOR, 1, 1, DESC_V, DESC(CTXT_))
       ENDIF

       ! bcast eigenvector LAMBDA to all other ranks
       CALL M_bcastbse_g_from(WHF%WDES%COMM_INTER, BSEVECTOR, SIZE(BSEVECTOR), 1)

       DO ISP=ISP_LOW,ISP_HIGH
          BD=BDS(ISP)

          DO WHILE (NKG_FIRST <= WHF%WDES%NKPTS)
              CW = zzero

              ! gather orbitals for k-points
              ! loop from the first global k-point to the last global k-point
              DO NK = NKG_FIRST,MIN(WHF%WDES%NKPTS,NKG_FIRST + WHF%WDES%NB_PAR - 1)
                 CALL SETWDES(WHF%WDES,WDESK1,NK)
                 CALL SETWDES(WHF%WDES,WDESK3,NK)

                 ! set WDESK for all W1 and W3
                 FORALL (N=1:SIZE(BSE_INDEX%INDEX,1)) W1(N)%WDES1 => WDESK1
                 FORALL (N=1:SIZE(BSE_INDEX%INDEX,2)) W3(N)%WDES1 => WDESK3

                 ! node where this k-point is computed
                 NODE = MOD(NK-1,WHF%WDES%NB_PAR)+1

                 ! only need WF of fixed particle on real grid
                 IF(LFIXHOLE) THEN
                    CALL W1_GATHER_KSEL_NOFFT(WHF, BD%CBMIN, BD%CBMAX, ISP, W3, NODE, NK)
                    CALL W1_GATHER_KSEL(WHF, BD%VBMIN, BD%VBMAX, ISP, W1, NK)
                 ELSE
                    CALL W1_GATHER_KSEL_NOFFT(WHF, BD%VBMIN, BD%VBMAX, ISP, W1, NODE, NK)
                    CALL W1_GATHER_KSEL(WHF, BD%CBMIN, BD%CBMAX, ISP, W3, NK)
                 ENDIF
              ENDDO

              ! k-point copmuted locally
              NK = NKL

              IF (NKL <= WHF%WDES%NKPTS) THEN
                 ! phi=2pi*k*r
                 PHI = CITPI*DOT_PRODUCT(WHF%WDES%VKPT(:,NK),R)
#ifndef gammareal
                 ! exp(ikr)
                 IF(LFIXHOLE) THEN
                    EPHI = EXP(-PHI)
                 ELSE
                    EPHI = EXP(PHI)
                 ENDIF
#else
                 EPHI = one
#endif
                 DO NV_INTO_TOT =  BD%VBMIN, BD%VBMAXpCORE
                    NV = NV_INTO_TOT - BD%VBMIN + 1
                    DO NC_INTO_TOT = BD%CBMIN, BD%CBMAX
                       NC = NC_INTO_TOT - BD%CBMIN + 1

                       IF (REAL(WHF%CELTOT(NC,NK,ISP)-WHF%CELTOT(NV,NK,ISP),q) < OMEGAMAX ) THEN
                         DO N=1,S%NGVECTOR(NK)
                            ! PW coefficient in supercell
                            NS=S%INDEX(N,NK)
                            NCV13 = BSE_INDEX%INDEX(NV,NC,NK,ISP)

                            IF(LFIXHOLE) THEN
                               IF(NV_INTO_TOT > BD%VBMAX) THEN
                                  CW(NS)=CW(NS)+W3(NC)%CW(N)*EPHI*BSEVECTOR(NCV13)
                               ELSE
                                  CW(NS)=CW(NS)+W3(NC)%CW(N)*CONJG(W1(NV)%CR(NR))*EPHI*BSEVECTOR(NCV13)
                               ENDIF
                            ELSE
                               CW(NS)=CW(NS)+CONJG(W3(NC)%CR(NR))*W1(NV)%CW(N)*EPHI*BSEVECTOR(NCV13)
                            ENDIF
                         ENDDO
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF

              ! next global k-point
              NKG_FIRST = NKG_FIRST + WHF%WDES%NB_PAR
              ! next local k-point
              NKL = NKL + WHF%WDES%NB_PAR

              ! add PWs for all k-points
              CALL M_sum_z8(WHF%WDES%COMM_INTER, CW, SIZE(CW,KIND=qi8))
              EWF%CW = EWF%CW + CW
          ENDDO ! k-points
       ENDDO ! spin

       ! FFT WF in supercell to real space
       CALL FFTWAV_W1(EWF)

       ! plane wave contribution
       CALL PW_CHARGE_TRACE(EWF%WDES1, GWORK(1), EWF%CR(1), EWF%CR(1))

       ! compute augmentation charge, false by default
       IF(LBSEWFAUG) THEN

          ! setup non-local projection operators
          IF(LBSEWFAUG) THEN
             ALLOCATE(S%NONL_S)
             CALL NONL_ALLOC(S%NONL_S,S%T_INFO,P,S%WDES, .FALSE.)
             ! set phase once and forever
             CALL SPHER(S%GRID, S%NONL_S, P, S%WDES, S%LATT_CUR,  1)
             CALL PHASE(S%WDES,S%NONL_S,1)
          ENDIF

          ! project psi onto all projectors
          CALL PROJ1(S%NONL_S,S%WDES1,EWF)

          ! add augmentation part
          IF (EWF%WDES1%LOVERL) THEN
             CALL DEPSUM_TWO_BANDS_RHOLM_TRACE(EWF%CPROJ(:),EWF%CPROJ(:), EWF%WDES1, S%AUG_DES, &
                  S%TRANS_MATRIX, CRHOLM, 1._q, EWF%WDES1%LOVERL)

             S%AUG_DES%RINPL=1._q ! multiplicator used by RACC0
             CALL RACC0_HF(S%FAST_AUG, S%AUG_DES, CRHOLM(1), GWORK(1))
          ENDIF
       ENDIF

       CHG = REAL(GWORK)

       ! writing exciton charge density into CHGCAR files fith suffix LAMBDA
#ifdef VASP_HDF5
       IF (IO%WRT_DENSITY%LCHARGH5) THEN
           CALL VH5_WRITE_EXCITON_STEP(IH5OUTFILEID, "/results/linear_response", &
              S%GRID%NGPTAR(1),S%GRID%NGPTAR(2),S%GRID%NGPTAR(3), LAMBDA, CWORK)
       END IF
#endif

       IF (IO%IU6>=0 .AND. .NOT. IO%WRT_DENSITY%LCHARGH5) THEN
          CALL DUMP_CHARGE_RGRID("CHG."//str(LAMBDA), S%GRID,CHG,IO%IU6)
       ENDIF

    ENDDO ! eigenvector

#ifdef VASP_HDF5
    IF (IO%WRT_DENSITY%LCHARGH5) THEN
       CALL VH5_WRITE_LATTICE_ION(IH5OUTFILEID, INFO%SZNAM1, S%T_INFO, LATT_CUR%SCALE, S%LATT_CUR%A, &
          .FALSE., S%T_INFO%POSION, GROUP="/results", SUBGROUP="supercell")
    END IF
#endif

    ! clean up
    DEALLOCATE(BSEVECTOR,CHG,CW,CWORK)

    CALL DEALLOCATE_SUPERCELL(S)
    CALL DELWAV(EWF,.TRUE.)


    DO N=1,VBANDS
       CALL DELWAV(W1(N), .TRUE.)
    ENDDO

    DO N=1,CBANDS
       CALL DELWAV(W3(N), .TRUE.)
    ENDDO
    CALL BLACS_GRIDINFO(DESC_V(CTXT_), NPROW, NPCOL, MYROW, MYCOL)
    IF(MYROW /= -1 .OR. MYCOL /= -1) CALL BLACS_GRIDEXIT(ICTXC)

    PROFILING_STOP('exciton_wf')

    CONTAINS

    SUBROUTINE DUMP_CHARGE_RGRID( FILEN, GRID, CHG, IU0)
       CHARACTER(LEN=*) :: FILEN
       TYPE (grid_3d)   :: GRID
       RGRID            :: CHG(:)
       !local
       INTEGER          :: ISP, IU0
       IF( IU0>=0 ) THEN
          OPEN(UNIT=99,FILE=FILEN,STATUS='UNKNOWN')
          CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,S%T_INFO,LATT_CUR%SCALE,S%LATT_CUR%A, &
         &             .FALSE., S%T_INFO%POSION)
       ENDIF
       CALL OUTCHG_RGRID(GRID,99,.FALSE.,CHG(:))
       IF( IU0>=0 ) CLOSE(99)
    END SUBROUTINE DUMP_CHARGE_RGRID

  END SUBROUTINE EXCITON_WF

  SUBROUTINE CHECK_BSEMATRIX(BSEMATRIX,NDIM,IFAIL)
    IMPLICIT NONE
#ifdef single_prec_bse
    GDEFS, TARGET :: BSEMATRIX(:)
#endif
#ifdef double_prec_bse
    GDEF, TARGET :: BSEMATRIX(:)
#endif
    INTEGER, INTENT(IN) :: NDIM
    INTEGER, INTENT(INOUT) :: IFAIL
    INTEGER :: I

    !$ACC ENTER DATA COPYIN(IFAIL) __IF_ASYNC__
    !$ACC PARALLEL LOOP PRESENT(BSEMATRIX) REDUCTION(+:IFAIL) __IF_ASYNC__
    DO I=1,SIZE(BSEMATRIX)
      IF(BSEMATRIX(I) == MAGIC) IFAIL=IFAIL+1
    ENDDO
    !$ACC EXIT DATA COPYOUT(IFAIL) __IF_ASYNC__
    !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)


  END SUBROUTINE CHECK_BSEMATRIX

!**********************************************************************
!
!> Writes out BSE eigenfunctions and eigenvalues
!> 			M.Bokdam Nov 2014
!**********************************************************************

  SUBROUTINE WRITE_BSE(WHF,NCV,BSE_INDEX,EIGVECT,R,NBSEEIG)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    TYPE (latt) LATT_CUR
    TYPE(bse_matrix_index) :: BSE_INDEX
    INTEGER :: NCV
    INTEGER :: NBSEEIG
    COMPLEX(q) :: EIGVECT(:)
    REAL(q) :: R(:)
 ! local
    INTEGER K1,K2,K3,K4,NPOS1,NPOS2,MINL
    INTEGER I, J
    INTEGER     :: LAMBDA, ISP
    REAL(q) :: RIP(NCV),MINV 
#ifdef VASP_HDF5
    COMPLEX(q) :: FATBANDS(NCV,NBSEEIG)
#endif

      K1=1
      DO NPOS2=1,SIZE(R)
        K2=BSE_INDEX%NK(NPOS2)
        K3=BSE_INDEX%N1(NPOS2)
        K4=BSE_INDEX%N3(NPOS2)
        ISP=BSE_INDEX%ISP(NPOS2)
        RIP(NPOS2)=REAL(WHF%CELTOT(K4,K2,ISP)-WHF%CELTOT(K3,K2,ISP))
        K1=K1+1
      ENDDO
!MB very inefficient SORTING routine, but ok for now
      DO NPOS1=1,NBSEEIG
        MINV=MINVAL(RIP(NPOS1:NCV)) !locate
        MINL=INT(MINLOC(RIP(NPOS1:NCV),1))+NPOS1-1
        RIP(MINL)=RIP(NPOS1)        !swap
        RIP(NPOS1)=MINV
      ENDDO

      OPEN(unit=78, file='BSEFATBAND', status="UNKNOWN")
      WRITE(78,'(2I18)') SIZE(R),NBSEEIG
      K1=1
      DO NPOS1=1,NBSEEIG
        WRITE(78,'(I6,A14,F14.8,A20,F14.8)') NPOS1,"BSE eigenvalue:",R(NPOS1), "IP-eigenvalue:",RIP(NPOS1)
        DO NPOS2=1,SIZE(R)
          K2=BSE_INDEX%NK(NPOS2)
          K3=BSE_INDEX%N1(NPOS2)
          K4=BSE_INDEX%N3(NPOS2)
          ISP=BSE_INDEX%ISP(NPOS2)
          WRITE(78,'(3F9.5,3F14.7,2I6,F14.6,A5,F14.6)') WHF%WDES%VKPT(:,K2),REAL(WHF%CELTOT(K3,K2,ISP)), &
          & REAL(WHF%CELTOT(K4,K2,ISP)),ABS(EIGVECT(K1))/WHF%WDES%WTKPT(K2),K3,K4,REAL(EIGVECT(K1)), &
          & " +i* ",AIMAG(EIGVECT(K1))

#ifdef VASP_HDF5
          FATBANDS(NPOS2,NPOS1) = EIGVECT(K1) / WHF%WDES%WTKPT(K2)
#endif
          K1=K1+1
        ENDDO
      ENDDO
      CLOSE(78)

#ifdef VASP_HDF5
      CALL VH5_WRITE_BSEFATBAND(IH5OUTFILEID, FATBANDS)
#endif

  END SUBROUTINE WRITE_BSE


!**********************************************************************
!
!> Writes out BSE eigenfunctions and eigenvalues
!> It is a scaLAPACK-aware version
!> based on #CALCULATE_BSE_OSZI_SCALA
!> 			M.Bokdam April  2016, Aug 2018
!**********************************************************************

  SUBROUTINE WRITE_BSE_SCALA(WHF, ISP_LOW, ISP_HIGH,  &
       BD, BSE_INDEX, R, AMAT, NCV, NBSEEIG, IU6, DESC)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP_LOW, ISP_HIGH
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: R(:)
#ifdef single_prec_bse
    GDEFS :: AMAT(:)
#else
    GDEF  :: AMAT(:)
#endif
    INTEGER :: NCV
    INTEGER :: ANTIRES
    INTEGER :: IU6
 ! local
    INTEGER :: NCV13
    INTEGER K1, ISP
    INTEGER N1, N3
    INTEGER K2,K3,K4,NPOS1,NPOS2 !MB
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    INTEGER I, J, IERR
    INTEGER     :: LAMBDA
    INTEGER, OPTIONAL :: DESC(DLEN_)
!MB
    COMPLEX(q), ALLOCATABLE :: EIGVECT(:)
    INTEGER :: NBSEEIG

! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP,NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2

!Large NBSEEIG value can overrun the memory of a local node.
    ALLOCATE(EIGVECT(NBSEEIG*NCV))
    EIGVECT=0.0_q

    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    IF(MYCOL /= -1 .OR. MYROW /= -1) THEN
       ! loop over column indices (eigenvectors)
       DO LAMBDA=1, NBSEEIG
          ! local storage index in block (J1) and local index of block (J2)
          J1=   ((LAMBDA-DESC(NB_)*MYCOL-1)/(NPCOL*DESC(NB_)))*DESC(NB_)+1
          J2=MOD((LAMBDA-DESC(NB_)*MYCOL-1) ,NPCOL*DESC(NB_))+1

          JCOL=J1+J2-1
          IF (J2 <= 0 .OR. J2 > DESC(NB_)) THEN
             JCOL  =-1
          ELSE
             IF (LAMBDA /= DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2) THEN
                CALL vtutor%bug("internal error in WRITE_BSE_SCALA: " // str(LAMBDA) // " " // &
                   str(DESC(NB_)) // " " // str(NPCOL*(J1-1)) // " " // str(J2), __FILE__, __LINE__)
             ENDIF
          ENDIF

          IF (JCOL==-1) THEN
             ! no data on local node for eigenvector/eigenvalue pair LAMBDA_DONE
          ELSE
             IROW=0
             DO I1=1,NP,DESC(MB_)
                I1RES=MIN(DESC(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1

                   NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                   ISP= BSE_INDEX%ISP(NCV13)
                   K1 = BSE_INDEX%NK(NCV13)

                   IF (WHF%WDES%WTKPT(K1)==0) CYCLE

                   NB1_INTO_TOT=BSE_INDEX%N1(NCV13)
                   N1=NB1_INTO_TOT-BD(ISP)%VBMIN+1

                   NB3_INTO_TOT=BSE_INDEX%N3(NCV13)
                   N3=NB3_INTO_TOT-BD(ISP)%CBMIN+1

                   IF (NCV13/=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)) THEN
                      CALL vtutor%bug("internal error in WRITE_BSE_SCALA: " // str(NCV13) // " " // &
                         str(NB1_INTO_TOT-BD(ISP)%VBMIN+1) // " " // str(NB3_INTO_TOT) // " " // str(K1) // " "&
                          // str(ISP) // " " // str(BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1,&
                         NB3_INTO_TOT-BD(ISP)%CBMIN+1,K1,ISP -ISP_LOW+1)), __FILE__, __LINE__)
                   ENDIF

                   IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K1,ISP))) CYCLE
                   !IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)))  CYCLE

                   EIGVECT((LAMBDA-1)*NCV+NCV13)=AMAT(IROW+(JCOL-1)*DESC(LLD_))
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDIF

!Gather all parts of the eigenvectors stored on the different nodes
    CALLMPI( M_sum_z(WHF%WDES%COMM, EIGVECT, NBSEEIG*NCV))
    IF (IU6>=0) CALL WRITE_BSE(WHF,NCV,BSE_INDEX,EIGVECT,R,NBSEEIG)
      
    DEALLOCATE(EIGVECT)

  END SUBROUTINE WRITE_BSE_SCALA

!****************** subroutine hermiticity    *********************************
!
!> enforces hermiticity of a matrix by removing the non-zero imaginary part in the diagonal elements
!
!******************************************************************************
  subroutine hermicity(n, a, desca ) 

      implicit none
      integer,parameter :: ctxt_=2 !< the blocking factor used to distribute the rows of the array.
      integer,parameter :: mb_=5   !< the blocking factor used to distribute the columns of the array.
      integer,parameter :: nb_=6   !< the leading dimension of the local array.  lld_a >= max(1,locr(m_a)).
      integer,parameter :: lld_=9  !< dimension of desca
      integer,parameter :: dlen_=9

      integer desca( dlen_ ) !< distributed matrix descriptor array
      complex(mq) a(:)       !< distributed matrix
      integer n              !< dimension of the matrix
      ! local
      integer np,nq
      integer i1res,j1res,jcol
      integer myrow, mycol, nprow, npcol
      integer i1,i2,j1,j2

      intrinsic  min
      integer    numroc,irow_jcol
      external   numroc,blacs_gridinfo


      call blacs_gridinfo(desca(ctxt_),nprow,npcol,myrow,mycol)
      np = numroc(n,desca(mb_),myrow,0,nprow)
      nq = numroc(n,desca(nb_),mycol,0,npcol)

      jcol = 0
      ! loop over local columns
      do j1=1,nq,desca(nb_)
        j1res=min(desca(nb_),nq-j1+1)
        do j2=1,j1res
          irow_jcol = jcol

          ! loop over local rows
          do i1=1,np,desca(mb_)
            i1res=min(desca(mb_),np-i1+1)
            do i2=1,i1res
              irow_jcol = irow_jcol+1
              if(  (desca(mb_)*myrow+nprow*(i1-1)+i2) == &
                   (desca(nb_)*mycol+npcol*(j1-1)+j2) ) then
                a(irow_jcol) = real(a(irow_jcol),mq)
                goto 100 ! next column, since we have found the right row
              endif
            enddo
          enddo
100       jcol = jcol + desca(lld_)
        enddo
      enddo

      return

  end subroutine hermicity


!****************** subroutine matmul_kahan    ******************************
!
!> performs matrix-matrix multiplication, C=AB. Uses the Kahan method to keep track of errors due to numerical precision limits
!> has both a distribute and non-distribute option, though I do not know why one would choose the non-distributed version.
!
!*************************************************************************************
  subroutine matmul_kahan(comm, matrix1, matrix2, matrix3, nrows_1, ncols_1, nrows_2, ncols_2, nrows_3, ncols_3, distribute, matdim_rank)

    implicit none
    type(communic) :: comm      !< communicator to be used (set by vasp).
    integer     :: nrows_1      !< number of rows of matrix1, A 
    integer     :: ncols_1      !< number of columns of matrix1, A
    integer     :: nrows_2      !< number of rows of matrix2, B
    integer     :: ncols_2      !< number of columns of matrix2, B
    integer     :: nrows_3      !< number of rows of matrix3, C
    integer     :: ncols_3      !< number of columns of matrix3, C
    character(1) :: distribute                                        !< split over MPI tasks or not
    integer     :: matdim_rank(:)
#ifdef gammareal
    real(mq) :: matrix1(:,:)    !< matrix A, as defined in the header
    real(mq) :: matrix2(:,:)    !< matrix B, as defined in the header
    real(mq) :: matrix3(:,:)    !< matrix C, as defined in the header
    real(mq)              :: tmp(1:nrows_1,1:ncols_3)                    ! temporary array
#else
    complex(mq) :: matrix1(:,:) !< matrix A, as defined in the header
    complex(mq) :: matrix2(:,:) !< matrix B, as defined in the header
    complex(mq) :: matrix3(:,:) !< matrix C, as defined in the header
    complex(mq)              :: tmp(1:nrows_1,1:ncols_3)                 ! temporary array
    ! Local variables for Kahan algorithm
    real(mq) :: sum_i(1:nrows_1,1:ncols_3), c_i(1:nrows_1,1:ncols_3)  ! Sum and compensation for the Imaginary part
    real(mq) :: y_i(1:nrows_1,1:ncols_3), t_i(1:nrows_1,1:ncols_3)    ! Temporary variables
#endif
    ! Local variables for Kahan algorithm
    real(mq) :: sum_r(1:nrows_1,1:ncols_3), c_r(1:nrows_1,1:ncols_3)  ! Sum and compensation for the Real part
    real(mq) :: y_r(1:nrows_1,1:ncols_3), t_r(1:nrows_1,1:ncols_3)    ! Temporary variables
    integer  :: i, n, m, start, finish                                ! loop indices

    ! check that dimensions are compatible in C_nm = A_nk*B_km operation, sum over k
    if (nrows_1 /= nrows_3) call vtutor%error("Error in Kahan routine for matrix-matrix multiplication, C=AB. Number of rows of C &
and A differ. You should not be encontering this error, please contact the developers' team.")
    if (nrows_2 /= ncols_1) call vtutor%error("Error in Kahan routine for matrix-matrix multiplication, C=AB. Number of columns of A &
and number of rows of B differ. You should not be encontering this error, please contact the developers' team.")
    if (ncols_2 /= ncols_3) call vtutor%error("Error in Kahan routine for matrix-matrix multiplication, C=AB. Number of columns of C &
and B differ. You should not be encontering this error, please contact the developers' team.")

    ! Initialization
    sum_r = 0.0_mq
    c_r   = 0.0_mq
#ifndef gammareal
    sum_i = 0.0_mq
    c_i   = 0.0_mq
#endif
    ! safe zeroing of matrix3
    matrix3 = zerobse

    start = 0; finish = 0

    if (distribute == 't') then
      !fancier method that uses the information provided by the descriptor vector to know how many elements 
      !are store in each row 
      !note that the size of the array that's stored in each rank is saved in matdim_rank(comm%ncpu+1)
      if (matdim_rank(comm%node_me+1)>0) then
        start  = sum(matdim_rank(1:comm%node_me)) + 1
        finish = sum(matdim_rank(1:comm%node_me + 1))
      endif
      !dumb method that just divides the number of elements in the row by the number of ranks
!     start  = int((comm%node_me-1)*nrows_2/comm%ncpu) + 1
!     finish = int(comm%node_me*nrows_2/comm%ncpu)
    else
      start = 1
      finish = nrows_2
    endif
    ! Kahan Summation Loop
    if (matdim_rank(comm%node_me+1)>0) then
      do n = 1, ncols_3
        do i = start, finish
          tmp(1:nrows_1,n) = matrix1(1:nrows_1,i)*matrix2(i,n)
          ! Real Part Compensation
          ! y = input_real - compensation_real
          y_r(1:nrows_1,n) = real(tmp(1:nrows_1,n), mq) - c_r(1:nrows_1,n)

          ! t = sum_real + y
          t_r(1:nrows_1,n) = sum_r(1:nrows_1,n) + y_r(1:nrows_1,n)

          ! c = (t - sum_real) - y
          ! This calculates the lost low-order bits from the previous Sum_R + Y_R
          c_r(1:nrows_1,n) = (t_r(1:nrows_1,n) - sum_r(1:nrows_1,n)) - y_r(1:nrows_1,n)

          ! sum_real = t
          sum_r(1:nrows_1,n) = t_r(1:nrows_1,n)
#ifndef gammareal
          ! Imaginary Part Compensation ---
          ! y = input_imaginary - compensation_imaginary
          y_i(1:nrows_1,n) = aimag(tmp(1:nrows_1,n)) - c_i(1:nrows_1,n)

          ! t = sum_imaginary + y
          t_i(1:nrows_1,n) = sum_i(1:nrows_1,n) + y_i(1:nrows_1,n)

          ! c = (t - sum_imaginary) - y
          c_i(1:nrows_1,n) = (t_i(1:nrows_1,n) - sum_i(1:nrows_1,n)) - y_i(1:nrows_1,n)

          ! sum_imaginary = t
          sum_i(1:nrows_1,n) = t_i(1:nrows_1,n)  
#endif
        enddo
      enddo
      ! final result construction
    endif
#ifdef gammareal
    matrix3 = sum_r
#else
    matrix3 = cmplx(sum_r, sum_i, mq)
#endif

    if (distribute == 't')  call M_sumbse_v(comm,matrix3,36)

  end subroutine matmul_kahan

#endif
END MODULE bse_base
