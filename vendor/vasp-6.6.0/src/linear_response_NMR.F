#include "symbol.inc"

#define use_davidson
! #define apply_greens_function
! #define provide_w1
! #define exact_diagonalisation

      MODULE mlr_main_nmr
      USE prec
      USE vaspxml
      IMPLICIT NONE

      PRIVATE

      INTEGER, PARAMETER :: NSTEPS_DAVIDSON=15
      INTEGER :: NKPTS_ORIG

      COMPLEX(q) :: QQ(3,3,2)
      REAL(q) :: CHI_MAG(3,3)

      !< TOTAL CHEMICAL SHIFT/ABSOLUTE CHEMICAL SHIFT
      COMPLEX(q), ALLOCATABLE :: ABS_CHEM_SHIFT(:,:,:)       !< total chemical shift
      !< CONTRIBUTIONS
      COMPLEX(q), ALLOCATABLE :: B_OUT_PARA_PS(:,:,:)        !> paramagnetic planewave field
      REAL(q), ALLOCATABLE :: B_PARA_ONE_CENTR(:,:,:)        !> paramagnetic one centre field (star)
      REAL(q), ALLOCATABLE :: B_DIA_ONE_CENTR(:,:,:)         !> diamagnetic one centre field (center)
!< SOC_
      REAL(q), ALLOCATABLE :: B_SPIN_PS(:,:,:)               !> spin planewave field                         SOC
      REAL(q), ALLOCATABLE :: B_SPIN_ONE_CENTR_FC(:,:,:)     !> spin one centre fermi fontact field (star)   SOC
      REAL(q), ALLOCATABLE :: B_SPIN_ONE_CENTR_DIP(:,:,:)    !> spin one centre dipolar field (star)         SOC
      REAL(q), ALLOCATABLE :: B_SPIN_ONE_CENTR_LQ_FC(:,:,:)  !> spin one centre fermi fontact field (center) SOC
      REAL(q), ALLOCATABLE :: B_SPIN_ONE_CENTR_LQ_DIP(:,:,:) !> spin one centre dipolar field (center)       SOC
!> SOC_

!< NICS_
      REAL(q), ALLOCATABLE :: POSROT(:,:)
      INTEGER :: NPOSROT=1
!> NICS_

      PUBLIC :: MLR_SKELETON, EIJK, NMR_MAGRES_TENS

      CONTAINS

!> ************************* SUBROUTINE MLR_SKELETON **************************
!>
!> Wrapper routine for calculating chemical shielding using GIPAW.
!> VASP Literature:
!> J. Chem. Phys. 146, 064115 (2017)
!> LLRAUG - J. Chem. Phys. 155, 234101 (2021)
!> SOC extension (+ review) - arXiv:2504.17380 [cond-mat.mtrl-sci]
!>
!> Background Literature:
!> GIPAW - Pickard and Mauri (PM), Phys. Rev. B, 63 245101 (2001)
!> Ultrasoft - Yates, Pickard, and Mauri (YPM), Phys. Rev. B 76, 024401 (2007)
!> ZORA (scalar-relativistic) - J. Chem. Phys. 118, 5746 (2003)
!>
!> Wrapper routine:
!> - Reduces kpoint symmetry if necessary (LNMR_SYM_RED)
!> - Calls appropriate subroutines
!>   (mainly MLR_B_MAIN within b-loop,
!>    i.e. loop over magnetic field directions)
!> - Writes output (CSA Tensor)
!>
!> ****************************************************************************
      SUBROUTINE MLR_SKELETON( &
     &   HAMILTONIAN,W,WDES,GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,KPOINTS,LATT_CUR,LATT_INI, &
     &   T_INFO,DYN,SYMM,P,XC,NONL_S,NONLR_S,LMDIM,CDIJ,CQIJ,SV,E,INFO,IRDMAX,IO)
      USE prec
      USE base
      USE pseudo
      USE poscar
      USE wave_high
      USE nonl_high
      USE hamil_struct_def
      USE mgrid
      USE lattice
      USE constant
      USE kpoints_change
      USE msymmetry
      USE david
      USE subrot
      USE paw
      USE pead, ONLY : LUSEPEAD,PEAD_RESETUP_WDES,PEAD_DPSI_DK_ALL,LPEAD_SYM_RED
      USE morbitalmag
      USE choleski
      USE setexm_struct_def, ONLY : xc_info
      USE fileio,            ONLY: WRITE_GRID_QUANTITY
#ifdef VASP_HDF5
      USE vhdf5
#endif
      IMPLICIT NONE
      TYPE (ham_handle) HAMILTONIAN
      TYPE (wavespin) W
      TYPE (wavedes) WDES
      TYPE (grid_3d) GRID
      TYPE (grid_3d) GRID_SOFT               ! soft grid for pseudized potentials/ charge etc.
      TYPE (grid_3d) GRIDC                   ! full, fine grid
      TYPE (grid_3d) GRIDUS
      TYPE (transit) C_TO_US
      TYPE (transit) SOFT_TO_C
      TYPE (kpoints_struct) KPOINTS
      TYPE (latt) LATT_CUR
      TYPE (latt) LATT_INI
      TYPE (type_info) T_INFO
      TYPE (dynamics) DYN
      TYPE (symmetry) SYMM
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct)NONLR_S
      TYPE (energy) E
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (potcar), TARGET :: P(:)
      TYPE (xc_info) XC
      INTEGER LMDIM,IRDMAX
      OVERLAP CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      OVERLAP CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)

      !< local variables
      INTEGER I,J,NI
      INTEGER BDIR !< cartesian component of B-field, 1=x, 2=y, 3=z

      COMPLEX(qs), ALLOCATABLE :: RPHI(:,:,:,:,:)
      GDEFS, ALLOCATABLE ::  RPHI_CPROJ(:,:,:,:,:)
      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
           GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      REAL(q) TMP(3,3), TMP_TENS_AT(3,3), TMP_MACRO(3,3)
      REAL(q) TMP_TENS(3,3,T_INFO%NIOND)
      REAL(q) ISO1,SPAN1,SKEW1
      REAL(q) ISO2,SPAN2,SKEW2
      REAL(q), ALLOCATABLE :: CORE_SHIFT(:)
      REAL(q) CORE_CHI, CHI_FACT
      INTEGER K1,K2
      INTEGER ISUSP
      CHARACTER*5 MODUS

      TYPE (skpoints_trans) :: KPOINTS_TRANS
      INTEGER NKORIG
#ifndef use_davidson
      REAL(q) DESUM
#else
      INTEGER NELM
      INTEGER ICOUEV,NSIM
      REAL(q) RMS,DESUM1,TOTEN
#endif

!< NICS_
      INTEGER NPOSROTMAX,NSPTMP,NSXTMP
      INTEGER, DIMENSION(1):: NANICS
      INTEGER, ALLOCATABLE :: MAPNICS(:,:,:,:)
      REAL(q), ALLOCATABLE :: TRNICS(:,:),TAUWRK(:,:)
      COMPLEX(q), ALLOCATABLE :: NICSDATA(:,:,:)
      REAL(q), ALLOCATABLE :: NICSDATA_UNSYM(:,:,:)
      REAL(q), ALLOCATABLE :: NICSDATA_SYM(:,:,:)
!> NICS_
      REAL(q) CHI_MAG_DIAG(3,3),EIG_CHI(3)
      INTEGER :: IDIR1,IDIR2
      INTEGER, PARAMETER :: LWORK=6
      REAL(q) WORK(3*LWORK)
      INTEGER IFAIL

      !< ALLOCATE ARRAYS FOR RESULTS
      ALLOCATE(ABS_CHEM_SHIFT(3,3,T_INFO%NIONS))

      !< contribs (note awkard differences in dimensions, thank MPI)
      ALLOCATE(B_OUT_PARA_PS(3,T_INFO%NIONS,3))
      ALLOCATE(B_PARA_ONE_CENTR(3,T_INFO%NIONS,3))
      ALLOCATE(B_DIA_ONE_CENTR(3,T_INFO%NIONS,3))
      IF (LSOSHIFT) THEN
         ALLOCATE(B_SPIN_PS(3,T_INFO%NIONS,3))
         ALLOCATE(B_SPIN_ONE_CENTR_FC(3,T_INFO%NIONS,3))
         ALLOCATE(B_SPIN_ONE_CENTR_DIP(3,T_INFO%NIONS,3))
         ALLOCATE(B_SPIN_ONE_CENTR_LQ_FC(3,T_INFO%NIONS,3))
         ALLOCATE(B_SPIN_ONE_CENTR_LQ_DIP(3,T_INFO%NIONS,3))
      ENDIF

      !< initialized fields
      B_OUT_PARA_PS=0._q
      B_PARA_ONE_CENTR=0._q
      B_DIA_ONE_CENTR=0._q
      IF (LSOSHIFT) THEN
         B_SPIN_PS = 0._q
         B_SPIN_ONE_CENTR_FC = 0._q
         B_SPIN_ONE_CENTR_DIP = 0._q
         B_SPIN_ONE_CENTR_LQ_FC = 0._q
         B_SPIN_ONE_CENTR_LQ_DIP = 0._q
      ENDIF

#ifdef use_davidson
      NSIM=W%WDES%NSIM*2
#ifdef MPI
      NSIM=((W%WDES%NSIM*2+W%WDES%COMM_INTER%NCPU-1)/W%WDES%COMM_INTER%NCPU)*W%WDES%COMM_INTER%NCPU
#endif
#endif

      !< setup one-centre partial wave matrix elements for current augmentation
      IF (LLRAUG) CALL CALC_JPAW(T_INFO, P)

      QQ = 0
      CHI_MAG = 0

      !> possibly we want to reduce the symmetry to strictly
      !> keep the q-stars mapped on themselves for all symmetry operations
      IF (LNMR_SYM_RED.AND.SYMM%ISYM>0) THEN
         NKORIG=W%WDES%NKPTS

         CALL MLR_SYM_REDUCE(LATT_CUR,T_INFO,SYMM,IO%IU6)
#ifdef oldsym
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
              SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
              T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#else
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
              SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
              T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#endif
         CALL KPAR_SYNC_ALL(WDES,W)
         CALL RE_GEN_LAYOUT(GRID,WDES,KPOINTS,LATT_CUR,LATT_INI,IO%IU6,IO%IU0)
         IF (LUSEPEAD()) CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)
         CALL REALLOCATE_WAVE(W,GRID,WDES,NONL_S,T_INFO,P,LATT_CUR,KPOINTS_TRANS)

         CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,W)
#ifndef use_davidson
         !< Loewdin perturbation theory to improve states at added k-points
         CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
         &     LMDIM,CDIJ,CQIJ,4,SV,T_INFO,P,IO%IU0,DESUM,NKSTART=NKORIG+1)
         CALL KPAR_SYNC_ALL(WDES,W)
#else
         !< optimize the states at added k-points with Davidson optimization
         TOTEN=0; INFO%IALGO=8
         DO NELM=1,NSTEPS_DAVIDSON
            CALL EDDAV(HAMILTONIAN,P,GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,NSIM, &
           &     LMDIM,CDIJ,CQIJ,RMS,DESUM1,ICOUEV,SV,E%EXHF,IO%IU6,IO%IU0, &
           &     LDELAY=.FALSE.,LSUBROTI=.TRUE.,LEMPTY=.FALSE.,LHF=.TRUE., &
           &     NKSTART=NKORIG+1)

            E%EBANDSTR=BANDSTRUCTURE_ENERGY(W%WDES,W)

            IF (IO%IU0>=0) WRITE(17, 200)      NELM,E%EBANDSTR,E%EBANDSTR-TOTEN,DESUM1,ICOUEV,RMS
            IF (IO%IU0>=0) WRITE(IO%IU0, 200)  NELM,E%EBANDSTR,E%EBANDSTR-TOTEN,DESUM1,ICOUEV,RMS
            TOTEN=E%EBANDSTR

 200        FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
                 &       I6,'  ',E10.3)

            CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)

            IF (ABS(DESUM1) < ABS(INFO%EDIFF) .AND. NELM>1) EXIT
         ENDDO
         CALL KPAR_SYNC_ALL(W%WDES,W)
#endif
      ENDIF

!< NICS_
      IF(NUCIND.AND..NOT.LNICSALL) THEN
        NPOSROTMAX=NROTK*NPCELL*NPOS_NICS
        ALLOCATE(TAUWRK(NPOSROTMAX,3))
        ! set up all symmetry equivalent positions of sites read from POSNICS,
        ! non-primitive translations are included, but not needed.
        CALL PSYMMAKE(TAUNICS,NPOS_NICS,ISYMOP,GTRANS,SYMM%PTRANS,NROTK,NPCELL,TAUWRK,NPOSROT,NPOSROTMAX)
        ! POSROT starts with NPOS_NICS sites from POSNICS and has all symmetry equivalent positions. Total: NPOSROT
        ALLOCATE(POSROT(NPOSROT,3))
        DO K1=1,NPOSROT
           POSROT(K1,:)=TAUWRK(K1,:)
        ENDDO
        DEALLOCATE(TAUWRK)
        NANICS(1)=NPOSROT
        ALLOCATE(MAPNICS(NPOSROT,1,NROTK,NPCELL))
        ALLOCATE(TRNICS(NPOSROT,3))
        NSPTMP=1
        NSXTMP=1
        CALL POSMAP(POSROT,ISYMOP,GTRANS,NROTK,SYMM%PTRANS,NPCELL,T_INFO%NIOND,NSXTMP,NSPTMP,NPOSROT,NANICS,MAPNICS,TRNICS)
        DEALLOCATE(TRNICS)
        ALLOCATE(NICSDATA(NPOSROT,3,3))
        ALLOCATE(NICSDATA_UNSYM(3,3,NPOSROT))
      ELSEIF (NUCIND.AND.LNICSALL) THEN
        ALLOCATE(NICSDATA(GRIDC%MPLWV,3,3))
      ELSE
        ALLOCATE(NICSDATA(1,3,3))
      ENDIF
!> NICS_

      !< Primary loop for calculating chemical shielding tensor.
      bloop: DO BDIR=1,3 !< loop over magnetic field directions

         IF (LUSEPEAD()) THEN
            ALLOCATE(RPHI(W%WDES%NRPLWV,W%WDES%NBANDS,W%WDES%NKPTS,W%WDES%ISPIN,3), &
           &   RPHI_CPROJ(W%WDES%NPROD,W%WDES%NBANDS,W%WDES%NKPTS,W%WDES%ISPIN,3))
            RPHI=0; RPHI_CPROJ=0
            CALL PEAD_DPSI_DK_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO,RPHI,RPHI_CPROJ)
            !< calculate induced field
            CALL MLR_B_MAIN( &
           &   HAMILTONIAN,W,GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,KPOINTS,LATT_CUR,LATT_INI,T_INFO,SYMM,P, &
           &   NONL_S,NONLR_S,LMDIM,CDIJ,CQIJ,SV,E,INFO,IRDMAX,IO,BDIR,XC,NICSDATA(:,:,BDIR),RPHI)
            DEALLOCATE(RPHI,RPHI_CPROJ)
         ELSE
            CALL MLR_B_MAIN( &
           &   HAMILTONIAN,W,GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,KPOINTS,LATT_CUR,LATT_INI,T_INFO,SYMM,P, &
           &   NONL_S,NONLR_S,LMDIM,CDIJ,CQIJ,SV,E,INFO,IRDMAX,IO,BDIR,XC,NICSDATA(:,:,BDIR))
         ENDIF
      ENDDO bloop

!< NICS_
      IF (NUCIND.AND.LNICSALL) THEN
         ! symmetrization of NICS in rec space
         IF (SYMM%ISYM>0) THEN
            CALL SYMTENS(NICSDATA,GRIDC,SYMM%PTRANS,T_INFO%NIOND,LATT_CUR,[0._q,0._q,1._q],.FALSE.)
         ENDIF
         CALL WRITE_GRID_QUANTITY("NICS",W%WDES,T_INFO%POSION,T_INFO,LATT_CUR,GRIDC,&
            &    RESHAPE(NICSDATA,[SIZE(NICSDATA,1),SIZE(NICSDATA,2)*SIZE(NICSDATA,3)]),&
            &    LFFT=.TRUE.)
      ENDIF
!> NICS_

      CALLMPI( M_sum_z(W%WDES%COMM,QQ,18))

      DO ISUSP=1,2   !< loop over susceptibilty formulas: 1=pGv, 2=vGv
      DO I=1,3
      DO J=1,3
         QQ(I,J,ISUSP) = QQ(I,J,ISUSP) /DQ/DQ/TPI/TPI/LATT_CUR%OMEGA
         IF (I.NE.J) QQ(I,J,ISUSP) = 2.*QQ(I,J,ISUSP)
         QQ(I,J,ISUSP) = QQ(I,J,ISUSP) * 8._q * PI/3._q
         QQ(I,J,ISUSP) = QQ(I,J,ISUSP) * 1E6 *MOMTOMOM
      ENDDO
      ENDDO
      ENDDO

#ifdef VASP_HDF5
      ! meta-data for WRT_NMRCUR
      IF (WRT_NMRCUR>0) CALL VH5_WRITE_NMRCUR_META(IH5OUTFILEID,WRT_NMRCUR)
#endif

      IF (IO%IU6>=0) THEN
         IF (LNMRLEG) THEN !< Legacy output
         WRITE (IO%IU6,'(/A)') ' -------------------------------------------------------------'
         WRITE (IO%IU6,*)      ' Absolute Chemical Shift tensors'
         WRITE (IO%IU6,*)      '-------------------------------------------------------------'
         TMP_TENS=0._q
         DO NI=1,T_INFO%NIONS

            TMP_TENS(1:3,1:3,NI)= SHIFTSHIELD * REAL(ABS_CHEM_SHIFT(1:3,1:3,NI),q)

            IF (NI.GT.T_INFO%NIOND) THEN
               CALL vtutor%error("ERROR, STOP")
            ENDIF
         ENDDO
         WRITE (IO%IU6,*) ' UNSYMMETRIZED TENSORS '
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(" ion ",I4)') NI
            DO BDIR=1,3
               WRITE(IO%IU6,'(3F18.6)') TMP_TENS(:,BDIR,NI)
            ENDDO
         ENDDO
!        CALL TSYM(TMP,ISYMOP,NROTK,LATT_CUR%A)
         CALL TENSYM(TMP_TENS, &
      &     SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,LATT_CUR%A)
         WRITE (IO%IU6,*) ' SYMMETRIZED TENSORS '
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(" ion ",I4)') NI
            DO BDIR=1,3
               WRITE(IO%IU6,'(3F18.6)') TMP_TENS(:,BDIR,NI)
            ENDDO
         ENDDO
         WRITE (IO%IU6,'(A/)') ' -------------------------------------------------------------'
         ELSE !< New output
            WRITE (IO%IU6,'(/A)') ' --------------------------------------------------------------'

            IF (LNMRSHIELD) THEN
               WRITE (IO%IU6,*)      ' Absolute Chemical Shielding Tensors'
            ELSE
               WRITE (IO%IU6,*)      ' Absolute Chemical Shift Tensors (= Negative Chemical Shielding)'
            ENDIF

            WRITE (IO%IU6,*)      ' VALENCE ONLY (Space Group Symmetrized)'
            WRITE (IO%IU6,*)
            WRITE (IO%IU6,*)      ' BDIR labels direction of applied magnetic field (i.e. B0)'
            WRITE (IO%IU6,*)      ' For BDIR==1, B0 along x-axis, etc.'
            WRITE (IO%IU6,*)      ' Induced field listed along cartesian directions for each BDIR'
            WRITE (IO%IU6,'(/A)') ' --------------------------------------------------------------'
            WRITE (IO%IU6,'(/A)') ' ORBITAL CONTRIBUTION (DIAMAGNETIC + PARAMAGNETIC)'
            TMP_TENS = 0._q
            DO NI=1, T_INFO%NIONS
               TMP_TENS(1:3,1:3,NI) = SHIFTSHIELD * (REAL(B_OUT_PARA_PS(1:3,NI,1:3),q)+B_PARA_ONE_CENTR(1:3,NI,1:3)+B_DIA_ONE_CENTR(1:3,NI,1:3))

               IF (NI.GT.T_INFO%NIOND) THEN
                  CALL vtutor%error("ERROR, STOP")
               ENDIF
            ENDDO

            !<Space group symmetrize
            CALL TENSYM(TMP_TENS, &
      &     SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,LATT_CUR%A)

            !< Write tensors per ion
            CALL NMR_WRITE_TENS(T_INFO%NIONS, T_INFO%NIOND, TMP_TENS, IO)

            IF (LSOSHIFT) THEN
               WRITE (IO%IU6,'(/A)') ' SPIN CONTRIBUTION'
               TMP_TENS = 0._q
               DO NI=1, T_INFO%NIONS
                  TMP_TENS(1:3,1:3,NI) = SHIFTSHIELD * (B_SPIN_PS(1:3,NI,1:3)+B_SPIN_ONE_CENTR_FC(1:3,NI,1:3)+B_SPIN_ONE_CENTR_DIP(1:3,NI,1:3) &
  &                                                           +B_SPIN_ONE_CENTR_LQ_FC(1:3,NI,1:3)+B_SPIN_ONE_CENTR_LQ_DIP(1:3,NI,1:3))

                  IF (NI.GT.T_INFO%NIOND) THEN
                     CALL vtutor%error("ERROR, STOP")
                  ENDIF
               ENDDO

               !< Space group symmetrize
               CALL TENSYM(TMP_TENS, &
         &     SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,LATT_CUR%A)

               !< Write tensors per ion
               CALL NMR_WRITE_TENS(T_INFO%NIONS, T_INFO%NIOND, TMP_TENS, IO)
            ENDIF

            WRITE (IO%IU6,'(/A)') ' TOTAL'
            TMP_TENS = 0._q
            DO NI=1, T_INFO%NIONS
               TMP_TENS(1:3,1:3,NI)= SHIFTSHIELD * REAL(ABS_CHEM_SHIFT(1:3,1:3,NI),q)

               IF (NI.GT.T_INFO%NIOND) THEN
                  CALL vtutor%error("ERROR, STOP")
               ENDIF
            ENDDO

            !< Space group symmetrize
            CALL TENSYM(TMP_TENS, &
      &     SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,LATT_CUR%A)

            !< Write tensors per ion
            CALL NMR_WRITE_TENS(T_INFO%NIONS, T_INFO%NIOND, TMP_TENS, IO)

            WRITE (IO%IU6,'(A/)') ' --------------------------------------------------------------'
         ENDIF
!< NICS_
         IF (NUCIND.AND..NOT.LNICSALL) THEN ! POSNICS related output
            ! legacy output prints the unsymmetrized nics: create backup
            ! and reorganize data to call tensym
            DO IDIR1=1,3
              DO IDIR2=1,3
                DO I=1,NPOSROT
                  NICSDATA_UNSYM(IDIR1,IDIR2,I)=REAL(NICSDATA(I,IDIR1,IDIR2))
                ENDDO
              ENDDO
            ENDDO
            DEALLOCATE(NICSDATA)
            ! symmetrization
            ALLOCATE(NICSDATA_SYM, SOURCE=NICSDATA_UNSYM)
            CALL TENSYM(NICSDATA_SYM(1,1,1),MAPNICS,1,NANICS,NPOSROT,LATT_CUR%A)
            ! write OUTCAR and vaspout.h5
            CALL WRITE_NICS_FOR_POSNICS(IO,LATT_CUR,T_INFO,TAUNICS,SHIFTSHIELD*NICSDATA_SYM,NPOS_NICS,&
              & NICSDATA_UNSYM,NPOSROT,LNMRLEG,LNMRSHIELD)
         ENDIF
!> NICS_

         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) ' ORBITAL MAGNETIC SUSCEPTIBILITY, excluding core contribution'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO ISUSP=1,2
            TMP_MACRO(1:3,1:3)=REAL(QQ(1:3,1:3,ISUSP),q)
            CALL TSYM(TMP_MACRO,ISYMOP,NROTK,LATT_CUR%A)
            CHI_FACT= 3._q/8._q/PI*LATT_CUR%OMEGA*6.022142E23/1.E24 !< N_Avogadro * Ang_to_cm^3
            CHI_MAG(1:3,1:3)= TMP_MACRO(1:3,1:3) * CHI_FACT

            ! These if statements derive functionality from the ISUSP loop.
            IF     ((     LVGVAPPL).AND.(ISUSP.EQ.2)) THEN ! use vGv
               TMP(1:3,1:3)=TMP_MACRO(1:3,1:3)
            ELSEIF ((.NOT.LVGVAPPL).AND.(ISUSP.EQ.1)) THEN ! use pGv
               TMP(1:3,1:3)=TMP_MACRO(1:3,1:3)      
            ENDIF

            IF (ISUSP.EQ.1) WRITE (IO%IU6,*) ' Approximate magnetic susceptibility, pGv (10^-6 cm^3/mole) '
            IF (ISUSP.EQ.2) WRITE (IO%IU6,*) ' Approximate magnetic susceptibility, vGv (10^-6 cm^3/mole) '
            DO BDIR=1,3
               WRITE (IO%IU6,'(I6,6F18.6)') BDIR,CHI_MAG(1,BDIR),CHI_MAG(2,BDIR),CHI_MAG(3,BDIR)
            ENDDO
            IF (ISUSP.EQ.2) THEN
               !< squeeze out any remaining asymmetry
               DO IDIR1=1,3
               DO IDIR2=1,3
                 CHI_MAG_DIAG(IDIR1,IDIR2)=(CHI_MAG(IDIR1,IDIR2)+CHI_MAG(IDIR2,IDIR1))/2._q
               ENDDO
               ENDDO
               !< Get eigenvalues and eigenvectors
               CALL DSYEV('V','U',3,CHI_MAG_DIAG,3,EIG_CHI,WORK,3*LWORK,IFAIL)
               WRITE (IO%IU6,*)
               WRITE (IO%IU6,*) '        principal value                      axis'
               WRITE (IO%IU6,*) '      (10^-6 cm^3/mole)           x,          y,          z'
               WRITE (IO%IU6,*) '     --------------------------------------------------------'
               WRITE (IO%IU6,'(6X,F18.6,3F12.4)') EIG_CHI(1),CHI_MAG_DIAG(1,1),CHI_MAG_DIAG(2,1),CHI_MAG_DIAG(3,1)
               WRITE (IO%IU6,'(6X,F18.6,3F12.4)') EIG_CHI(2),CHI_MAG_DIAG(1,2),CHI_MAG_DIAG(2,2),CHI_MAG_DIAG(3,2)
               WRITE (IO%IU6,'(6X,F18.6,3F12.4)') EIG_CHI(3),CHI_MAG_DIAG(1,3),CHI_MAG_DIAG(2,3),CHI_MAG_DIAG(3,3)
            ENDIF
            WRITE (IO%IU6,*) '-------------------------------------------------------------'
            IF (.NOT.LVGVCALC) EXIT
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         IF (LNMRSHIELD) THEN
            WRITE (IO%IU6,*) ' G=0 CONTRIBUTION TO CHEMICAL SHIELDING (field along BDIR) '
            TMP(:,:) = SHIFTSHIELD * TMP(:,:)
         ELSE
            WRITE (IO%IU6,*) ' G=0 CONTRIBUTION TO CHEMICAL SHIFT (field along BDIR) '
         ENDIF

         IF (LVGVAPPL) THEN
            WRITE (IO%IU6,*) ' using vGv susceptibility, excluding core contribution'
         ELSE
            WRITE (IO%IU6,*) ' using pGv susceptibility, excluding core contribution'
         ENDIF
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  BDIR                X                 Y                 Z'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO BDIR=1,3
            WRITE (IO%IU6,'(I6,6F18.6)')  BDIR,TMP(1,BDIR),TMP(2,BDIR),TMP(3,BDIR)
         ENDDO
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*)

      ENDIF

      !< restore the original symmetry
      IF (LNMR_SYM_RED.AND.SYMM%ISYM>0) THEN
         CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,DYN%EFOR,T_INFO%LSFOR, &
              T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
              SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
              SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
              SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
              T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#else
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
              SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
              T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#endif
         CALL RE_GEN_LAYOUT(GRID,WDES,KPOINTS,LATT_CUR,LATT_INI,IO%IU6,IO%IU0)
         IF (LUSEPEAD()) CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)
         CALL REALLOCATE_WAVE(W,GRID,WDES,NONL_S,T_INFO,P,LATT_CUR)

         CALL DEALLOCATE_KPOINTS_TRANS(KPOINTS_TRANS)
      ENDIF

      ! Do not deallocate ABS_CHEM_SHIFT, needed for Magres writing
      !DEALLOCATE(ABS_CHEM_SHIFT)
      IF (NUCIND.AND..NOT.LNICSALL) THEN
        DEALLOCATE(POSROT,MAPNICS)
      ENDIF

      ALLOCATE(CORE_SHIFT(T_INFO%NTYP))

      CALL MLR_CORE_SHIFTS(P,T_INFO,LATT_CUR,XC,IO%IU6,CORE_SHIFT,CORE_CHI)

      IF (IO%IU6>=0) THEN
         ! write csa info
         IF (LNMRSHIELD) THEN
            MODUS='SHIEL'
         ELSE
            MODUS='SHIFT'
         ENDIF

         IF (LNMRLEG.OR.(INMRPRINT.LT.0)) THEN
            WRITE (IO%IU6,'(A)') ' ---------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  CSA tensor (J. Mason, Solid State Nucl. Magn. Reson. 2, 285 (1993))             '
            WRITE (IO%IU6,'(A)') ' ---------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '             EXCLUDING G=0 CONTRIBUTION             INCLUDING G=0 CONTRIBUTION    '
            WRITE (IO%IU6,'(A)') '         -----------------------------------   -----------------------------------'
            IF (LNMRSHIELD) THEN
               WRITE (IO%IU6,'(A)') '  ATOM   ISO_SHIELD        SPAN        SKEW    ISO_SHIELD        SPAN        SKEW '
            ELSE
               WRITE (IO%IU6,'(A)') '  ATOM    ISO_SHIFT        SPAN        SKEW     ISO_SHIFT        SPAN        SKEW '
            ENDIF
            WRITE (IO%IU6,'(A)') ' ---------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  (absolute, valence only)                                                        '

            DO NI=1,T_INFO%NIONS
               TMP_TENS_AT(:,:) = TMP_TENS(:,:,NI)
               CALL DIAGSHIFT(TMP_TENS_AT(1,1),MODUS,ISO1,SPAN1,SKEW1)
               TMP_TENS_AT(:,:) = TMP_TENS(:,:,NI) + TMP(:,:) ! add suscept, already SHIFTSHIELD converted above
               CALL DIAGSHIFT(TMP_TENS_AT(1,1),MODUS,ISO2,SPAN2,SKEW2)
               IF (SPAN1.LT.0.00005_q) SKEW1=0._q
               IF (SPAN2.LT.0.00005_q) SKEW2=0._q
               WRITE (IO%IU6,'(I6,X,3F12.4,2X,3F12.4)') NI, ISO1, SPAN1, SKEW1, ISO2, SPAN2, SKEW2
            ENDDO

            WRITE (IO%IU6,'(A)') ' ---------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  (absolute, valence and core)                                                    '

            DO NI=1,T_INFO%NIONS
               TMP_TENS_AT(:,:) = TMP_TENS(:,:,NI)
               TMP_TENS_AT(1,1) = TMP_TENS_AT(1,1) + SHIFTSHIELD * CORE_SHIFT(T_INFO%ITYP(NI))
               TMP_TENS_AT(2,2) = TMP_TENS_AT(2,2) + SHIFTSHIELD * CORE_SHIFT(T_INFO%ITYP(NI))
               TMP_TENS_AT(3,3) = TMP_TENS_AT(3,3) + SHIFTSHIELD * CORE_SHIFT(T_INFO%ITYP(NI))
               CALL DIAGSHIFT(TMP_TENS_AT(1,1),MODUS,ISO1,SPAN1,SKEW1)
               
               ! Overwrite ABS_CHEM_SHIFT and store core + valence + suscept for magres writing
               ABS_CHEM_SHIFT(:,:,NI) = TMP_TENS(:,:,NI) + TMP(:,:) ! add suscept, already SHIFTSHIELD converted above
               ABS_CHEM_SHIFT(1,1,NI) = ABS_CHEM_SHIFT(1,1,NI) + SHIFTSHIELD * (CORE_SHIFT(T_INFO%ITYP(NI)) + CORE_CHI/CHI_FACT)
               ABS_CHEM_SHIFT(2,2,NI) = ABS_CHEM_SHIFT(2,2,NI) + SHIFTSHIELD * (CORE_SHIFT(T_INFO%ITYP(NI)) + CORE_CHI/CHI_FACT)
               ABS_CHEM_SHIFT(3,3,NI) = ABS_CHEM_SHIFT(3,3,NI) + SHIFTSHIELD * (CORE_SHIFT(T_INFO%ITYP(NI)) + CORE_CHI/CHI_FACT)
               TMP_TENS_AT(:,:)=ABS_CHEM_SHIFT(:,:,NI)
               CALL DIAGSHIFT(TMP_TENS_AT(1,1),MODUS,ISO2,SPAN2,SKEW2)
               
               IF (SPAN1.LT.0.00005_q) SKEW1=0._q
               IF (SPAN2.LT.0.00005_q) SKEW2=0._q
               WRITE (IO%IU6,'(I6,X,3F12.4,2X,3F12.4)') NI, ISO1, SPAN1, SKEW1, ISO2, SPAN2, SKEW2
            ENDDO

            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  IF SPAN.EQ.0, THEN SKEW IS ILL-DEFINED                                          '
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
         ELSE
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  CSA tensor (J. Mason, Solid State Nucl. Magn. Reson. 2, 285 (1993))'
            IF (LNMRSHIELD) THEN
               WRITE (IO%IU6,'(A)') '  for chemical shielding  (including the isotropic core contribution)'
            ELSE
               WRITE (IO%IU6,'(A)') '  (including the isotropic core contribution)'
            ENDIF
            WRITE (IO%IU6,*)
            WRITE (IO%IU6,*)     ' BDIR labels direction of applied magnetic field (i.e. B0)'
            WRITE (IO%IU6,*)     ' For BDIR==1, B0 along x-axis, etc.'
            WRITE (IO%IU6,*)     ' Induced field listed along cartesian directions for each BDIR'
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '                 EXCLUDING G=0 CONTRIBUTION              INCLUDING G=0 CONTRIBUTION        '
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  ion  BDIR            X           Y           Z                X           Y           Z'
            IF (LNMRSHIELD) THEN
               WRITE (IO%IU6,'(A)') '          (   iso_shield        span        skew ) (   iso_shield        span        skew )'
            ELSE
               WRITE (IO%IU6,'(A)') '          (    iso_shift        span        skew ) (    iso_shift        span        skew )'
            ENDIF
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'

            DO NI=1, T_INFO%NIONS
               ! Add core contrib to tensor (permanently changes TMP_TENS(1,1,NI))
               TMP_TENS(1,1,NI) = TMP_TENS(1,1,NI) + SHIFTSHIELD * CORE_SHIFT(T_INFO%ITYP(NI))
               TMP_TENS(2,2,NI) = TMP_TENS(2,2,NI) + SHIFTSHIELD * CORE_SHIFT(T_INFO%ITYP(NI))
               TMP_TENS(3,3,NI) = TMP_TENS(3,3,NI) + SHIFTSHIELD * CORE_SHIFT(T_INFO%ITYP(NI))

               ! Excl. G=0
               TMP_TENS_AT(:,:) = TMP_TENS(:,:,NI)
               CALL DIAGSHIFT(TMP_TENS_AT(1,1),MODUS,ISO1,SPAN1,SKEW1)
               IF (SPAN1.LT.0.00005_q) SKEW1=0._q

               ! Incl. G=0
               ! Overwrite ABS_CHEM_SHIFT and store core + valence + suscept for magres writing
               ABS_CHEM_SHIFT(:,:,NI) = TMP_TENS(:,:,NI) + TMP(:,:) ! add suscept, already SHIFTSHIELD converted above
               ABS_CHEM_SHIFT(1,1,NI) = ABS_CHEM_SHIFT(1,1,NI) + SHIFTSHIELD * (CORE_CHI/CHI_FACT)
               ABS_CHEM_SHIFT(2,2,NI) = ABS_CHEM_SHIFT(2,2,NI) + SHIFTSHIELD * (CORE_CHI/CHI_FACT)
               ABS_CHEM_SHIFT(3,3,NI) = ABS_CHEM_SHIFT(3,3,NI) + SHIFTSHIELD * (CORE_CHI/CHI_FACT)
               TMP_TENS_AT(:,:)=ABS_CHEM_SHIFT(:,:,NI)

               CALL DIAGSHIFT(TMP_TENS_AT(1,1),MODUS,ISO2,SPAN2,SKEW2)
               IF (SPAN2.LT.0.00005_q) SKEW2=0._q

               IF (NI.NE.1) THEN
                  WRITE (IO%IU6,*)
               ENDIF
                                                             !ion bdir  x y z
               WRITE (IO%IU6,'(I5,2X,I4,1X,3F12.4,5X,3F12.4)') NI,   1, TMP_TENS(:,1,NI), TMP_TENS_AT(:,1)
               WRITE (IO%IU6,'(5X,2X,I4,1X,3F12.4,5X,3F12.4)')       2, TMP_TENS(:,2,NI), TMP_TENS_AT(:,2)
               WRITE (IO%IU6,'(5X,2X,I4,1X,3F12.4,5X,3F12.4)')       3, TMP_TENS(:,3,NI), TMP_TENS_AT(:,3)
               WRITE (IO%IU6,'(10X,"( ",3F12.4," ) ( ",3F12.4," )")') ISO1, SPAN1, SKEW1, ISO2, SPAN2, SKEW2

               ! compose susceptibility tensor for Magres writing
               CHI_MAG(1,1) = CHI_MAG(1,1) + CORE_CHI
               CHI_MAG(2,2) = CHI_MAG(2,2) + CORE_CHI
               CHI_MAG(3,3) = CHI_MAG(3,3) + CORE_CHI

            ENDDO
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
            WRITE (IO%IU6,'(A)') '  IF SPAN.EQ.0, THEN SKEW IS ILL-DEFINED                                                   '
            WRITE (IO%IU6,'(A)') ' ------------------------------------------------------------------------------------------'
         ENDIF
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*)
      ENDIF

      DEALLOCATE(CORE_SHIFT)

      RETURN
      END SUBROUTINE MLR_SKELETON

      SUBROUTINE WRITE_NICS_FOR_POSNICS(IO,LATT_CUR,T_INFO,TAUNICS,NICSDATA_SYM,NPOS_NICS,&
         & NICSDATA_UNSYM,NPOSROT,LNMRLEG,LNMRSHIELD)
        USE base, ONLY: in_struct
        USE poscar, ONLY: latt, type_info
#ifdef VASP_HDF5
        USE vhdf5
#endif
        IMPLICIT NONE

        TYPE (in_struct) :: IO
        TYPE (latt) :: LATT_CUR
        TYPE (type_info) :: T_INFO
        REAL(q) :: NICSDATA_SYM(:,:,:) !< symmetrized NICS at positions in POSNICS
        REAL(q) :: NICSDATA_UNSYM(:,:,:) !< unsymmetrized NICS
        REAL(q) :: TAUNICS(:,:) !< positions provided in the POSNICS file
        LOGICAL :: LNMRLEG !< write legacy NMR output
        LOGICAL :: LNMRSHIELD !< print shielding, F: or shift
        INTEGER :: NPOS_NICS !< no of positions provided on POSNICS file
        INTEGER :: NPOSROT !< no of symmetry-related positions; dimension of NICSDATA_UNSYM

        ! local
        INTEGER NI, BDIR
#ifdef VASP_HDF5
        CHARACTER(LEN=40) :: GROUP_
#endif

! - OUTCAR writing
        WRITE (IO%IU6,*) '-------------------------------------------------------------'
        IF (LNMRSHIELD) THEN
          WRITE (IO%IU6,*) ' Nucleus Independent Chemical Shielding, NICS (absolute) '
        ELSE
          WRITE (IO%IU6,*) ' Nucleus Independent Chemical Shift, NICS (absolute) '
        ENDIF
        WRITE (IO%IU6,*) '-------------------------------------------------------------'
! legacy output (unsymmetrized)
        IF (LNMRLEG) THEN
          WRITE (IO%IU6,*) ' UNSYMMETRIZED TENSORS '
          DO NI=1,NPOSROT   ! .NE. NPOS_NICS IN GENERAL
            WRITE (IO%IU6,'(" nics ",I4)') NI
            DO BDIR=1,3
              WRITE(IO%IU6,'(3F18.6)') NICSDATA_UNSYM(:,BDIR,NI)
            ENDDO
          ENDDO
        ENDIF

! symmetriezed NICS to OUTCAR
        WRITE (IO%IU6,*) ' SYMMETRIZED TENSORS '
        DO NI=1,NPOS_NICS
          WRITE (IO%IU6,'(" nics ",I4)') NI
          DO BDIR=1,3
            WRITE(IO%IU6,'(3F18.6)') NICSDATA_SYM(:,BDIR,NI)
          ENDDO
        ENDDO
        WRITE (IO%IU6,'(A/)') &
         &         '-------------------------------------------------------------'

! - vaspout.h5 writing
#ifdef VASP_HDF5
        GROUP_="/posnics"
        CALL VH5_WRITE_NICS_FOR_POSNICS(IH5OUTFILEID,GRP_RESULTS//trim(GROUP_),NICSDATA_SYM,TAUNICS,NPOS_NICS,LNMRSHIELD)
#endif

      END SUBROUTINE WRITE_NICS_FOR_POSNICS

!> *********************** SUBROUTINE MLR_CORE_SHIFTS *************************
!>
!> Calculation of the core contribution to the chemical shielding.
!>
!> ****************************************************************************
      SUBROUTINE MLR_CORE_SHIFTS(P,T_INFO,LATT_CUR,XC,IU6,CORE_SHIFT,CORE_CHI)
      USE pseudo
      USE poscar
      USE rhfatm
      USE radial
      USE lattice
      USE constant
      USE morbitalmag
      USE setexm_struct_def, ONLY : xc_info
      TYPE (potcar), TARGET :: P(:)
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (xc_info) XC
      INTEGER IU6
      REAL(q) CORE_SHIFT(*),CORE_CHI
      ! local variables
      TYPE (potcar), POINTER :: PP
      TYPE (atomic) :: ATOM
      INTEGER NT,I
      REAL(q), ALLOCATABLE :: RHOC(:)
      REAL(q) RES

      REAL(q), PARAMETER :: C=137.037_q

      IF (IU6<0) RETURN

      DO NT=1,T_INFO%NTYP
         CORE_SHIFT(NT)=0._q
      ENDDO
      CORE_CHI=0._q

      DO NT=1,T_INFO%NTYP
         PP=>P(NT)
         IF (.NOT.ASSOCIATED(PP%ATOMIC_N)) THEN
            WRITE(IU6,'(/A,I4)') ' WARNING: MLR_CORE_SHIFTS: can not compute core contributions for atomic type:',NT
            WRITE(IU6,'(A/)')    '          CORE CONTRIBUTION TO SHIELDING AND SUSCEPTIBILITY ARE INCORRECT'
            CYCLE
         ENDIF
         CALL RHFATM_DFATOM(PP,ATOM,XC)
         ALLOCATE(RHOC(ATOM%R%NMAX))
         RHOC=0._q
         DO I=1,ATOM%NSCORE
!           RHOC(:)=RHOC(:)+2*(2*ATOM%L(I)+1)*ATOM%W(:,I)**2/ATOM%R%R(:)
            RHOC(:)=RHOC(:)+2*(2*ATOM%L(I)+1)*(ATOM%A(:,I)**2+ATOM%B(:,I)**2/C/C)/ATOM%R%R(:)
         ENDDO
         CALL SIMPI(ATOM%R,RHOC,RES)
         CORE_SHIFT(NT)=-RES/2._q/C*2._q/3._q/C*AUTOA

         RHOC=0._q
         DO I=1,ATOM%NSCORE
            RHOC(:)=RHOC(:)+2*(2*ATOM%L(I)+1)*ATOM%W(:,I)**2*ATOM%R%R(:)**2
         ENDDO
         CALL SIMPI(ATOM%R,RHOC,RES)
         CORE_CHI=CORE_CHI-T_INFO%NITYP(NT)*RES/C/C*AUTOA

         DEALLOCATE(RHOC)
         CALL RHFATM_DEALLOCW(ATOM)
      ENDDO

      WRITE(IU6,'(A)')  ' --------------------------------------------------------------------------'
      WRITE(IU6,'(A)')  '  Core NMR properties        '

      IF (LNMRSHIELD) THEN
         WRITE(IU6,'(/A)') '  typ  El   Core shielding (ppm)'
      ELSE
         WRITE(IU6,'(/A)') '  typ  El   Core shift (ppm)'
      ENDIF
      WRITE(IU6,'(A)')  ' ----------------------------'
      DO NT=1,T_INFO%NTYP
         WRITE(IU6,'(X,I4,2X,A,2X,F14.7)') NT,P(NT)%ELEMENT,CORE_SHIFT(NT)*1.E6_q*SHIFTSHIELD
      ENDDO
      WRITE(IU6,'(A)')  ' ----------------------------'

      WRITE(IU6,'(/A,F10.2,X,A)') '  Core contribution to magnetic susceptibility:',CORE_CHI*1.E5_q,' 10^-6 cm^3/mole'
      WRITE(IU6,'(A)')  ' --------------------------------------------------------------------------'
      WRITE(IU6,*)
      WRITE(IU6,*)

      DO NT=1,T_INFO%NTYP
         CORE_SHIFT(NT) = CORE_SHIFT(NT)*1.E6_q
      ENDDO
      CORE_CHI = CORE_CHI*1.E5_q

      RETURN
      END SUBROUTINE MLR_CORE_SHIFTS

!> ********************** SUBROUTINE MLR_CORE_SHIFTS_FC ***********************
!>
!> Experimental subroutine for the calculation of the hyperfine contribution
!> to the core chemical shielding...
!>
!> Borrowed corepol parts from hyperfine_rad
!> Subroutine works for the center parts of the star (i.e. "LQ")
!> For the on-star stuff, it does not work so well,
!> problems with convergence of taylor series in hyperfine_FC_ZORA
!>
!> Could be merged with MLR_CORE_SHIFTS_FCmy like the hyperfine routines are
!>
!> ****************************************************************************
      SUBROUTINE MLR_CORE_SHIFTS_FC(T_INFO,NI,PP,XC,RHOAE,RHOPS,B_OUT)
      USE poscar
      USE pseudo
      USE constant
      USE radial
      USE rhfatm
      USE setexm
      USE setexm_struct_def, ONLY: xc_info
      USE ldalib, ONLY: EXCOR_DER
      USE morbitalmag
      USE hyperfine

      IMPLICIT NONE

      TYPE (potcar), POINTER :: PP
      TYPE (type_info) T_INFO
      TYPE (xc_info) XC

      REAL(q) RHOAE(:,:,:), RHOPS(:,:,:)
      INTEGER NI
      REAL(q) B_OUT

      !local variables
!#ifdef corepol
!      TYPE (atomic), ALLOCATABLE, SAVE :: ATOM (:)
      TYPE(atomic) :: ATOM
      REAL(q), ALLOCATABLE :: C(:,:)
      REAL(q) DV(PP%R%NMAX),RHOC(PP%R%NMAX)
      REAL(q) RHOUP,RHODN
      REAL(q) EXC,EXCD(2),EXCDD(2,2),EXCDDD(2,2,2)
!#endif

      REAL(q), ALLOCATABLE :: RHOAE_(:,:,:), RHOPS_(:,:,:)

      INTEGER I,J,K,NT
      REAL(q) SCALE,SUM(5),RHOW(PP%R%NMAX),RES
      REAL(q) TMP_AE

!      INTEGER IDIR
!      INTEGER INDX(3)
!      INDX(1)=2
!      INDX(2)=3
!      INDX(3)=4

      ALLOCATE(RHOAE_(PP%R%NMAX,(PP%LMAX+1)*(PP%LMAX+1),2))
      ALLOCATE(RHOPS_(PP%R%NMAX,(PP%LMAX+1)*(PP%LMAX+1),2))
      RHOAE_ = 0._q
      RHOPS_ = 0._q
      CALL RAD_MAG_DENSITY(RHOAE,RHOAE_,PP%LMAX,PP%R)
      CALL RAD_MAG_DENSITY(RHOPS,RHOPS_,PP%LMAX,PP%R)

      !Y00
      SCALE=1/(2._q*SQRT(PI))

!#ifdef corepol
!      sdir: DO IDIR=1,3
         !need to do RHF on atom... only per atomtype?
!         IF (.NOT.ALLOCATED(ATOM)) THEN
!            ALLOCATE(ATOM(T_INFO%NTYP)); ATOM(:)%Z=-1
!         ENDIF
!         NT=T_INFO%ITYP(NI)
!         IF (ATOM(NT)%Z<0) CALL RHFATM_DFATOM(PP,ATOM(NT),XC)
         CALL RHFATM_DFATOM(PP,ATOM,XC)

         DV=0._q

         !setup spin polarization of valence
         DO K=1, PP%R%NMAX
!            RHOUP=(RHOAE(K,1,1)+RHOAE(K,1,INDX(IDIR))+PP%RHOAE(K))*SCALE/2._q/PP%R%R(K)**2
!            RHODN=(RHOAE(K,1,1)-RHOAE(K,1,INDX(IDIR))+PP%RHOAE(K))*SCALE/2._q/PP%R%R(K)**2
            RHOUP=(RHOAE_(K,1,1)+RHOAE_(K,1,2)+PP%RHOAE(K))*SCALE/2._q/PP%R%R(K)**2
            RHODN=(RHOAE_(K,1,1)-RHOAE_(K,1,2)+PP%RHOAE(K))*SCALE/2._q/PP%R%R(K)**2
            CALL EXCOR_DER(XC,RHOUP,RHODN,0,EXC,EXCD,EXCDD,EXCDDD,.FALSE.)
            DV(K)=(EXCD(1)-EXCD(2))/SCALE/2._q
         ENDDO

         DEALLOCATE(RHOAE_, RHOPS_)

         !combine valence spin polarization with core...
         ALLOCATE(C(ATOM%NSCORE,ATOM%NS)); C=0._q
         RHOC=0._q

         DO I=1,ATOM%NSCORE
            IF (ATOM%L(I)/=0) CYCLE
            DO J=ATOM%NSCORE+1,ATOM%NS
               IF (ATOM%L(J)/=0) CYCLE
               RHOW=0._q
               DO K=1,PP%R%NMAX
   !              IF (PP%RDEP>0.AND.PP%R%R(K)-PP%RDEP>-5E-3) EXIT
                  RHOW(K)=ATOM%W(K,J)*DV(K)*ATOM%W(K,I)
               ENDDO
               CALL SIMPI(PP%R,RHOW,RES)
               C(I,J)=RES/(ATOM%E(I)-ATOM%E(J))
               DO K=1,PP%R%NMAX
                  RHOC(K)=RHOC(K)+4._q*C(I,J)*ATOM%W(K,J)*ATOM%W(K,I)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(C)

         RHOC(:)=RHOC(:)*SCALE
!#endif
!         CALL HYPERFINE_FC_ZORA(PP,RHOC,B_OUT(IDIR))
!      ENDDO sdir
         CALL HYPERFINE_FC_ZORA(PP,RHOC,B_OUT)

      RETURN
      END SUBROUTINE MLR_CORE_SHIFTS_FC

!> ********************* SUBROUTINE MLR_CORE_SHIFTS_FCmy **********************
!>
!> Experimental subroutine for the calculation of the hyperfine contribution
!> to the core chemical shielding...
!>
!> Borrowed corepol_my parts from hyperfine_rad
!> Subroutine works for the center parts of the star (i.e. "LQ")
!> For the on-star stuff, it does not work so well,
!> problems with convergence of taylor series in hyperfine_FC_ZORA
!>
!> Could be merged with MLR_CORE_SHIFTS_FC like the hyperfine routines are
!>
!> ****************************************************************************
      SUBROUTINE MLR_CORE_SHIFTS_FCmy(T_INFO,NI,NDIM,LYMAX,LMMAX,PP,XC,RHOAE,RHOPS,B_OUT)
      USE us
      USE poscar
      USE pseudo
      USE constant
      USE radial
      USE cl
      USE setexm
      USE setexm_struct_def, ONLY: xc_info
      USE morbitalmag
      USE hyperfine

      IMPLICIT NONE

      TYPE (potcar), POINTER :: PP
      TYPE (type_info) T_INFO
      TYPE (xc_info) XC

      REAL(q) RHOAE(:,:,:), RHOPS(:,:,:)
      INTEGER NI,NDIM,LYMAX,LMMAX
      REAL(q) B_OUT

      !local variables
!#ifdef corepol_my
      REAL(q) RHOC(PP%R%NMAX),RHO(PP%R%NMAX,1,2), POT(PP%R%NMAX,1,2),DMY(PP%R%NMAX),RHOP(PP%R%NMAX,2),DRHO
      REAL(q), ALLOCATABLE :: POTAE(:,:,:)
      REAL(q), ALLOCATABLE :: WUP(:,:),WDN(:,:),EUP(:),EDN(:)
      INTEGER, ALLOCATABLE :: N(:),LC(:)
      INTEGER MAXNL

      REAL(q), PARAMETER :: MIX=0.2_q, RCONV=1.E-8_q
      INTEGER, PARAMETER :: NSCF=100

      REAL(q) DV(NDIM)
      REAL(q) RHOUP,RHODN
      REAL(q) EXC,EXCD(2),EXCDD(2,2),EXCDDD(2,2,2)
!#endif

      REAL(q), ALLOCATABLE :: RHOAE_(:,:,:), RHOPS_(:,:,:)

      INTEGER I,J,K,NT
      REAL(q) SCALE,SUM(5),RHOW(NDIM),RES
      REAL(q) TMP_AE, DUMMY

      !Y00
      SCALE=1/(2._q*SQRT(PI))

      ALLOCATE(POTAE(NDIM, LMMAX, 4))

      !WRITE (*,*) "TOM RAD_POT START"
      !Hard-coded LASPH
      CALL RAD_POT(PP%R, XC, 2, LYMAX, PP%LMAX_CALC, .TRUE., RHOAE, PP%RHOAE, PP%POTAE, POTAE, DUMMY, DUMMY)
      !WRITE (*,*) "TOM RAD_POT DONE"

      !ALTHOUGH THE 4's below should be 2's, seg fault on 2...
      !WRITE (*,*) PP%R%NMAX, (PP%LMAX+1)*(PP%LMAX+1)
      !WRITE (*,*) NDIM, LMMAX

      ALLOCATE(RHOAE_(PP%R%NMAX, (PP%LMAX+1)*(PP%LMAX+1), 2))
      ALLOCATE(RHOPS_(PP%R%NMAX, (PP%LMAX+1)*(PP%LMAX+1), 2))
      RHOAE_ = 0._q
      RHOPS_ = 0._q
      !WRITE (*,*) "TOM RAD_MAG_DENSITY"
      CALL RAD_MAG_DENSITY(RHOAE,RHOAE_,LYMAX,PP%R)
      CALL RAD_MAG_DENSITY(RHOPS,RHOPS_,LYMAX,PP%R)
      !WRITE (*,*) "TOM RAD_MAG_DENSITY"

      CALL PUSH_XC_TYPE(XC,XC_DATA_PP)

      !WRITE (*,*) "TOM WF START"

      CALL CL_INIT_CORE_CONF(PP,MAXNL)
      ALLOCATE(WUP(PP%R%NMAX,MAXNL),WDN(PP%R%NMAX,MAXNL),N(MAXNL),LC(MAXNL),EUP(MAXNL),EDN(MAXNL))

      CALL SET_CORE_WF(PP%RHOAE,PP%POTAE,PP%R,PP%ZCORE,PP%ZVALF_ORIG, &
     &   WUP,N,LC,EUP,(POTAE(1:PP%R%NMAX,1,1)+POTAE(1:PP%R%NMAX,1,2))/2._q)

      !WRITE (*,*) "TOM WF DONE"

      RHOC=0._q
      DO I=1,MAXNL
         DO K=1,PP%R%NMAX
            RHOC(K)=RHOC(K)+2._q*(WUP(K,I))*REAL(2*LC(I)+1,KIND=q)
         ENDDO
      ENDDO
      RHOC=RHOC*SCALE

      DO I=1,PP%R%NMAX
         RHO(I,1,1)=RHOAE_(I,1,1)+RHOC(I)
      ENDDO
!     RHO(:,1,1)=RHOAE(:,1,1)+PP%RHOAE(:)
      DO I=1,PP%R%NMAX
         RHO(I,1,2)=RHOAE_(I,1,2)
      ENDDO

      DMY=0._q

      ! TOM: Loop does not include indication of not converging
      ! start scf cycle here
      RHOP(:,1)=RHOC(:); RHOP(:,2)=0._q
      scf: DO J=1,NSCF

         CALL RAD_POT(PP%R,XC,2,0,0,.FALSE.,RHO,DMY,DMY,POT,DUMMY,DUMMY)

         ! spin up
         CALL SET_CORE_WF(DMY,DMY,PP%R,PP%ZCORE,PP%ZVALF_ORIG, &
        &   WUP,N,LC,EUP,POT(1:PP%R%NMAX,1,1)) !,NMAX= RNMAX_CL)
         ! spin down
         CALL SET_CORE_WF(DMY,DMY,PP%R,PP%ZCORE,PP%ZVALF_ORIG, &
        &   WDN,N,LC,EDN,POT(1:PP%R%NMAX,1,2)) !,NMAX= RNMAX_CL)

#ifdef debug
         DO I=1,MAXNL
            WRITE(*,'(a,i4,x,a,i4,x,a,f14.7,x,a,f14.7)') 'n:',n(i),'l:',lc(i),'eup:',eup(i),'edn:',edn(i)
         ENDDO
#endif
         RHO=0._q
         DO I=1,MAXNL
            DO K=1,PP%R%NMAX
               RHO(K,1,1)=RHO(K,1,1)+(WUP(K,I)+WDN(K,I))*REAL(2*LC(I)+1,KIND=q)
               RHO(K,1,2)=RHO(K,1,2)+(WUP(K,I)-WDN(K,I))*REAL(2*LC(I)+1,KIND=q)
            ENDDO
         ENDDO
         RHO=RHO*SCALE

         ! check for convergence
         DRHO=0._q
         DO K=1,PP%R%NMAX
            DRHO=DRHO+ABS(RHO(K,1,2)-RHOP(K,2))
         ENDDO
         !WRITE (*,*) "TOM FCmy", J, DRHO
         IF (DRHO<RCONV) EXIT scf

         RHO(:,1,1)=RHOP(:,1)+MIX*(RHO(:,1,1)-RHOP(:,1))
         RHO(:,1,2)=RHOP(:,2)+MIX*(RHO(:,1,2)-RHOP(:,2))

         RHOP(:,1)=RHO(:,1,1)
         RHOP(:,2)=RHO(:,1,2)

         DO I=1,PP%R%NMAX
            RHO(I,1,1)=RHO(I,1,1)+RHOAE_(I,1,1)
         ENDDO
         DO I=1,PP%R%NMAX
            RHO(I,1,2)=RHO(I,1,2)+RHOAE_(I,1,2)
         ENDDO

      ENDDO scf
#ifdef debug
      WRITE(*,*) 'iscf:',J
      DO I=1,MAXNL
         WRITE(*,'(a,i4,x,a,i4,x,a,f14.7,x,a,f14.7)') 'n:',n(i),'l:',lc(i),'eup:',eup(i),'edn:',edn(i)
      ENDDO
#endif

      RHOC=0._q
      DO I=1,MAXNL
         DO K=1,PP%R%NMAX
            RHOC(K)=RHOC(K)+(WUP(K,I)-WDN(K,I))*REAL(2*LC(I)+1,KIND=q)
         ENDDO
      ENDDO
      RHOC=RHOC*SCALE

      DEALLOCATE(WUP,WDN,N,LC,EUP,EDN)
      DEALLOCATE(POTAE, RHOAE_, RHOPS_)

      CALL CL_CLEAR_CORE_CONF

      ! Switch back to requested xc-type
      CALL POP_XC_TYPE(XC)

      CALL HYPERFINE_FC_ZORA(PP,RHOC,B_OUT)
      RETURN
      END SUBROUTINE MLR_CORE_SHIFTS_FCmy

!> ************************** SUBROUTINE MLR_B_MAIN ***************************
!>
!> Called from within loop over magnetic field direction (bloop) in
!> MLR_SKELETON. Treats response for each cartesian direction under the
!> applied field (iloop).
!>
!> ****************************************************************************
      SUBROUTINE MLR_B_MAIN( &
     &   HAMILTONIAN,W,GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,KPOINTS,LATT_CUR,LATT_INI,T_INFO,SYMM,P, &
     &   NONL_S,NONLR_S,LMDIM,CDIJ,CQIJ,SV,E,INFO,IRDMAX,IO,BDIR,XC,NICSDATA,RPHI)
      USE prec
      USE base
      USE pseudo
      USE poscar
      USE wave_high
      USE nonl_high
      USE hamil_struct_def
      USE mgrid
      USE lattice
      USE constant
      USE kpoints_change
      USE david
      USE subrot
      USE paw
      USE choleski
      USE main_mpi
      USE morbitalmag
      USE charge                   
      USE setexm_struct_def, ONLY : xc_info

      USE fock_ace, ONLY : LFOCK_ACE,FOCK_ACE_ALLOCW,WACE,FOCK_ACE_CONSTRUCT,FOCK_ACE_DEALLOCW, WACE
      USE fock, ONLY: LHFCALC, RESETUP_FOCK_WDES
      USE fock_dbl, ONLY : FOCK_ALL_DBLBUF

      IMPLICIT NONE
      TYPE (ham_handle) HAMILTONIAN
      TYPE (wavespin) W
      TYPE (grid_3d) GRID
      TYPE (grid_3d) GRID_SOFT               !< soft grid for pseudized potentials/ charge etc.
      TYPE (grid_3d) GRIDC                   !< full, fine grid
      TYPE (grid_3d) GRIDUS
      TYPE (transit) C_TO_US
      TYPE (transit) SOFT_TO_C
      TYPE (kpoints_struct) KPOINTS
      TYPE (latt) LATT_CUR
      TYPE (latt) LATT_INI
      TYPE (type_info) T_INFO
      TYPE (symmetry) SYMM
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct)NONLR_S
      TYPE (energy) E
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (potcar), TARGET :: P(:)
      INTEGER LMDIM,IRDMAX
      OVERLAP CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      OVERLAP CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)
      INTEGER BDIR !< cartesian component of B-field, 1=x, 2=y, 3=z
      TYPE (xc_info) XC
      COMPLEX(q), INTENT(OUT)  :: NICSDATA(:,:)
      COMPLEX(qs), OPTIONAL :: RPHI(:,:,:,:,:)

      RGRID :: S_BARE(W%WDES%GRID%RL%NP,3)   !< three component vector storing current density in real space,
                                             !> we should be careful allocating this
      RGRID :: S_DENS(W%WDES%GRID%RL%NP,3)   !< three component vector storing spin    density in real space, GIL SOC
                                             !> we should really allocate this, and only if LSOSHIFT GIL SOC

      ! local variables
      TYPE (wavespin) WGKPQ !< contains G u_i v_k+q, k| u>
      TYPE (wavespin) WXI
      TYPE (wavespin) WXI_1   !< vGv
      TYPE (wavespin) WXI_2   !< vGv
      TYPE (wavespin) WXI_3   !< vGv

      TYPE (wavedes), TARGET  :: WDES_kpq

      TYPE (nonl_struct) NONL_D
      TYPE (nonlr_struct), TARGET :: NONLR_D

      TYPE (nonl_struct) NONL_S_tmp

      INTEGER IDIR !< cartesian component of k-point shift
      INTEGER VDIR !< cartesian component of velocity operator, 1=x, 2=y, 3=z
      INTEGER NK,ISP,IQ
      INTEGER N,NELM
      INTEGER NSIM
      REAL(q) DISPL(3),DISPL_CART(3)
      INTEGER ICOUEV
      REAL(q) RMS,DESUM1,TOTEN
      REAL(q) WSCAL

      REAL(q) B_TMP(3)     !< helper array
      REAL(q) B_TMP_FC(3)  !< helper array SOC Fermi Contact
      REAL(q) B_TMP_DP(3)  !< helper array SOC Dipolar
      REAL(q) MAXTMP
      TYPE (potcar), POINTER :: PP
      COMPLEX(q) CRHODE(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      COMPLEX(q) CRHODE_AUG(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ) !< augmentation strengths
      COMPLEX(q), ALLOCATABLE :: CRHODE_SYM(:,:)
      !> because the current augmentation has to be done in PS_BIOT_SAVART_DR the strengths have to be
      !> processed and carried around separately. alternatively the calls to CALC_B_MAGATOM_LR might be
      !> merged so that CRHODE and CRHODE_AUG can be identical arrays, however, this makes for more
      !> complications
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      COMPLEX(q), ALLOCATABLE :: S_DENS_RCP(:,:)                   !< for magnetization density SOC
      COMPLEX(q), ALLOCATABLE :: S_DENS_RCP_FINE(:,:)              !< for magnetization density SOC
      OVERLAP,ALLOCATABLE :: COCC(:,:,:),COCC_IM(:,:)
      REAL(q), ALLOCATABLE :: RHOAE(:,:,:), RHOPS(:,:,:)
      REAL(q), ALLOCATABLE :: RHOLM(:,:),RHOLM_HELP(:)
      INTEGER NI, NIP, NT, IDP, I, J, ISUSP
      INTEGER M, NDIM, LYMAX, LMMAX
      INTEGER, EXTERNAL :: MAXL_AUG
      COMPLEX(q) TMP(3)

      COMPLEX(q) NAB_QQ(3)
      REAL(q) RTMP
      INTEGER ICART_K_DIR
      INTEGER IQDIR, JQDIR

      REAL(q) STENCIL(7,3)
      DATA STENCIL /  1.0_q,            -2.0_q,           1.0_q,  0.0_q,           0.0_q,             0.0_q,  0.0_q, &
                     -0.083333333333_q,  1.3333333333_q, -2.5_q,  1.3333333333_q, -0.083333333333_q,  0.0_q,  0.0_q, &
                      0.011111111111_q, -0.15_q,          1.5_q, -2.7222222222_q,  1.5_q,            -0.15_q, 0.011111111111_q /

      INTEGER CHI_BARE_STENCIL

      ! To do this check here is not very nice because the scf calculation
      ! is already done. But in the ORBITALMAG_READER we don't know how
      ! many k points the user wants to use. Note: We must do this check
      ! before doubling the k-point grid for the linear response calculation.
      IF (WRT_NMRCUR>0.AND.SYMM%ISYM>0.AND..NOT.((W%WDES%NKPTS==1.AND.W%WDES%VKPT(1,1)**2+W%WDES%VKPT(2,1)**2+W%WDES%VKPT(3,1)**2<1E-10_q))) THEN
         CALL vtutor%alert( "Setting WRT_NMRCUR>0 is incompatible"//&
             " with ISYM>0 in bulk systems (that inlude more than a"//&
             " single k point at 0). Hence, writing the"//&
             " nmr current density will be skipped. \n"//&
             "Set ISYM<=0 to write the nmr current density for"//&
             " a bulk system!")
         WRT_NMRCUR=0
      ENDIF
#ifdef exact_diagonalisation
      INTEGER IFLAG
#endif
      CHI_BARE_STENCIL=ICHIBARE

      IF (LSOSHIFT) THEN
         ALLOCATE(S_DENS_RCP(GRID_SOFT%MPLWV,3))
         ALLOCATE(S_DENS_RCP_FINE(GRIDC%MPLWV,3))
         NDIM=0
         DO NT=1,T_INFO%NTYP
            IF (ASSOCIATED(P(NT)%QPAW)) THEN
               NDIM=MAX(NDIM, P(NT)%R%NMAX)
            END IF
         ENDDO
         LYMAX =MAXL_AUG(T_INFO%NTYP,P)
         LMMAX=(LYMAX+1)**2
         ALLOCATE(RHOAE( NDIM, LMMAX, W%WDES%NCDIJ))
         ALLOCATE(RHOPS( NDIM, LMMAX, W%WDES%NCDIJ))
         ALLOCATE(RHOLM( LMDIM*LMDIM, W%WDES%NCDIJ))
         ALLOCATE(RHOLM_HELP( LMDIM*LMDIM))
         ALLOCATE(COCC(LMDIM,LMDIM,W%WDES%NCDIJ),COCC_IM(LMDIM,LMDIM))
         ALLOCATE(CRHODE_SYM(LMDIM,LMDIM))
      ENDIF

      NSIM=W%WDES%NSIM*2
#ifdef MPI
      NSIM=((W%WDES%NSIM*2+W%WDES%COMM_INTER%NCPU-1)/W%WDES%COMM_INTER%NCPU)*W%WDES%COMM_INTER%NCPU
#endif

      S_BARE = 0._q
      IF (LSOSHIFT) S_DENS = 0._q

      ! allocate before doubling the k-point grid
      CALL ALLOCW(W%WDES,WGKPQ)
      CALL ALLOCW(W%WDES,WXI)
      IF (LVGVCALC) THEN !< vGv
         CALL ALLOCW(W%WDES,WXI_1)
         CALL ALLOCW(W%WDES,WXI_2)
         CALL ALLOCW(W%WDES,WXI_3)
      ENDIF

      CALL CHECK_FULL_KPOINTS

      ! original number of k-points
      NKPTS_ORIG=W%WDES%NKPTS

      ! double the number of k-points
#ifdef oldsym
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LSORBIT.AND..NOT.W%WDES%LSPIRAL, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0,W%WDES%VKPT(:,1:NKPTS_ORIG))
#else
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0,W%WDES%VKPT(:,1:NKPTS_ORIG))
#endif
      CALL KPAR_SYNC_ALL(W%WDES,W)
      CALL RE_GEN_LAYOUT(GRID,W%WDES,KPOINTS,LATT_CUR,LATT_INI,IO%IU6,IO%IU0)
      CALL REALLOCATE_WAVE(W,GRID,W%WDES,NONL_S,T_INFO,P,LATT_CUR)

      IF (LHFCALC) CALL RESETUP_FOCK_WDES(W%WDES,LATT_CUR,LATT_INI,IO%IU6)

!     CALL SETWDES(W%WDES,WDES1,0)
!     CALL NEWWAV(W1, WDES1, .FALSE.)

      !> WGKPQ and WXI will live at k+q, their WDES remains the same except WDES%NKPTS, WDES%VKPT,
      !> WDES%DATAKE should point to the upper panel of NKPTS_ORIG k-points
      WDES_kpq=W%WDES
      WDES_kpq%NKPTS=NKPTS_ORIG
      WDES_kpq%VKPT=>W%WDES%VKPT(:,NKPTS_ORIG+1:W%WDES%NKPTS)
      WDES_kpq%DATAKE=>W%WDES%DATAKE(:,:,NKPTS_ORIG+1:W%WDES%NKPTS)

      WGKPQ%WDES=>WDES_kpq; WXI%WDES=>WDES_kpq
      IF (LVGVCALC) THEN !< vGv
         WXI_1%WDES=>WDES_kpq
         WXI_2%WDES=>WDES_kpq
         WXI_3%WDES=>WDES_kpq
      ENDIF

      CRHODE_AUG=0._q

130   FORMAT (5X, //, &
      &'------------------------------------ EXTERNAL FIELD DIRECTION ', &
      &I4,'  ------------------------------------'//)

      IF (IO%IU6>=0) THEN
         WRITE(IO%IU6,130) BDIR
      ENDIF

      iloop: DO IDIR=1,3
         ! skip the component along the B-field
         IF (IDIR==BDIR) CYCLE iloop

         ! the direction of the velocity operator we need here is \hat{B}_bdir x \hat{q}_idir
         VDIR=6-BDIR-IDIR

         ! setup  r | p_i >
         IF (INFO%LREAL) THEN
!           NONLR_D=NONLR_S
            CALL NONLR_ASSIGN(NONLR_D,NONLR_S)
            CALL NONLR_ALLOC_CRREXP(NONLR_D)
         ELSE
            CALL NONL_ALLOC_DER(NONL_S,NONL_D)
!           CALL SPHER_DER(GRID,NONL_D,P,W%WDES,LATT_CUR,VDIR)
         ENDIF

         DISPL_CART=0; DISPL_CART(IDIR)=DQ
         DISPL=DISPL_CART
         CALL KARDIR(1, DISPL(1), LATT_CUR%A)

         qloop: DO IQ=-CHI_BARE_STENCIL,CHI_BARE_STENCIL
            IF (IO%IU0>=0) WRITE(*,'("BDIR=",I2,2X,"IDIR=",I2,2X,"IQ=",I3)') BDIR,IDIR,IQ

            WGKPQ%CW=(0._q,0._q); WGKPQ%CPROJ=zero; WXI%CW=(0._q,0._q); WXI%CPROJ=zero
            IF (LVGVCALC) THEN !< vGv
               WXI_1%CW=(0._q,0._q); WXI_1%CPROJ=zero
               WXI_2%CW=(0._q,0._q); WXI_2%CPROJ=zero
               WXI_3%CW=(0._q,0._q); WXI_3%CPROJ=zero
            ENDIF

!> ==============================================================================
!>  Start with the outer-points of the k-point star: ``k+q''
!>  (See, e.g., eq. 67 of arXiv:2504.17380 [cond-mat.mtrl-sci], first term RHS.)
!>
!>  first half of array contains KPOINTS on center of star ``k''
!>  second half of array contains KPOINTS on displacement ``+q''
!>
!>  ``+q'' and ``-q'' are treated in separate iterations of the loop
!>  but added to the same array, which is processed after the loop.
!>
!> ==============================================================================

            ! set additional k-points to k+q
            DO NK=1,NKPTS_ORIG
               W%WDES%VKPT(:,NKPTS_ORIG+NK)=W%WDES%VKPT(:,NK)+IQ*DISPL(:)
!              W%WDES%VKPT(:,NKPTS_ORIG+NK)=W%WDES%VKPT(:,NK)
            ENDDO

            !BUG FIX FOR GNU COMPILER (IF-statement)
            IF (.NOT.INFO%LREAL) THEN                                   !FIX
               CALL SPHER_DER(GRID,NONL_D,P,W%WDES,LATT_CUR,VDIR)
            ENDIF                                                       !FIX

            ! fill additional wave functions
            IF (PRESENT(RPHI)) THEN
               DO ISP=1,W%WDES%ISPIN
                  DO NK=1,NKPTS_ORIG
                     DO N=1,W%WDES%NBANDS
                        W%CW(:,N,NKPTS_ORIG+NK,ISP)=W%CW(:,N,NK,ISP)+TPI*( &
                       &   RPHI(:,N,NK,ISP,1)*(0.0_q,1.0_q)*IQ*DISPL_CART(1)+ &
                       &   RPHI(:,N,NK,ISP,2)*(0.0_q,1.0_q)*IQ*DISPL_CART(2)+ &
                       &   RPHI(:,N,NK,ISP,3)*(0.0_q,1.0_q)*IQ*DISPL_CART(3))
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               DO ISP=1,W%WDES%ISPIN
                  DO NK=1,NKPTS_ORIG
                     DO N=1,W%WDES%NBANDS
                        W%CW(:,N,NKPTS_ORIG+NK,ISP)=W%CW(:,N,NK,ISP)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

            CALL SET_DATAKE(W%WDES,LATT_CUR%B)
            IF (NONLR_S%LREAL) THEN
               CALL RSPHER(GRID,NONLR_S,LATT_CUR)
            ELSE
               CALL SPHER(GRID,NONL_S,P,W%WDES,LATT_CUR,1)
            ENDIF
            NONL_S%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,W)

            CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)

!           IF (IQ/=0) THEN
#ifdef exact_diagonalisation
!              CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
!                          LMDIM,CDIJ,CQIJ,3,SV,T_INFO,P,IO%IU0,E%EXHF, &
!                          NKSTART=NKPTS_ORIG+1)
!               CALL EDDIAG_EXACT(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
!                                 LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IO%IU0,IO%IU6,E%EXHF)
               CALL EDDIAG_EXACT(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
                           LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IO%IU0,IO%IU6,E%EXHF,E%EXHF_ACFDT)
               CALL KPAR_SYNC_ALL(W%WDES,W)
#else
               ! let's do a bit of Davidson optimization
               W%FERTOT(:,NKPTS_ORIG+1:W%WDES%NKPTS,:)=W%FERTOT(:,1:NKPTS_ORIG,:)
               W%WDES%WTKPT(NKPTS_ORIG+1:W%WDES%NKPTS)=W%WDES%WTKPT(1:NKPTS_ORIG)

               IF (LFOCK_ACE()) CALL FOCK_ACE_ALLOCW(W%WDES)

               TOTEN=0; INFO%IALGO=8
               DO NELM=1,NSTEPS_DAVIDSON

                  IF (LFOCK_ACE()) THEN
                     ! temporarily set the weights for the orbitals beyond NKPOINTS_ORIG to zero
                     W%WDES%WTKPT(NKPTS_ORIG+1:W%WDES%NKPTS)=0
                     ! compute the action of the Fock potential on these orbitals
                     CALL FOCK_ALL_DBLBUF(W%WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ,EX=E%EXHF,NKSTART=NKPTS_ORIG+1,XI=WACE)
                     ! and construct the corresponding ACE representation of the Fock operator
                     CALL FOCK_ACE_CONSTRUCT(W%WDES,W,NKSTART=NKPTS_ORIG+1)
                     ! return the weights back to their appropriate values
                     W%WDES%WTKPT(NKPTS_ORIG+1:W%WDES%NKPTS)=W%WDES%WTKPT(1:NKPTS_ORIG)
                  ENDIF

                  CALL EDDAV(HAMILTONIAN,P,GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,NSIM, &
                 &     LMDIM,CDIJ,CQIJ,RMS,DESUM1,ICOUEV,SV,E%EXHF,IO%IU6,IO%IU0, &
                 &     LDELAY=.FALSE.,LSUBROTI=.TRUE.,LEMPTY=.FALSE.,LHF=.TRUE., &
                 &     NKSTART=NKPTS_ORIG+1)

                  E%EBANDSTR=BANDSTRUCTURE_ENERGY(W%WDES,W)

                  IF (IO%IU0>=0) WRITE(17, 200)      NELM,E%EBANDSTR,E%EBANDSTR-TOTEN,DESUM1,ICOUEV,RMS
                  IF (IO%IU0>=0) WRITE(IO%IU0, 200)  NELM,E%EBANDSTR,E%EBANDSTR-TOTEN,DESUM1,ICOUEV,RMS
                  TOTEN=E%EBANDSTR

 200              FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
                       &       I6,'  ',E10.3)

                  CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)

                  IF (ABS(DESUM1) < ABS(INFO%EDIFF) .AND. NELM>1) EXIT
               ENDDO

               CALL KPAR_SYNC_ALL(W%WDES,W)
#endif
!           ENDIF

            IF (LFOCK_ACE()) THEN
               ! temporarily set the weights for the orbitals beyond NKPOINTS_ORIG to zero
               W%WDES%WTKPT(NKPTS_ORIG+1:W%WDES%NKPTS)=0
               ! compute the action of the Fock potential on these orbitals
               CALL FOCK_ALL_DBLBUF(W%WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ,EX=E%EXHF,NKSTART=NKPTS_ORIG+1,XI=WACE)
               ! and construct the corresponding ACE representation of the Fock operator
               CALL FOCK_ACE_CONSTRUCT(W%WDES,W,NKSTART=NKPTS_ORIG+1)
               ! return the weights back to their appropriate values
               W%WDES%WTKPT(NKPTS_ORIG+1:W%WDES%NKPTS)=W%WDES%WTKPT(1:NKPTS_ORIG)
               ! and copy the upper panel into the lower one ...
               WACE%CW(:,:,1:NKPTS_ORIG,:)=WACE%CW(:,:,NKPTS_ORIG+1:W%WDES%NKPTS,:)
            ENDIF

            ! let's start the linear response stuff to get the elements | psi^(1e) >
            CALL MLR_PSI_RESPONSE_EMPTY(HAMILTONIAN,W,GRID,P,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,SV,INFO,IO,VDIR,WXI,WGKPQ)

            !< Susceptibility
            IF (LVGVCALC) THEN  !< vGv
               IF (.NOT.INFO%LREAL) CALL SPHER_DER(GRID,NONL_D,P,W%WDES,LATT_CUR,1)
               CALL MLR_COMMUTATOR(HAMILTONIAN,W,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,INFO,1,WXI_1)
               IF (.NOT.INFO%LREAL) CALL SPHER_DER(GRID,NONL_D,P,W%WDES,LATT_CUR,2)
               CALL MLR_COMMUTATOR(HAMILTONIAN,W,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,INFO,2,WXI_2)
               IF (.NOT.INFO%LREAL) CALL SPHER_DER(GRID,NONL_D,P,W%WDES,LATT_CUR,3)
               CALL MLR_COMMUTATOR(HAMILTONIAN,W,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,INFO,3,WXI_3)
               IF (.NOT.INFO%LREAL) CALL SPHER_DER(GRID,NONL_D,P,W%WDES,LATT_CUR,VDIR)
            ENDIF

            ! here we should get the contributions Q(q)
            DO ISUSP=1,2
               IF (ISUSP==1) CALL NABIJ_ASYM(W, WGKPQ, GRID, LATT_CUR, NKPTS_ORIG, INFO, NAB_QQ)
               IF (ISUSP==2) CALL NABIJ_VGV(W, WGKPQ, WXI_1, WXI_2, WXI_3,  GRID, LATT_CUR, NKPTS_ORIG, INFO, NAB_QQ) !     vGv
               DO ICART_K_DIR=1,3
                  IQDIR =6-IDIR-ICART_K_DIR
                  JQDIR =6-IDIR-VDIR
                  RTMP = EIJK(IQDIR,IDIR,ICART_K_DIR)*EIJK(JQDIR,IDIR,VDIR)
                  IF (ABS(RTMP)<.5_q) CYCLE
!                 IF ((IQ==-1).OR.(IQ==1)) QQ(IQDIR,JQDIR)=QQ(IQDIR,JQDIR)  &
!    &                                  + NAB_QQ(ICART_K_DIR)*EIJK(IQDIR,IDIR,ICART_K_DIR)*EIJK(JQDIR,IDIR,VDIR)
!                 IF (IQ==0)               QQ(IQDIR,JQDIR)=QQ(IQDIR,JQDIR)  &
!    &                           - 2._q * NAB_QQ(ICART_K_DIR)*EIJK(IQDIR,IDIR,ICART_K_DIR)*EIJK(JQDIR,IDIR,VDIR)
                  QQ(IQDIR,JQDIR,ISUSP)=QQ(IQDIR,JQDIR,ISUSP)+ &
                 &   NAB_QQ(ICART_K_DIR)*RTMP*STENCIL(IQ+CHI_BARE_STENCIL+1,CHI_BARE_STENCIL)
               ENDDO
               IF (.NOT.LVGVCALC) EXIT
            ENDDO
            !< Susceptibility

            IF (IQ==0.OR.ABS(IQ)>1) CYCLE qloop

            ! and the elements | \psi^(1o) >
            CALL MLR_PSI_RESPONSE_OCC(W,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CQIJ,IO,INFO,VDIR,WXI,WGKPQ)

            WGKPQ%CW=WGKPQ%CW*(0._q,-1._q)

            ! we should use the upper panel of NONL_S%QPROJ, i.e., the projectors at k+q
            NONL_S_tmp=NONL_S; NONL_S_tmp%NK=-1
            NONL_S_tmp%QPROJ=>NONL_S%QPROJ(:,:,:,NKPTS_ORIG+1:W%WDES%NKPTS,:)
            CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S_tmp,WGKPQ)

            ! the sign of \hat{B}_bdir x \hat{q}_idir is given by the Levi-Civita symbol
            WGKPQ%CW=EIJK(VDIR,BDIR,IDIR)*WGKPQ%CW; WGKPQ%CPROJ=EIJK(VDIR,BDIR,IDIR)*WGKPQ%CPROJ

            ! mupltiply with the sign of q
            WGKPQ%CW=WGKPQ%CW*IQ; WGKPQ%CPROJ=WGKPQ%CPROJ*IQ

            !< compute current contributions to S_bare(q). (if LSOSHIFT, also SDENS)
            !< eq. 68 arXiv:2504.17380 [cond-mat.mtrl-sci] (if LSOSHIFT, also eq. 76)
            CALL MLR_SBARE(W,WGKPQ,S_BARE,S_DENS,LSOSHIFT)

            ! compute current contributions to B_Dp.
            CALL DEPSUM_ASYM(W ,WGKPQ , W%WDES, LMDIM, CRHODE, W%WDES%LOVERL)   ! or cc?

            CALL US_FLIP_CMPLX(W%WDES, LMDIM, CRHODE, W%WDES%LOVERL, .FALSE.)
            ! now we should have charge and magnetization

            ! augmentation strenghts
            CRHODE_AUG(:,:,:,:)=CRHODE_AUG(:,:,:,:)+CRHODE(:,:,:,:)

            IDP=1 !< only paramagnetic contribution to field
                  !< (diamagnetic is a ``Q_R'' term.)

            DO NI=1,T_INFO%NIONS        !< index of ion

               NIP=NI_LOCAL(NI, W%WDES%COMM_INB)   !< local index of ion
               IF (NIP==0) CYCLE                   !< not stored locally return
               NT=T_INFO%ITYP(NI)                  !< type of ion
!              PP=> P(NT)
               PP=>PP_POINTER(P,NI,NT)

               IF ( ASSOCIATED(PP%QPAW)) THEN
                 IF (DO_LOCAL(NI)) THEN
                    !< first term RHS eq. 73 arXiv:2504.17380 [cond-mat.mtrl-sci]
                    CALL CALC_B_MAGATOM_LR(PP, CRHODE(:,:,NIP,1),BDIR,IDP,NIP,NI, &
                         B_TMP(:),LLRAUG)

                    B_TMP(:)=-B_TMP(:) !< for paramagnetic, -B_TMP
                    B_PARA_ONE_CENTR(:,NI,BDIR) = B_PARA_ONE_CENTR(:,NI,BDIR) + B_TMP(:)

! TOM 27-10-2023 fix for compiling vasp_std, but not nice fix
#ifndef realmode
                    IF (LSOSHIFT) THEN
                       RHOLM=0._q; RHOAE=0._q; RHOPS=0._q; COCC=0; COCC_IM=0
                       DO ISP=1,W%WDES%NCDIJ
                          CALL TRANS_RHOLM_ASYM_PREPARE(CRHODE(:,:,NIP,ISP),CRHODE_SYM(:,:), PP)     ! DANGER
                          CALL TRANS_RHOLM(CRHODE_SYM(:,:), RHOLM(:,ISP), PP)    ! pas op
                          CALL RAD_CHARGE(RHOAE(:,:,ISP), PP%R, RHOLM(:,ISP), PP%LMAX, PP%LPS, PP%WAE)
                          CALL RAD_CHARGE(RHOPS(:,:,ISP), PP%R, RHOLM(:,ISP), PP%LMAX, PP%LPS, PP%WPS)
                          CALL TRANS_RHOLM_IM( CRHODE_SYM(:,:), RHOLM_HELP, PP )
                          CALL TRANS_RHOLMI_IM(COCC_IM, RHOLM_HELP, PP)
                          COCC(:,:,ISP)=COCC(:,:,ISP)+(0._q,1._q)*COCC_IM(:,:)
                       ENDDO
                       !< eq. 78 arXiv:2504.17380 [cond-mat.mtrl-sci]
                       CALL CALC_B_MAGATOM_LR_FC_DP(PP, CRHODE(:,:,:,:), RHOPS, RHOAE, BDIR,NIP,NI, &
                            COCC,B_TMP_FC(:),B_TMP_DP(:))
                       B_SPIN_ONE_CENTR_FC(:,NI,BDIR)  = B_SPIN_ONE_CENTR_FC(:,NI,BDIR)  + B_TMP_FC(:)
                       B_SPIN_ONE_CENTR_DIP(:,NI,BDIR) = B_SPIN_ONE_CENTR_DIP(:,NI,BDIR) + B_TMP_DP(:)
                    ENDIF
#endif
                 ENDIF
               ENDIF
            ENDDO

            IF (LFOCK_ACE()) CALL FOCK_ACE_DEALLOCW

         ENDDO qloop

         ! some deallocation
         IF (INFO%LREAL) THEN
            CALL NONLR_DEALLOC_CRREXP(NONLR_D)
         ELSE
            CALL NONL_DEALLOC_DER(NONL_D)
         ENDIF

      ENDDO iloop

      ! here one has to divide by q, both soft and one centre
      DO I=1,W%WDES%GRID%RL%NP
         DO J=1,3
            S_BARE(I,J) = S_BARE(I,J)/DQ/2._q/TPI
         ENDDO
      ENDDO
      DO I=1,3
         DO J=1,T_INFO%NIONS
            B_PARA_ONE_CENTR(I,J,BDIR) = B_PARA_ONE_CENTR(I,J,BDIR)/DQ/2._q/TPI
         ENDDO
      ENDDO
      IF (LSOSHIFT) THEN
         DO I=1,W%WDES%GRID%RL%NP
            DO J=1,3
               S_DENS(I,J) = S_DENS(I,J)/DQ/2._q/TPI
            ENDDO
         ENDDO
         DO I=1,3
            DO J=1,T_INFO%NIONS
               B_SPIN_PS(I,J,BDIR)            = B_SPIN_PS(I,J,BDIR)           /DQ/2._q/TPI !TOM 16-01-2024 B_SPIN_PS has not yet been used?
               B_SPIN_ONE_CENTR_FC(I,J,BDIR)  = B_SPIN_ONE_CENTR_FC(I,J,BDIR) /DQ/2._q/TPI
               B_SPIN_ONE_CENTR_DIP(I,J,BDIR) = B_SPIN_ONE_CENTR_DIP(I,J,BDIR)/DQ/2._q/TPI
            ENDDO
         ENDDO
      ENDIF

      CRHODE_AUG(:,:,:,:)=CRHODE_AUG(:,:,:,:)/DQ/2._q/TPI

!> ============================================================================
!>
!>  Now, we calculate the terms on the points ``k'', i.e. the ``Q_R'' terms.
!>  (eq. 57 arXiv:2504.17380 [cond-mat.mtrl-sci])
!>
!>  First prepare the hamiltonian H^{(1)}_{Q_R} in array CLQIJ, subroutines:
!>  - MLR_CALC_CLQIJ
!>  - MLR_CALC_ZEEMAN_ONE
!>  - MLR_CALC_SPINORB_ONE_*
!>
!>  MLR_PSI_RESPONSE_LQ prepares for calculation of current contributions and
!>  includes Zeeman term of H^{(1)}
!>
!> ============================================================================
      !CLQIJ=0._q ! make 0 here, because next routines are subsequently going to add to this array
      !CLQIJ set to 0 inside subroutine

      ! here one has to add j_bare,Q_R and j_dp,Q_R
      CALL MLR_CALC_CLQIJ(W%WDES,P,T_INFO,LMDIM,BDIR,CLQIJ) !< eq. 24 arXiv:2504.17380 [cond-mat.mtrl-sci]

      !GIL: Check CLQIJ is not set to 0 inside MLR_CALC_CLQIJ
      !Then it can be removed above, but must call routine before those
      !below, and we neednt modify it. Seems more elegant for ISPIN=1 and ISPIN=2

      !TOM: keep CLQIJ=0 inside MLR_CALC_CLQIJ,
      !need to call MLR_CALC_CLQIJ before the routines below and needn't
      !be modified (seems more elegant for ISPIN=1 and ISPIN=2)

      IF (LSOSHIFT) THEN
          ! one has to add to CLQIJ the one-center Zeeman and the one-center spin-orbit due to the vectorpotential
          ! to get one-centre Zeeman perturbation, switch on either MLR_CALC_ZEEMAN_ONE or MLR_CALC_ZEEMAN_ONE_SIMPLE
          CALL MLR_CALC_ZEEMAN_ONE(W%WDES,P,T_INFO,LMDIM,BDIR,CLQIJ) !< eq. 23 arXiv:2504.17380 [cond-mat.mtrl-sci]
          !CALL MLR_CALC_ZEEMAN_ONE_SIMPLE(W%WDES,P,T_INFO,LMDIM,BDIR,CLQIJ)
          ! this is a spin-orbit gauge term AND MAYBE COULD BE REMOVED (very small contribution)
          !CALL MLR_CALC_SPINORB_ONE(W%WDES,P,T_INFO,LMDIM,BDIR,CLQIJ)
          CALL MLR_CALC_SPINORB_ONE_ALTS(W%WDES,P,T_INFO,LMDIM,BDIR,CLQIJ) !< eq. 25 arXiv:2504.17380 [cond-mat.mtrl-sci]
          CALL MLR_CALC_SPINORB_ONE_ALTD(W%WDES,P,T_INFO,LMDIM,BDIR,CLQIJ) !< eq. 25 arXiv:2504.17380 [cond-mat.mtrl-sci]
          ! MLR_CALC_CLQIJ, MLR_CALC_ZEEMAN_ONE and MLR_CALC_SPINORB_ONE_* could be made into a single subroutine.
      ENDIF

! test
!     CALL DUMP_DLLMM( "CLQIJ",CLQIJ(:,:,1,1), P(1))
!     CALL DUMP_DLLMM( "CQIJ",  CQIJ(:,:,1,1), P(1))
! test

      ! everything here will be at k
      W%WDES%VKPT(:,NKPTS_ORIG+1:W%WDES%NKPTS)=W%WDES%VKPT(:,1:NKPTS_ORIG)
      W%CW(:,:,NKPTS_ORIG+1:W%WDES%NKPTS,:)=W%CW(:,:,1:NKPTS_ORIG,:)

      CALL SET_DATAKE(W%WDES,LATT_CUR%B)
      IF (NONLR_S%LREAL) THEN
         CALL RSPHER(GRID,NONLR_S,LATT_CUR)
      ELSE
         CALL SPHER(GRID,NONL_S,P,W%WDES,LATT_CUR,1)
      ENDIF
      NONL_S%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,W)

      CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)

#ifdef exact_diagonalisation
      !CALL EDDIAG_EXACT(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
      !            LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IO%IU0,IO%IU6,E%EXHF)
      CALL EDDIAG_EXACT(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,W%WDES,SYMM, &
                  LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IO%IU0,IO%IU6,E%EXHF,E%EXHF_ACFDT)
#endif
      CALL KPAR_SYNC_ALL(W%WDES,W)

      IF (LFOCK_ACE()) THEN
         CALL FOCK_ACE_ALLOCW(W%WDES)
         ! compute the action of the Fock potential on orbitals NK <= NKPOINTS_ORIG
         CALL FOCK_ALL_DBLBUF(W%WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ,EX=E%EXHF,NKSTOP=NKPTS_ORIG,XI=WACE)
         ! and construct the corresponding ACE representation of the Fock operator
         CALL FOCK_ACE_CONSTRUCT(W%WDES,W,NKSTOP=NKPTS_ORIG)
      ENDIF

      WGKPQ%WDES%VKPT=>W%WDES%VKPT(:,1:NKPTS_ORIG); WXI%WDES%VKPT=>WGKPQ%WDES%VKPT
      WGKPQ%WDES%DATAKE=>W%WDES%DATAKE(:,:,1:NKPTS_ORIG); WXI%WDES%DATAKE=>WGKPQ%WDES%DATAKE

      WGKPQ%CW=(0._q,0._q); WGKPQ%CPROJ=zero; WXI%CW=(0._q,0._q); WXI%CPROJ=zero

      CALL MLR_PSI_RESPONSE_LQ(HAMILTONIAN,W,GRID,P,NONL_S,NONLR_S,LATT_CUR,LMDIM,CDIJ,CQIJ,CLQIJ,BDIR,SV,LSOSHIFT,INFO,IO,WXI,WGKPQ)

      IF (LFOCK_ACE()) CALL FOCK_ACE_DEALLOCW

      NONL_S%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,WGKPQ)

      ! compute current contributions to S_bare(q).
      ! note that WGKPQ now contains result of MLR_PSI_RESPONSE_LQ, and not the points k+q
      !< eq. 58 arXiv:2504.17380 [cond-mat.mtrl-sci] (if LSOSHIFT also eq. 60)
      CALL MLR_SBARE(W,WGKPQ,S_BARE,S_DENS,LSOSHIFT)

      ! S_DENS MERGING AND WRITING AND FOURIER TRANSFORM
      IF (LSOSHIFT) THEN
         ! merge the spin density from all nodes
         DO I=1,3
#ifdef realmode
            CALLMPI( M_sum_d(W%WDES%COMM_INTER, S_DENS(1,I), W%WDES%GRID%RL%NP))
            CALLMPI( M_sum_d(W%WDES%COMM_KINTER,S_DENS(1,I), W%WDES%GRID%RL%NP))
#else
            CALLMPI( M_sum_z(W%WDES%COMM_INTER, S_DENS(1,I), W%WDES%GRID%RL%NP))
            CALLMPI( M_sum_z(W%WDES%COMM_KINTER,S_DENS(1,I), W%WDES%GRID%RL%NP))
#endif
         ENDDO
         ! this is all complex
         DO M=1,W%WDES%GRID%RL%NP
            DO I=1,3
               TMP(I)=S_DENS(M,I)
            ENDDO
            CALL CDIRKAR(1,TMP,LATT_CUR%B)
            DO I=1,3
               S_DENS(M,I)=TMP(I)
            ENDDO
         ENDDO
         ! write spin density, for debugging propose
!        IF (IO%IU6>=0) THEN
!           OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'SCAR'//DIR_TEXT(BDIR),STATUS='UNKNOWN')
!           DO I=1,3
!              MAXTMP=0._q
!              DO J=1,GRID_SOFT%NGX*GRID_SOFT%NGY*GRID_SOFT%NGZ
!                 IF (REAL(S_DENS(J,I),q).GT.MAXTMP) MAXTMP=REAL(S_DENS(J,I),q)
!              ENDDO
!              WRITE(99,'(3I5,1X,E17.11)') GRID_SOFT%NGX,GRID_SOFT%NGY,GRID_SOFT%NGZ,MAXTMP
!              WRITE(99,'(5(1X,E17.11))') &
!          &     (REAL(S_DENS(J,I),q),J=1,GRID_SOFT%NGX*GRID_SOFT%NGY*GRID_SOFT%NGZ)
!           ENDDO
!           DO I=1,3
!              MAXTMP=0._q
!              DO J=1,GRID_SOFT%NGX*GRID_SOFT%NGY*GRID_SOFT%NGZ
!                 IF (AIMAG(S_DENS(J,I)).GT.MAXTMP) MAXTMP=AIMAG(S_DENS(J,I))
!              ENDDO
!              WRITE(99,'(3I5,1X,E17.11)') GRID_SOFT%NGX,GRID_SOFT%NGY,GRID_SOFT%NGZ,MAXTMP
!              WRITE(99,'(5(1X,E17.11))') &
!          &     (AIMAG(S_DENS(J,I)),J=1,GRID_SOFT%NGX*GRID_SOFT%NGY*GRID_SOFT%NGZ)
!           ENDDO
!           CLOSE(99)
!        ENDIF
         ! bring the spin density to reciprocal space
         S_DENS_RCP=0._q
         S_DENS_RCP_FINE=0._q
         DO I=1,3
            CALL FFT_RC_SCALE(S_DENS(1,I),S_DENS_RCP(1,I),GRID_SOFT)
            ! set the spin-density of unbalanced lattic-vectors to 0
            CALL SETUNB(S_DENS_RCP(1,I),GRID_SOFT)
            ! bring to full grid
            IF (.NOT.W%WDES%LOVERL) THEN
               CALL RC_ADD(S_DENS_RCP(1,I),1.0_q,S_DENS_RCP(1,I),0.0_q,S_DENS_RCP_FINE(1,I),GRID_SOFT)
            ELSE
               CALL ADD_GRID(GRIDC, GRID_SOFT, SOFT_TO_C,S_DENS_RCP(1,I),S_DENS_RCP_FINE(1,I))
               CALL SETUNB_COMPAT(S_DENS_RCP_FINE(1,I),GRIDC)
            ENDIF
         ENDDO
         ! let's try to call this...
         CALL FERMICONT_DIPOLAR_PW_S_DENS(T_INFO,LATT_CUR,GRIDC,W%WDES%NCDIJ,S_DENS_RCP_FINE,B_SPIN_PS(:,:,BDIR))
      ENDIF

      ! compute current contributions to B_Dp.
      ! note that WGKPQ now contains result of MLR_PSI_RESPONSE_LQ, and not the points k+q
      CALL DEPSUM_ASYM(W ,WGKPQ , W%WDES, LMDIM, CRHODE, W%WDES%LOVERL)   ! or cc?

      CALL US_FLIP_CMPLX(W%WDES, LMDIM, CRHODE, W%WDES%LOVERL, .FALSE.)

      CRHODE_AUG(:,:,:,:)=CRHODE_AUG(:,:,:,:)+CRHODE(:,:,:,:)

      IDP=1                       !< only paramagnetic contribution to field

      DO NI=1,T_INFO%NIONS        !< index of ion

         NIP=NI_LOCAL(NI, W%WDES%COMM_INB)   !< local index of ion
         IF (NIP==0) CYCLE                   !< not stored locally return
         NT=T_INFO%ITYP(NI)                  !< type of ion
!        PP=> P(NT)
         PP=>PP_POINTER(P,NI,NT)

         IF ( ASSOCIATED(PP%QPAW)) THEN
           IF (DO_LOCAL(NI)) THEN
              !< eq. 59 arXiv:2504.17380 [cond-mat.mtrl-sci]
              CALL CALC_B_MAGATOM_LR(PP, CRHODE(:,:,NIP,1),BDIR,IDP,NIP,NI, &
                   B_TMP(:),LLRAUG)

              B_TMP(:)=-B_TMP(:) !<for paramagnetic, -B_TMP
              B_PARA_ONE_CENTR(:,NI,BDIR) = B_PARA_ONE_CENTR(:,NI,BDIR) + B_TMP(:)

!TOM 27-10-2023
#ifndef realmode
              IF (LSOSHIFT) THEN
                 RHOLM=0._q; RHOAE=0._q; RHOPS=0._q; COCC=0; COCC_IM=0
                 DO ISP=1,W%WDES%NCDIJ
                    ! either uncomment TRANS_RHOLOM_ASYM_PREPARE of TRANS_RHOLM_ASYM_PREPARE_SIMPLE.
                    ! convert "asymmetric" CRHODE by adding the transpose of the c.c. and deviding by 2.
                    ! 01-11-2023; changed to be equal to above 10-01-2023: note that above we choose the other
                    CALL TRANS_RHOLM_ASYM_PREPARE(CRHODE(:,:,NIP,ISP),CRHODE_SYM(:,:), PP)            ! danger
!                    CALL TRANS_RHOLM_ASYM_PREPARE_SIMPLE(CRHODE(:,:,NIP,ISP),CRHODE_SYM(:,:), PP)     ! danger
                    ! for debugging
!                   if (ni.eq.1) then
!                      write (*,*)
!                      write (*,*) 'NI,ISP(ISP=1,NCDIJ),BDIR',ni,isp,bdir
!                      write (*,*) 'CRHODE real'
!                      do i11=1,lmdim
!                         write (*,'(i3,100f12.5)') i11,(REAL(CRHODE(i11,i22,nip,isp)),i22=1,lmdim)
!                      enddo
!                      write (*,*) 'CRHODE_SYM real'
!                      do i11=1,lmdim
!                         write (*,'(i3,100f12.5)') i11,(REAL(CRHODE_SYM(i11,i22)),i22=1,lmdim)
!                      enddo
!                      write (*,*) 'CRHODE imag'
!                      do i11=1,lmdim
!                         write (*,'(i3,100f12.5)') i11,(AIMAG(CRHODE(i11,i22,nip,isp)),i22=1,lmdim)
!                      enddo
!                      write (*,*) 'CRHODE_SYM imag'
!                      do i11=1,lmdim
!                         write (*,'(i3,100f12.5)') i11,(AIMAG(CRHODE_SYM(i11,i22)),i22=1,lmdim)
!                      enddo
!                      write (*,*)
!                   endif
                    CALL TRANS_RHOLM(CRHODE_SYM(:,:),  RHOLM(:,ISP), PP)    ! pas op
                    CALL RAD_CHARGE(RHOAE(:,:,ISP), PP%R, RHOLM(:,ISP), PP%LMAX, PP%LPS, PP%WAE)
                    CALL RAD_CHARGE(RHOPS(:,:,ISP), PP%R, RHOLM(:,ISP), PP%LMAX, PP%LPS, PP%WPS)
                    CALL TRANS_RHOLM_IM( CRHODE_SYM(:,:), RHOLM_HELP, PP )
                    CALL TRANS_RHOLMI_IM(COCC_IM, RHOLM_HELP, PP)
                    COCC(:,:,ISP)=COCC(:,:,ISP)+(0._q,1._q)*COCC_IM(:,:)
                 ENDDO
                 !< eq. 61 arXiv:2504.17380 [cond-mat.mtrl-sci]
                 CALL CALC_B_MAGATOM_LR_FC_DP(PP, CRHODE(:,:,:,:), RHOPS, RHOAE, BDIR,NIP,NI, &
                      COCC,B_TMP_FC(:),B_TMP_DP(:))
                 B_SPIN_ONE_CENTR_LQ_FC(:,NI,BDIR)  = B_SPIN_ONE_CENTR_LQ_FC(:,NI,BDIR)  + B_TMP_FC(:)
                 B_SPIN_ONE_CENTR_LQ_DIP(:,NI,BDIR) = B_SPIN_ONE_CENTR_LQ_DIP(:,NI,BDIR) + B_TMP_DP(:)

                  !WRITE(*,*) "TOM_CORE", BDIR, PP%ELEMENT

                  !MAXTMP=0._q
                  !CALL MLR_CORE_SHIFTS_FC(T_INFO,NI,PP,XC,RHOAE,RHOPS,MAXTMP)
                  !CALL MLR_CORE_SHIFTS_FCmy(T_INFO,NI,NDIM,LYMAX,LMMAX,PP,XC,RHOAE,RHOPS,MAXTMP)
                  !WRITE (*,*) "TOM_CORE", BDIR, PP%ELEMENT, MAXTMP
              ENDIF
#endif
           ENDIF
         ENDIF
      ENDDO

      ! here one has to call Biot-Savart for the PS current
      CALL PS_BIOT_SAVART_DR(S_BARE,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,&
        SOFT_TO_C,T_INFO,LATT_CUR,W%WDES,B_OUT_PARA_PS(:,:,BDIR),&
        CRHODE_AUG,LMDIM,IRDMAX,P,IO%IU0,IO%IU6,&
        BDIR,NICSDATA(:,:),SYMM%ISYM)

      ! calculate the diamagnetic one-centre field
      ! reset CRHODE appropriately for diamagnetic
      CALL DEPSUM_ASYM(W, W, W%WDES, LMDIM, CRHODE, W%WDES%LOVERL)

      CALL US_FLIP_CMPLX(W%WDES, LMDIM, CRHODE, W%WDES%LOVERL, .FALSE.)

      IDP=2                       ! only diamagnetic contribution to field

      DO NI=1,T_INFO%NIONS        !< index of ion

         NIP=NI_LOCAL(NI, W%WDES%COMM_INB)   !< local index of ion
         IF (NIP==0) CYCLE                   !< not stored locally return
         NT=T_INFO%ITYP(NI)                  !< type of ion
!        PP=> P(NT)
         PP=>PP_POINTER(P,NI,NT)

         IF ( ASSOCIATED(PP%QPAW)) THEN
           IF (DO_LOCAL(NI)) THEN
              !< eq. 50 arXiv:2504.17380 [cond-mat.mtrl-sci]
              CALL CALC_B_MAGATOM_LR(PP, CRHODE(:,:,NIP,1),BDIR,IDP,NIP,NI, &
                   B_TMP(:),LLRAUG)
              B_DIA_ONE_CENTR(:,NI,BDIR) = B_TMP(:)
           ENDIF
         ENDIF
      ENDDO

!> ============================================================================
!>
!> Shielding calculation finished, now summarize/wrap-up.
!> (Below covers legacy output)
!>
!> ============================================================================
      ! sum contributions to one-centre magnetic fields along BDIR
      CALLMPI( M_sum_d(W%WDES%COMM_KIN, B_PARA_ONE_CENTR(1,1,BDIR), 3*T_INFO%NIONS))
      CALLMPI( M_sum_d(W%WDES%COMM_KIN, B_DIA_ONE_CENTR(1,1,BDIR), 3*T_INFO%NIONS))

      IF (LSOSHIFT) THEN
         CALLMPI( M_sum_d(W%WDES%COMM_KIN, B_SPIN_ONE_CENTR_FC(1,1,BDIR), 3*T_INFO%NIONS))
         CALLMPI( M_sum_d(W%WDES%COMM_KIN, B_SPIN_ONE_CENTR_DIP(1,1,BDIR), 3*T_INFO%NIONS))
         CALLMPI( M_sum_d(W%WDES%COMM_KIN, B_SPIN_ONE_CENTR_LQ_FC(1,1,BDIR), 3*T_INFO%NIONS))
         CALLMPI( M_sum_d(W%WDES%COMM_KIN, B_SPIN_ONE_CENTR_LQ_DIP(1,1,BDIR), 3*T_INFO%NIONS))
      ENDIF

      ! write out one column of the tensor, should maybe go elsewhere
      B_OUT_PARA_PS(:,:,BDIR)     =B_OUT_PARA_PS(:,:,BDIR)   *1E6 *MOMTOMOM
      B_PARA_ONE_CENTR(:,:,BDIR)  =B_PARA_ONE_CENTR(:,:,BDIR)*1E6 *MOMTOMOM
      B_DIA_ONE_CENTR(:,:,BDIR)   =B_DIA_ONE_CENTR(:,:,BDIR) *1E6 !MOMTOMOM is in CALC_B_MAGATOM_LR for IDP=2

      DO NI=1,T_INFO%NIONS
         ABS_CHEM_SHIFT(:,BDIR,NI)=B_OUT_PARA_PS(:,NI,BDIR)+B_PARA_ONE_CENTR(:,NI,BDIR)+B_DIA_ONE_CENTR(:,NI,BDIR)
      ENDDO
      IF (LSOSHIFT) THEN
         B_SPIN_PS(:,:,BDIR)               =B_SPIN_PS(:,:,BDIR)              *1E6 *MOMTOMOM
         B_SPIN_ONE_CENTR_FC(:,:,BDIR)     =B_SPIN_ONE_CENTR_FC(:,:,BDIR)    *1E6 *MOMTOMOM
         B_SPIN_ONE_CENTR_DIP(:,:,BDIR)    =B_SPIN_ONE_CENTR_DIP(:,:,BDIR)   *1E6 *MOMTOMOM
         B_SPIN_ONE_CENTR_LQ_FC(:,:,BDIR)  =B_SPIN_ONE_CENTR_LQ_FC(:,:,BDIR) *1E6 *MOMTOMOM
         B_SPIN_ONE_CENTR_LQ_DIP(:,:,BDIR) =B_SPIN_ONE_CENTR_LQ_DIP(:,:,BDIR)*1E6 *MOMTOMOM

         DO NI=1,T_INFO%NIONS
            ABS_CHEM_SHIFT(:,BDIR,NI)=ABS_CHEM_SHIFT(:,BDIR,NI)   &
     &                               +B_SPIN_PS(:,NI,BDIR)+B_SPIN_ONE_CENTR_FC(:,NI,BDIR)   +B_SPIN_ONE_CENTR_DIP(:,NI,BDIR)  &
     &                                               +B_SPIN_ONE_CENTR_LQ_FC(:,NI,BDIR)+B_SPIN_ONE_CENTR_LQ_DIP(:,NI,BDIR)
         ENDDO
      ENDIF

      IF (LNMRLEG.OR.(INMRPRINT.GT.0)) THEN !< legacy output for single Bext direction
      IF (IO%IU6>=0) THEN
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) ' EXTERNAL FIELD APPLIED IN DIRECTION',BDIR
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '                       MAGNETIC FIELD (mu_B/Angstrom3/1e6)'
         WRITE (IO%IU6,*) ' ATOM                 X                 Y                 Z'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'

         IF (.NOT.LSOSHIFT) THEN
            DO NI=1,T_INFO%NIONS
               WRITE (IO%IU6,'(I6,3E18.6)') NI, &
     &         REAL(B_OUT_PARA_PS(1,NI,BDIR),q)+B_PARA_ONE_CENTR(1,NI,BDIR)+B_DIA_ONE_CENTR(1,NI,BDIR), &
     &         REAL(B_OUT_PARA_PS(2,NI,BDIR),q)+B_PARA_ONE_CENTR(2,NI,BDIR)+B_DIA_ONE_CENTR(2,NI,BDIR), &
     &         REAL(B_OUT_PARA_PS(3,NI,BDIR),q)+B_PARA_ONE_CENTR(3,NI,BDIR)+B_DIA_ONE_CENTR(3,NI,BDIR)
            ENDDO
         ELSE
            DO NI=1,T_INFO%NIONS
               WRITE (IO%IU6,'(I6,3E18.6)') NI, &
     &         REAL(B_OUT_PARA_PS(1,NI,BDIR),q)+B_PARA_ONE_CENTR(1,NI,BDIR)+B_DIA_ONE_CENTR(1,NI,BDIR) &
     &        +B_SPIN_PS(1,NI,BDIR)+B_SPIN_ONE_CENTR_FC(1,NI,BDIR)+B_SPIN_ONE_CENTR_DIP(1,NI,BDIR)     &
     &        +B_SPIN_ONE_CENTR_LQ_FC(1,NI,BDIR)+B_SPIN_ONE_CENTR_LQ_DIP(1,NI,BDIR)              ,&
     &         REAL(B_OUT_PARA_PS(2,NI,BDIR),q)+B_PARA_ONE_CENTR(2,NI,BDIR)+B_DIA_ONE_CENTR(2,NI,BDIR) &
     &        +B_SPIN_PS(2,NI,BDIR)+B_SPIN_ONE_CENTR_FC(2,NI,BDIR)+B_SPIN_ONE_CENTR_DIP(2,NI,BDIR)     &
     &        +B_SPIN_ONE_CENTR_LQ_FC(2,NI,BDIR)+B_SPIN_ONE_CENTR_LQ_DIP(2,BDIR,NI)              ,&
     &         REAL(B_OUT_PARA_PS(3,NI,BDIR),q)+B_PARA_ONE_CENTR(3,NI,BDIR)+B_DIA_ONE_CENTR(3,NI,BDIR) &
     &        +B_SPIN_PS(3,NI,BDIR)+B_SPIN_ONE_CENTR_FC(3,NI,BDIR)+B_SPIN_ONE_CENTR_DIP(3,NI,BDIR)     &
     &        +B_SPIN_ONE_CENTR_LQ_FC(3,NI,BDIR)+B_SPIN_ONE_CENTR_LQ_DIP(3,NI,BDIR)
            ENDDO
         ENDIF

         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  plane wave contribution'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,6E18.6)') NI,REAL(B_OUT_PARA_PS(1,NI,BDIR),q), &
     &                                   REAL(B_OUT_PARA_PS(2,NI,BDIR),q), &
     &                                   REAL(B_OUT_PARA_PS(3,NI,BDIR),q)
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  one center paramagnetic contribution'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_PARA_ONE_CENTR(1,NI,BDIR), &
     &                                   B_PARA_ONE_CENTR(2,NI,BDIR), &
     &                                   B_PARA_ONE_CENTR(3,NI,BDIR)
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  one center diamagnetic contribution'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_DIA_ONE_CENTR(1,NI,BDIR), &
     &                                   B_DIA_ONE_CENTR(2,NI,BDIR), &
     &                                   B_DIA_ONE_CENTR(3,NI,BDIR)
         ENDDO
      IF (LSOSHIFT) THEN
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  plane wave spin current contribution'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_SPIN_PS(1,NI,BDIR), &
     &                                   B_SPIN_PS(2,NI,BDIR), &
     &                                   B_SPIN_PS(3,NI,BDIR)
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  one center spin current contribution (Fermi contact) star'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_SPIN_ONE_CENTR_FC(1,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_FC(2,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_FC(3,NI,BDIR)
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  one center spin current contribution (dipolar) star'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_SPIN_ONE_CENTR_DIP(1,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_DIP(2,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_DIP(3,NI,BDIR)
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  one center spin current contribution (Fermi contact) LQ'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_SPIN_ONE_CENTR_LQ_FC(1,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_LQ_FC(2,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_LQ_FC(3,NI,BDIR)
         ENDDO
         WRITE (IO%IU6,*)
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*) '  one center spin current contribution (dipolar) LQ'
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         DO NI=1,T_INFO%NIONS
            WRITE (IO%IU6,'(I6,3E18.6)') NI,B_SPIN_ONE_CENTR_LQ_DIP(1,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_LQ_DIP(2,NI,BDIR), &
     &                                   B_SPIN_ONE_CENTR_LQ_DIP(3,NI,BDIR)
         ENDDO
      ENDIF

         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         IF (NUCIND.AND..NOT.LNICSALL) THEN
            WRITE (IO%IU6,*)
            WRITE (IO%IU6,*) '-------------------------------------------------------------'
            WRITE (IO%IU6,*) '  NICS at the position given in file POSNICS'
            WRITE (IO%IU6,*) '-------------------------------------------------------------'
            DO NI=1,NPOS_NICS
               WRITE (IO%IU6,'(I6,3E18.6)') NI,NICSDATA(NI,1), &
     &                                         NICSDATA(NI,2), &
     &                                         NICSDATA(NI,3)
            ENDDO
         ENDIF
         WRITE (IO%IU6,*) '-------------------------------------------------------------'
         WRITE (IO%IU6,*)
      ENDIF
      ENDIF

      CALL DEALLOCW(WGKPQ)
      CALL DEALLOCW(WXI)
      IF (LVGVCALC) THEN
         CALL DEALLOCW(WXI_1)
         CALL DEALLOCW(WXI_2)
         CALL DEALLOCW(WXI_3)
      ENDIF

      ! restore to the original number of k-points
#ifdef oldsym
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LSORBIT.AND..NOT.W%WDES%LSPIRAL, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,-1)
#else
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
           SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, &
           T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,-1)
#endif
      CALL RE_GEN_LAYOUT(GRID,W%WDES,KPOINTS,LATT_CUR,LATT_INI,IO%IU6,IO%IU0)
      CALL REALLOCATE_WAVE(W,GRID,W%WDES,NONL_S,T_INFO,P,LATT_CUR)

      IF (LSOSHIFT) THEN
         DEALLOCATE(S_DENS_RCP)
         DEALLOCATE(S_DENS_RCP_FINE)
         DEALLOCATE(RHOAE)
         DEALLOCATE(RHOPS)
         DEALLOCATE(RHOLM)
         DEALLOCATE(RHOLM_HELP)
         DEALLOCATE(COCC,COCC_IM)
         DEALLOCATE(CRHODE_SYM)
      ENDIF

      RETURN
      END SUBROUTINE MLR_B_MAIN

!> ******************** SUBROUTINE MLR_PSI_RESPONSE_EMPTY *********************
!>
!> YPM, Phys. Rev. B 76, 024401 (2017); Appendix B (B10)
!> Response from first order wavefunction can be split into
!> projection into empty and occupied subspaces...
!>
!> ****************************************************************************
      SUBROUTINE MLR_PSI_RESPONSE_EMPTY( &
     &   HAMILTONIAN,W0,GRID,P,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,SV,INFO,IO,IDIR,WXI,W1)
      USE prec
      USE base
      USE wave_high
      USE nonl_high
      USE lattice
      USE mgrid
      USE pseudo
      USE rmm_diis_mlr
      USE hamil_struct_def
      USE hamil
      USE vaspxml
      USE ini
      IMPLICIT NONE
      TYPE (ham_handle) HAMILTONIAN
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (potcar) P(:)
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (nonl_struct) NONL_D            ! k derivative of projector in real space
      TYPE (nonlr_struct) NONLR_D          ! projector times r (only CRREXP distinct)
      TYPE (wavespin) WXI                  ! stores H(1) - epsilon S(1) | phi_0>
      TYPE (wavespin) W1                   ! first order change of wavefunction
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (latt) LATT_CUR
      INTEGER LMDIM
      OVERLAP CDIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ), &
              CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),W0%WDES%NCDIJ)
      INTEGER IDIR !< component of velocity operator 1=x, 2=y, 3=z
      ! local variables
      TYPE (wavespin) W0_tmp
      TYPE (nonl_struct) NONL_S_tmp

      COMPLEX(q), POINTER :: SW0(:,:,:,:)=>NULL()

      INTEGER N
      REAL(q) TOTENL,TOTEN,DE

      INTEGER ICOUEV,IERROR
      REAL(q) DESUM,RMS(INFO%NELM),CSHIFT
      REAL(q) EBREAK_STORE
      LOGICAL LRESET
      REAL(q) FBREAK_MEAN

      ! calculate  v_k+q,k | u_nk >
      CALL MLR_COMMUTATOR(HAMILTONIAN,W0,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,INFO,IDIR,WXI)
      CALL MLR_COMMUTATOR_FOCK(W0,GRID,P,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,IDIR,WXI)
#ifdef apply_greens_function
      CALL MLR_APPLY_GREENS_FUNC(W0,WXI,W1)
#else
#ifdef provide_w1
      CALL MLR_APPLY_GREENS_FUNC(W0,WXI,W1)
#endif
      ! project onto empty subspace | xi_k+q > = P_e v_k+q,k | u_nk >
!     CALL START_TIMING("PROJ")
!     CALL MLR_PROJ_IN_EMPTY_SUBSPACE(W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)
      CALL MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS(W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)
!     CALL STOP_TIMING("PROJ",IO%IU6)

!     we need to trick LINEAR_RESPONSE_DIIS to work with the second panel of NKPTS_ORIG k-points
!     as far as some quantities are concerned (WDES%DATAKE, WDES%VKPT, W0%CW and NONL_S), but for
!     instance with eigenvalues from the lower panel (e(0)_k S_k+q - H_k+q)|u(1)_k+q> = |xi_k+q>,
!     (W0%CELTOT, and W0%CELEN).
      W0_tmp=W0
!     W0_tmp%WDES%NKPTS=NKPTS_ORIG
!     W0_tmp%WDES%VKPT=>W0%WDES%VKPT(:,NKPTS_ORIG+1:W0%WDES%NKPTS)
!     W0_tmp%WDES%DATAKE=>W0%WDES%DATAKE(:,:,NKPTS_ORIG+1:W0%WDES%NKPTS)
      W0_tmp%WDES=>W1%WDES

      W0_tmp%CW=>W0%CW(:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
      W0_tmp%CPROJ=>W0_tmp%CPROJ(:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)

      NONL_S_tmp=NONL_S; NONL_S_tmp%NK=-1
      NONL_S_tmp%QPROJ=>NONL_S%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)

      NONL_S_tmp%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S_tmp,WXI)

      EBREAK_STORE=INFO%EBREAK; INFO%EBREAK=INFO%EDIFF*0.25_q
      W0_tmp%OVER_BAND=.FALSE.

      IF (IO%IU0>=0) THEN
         WRITE(IO%IU0,142)
         WRITE(17,142)
      ENDIF
  142 FORMAT('       N       E                     dE             ' &
           ,'d eps       ncg     rms')

!  140 FORMAT (5X, //, &
!       &'----------------------------------------- Iteration ', &
!       &I4,'(',I4,')  ---------------------------------------'//)

!TOM: Although variable is IDIR, it actually contains "VDIR".
!  140 FORMAT (5X, //, &
!       &'--------------------------------------- IDIR', &
!       &I3,' Iteration (',I4,') ---------------------------------------'//)
  140 FORMAT (5X, //, &
       &'--------------------------------------- VDIR', &
       &I3,' Iteration (',I4,') ---------------------------------------'//)

!> ============================================================================
!> Main selfconsistent loop for linear response
!> ============================================================================
      TOTEN=0

      CSHIFT=0._q
      LRESET=.TRUE.
#ifdef provide_w1
      NONL_S_tmp%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S_tmp,W1)
      LRESET=.FALSE.; W1%CW=-W1%CW; W1%CPROJ=-W1%CPROJ
#endif
      FBREAK_MEAN=0.0
      iter: DO N=1,INFO%NELM
         CALL START_TIMING("LOOP")

         IF(IO%IU6>=0) WRITE(IO%IU6,140) IDIR,N
         CALL XML_TAG("scstep")

         CALL LINEAR_RESPONSE_DIIS(HAMILTONIAN,GRID,P,INFO,LATT_CUR,NONLR_S,NONL_S_tmp,W1,WXI,W0_tmp,W0_tmp%WDES, &
        &   LMDIM,CDIJ,CQIJ,RMS(N),DESUM,ICOUEV,SV,CSHIFT,IO%IU6,IO%IU0,LRESET,IERROR,SW0, FBREAK_MEAN)
!        CALL LINEAR_RESPONSE_DIIS(HAMILTONIAN,GRID,INFO,LATT_CUR,NONLR_S,NONL_S_tmp,W1,WXI,W0_tmp,W0_tmp%WDES, &
!       &   LMDIM,CDIJ,CQIJ,RMS(N),DESUM,ICOUEV,SV,CSHIFT,IO%IU6,IO%IU0,LRESET,IERROR)

         CALL MRG_CEL(W1%WDES,W1)
         LRESET=.FALSE.
         ! break criterion?
         W1%FERTOT(:,1:NKPTS_ORIG,:)=W0%FERTOT(:,1:NKPTS_ORIG,:)
         TOTENL=TOTEN
         TOTEN =BANDSTRUCTURE_ENERGY(W1%WDES,W1)
         DE= (TOTEN-TOTENL)

         IF (IO%IU0>=0) THEN
            WRITE(IO%IU0,200) N, TOTEN, DE, DESUM, ICOUEV, RMS(N)
            WRITE(17,200)     N, TOTEN, DE, DESUM, ICOUEV, RMS(N)
         ENDIF
 200     FORMAT('RMM: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5,I6,'  ',E10.3)
!> ============================================================================
!> total time used for this step
!> ============================================================================
         CALL SEPERATOR_TIMING(IO%IU6)
         CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')

         INFO%LABORT=.FALSE.
         IF(ABS(DESUM)<INFO%EDIFF.AND.ABS(DE)<INFO%EDIFF) INFO%LABORT=.TRUE.
!        IF(N>=3) THEN
!           IF (ABS((RMS(N)-RMS(N-1))/RMS(N))<1E-1 .AND. ABS((RMS(N)-RMS(N-2))/RMS(N))<1E-1) INFO%LABORT=.TRUE.
!        ENDIF

         IF (IO%IU6>=0)  THEN
            WRITE(IO%IU6,210) TOTEN

210   FORMAT(/ &
             '  free energy    TOTEN  = ',F18.8,' eV'/ &
             '  ---------------------------------------------------')

            IF (IO%LOPEN) CALL WFORCE(IO%IU6)
            IF (IO%LOPEN) CALL WFORCE(17)
         ENDIF

!> ============================================================================
!> xml related output
!> ============================================================================
         CALL XML_TAG("energy")
         IF (INFO%LABORT .OR. N==1) THEN
            CALL XML_ENERGY(TOTEN, TOTEN, TOTEN)
         ELSE
            CALL XML_ENERGY(TOTEN, TOTEN, TOTEN)
         ENDIF
         CALL XML_CLOSE_TAG

         CALL XML_CLOSE_TAG("scstep")


         IF (INFO%LABORT) THEN
            IF (IO%IU6>=0)  THEN

               WRITE(IO%IU6,131)
131       FORMAT (5X, //, &
       &  '------------------------ aborting loop because EDIFF', &
       &  ' is reached ----------------------------------------'//)
            ENDIF
            EXIT  iter
         ENDIF

      ENDDO iter

      CALL REDIS_PW_OVER_BANDS(W0_tmp%WDES,W0_tmp)
      DEALLOCATE(SW0); NULLIFY(SW0)

      INFO%EBREAK=EBREAK_STORE

      ! a bit heuristically
      W1%CW=-W1%CW
#endif

      RETURN
      END SUBROUTINE MLR_PSI_RESPONSE_EMPTY

!> ********************* SUBROUTINE MLR_PSI_RESPONSE_OCC **********************
!>
!> YPM, Phys. Rev. B 76, 024401 (2017); Appendix B (B15)
!> Response from first order wavefunction can be split into
!> projection into empty and occupied subspaces...
!>
!> ****************************************************************************
      SUBROUTINE MLR_PSI_RESPONSE_OCC( &
     & W0,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CQIJ,IO,INFO,IDIR,WXI,W1)
      USE prec
      USE base
      USE lattice
      USE wave_high
      USE nonl_high
      USE mgrid
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (nonl_struct) NONL_D            !< k derivative of projector in real space
      TYPE (nonlr_struct) NONLR_D          !< projector times r (only CRREXP distinct)
      TYPE (wavespin) WXI                  !< stores -i [r_idir,S_k+q,k] |u(0)_nk>
      TYPE (wavespin) W1
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (latt) LATT_CUR
      INTEGER LMDIM
      OVERLAP CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      INTEGER IDIR !< component of velocity operator 1=x, 2=y, 3=z
      ! local variables

      ! compute (1/i) [r_idir,S_k+q,k] |u(0)_nk> and store it in WXI
      CALL MLR_COMMUTATOR_RQ(W0,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CQIJ,INFO,IDIR,WXI)

      ! add - \sum_j | u(0)_j,k+q > < u(0)_j,k+q | (1/i) [r_idir,S_k+q,k] |u(0)_nk> to  | u(1)_n,k+q >
      ! where j runs over all occupied states
!     CALL MLR_PROJ_IN_OCC_SUBSPACE(W0,WXI,W1)
      CALL MLR_PROJ_IN_OCC_SUBSPACE_BLAS(W0,WXI,W1)

      RETURN
      END SUBROUTINE MLR_PSI_RESPONSE_OCC

!> ********************* SUBROUTINE MLR_PSI_RESPONSE_LQ ***********************
!>
!> YPM, Phys. Rev. B 76, 024401 (2017);
!> eq. 29 part to the right of current operator
!>
!> ****************************************************************************
      SUBROUTINE MLR_PSI_RESPONSE_LQ( &
     &   HAMILTONIAN,W0,GRID,P,NONL_S,NONLR_S,LATT_CUR,LMDIM,CDIJ,CQIJ,CLQIJ,BDIR,SV,LSOSHIFT,INFO,IO,WXI,W1)
      USE prec
      USE base
      USE mgrid
      USE pseudo
      USE wave_high
      USE nonl_high
      USE lattice
      USE rmm_diis_mlr
      USE hamil_struct_def
      USE hamil
      USE vaspxml
      USE ini
      USE constant
      IMPLICIT NONE
      TYPE (ham_handle) HAMILTONIAN
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (potcar) P(:)
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (wavespin) WXI
      TYPE (wavespin) W1
      INTEGER LMDIM,BDIR
      OVERLAP CDIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ), &
              CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),W0%WDES%NCDIJ)
      LOGICAL LSOSHIFT
      ! local variables
      TYPE (wavespin) W0_tmp

      COMPLEX(q), POINTER :: SW0(:,:,:,:)=>NULL()

      INTEGER N
      REAL(q) TOTENL,TOTEN,DE

      INTEGER ICOUEV,IERROR
      REAL(q) DESUM,RMS(INFO%NELM),CSHIFT
      REAL(q) EBREAK_STORE
      LOGICAL LRESET
      REAL(q) FBREAK_MEAN

      !< calculate L_R Q_R | u_nk >
      !< essentially applies previously constructed Hamiltonian in CLQIJ
      !< (eq. 57 of arXiv:2504.17380 [cond-mat.mtrl-sci])
      CALL MLR_LQ_PSI(W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CLQIJ,INFO,WXI)

      ! calculate PW Zeeman on | u_nk > and add to WXI
      IF (LSOSHIFT) THEN
         ! to add plane wave part of Zeeman perturbation, uncomment either MLR_ZEEMAN_PSI or MLR_ZEEMAN_PSI_SIMPLE
         CALL MLR_ZEEMAN_PSI(W0,GRID,BDIR,INFO,WXI)
         !CALL MLR_ZEEMAN_PSI_SIMPLE(W0,GRID,BDIR,INFO,WXI)
      ENDIF

! unit conversion: perturbation entering the Greens function has same units as with S
      WXI%CW=WXI%CW*AUTOA*AUTOA*RYTOEV*2._q

#ifdef apply_greens_function
      CALL MLR_APPLY_GREENS_FUNC(W0,WXI,W1)
#else
#ifdef provide_w1
      CALL MLR_APPLY_GREENS_FUNC(W0,WXI,W1)
#endif
      ! project onto empty subspace | xi_k > = P_e L_R Q_R | u_nk >
!     CALL MLR_PROJ_IN_EMPTY_SUBSPACE(W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)
      CALL MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS(W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)

      W0_tmp=W0
!     W0_tmp%WDES%NKPTS=W1%WDES%NKPTS
      W0_tmp%WDES=>W1%WDES

      NONL_S%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,WXI)

      EBREAK_STORE=INFO%EBREAK; INFO%EBREAK=INFO%EDIFF*0.25_q
      W0_tmp%OVER_BAND=.FALSE.

      IF (IO%IU0>=0) THEN
         WRITE(IO%IU0,142)
         WRITE(17,142)
      ENDIF
  142 FORMAT('       N       E                     dE             ' &
           ,'d eps       ncg     rms')


  130 FORMAT (5X, //, &
       &'----------------------------------------------------', &
       &'----------------------------------------------------'//)

  140 FORMAT (5X, //, &
       &'----------------------------------------- Iteration ', &
       &'(',I4,')  ---------------------------------------'//)

      TOTEN=0

      CSHIFT=0._q
      LRESET=.TRUE.
#ifdef provide_w1
      NONL_S%NK=-1; CALL PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,W1)
      LRESET=.FALSE.; W1%CW=-W1%CW; W1%CPROJ=-W1%CPROJ
#endif
      FBREAK_MEAN=0.0
      iter: DO N=1,INFO%NELM
         CALL START_TIMING("LOOP")

         IF(IO%IU6>=0) WRITE(IO%IU6,140) N
         CALL XML_TAG("scstep")

         CALL LINEAR_RESPONSE_DIIS(HAMILTONIAN,GRID,P,INFO,LATT_CUR,NONLR_S,NONL_S,W1,WXI,W0_tmp,W0_tmp%WDES, &
        &   LMDIM,CDIJ,CQIJ,RMS(N),DESUM,ICOUEV,SV,CSHIFT,IO%IU6,IO%IU0,LRESET,IERROR,SW0, FBREAK_MEAN)

         CALL MRG_CEL(W1%WDES,W1)

         LRESET=.FALSE.
         ! break criterion?
         W1%FERTOT(:,1:NKPTS_ORIG,:)=W0%FERTOT(:,1:NKPTS_ORIG,:)
         TOTENL=TOTEN
         TOTEN =BANDSTRUCTURE_ENERGY(W1%WDES,W1)
         DE= (TOTEN-TOTENL)

         IF (IO%IU0>=0) THEN
            WRITE(IO%IU0,200) N, TOTEN, DE, DESUM, ICOUEV, RMS(N)
            WRITE(17,200)     N, TOTEN, DE, DESUM, ICOUEV, RMS(N)
         ENDIF
 200     FORMAT('RLQ: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5,I6,'  ',E10.3)
!> ============================================================================
!> total time used for this step
!> ============================================================================
         CALL SEPERATOR_TIMING(IO%IU6)
         CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')

         INFO%LABORT=.FALSE.
         IF(ABS(DESUM)<INFO%EDIFF.AND.ABS(DE)<INFO%EDIFF) INFO%LABORT=.TRUE.
!        IF(N>=3) THEN
!           IF (ABS((RMS(N)-RMS(N-1))/RMS(N))<1E-1 .AND. ABS((RMS(N)-RMS(N-2))/RMS(N))<1E-1) INFO%LABORT=.TRUE.
!        ENDIF

         IF (IO%IU6>=0)  THEN
            WRITE(IO%IU6,210) TOTEN

210   FORMAT(/ &
             '  free energy    TOTEN  = ',F18.8,' eV'/ &
             '  ---------------------------------------------------')

            IF (IO%LOPEN) CALL WFORCE(IO%IU6)
            IF (IO%LOPEN) CALL WFORCE(17)
         ENDIF

!> ============================================================================
!> xml related output
!> ============================================================================
         CALL XML_TAG("energy")
         IF (INFO%LABORT .OR. N==1) THEN
            CALL XML_ENERGY(TOTEN, TOTEN, TOTEN)
         ELSE
            CALL XML_ENERGY(TOTEN, TOTEN, TOTEN)
         ENDIF
         CALL XML_CLOSE_TAG

         CALL XML_CLOSE_TAG("scstep")


         IF (INFO%LABORT) THEN
            IF (IO%IU6>=0)  THEN

               WRITE(IO%IU6,131)
131       FORMAT (5X, //, &
       &  '------------------------ aborting loop because EDIFF', &
       &  ' is reached ----------------------------------------'//)
            ENDIF
            EXIT  iter
         ENDIF

      ENDDO iter

      CALL REDIS_PW_OVER_BANDS(W0_tmp%WDES,W0_tmp)
      DEALLOCATE(SW0); NULLIFY(SW0)

      INFO%EBREAK=EBREAK_STORE

      ! a bit heuristically
      W1%CW=-W1%CW
#endif

      RETURN
      END SUBROUTINE MLR_PSI_RESPONSE_LQ

!> ********************* SUBROUTINE MLR_CALC_CLQIJ ****************************
!>
!> Prepares Hamiltonian for calculation of Q_R components, V^{nl}_{Q_R, R}
!> (eq. 24 arXiv:2504.17380 [cond-mat.mtrl-sci])
!>
!> ****************************************************************************
      SUBROUTINE MLR_CALC_CLQIJ(WDES,P,T_INFO,LMDIM,IDIR,CLQIJ)
      USE prec
      USE pseudo
      USE wave_high
      USE poscar
      USE relativistic
      USE radial
      USE constant
      USE morbitalmag, ONLY : LZORA, LSOSHIFT
      IMPLICIT NONE
      TYPE (wavedes) WDES
      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      INTEGER LMDIM
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      INTEGER IDIR
      ! local variables
      TYPE (potcar), POINTER :: PP
      INTEGER NI,NT,NIP
      INTEGER LYMAX,LDIMP,L,LL,LLP
      INTEGER LM,LMP,M,MP,CH1,CH2
      COMPLEX(q), ALLOCATABLE :: L_OP_R(:,:,:,:),DUMMY(:,:,:)

      INTEGER, EXTERNAL :: MAXL1

      REAL(q), ALLOCATABLE :: POT(:),KQ(:)
      REAL(q) KQION

      CLQIJ = 0._q

      ion: DO NI=1,T_INFO%NIONS
         NT=T_INFO%ITYP(NI)

         ! cycle if this ion is not treated by this node
         NIP=NI_LOCAL(NI,WDES%COMM_INB)
         IF (NIP==0) CYCLE ion

         PP=>PP_POINTER(P, NI, NT)

         IF (LZORA.OR.LSOSHIFT) THEN
            ALLOCATE(POT(PP%R%NMAX),KQ(PP%R%NMAX))
            CALL NMR_ZORA_POT(PP, POT)
         ENDIF

         ! cycle if this ion has no depletion charge
         IF (PP%PSDMAX==0) CYCLE ion
         LYMAX=MAXL1(PP)

         LDIMP=3
         IF (LYMAX>2) LDIMP=LYMAX+1

         ALLOCATE(L_OP_R(2*LDIMP+1,2*LDIMP+1,3,LDIMP),DUMMY(2*LDIMP+1,2*LDIMP+1,4))
         L_OP_R=(0._q,0._q)
         DO L=1,LDIMP
            CALL SETUP_LS(L,0._q,0._q,L_OP_R(1:2*L+1,1:2*L+1,1:3,L),DUMMY(1:2*L+1,1:2*L+1,1:4))
         ENDDO

         LM=1
         DO CH1=1,PP%LMAX
         LMP=1
         DO CH2=1,PP%LMAX
!           IF (PP%NDEP(CH1,CH2)==0) GOTO 510
            LL =PP%LPS(CH1)
            LLP=PP%LPS(CH2)
            IF (LL == LLP .AND. LL>0 .AND. LL<=LDIMP ) THEN
               IF (LZORA.OR.LSOSHIFT) THEN
                  ! ZORA: radially variant relativistic mass correction 1/(1-V(r)/2/c^2) in Hartree atomic units,
                  ! in accordance with Eq. A13 of Yates et al., J. Chem. Phys. 118, 5746 (2003). See NMR_ZORA_K.
                  KQ(:)=NMR_ZORA_K(PP%R%NMAX,POT)*PP%WAE(:,CH1)*PP%WAE(:,CH2)-PP%WPS(:,CH1)*PP%WPS(:,CH2)
                  CALL SIMPI(PP%R,KQ,KQION)
!                 WRITE(*,*) 'QION=',PP%QION(CH1,CH2),'KQION=',KQION
               ELSE
                  KQION=PP%QION(CH1,CH2)
               ENDIF

               DO M =1,2*LL+1
               DO MP=1,2*LLP+1
                  CLQIJ(LMP+MP-1,LM+M-1,NIP,1)=CLQIJ(LMP+MP-1,LM+M-1,NIP,1)+L_OP_R(M,MP,IDIR,LL)*KQION
               ENDDO
               ENDDO
            ENDIF
 510        LMP=LMP+(2*LLP+1)
         ENDDO
         LM=LM+(2*LL+1)
         ENDDO

         DEALLOCATE(L_OP_R,DUMMY)
         IF (LZORA.OR.LSOSHIFT) DEALLOCATE(POT,KQ)
      ENDDO ion

      IF (WDES%ISPIN==2) CLQIJ(:,:,:,2)=CLQIJ(:,:,:,1)
      IF (WDES%NCDIJ==4) CLQIJ(:,:,:,4)=CLQIJ(:,:,:,1)

      RETURN
      END SUBROUTINE MLR_CALC_CLQIJ

!> ********************* SUBROUTINE MLR_CALC_ZEEMAN_ONE ***********************
!>
!> Calculates V^{nl}_{\sigma, R}
!> (eq. 23 arXiv:2504.17380 [cond-mat.mtrl-sci])
!>
!> ****************************************************************************
      SUBROUTINE MLR_CALC_ZEEMAN_ONE(WDES,P,T_INFO,LMDIM,IDIR,CLQIJ)
      USE prec
      USE pseudo
      USE wave_high
      USE poscar
      USE relativistic
      USE radial
      USE constant
      USE morbitalmag, ONLY : LZORA, LSOSHIFT
      IMPLICIT NONE
      TYPE (wavedes) WDES
      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      INTEGER LMDIM
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      INTEGER IDIR
      ! local variables
      TYPE (potcar), POINTER :: PP
      INTEGER NI,NT,NIP,I,J
      INTEGER LYMAX,LDIMP,L,LL,LLP
      INTEGER LM,LMP,M,MP,CH1,CH2
      COMPLEX(q), ALLOCATABLE :: L_OP_R(:,:,:,:),DUMMY(:,:,:)

      INTEGER, EXTERNAL :: MAXL1

      REAL(q), ALLOCATABLE :: POT(:),KQ(:)
      REAL(q) RDUM,KQION
      REAL(q) BDIR_ONE(3)
      REAL(q) BSTRENGTH
      COMPLEX(q) BDOTSIG(4)

      BDIR_ONE=0._q
      BDIR_ONE(IDIR)=1._q

      ion: DO NI=1,T_INFO%NIONS
         NT=T_INFO%ITYP(NI)

         ! cycle if this ion is not treated by this node
         NIP=NI_LOCAL(NI,WDES%COMM_INB)
         IF (NIP==0) CYCLE ion

         PP=>PP_POINTER(P, NI, NT)

         IF (LZORA.OR.LSOSHIFT) THEN
            ALLOCATE(POT(PP%R%NMAX),KQ(PP%R%NMAX))
            CALL NMR_ZORA_POT(PP, POT)
         ENDIF

!        ! cycle if this ion has no depletion charge
!        IF (PP%PSDMAX==0) CYCLE ion
!        LYMAX=MAXL1(PP)

         LDIMP=3                      ! ??
         IF (LYMAX>2) LDIMP=LYMAX+1   ! ??

! Modelled on relativistic.F
! <up| B \cdot \sigma |up>
         BDOTSIG(1)=  (1._q,0._q)*BDIR_ONE(3)
! <up| B \cdot \sigma |down>
         BDOTSIG(2)=  (1._q,0._q)*BDIR_ONE(1) + (0._q,1._q)*BDIR_ONE(2)
! <down| B \cdot \sigma |up>
         BDOTSIG(3)=  (1._q,0._q)*BDIR_ONE(1) - (0._q,1._q)*BDIR_ONE(2)
! <down| B \cdot \sigma |down>
         BDOTSIG(4)= -(1._q,0._q)*BDIR_ONE(3)

         !BSTRENGTH=-1._q/2._q
         !BSTRENGTH = HSQDTM !AUTOA*AUTOA*RYTOEV=HSQDTM !TOM 03-12-2024 not this
         BSTRENGTH = 1._q !TOM 01-11-2023
         BDOTSIG(:)=BDOTSIG(:)*BSTRENGTH

         LM=1
         DO CH1=1,PP%LMAX
         LMP=1
         DO CH2=1,PP%LMAX
!           IF (PP%NDEP(CH1,CH2)==0) GOTO 510
            LL =PP%LPS(CH1)
            LLP=PP%LPS(CH2)
            IF (LL == LLP .AND. LL<=LDIMP ) THEN
               IF (LZORA.OR.LSOSHIFT) THEN
                  ! ZORA: radially variant relativistic mass correction 1/(1-V(r)/2/c^2) in Hartree atomic units,
                  ! in accordance with Eq. A13 of Yates et al., J. Chem. Phys. 118, 5746 (2003).
                  KQ(:)=NMR_ZORA_K(PP%R%NMAX,POT)*PP%WAE(:,CH1)*PP%WAE(:,CH2)-PP%WPS(:,CH1)*PP%WPS(:,CH2)
                  CALL SIMPI(PP%R,KQ,KQION)
               ELSE
                  KQION=PP%QION(CH1,CH2)
               ENDIF

               ! modelled on relativistic.F
               DO I=0,1
               DO J=0,1
                  DO MP=1,2*LLP+1
                     CLQIJ(LMP+MP-1,LM+MP-1,NIP,J+2*I+1)=CLQIJ(LMP+MP-1,LM+MP-1,NIP,J+2*I+1)+KQION*BDOTSIG(I+2*J+1)
                  ENDDO
               ENDDO
               ENDDO
            ENDIF
 510        LMP=LMP+(2*LLP+1)
         ENDDO
         LM=LM+(2*LL+1)
         ENDDO

         IF (LZORA.OR.LSOSHIFT) DEALLOCATE(POT,KQ)
      ENDDO ion

      RETURN
      END SUBROUTINE MLR_CALC_ZEEMAN_ONE

!> ************************** SUBROUTINE MLR_LQ_PSI ***************************
!>
!> Supposedly calculates L_R Q_R | u_nk >
!>
!> ****************************************************************************
      SUBROUTINE MLR_LQ_PSI(W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CLQIJ,INFO,WXI)
      USE prec
      USE base
      USE mgrid
      USE wave_high
      USE nonl_high
      USE lattice
      USE dfast
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      TYPE (wavespin) WXI
      INTEGER LMDIM
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      ! local variables
      TYPE (wavefun1) W0_1(W0%WDES%NSIM)  !< current wavefunction
      TYPE (wavedes1) WDES1

      COMPLEX(q),ALLOCATABLE :: CH(:,:),CWORK1(:,:)

      INTEGER ISP,NK,ISPINOR,M,MM
      INTEGER NSIM,NGVECTOR
      INTEGER NB(W0%WDES%NSIM),NB_DONE,NP,N
      LOGICAL LDO(W0%WDES%NSIM)
      REAL(q) DUMMY(W0%WDES%NSIM)

      LOGICAL LSTOP

      NSIM=W0%WDES%NSIM
      DUMMY=0._q

      ALLOCATE(CH(W0%WDES%NRPLWV,NSIM))

      IF (INFO%LREAL) THEN
         DO NP=1,NSIM
            ALLOCATE(W0_1(NP)%CR(GRID%MPLWV*W0%WDES%NRSPINORS))
         ENDDO
         ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,NSIM))
      ENDIF

      WXI%CPROJ=0

!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,WXI%WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
         CALL SETWDES(W0%WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)

         NGVECTOR=WDES1%NGVECTOR

         IF (INFO%LREAL) THEN
            CALL PHASER (GRID,LATT_CUR,NONLR_S,NK,W0%WDES)
         ELSE
            NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK)
         ENDIF

         NB_DONE=0  !< index of the bands allready optimised
         bands: DO
            NB=0    !< empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
            newband: DO NP=1,NSIM
               IF (NB_DONE < WXI%WDES%NBANDS ) THEN
                  NB_DONE=NB_DONE+1
                  N     =NB_DONE
                  NB(NP)=NB_DONE

                  CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

                  IF (INFO%LREAL) THEN
                     ! fft to real space
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO newband
!> ============================================================================
!> if the NB list is now empty end the bands DO loop
!> ============================================================================
            LSTOP=.TRUE.
            LDO  =.FALSE.
            DO NP=1,NSIM
               IF (NB(NP)/=0) THEN
                  LSTOP  =.FALSE.
                  LDO(NP)=.TRUE.     ! band not finished yet
               ENDIF
            ENDDO
            IF (LSTOP) EXIT bands

            CH=(0._q,0._q)

            ! contribution | p_i,k > LQ_ij < p_j,k | u_nk >
            IF (INFO%LREAL) THEN
               CWORK1=0
               CALL RACCMU_CCDIJ(NONLR_S,WDES1,W0_1,LMDIM,CLQIJ,CLQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     CALL VNLACC_ADD_CCDIJ(NONL_S,W0_1(NP),CLQIJ,CLQIJ,1,DUMMY(NP),CH(:,NP))
                  ENDIF
               ENDDO
            ENDIF

            ! store in WXI
            DO NP=1,NSIM
               N=NB(NP); IF (.NOT.LDO(NP)) CYCLE
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,NGVECTOR
                     MM=M+ISPINOR*NGVECTOR
                     WXI%CW(MM,N,NK,ISP)=CH(MM,NP)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO bands
!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      DEALLOCATE (CH)
      IF (ALLOCATED(CWORK1)) DEALLOCATE(CWORK1)

      RETURN
      END SUBROUTINE MLR_LQ_PSI

!> ************************ SUBROUTINE MLR_ZEEMAN_PSI *************************
!>
!> Calculates the \sigma \cdot B term of H^{(1)}
!>
!> ****************************************************************************
      SUBROUTINE MLR_ZEEMAN_PSI(W0,GRID,BDIR,INFO,WXI)
      USE prec
      USE base
      USE mgrid
      USE wave_high
      USE nonl_high
      USE lattice
      USE dfast
      USE constant
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (info_struct) INFO
      TYPE (wavespin) WXI
      INTEGER BDIR
      ! local variables
      TYPE (wavefun1) W0_1(W0%WDES%NSIM)  ! current wavefunction
      TYPE (wavedes1) WDES1

      COMPLEX(q),ALLOCATABLE :: CH(:,:)

      INTEGER ISP,NK,ISPINOR,M,MM
      INTEGER NSIM,NGVECTOR
      INTEGER NB(W0%WDES%NSIM),NB_DONE,NP,N
      LOGICAL LDO(W0%WDES%NSIM)
      REAL(q) DUMMY(W0%WDES%NSIM)
      REAL(q) BSTRENGTH
      COMPLEX(q) IC

      LOGICAL LSTOP

      NSIM=W0%WDES%NSIM
      DUMMY=0._q

      ALLOCATE(CH(W0%WDES%NRPLWV,NSIM))

!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,WXI%WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
         CALL SETWDES(W0%WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)

         NGVECTOR=WDES1%NGVECTOR

         NB_DONE=0  ! index of the bands allready optimised
         bands: DO
            NB=0    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
            newband: DO NP=1,NSIM
               IF (NB_DONE < WXI%WDES%NBANDS ) THEN
                  NB_DONE=NB_DONE+1
                  N     =NB_DONE
                  NB(NP)=NB_DONE

                  CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

               ENDIF
            ENDDO newband

!> ============================================================================
!> if the NB list is now empty end the bands DO loop
!> ============================================================================
            LSTOP=.TRUE.
            LDO  =.FALSE.
            DO NP=1,NSIM
               IF (NB(NP)/=0) THEN
                  LSTOP  =.FALSE.
                  LDO(NP)=.TRUE.     ! band not finished yet
               ENDIF
            ENDDO
            IF (LSTOP) EXIT bands

            CH=(0._q,0._q)

            IC=CMPLX(0._q,1._q, q)
            !BSTRENGTH=-1._q/2._q
            !BSTRENGTH = HSQDTM !AUTOA*AUTOA*RYTOEV=HSQDTM !TOM 03-12-2024 not this
            BSTRENGTH = 1._q !TOM 01-11-2023
            DO NP=1,NSIM
               IF (LDO(NP)) THEN
                  IF (BDIR.EQ.1) THEN
                     DO M=1,NGVECTOR
                        MM=M+NGVECTOR
                        CH( M,NP)=CH( M,NP)+W0_1(NP)%CW(MM)*BSTRENGTH
                        CH(MM,NP)=CH(MM,NP)+W0_1(NP)%CW( M)*BSTRENGTH
                     ENDDO
                  ELSEIF (BDIR.EQ.2) THEN
                     DO M=1,NGVECTOR
                        MM=M+NGVECTOR
                        CH( M,NP)=CH( M,NP)-W0_1(NP)%CW(MM)*BSTRENGTH*IC
                        CH(MM,NP)=CH(MM,NP)+W0_1(NP)%CW( M)*BSTRENGTH*IC
                     ENDDO
                  ELSEIF (BDIR.EQ.3) THEN
                     DO M=1,NGVECTOR
                        MM=M+NGVECTOR
                        CH( M,NP)=CH( M,NP)+W0_1(NP)%CW( M)*BSTRENGTH
                        CH(MM,NP)=CH(MM,NP)-W0_1(NP)%CW(MM)*BSTRENGTH
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO

            ! add to WXI
            DO NP=1,NSIM
               N=NB(NP); IF (.NOT.LDO(NP)) CYCLE
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,NGVECTOR
                     MM=M+ISPINOR*NGVECTOR
                     WXI%CW(MM,N,NK,ISP)=WXI%CW(MM,N,NK,ISP)+CH(MM,NP)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO bands
!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      DEALLOCATE (CH)

      RETURN
      END SUBROUTINE MLR_ZEEMAN_PSI

!> ************************ SUBROUTINE MLR_COMMUTATOR *************************
!>
!> Commutator for vGv calcs...
!>
!> ****************************************************************************
      SUBROUTINE MLR_COMMUTATOR( &
     &   HAMILTONIAN,W0,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CDIJ,CQIJ,INFO,IDIR,WXI)
      USE prec
      USE base
      USE wave_high
      USE nonl_high
      USE lattice
      USE mgrid
      USE constant
      USE hamil_struct_def
      IMPLICIT NONE
      TYPE (ham_handle) HAMILTONIAN
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_D
      TYPE (nonlr_struct) NONLR_D          ! projector times r (only CRREXP distinct)
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (latt) LATT_CUR
      INTEGER LMDIM
      OVERLAP CDIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ), &
              CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      INTEGER IDIR ! component of velocity operator 1=x, 2=y, 3=z
      TYPE (wavespin) WXI                  ! stores H(1) - epsilon S(1) | phi_0>
      ! local variables
      TYPE (wavedes1)  WDES1               ! descriptor for one k-point
      TYPE (wavefun1)  W0_1(W0%WDES%NSIM)  ! current wavefunction
      TYPE (nonl_struct) NONL_S_tmp
      TYPE (nonl_struct) NONL_D_tmp

      REAL(q), ALLOCATABLE :: GC(:)
      COMPLEX(q), ALLOCATABLE :: CH(:,:),CWORK1(:,:)
      COMPLEX(q), ALLOCATABLE :: NABLACW(:), NABLACR(:)
      GDEF, TARGET, ALLOCATABLE :: CPROJ(:,:)

      INTEGER ISP,NK,M,MM,ISPINOR,ISPINOR_
      INTEGER NSIM,NGVECTOR
      INTEGER NB(W0%WDES%NSIM),NB_DONE,NP,N
      LOGICAL LDO(W0%WDES%NSIM)
      LOGICAL LSTOP

      REAL(q) EVALUE0(W0%WDES%NSIM)
      REAL(q) G1,G2,G3,RINPLW

      NSIM=W0%WDES%NSIM

      ALLOCATE(CH(W0%WDES%NRPLWV,NSIM),GC(W0%WDES%NRPLWV),CPROJ(W0%WDES%NPROD,NSIM))

      IF (INFO%LREAL.OR.ASSOCIATED(HAMILTONIAN%MU)) THEN
         DO NP=1,NSIM
            ALLOCATE(W0_1(NP)%CR(GRID%MPLWV*W0%WDES%NRSPINORS))
         ENDDO
         N=NSIM; IF (ASSOCIATED(HAMILTONIAN%MU)) N=MAX(N,2)
         ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,N))
      ENDIF

      WXI%CPROJ=0

      spin:    DO ISP=1,W0%WDES%ISPIN
! We have wave functions at k and k+q, the first set have the k-point index running from 1,..,NKPTS/2,
! the latter from NKPTS/2+1,2*NKPTS
      kpoints: DO NK=1,NKPTS_ORIG  ! loop over k
#ifdef MPI
         IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         CALL SETWDES(W0%WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)

         NGVECTOR=WDES1%NGVECTOR

         IF (INFO%LREAL) THEN
            CALL PHASER (GRID,LATT_CUR,NONLR_S,NK+NKPTS_ORIG,W0%WDES); NONLR_S%NK=NK
            CALL PHASERR(GRID,LATT_CUR,NONLR_D,NK,W0%WDES,IDIR)
         ELSE
            NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK)
            NONL_D%NK=NONL_S%NK        ! uses the same phasefactor array
            NONL_S_tmp=NONL_S
            NONL_S_tmp%QPROJ=>NONL_S%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
         ENDIF

         NB_DONE=0  ! index of the bands allready optimised
         bands: DO
            NB=0    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
            newband: DO NP=1,NSIM
               IF (NB_DONE < WXI%WDES%NBANDS ) THEN
                  NB_DONE=NB_DONE+1
                  N     =NB_DONE
                  NB(NP)=NB_DONE

                  CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

                  IF (INFO%LREAL.OR.ASSOCIATED(HAMILTONIAN%MU)) THEN
                     ! fft to real space
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
                     ENDDO
                  ENDIF
                  ! W0_1 points to elements of W0, we don't want to destroy them so we
                  ! create some space to hold the projection of the wave functions onto
                  ! the derivative of the projectors w.r.t. k, as used below.
                  W0_1(NP)%CPROJ => CPROJ(:,NP)

                  EVALUE0(NP)=W0_1(NP)%CELEN
               ENDIF
            ENDDO newband
!> ============================================================================
!> if the NB list is now empty end the bands DO loop
!> ============================================================================
            LSTOP=.TRUE.
            LDO  =.FALSE.
            DO NP=1,NSIM
               IF (NB(NP)/=0) THEN
                  LSTOP  =.FALSE.
                  LDO(NP)=.TRUE.     ! band not finished yet
               ENDIF
            ENDDO
            IF (LSTOP) EXIT bands

            !  -i\nabla + k | u_nk > = hbar^2/ m_e (G+k) C_G | u_nk >
            DO M=1,NGVECTOR
               G1=WDES1%IGX(M)+W0%WDES%VKPT(1,NK)
               G2=WDES1%IGY(M)+W0%WDES%VKPT(2,NK)
               G3=WDES1%IGZ(M)+W0%WDES%VKPT(3,NK)
               GC(M)=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))*(TPI*(2*HSQDTM))
            ENDDO

            CH=(0._q,0._q)

            DO NP=1,NSIM
               IF (LDO(NP)) THEN
                  DO ISPINOR =0,WDES1%NRSPINORS-1
                     DO M=1,NGVECTOR
                        MM=M+ISPINOR*NGVECTOR
                        ! mind the factor -i here, it is removed lateron
                        CH(MM,NP)=CH(MM,NP)+(W0_1(NP)%CW(MM)*GC(M))*(0._q,-1._q)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO

            ! add term [- nabla . mu(r) nabla , r] = - mu(r) nabla - nabla mu(r)
            IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
               ALLOCATE(NABLACR(GRID%MPLWV))
               ALLOCATE(NABLACW(GRID%MPLWV))
               RINPLW=1._q/GRID%NPLWV
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        ! apply nabla | psi >
                        DO M=1,NGVECTOR
                           MM=M+ISPINOR*NGVECTOR
                           NABLACW(MM)=W0_1(NP)%CW(MM)*GC(M)*(0._q,-1._q)/(2*HSQDTM)
                        ENDDO
                        ! fft to real space
                        CALL FFTWAV(NGVECTOR,WDES1%NINDPW,NABLACR,NABLACW,GRID)
                        ! apply mu to ( | psi > + nabla | psi > )
                        CWORK1 = 0
                        DO ISPINOR_=0,WDES1%NRSPINORS-1
                           DO M=1,GRID%RL%NP
                              MM=M+ISPINOR*GRID%MPLWV
                              ! apply mu to | psi >
                              CWORK1(M,1)=CWORK1(M,1)+RINPLW*HAMILTONIAN%MU(M,1+ISPINOR_+2*ISPINOR)*W0_1(NP)%CR(MM)
                              ! apply mu to nabla | psi>
                              CWORK1(M,2)=CWORK1(M,2)+RINPLW*HAMILTONIAN%MU(M,1+ISPINOR_+2*ISPINOR)*NABLACR(M)
                           ENDDO
                        ENDDO
                        ! fft mu |psi> to reciprocal space
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1,1),NABLACW(1),GRID,.FALSE.)
                        ! fft mu nabla | psi > to reciprocal space and add to CH
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1,2),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
                        ! apply nabla  to mu |psi>
                        DO M=1,NGVECTOR
                           MM=M+ISPINOR*NGVECTOR
                           ! nabla mu | psi >
                           CH(MM,NP)=CH(MM,NP)+GC(M)*(0._q,-1._q)*NABLACW(M)/(2*HSQDTM)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               DEALLOCATE(NABLACW)
               DEALLOCATE(NABLACR)
            ENDIF

            ! contribution   | p_i,k+q > [D_ij - e(0)Q_ij] < p_j,k|  r_idir | u_nk >
            IF (INFO%LREAL) THEN
               CWORK1=0

               CALL RPROMU(NONLR_D,WDES1,W0_1,NSIM,LDO)
               CALL RACCMU(NONLR_S,WDES1,W0_1,LMDIM,CDIJ,CQIJ,EVALUE0,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     CALL PROJ1(NONL_D,WDES1,W0_1(NP))
                     CALL VNLACC_ADD(NONL_S_tmp,W0_1(NP),CDIJ,CQIJ,1,EVALUE0(NP),CH(:,NP))
                  ENDIF
               ENDDO
            ENDIF

            ! store in WXI
            DO NP=1,NSIM
               N=NB(NP); IF (.NOT.LDO(NP)) CYCLE
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,NGVECTOR
                     MM=M+ISPINOR*NGVECTOR
                     WXI%CW(MM,N,NK,ISP)=CH(MM,NP)*(0._q,1._q)
                  ENDDO
                  ! we store  i< p_j,k | r_idir | u_nk >  in WXI%CPROJ
                  WXI%CPROJ(:,N,NK,ISP)=W0_1(NP)%CPROJ*(0._q,1._q)
               ENDDO
            ENDDO
         ENDDO bands

         IF (INFO%LREAL) THEN
            CALL PHASER (GRID,LATT_CUR,NONLR_S,NK,W0%WDES)
            CALL PHASERR(GRID,LATT_CUR,NONLR_D,NK+NKPTS_ORIG,W0%WDES,IDIR); NONLR_D%NK=NK
         ELSE
! Gilles and Martijn think this does not need to be recomputed
!           CALL PHASE(WDES,NONL_S,NK)
!           NONL_D%NK=NONL_S%NK        ! uses the same phasefactor array
            NONL_D_tmp=NONL_D
            NONL_D_tmp%QPROJ=>NONL_D%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
         ENDIF


         NB_DONE=0  ! index of the bands allready optimised
         bands2: DO
            NB=0    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
            newband2: DO NP=1,NSIM
               IF (NB_DONE < WXI%WDES%NBANDS ) THEN
                  NB_DONE=NB_DONE+1
                  N     =NB_DONE
                  NB(NP)=NB_DONE

                  CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

                  IF (INFO%LREAL) THEN
                     ! fft to real space
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
                     ENDDO
                  ENDIF

                  EVALUE0(NP)=W0_1(NP)%CELEN
               ENDIF
            ENDDO newband2
!> ============================================================================
!> if the NB list is now empty end the bands DO loop
!> ============================================================================
            LSTOP=.TRUE.
            LDO  =.FALSE.
            DO NP=1,NSIM
               IF (NB(NP)/=0) THEN
                  LSTOP  =.FALSE.
                  LDO(NP)=.TRUE.     ! band not finished yet
               ENDIF
            ENDDO
            IF (LSTOP) EXIT bands2

            CH=(0._q,0._q)

            ! contribution - r_idir | p_i,k+q > [D_ij - e(0)Q_ij] < p_j,k | u_nk >
            IF (INFO%LREAL) THEN
               CWORK1=0
               CALL RACCMU(NONLR_D,WDES1,W0_1,LMDIM,CDIJ,CQIJ,EVALUE0,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     CALL VNLACC_ADD(NONL_D_tmp,W0_1(NP),CDIJ,CQIJ,1,EVALUE0(NP),CH(:,NP))
                  ENDIF
               ENDDO
            ENDIF

            ! store in WXI
            DO NP=1,NSIM
               N=NB(NP); IF (.NOT.LDO(NP)) CYCLE
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,NGVECTOR
                     MM=M+ISPINOR*NGVECTOR
                     WXI%CW(MM,N,NK,ISP)=WXI%CW(MM,N,NK,ISP)-CH(MM,NP)*(0._q,1._q)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO bands2
      ENDDO kpoints
      ENDDO spin

      IF (INFO%LREAL.OR.ASSOCIATED(HAMILTONIAN%MU)) THEN
         DO NP=1,NSIM
            DEALLOCATE(W0_1(NP)%CR)
         ENDDO
         DEALLOCATE(CWORK1)
      ENDIF

      DEALLOCATE(CH,GC,CPROJ)

      RETURN
      END SUBROUTINE MLR_COMMUTATOR

!> *********************** SUBROUTINE MLR_COMMUTATOR_RQ ***********************
!>
!> Calculate commutator for MLR_PSI_RESPONSE_OCC
!>
!> ****************************************************************************
      SUBROUTINE MLR_COMMUTATOR_RQ( &
     &   W0,GRID,NONL_S,NONLR_S,NONL_D,NONLR_D,LATT_CUR,LMDIM,CQIJ,INFO,IDIR,WXI)
      USE prec
      USE base
      USE wave_high
      USE nonl_high
      USE lattice
      USE mgrid
      USE constant
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (grid_3d) GRID
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_D
      TYPE (nonlr_struct) NONLR_D          ! projector times r (only CRREXP distinct)
      TYPE (info_struct) INFO
      TYPE (in_struct) IO
      TYPE (latt) LATT_CUR
      INTEGER LMDIM
      OVERLAP CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      INTEGER IDIR ! component of velocity operator 1=x, 2=y, 3=z
      TYPE (wavespin) WXI                  ! S(1) | phi_0>
      ! local variables
      TYPE (wavedes1)  WDES1               ! descriptor for one k-point
      TYPE (wavefun1)  W0_1(W0%WDES%NSIM)  ! current wavefunction
      TYPE (nonl_struct) NONL_S_tmp
      TYPE (nonl_struct) NONL_D_tmp

      COMPLEX(q),ALLOCATABLE :: CH(:,:),CWORK1(:,:)
      GDEF, TARGET, ALLOCATABLE ::  CPROJ(:,:)

      INTEGER ISP,NK,M,MM,ISPINOR
      INTEGER NSIM,NGVECTOR
      INTEGER NB(W0%WDES%NSIM),NB_DONE,NP,N
      LOGICAL LDO(W0%WDES%NSIM)
      LOGICAL LSTOP

      REAL(q) DUMMY(W0%WDES%NSIM)

      NSIM=W0%WDES%NSIM
      DUMMY=0._q

      ALLOCATE(CH(W0%WDES%NRPLWV,NSIM),CPROJ(W0%WDES%NPROD,NSIM))

      IF (INFO%LREAL) THEN
         DO NP=1,NSIM
            ALLOCATE(W0_1(NP)%CR(GRID%MPLWV*W0%WDES%NRSPINORS))
         ENDDO
         ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,NSIM))
      ENDIF

      WXI%CPROJ=0

      spin:    DO ISP=1,W0%WDES%ISPIN
! We wave functions at k and k+q, the first set have the k-point index running from 1,..,NKPTS/2,
! the latter from NKPTS/2+1,2*NKPTS
      kpoints: DO NK=1,NKPTS_ORIG  ! loop over k
#ifdef MPI
         IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         CALL SETWDES(W0%WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)

         NGVECTOR=WDES1%NGVECTOR

         IF (INFO%LREAL) THEN
            CALL PHASER (GRID,LATT_CUR,NONLR_S,NK+NKPTS_ORIG,W0%WDES); NONLR_S%NK=NK
            CALL PHASERR(GRID,LATT_CUR,NONLR_D,NK,W0%WDES,IDIR)
         ELSE
            NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK)
            NONL_D%NK=NONL_S%NK        ! uses the same phasefactor array
            NONL_S_tmp=NONL_S
            NONL_S_tmp%QPROJ=>NONL_S%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
         ENDIF

         NB_DONE=0  ! index of the bands allready optimised
         bands: DO
            NB=0    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
            newband: DO NP=1,NSIM
               IF (NB_DONE < WXI%WDES%NBANDS ) THEN
                  NB_DONE=NB_DONE+1
                  N     =NB_DONE
                  NB(NP)=NB_DONE

                  CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

                  IF (INFO%LREAL) THEN
                     ! fft to real space
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
                     ENDDO
                  ENDIF
                  ! W0_1 points to elements of W0, we don't want to destroy them so we
                  ! create some space to hold the projection of the wave functions onto
                  ! the derivative of the projectors w.r.t. k, as used below.
                  W0_1(NP)%CPROJ => CPROJ(:,NP)
               ENDIF
            ENDDO newband
!> ============================================================================
!>  if the NB list is now empty end the bands DO loop
!> ============================================================================
            LSTOP=.TRUE.
            LDO  =.FALSE.
            DO NP=1,NSIM
               IF (NB(NP)/=0) THEN
                  LSTOP  =.FALSE.
                  LDO(NP)=.TRUE.     ! band not finished yet
               ENDIF
            ENDDO
            IF (LSTOP) EXIT bands

            CH=(0._q,0._q)

            ! contribution   | p_i,k+q > Q_ij < p_j,k|  r_idir | u_nk >
            IF (INFO%LREAL) THEN
               CWORK1=0
               CALL RPROMU(NONLR_D,WDES1,W0_1,NSIM,LDO)
               CALL RACCMU(NONLR_S,WDES1,W0_1,LMDIM,CQIJ,CQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     CALL PROJ1(NONL_D,WDES1,W0_1(NP))
                     CALL VNLACC_ADD(NONL_S_tmp,W0_1(NP),CQIJ,CQIJ,1,DUMMY(NP),CH(:,NP))
                  ENDIF
               ENDDO
            ENDIF

            ! store in WXI
            DO NP=1,NSIM
               N=NB(NP); IF (.NOT.LDO(NP)) CYCLE
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,NGVECTOR
                     MM=M+ISPINOR*NGVECTOR
                     WXI%CW(MM,N,NK,ISP)=CH(MM,NP)*(0._q,1._q)
                  ENDDO
                  ! we store  i< p_j,k | r_idir | u_nk >  in WXI%CPROJ
                  WXI%CPROJ(:,N,NK,ISP)=W0_1(NP)%CPROJ*(0._q,1._q)
               ENDDO
            ENDDO
         ENDDO bands

         IF (INFO%LREAL) THEN
            CALL PHASER (GRID,LATT_CUR,NONLR_S,NK,W0%WDES)
            CALL PHASERR(GRID,LATT_CUR,NONLR_D,NK+NKPTS_ORIG,W0%WDES,IDIR); NONLR_D%NK=NK
         ELSE
! Gilles and Martijn think this does not need to be recomputed
!           CALL PHASE(WDES,NONL_S,NK)
!           NONL_D%NK=NONL_S%NK        ! uses the same phasefactor array
            NONL_D_tmp=NONL_D
            NONL_D_tmp%QPROJ=>NONL_D%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
         ENDIF


         NB_DONE=0  ! index of the bands allready optimised
         bands2: DO
            NB=0    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
            newband2: DO NP=1,NSIM
               IF (NB_DONE < WXI%WDES%NBANDS ) THEN
                  NB_DONE=NB_DONE+1
                  N     =NB_DONE
                  NB(NP)=NB_DONE

                  CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

                  IF (INFO%LREAL) THEN
                     ! fft to real space
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO newband2
!> ============================================================================
!> if the NB list is now empty end the bands DO loop
!> ============================================================================
            LSTOP=.TRUE.
            LDO  =.FALSE.
            DO NP=1,NSIM
               IF (NB(NP)/=0) THEN
                  LSTOP  =.FALSE.
                  LDO(NP)=.TRUE.     ! band not finished yet
               ENDIF
            ENDDO
            IF (LSTOP) EXIT bands2

            CH=(0._q,0._q)

            ! contribution - r_idir | p_i,k+q > Q_ij < p_j,k | u_nk >
            IF (INFO%LREAL) THEN
               CWORK1=0
               CALL RACCMU(NONLR_D,WDES1,W0_1,LMDIM,CQIJ,CQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     DO ISPINOR=0,WDES1%NRSPINORS-1
                        CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO NP=1,NSIM
                  IF (LDO(NP)) THEN
                     CALL VNLACC_ADD(NONL_D_tmp,W0_1(NP),CQIJ,CQIJ,1,DUMMY(NP),CH(:,NP))
                  ENDIF
               ENDDO
            ENDIF

            ! store in WXI
            DO NP=1,NSIM
               N=NB(NP); IF (.NOT.LDO(NP)) CYCLE
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,NGVECTOR
                     MM=M+ISPINOR*NGVECTOR
                     WXI%CW(MM,N,NK,ISP)=WXI%CW(MM,N,NK,ISP)-CH(MM,NP)*(0._q,1._q)
                  ENDDO
               ENDDO
            ENDDO

         ENDDO bands2
      ENDDO kpoints
      ENDDO spin
! test
!     WXI%CW=WXI%CW*AUTOA
      WXI%CW=WXI%CW/2._q
! test

      IF (INFO%LREAL) THEN
         DO NP=1,NSIM
            DEALLOCATE(W0_1(NP)%CR)
         ENDDO
         DEALLOCATE(CWORK1)
      ENDIF

      DEALLOCATE(CH,CPROJ)

      RETURN
      END SUBROUTINE MLR_COMMUTATOR_RQ

!> ************************ SUBROUTINE MLR_COMMUTATOR_FOCK ********************
!>
!> ****************************************************************************
      SUBROUTINE MLR_COMMUTATOR_FOCK( &
         W0, GRID, P, NONL_S, NONLR_S, NONL_D, NONLR_D, LATT_CUR, LMDIM, IDIR, WXI)

      USE prec
      USE wave_high
      USE nonl_high
      USE lattice
      USE pseudo
      USE mgrid
      USE fock_glb, ONLY : LHFCALC, AEXX, WDES_FOCK, GRID_FOCK, GRIDHF, &
                           AUG_DES, FAST_AUG_FOCK, TRANS_MATRIX_FOCK, &
                           HFKIDENT, FSG_STORE, NBANDSGWLOW_FOCK
      USE fock, ONLY : SKIP_THIS_KPOINT_IN_FOCK, SET_GFAC, SET_GFAC_QDER,&
                       FOCK_CHARGE, FOCK_CHARGE_AUG
      USE augfast, ONLY : CALC_DLLMM_TRANS
      USE mkpoints_struct_def, ONLY : KPOINTS_FULL, ROTATION_HANDLE
      USE full_kpoints, ONLY : CHECK_FULL_KPOINTS, DEALLOCATE_ROT_HANDLE
#ifndef gammareal
      USE full_kpoints, ONLY : W1_ROTATE_AND_FFT
#endif
      USE dfast, ONLY: NSTRIP_STANDARD
      USE sym_prec, ONLY : TINY

      TYPE (wavespin) :: W0
      TYPE (grid_3d) :: GRID
      TYPE (nonl_struct) :: NONL_S
      TYPE (nonlr_struct) :: NONLR_S
      TYPE (nonl_struct) :: NONL_D
      TYPE (nonlr_struct) :: NONLR_D
      TYPE (potcar) :: P(:)
      TYPE (latt) :: LATT_CUR
      TYPE (wavespin) :: WXI

      INTEGER :: LMDIM
      INTEGER :: IDIR

      ! local variables
      TYPE (wavespin) :: WHF

      TYPE( rotation_handle), POINTER :: ROT_HANDLE

      TYPE (wavedes1) WDESK0, WDESK, WDESQ, WDESQ_IRZ

      TYPE (wavefun1), ALLOCATABLE :: WIN(:), WTMP(:)
      TYPE (wavefun1) :: W1, WQ

      TYPE (nonl_struct) :: NONL_S_KPQ, NONL_D_KPQ
      TYPE (nonlr_struct) :: NONLR_S_KPQ, NONLR_D_KPQ

      COMPLEX(q), ALLOCATABLE :: CXI(:,:), CH(:)
      GDEF, ALLOCATABLE :: CKAPPA1(:,:), CKAPPA2(:,:)
      GDEF, ALLOCATABLE :: CRHOLM(:), CDIJ(:,:,:,:)
      GDEF, ALLOCATABLE, TARGET :: CDLM(:)
      REAL(q), ALLOCATABLE :: POTFAK(:), POTFAKQ(:)

      GDEF :: GWORK (m_ GRIDHF%MPLWV,W0%WDES%NRSPINORS)
      GDEF :: GWORK2(m_ GRIDHF%MPLWV,W0%WDES%NRSPINORS)

      COMPLEX(q) :: CWORK(W0%WDES%GRID%MPLWV*W0%WDES%NRSPINORS)

      LOGICAL, ALLOCATABLE :: LDO(:)

      REAL(q) :: FSG, WEIGHT

      INTEGER :: ISP, ISP_IRZ, NK, NPOS, NSTRIP, NSTRIP_ACT, NGLB, NGLB_ACT, N, NQ, MQ, NP, ISPINOR

      LOGICAL :: LSHIFT

      IF ((.NOT.LHFCALC) .OR. (AEXX==0)) RETURN

      WHF=W0 ; WHF%WDES => WDES_FOCK

      CALL CHECK_FULL_KPOINTS

      NULLIFY(ROT_HANDLE)

      NSTRIP= MIN(NSTRIP_STANDARD, W0%WDES%NBANDS)
      NGLB = NSTRIP*W0%WDES%NB_PAR

      ALLOCATE( &
         CXI(GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS,NGLB), &
         CKAPPA1(WHF%WDES%NPROD,NGLB), &
         CKAPPA2(WHF%WDES%NPROD,NGLB), &
         CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS), &
         CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS), &
         POTFAK(GRIDHF%MPLWV), &
         POTFAKQ(GRIDHF%MPLWV), &
         CH(W0%WDES%NRPLWV), &
         LDO(NGLB))

      CALL SETWDES(WHF%WDES, WDESK, 0)
      ALLOCATE(WIN(NGLB))
      DO N=1,NGLB
         CALL NEWWAV(WIN(N), WDESK, .TRUE.)
      ENDDO

      CALL SETWDES(W0%WDES, WDESK0, 0)
      ALLOCATE(WTMP(NGLB))
      DO N=1,NGLB
!!         CALL NEWWAV(WTMP(N), WDESK, .TRUE.)
         CALL NEWWAV_NOCW(WTMP(N), WDESK0)
         WTMP(N)%CW => WIN(N)%CW
      ENDDO


      CALL SETWDES(WHF%WDES, WDESQ, 0)
      CALL NEWWAV(WQ, WDESQ, .TRUE.)

      IF (NONLR_S%LREAL) THEN
         CALL NONLR_ASSIGN(NONLR_S_KPQ, NONLR_S)
         CALL NONLR_ALLOC_CRREXP(NONLR_S_KPQ)

         CALL NONLR_ASSIGN(NONLR_D_KPQ, NONLR_S)
         CALL NONLR_ALLOC_CRREXP(NONLR_D_KPQ)
      ENDIF

      spin: DO ISP = 1, W0%WDES%ISPIN
         kpoints: DO NK = 1, NKPTS_ORIG
#ifdef MPI
            IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE kpoints
#endif
            CALL SETWDES(WHF%WDES, WDESK, NK)
            CALL SETWDES(W0%WDES, WDESK0, NK)

            IF (NONLR_S%LREAL) THEN
               CALL PHASER (GRID, LATT_CUR, NONLR_S, NK, W0%WDES)
               CALL PHASERR(GRID, LATT_CUR, NONLR_D, NK, W0%WDES, IDIR)

               CALL PHASER (GRID, LATT_CUR, NONLR_S_KPQ, NK+NKPTS_ORIG, W0%WDES) ; NONLR_S_KPQ%NK = NK
               CALL PHASERR(GRID, LATT_CUR, NONLR_D_KPQ, NK+NKPTS_ORIG, W0%WDES, IDIR) ; NONLR_D_KPQ%NK = NK
            ELSE
               CALL PHASE(W0%WDES, NONL_S, NK)
               NONL_D%NK = NONL_S%NK

               NONL_S_KPQ = NONL_S ; NONL_S_KPQ%QPROJ => NONL_S%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
               NONL_D_KPQ = NONL_D ; NONL_D_KPQ%QPROJ => NONL_D%QPROJ(:,:,:,NKPTS_ORIG+1:W0%WDES%NKPTS,:)
            ENDIF

            FSG=FSG_STORE(NK)

            strips: DO NPOS=1,W0%WDES%NBANDS,NSTRIP

               NSTRIP_ACT=MIN(W0%WDES%NBANDS+1-NPOS,NSTRIP)
               NGLB_ACT = NSTRIP_ACT*W0%WDES%NB_PAR

               CDIJ=0; CDLM=0; CXI=0; CKAPPA1=0; CKAPPA2=0

               CALL W1_GATHER(WHF, NPOS, NPOS+NSTRIP_ACT-1, ISP, WIN)

               IF (NONLR_S%LREAL) THEN
!!                  CALL W1_GATHER(W0, NPOS, NPOS+NSTRIP_ACT-1, ISP, WTMP)
                  DO N = 1, NGLB_ACT
                     CALL FFTWAV_W1(WTMP(N))
                  ENDDO
                  LDO = .TRUE. ; CALL RPROMU(NONLR_D, WDESK0, WTMP, NGLB_ACT, LDO)
               ELSE
                  DO N = 1, NGLB_ACT
                     CALL PROJ1(NONL_D, WDESK0, WTMP(N))
                  ENDDO
               ENDIF

               qpoints: DO NQ=1,KPOINTS_FULL%NKPTS

                  IF((KPOINTS_FULL%WTKPT(NQ)==0) .OR. &
                     (HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
                     (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) CYCLE qpoints

                  CALL SETWDES(WHF%WDES, WDESQ, NQ)
                  CALL SETWDES(WHF%WDES, WDESQ_IRZ, KPOINTS_FULL%NEQUIV(NQ))

                  ISP_IRZ = ISP
                  IF (KPOINTS_FULL%SPINFLIP(NQ)==1) ISP_IRZ = 3-ISP

                  CALL SET_GFAC(GRIDHF, LATT_CUR, NK, NQ, FSG, POTFAK)
                  CALL SET_GFAC_QDER(GRIDHF, LATT_CUR, NK, NQ, 0.0_q, POTFAKQ, IDIR)

                  mband: DO MQ = 1, WHF%WDES%NBANDS

                     IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10_q) CYCLE mband
                     IF ((MQ-1)*W0%WDES%NB_PAR+W0%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband

                     IF (NQ<=WHF%WDES%NKPTS) THEN
                        CALL W1_COPY(ELEMENT(WHF, WDESQ, MQ, ISP), WQ)
                        CALL FFTWAV_W1(WQ)
                     ELSE
#ifndef gammareal
                        LSHIFT=.FALSE.
                        IF ((ABS(KPOINTS_FULL%TRANS(1,NQ))>TINY) .OR. &
                            (ABS(KPOINTS_FULL%TRANS(2,NQ))>TINY) .OR. &
                            (ABS(KPOINTS_FULL%TRANS(3,NQ))>TINY)) LSHIFT=.TRUE.
                        CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)
#endif
                     ENDIF

                     nband: DO N = 1, NGLB_ACT

                        IF (WHF%WDES%LOVERL) THEN
                           CALL FOCK_CHARGE(WIN(N), WQ, GWORK(:,1), CRHOLM)

                           CALL FFT3D(GWORK(1,1), GRIDHF, -1)
                           GWORK2 = GWORK

                           CALL APPLY_GFAC(GRIDHF, GWORK(1,1), POTFAK(1))

                           CALL FFT3D(GWORK(1,1), GRIDHF, 1)

                           WEIGHT = WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV

                           W1%CPROJ => CDLM(:) ; AUG_DES%RINPL = WEIGHT
                           CALL RPRO1_HF(FAST_AUG_FOCK, AUG_DES, W1, GWORK(:,1))

                           IF (WHF%WDES%NRSPINORS==2) CDLM(AUG_DES%NPRO+1:AUG_DES%NPRO*2) = CDLM(1:AUG_DES%NPRO)
                           CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ, CDLM)

                           CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), CKAPPA1(1,N), .TRUE.)
                        ELSE
                           GWORK2 = 0
                        ENDIF

                        GWORK = 0
                        CALL FOCK_CHARGE_AUG( WTMP(N), WQ, GWORK(:,1), CRHOLM)

                        CALL FFT3D(GWORK(1,1), GRIDHF, -1)

                        CALL APPLY_GFAC(GRIDHF, GWORK(1,1), POTFAK(1))

                        CALL APPLY_GFAC_QDER(GRIDHF, GWORK2(1,1),  GWORK(1,1), POTFAKQ(1))

                        CALL FFT3D(GWORK(1,1), GRIDHF, 1)

                        WEIGHT = WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV

                        CALL VHAMIL_TRACE(WDESK, GRID_FOCK, GWORK(1,1), WQ%CR(1), CXI(1,N), WEIGHT)

                        IF (WHF%WDES%LOVERL) THEN
                           W1%CPROJ => CDLM(:) ; AUG_DES%RINPL = WEIGHT
                           CALL RPRO1_HF(FAST_AUG_FOCK,AUG_DES, W1, GWORK(:,1))

                           IF (WHF%WDES%NRSPINORS==2) CDLM(AUG_DES%NPRO+1:AUG_DES%NPRO*2) = CDLM(1:AUG_DES%NPRO)
                           CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ, CDLM)

                           CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), CKAPPA2(1,N), .TRUE.)
                        ENDIF

                     ENDDO nband
                  ENDDO mband
               ENDDO qpoints

               IF (WHF%WDES%DO_REDIS) THEN
                  CALLMPI( M_sum_z(WDESK%COMM_INTER, CXI(1,1), GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS*NGLB) )
                  CALLMPI( M_sum_g(WDESK%COMM_INTER, CKAPPA1(1,1), WHF%WDES%NPROD*NGLB) )
                  CALLMPI( M_sum_g(WDESK%COMM_INTER, CKAPPA2(1,1), WHF%WDES%NPROD*NGLB) )
               ENDIF


               fft_back: DO N = W0%WDES%NB_LOW, NGLB_ACT, W0%WDES%NB_PAR
                  NP = (N-1)/W0%WDES%NB_PAR+NPOS

                  CH = 0

                  IF (WHF%WDES%LOVERL) THEN
                     IF (NONLR_S%LREAL) THEN

                        CWORK = 0
                        CALL RACC0(NONLR_D_KPQ, WDESK0, CKAPPA1(1,N), CWORK(1))

                        CWORK=-CWORK
                        CALL RACC0(NONLR_S_KPQ, WDESK0, CKAPPA2(1,N), CWORK(1))

                        DO ISPINOR = 0, WDESK0%NRSPINORS-1
                           CALL FFTEXT(WDESK0%NGVECTOR, WDESK0%NINDPW(1), CWORK(1+ISPINOR*WDESK0%GRID%MPLWV), &
                                       CH(1+ISPINOR*WDESK0%NGVECTOR), WDESK0%GRID, .TRUE.)
                        ENDDO
                     ELSE
                        CALL VNLAC0(NONL_D_KPQ, WDESK, CKAPPA1(1,N), CH(1))

                        CH = -CH
                        CALL VNLAC0(NONL_S_KPQ, WDESK, CKAPPA2(1,N), CH(1))
                     ENDIF
                  ENDIF

                  DO ISPINOR = 0, WDESK%NRSPINORS-1
                     CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW(1), CXI(1+ISPINOR*GRID_FOCK%MPLWV,N), &
                                 CH(1+ISPINOR*WDESK%NGVECTOR), GRID_FOCK, .TRUE.)
                  ENDDO

!!                  WXI%CW(:,NP,NK,ISP) = WXI%CW(:,NP,NK,ISP) + CH(:)*(0._q,1._q)
                  WXI%CW(:,NP,NK,ISP) = WXI%CW(:,NP,NK,ISP) - CH(:)*(0._q,1._q)

               ENDDO fft_back

            ENDDO strips

         ENDDO kpoints
      ENDDO spin

      DEALLOCATE(CXI,CKAPPA1,CKAPPA2,CRHOLM,CDIJ,CDLM,POTFAK,POTFAKQ,CH,LDO)

      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)

      IF (NONLR_S%LREAL) THEN
         CALL NONLR_DEALLOC_CRREXP(NONLR_S_KPQ)
         CALL NONLR_DEALLOC_CRREXP(NONLR_D_KPQ)
      ENDIF

      CALL DELWAV(WQ, .TRUE.)

      DO N = 1, NGLB
         CALL DELWAV(WIN(N), .TRUE.)
!!         CALL DELWAV(WTMP(N), .TRUE.)
         CALL DELWAV_NOCW(WTMP(N))
      ENDDO
      DEALLOCATE(WIN,WTMP)

      RETURN
      END SUBROUTINE MLR_COMMUTATOR_FOCK

!> ********************** SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE ***************
!>
!> ****************************************************************************
      SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE(&
     & W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)
      USE prec
      USE base
      USE dfast
      USE lattice
      USE mgrid
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (wavespin) WXI
      TYPE (grid_3d) GRID
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      INTEGER LMDIM
      OVERLAP CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      ! local variables
      TYPE (wavedes1) WDES1
      TYPE (wavefun1), ALLOCATABLE :: WTMP(:)

      INTEGER ISP,NK,NPOS,N,NB,M,MM,ISPINOR
      INTEGER NB_TOT,NBANDS,NB_START,NB_STOP,NSTRIP,NSTRIP_ACT

      REAL(q), ALLOCATABLE :: DUMMY(:)
      LOGICAL, ALLOCATABLE :: LDO(:)

      GDEF COVL
      COMPLEX(q), ALLOCATABLE :: CH(:,:),CWORK1(:,:)

      NB_TOT=W0%WDES%NB_TOT

      ! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

      CALL SETWDES(W0%WDES,WDES1,0)
      ALLOCATE(WTMP(NSTRIP))
      DO N=1,NSTRIP
         CALL NEWWAV(WTMP(N),WDES1,INFO%LREAL)
      ENDDO

      ALLOCATE(DUMMY(NSTRIP),LDO(NSTRIP))

      ALLOCATE(CH(W0%WDES%NRPLWV,NSTRIP))
      IF (INFO%LREAL) THEN
         ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,NSTRIP))
      ENDIF

!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,NKPTS_ORIG
#ifdef MPI
      IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
      CALL SETWDES(W0%WDES,WDES1,NK+NKPTS_ORIG)

!> ----------------------------------------------------------------------------
!> calculate
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!> |phi(1)_i> - \sum_j S|phi(0)_j><phi(0)_j|phi(1)_i>
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> where j runs over occupied orbitals
!> ----------------------------------------------------------------------------

      IF (INFO%LREAL) THEN
         CALL PHASER(GRID,LATT_CUR,NONLR_S,NK+NKPTS_ORIG,W0%WDES)
      ELSE
         NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK+NKPTS_ORIG)
      ENDIF

      DO NB_START=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NSTRIP,NB_TOT-NB_START+1)
         NB_STOP=MIN(NB_START+NSTRIP-1,NB_TOT)
         ! collect the wave functions from NB_START to NB_STOP (global band index)
         ! and store into WTMP (for LREAL=.TRUE. WTMP(:)%CR is calculated as well).
         IF (INFO%LREAL) THEN
            CALL W1_GATHER_GLB(W0,NB_START,NB_STOP,ISP,WTMP)
         ELSE
            CALL W1_GATHER_GLB_NOCR(W0,NB_START,NB_STOP,ISP,WTMP)
         ENDIF

         CH=(0._q,0._q)
         DO N=1,NSTRIP_ACT
            CH(1:WDES1%NRSPINORS*WDES1%NGVECTOR,N)=WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
         ENDDO

         ! compute S|phi(0)> and store it in CH
         IF (INFO%LREAL) THEN
            CWORK1=0; LDO=.FALSE.; LDO(1:NSTRIP_ACT)=.TRUE.; DUMMY=0
            CALL RACCMU(NONLR_S,WDES1,WTMP,LMDIM,CQIJ,CQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSTRIP,LDO)
            DO N=1,NSTRIP_ACT
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,N),CH(1+ISPINOR*WDES1%NGVECTOR,N),GRID,.TRUE.)
               ENDDO
            ENDDO
         ELSE
            DO N=1,NSTRIP_ACT
               CALL VNLACC_ADD(NONL_S,WTMP(N),CQIJ,CQIJ,1,0._q,CH(:,N))
            ENDDO
         ENDIF

         ! subtract S|phi(0)_j><phi(0)_j|phi(1)_i> from |phi(1)_i>
         ! the latter are distributed over nodes
         DO NB=1,WXI%WDES%NBANDS
            DO N=1,NSTRIP_ACT
               IF (W0%FERTOT(NB_START+N-1,NK,ISP)<1E-4_q) CYCLE
               COVL=zero
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,WDES1%NGVECTOR
                     MM=M+WDES1%NGVECTOR*ISPINOR
                     COVL=COVL+GCONJG(WTMP(N)%CW(MM))*WXI%CW(MM,NB,NK,ISP)
                  ENDDO
               ENDDO
! test
!              WRITE(70+W0%WDES%COMM%NODE_ME,'(2I4,2F14.7)') N+NB_START-1,W0%WDES%NB_LOW+(NB-1)*W0%WDES%NB_PAR,COVL
! test
               WXI%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)=WXI%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)- &
              &   COVL*CH(1:WDES1%NRSPINORS*WDES1%NGVECTOR,N)
            ENDDO
         ENDDO
      ENDDO
!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      DO N=1,NSTRIP
         CALL DELWAV(WTMP(N),INFO%LREAL)
      ENDDO
      DEALLOCATE(WTMP,DUMMY,LDO,CH)
      IF (INFO%LREAL) DEALLOCATE(CWORK1)

      RETURN
      END SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE

!> *************** SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS *****************
!>
!> ****************************************************************************
      SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS(&
     & W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)
      USE prec
      USE base
      USE dfast
      USE lattice
      USE mgrid
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (wavespin) WXI
      TYPE (grid_3d) GRID
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      INTEGER LMDIM
      OVERLAP CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      ! local variables
      TYPE (wavedes1) WDES1
      TYPE (wavefun1), ALLOCATABLE :: WTMP(:)
      ! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW0_RED(:,:),CW1_RED(:,:),CH_RED(:,:)
      INTEGER NCPU,NRPLWV_RED,NPL,NPRO
      LOGICAL DO_REDIS

      INTEGER ISP,NK,ISPINOR
      INTEGER NB_TOT,NBANDS,NSTRIP,NSTRIP_ACT,NB_START,N
      INTEGER NPOSB1,NPOSB2,NSTRPB,NPOSPL,NSTRPL

      REAL(q), ALLOCATABLE :: DUMMY(:)
      LOGICAL, ALLOCATABLE :: LDO(:)

      COMPLEX(q), ALLOCATABLE, TARGET :: CH(:,:)
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)

      GDEF,ALLOCATABLE,TARGET :: COVL(:,:)

#ifdef MPI
      ! number of procs involved in band distribution
      NCPU=W0%WDES%COMM_INTER%NCPU
#else
      NCPU=1
#endif
      IF (NCPU/=1) THEN
         DO_REDIS=.TRUE.
         NRPLWV_RED=W0%WDES%NRPLWV/NCPU
      ELSE
         DO_REDIS=.FALSE.
         NRPLWV_RED=W0%WDES%NRPLWV
      ENDIF

      NB_TOT=W0%WDES%NB_TOT
      NBANDS=W0%WDES%NBANDS

!     ! set NSTRIP between [1 and 32]
!     NSTRIP=NSTRIP_STANDARD_GLOBAL
      ! the above was much too large
      NSTRIP=MIN(4,NBANDS)

      CALL SETWDES(W0%WDES,WDES1,0)
      ALLOCATE(WTMP(NSTRIP))
      DO N=1,NSTRIP
         CALL NEWWAV(WTMP(N),WDES1,INFO%LREAL)
      ENDDO

      ALLOCATE(DUMMY(NSTRIP),LDO(NSTRIP))

      IF (INFO%LREAL) THEN
         ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,NSTRIP))
      ENDIF

      ALLOCATE(CH(W0%WDES%NRPLWV,NBANDS))
      ! get pointers for redistribution to "over-plane-wave coefficients"
      IF (NCPU/=1) THEN
         CALL SET_WPOINTER( CH_RED, NRPLWV_RED, NB_TOT, CH(1,1))
      ELSE
         CH_RED => CH(:,:)
      ENDIF

      ALLOCATE(COVL(NB_TOT,NB_TOT))
!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,NKPTS_ORIG
#ifdef MPI
      IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
      CALL SETWDES(W0%WDES,WDES1,NK+NKPTS_ORIG)

      IF (INFO%LREAL) THEN
         CALL PHASER(GRID,LATT_CUR,NONLR_S,NK+NKPTS_ORIG,W0%WDES)
      ELSE
         NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK+NKPTS_ORIG)
      ENDIF

      ! compute S|phi(0)> and store it in CH
      CH=(0._q,0._q)
      DO NB_START=1,NBANDS,NSTRIP
         NSTRIP_ACT=MIN(NSTRIP,NBANDS-NB_START+1)

         DO N=1,NSTRIP_ACT
            CALL W1_COPY(ELEMENT(W0,WDES1,NB_START+N-1,ISP),WTMP(N)); IF (INFO%LREAL) CALL FFTWAV_W1(WTMP(N))
            CH(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB_START+N-1)=WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
         ENDDO

         IF (INFO%LREAL) THEN
            CWORK1=0; LDO=.FALSE.; LDO(1:NSTRIP_ACT)=.TRUE.; DUMMY=0
            CALL RACCMU(NONLR_S,WDES1,WTMP,LMDIM,CQIJ,CQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSTRIP,LDO)
            DO N=1,NSTRIP_ACT
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,N),CH(1+ISPINOR*WDES1%NGVECTOR,NB_START+N-1),GRID,.TRUE.)
               ENDDO
            ENDDO
         ELSE
            DO N=1,NSTRIP_ACT
               CALL VNLACC_ADD(NONL_S,WTMP(N),CQIJ,CQIJ,1,0._q,CH(:,NB_START+N-1))
            ENDDO
         ENDIF
      ENDDO

      ! get pointers for redistribution to "over-plane-wave coefficients"
      IF (NCPU/=1) THEN
         CALL SET_WPOINTER(CW0_RED, NRPLWV_RED, NB_TOT,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL SET_WPOINTER(CW1_RED, NRPLWV_RED, NB_TOT, WXI%CW(1,1,NK,ISP))
      ELSE
         CW0_RED =>  W0%CW(:,:,NK+NKPTS_ORIG,ISP)
         CW1_RED => WXI%CW(:,:,NK,ISP)
      ENDIF

      ! set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO
      CALLMPI(SET_NPL_NPRO(WDES1, NPL, NPRO))

      ! redistribute to "over-plane-wave coefficients"
      IF (DO_REDIS) THEN
         CALL REDIS_PW(WDES1, NBANDS,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL REDIS_PW(WDES1, NBANDS, WXI%CW(1,1,NK,ISP))
         CALL REDIS_PW(WDES1, NBANDS,     CH(1,1)       )
      ENDIF

      ! calculate COVL(i,j)=<phi(0)_i|phi(1)_j>
      COVL=zero
      IF (NPL/=0) THEN
!     CALL GGEMM('C','N',NB_TOT,NB_TOT,m_ NPL,one,CW0_RED(1,1),m_ NRPLWV_RED,CW1_RED(1,1),m_ NRPLWV_RED,zero,COVL(1,1),NB_TOT)
!     NSTRPB=NSTRIP
      NSTRPB=NB_TOT
      DO NPOSB1=1,NB_TOT,NSTRPB; DO NPOSB2=1,NB_TOT,NSTRPB
      CALL GGEMM('C','N',MIN(NSTRPB,NB_TOT-NPOSB2+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),m_ NPL,one,CW0_RED(1,NPOSB2),m_ NRPLWV_RED,CW1_RED(1,NPOSB1),m_ NRPLWV_RED,zero,COVL(NPOSB2,NPOSB1),NB_TOT)
      ENDDO; ENDDO
      ENDIF

      CALLMPI( M_sum_g(W0%WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))

! test
!     IF (W0%WDES%COMM%IONODE==W0%WDES%COMM%NODE_ME) THEN
!        DO N=1,W0%WDES%NB_TOT
!           IF (W0%FERTOT(N,NK,ISP)<1E-4_q) COVL(N,:)=zero
!           WRITE(80,'(32(2F14.7,2X))') COVL(N,:)
!        ENDDO
!     ENDIF
! test
      ! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
      DO N=1,W0%WDES%NB_TOT
         IF (W0%FERTOT(N,NK,ISP)<1E-4_q) COVL(N,:)=zero
      ENDDO

      ! subtract S|phi(0)_j><phi(0)_j|phi(1)_i> from |phi(1)_i>
!     CALL GGEMM('N','N',m_ NPL,NB_TOT,NB_TOT,-one,CH_RED(1,1),m_ NRPLWV_RED,COVL(1,1),NB_TOT,one,CW1_RED(1,1),m_ NRPLWV_RED)
!     NSTRPL=m_ NPL
!     DO NPOSPL=1,m_ NPL,NSTRPL
!     CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),NB_TOT,NB_TOT,-one,CH_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,1),NB_TOT,one,CW1_RED(NPOSPL,1),m_ NRPLWV_RED)
!     ENDDO
      IF (NPL/=0) THEN
      NSTRPB=NB_TOT; NSTRPL=m_ NPL
      DO NPOSB1=1,NB_TOT,NSTRPB; DO NPOSPL=1,m_ NPL,NSTRPL
      CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),NB_TOT,-one,CH_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,NPOSB1),NB_TOT,one,CW1_RED(NPOSPL,NPOSB1),m_ NRPLWV_RED)
      ENDDO; ENDDO
      ENDIF

      ! redistribute back to original data distribution
      IF (DO_REDIS) THEN
         CALL REDIS_PW(WDES1, NBANDS,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL REDIS_PW(WDES1, NBANDS, WXI%CW(1,1,NK,ISP))
      ENDIF

!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      ! deallocation
      IF (NCPU/=1) NULLIFY(CW0_RED,CW1_RED,CH_RED)

      DO N=1,NSTRIP
         CALL DELWAV(WTMP(N),INFO%LREAL)
      ENDDO

      DEALLOCATE(WTMP,DUMMY,LDO,CH,COVL)

      IF (INFO%LREAL) DEALLOCATE(CWORK1)

      RETURN
      END SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS

!> *************** SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS_ ****************
!>
!> ****************************************************************************
      SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS_(&
     & W0,GRID,NONL_S,NONLR_S,LATT_CUR,LMDIM,CQIJ,INFO,WXI)
      USE prec
      USE base
      USE dfast
      USE lattice
      USE mgrid
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (wavespin) WXI
      TYPE (grid_3d) GRID
      TYPE (nonl_struct) NONL_S
      TYPE (nonlr_struct) NONLR_S
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      INTEGER LMDIM
      OVERLAP CQIJ(LMDIM,LMDIM,W0%WDES%NIONS,W0%WDES%NCDIJ)
      ! local variables
      TYPE (wavedes1) WDES1
      TYPE (wavefun1), ALLOCATABLE :: WTMP(:)
      ! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW0_RED(:,:),CW1_RED(:,:),CH_RED(:,:)
      INTEGER NCPU,NRPLWV_RED,NPL,NPRO
      LOGICAL DO_REDIS

      INTEGER ISP,NK,ISPINOR
      INTEGER NB_TOT,NBANDS,NSTRIP,NSTRIP_ACT,NB_START,N
      INTEGER NPOSB1,NPOSB2,NSTRPB,NPOSPL,NSTRPL

      REAL(q), ALLOCATABLE :: DUMMY(:)
      LOGICAL, ALLOCATABLE :: LDO(:)

      COMPLEX(q), ALLOCATABLE, TARGET :: CH(:,:)
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)

      GDEF,ALLOCATABLE,TARGET :: COVL(:,:)

#ifdef MPI
      ! number of procs involved in band distribution
      NCPU=W0%WDES%COMM_INTER%NCPU
#else
      NCPU=1
#endif
      IF (NCPU/=1) THEN
         DO_REDIS=.TRUE.
         NRPLWV_RED=W0%WDES%NRPLWV/NCPU
      ELSE
         DO_REDIS=.FALSE.
         NRPLWV_RED=W0%WDES%NRPLWV
      ENDIF

      NB_TOT=W0%WDES%NB_TOT
      NBANDS=W0%WDES%NBANDS

!     ! set NSTRIP between [1 and 32]
!     NSTRIP=NSTRIP_STANDARD_GLOBAL
      ! the above was much too large
      NSTRIP=MIN(4,NBANDS)

      CALL SETWDES(W0%WDES,WDES1,0)
      ALLOCATE(WTMP(NSTRIP))
      DO N=1,NSTRIP
         CALL NEWWAV(WTMP(N),WDES1,INFO%LREAL)
      ENDDO

      ALLOCATE(DUMMY(NSTRIP),LDO(NSTRIP))

      IF (INFO%LREAL) THEN
         ALLOCATE(CWORK1(GRID%MPLWV*W0%WDES%NRSPINORS,NSTRIP))
      ENDIF

      ALLOCATE(CH(W0%WDES%NRPLWV,NBANDS))
      ! get pointers for redistribution to "over-plane-wave coefficients"
      IF (NCPU/=1) THEN
         CALL SET_WPOINTER( CH_RED, NRPLWV_RED, NB_TOT, CH(1,1))
      ELSE
         CH_RED => CH(:,:)
      ENDIF

      ALLOCATE(COVL(NB_TOT,NB_TOT))

!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,NKPTS_ORIG
#ifdef MPI
      IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
      CALL SETWDES(W0%WDES,WDES1,NK+NKPTS_ORIG)

      IF (INFO%LREAL) THEN
         CALL PHASER(GRID,LATT_CUR,NONLR_S,NK+NKPTS_ORIG,W0%WDES)
      ELSE
         NONL_S%NK=-1; CALL PHASE(W0%WDES,NONL_S,NK+NKPTS_ORIG)
      ENDIF

      ! compute S|phi(0)> and store it in CH
      CH=(0._q,0._q)
      DO NB_START=1,NBANDS,NSTRIP
         NSTRIP_ACT=MIN(NSTRIP,NBANDS-NB_START+1)

         DO N=1,NSTRIP_ACT
            CALL W1_COPY(ELEMENT(W0,WDES1,NB_START+N-1,ISP),WTMP(N)); IF (INFO%LREAL) CALL FFTWAV_W1(WTMP(N))
            CH(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB_START+N-1)=WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
         ENDDO

         IF (INFO%LREAL) THEN
            CWORK1=0; LDO=.FALSE.; LDO(1:NSTRIP_ACT)=.TRUE.; DUMMY=0
            CALL RACCMU(NONLR_S,WDES1,WTMP,LMDIM,CQIJ,CQIJ,DUMMY,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSTRIP,LDO)
            DO N=1,NSTRIP_ACT
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,N),CH(1+ISPINOR*WDES1%NGVECTOR,NB_START+N-1),GRID,.TRUE.)
               ENDDO
            ENDDO
         ELSE
            DO N=1,NSTRIP_ACT
               CALL VNLACC_ADD(NONL_S,WTMP(N),CQIJ,CQIJ,1,0._q,CH(:,NB_START+N-1))
            ENDDO
         ENDIF
      ENDDO

      ! get pointers for redistribution to "over-plane-wave coefficients"
      IF (NCPU/=1) THEN
         CALL SET_WPOINTER(CW0_RED, NRPLWV_RED, NB_TOT,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL SET_WPOINTER(CW1_RED, NRPLWV_RED, NB_TOT, WXI%CW(1,1,NK,ISP))
      ELSE
         CW0_RED =>  W0%CW(:,:,NK+NKPTS_ORIG,ISP)
         CW1_RED => WXI%CW(:,:,NK,ISP)
      ENDIF

      ! set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO
      CALLMPI(SET_NPL_NPRO(WDES1, NPL, NPRO))

      ! redistribute to "over-plane-wave coefficients"
      IF (DO_REDIS) THEN
         CALL REDIS_PW(WDES1, NBANDS,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL REDIS_PW(WDES1, NBANDS, WXI%CW(1,1,NK,ISP))
         CALL REDIS_PW(WDES1, NBANDS,     CH(1,1)       )
      ENDIF

      ! calculate COVL(i,j)=<phi(0)_i|S|phi(1)_j>
      COVL=zero
!     CALL GGEMM('C','N',NB_TOT,NB_TOT,m_ NPL,one,CW0_RED(1,1),m_ NRPLWV_RED,CW1_RED(1,1),m_ NRPLWV_RED,zero,COVL(1,1),NB_TOT)
!     NSTRPB=NSTRIP
      NSTRPB=NB_TOT
      IF (NPL/=0) THEN
      DO NPOSB1=1,NB_TOT,NSTRPB; DO NPOSB2=1,NB_TOT,NSTRPB
      CALL GGEMM('C','N',MIN(NSTRPB,NB_TOT-NPOSB2+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),m_ NPL,one,CH_RED(1,NPOSB2),m_ NRPLWV_RED,CW1_RED(1,NPOSB1),m_ NRPLWV_RED,zero,COVL(NPOSB2,NPOSB1),NB_TOT)
      ENDDO; ENDDO
      ENDIF

      CALLMPI( M_sum_g(W0%WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))

! test
!     IF (W0%WDES%COMM%IONODE==W0%WDES%COMM%NODE_ME) THEN
!        DO N=1,W0%WDES%NB_TOT
!           IF (W0%FERTOT(N,NK,ISP)<1E-4_q) COVL(N,:)=zero
!           WRITE(80,'(32(2F14.7,2X))') COVL(N,:)
!        ENDDO
!     ENDIF
! test
      ! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
      DO N=1,W0%WDES%NB_TOT
         IF (W0%FERTOT(N,NK,ISP)<1E-4_q) COVL(N,:)=zero
      ENDDO

      ! subtract |phi(0)_j><phi(0)_j|S|phi(1)_i> from |phi(1)_i>
!     CALL GGEMM('N','N',m_ NPL,NB_TOT,NB_TOT,-one,CH_RED(1,1),m_ NRPLWV_RED,COVL(1,1),NB_TOT,one,CW1_RED(1,1),m_ NRPLWV_RED)
!     NSTRPL=m_ NPL
!     DO NPOSPL=1,m_ NPL,NSTRPL
!     CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),NB_TOT,NB_TOT,-one,CH_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,1),NB_TOT,one,CW1_RED(NPOSPL,1),m_ NRPLWV_RED)
!     ENDDO
      NSTRPB=NB_TOT; NSTRPL=m_ NPL
      IF (NPL/=0) THEN
      DO NPOSB1=1,NB_TOT,NSTRPB; DO NPOSPL=1,m_ NPL,NSTRPL
      CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),NB_TOT,-one,CW0_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,NPOSB1),NB_TOT,one,CW1_RED(NPOSPL,NPOSB1),m_ NRPLWV_RED)
      ENDDO; ENDDO
      ENDIF

      ! redistribute back to original data distribution
      IF (DO_REDIS) THEN
         CALL REDIS_PW(WDES1, NBANDS,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL REDIS_PW(WDES1, NBANDS, WXI%CW(1,1,NK,ISP))
      ENDIF

!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      ! deallocation
      IF (NCPU/=1) NULLIFY(CW0_RED,CW1_RED,CH_RED)

      DO N=1,NSTRIP
         CALL DELWAV(WTMP(N),INFO%LREAL)
      ENDDO

      DEALLOCATE(WTMP,DUMMY,LDO,CH,COVL)

      IF (INFO%LREAL) DEALLOCATE(CWORK1)

      RETURN
      END SUBROUTINE MLR_PROJ_IN_EMPTY_SUBSPACE_BLAS_

!> ******************* SUBROUTINE MLR_PROJ_IN_OCC_SUBSPACE ********************
!>
!> ****************************************************************************
      SUBROUTINE MLR_PROJ_IN_OCC_SUBSPACE(W0,WXI,W1)
      USE prec
      USE dfast
      USE wave_high
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (wavespin) W1
      TYPE (wavespin) WXI
      ! local variables
      TYPE (wavedes1) WDES1
      TYPE (wavefun1), ALLOCATABLE :: WTMP(:)

      INTEGER ISP,NK,NB,N,ISPINOR,M,MM
      INTEGER NB_TOT,NB_START,NB_STOP,NSTRIP,NSTRIP_ACT

      GDEF COVL

      NB_TOT=W0%WDES%NB_TOT

      ! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

      CALL SETWDES(W0%WDES,WDES1,0)
      ALLOCATE(WTMP(NSTRIP))
      DO N=1,NSTRIP
         CALL NEWWAV(WTMP(N),WDES1,.FALSE.)
      ENDDO

!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,NKPTS_ORIG
#ifdef MPI
      IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
      CALL SETWDES(W0%WDES,WDES1,NK+NKPTS_ORIG)

!> ----------------------------------------------------------------------------
!> calculate
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!>  |u(1)_k+q,n>
!>   - \sum_j | u(0)_j,k+q >< u(0)_j,k+q | (1/i) [r_idir,S_k+q,k] |u(0)_nk>
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> where j runs over occupied orbitals
!> ----------------------------------------------------------------------------

      DO NB_START=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NSTRIP,NB_TOT-NB_START+1)
         NB_STOP=MIN(NB_START+NSTRIP-1,NB_TOT)
         ! collect the wave functions from NB_START to NB_STOP (global band index) into WTMP
         CALL W1_GATHER_GLB_NOCR(W0,NB_START,NB_STOP,ISP,WTMP)

         DO NB=1,W1%WDES%NBANDS
            DO N=1,NSTRIP_ACT
               IF (W0%FERTOT(NB_START+N-1,NK,ISP)<1E-4_q) CYCLE
               COVL=zero
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,WDES1%NGVECTOR
                     MM=M+WDES1%NGVECTOR*ISPINOR
                     COVL=COVL+GCONJG(WTMP(N)%CW(MM))*WXI%CW(MM,NB,NK,ISP)
                  ENDDO
               ENDDO
! test
!              W1%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)=W1%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)- &
!             &   COVL*WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
               W1%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)=W1%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)+ &
              &   COVL*WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
! test
            ENDDO
         ENDDO
      ENDDO
!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      DO N=1,NSTRIP
         CALL DELWAV(WTMP(N),.FALSE.)
      ENDDO
      DEALLOCATE(WTMP)

      RETURN
      END SUBROUTINE MLR_PROJ_IN_OCC_SUBSPACE

!> ***************** SUBROUTINE MLR_PROJ_IN_OCC_SUBSPACE_BLAS *****************
!>
!> ****************************************************************************
      SUBROUTINE MLR_PROJ_IN_OCC_SUBSPACE_BLAS(&
     &   W0,WXI,W1)
      USE prec
      USE base
      USE dfast
      USE lattice
      USE mgrid
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (wavespin) WXI
      TYPE (wavespin) W1
      ! local variables
      TYPE (wavedes1) WDES1
      ! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW0_RED(:,:),CW1_RED(:,:),CWX_RED(:,:)
      INTEGER NCPU,NRPLWV_RED,NPL,NPRO
      LOGICAL DO_REDIS

      INTEGER ISP,NK
      INTEGER NB_TOT,NBANDS,N
      INTEGER NPOSB1,NPOSB2,NSTRPB,NPOSPL,NSTRPL

      GDEF,ALLOCATABLE,TARGET :: COVL(:,:)

#ifdef MPI
      ! number of procs involved in band distribution
      NCPU=W0%WDES%COMM_INTER%NCPU
#else
      NCPU=1
#endif
      IF (NCPU/=1) THEN
         DO_REDIS=.TRUE.
         NRPLWV_RED=W0%WDES%NRPLWV/NCPU
      ELSE
         DO_REDIS=.FALSE.
         NRPLWV_RED=W0%WDES%NRPLWV
      ENDIF

      NB_TOT=W0%WDES%NB_TOT
      NBANDS=W0%WDES%NBANDS

      ALLOCATE(COVL(NB_TOT,NB_TOT))
!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,NKPTS_ORIG
#ifdef MPI
      IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
      CALL SETWDES(W0%WDES,WDES1,NK+NKPTS_ORIG)

      ! get pointers for redistribution to "over-plane-wave coefficients"
      IF (NCPU/=1) THEN
         CALL SET_WPOINTER(CW0_RED, NRPLWV_RED, NB_TOT,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL SET_WPOINTER(CW1_RED, NRPLWV_RED, NB_TOT,  W1%CW(1,1,NK,ISP))
         CALL SET_WPOINTER(CWX_RED, NRPLWV_RED, NB_TOT, WXI%CW(1,1,NK,ISP))
      ELSE
         CW0_RED =>  W0%CW(:,:,NK+NKPTS_ORIG,ISP)
         CW1_RED =>  W1%CW(:,:,NK,ISP)
         CWX_RED => WXI%CW(:,:,NK,ISP)
      ENDIF

      ! set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO
      CALLMPI(SET_NPL_NPRO(WDES1, NPL, NPRO))

      ! redistribute to "over-plane-wave coefficients"
      IF (DO_REDIS) THEN
         CALL REDIS_PW(WDES1, NBANDS,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL REDIS_PW(WDES1, NBANDS,  W1%CW(1,1,NK,ISP))
         CALL REDIS_PW(WDES1, NBANDS, WXI%CW(1,1,NK,ISP))
      ENDIF

      ! calculate COVL(i,j)=<phi(0)_i|WXI_j>
      COVL=zero
      IF (NPL/=0) THEN
      NSTRPB=NB_TOT
      DO NPOSB1=1,NB_TOT,NSTRPB; DO NPOSB2=1,NB_TOT,NSTRPB
      CALL GGEMM('C','N',MIN(NSTRPB,NB_TOT-NPOSB2+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),m_ NPL,one,CW0_RED(1,NPOSB2),m_ NRPLWV_RED,CWX_RED(1,NPOSB1),m_ NRPLWV_RED,zero,COVL(NPOSB2,NPOSB1),NB_TOT)
      ENDDO; ENDDO
      ENDIF

      CALLMPI( M_sum_g(W0%WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))

      ! set COVL(i,j)=0 if |phi(0)_i> is an empty state.
      DO N=1,W0%WDES%NB_TOT
         IF (W0%FERTOT(N,NK,ISP)<1E-4_q) COVL(N,:)=zero
      ENDDO

      ! add |phi(0)_j><phi(0)_j|WXI_i> to |phi(1)_i>
      IF (NPL/=0) THEN
      NSTRPB=NB_TOT; NSTRPL=m_ NPL
      DO NPOSB1=1,NB_TOT,NSTRPB; DO NPOSPL=1,m_ NPL,NSTRPL
      CALL GGEMM('N','N',MIN(NSTRPL,m_ NPL-NPOSPL+1),MIN(NSTRPB,NB_TOT-NPOSB1+1),NB_TOT,one,CW0_RED(NPOSPL,1),m_ NRPLWV_RED,COVL(1,NPOSB1),NB_TOT,one,CW1_RED(NPOSPL,NPOSB1),m_ NRPLWV_RED)
      ENDDO; ENDDO
      ENDIF

      ! redistribute back to original data distribution
      IF (DO_REDIS) THEN
         CALL REDIS_PW(WDES1, NBANDS,  W0%CW(1,1,NK+NKPTS_ORIG,ISP))
         CALL REDIS_PW(WDES1, NBANDS,  W1%CW(1,1,NK,ISP))
         CALL REDIS_PW(WDES1, NBANDS, WXI%CW(1,1,NK,ISP))
      ENDIF

!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      ! deallocation
      IF (NCPU/=1) NULLIFY(CW0_RED,CW1_RED,CWX_RED)
      DEALLOCATE(COVL)

      RETURN
      END SUBROUTINE MLR_PROJ_IN_OCC_SUBSPACE_BLAS

!> ********************* SUBROUTINE MLR_APPLY_GREENS_FUNC *********************
!>
!> ****************************************************************************
      SUBROUTINE MLR_APPLY_GREENS_FUNC(W0,WXI,W1)
      USE prec
      USE dfast
      USE wave_high
      IMPLICIT NONE
      TYPE (wavespin) W0
      TYPE (wavespin) W1
      TYPE (wavespin) WXI
      ! local variables
      TYPE (wavedes1) WDES1
      TYPE (wavefun1), ALLOCATABLE :: WTMP(:)

      INTEGER ISP,NK,NB,N,ISPINOR,M,MM
      INTEGER NB_TOT,NB_START,NB_STOP,NSTRIP,NSTRIP_ACT

      REAL(q), PARAMETER :: TINY=1.E-6_q

      GDEF COVL

      NB_TOT=W0%WDES%NB_TOT

      ! set NSTRIP between [1 and 32]
      NSTRIP=NSTRIP_STANDARD_GLOBAL

      CALL SETWDES(W0%WDES,WDES1,0)
      ALLOCATE(WTMP(NSTRIP))
      DO N=1,NSTRIP
         CALL NEWWAV(WTMP(N),WDES1,.FALSE.)
      ENDDO

!> ============================================================================
      spin:  DO ISP=1,W0%WDES%ISPIN
      kpoints: DO NK=1,NKPTS_ORIG
#ifdef MPI
      IF (MOD(NK-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!> ============================================================================
      CALL SETWDES(W0%WDES,WDES1,NK+NKPTS_ORIG)

!> ----------------------------------------------------------------------------
!> calculate
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!>  |u(1)_k+q,n>
!>            | u(0)_j,k+q >< u(0)_j,k+q |
!>   + \sum_j ---------------------------- |Xi(0)_nk+q>
!>                e(0)_n,k - e(0)_j,k+q
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> where j runs over occupied orbitals
!> ----------------------------------------------------------------------------

      DO NB_START=1,NB_TOT,NSTRIP
         NSTRIP_ACT=MIN(NSTRIP,NB_TOT-NB_START+1)
         NB_STOP=MIN(NB_START+NSTRIP-1,NB_TOT)
         ! collect the wave functions from NB_START to NB_STOP (global band index) into WTMP
         CALL W1_GATHER_GLB_NOCR(W0,NB_START,NB_STOP,ISP,WTMP)

         DO NB=1,W1%WDES%NBANDS
            DO N=1,NSTRIP_ACT
               IF (W0%FERTOT(NB_START+N-1,NK,ISP)>1E-4_q) CYCLE
               COVL=zero
               DO ISPINOR=0,WDES1%NRSPINORS-1
                  DO M=1,WDES1%NGVECTOR
                     MM=M+WDES1%NGVECTOR*ISPINOR
                     COVL=COVL+GCONJG(WTMP(N)%CW(MM))*WXI%CW(MM,NB,NK,ISP)
                  ENDDO
               ENDDO
               IF (ABS(REAL(W0%CELEN(NB,NK,ISP)-W0%CELTOT(NB_START+N-1,NK+NKPTS_ORIG,ISP),q))>TINY) THEN
                  W1%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)=W1%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR,NB,NK,ISP)- &
                 &   COVL/REAL(W0%CELEN(NB,NK,ISP)-W0%CELTOT(NB_START+N-1,NK+NKPTS_ORIG,ISP),q) &
                 &   *WTMP(N)%CW(1:WDES1%NRSPINORS*WDES1%NGVECTOR)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!> ============================================================================
      ENDDO kpoints
      ENDDO spin
!> ============================================================================

      DO N=1,NSTRIP
         CALL DELWAV(WTMP(N),.FALSE.)
      ENDDO
      DEALLOCATE(WTMP)

      RETURN
      END SUBROUTINE MLR_APPLY_GREENS_FUNC

!> *************************** SUBROUTINE MLR_SBARE ***************************
!>
!> Calculates S_{bare} as per the theory
!> Orbital contribution contained in variable S_BARE
!> Spin contribution, denoted M_{bare} in the theory, contained in S_DENS
!>
!> first term RHS eq. 67 arXiv:2504.17380 [cond-mat.mtrl-sci]
!> (if LSOSHIFT, also first term RHS eq. 75)
!>
!> ****************************************************************************
      SUBROUTINE MLR_SBARE(W,WP,S_BARE,S_DENS,LSOSHIFT)
      USE prec
      USE wave_high
      IMPLICIT NONE
      TYPE (wavespin) W,WP
      RGRID :: S_BARE(W%WDES%GRID%RL%NP,3)
      RGRID :: S_DENS(W%WDES%GRID%RL%NP,3)
      LOGICAL LSOSHIFT
      ! local variables
      TYPE (wavedes1) WDES1,WDESP1
      TYPE (wavefun1) W1,WP1
      COMPLEX(q), ALLOCATABLE :: CWORK1(:)
      COMPLEX(q), ALLOCATABLE :: CWORK2(:)
      INTEGER ISP,NK,N
      REAL(q) WEIGHT

      ALLOCATE(CWORK1(W%WDES%NRPLWV))
      ALLOCATE(CWORK2(W%WDES%GRID%MPLWV))

      CALL SETWDES(W%WDES ,WDES1 ,0)
      CALL SETWDES(WP%WDES,WDESP1,0)
      CALL NEWWAV(W1 ,WDES1, .TRUE.)
      CALL NEWWAV(WP1,WDESP1,.TRUE.)

      spin:  DO ISP=1,W%WDES%ISPIN
      kpoint: DO NK=1,WP%WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         CALL SETWDES(W%WDES ,WDES1 ,NK)
         CALL SETWDES(WP%WDES,WDESP1,NK)  ! but this one is shifted with q

         DO N=1,W%WDES%NBANDS
            CALL W1_COPY( ELEMENT( W , WDES1 , N, ISP), W1 )
            CALL W1_COPY( ELEMENT( WP, WDESP1, N, ISP), WP1)
            CALL FFTWAV_W1(W1 )
            CALL FFTWAV_W1(WP1)
            WEIGHT=WDES1%RSPIN*W%FERWE(N,NK,ISP)*W%WDES%WTKPT(NK)     ! occupancies assumed identical

            CALL PSCURRENT_PARA_TWO_KPT(WDES1, WDES1%GRID, W1%CR(1), WP1%CR(1), WDES1%IGX(1), &
           &     WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT, WDESP1%VKPT, WEIGHT, CWORK1(1), CWORK2(1), &
           &     W1%CW(1), WP1%CW(1), S_BARE(1,1), S_DENS(1,1), LSOSHIFT)

         ENDDO
      ENDDO kpoint
      ENDDO spin

      CALL DELWAV(W1  ,.TRUE.)
      CALL DELWAV(WP1,.TRUE.)

      DEALLOCATE(CWORK1)
      DEALLOCATE(CWORK2)

      RETURN
      END SUBROUTINE MLR_SBARE

!> **************************** CALC_B_MAGATOM_LR *****************************
!>
!> Taken from CALC_B_MAGATOM, augmentation removed
!>
!> To calculate one-centre contribution to magnetic field at atomic nuclei.
!> First paramagnetic and diamagnetic "currents" are calculated:
!>
!> COMMENT: change from p -> -i nabla
!> COMMENT: and  A to M x r / r^3
!>
!> para:
!> ~~~~~~~~~~~~~~~~~~~~~
!> J_LM^alpha(r) =   Sum_ij \int \phi_i ((-i nabla)^* + -i nabla )/2 Y_LM \phi_j dOmega \rho_ji
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> dia:
!> ~~~~~~~~~~~~~~~~~~~~~
!> J_LM^alpha(r) =   - Sum_ij \int \phi_i  M x r / r^3   Y_LM \phi_j dOmega \rho_ji
!>                   * MOMTOMOM
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> These are collected into one array
!>
!> Next the augmentation currents are NOT subracted
!>
!> Finally Biot-Savart is applied to obtain the field:
!> ~~~~~~~~~~~~~~~~~~~~~
!>   B = - r x J /|r^3|
!> ~~~~~~~~~~~~~~~~~~~~~
!> ~~~~~~~~~~~~~~~~~~~~~
!> A=(Bxr)/2*MAGMOMTOENERGY
!> ~~~~~~~~~~~~~~~~~~~~~
!> This routine is similar to and derived from CALC_MU_MAGATOM
!> It is different though, as the augmentation current has to be subtracted in
!> a "different" place.
!>
!> Note for SOC/ZORA:
!> LBONE in this subroutine reconstructs the small-component in a ZORA-like
!> fashion (see J.Chem.Phys. 146, 054115 (2017))
!> this is not equal to what we intend for LZORA, hence both should be implemented
!> LBONE and LZORA should not be used simultaneously
!>
!> ****************************************************************************

      SUBROUTINE CALC_B_MAGATOM_LR(PP,CRHODE,BDIR,IDP,NIP,NI,B_OUT,LLRAUG)
      USE pseudo
      USE constant
      USE asa
      USE radial
      USE wave
      USE paw
      USE poscar
      USE prec
      USE morbitalmag, ONLY : LBONE, LZORA

      IMPLICIT NONE

      TYPE (potcar), POINTER :: PP !< pseudopotential descriptor
      COMPLEX(q) :: CRHODE(:,:) !< one center occupancy matrix
      REAL(q) :: B_OUT(3)
      INTEGER :: NIP,NI,BDIR
      LOGICAL LLRAUG

      INTEGER, EXTERNAL :: MAXL1

      !local variables
      INTEGER CHANNELS          !< number of channels
      INTEGER NMAX,I,NT
      INTEGER I0,I1
      INTEGER L0,L1,L2,L3
      INTEGER M0,M1,M3
      INTEGER LM0,LM1,LM2,LM3,LM2MAX
      INTEGER ISTART1,IEND1,LMIND1,IC1
      INTEGER LMAXNABLA,LMAX,K,IDP
      INTEGER II1,II2,II3,LL,M
      INTEGER LRHO,LMRHO,MRHO,L2MAX
      INTEGER IRMAX
      REAL(q) :: CGOR
      REAL(q) :: SUM
      REAL(q), ALLOCATABLE :: WTMP(:),DWAE(:),DWPS(:),TMP(:)
      REAL(q), ALLOCATABLE :: WTMP2(:),WTMP3(:),WBAE(:,:),DWBAE(:)
      REAL(q), ALLOCATABLE :: YLM_NABLA_YLM(:,:,:),YLM_X_YLM(:,:,:),YLM_BxR_YLM(:,:,:)
      REAL(q), ALLOCATABLE :: JPAWRAD(:,:,:)
      REAL(q) :: RHOLM(256)
      COMPLEX(q) :: JLM(256,3)
      REAL(q) :: OSUM
      REAL(q) :: OSUM1,OSUM2
      REAL(q), ALLOCATABLE :: OTMP(:)
      COMPLEX(q) :: BF(3)
      REAL(q) :: BCONST_LR(3)

      REAL(q), PARAMETER :: C = 137.037  !< speed of light in a.u.
      REAL(q), PARAMETER :: INVMC2=7.45596E-6 !< hbar^2/2(m_e c)^2 in A^2

      REAL(q), ALLOCATABLE :: POT(:),KQ(:)
      REAL(q), ALLOCATABLE :: WEIGHT_WAAE(:,:),WEIGHT_WBAE(:,:)
      REAL(q) RDUM
      REAL(q) AAE_CHAR,BAE_CHAR,WEIGHT_TMP

      IF (.NOT. ASSOCIATED(PP%JPAW)) RETURN
      BCONST_LR=0._q
      IF (IDP.EQ.2) THEN  ! diamagnetic
         BCONST_LR(BDIR)=1._q
      ENDIF

      ! allocate
      NMAX=PP%R%NMAX
      ALLOCATE(WTMP(NMAX),DWAE(NMAX),DWPS(NMAX),TMP(NMAX),OTMP(NMAX))
      ALLOCATE(WTMP2(NMAX),WTMP3(NMAX),WBAE(NMAX,PP%LMAX),DWBAE(NMAX))
      ALLOCATE(POT(NMAX),KQ(NMAX),WEIGHT_WAAE(NMAX,PP%LMAX),WEIGHT_WBAE(NMAX,PP%LMAX))

      CHANNELS = PP%LMAX
      LMAX=MAXL1(PP)
      LMAXNABLA=MAXL1(PP)+1   !< nabla Y_lm has components up to LMAX+1
      ALLOCATE(JPAWRAD(PP%R%NMAX,(LMAXNABLA+LMAX+1)*(LMAXNABLA+LMAX+1),3)) ! beware

      ALLOCATE(YLM_NABLA_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAXNABLA+2)*(2*LMAXNABLA+2),0:3))
      ALLOCATE(YLM_X_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAX+2)*(2*LMAX+2),0:3))
      ALLOCATE(YLM_BxR_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAX+2)*(2*LMAX+2),0:3))

      YLM_NABLA_YLM=0._q
      YLM_X_YLM=0._q
      YLM_BxR_YLM=0._q

!     CALL SETYLM_NABLA_YLM(LMAXNABLA,YLM_NABLA_YLM,YLM_X_YLM)
! in principle we should go up to LMAX+LMAXNABLA+1 = 2 LMAX + 2
! see below, but that does not work for f electrons
      CALL SETYLM_NABLA_YLM(2*LMAX+1,YLM_NABLA_YLM,YLM_X_YLM)

      YLM_BxR_YLM(:,:,1)=BCONST_LR(2)*YLM_X_YLM(:,:,3)-BCONST_LR(3)*YLM_X_YLM(:,:,2)
      YLM_BxR_YLM(:,:,2)=BCONST_LR(3)*YLM_X_YLM(:,:,1)-BCONST_LR(1)*YLM_X_YLM(:,:,3)
      YLM_BxR_YLM(:,:,3)=BCONST_LR(1)*YLM_X_YLM(:,:,2)-BCONST_LR(2)*YLM_X_YLM(:,:,1)
      YLM_BxR_YLM(:,:,:)=YLM_BxR_YLM(:,:,:)/2. !*MAGMOMTOENERGY

      JPAWRAD=0._q
      B_OUT=0._q

      IF (LLRAUG) THEN
         DO I=1,3
            CALL CALC_RHOLM_JVEC_CMPLX(  CRHODE, JLM(:,I), PP, I)
         ENDDO
      ENDIF

      !> J.Chem.Phys. 146, 054115 (2017)
      !> LBONE reconstructs the small-component in a ``ZORA-like'' fashion
      !> (Served for comparison to Wien2K)
      !>
      !> LBONE NOT the same as ZORA K-factor (LZORA implemented separately)
      !> LBONE and LZORA should not be used simultaneously
      IF (LBONE) THEN
         POT=0._q
         ! compute V_H[n^a_c]
         CALL RAD_POT_HAR(0,PP%R,POT,PP%RHOAE,RDUM)
         ! add -Z/r
         POT=POT/(2._q*SQRT(PI))-(PP%ZCORE+PP%ZVALF_ORIG)*FELECT/PP%R%R
         ! add V_H[n^a_v] + V_xc[n^a_v+n^a_c]
         POT=POT-PP%POTAE
         ! and make sure this AE-"atomic reference"-potential lines up exactly with
         ! the potential that was constructed by the pseudopotential generation code.
         POT=POT+(PP%VPSRMAX-POT(NMAX))

         KQ(:)=1._q/(1._q-POT(:)/(2._q*RYTOEV)/CLIGHT/CLIGHT/2._q)
         KQ(:)=KQ(:)*(INVMC2**0.5_q)

! Alternative way to make LBONE using the NMR ZORA routines
! prefer to keep the original above approach to separate the two clearly
! additionally guarantees LBONE is not accidentally affected.
!         CALL NMR_ZORA_POT(PP, POT)
!
!!         KQ(:)=1._q/(1._q-POT(:)/(2._q*RYTOEV)/CLIGHT/CLIGHT/2._q)
!!         KQ(:)=KQ(:)*(INVMC2**0.5_q)
!         KQ(:)=NMR_ZORA_K(PP%R%NMAX, POT)*(INVMC2**0.5_q)

         DO I=1,NMAX-1
            IF ((PP%RDEP>0).AND.(PP%R%R(I)-PP%RDEP>-5E-3)) EXIT
         ENDDO
!        WRITE(*,*) 'outermost point',PP%R%R(I),PP%RDEP
         IRMAX=I
      ELSE
         KQ(:)=1._q
      ENDIF

!     OPEN (UNIT=99,FILE='ABTMPAE',STATUS="UNKNOWN")
      DO I1=1,PP%LMAX   !CHANNELS
         WTMP2(:)=PP%WAE(:,I1)/PP%R%R(:)
         CALL GRAD_(PP%R,WTMP2,WTMP3)
         WTMP3(:)=WTMP3(:)*KQ(:)/SQRT(2._q)
!        WTMP3(:)=WTMP3(:)/(2._q**0.5_q)
         WBAE(:,I1)=WTMP3(:)*PP%R%R(:)
!        WRITE (99,*) I1
!        DO I=1,NMAX
!           WRITE (99,*) PP%R%R(I),PP%WAE(I,I1),WBAE(I,I1)
!        ENDDO
      ENDDO
!     CLOSE (99)
!     OPEN (UNIT=99,FILE='ABTMPPS',STATUS="UNKNOWN")
!     DO I1=1,PP%LMAX   !CHANNELS
!        WTMP2(:)=PP%WPS(:,I1)/PP%R%R
!        CALL GRAD_(PP%R,WTMP2,WTMP3)
!        WTMP2=WTMP3*PP%R%R
!        WRITE (99,*) I1
!        DO I=1,NMAX
!           WRITE (99,*) PP%R%R(I),PP%WPS(I,I1),WTMP2(I)
!        ENDDO
!     ENDDO
!     CLOSE (99)

      WEIGHT_WAAE(:,:)=1._q
      WEIGHT_WBAE(:,:)=0._q
      IF (LBONE) THEN
         DO I1=1,PP%LMAX   !CHANNELS
            WTMP2(:)=0._q
            DO I=1,IRMAX
               WTMP2(I)=PP%WAE(I,I1)*PP%WAE(I,I1)
            ENDDO
            CALL SIMPI(PP%R,WTMP2,AAE_CHAR)
            WTMP2(:)=0._q
            DO I=1,IRMAX
               WTMP2(I)=WBAE(I,I1)*WBAE(I,I1)
            ENDDO
            CALL SIMPI(PP%R,WTMP2,BAE_CHAR)
            WEIGHT_TMP=SQRT(AAE_CHAR/(AAE_CHAR+BAE_CHAR))
!           WRITE (*,*) 'charge',AAE_CHAR,BAE_CHAR,WEIGHT_TMP
            DO I=1,IRMAX
               WEIGHT_WAAE(I,I1)=WEIGHT_TMP
               WEIGHT_WBAE(I,I1)=WEIGHT_TMP
            ENDDO
         ENDDO
      ENDIF

! diapara: DO IDP=1,2

! if IDP.eq.1 paramagnetic contribution
! if IDP.eq.2 diamagnetic contribution

!> ----------------------------------------------------------------------------
!> loop i0,i1
!> i0 corresponds to nl
!> i1 corresponds to n'l'
!> LMI combined index (n,l,m)
!> ----------------------------------------------------------------------------

      IF (LZORA) THEN !TOM SOC ZORA>
         !might be duplicate to reset these values, but just to be sure
         WEIGHT_WAAE(:,:)=1._q
         WEIGHT_WBAE(:,:)=0._q

         CALL NMR_ZORA_POT(PP, POT)
         KQ(:)=NMR_ZORA_K(PP%R%NMAX, POT)
      ENDIF !TOM SOC ZORA<

      LM2MAX=0
      L2MAX=0

      LM0=1
      DO I0=1,PP%LMAX   !CHANNELS
      LM1=1
      DO I1=1,PP%LMAX   !CHANNELS

         L0 = PP%LPS(I0)
         L1 = PP%LPS(I1)

         WTMP(:)=PP%WAE(:,I1)/PP%R%R(:)
! perferable maybe to use y d/dr [r*psi(r)] = r * [d/dr psi(r) + psi(r)/r]
         CALL GRAD_(PP%R,WTMP,DWAE)
         DWAE=DWAE*PP%R%R

         WTMP(:)=WBAE(:,I1)/PP%R%R(:)
         CALL GRAD_(PP%R,WTMP,DWBAE)
         DWBAE=DWBAE*PP%R%R

         WTMP(:)=PP%WPS(:,I1)/PP%R%R(:)
         CALL GRAD_(PP%R,WTMP,DWPS)
         DWPS=DWPS*PP%R%R

         DO M1=1,2*L1+1

            LM3=1
            DO L3=0,LMAXNABLA

            CALL YLM3LOOKUP(L0,L3,LMIND1)

            DO M0=1,2*L0+1
            DO M3=1,2*L3+1

               LMIND1 = LMIND1+1
               ISTART1=INDCG(LMIND1)
               IEND1=INDCG(LMIND1+1)

               DO IC1=ISTART1,IEND1-1

                  CGOR=YLM3(IC1)
                  LM2 =JS(IC1)
                  L2  =JL(IC1)
                  IF (LM2 <=0 .OR. LM2 > SIZE(JPAWRAD,2)) THEN
                     CALL vtutor%bug("internal error in CALC_B_MAGATOM_LR: index 2 into JPAWRAD " &
                        // "exceeds bounds", __FILE__, __LINE__)
                  ENDIF
                  LM2MAX=MAX(LM2, LM2MAX)
                  L2MAX=MAX(L2, L2MAX)

                  DO I=1,3

!> ----------------------------------------------------------------------------
!> paramagnetic contribution
!> ----------------------------------------------------------------------------
                     IF (IDP.EQ.1) THEN
                     IF (ABS(YLM_X_YLM(LM3+M3-1,L1*L1+M1,I))>1E-8 .OR. ABS(YLM_NABLA_YLM(LM3+M3-1,L1*L1+M1,I))>1E-8) THEN
                        TMP=0._q
                        SUM=0._q
                        DO K=1,PP%R%NMAX
                           TMP(K) = CGOR * ( &
        &   WEIGHT_WAAE(K,I0) *                                                     WEIGHT_WAAE(K,I1)   * &
        &        PP%WAE(K,I0) * KQ(K) * ( DWAE(K) * YLM_X_YLM(LM3+M3-1,L1*L1+M1,I)                        & !TOM SOC ZORA; note KQ(K)
        &                              + YLM_NABLA_YLM(LM3+M3-1,L1*L1+M1,I)/PP%R%R(K) *  PP%WAE(K,I1) )   &
        & + WEIGHT_WBAE(K,I0) *                                                     WEIGHT_WBAE(K,I1)   * &
        &          WBAE(K,I0) * ( DWBAE(K) * YLM_X_YLM(LM3+M3-1,L1*L1+M1,I)                               &
        &                              + YLM_NABLA_YLM(LM3+M3-1,L1*L1+M1,I)/PP%R%R(K) *    WBAE(K,I1) )   &
        & -    ( PP%WPS(K,I0) * (  DWPS(K) * YLM_X_YLM(LM3+M3-1,L1*L1+M1,I)                               &
        &                              + YLM_NABLA_YLM(LM3+M3-1,L1*L1+M1,I)/PP%R%R(K) *  PP%WPS(K,I1) ))  &
        &                   )
                        ENDDO

!                       IF (LM0+M0-1 <=0 .OR. LM0+M0-1 > SIZE(PP%JPAW,1)) THEN
!                          WRITE(0,*) 'internal error in CALC_JPAW: index 1 into JPAR exceeds bounds'
!                          STOP
!                       ENDIF
!                       IF (LM1+M1-1 <=0 .OR. LM1+M1-1 > SIZE(PP%JPAW,2)) THEN
!                          WRITE(0,*) 'internal error in CALC_JPAW: index 2 into JPAR exceeds bounds'
!                          STOP
!                       ENDIF
!                       PP%JPAW(LM1+M1-1,LM0+M0-1,LM2,I) = PP%JPAW(LM1+M1-1,LM0+M0-1,LM2,I) + SUM

                        !multipy by -i and make symmetric (!)
                        DO K=1,PP%R%NMAX
#ifdef realmode
                           JPAWRAD(K,LM2,I)=JPAWRAD(K,LM2,I)-AIMAG(TMP(K)*(   CONJG(CRHODE(LM1+M1-1,LM0+M0-1))+CRHODE(LM0+M0-1,LM1+M1-1))/2)
#else
                           JPAWRAD(K,LM2,I)=JPAWRAD(K,LM2,I)-AIMAG(TMP(K)*(CHGCONJG(CRHODE(LM1+M1-1,LM0+M0-1))+CRHODE(LM0+M0-1,LM1+M1-1))/2)
#endif
! note that one can save 4 lines above
                        ENDDO
                     ENDIF
                     ENDIF
!> ----------------------------------------------------------------------------
!> diamagnetic contribution
!> ----------------------------------------------------------------------------
                     IF (IDP.EQ.2) THEN
                     IF (ABS(YLM_BxR_YLM(LM3+M3-1,L1*L1+M1,I))>1E-8) THEN
                        TMP=0._q
                        SUM=0._q
!                       |r| from Y*x/y/z*Y
                        DO K=1,PP%R%NMAX
                           TMP(K) = CGOR * PP%R%R(K) * ( &
                               WEIGHT_WAAE(K,I0) * WEIGHT_WAAE(K,I1)                                      * &
                                    PP%WAE(K,I0) *      PP%WAE(K,I1) *( KQ(K) * YLM_BxR_YLM(LM3+M3-1,L1*L1+M1,I) )  &  !TOM SOC ZORA; note KQ(K)
                             + WEIGHT_WBAE(K,I0) * WEIGHT_WBAE(K,I1)                                      * &
                                      WBAE(K,I0) *        WBAE(K,I1) *( YLM_BxR_YLM(LM3+M3-1,L1*L1+M1,I) )  &
                             -      PP%WPS(K,I0) *      PP%WPS(K,I1) *( YLM_BxR_YLM(LM3+M3-1,L1*L1+M1,I) )  &
                            )
                        ENDDO
!                       WRITE (*,'(3I3,2f16.10)') I1,I0,L2,PP%QPAW(I0,I1,L2)

                        DO K=1,PP%R%NMAX
!gK change
!                           JPAWRAD(K,LM2,I)=JPAWRAD(K,LM2,I)+TMP(K)*MOMTOMOM*(-1._q)*CRHODE(LM1+M1-1,LM0+M0-1)
                           JPAWRAD(K,LM2,I)=JPAWRAD(K,LM2,I)+TMP(K)*MOMTOMOM*(-1._q)*CRHODE(LM0+M0-1,LM1+M1-1)
                        ENDDO
                     ENDIF
                     ENDIF
                  ENDDO  ! I

               ENDDO  ! IC1

            ENDDO ! M3
            ENDDO ! M0

            LM3=LM3+2*L3+1
            ENDDO ! L3

         ENDDO ! M1

      LM1=LM1+2*L1+1
      ENDDO ! I1
      LM0=LM0+2*L0+1
      ENDDO ! I0

      IF (IDP==1) THEN
!> ----------------------------------------------------------------------------                                              !NO AUG
!> paramagnetic current augmentation contribution (matches                                                                   !NO AUG
!>  CURRENT_AUGMENTATION in us.F)                                                                                            !NO AUG
!> i.e. depends on whether                                                                                                   !NO AUG
!>                                                                                                                           !NO AUG
!> ~~~~~~~~~~~~~~~~~~~~~                                                                                                     !NO AUG
!>   j_para = e hbar/m_e sum_n Im(psi nabla psi)                                                                             !NO AUG
!> ~~~~~~~~~~~~~~~~~~~~~                                                                                                     !NO AUG
!>                                                                                                                           !NO AUG
!> is augmented with augmentation currents or not                                                                            !NO AUG
!> ----------------------------------------------------------------------------                                              !NO AUG
      IF (LLRAUG) THEN
      LMRHO=1
      DO LRHO=0,2*LMAX
         DO MRHO=1,2*LRHO+1
            DO I=1,3
               DO K=1,PP%R%NMAX
                  JPAWRAD(K,LMRHO+MRHO-1,I)=JPAWRAD(K,LMRHO+MRHO-1,I) +  & ! -1 as above -1 for augmentation correction
                       PP%AUG(K,LRHO)*AIMAG((1.0_q,0.0_q)*JLM(LMRHO+MRHO-1,I)) ! R^2 in PP%AUG ?
               ENDDO
            ENDDO
         ENDDO
         LMRHO=LMRHO+2*LRHO+1
      ENDDO
      ENDIF
!> ----------------------------------------------------------------------------                                              !NO AUG
!> diagmagnetic augmentation correction                                                                                      !NO AUG
!>  depends on whether augmented charge density or plane wave                                                                !NO AUG
!>  contribution only is used when                                                                                           !NO AUG
!>                                                                                                                           !NO AUG
!> ~~~~~~~~~~~~~~~~~~~~~                                                                                                     !NO AUG
!>   j_dia =  e^2/m_e c psi^2 A                                                                                              !NO AUG
!> ~~~~~~~~~~~~~~~~~~~~~                                                                                                     !NO AUG
!>                                                                                                                           !NO AUG
!> is calculated                                                                                                             !NO AUG
!> ----------------------------------------------------------------------------                                              !NO AUG
      ELSE IF (IDP==2) THEN
!NO AUG #ifdef current_augmentation
!NO AUG       ! loop over all L,M pairs in the current density
!NO AUG !      DO LM2=1,(LMAX+LMAXNABLA+1)*(LMAX+LMAXNABLA+1)  ! this does not yet work
!NO AUG       DO LM2=1,(LMAX+LMAX+1)*(LMAX+LMAX+1)
!NO AUG          LMRHO=1
!NO AUG          DO LRHO=0,2*LMAX
!NO AUG          DO MRHO=1,2*LRHO+1
!NO AUG             DO I=1,3
!NO AUG                IF (LM2 <1 .OR. LM2 > SIZE(JPAWRAD,2)) THEN
!NO AUG                   WRITE(0,*) 'internal error CALC_B_MAGATOM, JPAWRAD',LM2,SIZE(JPAWRAD,2)
!NO AUG                   STOP
!NO AUG                ENDIF
!NO AUG                IF (LM2 <1 .OR. LM2 > SIZE(YLM_BxR_YLM,1)) THEN
!NO AUG                   WRITE(0,*) 'internal error CALC_B_MAGATOM, LM2',LM2,SIZE(YLM_BxR_YLM,1)
!NO AUG                   STOP
!NO AUG                ENDIF
!NO AUG                IF (LMRHO+MRHO-1 <1 .OR. LMRHO+MRHO-1 > SIZE(YLM_BxR_YLM,2)) THEN
!NO AUG                   WRITE(0,*) 'internal error CALC_B_MAGATOM, YLM_MxR_YLM, LMRHO+MRHO-1',LMRHO, MRHO
!NO AUG                   STOP
!NO AUG                ENDIF
!NO AUG                IF (LMRHO+MRHO-1 <1 .OR. LMRHO+MRHO-1 > SIZE(RHOLM,1)) THEN
!NO AUG                   WRITE(0,*) 'internal error CALC_B_MAGATOM, RHOLM, LMRHO+MRHO-1',LMRHO, MRHO
!NO AUG                   STOP
!NO AUG                ENDIF
!NO AUG                IF (LRHO <0 .OR. LRHO > SIZE(PP%AUG,2)-1) THEN
!NO AUG                   WRITE(0,*) 'internal error CALC_B_MAGATOM,PPAUG',LRHO,SIZE(PP%AUG,2)
!NO AUG                   STOP
!NO AUG                ENDIF
!NO AUG
!NO AUG                DO K=1,PP%R%NMAX
!NO AUG                   JPAWRAD(K,LM2,I)=JPAWRAD(K,LM2,I) +  &
!NO AUG                      MOMTOMOM*  &                                             ! -1 for units, -1 for augmentation correction
!NO AUG                      PP%AUG(K,LRHO)*RHOLM(LMRHO+MRHO-1)*YLM_BxR_YLM(LM2,LMRHO+MRHO-1,I)*PP%R%R(K)
!NO AUG                ENDDO
!NO AUG             ENDDO
!NO AUG          ENDDO
!NO AUG          LMRHO=LMRHO+2*LRHO+1
!NO AUG          ENDDO
!NO AUG       ENDDO
!NO AUG #endif
       ENDIF

! calculate moment
! \int r \times J_LM^alpha Y_LM dV
! x Y_LM = sqrt(4pi) Y_00 x Y_LM = sqrt(4pi) |r|/|r| (Y_00 x Y_LM) = sqrt(4pi) |r|/|r| YLM_X_YLM(LM=0,LM,x)

      BF=0._q

      DO LM2=1,LM2MAX
         DO I=1,3
            ! please put in some guard against going over bounds

            DO K=1,PP%R%NMAX
               IF (I.EQ.1) OTMP(K)=(YLM_X_YLM(1,LM2,2)*JPAWRAD(K,LM2,3)-YLM_X_YLM(1,LM2,3)*JPAWRAD(K,LM2,2)) *PP%R%R(K)
               IF (I.EQ.2) OTMP(K)=(YLM_X_YLM(1,LM2,3)*JPAWRAD(K,LM2,1)-YLM_X_YLM(1,LM2,1)*JPAWRAD(K,LM2,3)) *PP%R%R(K)
               IF (I.EQ.3) OTMP(K)=(YLM_X_YLM(1,LM2,1)*JPAWRAD(K,LM2,2)-YLM_X_YLM(1,LM2,2)*JPAWRAD(K,LM2,1)) *PP%R%R(K)
               OTMP(K)=OTMP(K)/PP%R%R(K)**3
               !TOM SOC ZORA>
               !Currently KQ(K) is implemented above, which is more inline
               !with the developed formulas.
               !Alternatively, can be included here at OTMP,
               !which is more flexible but less recognizable
               !TOM SOC ZORA<
               !OTMP(K)= KQ(K) * OTMP(K)/PP%R%R(K)**3
            ENDDO

            !ESSENTIALLY A SIMPI
            OSUM=0._q
            DO K=1,PP%R%NMAX
               OSUM=OSUM+OTMP(K)*PP%R%SI(K)
            ENDDO
            BF(I)=BF(I)+OSUM *SQRT(4._q*PI) *2._q !normalisation Y00 * 2

         ENDDO
      ENDDO

!      WRITE (*,'("one center magnetic field*1E6",I3,6E16.7)') NI,REAL(BF,q)*1E6

! default here: calculate dia and paramagnetic contr. seperately
        JPAWRAD=0

! make one-centre contributions to field real
        B_OUT(1)=REAL(BF(1),q)
        B_OUT(2)=REAL(BF(2),q)
        B_OUT(3)=REAL(BF(3),q)

!     ENDDO diapara

      DEALLOCATE(YLM_NABLA_YLM, YLM_X_YLM, YLM_BxR_YLM)
      DEALLOCATE(WTMP,DWAE,DWPS,TMP,OTMP)
      DEALLOCATE(WTMP2,WTMP3,WBAE,DWBAE)
      DEALLOCATE(JPAWRAD)
      DEALLOCATE(POT,KQ)
      DEALLOCATE(WEIGHT_WAAE,WEIGHT_WBAE)

      RETURN
      END SUBROUTINE CALC_B_MAGATOM_LR

!> ************************** SUBROUTINE DEPSUM_ASYM **************************
!>
!> modelled on DEPSUM from us.F, left-right wave functions are different
!> and from different k-points
!>
!> this subroutine calculates  the total one center (on site) "occupancy"
!> matrix of each ll'mm'augmentation channel from the FERMI weights,
!> the weight of each k-point WTKPT and  the wavefunction character
!> of all bands
!> result is stored in CRHODE (thesis gK 10.32)
!>
!> ****************************************************************************
      SUBROUTINE DEPSUM_ASYM(W_LEFT, W_RIGHT ,WDES, LMDIM, CRHODE, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavespin)    W_LEFT
      TYPE (wavespin)    W_RIGHT
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      COMPLEX(q) CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      INTEGER ISP, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      REAL(q) WEIGHT

      IF (.NOT.LOVERL) RETURN
!> ============================================================================
!> initialise to zero
!> ============================================================================
      CRHODE=0
!> ============================================================================
!> loop over all bands and k-points
!> ============================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,NKPTS_ORIG    ! beware original number of k-points
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      band:   DO N=1,WDES%NBANDS

      WEIGHT=WDES%RSPIN*W_LEFT%FERWE(N,NK,ISP)*WDES%WTKPT(NK) ! assumed identical occupancies for left and right channel

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)
      LMBASE_=ISPINOR_*(WDES%NPRO/2)

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 210

      ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

!DIR$ IVDEP
!OCL NOVREC
        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
        DO LP=1,LMMAXC
           CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
              WEIGHT*W_RIGHT%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W_LEFT%CPROJ(LP+LMBASE_,N,NK,ISP))
        ENDDO
        ENDDO

      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_

      ENDDO ion

  210 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
!#ifdef realmode
!      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
!      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
!#else
      CALLMPI( M_sum_z(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_z(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
!#endif

      RETURN
      END SUBROUTINE DEPSUM_ASYM

!> ************************ SUBROUTINE DEPSUM_SYM_ASYM ************************
!>
!> modelled on DEPSUM from us.F, left-right wave functions are different
!> and from different k-points
!>
!> this subroutine calculates  the total one center (on site) "occupancy"
!> applying
!> symmetry to determine the wavefunctions at each k-point
!> it should be fully compatible with the Fock like routines but
!> of course is slightly slower than DEPSUM
!> if this routine is used to construct CRHODE no
!> a posteriori symmetrization of the charge density should be required
!>
!> ****************************************************************************
      SUBROUTINE DEPSUM_SYM_ASYM(W_LEFT, W_RIGHT, WDES, LMDIM, CRHODE, LOVERL, P, LATT_CUR)
      USE prec
      USE sym_prec
      USE wave_high
      USE full_kpoints
      USE pseudo
      USE lattice
      IMPLICIT NONE

      TYPE (wavespin)    W_LEFT, W_RIGHT
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      COMPLEX(q) CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      TYPE (potcar)       P(:)
      TYPE (latt) LATT_CUR

! local
      INTEGER ISP, ISP_IRZ, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      LOGICAL :: LSHIFT
      TYPE( rotation_handle) :: ROT_HANDLE
      TYPE (wavedes1)     WDES1, WDES1_IRZ
      TYPE (wavefun1)     W1_LEFT
      TYPE (wavefun1)     W1_RIGHT
      REAL(q) WEIGHT

      IF (.NOT.LOVERL) RETURN

      CALL CHECK_FULL_KPOINTS
      ROT_HANDLE%LSETUP = .FALSE.

      CALL SETWDES(WDES, WDES1, 0)
      CALL NEWWAV(W1_LEFT  , WDES1, .FALSE.)
      CALL NEWWAV(W1_RIGHT , WDES1, .FALSE.)
!> ============================================================================
!> initialise to zero
!> ============================================================================
      CRHODE=0
!> ============================================================================
!> loop over all bands and k-points
!> ============================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,KPOINTS_FULL%NKPTS       ! BE CAREFUL
         CALL SETWDES(WDES, WDES1, NK)
         CALL SETWDES(WDES, WDES1_IRZ, KPOINTS_FULL%NEQUIV(NK))

         ISP_IRZ=ISP
         IF (KPOINTS_FULL%SPINFLIP(NK)==1) THEN
            ISP_IRZ=3-ISP
         ENDIF

         band:   DO N=1,WDES%NBANDS

            WEIGHT=WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*W_LEFT%FERWE(N,KPOINTS_FULL%NEQUIV(NK),ISP)
                                                     ! assume identical occupancies left and right

            IF (NK <= WDES%NKPTS) THEN
               CALL W1_COPY(ELEMENT(W_LEFT,  WDES1, N, ISP), W1_LEFT )
               CALL W1_COPY(ELEMENT(W_RIGHT, WDES1, N, ISP), W1_RIGHT)
            ELSE
#ifndef gammareal
           !
           ! use symmetry to construct wave function charakter at k
               CALL ROTATE_WAVE_CHARACTER(ROT_HANDLE, P, LATT_CUR, WDES1, ELEMENT(W_LEFT, WDES1_IRZ, N, ISP_IRZ), W1_LEFT )
               CALL ROTATE_WAVE_CHARACTER(ROT_HANDLE, P, LATT_CUR, WDES1, ELEMENT(W_RIGHT, WDES1_IRZ, N, ISP_IRZ), W1_RIGHT)
#endif
            ENDIF

            spinor: DO ISPINOR =0,WDES%NRSPINORS-1
               DO ISPINOR_=0,WDES%NRSPINORS-1

                  LMBASE =ISPINOR *(WDES%NPRO/2)
                  LMBASE_=ISPINOR_*(WDES%NPRO/2)

                  NIS   =1
                  typ:  DO NT=1,WDES%NTYP
                     LMMAXC=WDES%LMMAX(NT)
                     IF (LMMAXC==0) GOTO 210

                     ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
!DIR$ IVDEP
!OCL NOVREC
                        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                           DO LP=1,LMMAXC
                              CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
                                   WEIGHT*W1_RIGHT%CPROJ(L+LMBASE)*GCONJG(W1_LEFT%CPROJ(LP+LMBASE_))
                           ENDDO
                        ENDDO

                        LMBASE = LMMAXC+LMBASE
                        LMBASE_= LMMAXC+LMBASE_

                     ENDDO ion

210                  NIS = NIS+WDES%NITYP(NT)
                  ENDDO typ
               ENDDO
            ENDDO spinor

         ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ*2))
#else
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ*2))
#endif
! note that one can save 4 lines above

      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
      CALL DELWAV(W1_LEFT,  .FALSE.)
      CALL DELWAV(W1_RIGHT, .FALSE.)

      RETURN
      END SUBROUTINE DEPSUM_SYM_ASYM

!> *********************** SUBROUTINE PS_BIOT_SAVART_DR ***********************
!>
!> changes representation current density from reciprocal to direct,
!> brings current to fine grid, entirely unnecessary, but nice for cosmetics,
!> call to Biot-Savart in reciprocal space
!>
!> routine destroys contents input current JVEC_IN,
!> routine returns magnetic fields on atoms B_OUT_PARA
!>
!> modelled from SUBROUTINE CURRENT
!>
!> ****************************************************************************
      SUBROUTINE PS_BIOT_SAVART_DR(JVEC_IN,GRID_SOFT,GRIDC,GRIDUS,&
          C_TO_US,SOFT_TO_C,T_INFO,LATT_CUR,WDES,B_OUT_PARA,&
          CRHODE_AUG,LMDIM,IRDMAX,P,IU0,IU6,&
          BDIR,NICSDATA,ISYM)
      USE prec
      USE lattice
      USE morbitalmag
      USE mgrid
      USE poscar
      USE wave_high
      USE fileio,      only: WRITE_GRID_QUANTITY
      USE constant,    only: TPI, HSQDTM, JDENTOAU
      IMPLICIT NONE
      TYPE (grid_3d)     GRIDC         !< full, fine grid
      TYPE (grid_3d)     GRID_SOFT     !< soft grid for pseudized potentials/ charge etc.
      TYPE (grid_3d) GRIDUS
      TYPE (transit) C_TO_US
      TYPE (transit)     SOFT_TO_C
      TYPE (latt)        LATT_CUR      !< lattice
      TYPE (type_info):: T_INFO        !< type information
      TYPE (wavedes)     WDES          !< descriptor for wavefunction
      TYPE (potcar), TARGET :: P(:)
      COMPLEX(q) :: B_OUT_PARA(3,T_INFO%NIONS)   !< output field
      COMPLEX(q) CRHODE_AUG(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) !< augmentation input
      RGRID :: JVEC_IN(WDES%GRID%RL%NP,3)
      INTEGER LMDIM,IRDMAX
      INTEGER IU0, IU6
      INTEGER BDIR
      COMPLEX(q) :: NICSDATA(:,:)
      INTEGER ISYM
      ! local
      COMPLEX(q) TMP(3)
      COMPLEX(q), ALLOCATABLE :: JVEC(:)
      COMPLEX(q), ALLOCATABLE :: JTOT(:,:)
      REAL(q) DISPL_NIONS(3,T_INFO%NIONS)
      INTEGER M, I, J
      CHARACTER(2) :: DIR_TEXT(3)=(/"BX","BY","BZ"/)

!> ----------------------------------------------------------------------------
!> change from reciprocal space presentation to cartesian coordinates
!> sign change because electronic charge is negative
!> ----------------------------------------------------------------------------
      DO M=1,WDES%GRID%RL%NP
         DO I=1,3
            TMP(I)=JVEC_IN(M,I) ! make real from ... well
         ENDDO
         CALL CDIRKAR(1,TMP,LATT_CUR%B)
         DO I=1,3
            JVEC_IN(M,I)=TMP(I)
         ENDDO
      ENDDO

!> ----------------------------------------------------------------------------
!> bring the current density to the fine grid
!> ----------------------------------------------------------------------------

      ALLOCATE(JVEC(GRID_SOFT%MPLWV))
      ALLOCATE(JTOT(GRIDC%MPLWV,3))

      JTOT=0

      DO I=1,3
! now merge the current from all nodes
#ifdef realmode
         CALLMPI( M_sum_d(WDES%COMM_INTER, JVEC_IN(1,I), WDES%GRID%RL%NP))
         CALLMPI( M_sum_d(WDES%COMM_KINTER,JVEC_IN(1,I), WDES%GRID%RL%NP))
#else
         CALLMPI( M_sum_z(WDES%COMM_INTER, JVEC_IN(1,I), WDES%GRID%RL%NP))
         CALLMPI( M_sum_z(WDES%COMM_KINTER,JVEC_IN(1,I), WDES%GRID%RL%NP))
#endif

         CALL FFT_RC_SCALE(JVEC_IN(1,I),JVEC(1),GRID_SOFT)
         ! set the current density of unbalanced lattic-vectors to 0
         CALL SETUNB(JVEC(1),GRID_SOFT)

         ! bring to full grid
         IF (.NOT.WDES%LOVERL) THEN
            CALL RC_ADD(JVEC(1),1.0_q,JVEC(1),0.0_q,JTOT(1,I),GRID_SOFT)
         ELSE
           CALL ADD_GRID(GRIDC, GRID_SOFT, SOFT_TO_C,JVEC(1),JTOT(1,I))
           CALL SETUNB_COMPAT(JTOT(1,I),GRIDC)
        ENDIF
      ENDDO

      IF (LLRAUG) THEN
         DISPL_NIONS=0
         CRHODE_AUG(:,:,:,:)=-CRHODE_AUG(:,:,:,:) !GIL empirical test
         CALL CURRENT_AUGMENTATION_CMPLX( &
              WDES, GRIDC, GRIDUS, C_TO_US, LATT_CUR, P, T_INFO, &
              LMDIM, CRHODE_AUG(1,1,1,1), JTOT, IRDMAX, DISPL_NIONS, .TRUE.)
      ENDIF

      ! bring total pseudo current and total charge to real space
      DO I=1,3
         CALL FFT3D(JTOT(1,I),GRIDC,1)
      ENDDO

      IF (WRT_NMRCUR==BDIR.OR.WRT_NMRCUR==4) THEN
         ! convert JTOT to a.u. and write to file
         CALL WRITE_GRID_QUANTITY("NMRCUR"//DIR_TEXT(BDIR),WDES,&
          T_INFO%POSION,T_INFO,LATT_CUR,GRIDC,JDENTOAU*HSQDTM*2_q*TPI*JTOT/LATT_CUR%OMEGA)
      ENDIF

!> ----------------------------------------------------------------------------
!> call Biot-Savart
!> ----------------------------------------------------------------------------

      CALL PS_BFIELD_C(JTOT(1,1),GRIDC,LATT_CUR,T_INFO,NPOSROT,POSROT,&
          B_OUT_PARA,IU0,IU6,NICSDATA)

      DEALLOCATE(JTOT)
      DEALLOCATE(JVEC)

      END SUBROUTINE PS_BIOT_SAVART_DR

!>******************************* FUNCTION EIJK *******************************
!>
!>*****************************************************************************
      FUNCTION EIJK(I,J,K)
      USE prec
      INTEGER I,J,K
      REAL(q) EIJK
      EIJK=REAL((J-I)*(K-I)*(K-J)/2,q)
      END FUNCTION EIJK

!> *************************** SUBROUTINE NABIJ_VGV ***************************
!>
!> Modelled from NABIJ_SOFT from optics.F
!>
!> ****************************************************************************
      SUBROUTINE NABIJ_VGV(W, WGKPQ, W_1, W_2, W_3, GRID, LATT_CUR, NKPT, INFO, NAB_QQ)
      USE prec
      USE lattice
      USE mgrid
      USE wave
      USE constant
      USE base
      IMPLICIT NONE

      TYPE (wavespin) W
      TYPE (wavespin) WGKPQ ! G u_i v_k+q, k| u>
      TYPE (wavespin) W_1
      TYPE (wavespin) W_2
      TYPE (wavespin) W_3
      TYPE (grid_3d) GRID
      TYPE (latt)    LATT_CUR
      TYPE (info_struct) INFO
      COMPLEX(q)  NAB_QQ(3)     !not GDEF
      INTEGER NKPT

      !local
      INTEGER I, NB, NBP, IDIR, NK, ISP
      REAL(q) G1, G2, G3, GC, WEIGHT
      COMPLEX(q) CSUM

      NAB_QQ = 0

! loop over all (current) k-points, spin component and cartesian directions
      DO IDIR=1,3

      DO NK=1,NKPT
#ifdef MPI
         IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         DO ISP=1,INFO%ISPIN
            DO NB =1,W%WDES%NBANDS
               WEIGHT=W%WDES%RSPIN*W%FERWE(NB,NK,ISP)*W%WDES%WTKPT(NK)     ! occupancies assumed identical
               IF (IDIR.EQ.1) THEN   ! avoid IF statement in next loop
                  DO I  =1,W%WDES%NPLWKP(NK)
!                    G1=W%WDES%IGX(I,NK)+W%WDES%VKPT(1,NK)   ! beware take k-points from W.
!                    G2=W%WDES%IGY(I,NK)+W%WDES%VKPT(2,NK)
!                    G3=W%WDES%IGZ(I,NK)+W%WDES%VKPT(3,NK)
!                    GC=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))*TPI

                     CSUM=WGKPQ%CW(I,NB,NK,ISP)*CONJG(W_1%CW(I,NB,NK,ISP))
!                    CSUM=WGKPQ%CW(I,NB,NK,ISP)*GC*CONJG(W%CW(I,NB,NK,ISP)) !*(0._q,1._q)  ! Gilles doesn't understand the i
                     CSUM = CSUM /(2*HSQDTM)   ! heuristic, but should be correct for units, see MLR_COMMUTATOR
                     NAB_QQ(IDIR)=NAB_QQ(IDIR)+CSUM*WEIGHT
                  ENDDO
               ENDIF
               IF (IDIR.EQ.2) THEN
                  DO I  =1,W%WDES%NPLWKP(NK)
                     CSUM=WGKPQ%CW(I,NB,NK,ISP)*CONJG(W_2%CW(I,NB,NK,ISP))
                     CSUM = CSUM /(2*HSQDTM)   ! heuristic, but should be correct for units, see MLR_COMMUTATOR
                     NAB_QQ(IDIR)=NAB_QQ(IDIR)+CSUM*WEIGHT
                  ENDDO
               ENDIF
               IF (IDIR.EQ.3) THEN
                  DO I  =1,W%WDES%NPLWKP(NK)
                     CSUM=WGKPQ%CW(I,NB,NK,ISP)*CONJG(W_3%CW(I,NB,NK,ISP))
                     CSUM = CSUM /(2*HSQDTM)   ! heuristic, but should be correct for units, see MLR_COMMUTATOR
                     NAB_QQ(IDIR)=NAB_QQ(IDIR)+CSUM*WEIGHT
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ENDDO

!     CALLMPI(M_sum_z(W%WDES%COMM_KINTER,NAB_QQ(1),3))

      RETURN
      END SUBROUTINE NABIJ_VGV

!> ************************** SUBROUTINE NABIJ_ASYM ***************************
!>
!>  Modelled from NABIJ_SOFT from optics.F
!>
!> ****************************************************************************
      SUBROUTINE NABIJ_ASYM(W, WGKPQ, GRID, LATT_CUR, NKPT, INFO, NAB_QQ)
      USE prec
      USE lattice
      USE mgrid
      USE wave
      USE constant
      USE base
      IMPLICIT NONE

      TYPE (wavespin) W
      TYPE (wavespin) WGKPQ ! G u_i v_k+q, k| u>
      TYPE (grid_3d) GRID
      TYPE (latt)    LATT_CUR
      TYPE (info_struct) INFO
      COMPLEX(q)  NAB_QQ(3)     !not GDEF
      INTEGER NKPT

      !local
      INTEGER I, NB, NBP, IDIR, NK, ISP
      INTEGER ISPINOR, II
      REAL(q) G1, G2, G3, GC, WEIGHT
      COMPLEX(q) CSUM

      NAB_QQ = 0

! loop over all (current) k-points, spin component and cartesian directions
      DO IDIR=1,3

      DO NK=1,NKPT
#ifdef MPI
         IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         DO ISP=1,INFO%ISPIN
            DO NB =1,W%WDES%NBANDS
               WEIGHT=W%WDES%RSPIN*W%FERWE(NB,NK,ISP)*W%WDES%WTKPT(NK)     ! occupancies assumed identical
               DO I  =1,W%WDES%NGVECTOR(NK)
                  G1=W%WDES%IGX(I,NK)+W%WDES%VKPT(1,NK)   ! beware take k-points from W.
                  G2=W%WDES%IGY(I,NK)+W%WDES%VKPT(2,NK)
                  G3=W%WDES%IGZ(I,NK)+W%WDES%VKPT(3,NK)
                  GC=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))*TPI
                  DO ISPINOR=0,W%WDES%NRSPINORS-1
                     II=I+W%WDES%NGVECTOR(NK)*ISPINOR
                     CSUM=WGKPQ%CW(II,NB,NK,ISP)*GC*CONJG(W%CW(II,NB,NK,ISP))  !*(0._q,1._q)  ! Gilles doesn't understand the i
                     NAB_QQ(IDIR)=NAB_QQ(IDIR)+CSUM*WEIGHT
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ENDDO

!     CALLMPI(M_sum_z(W%WDES%COMM_KINTER,NAB_QQ(1),3))

      RETURN
      END SUBROUTINE NABIJ_ASYM

!> *************************** MLR_CALC_SPINORB_ONE ***************************
!>
!> Calculates V^{nl}_{\perp, R}
!> Taken from CALC_B_MAGATOM_LR, but largely rewritten and much simplified
!>
!> This is all terribly tricky: the routine gives a quite small contribution.
!> No way to check it is correct. Moreover, units are brain damaging, see below.
!>
!> Subroutine is replace by _ALTS and _ALTD, but contains valuable comments
!>
!> ****************************************************************************
      SUBROUTINE MLR_CALC_SPINORB_ONE(WDES,P,T_INFO,LMDIM,BDIR,CLQIJ)
      USE pseudo
      USE constant
      USE asa
      USE radial
      USE wave
      USE paw
      USE poscar
      USE prec
      USE morbitalmag, ONLY : LBONE !TOM 20240226, test whether this can be removed

      IMPLICIT NONE

      TYPE (wavedes) WDES
      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
!     COMPLEX(q) :: CRHODE(:,:)         ! one center occupancy matrix
      INTEGER :: BDIR,LMDIM
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      INTEGER, EXTERNAL :: MAXL1

      !local variables
      TYPE (potcar), POINTER :: PP   ! pseudopotential descriptor
      INTEGER :: NIP,NI
      INTEGER CHANNELS          ! number of channels
      INTEGER NMAX,I,NT
      INTEGER I0,I1
      INTEGER L0,L1
      INTEGER M0,M1
      INTEGER LM0,LM1
      INTEGER ISTART1,IEND1,LMIND1,IC1
      INTEGER LMAXNABLA,LMAX,K !,IDP
      INTEGER II1,II2,II3,LL,M
      INTEGER LRHO,LMRHO,MRHO
      INTEGER IRMAX
      INTEGER CH1,CH2,LM,LMP,LLP,J,MP
      INTEGER LMCOMB2
      REAL(q) :: CGOR
      REAL(q) :: SUM
      REAL(q), ALLOCATABLE :: TMP(:)
      REAL(q), ALLOCATABLE :: YLM_NABLA_YLM(:,:,:),YLM_X_YLM(:,:,:),YLM_BxR_YLM(:,:,:)
      REAL(q) :: BCONST_LR(3)

      REAL(q), PARAMETER :: C = 137.037  ! speed of light in a.u.
      REAL(q), PARAMETER :: INVMC2=7.45596E-6
!>                          INVMC2=hbar^2/2(m_e c)^2 in A^2
!>                          appears to be AUTOA**2/(2.d0*(CLIGHT)**2),
!>                          that yields: 7.45584116E-6
      REAL(q) UNITCONV
      REAL(q), ALLOCATABLE :: POT(:),KQ(:),KSI(:),APOT(:),DPOT(:)
      REAL(q) RDUM,SCALE
      REAL(q) VTMP(3)
      COMPLEX(q) SIGDOTVEC(4)

!> ----------------------------------------------------------------------------
!> ~~~~~~~~~~~~~~~~~~~~~
!> hbar/m^2  K^2/(4c^2 r) dV/dr r \times (sigma \times r)
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!>  all the r's cancel --> 1
!>  K is dimensionless --> 1
!>  we split off hbar/m (so we have the same units as in the perturbations on the k-point (L_R Q_R, etc.)):
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!> hbar/m [ V/(4 m c^2) ]
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!>  the quantity in [brackets] is dimensionless, but we have to consider V in eV. Hence:
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!> [...] = 1/(2 rytoev)/4/c^2
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!>  it can be checked with:
!>
!>  m_e =  9.10938356d-31 !kg
!>  c = 299792458.d0 ! m/sec
!>  ev = 1.602176635d-19 ! J = kg m^2/s^2
!>  [...] =  ev/m_e/c/c/4.d0 = 4.8923780379836913E-007 ! kg m^2/s^2 / kg / (m^2/s^2)
!> ----------------------------------------------------------------------------
      UNITCONV = 1._q/RYTOEV/2._q/4._q/CLIGHT/CLIGHT

!     IF (.NOT. ASSOCIATED(PP%JPAW)) RETURN
      BCONST_LR=0._q
      BCONST_LR(BDIR)=1._q

      ion: DO NI=1,T_INFO%NIONS
         NT=T_INFO%ITYP(NI)

         ! cycle if this ion is not treated by this node
         NIP=NI_LOCAL(NI,WDES%COMM_INB)
         IF (NIP==0) CYCLE ion

         PP=>PP_POINTER(P, NI, NT)

         ! allocate
         NMAX=PP%R%NMAX
         ALLOCATE(TMP(NMAX))
         ALLOCATE(KSI(NMAX),APOT(NMAX),DPOT(NMAX))
         ALLOCATE(POT(NMAX),KQ(NMAX))

         CHANNELS = PP%LMAX
         LMAX=MAXL1(PP)
         LMAXNABLA=MAXL1(PP)+1   ! nabla Y_lm has components up to LMAX+1

         ALLOCATE(YLM_NABLA_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAXNABLA+2)*(2*LMAXNABLA+2),0:3))
         ALLOCATE(YLM_X_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAX+2)*(2*LMAX+2),0:3))
         ALLOCATE(YLM_BxR_YLM((2*LMAX+2)*(2*LMAX+2),(2*LMAX+2)*(2*LMAX+2),0:3))

         YLM_NABLA_YLM=0._q
         YLM_X_YLM=0._q
         YLM_BxR_YLM=0._q

!GIL CALL SETYLM_NABLA_YLM(LMAXNABLA,YLM_NABLA_YLM,YLM_X_YLM)
!> in principle we should go up to LMAX+LMAXNABLA+1 = 2 LMAX + 2
!> see below, but that does not work for f electrons
         CALL SETYLM_NABLA_YLM(2*LMAX+1,YLM_NABLA_YLM,YLM_X_YLM)

         YLM_BxR_YLM(:,:,1)=BCONST_LR(2)*YLM_X_YLM(:,:,3)-BCONST_LR(3)*YLM_X_YLM(:,:,2)
         YLM_BxR_YLM(:,:,2)=BCONST_LR(3)*YLM_X_YLM(:,:,1)-BCONST_LR(1)*YLM_X_YLM(:,:,3)
         YLM_BxR_YLM(:,:,3)=BCONST_LR(1)*YLM_X_YLM(:,:,2)-BCONST_LR(2)*YLM_X_YLM(:,:,1)
         YLM_BxR_YLM(:,:,:)=YLM_BxR_YLM(:,:,:)/2. !*MAGMOMTOENERGY

!        IF (LLRAUG) THEN
!           DO I=1,3
!              CALL CALC_RHOLM_JVEC_CMPLX(  CRHODE, JLM(:,I), PP, I)
!           ENDDO
!        ENDIF

!> The radial ZORA factor from relativistic.F,
!> for comments, see relativstic.F and LBONE in linear_response_NMR.F.
         SCALE=1/(2*SQRT(PI))
!> Hartree potential of core electrons
         CALL RAD_POT_HAR(0, PP%R, APOT, PP%RHOAE, SUM)
!> add potential of total ionic core
         APOT=APOT*SCALE - (PP%ZCORE+PP%ZVALF_ORIG)*FELECT/PP%R%R    ! how do we know FELECT?

!> In paw.F valence Hartree and total XC potential are added to PP%POTAE_XCUPDATED and put in POTAE.
!> In order to add the valence Hartree and XC potential (of allelectrons) one has to
!> rescale POTAE with 1/2/sqrt(2) and subtract PP%POTAE_XCUPDATED in the
!> call to SPINORB_STRENGTH. That would give here:
!>        APOT=APOT-PP%POTAE_XCUPDATED
!>        APOT=APOT+PP%POTAE(1:SIZE(APOT),1,1)* SCALE
!> But here we have a different situation.
!>
!> Let's assume we have what is set in SET_PAW_ATOM_POT, then the scaling
!> is already done, and the sign has been inverted. When then do:
         POT=POT-PP%POTAE
!> That has added the Hartree potential of valence electrons and xc potential of all electrons.

!> Next is the alignment.
         POT=POT+(PP%VPSRMAX-POT(NMAX))
!> Now the setup of POT is exactly as with LBONE.

!> TOM TESTING: in principle, following two lines achieve the same as the above
!>         APOT = 0._q
!>         CALL NMR_ZORA_POT(PP, APOT)
!>TOM TESTING

         CALL GRAD(PP%R,APOT,DPOT)
         DO I=1,PP%R%NMAX
                    ! option 1: as in /vol/thchem/dewijs/vasp_src/vasp.6/check/vasp6_guppy_clean/build/ncl
                   KSI(I)=UNITCONV*(1._q/(1._q-APOT(I)/(2._q*RYTOEV)/C/C/2._q)**2)*DPOT(I)/PP%R%R(I)
                    ! option 2: as in /vol/thchem/dewijs/vasp_src/vasp.6/vasp6_guppy_clean/build/ncl/linear_response_NMR.F-vanbackup
!                    KSI(I)=INVMC2*(1._q/(1._q-APOT(I)/(2._q*RYTOEV)/C/C/2._q)**2)*DPOT(I)/PP%R%R(I)
!           write (*,*) 'KSI',I,KSI(I)
         ENDDO

!> ----------------------------------------------------------------------------
!> loop i0,i1
!> i0 corresponds to nl
!> i1 corresponds to n'l'
!> LMI combined index (n,l,m)
!> ----------------------------------------------------------------------------
         LM0=1
         DO I0=1,PP%LMAX !< CHANNELS
         LM1=1
         DO I1=1,PP%LMAX !< CHANNELS

            L0 = PP%LPS(I0)
            L1 = PP%LPS(I1)

            DO M1=1,2*L1+1
            DO M0=1,2*L0+1

               VTMP(:)=0._q
               DO LMCOMB2=1,(LMAX+1)*(LMAX+1) !,2*LMAX+1
                  VTMP(1)=VTMP(1)+YLM_X_YLM(L0*L0+M0,LMCOMB2,2)*YLM_BxR_YLM(LMCOMB2,L1*L1+M1,3)-YLM_X_YLM(L0*L0+M0,LMCOMB2,3)*YLM_BxR_YLM(LMCOMB2,L1*L1+M1,2)
                  VTMP(2)=VTMP(2)+YLM_X_YLM(L0*L0+M0,LMCOMB2,3)*YLM_BxR_YLM(LMCOMB2,L1*L1+M1,1)-YLM_X_YLM(L0*L0+M0,LMCOMB2,1)*YLM_BxR_YLM(LMCOMB2,L1*L1+M1,3)
                  VTMP(3)=VTMP(3)+YLM_X_YLM(L0*L0+M0,LMCOMB2,1)*YLM_BxR_YLM(LMCOMB2,L1*L1+M1,2)-YLM_X_YLM(L0*L0+M0,LMCOMB2,2)*YLM_BxR_YLM(LMCOMB2,L1*L1+M1,1)
               ENDDO

            !  IF (ABS(YLM_BxR_YLM(LM3+M3-1,L1*L1+M1,I))>1E-8) THEN
               TMP=0._q
               SUM=0._q
!              twice |r| from Y*x/y/z*Y
               DO K=1,PP%R%NMAX
                  TMP(K) = PP%R%R(K)**2 * KSI(K) * PP%WAE(K,I0) * PP%WAE(K,I1)  ! -  PP%WPS(K,I0) * PP%WPS(K,I1) !TOM 23-02-2024, no PS here indeed
               ENDDO
!              DO K=1,PP%R%NMAX
!                 JPAWRAD(K,LM2,I)=TMP(K)*(-1._q) !*MOMTOMOM*(-1._q)          !*CRHODE(LM0+M0-1,LM1+M1-1)
!              ENDDO
!              ENDIF
               DO K=1,PP%R%NMAX
                  SUM=SUM+TMP(K)*PP%R%SI(K)    ! pas op: indices
               ENDDO
               VTMP(:)=VTMP(:)*SUM

               ! modelled on relativistic.F
               SIGDOTVEC(1)=  (1._q,0._q)*VTMP(3)
               SIGDOTVEC(2)=  (1._q,0._q)*VTMP(1) + (0._q,1._q)*VTMP(2)  !,LM0+M0-1,LM1+M1-1)
               SIGDOTVEC(3)=  (1._q,0._q)*VTMP(1) - (0._q,1._q)*VTMP(2)  !,LM0+M0-1,LM1+M1-1)
               SIGDOTVEC(4)= -(1._q,0._q)*VTMP(3)

               DO I=0,1
               DO J=0,1
                  ! option 1: as in /vol/thchem/dewijs/vasp_src/vasp.6/check/vasp6_guppy_clean/build/ncl
!                 CLQIJ(LM0+M0-1,LM1+M1-1,NIP,J+2*I+1)=CLQIJ(LM0+M0-1,LM1+M1-1,NIP,J+2*I+1) + SIGDOTVEC(I+2*J+1) *(-2._q)
                  ! option 2: as in /vol/thchem/dewijs/vasp_src/vasp.6/vasp6_guppy_clean/build/ncl/linear_response_NMR.F-vanbackup
                  CLQIJ(LM0+M0-1,LM1+M1-1,NIP,J+2*I+1)=CLQIJ(LM0+M0-1,LM1+M1-1,NIP,J+2*I+1) + SIGDOTVEC(I+2*J+1) !GILT *(-2._q)/2._q !GILT    !deb
               ENDDO
               ENDDO

            ENDDO ! M1
            ENDDO ! M0

         LM1=LM1+2*L1+1
         ENDDO ! I1
         LM0=LM0+2*L0+1
         ENDDO ! I0

!        DO CH1=1,PP%LMAX
!           LL =PP%LPS(CH1)
!           WRITE (*,*) CH1,LL,'CH1'
!        ENDDO

         DEALLOCATE(YLM_NABLA_YLM, YLM_X_YLM, YLM_BxR_YLM)
         DEALLOCATE(TMP)
         DEALLOCATE(POT,KQ,KSI,APOT,DPOT)

      ENDDO ion

      RETURN
      END SUBROUTINE MLR_CALC_SPINORB_ONE

!> ******************* SUBROUTINE MLR_CALC_SPINORB_ONE_ALTS *******************
!>
!> Routine based on J.Chem.Phys., Vol. 110, No. 16, 22 april 1999 eq. 18
!> (i.e. ADF, Wolff ZORA)
!> routine covers first of the last two terms (the first order in B SOC terms)
!> (eq. 28 arXiv:2504.17380 [cond-mat.mtrl-sci])
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!> \sigma \cdot B ( r \cdot \nabla( (K-1)/4c ) )
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> together with MLR_CALC_SPINORB_ONE_ALTD, this replaces MLR_CALC_SPINORB_ONE
!>
!> difficult subroutine to verify because of its small effect.
!>
!> ****************************************************************************
      SUBROUTINE MLR_CALC_SPINORB_ONE_ALTS(WDES,P,T_INFO,LMDIM,IDIR,CLQIJ)
      USE pseudo
      USE constant
      USE asa
      USE radial
      USE wave
      USE paw
      USE poscar
      USE prec

      IMPLICIT NONE

      !< INPUT/OUTPUT
      TYPE (wavedes) WDES
      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      INTEGER LMDIM, IDIR
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      !< LOCAL VARIABLES
      REAL(q) BDIR_ONE(3)
      TYPE (potcar), POINTER :: PP
      INTEGER I,J,K
      INTEGER NI,NT,NIP
      INTEGER LMAX,LMAXNABLA,LYMAX,LDIMP,L,LL,LLP
      INTEGER LM,LMP,M,MP,CH1,CH2
      REAL(q) SUM

      REAL(q), ALLOCATABLE :: ZPOT(:),KPOT(:),KSI(:)
      COMPLEX(q) BDOTSIG(4), TMP(4)

      REAL(q) UNITCONV
      REAL(q), PARAMETER :: C = 137.037  ! speed of light in a.u.
      REAL(q), PARAMETER :: INVMC2=7.45596E-6

      INTEGER, EXTERNAL :: MAXL1

      UNITCONV = 1._q/RYTOEV/2._q/4._q/CLIGHT/CLIGHT

      BDIR_ONE=0._q
      BDIR_ONE(IDIR)=1._q

      ion: DO NI=1,T_INFO%NIONS
         NT=T_INFO%ITYP(NI)

         !cycle if ion not treated by node
         NIP=NI_LOCAL(NI,WDES%COMM_INB)
         IF (NIP==0) CYCLE ion

         PP => PP_POINTER(P,NI,T_INFO%ITYP(NI))

         LDIMP=3                      ! ??
         IF (LYMAX>2) LDIMP=LYMAX+1   ! ??

         !> PREP BDOTSIG
         !> Modelled on relativistic.F
         !< <up| B \cdot \sigma |up>
         BDOTSIG(1)=  (1._q,0._q)*BDIR_ONE(3)
         !< <up| B \cdot \sigma |down>
         BDOTSIG(2)=  (1._q,0._q)*BDIR_ONE(1) + (0._q,1._q)*BDIR_ONE(2)
         !< <down| B \cdot \sigma |up>
         BDOTSIG(3)=  (1._q,0._q)*BDIR_ONE(1) - (0._q,1._q)*BDIR_ONE(2)
         !< <down| B \cdot \sigma |down>
         BDOTSIG(4)= -(1._q,0._q)*BDIR_ONE(3)

         LMAX=MAXL1(PP)
         LMAXNABLA=MAXL1(PP)+1   !< nabla Y_lm has components up to LMAX+1

         !< CALCULATE GRADIENT OF ZORA K
         ALLOCATE(ZPOT(PP%R%NMAX),KPOT(PP%R%NMAX),KSI(PP%R%NMAX))
         ZPOT = 0._q
         CALL NMR_ZORA_POT(PP, ZPOT)

!> Option 1: First take the gradient of the potential, then build K
!         CALL GRAD(PP%R,ZPOT,KPOT)
!
!         DO I=1,PP%R%NMAX
!            ! option 1: as in /vol/thchem/dewijs/vasp_src/vasp.6/check/vasp6_guppy_clean/build/ncl
!!           KSI(I)=UNITCONV*(1._q/(1._q-ZPOT(I)/(2._q*RYTOEV)/C/C/2._q)**2)*KPOT(I)/PP%R%R(I)
            ! option 2: as in /vol/thchem/dewijs/vasp_src/vasp.6/vasp6_guppy_clean/build/ncl/linear_response_NMR.F-vanbackup
!            KSI(I)=INVMC2*(1._q/(1._q-ZPOT(I)/(2._q*RYTOEV)/C/C/2._q))**2/PP%R%R(I)*KPOT(I)
!!           write (*,*) 'KSI',I,KSI(I)
!         ENDDO
!> Option 2: First build K, then take gradient (of both K and potential).
         KPOT(:) = NMR_ZORA_K(PP%R%NMAX, ZPOT)/2._q
         CALL GRAD(PP%R,KPOT,KSI)

         LM=1
         DO CH1=1,PP%LMAX
         LMP=1
         DO CH2=1,PP%LMAX

            !integrate "SIMPI"
            SUM = 0._q
            DO K=1, PP%R%NMAX
               !here we multiply with PP%R%R to make "r \cdot \nabla K"
               SUM = SUM + PP%R%SI(K) * PP%R%R(K) * KSI(K) * PP%WAE(K,CH1) * PP%WAE(K,CH2)
            ENDDO

            !< \sigma \cdot B (r \cdot \nabla K)
            TMP(:)=BDOTSIG(:)*SUM

            LL=PP%LPS(CH1)
            LLP=PP%LPS(CH2)
            IF (LL == LLP .AND. LL<=LDIMP) THEN
               DO I=0,1
               DO J=0,1
                  DO MP=1,2*LLP+1
                     CLQIJ(LMP+MP-1,LM+MP-1,NIP,J+2*I+1)=CLQIJ(LMP+MP-1,LM+MP-1,NIP,J+2*I+1)+TMP(I+2*J+1)
                  ENDDO
               ENDDO
               ENDDO
            ENDIF

!            IF (LL == LLP .AND. LL>0 .AND. LL<=LMAX ) THEN
!               DO M=1,2*LL+1
!               DO MP=1,2*LLP+1
!
!                  !modelled on relativistic.F (loop for sigma)
!                  DO I=0,1
!                  DO J=0,1
!                     CLQIJ(LMP+MP-1,LM+M-1,NIP,J+2*I+1)=CLQIJ(LMP+MP-1,LM+M-1,NIP,J+2*I+1)+TMP(I+2*J+1)
!                  ENDDO
!                  ENDDO
!
!               ENDDO
!               ENDDO
!            ENDIF

         LMP=LMP+(2*LLP+1)
         ENDDO
         LM= LM+ (2*LL+1)
         ENDDO

         DEALLOCATE(ZPOT,KPOT,KSI)
      ENDDO ion

      RETURN
      END SUBROUTINE MLR_CALC_SPINORB_ONE_ALTS

!> ******************* SUBROUTINE MLR_CALC_SPINORB_ONE_ALTD *******************
!>
!> Routine based on J.Chem.Phys., Vol. 110, No. 16, 22 april 1999 eq. 18
!> (i.e. ADF, Wolff ZORA)
!> routine covers second of the last two terms (the first order in B SOC terms)
!> (eq. 29 arXiv:2504.17380 [cond-mat.mtrl-sci])
!>
!> ~~~~~~~~~~~~~~~~~~~~~
!> - \sigma \cdot r ( B \cdot \nabla( (K-1)/4c ) )
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> together with MLR_CALC_SPINORB_ONE_ALTS, this replaces MLR_CALC_SPINORB_ONE
!>
!> difficult subroutine to verify because of its small effect.
!>
!> ****************************************************************************
      SUBROUTINE MLR_CALC_SPINORB_ONE_ALTD(WDES,P,T_INFO,LMDIM,IDIR,CLQIJ)
      USE pseudo
      USE constant
      USE asa
      USE radial
      USE wave
      USE paw
      USE poscar
      USE prec

      IMPLICIT NONE

      !< INPUT/OUTPUT
      TYPE (wavedes) WDES
      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      INTEGER LMDIM, IDIR
      COMPLEX(q) CLQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      !< LOCAL VARIABLES
      REAL(q) BDIR_ONE(3)
      TYPE (potcar), POINTER :: PP
      INTEGER I,J,K
      INTEGER NI,NT,NIP
      INTEGER LMAX,LMMAX,LYMAX,LDIMP,L,LL,LLP
      INTEGER LM,LMP,M,MP,CH1,CH2
      REAL(q) SUM, BDK

      REAL(q), ALLOCATABLE :: YLM_XIXJ_YLM(:,:,:)

      REAL(q), ALLOCATABLE :: ZPOT(:),KPOT(:),KSI(:)
      COMPLEX(q) SIG(4)

      REAL(q) UNITCONV
      REAL(q), PARAMETER :: C = 137.037  ! speed of light in a.u.
      REAL(q), PARAMETER :: INVMC2=7.45596E-6

      INTEGER, EXTERNAL :: MAXL1

      UNITCONV = 1._q/RYTOEV/2._q/4._q/CLIGHT/CLIGHT

      BDIR_ONE=0._q
      BDIR_ONE(IDIR)=1._q

      ion: DO NI=1,T_INFO%NIONS
         NT=T_INFO%ITYP(NI)

         !cycle if ion not treated by node
         NIP=NI_LOCAL(NI,WDES%COMM_INB)
         IF (NIP==0) CYCLE ion

         PP => PP_POINTER(P,NI,T_INFO%ITYP(NI))

         LDIMP=3                      ! ??
         IF (LYMAX>2) LDIMP=LYMAX+1   ! ??

         !SEE DAB_HEAVYSIDE (dmatrix.F)
         LMAX=MAXL1(PP)
         LMMAX=(LMAX+1)*(LMAX+1)

         ALLOCATE(YLM_XIXJ_YLM(LMMAX,LMMAX,6))
         YLM_XIXJ_YLM(:,:,:) = 0._q
         CALL SETYLM_XIXJ_YLM(LMAX,YLM_XIXJ_YLM)

         !< CALCULATE GRADIENT OF ZORA K
         ALLOCATE(ZPOT(PP%R%NMAX),KPOT(PP%R%NMAX),KSI(PP%R%NMAX))
         ZPOT = 0._q
         CALL NMR_ZORA_POT(PP, ZPOT)

!> Option 1: First take the gradient of the potential, then build K
!         CALL GRAD(PP%R,ZPOT,KPOT)
!
!         DO I=1,PP%R%NMAX
!            ! option 1: as in /vol/thchem/dewijs/vasp_src/vasp.6/check/vasp6_guppy_clean/build/ncl
!!           KSI(I)=UNITCONV*(1._q/(1._q-ZPOT(I)/(2._q*RYTOEV)/C/C/2._q)**2)*KPOT(I)/PP%R%R(I)
            ! option 2: as in /vol/thchem/dewijs/vasp_src/vasp.6/vasp6_guppy_clean/build/ncl/linear_response_NMR.F-vanbackup
!            KSI(I)=INVMC2*(1._q/(1._q-ZPOT(I)/(2._q*RYTOEV)/C/C/2._q))**2/PP%R%R(I)*KPOT(I)
!!           write (*,*) 'KSI',I,KSI(I)
!         ENDDO
!> Option 2: First build K, then take gradient (of both K and potential).
         KPOT(:) = NMR_ZORA_K(PP%R%NMAX, ZPOT)/2._q !should pick more suitable array name
         CALL GRAD(PP%R,KPOT,KSI)

         LM=1
         DO CH1=1,PP%LMAX
         LMP=1
         DO CH2=1,PP%LMAX

            !< integrate "SIMPI"
            SUM = 0._q
            DO K=1, PP%R%NMAX
               SUM = SUM + PP%R%SI(K) * PP%R%R(K) * KSI(K) * PP%WAE(K,CH1) * PP%WAE(K,CH2)
            ENDDO

            !< BDIR_ONE now contains essentially "B r dK/dr"
            BDIR_ONE(:)=BDIR_ONE(:)*SUM

            LL=PP%LPS(CH1)
            LLP=PP%LPS(CH2)
            !IF (LL == LLP .AND. LL<=LDIMP) THEN
            !IF (LL == LLP .AND. LL>0 .AND. LL<=LMAX ) THEN
            DO M=1,2*LL+1
            DO MP=1,2*LLP+1

!> ----------------------------------------------------------------------------
!>
!> Remember we are calculating:
!> ~~~~~~~~~~~~~~~~~~~~~
!> - \sigma \cdot r ( B \cdot \nabla( (K-1)/4c ) )
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> ``B \cdot \nabla K'' yields terms:
!> ~~~~~~~~~~~~~~~~~~~~~
!> B_i * i/r * dK/dr
!> ~~~~~~~~~~~~~~~~~~~~~
!> ``\sigma \cdot r'' can be written as
!> ~~~~~~~~~~~~~~~~~~~~~
!> \sigma_j * j
!> ~~~~~~~~~~~~~~~~~~~~~
!> in turn we can write terms:
!> ~~~~~~~~~~~~~~~~~~~~~
!> \sigma_j * j/r
!> ~~~~~~~~~~~~~~~~~~~~~
!> if we transform the other side to
!> ~~~~~~~~~~~~~~~~~~~~~
!> B_i * i/r * r * dK/dr
!> ~~~~~~~~~~~~~~~~~~~~~
!> writing out the dot product then gives terms like:
!> ~~~~~~~~~~~~~~~~~~~~~
!> x_i x_j / |x|^2
!> ~~~~~~~~~~~~~~~~~~~~~
!>
!> YLM_XIXJ_YLM contains terms: x_i x_j / |x|^2
!> (array index 1=xx, 2=yy, 3=zz, 4=xy, 5=xz, 6=yz)
!>
!> Ultimately, for this term (contained in array SIG) we need:
!> SIG(1) = < up | \sigma \cdot ... | up >
!>          = (1._q,0._q)* ``z''
!> SIG(2) = < up | \sigma \cdot ... | down >
!>          =  (1._q,0._q)* ``x'' + (0._q,1._q)* ``y''
!> SIG(3) = < down | \sigma \cdot ... | up >
!>          =  (1._q,0._q)* ``x'' - (0._q,1._q)* ``y''
!> SIG(4) = < down | \sigma \cdot ... | down >
!>          = -(1._q,0._q)* ``z''
!>
!> where                (arr index)
!> ``x'' = xx + xy + xz (= 1, 4, 5)
!> ``y'' = yx + yy + yz (= 4, 2, 6)
!> ``z'' = zx + zy + zz (= 5, 6 ,3)
!>
!> The combination of the above leads to the lengthy procedures below.
!> ``r'' of r * dK/dr handled in the SIMPI like integration above.
!>
!> ----------------------------------------------------------------------------
               SIG(1)=(1._q,0._q)*(&
                        BDIR_ONE(1)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,5)&
                        +BDIR_ONE(2)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,6)&
                        +BDIR_ONE(3)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,3)&
                        )

               SIG(2)=(1._q,0._q)*(&
                        BDIR_ONE(1)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,1)&
                        +BDIR_ONE(2)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,4)&
                        +BDIR_ONE(3)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,5)&
                        )&
                        +(0._q,1._q)*(&
                        BDIR_ONE(1)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,4)&
                        +BDIR_ONE(2)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,2)&
                        +BDIR_ONE(3)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,6)&
                        )

               SIG(3)=(1._q,0._q)*(&
                        BDIR_ONE(1)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,1)&
                        +BDIR_ONE(2)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,4)&
                        +BDIR_ONE(3)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,5)&
                        )&
                        -(0._q,1._q)*(&
                        BDIR_ONE(1)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,4)&
                        +BDIR_ONE(2)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,2)&
                        +BDIR_ONE(3)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,6)&
                        )

               SIG(4)=-(1._q,0._q)*(&
                        BDIR_ONE(1)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,5)&
                        +BDIR_ONE(2)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,6)&
                        +BDIR_ONE(3)*YLM_XIXJ_YLM(LL*LL+M,LLP*LLP+MP,3)&
                        )

               !modelled on relativistic.F (loop for sigma)
               DO I=0,1
               DO J=0,1
                  CLQIJ(LMP+MP-1,LM+M-1,NIP,J+2*I+1)=CLQIJ(LMP+MP-1,LM+M-1,NIP,J+2*I+1)-SIG(I+2*J+1)
               ENDDO
               ENDDO

            ENDDO
            ENDDO

            !ENDIF
         LMP=LMP+(2*LLP+1)
         ENDDO
         LM= LM+ (2*LL+1)
         ENDDO

         DEALLOCATE(YLM_XIXJ_YLM)
         DEALLOCATE(ZPOT,KPOT,KSI)
      ENDDO ion

      RETURN
      END SUBROUTINE MLR_CALC_SPINORB_ONE_ALTD

!> ****************** SUBROUTINE FERMICONT_DIPOLAR_PW_S_DENS ******************
!>
!> N.B.: CHTOT must contain the magnetisation density
!>
!> adapted from from hyperfine.F,
!> SUBROUTINE HYPERFINE_PW(T_INFO,LATT_CUR,GRIDC,CHTOT)
!>
!> ****************************************************************************
      SUBROUTINE FERMICONT_DIPOLAR_PW_S_DENS(T_INFO,LATT_CUR,GRIDC,NCDIJ,CHTOT,B_SPIN_PS_)
      USE poscar
      USE lattice
      USE mgrid
      USE constant
      IMPLICIT NONE
      TYPE(type_info) T_INFO
      TYPE(latt) LATT_CUR
      TYPE(grid_3d) GRIDC
      COMPLEX(q) CHTOT(GRIDC%RC%NP,3)               !< only spin density
      REAL(q) B_SPIN_PS_(3,T_INFO%NIONS)
      INTEGER NCDIJ

      !< local variables
      INTEGER NI,I,N1,N2,N3,NC,INDX
      REAL(q) FACTM,GX,GY,GZ,GSQU
      COMPLEX(q) CGR
      DIMENSION INDX(3)

      REAL(q), ALLOCATABLE :: B_SPIN_PW_FC(:,:)     !< Fermi Contact
      REAL(q), ALLOCATABLE :: B_SPIN_PW_DIP(:,:)    !< Dipolar

      ALLOCATE(B_SPIN_PW_FC(3,T_INFO%NIONS))
      ALLOCATE(B_SPIN_PW_DIP(3,T_INFO%NIONS))
      B_SPIN_PW_FC=0._q
      B_SPIN_PW_DIP=0._q

      INDX(1)=1
      INDX(2)=2
      INDX(3)=3

      if (NCDIJ.ne.4) then
         write (*,*) 'error, fermicont_dipolar_pw, stop'
         stop
      endif

      I=0
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW

         I=I+1
         SETFACT1
         SETFACT

         !< x,y,z components and norm of G vector
         GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))*2*PI
         GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))*2*PI
         GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))*2*PI

         GSQU=GX**2+GY**2+GZ**2

         ions: DO NI=1,T_INFO%NIONS
            !< phase factor
            CGR=EXP(CITPI*(GRIDC%LPCTX(N1)*T_INFO%POSION(1,NI)+GRIDC%LPCTY(N2)*T_INFO%POSION(2,NI)+GRIDC%LPCTZ(N3)*T_INFO%POSION(3,NI)))

            !< isotropic (Fermi contact) contribution
            B_SPIN_PW_FC(1,NI)=B_SPIN_PW_FC(1,NI)+MULFACT CHTOT(I,INDX(1))*CGR ! +MULFACT 1._q/(1._q+GSQU)*CGR
            B_SPIN_PW_FC(2,NI)=B_SPIN_PW_FC(2,NI)+MULFACT CHTOT(I,INDX(2))*CGR ! +MULFACT 1._q/(1._q+GSQU)*CGR
            B_SPIN_PW_FC(3,NI)=B_SPIN_PW_FC(3,NI)+MULFACT CHTOT(I,INDX(3))*CGR ! +MULFACT 1._q/(1._q+GSQU)*CGR

            IF (GSQU<1.E-6) THEN
               !< we kill any G=0 is component. here it should vanish... and it should be included in the susceptibility.
            ELSE
               !< anisotropic (dipolar) contribution
               B_SPIN_PW_DIP(1,NI)=B_SPIN_PW_DIP(1,NI)+MULFACT CHTOT(I,INDX(1))*(GX*GX/GSQU-1._q/3._q)*CGR
               B_SPIN_PW_DIP(2,NI)=B_SPIN_PW_DIP(2,NI)+MULFACT CHTOT(I,INDX(2))*(GY*GY/GSQU-1._q/3._q)*CGR
               B_SPIN_PW_DIP(3,NI)=B_SPIN_PW_DIP(3,NI)+MULFACT CHTOT(I,INDX(3))*(GZ*GZ/GSQU-1._q/3._q)*CGR
               B_SPIN_PW_DIP(1,NI)=B_SPIN_PW_DIP(1,NI)+MULFACT CHTOT(I,INDX(2))*(GX*GY/GSQU)*CGR
               B_SPIN_PW_DIP(2,NI)=B_SPIN_PW_DIP(2,NI)+MULFACT CHTOT(I,INDX(1))*(GX*GY/GSQU)*CGR
               B_SPIN_PW_DIP(1,NI)=B_SPIN_PW_DIP(1,NI)+MULFACT CHTOT(I,INDX(3))*(GX*GZ/GSQU)*CGR
               B_SPIN_PW_DIP(3,NI)=B_SPIN_PW_DIP(3,NI)+MULFACT CHTOT(I,INDX(1))*(GX*GZ/GSQU)*CGR
               B_SPIN_PW_DIP(2,NI)=B_SPIN_PW_DIP(2,NI)+MULFACT CHTOT(I,INDX(3))*(GY*GZ/GSQU)*CGR
               B_SPIN_PW_DIP(3,NI)=B_SPIN_PW_DIP(3,NI)+MULFACT CHTOT(I,INDX(2))*(GY*GZ/GSQU)*CGR
            ENDIF
         ENDDO ions

      ENDDO row
      ENDDO col

      !< Bloechl, PRB 62, page 6177, right column
      B_SPIN_PW_FC = B_SPIN_PW_FC /LATT_CUR%OMEGA
      B_SPIN_PW_DIP=-B_SPIN_PW_DIP/LATT_CUR%OMEGA*4._q*PI

      !< the relative weight of the FC and DIP contribution
      B_SPIN_PW_FC = B_SPIN_PW_FC *8._q*PI/3._q
      B_SPIN_PW_DIP= B_SPIN_PW_DIP

      !< communicate
      CALLMPI( M_sum_d(GRIDC%COMM,B_SPIN_PW_FC, 3*T_INFO%NIONS) )
      CALLMPI( M_sum_d(GRIDC%COMM,B_SPIN_PW_DIP,3*T_INFO%NIONS) )

      DO NI=1,T_INFO%NIONS
      DO I=1,3
         B_SPIN_PS_(I,NI)=B_SPIN_PS_(I,NI)+B_SPIN_PW_FC(I,NI)+B_SPIN_PW_DIP(I,NI)
      ENDDO
      ENDDO

      DEALLOCATE(B_SPIN_PW_FC,B_SPIN_PW_DIP)

      RETURN
      END SUBROUTINE FERMICONT_DIPOLAR_PW_S_DENS

!> ******************** SUBROUTINE CALC_B_MAGATOM_LR_FC_DP ********************
!>
!>  based on and stolen from subroutine HYPERFINE_RAD in hyperfine.F
!>
!> ****************************************************************************
      SUBROUTINE CALC_B_MAGATOM_LR_FC_DP(PP,CRHODE,RHOPS,RHOAE,BDIR,NIP,NI,COCC,B_OUT_FC,B_OUT_DP)
!#define s_here
      USE pseudo
      USE constant
      USE asa
      USE radial
      USE wave
      USE paw
      USE poscar
      USE prec
      USE morbitalmag, ONLY : LZORA
      USE hyperfine

      IMPLICIT NONE

      TYPE (potcar), POINTER :: PP   !< pseudopotential descriptor
      COMPLEX(q) :: CRHODE(:,:,:,:)  !< one center occupancy matrix, BEWARE SPIN INDEX
      REAL(q) :: B_OUT_FC(3)
      REAL(q) :: B_OUT_DP(3)
      REAL(q) RHOPS(:,:,:),RHOAE(:,:,:)
      OVERLAP COCC(:,:,:)
      INTEGER :: NIP,NI,BDIR

      !< local variables
      COMPLEX(q) COCC_SPIN
      REAL(q) SCALE,SUM(5),RHOW(PP%R%NMAX)
      REAL(q), ALLOCATABLE :: RHOAE_(:),RHOPS_(:)
      INTEGER I,J,K,LM1,LM2,IDIR,LM
      INTEGER INDX(3)
      REAL(q) TMP(0:4),R(4),DUMMY
      REAL(q) TMP_AE,TMP_PS
      REAL(q) TMP_RAD_ANI_AE(6),TMP_RAD_ANI_PS(6)
      INTEGER, PARAMETER :: ISTRIDE=1
      REAL(q), ALLOCATABLE :: POT(:), KQ(:)

      INDX(1)=2
      INDX(2)=3
      INDX(3)=4

      ! Y00
      SCALE=1/(2._q*SQRT(PI))

      B_OUT_FC=0._q

      ALLOCATE(RHOAE_(PP%R%NMAX),RHOPS_(PP%R%NMAX))

      sdir: DO IDIR=1,3

      !< construct the s-wave spin density only, this is currently not
      !< used in hyperfine.F
      RHOAE_=0._q; RHOPS_=0._q
      !< s_here might yield small improvement, tests undecisive
#ifdef s_here
      LM1=0
      DO I=1,PP%LMAX
         LM1=LM1+(2*PP%LPS(I)+1)
!        IF (PP%LPS(I)/=0.OR.(.NOT.PP%QATO(I,I)/=0._q)) CYCLE
         IF (PP%LPS(I)/=0) CYCLE
         LM2=0
         DO J=1,PP%LMAX
            LM2=LM2+(2*PP%LPS(J)+1)
!           IF (PP%LPS(J)/=0.OR.(.NOT.PP%QATO(J,J)/=0._q)) CYCLE
            IF (PP%LPS(J)/=0) CYCLE
     !      COCC_SPIN=COCC(LM1,LM2,1)-COCC(LM1,LM2,2)
            COCC_SPIN=CRHODE(LM1,LM2,NIP,INDX(IDIR))
     !      COCC_SPIN=COCC(LM1,LM2,INDX(IDIR))
            !write (*,*) 'COCC',LM1,LM2,CRHODE(LM1,LM2,NIP,INDX(IDIR)),COCC(LM1,LM2,INDX(IDIR))
            DO K=1,PP%R%NMAX
               RHOAE_(K)=RHOAE_(K)+PP%WAE(K,I)*PP%WAE(K,J)*COCC_SPIN*SCALE!*PP%R%R(K)**2
               RHOPS_(K)=RHOPS_(K)+PP%WPS(K,I)*PP%WPS(K,J)*COCC_SPIN*SCALE!*PP%R%R(K)**2
            ENDDO
         ENDDO
      ENDDO
#else
      DO K=1,PP%R%NMAX
         RHOAE_(K)=RHOAE(K,1,INDX(IDIR))
         RHOPS_(K)=RHOPS(K,1,INDX(IDIR))
      ENDDO
#endif

      TMP=0._q
      TMP_AE=0._q
      TMP_PS=0._q
      R=0._q
      DUMMY=0._q

!< isotropic (Fermi contact) contribution

!< Below would not do anything, HYPERFINE_FC_ZORA takes priority
!>      DO J=1,4
!>         I=(J-1)*ISTRIDE+1
!>         R(J)=PP%R%R(I)
!>         TMP(J)=RHOAE_(I)/PP%R%R(I)**2
!>      ENDDO
!>      ! extrapolate TMP to R=0
!>      CALL POLINT(R(1:4),TMP(1:4),4,0._q,TMP(0),DUMMY)
!>      TMP_AE=TMP(0)*SCALE
!>

      !< HYPERFINE_FC_ZORA should do ZORA K-factor for Fermi Contact,
      !> should be included in the delta function.
      CALL HYPERFINE_FC_ZORA(PP,RHOAE_,TMP_AE)

      DO J=1,4
         I=(J-1)*ISTRIDE+1
         R(J)=PP%R%R(I)
         TMP(J)=RHOPS_(I)/PP%R%R(I)**2
      ENDDO
      ! extrapolate TMP to R=0
      CALL POLINT(R(1:4),TMP(1:4),4,0._q,TMP(0),DUMMY)
      TMP_PS=TMP(0)*SCALE

      B_OUT_FC(IDIR)=TMP_AE !GIL TEST 9-2-2023    - TMP_PS

      ENDDO sdir

      B_OUT_FC=B_OUT_FC*8._q*PI/3._q

      DEALLOCATE(RHOAE_,RHOPS_)

!< anisotropic (dipolar) contribution
      ALLOCATE(KQ(PP%R%NMAX))
      KQ(:)=1._q
      !> For future development;
      !> consider need for ZORA K-factor in dipolar terms critically
      !> (Most likely needs to be off...)
      IF (LZORA) THEN
         ALLOCATE(POT(PP%R%NMAX))
         CALL NMR_ZORA_POT(PP, POT)
         KQ(:)=NMR_ZORA_K(PP%R%NMAX, POT)
      ENDIF

      B_OUT_DP=0._q

      IF (SIZE(RHOAE,2)<5) RETURN !< no d channel

      SCALE=SQRT(12._q*PI/5._q)

      sdir2: DO IDIR=1,3

      DO LM=5,9
         DO I=1,PP%R%NMAX
            RHOW(I)=KQ(I)*RHOAE(I,LM,INDX(IDIR))/PP%R%R(I)**3
         ENDDO
         CALL SIMPI(PP%R,RHOW,SUM(LM-4))
      ENDDO

      TMP_RAD_ANI_AE(1)=SCALE*( SUM(5)-SUM(3)/SQRT(3._q)) !<xx
      TMP_RAD_ANI_AE(2)=SCALE*(-SUM(5)-SUM(3)/SQRT(3._q)) !<yy
      TMP_RAD_ANI_AE(3)=SCALE*( SUM(3)*SQRT(4._q/3._q))   !<zz
      TMP_RAD_ANI_AE(4)=SCALE*( SUM(1))                   !<xy
      TMP_RAD_ANI_AE(5)=SCALE*( SUM(4))                   !<xz
      TMP_RAD_ANI_AE(6)=SCALE*( SUM(2))                   !<yz

      DO LM=5,9
         DO I=1,PP%R%NMAX
            RHOW(I)=RHOPS(I,LM,INDX(IDIR))/PP%R%R(I)**3
         ENDDO
         CALL SIMPI(PP%R,RHOW,SUM(LM-4))
      ENDDO

      TMP_RAD_ANI_PS(1)=SCALE*( SUM(5)-SUM(3)/SQRT(3._q))
      TMP_RAD_ANI_PS(2)=SCALE*(-SUM(5)-SUM(3)/SQRT(3._q))
      TMP_RAD_ANI_PS(3)=SCALE*( SUM(3)*SQRT(4._q/3._q))
      TMP_RAD_ANI_PS(4)=SCALE*( SUM(1))
      TMP_RAD_ANI_PS(5)=SCALE*( SUM(4))
      TMP_RAD_ANI_PS(6)=SCALE*( SUM(2))

      !< GIL let's hope this is correct
      IF (IDIR.EQ.1) THEN
         B_OUT_DP(1) = B_OUT_DP(1) + TMP_RAD_ANI_AE(1) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(1)
         B_OUT_DP(2) = B_OUT_DP(2) + TMP_RAD_ANI_AE(4) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(4)
         B_OUT_DP(3) = B_OUT_DP(3) + TMP_RAD_ANI_AE(5) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(5)
      ENDIF

      IF (IDIR.EQ.2) THEN
         B_OUT_DP(1) = B_OUT_DP(1) + TMP_RAD_ANI_AE(4) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(4)
         B_OUT_DP(2) = B_OUT_DP(2) + TMP_RAD_ANI_AE(2) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(2)
         B_OUT_DP(3) = B_OUT_DP(3) + TMP_RAD_ANI_AE(6) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(6)
      ENDIF

      IF (IDIR.EQ.3) THEN
         B_OUT_DP(1) = B_OUT_DP(1) + TMP_RAD_ANI_AE(5) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(5)
         B_OUT_DP(2) = B_OUT_DP(2) + TMP_RAD_ANI_AE(6) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(6)
         B_OUT_DP(3) = B_OUT_DP(3) + TMP_RAD_ANI_AE(3) !GIL TEST 9-2-2023    - TMP_RAD_ANI_PS(3)
      ENDIF

      ENDDO sdir2

      DEALLOCATE(KQ)
      IF (LZORA) DEALLOCATE(POT)

      END SUBROUTINE CALC_B_MAGATOM_LR_FC_DP

!> ************************* SUBROUTINE NMR_ZORA_POT **************************
!>
!> Adapt potential to include ZORA effects
!> modelled from relativistic.F SPINORB_STRENGTH
!>
!> ****************************************************************************
      SUBROUTINE NMR_ZORA_POT(PP, POT)
      USE prec
      USE constant
      USE pseudo
      USE radial
      IMPLICIT NONE

      TYPE(potcar), POINTER :: PP
      REAL(q) :: POT(:), RDUM, SCALE

      SCALE=1/(2*SQRT(PI)) !< scale factor as in relativistic.F

      POT=0._q

      !< compute V_H[n^a_c]
      CALL RAD_POT_HAR(0,PP%R,POT,PP%RHOAE,RDUM)

      !< add -Z/r
      POT=POT*SCALE-(PP%ZCORE+PP%ZVALF_ORIG)*FELECT/PP%R%R

      !< add V_H[n^a_v] + V_xc[n^a_v+n^a_c] (contained in PP%POTAE)
      POT=POT-PP%POTAE

      !> and make sure this AE-"atomic reference"-potential lines up exactly with
      !> the potential that was constructed by the pseudopotential generation code.
      POT=POT+(PP%VPSRMAX-POT(PP%R%NMAX))
!      WRITE(1001,'(F14.7,4X,2F22.10)') (PP%R%R(M),POT(M),-PP%POTPS(M)+PP%POTPSC(M),M=1,SIZE(PP%POTAE))

!TESTING...
!      !WRITE GRID
!      WRITE(99, *) "NMR_ZORA_POT GRID"
!      WRITE(99, '(3I5)') PP%R%NMAX
!      WRITE(99, '(5(1x,E17.11))') (PP%R%R(I), I=1, PP%R%NMAX)
!      !WRITE RHOAE
!      WRITE(99, *) "NMR_ZORA_POT PP%RHOAE"
!      WRITE(99, '(3I5)') PP%R%NMAX
!      WRITE(99, '(5(1x,E17.11))') (PP%RHOAE(I), I=1, PP%R%NMAX)
!      !WRITE POTAE
!      WRITE(99, *) "NMR_ZORA_POT PP%POTAE"
!      WRITE(99, '(3I5)') PP%R%NMAX
!      WRITE(99, '(5(1x,E17.11))') (PP%POTAE(I), I=1, PP%R%NMAX)
!      !WRITE POT
!      WRITE(99, *) "NMR_ZORA_POT POT"
!      WRITE(99, '(3I5)') PP%R%NMAX
!      WRITE(99, '(5(1x,E17.11))') (POT(I), I=1, PP%R%NMAX)
!
!      CALL SLEEP(5)
!      STOP

      RETURN
      END SUBROUTINE NMR_ZORA_POT

!> *************************** FUNCTION NMR_ZORA_K ****************************
!>
!> ZORA: radially variant relativistic mass correction 1/(1-V(r)/2/c^2) in Hartree atomic units,
!> in accordance with Eq. A13 of Yates et al., J. Chem. Phys. 118, 5746 (2003).
!>
!> For VASP units, consider relativistic.F and Steiner et al. (Phys. Rev. B. 93, 224425 (2016))
!> K = (1-V/(2 m_e c^2))^-1 (Note, Steiner uses ^-2 which is not standard.)
!> Unit conversion of "me_2 c^2" to VASP leads to (2*RYTOEV*CLIGHT*CLIGHT)^-1.
!>
!> ****************************************************************************
      FUNCTION NMR_ZORA_K(NMAX,POT)
      USE prec
      USE constant
      IMPLICIT NONE

      INTEGER :: NMAX
      REAL(q) :: POT(:)
      REAL(q), dimension(NMAX) :: NMR_ZORA_K

      NMR_ZORA_K(:)=0._q
      NMR_ZORA_K(:)=1._q/(1._q-POT(:)/(2._q*RYTOEV)/CLIGHT/CLIGHT/2._q)

      END FUNCTION NMR_ZORA_K

! ************************ SUBROUTINE NMR_MAGRES_TENS *************************
!
!> Subroutine to write to magres file.
!
! ****************************************************************************
      SUBROUTINE NMR_MAGRES_TENS(P,T_INFO, IO)
         use prec
         use pseudo_struct_def
         use poscar
         use morbitalmag
         
         IMPLICIT NONE
   
         TYPE (potcar), TARGET :: P(:)
         TYPE (latt) LATT_CURR
         TYPE (type_info) T_INFO
         INTEGER IO, NI, I, J

         ! ABS_CHEM_SHIFT contains core + valence shielding and susceptibility,
         ! one could argue to exclude susceptibility and only write it separate.
         ! writing of comment not supported inside [magres][/magres] block.
         ! WRITE(IO, '(A74)') "#VASP ms includes susceptibility (G=0 component assuming spherical sample)"
         WRITE(IO, '(A12)') "units ms ppm"
         NI = 0
         DO I=1, T_INFO%NTYP
            DO J=1, T_INFO%NITYP(I)
               NI = NI + 1
               !ms atomtype #num xx xy xz yx yy yz zx zy zz
               WRITE(IO, '(A3,A2,I4,9(F12.4))') "ms ", P(T_INFO%ITYP(NI))%ELEMENT, J, REAL(ABS_CHEM_SHIFT(1:3,1:3,NI),q)
            ENDDO
         ENDDO

         WRITE(IO, '(A27)') "units sus 10^-6.cm^3.mol^-1"
         WRITE(IO, '(A4,9(F12.4))') "sus ", CHI_MAG(1:3,1:3)

      RETURN
      END SUBROUTINE NMR_MAGRES_TENS

      END MODULE mlr_main_nmr

!> ******************** SUBROUTINE PSCURRENT_PARA_TWO_KPT *********************
!>
!> gradient psi contribution to current: Im[psi^*(r) nabla psi(r)]
!> for a supplied weighted k-point and psi
!> with psi(r) = e(iGr)
!> we get contributions proportional to (iG) e(iGr)
!>
!> ****************************************************************************
      SUBROUTINE PSCURRENT_PARA_TWO_KPT(WDES1, GRID, CR_1, CR_2, IGX, IGY, IGZ, &
     &       VKPT_1, VKPT_2, WEIGHT, CWORK1, CWORK2, CW_1, CW_2, J, S_DENS, LSOSHIFT)

!> GIL 04-01-2023 In the guppy version of 2019 there are several more
!> GIL 04-01-2023 elaborate versions of the LSOSHIFT addition to this subroutine.
!> GIL 04-01-2023 It could well be possible that even the present !LSOSHIFT extension can
!> GIL 04-01-2023 be removed, because it might well be that the current density gives
!> GIL 04-01-2023 negligible contribution to the shielding. This needs to be checked.

      USE prec
      USE wave
      USE mgrid
      USE constant
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      TYPE (grid_3d)     GRID
      COMPLEX(q) :: CR_1(GRID%MPLWV*WDES1%NRSPINORS)  !< wavefunction in real space
      COMPLEX(q) :: CR_2(GRID%MPLWV*WDES1%NRSPINORS)  !< wavefunction in real space
      COMPLEX(q) :: CW_1(WDES1%NRPLWV)
      COMPLEX(q) :: CW_2(WDES1%NRPLWV)
      REAL(q)    :: VKPT_1(3)                         !< k-point in reciprocal lattice
      REAL(q)    :: VKPT_2(3)                         !< k-point in reciprocal lattice
      REAL(q)    :: WEIGHT            !< same for both k-points
      INTEGER    :: IGX(WDES1%NGDIM)  !< component of each G vector in direction of first rec lattice vector
      INTEGER    :: IGY(WDES1%NGDIM)  !< component of each G vector in direction of second rec lattice vector
      INTEGER    :: IGZ(WDES1%NGDIM)  !< component of each G vector in direction of third rec lattice vector
      !< we assume both states consist of the same set of G-vectors
      COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
      COMPLEX(q) :: CWORK2(GRID%MPLWV)
      RGRID      :: J(GRID%RL%NP,3)      !< three component vector storing current density in real space
      RGRID      :: S_DENS(GRID%RL%NP,3) !< three component vector storing spin density in real space GIL SOC
      LOGICAL LSOSHIFT

      !< local
      INTEGER    :: ISPINOR, M, MM
      COMPLEX(q) :: IC

      DO ISPINOR=0,WDES1%NRSPINORS-1
!> ----------------------------------------------------------------------------
!>  First do paramagnetic operator (which is jbare)
!> ----------------------------------------------------------------------------
      !< first we take the gradient of wave function 2, fft that to real space,
      !> and multiply with 1

         !< component along first axis
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CWORK1(M)=CW_2(MM)*(IGX(M)+VKPT_2(1))*TPI
         ENDDO
         !< FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            J(M,1)=J(M,1)+(CWORK2(M)*CONJG(CR_1(MM)))*WEIGHT / 2._q
         ENDDO

         !< component along second axis
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CWORK1(M)=CW_2(MM)*(IGY(M)+VKPT_2(2))*TPI
         ENDDO
         !< FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            J(M,2)=J(M,2)+(CWORK2(M)*CONJG(CR_1(MM)))*WEIGHT / 2._q
         ENDDO

         !< component along third axis
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CWORK1(M)=CW_2(MM)*(IGZ(M)+VKPT_2(3))*TPI
         ENDDO
         !< FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            J(M,3)=J(M,3)+(CWORK2(M)*CONJG(CR_1(MM)))*WEIGHT / 2._q
         ENDDO

      !> next we take the gradient of wave function 1, fft that to real space,
      !> and multiply with 2
      !> we keep the conjugation on the wave function 1

         !< component along first axis
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CWORK1(M)=CW_1(MM)*(IGX(M)+VKPT_1(1))*TPI
         ENDDO
         !< FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            J(M,1)=J(M,1)+(CONJG(CWORK2(M))*CR_2(MM))*WEIGHT / 2._q
         ENDDO

         !< component along second axis
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CWORK1(M)=CW_1(MM)*(IGY(M)+VKPT_1(2))*TPI
         ENDDO
         !< FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            J(M,2)=J(M,2)+(CONJG(CWORK2(M))*CR_2(MM))*WEIGHT / 2._q
         ENDDO

         !< component along third axis
         DO M=1,WDES1%NGVECTOR
            MM=M+ISPINOR*WDES1%NGVECTOR
            CWORK1(M)=CW_1(MM)*(IGZ(M)+VKPT_1(3))*TPI
         ENDDO
         !< FFT to real space
         CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *WDES1%GRID%MPLWV
            J(M,3)=J(M,3)+(CONJG(CWORK2(M))*CR_2(MM))*WEIGHT / 2._q
         ENDDO

      ENDDO

!> ----------------------------------------------------------------------------
!> If LSOSHIFT, do spin operator/pauli matrix (i.e. mbare)
!> ----------------------------------------------------------------------------

      IF (LSOSHIFT) THEN
         IC=CMPLX(0._q,1._q, q) !GIL this looks strange
         DO M=1,GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            !< the three rules have the pauli spin matrices x,y,z
            S_DENS(M,1)=S_DENS(M,1)+    ( CONJG(CR_1( M))*CR_2(MM)+CONJG(CR_1(MM))*CR_2( M) )*WEIGHT
            S_DENS(M,2)=S_DENS(M,2)+ IC*(-CONJG(CR_1( M))*CR_2(MM)+CONJG(CR_1(MM))*CR_2( M) )*WEIGHT
            S_DENS(M,3)=S_DENS(M,3)+    ( CONJG(CR_1( M))*CR_2( M)-CONJG(CR_1(MM))*CR_2(MM) )*WEIGHT
         ENDDO
      ENDIF

      END SUBROUTINE PSCURRENT_PARA_TWO_KPT

!> ************************ SUBROUTINE MLR_SYM_REDUCE *************************
!>
!> ****************************************************************************
!#define debug
      SUBROUTINE MLR_SYM_REDUCE(LATT_CUR,T_INFO,SYMM,IU6)
      USE prec
      USE base
      USE lattice
      USE poscar
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(type_info) T_INFO
      TYPE(symmetry) SYMM
      INTEGER IU6
      !< local variables
      INTEGER I,J,K,IS
      INTEGER INOT,INROT,ISORNOT,LIST(48)
      REAL(q) QST(3,3),RQST(3,3)
      REAL(q) V1,V2,V3
      REAL(q) MODQST,MODRQST,PRODMOD,DOT
      INTEGER ISYMOPTMP(3,3,48),IGRPOPTMP(3,3,48)
      REAL(q) GTRANSTMP(3,48)
      LOGICAL ISMAPPED(3)
      INTEGER IPGIND
      CHARACTER(4) GRPNAM
      REAL(q) TAU(T_INFO%NIONS,3),TAUROT(T_INFO%NIONS,3)
      !< symmetry common block
      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
           GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

#ifdef debug
      IF (IU6>0) THEN
      WRITE(*,'(A)') 'OLD ISYMOP'
      DO IS=1,NROTK
         WRITE(*,'(A,I5)') 'SYM. OP',IS
         DO I=1,3
            WRITE(*,'(7I5)') (ISYMOP(I,J,IS),J=1,3), (IGRPOP(I,J,IS),J=1,3),INVMAP(IS)
         ENDDO
         WRITE(*,*) 'OLD GTRANS'
         WRITE(*,'(3F8.5)') (GTRANS(J,IS),J=1,3)
      ENDDO
      WRITE(*,*) 'NROT=',NROT,'NROTK=',NROTK
      ENDIF
#endif

      INOT=0; INROT=0; LIST=1

      !< loop on sym. op.
      isym: DO IS=1,NROTK
!        WRITE (*,'(A,I5)') 'SYMOP=',IS
         !< q-star
         QST=0._q; QST(1,1)=1._q; QST(2,2)=1._q; QST(3,3)=1._q

         !< cartesian to direct
         DO I=1,3 !< 3 vectors
            V1=QST(1,I)*LATT_CUR%B(1,1)+QST(2,I)*LATT_CUR%B(2,1)+QST(3,I)*LATT_CUR%B(3,1)
            V2=QST(1,I)*LATT_CUR%B(1,2)+QST(2,I)*LATT_CUR%B(2,2)+QST(3,I)*LATT_CUR%B(3,2)
            V3=QST(1,I)*LATT_CUR%B(1,3)+QST(2,I)*LATT_CUR%B(2,3)+QST(3,I)*LATT_CUR%B(3,3)
            QST(1,I)=V1
            QST(2,I)=V2
            QST(3,I)=V3
         ENDDO

         !< rotate the q-star
         RQST=0._q
         DO K=1,3 !< x y z direction of IQST
            DO I=1,3 !< component of IRQST
               DO J=1,3 !< summation index (row column)
                  RQST(I,K)=RQST(I,K)+IGRPOP(I,J,IS)*QST(J,K)
               ENDDO
            ENDDO
         ENDDO

         ISORNOT = 0

         !< test colinearity abs ( u dot v ) == u*v
         idir: DO K=1,3 !< x y z direction for QST
            ISMAPPED=.FALSE.
            MODQST=QST(1,K)*QST(1,K)+QST(2,K)*QST(2,K)+QST(3,K)*QST(3,K)
            MODQST=SQRT(MODQST)
            DO I=1,3 !< x y z direction for RQST
               DOT = 0.0_q
               MODRQST = RQST(1,I)*RQST(1,I) + RQST(2,I)*RQST(2,I) + RQST(3,I)*RQST(3,I)
               MODRQST = SQRT(MODRQST)
               PRODMOD = MODQST * MODRQST
               ! dot product
               DO J=1,3
                  DOT=DOT+RQST(J,I)*QST(J,K)
               ENDDO
               IF (ABS(ABS(DOT)-PRODMOD).LT.1E-8) THEN
                  ISMAPPED(I)=.TRUE.
               ENDIF
#ifdef debug
              WRITE(*,'(A,3F9.5,A,L,3F9.5)') 'rotated ',RQST(:,I),' is it mapped ? ',ISMAPPED(I),DOT,PRODMOD,ABS(DABS(DOT)-PRODMOD)
#endif
            ENDDO

            IF ((.NOT.ISMAPPED(1)).AND.(.NOT.ISMAPPED(2)).AND.(.NOT.ISMAPPED(3)))  THEN
               ISORNOT=ISORNOT+1
            ENDIF
#ifdef debug
           WRITE(*,*) IS , K,ISORNOT,ISMAPPED(1),ISMAPPED(2),ISMAPPED(3)
#endif
         ENDDO idir


         IF (ISORNOT.NE.0) THEN
#ifdef debug
            IF (IU6>0) WRITE(*,'(A,I3,A)') 'SYM. OP.',IS,' is NOT SUITABLE FOR LR_NMR and set to zero'
#endif
            INOT=INOT+1
            ISYMOP(:,:,IS)=0
            IGRPOP(:,:,IS)=0
            GTRANS(:,IS)=0
            LIST(IS)=0
         ENDIF
      ENDDO isym

#ifdef debug
      IF (IU6>0) WRITE(*,*) 'NUMBER OF SYMMETRY OPERATIONS NOT SUITABLE',INOT
#endif

      !< reorganize the arrays
      K=0
      DO IS=1,NROTK
         IF (LIST(IS).EQ.1) THEN
            K=K+1
            ISYMOPTMP(:,:,K)=ISYMOP(:,:,IS)
            IGRPOPTMP(:,:,K)=IGRPOP(:,:,IS)
            GTRANSTMP(:,K)   =GTRANS(:,IS)
            IF ((GTRANS(1,IS).EQ.0).AND.(GTRANS(2,IS).EQ.0).AND.(GTRANS(3,IS).EQ.0)) INROT=INROT+1
         ENDIF
      ENDDO

      ISYMOP(:,:,:)=ISYMOPTMP(:,:,:)
      IGRPOP(:,:,:)=IGRPOPTMP(:,:,:)
      GTRANS(:,:)  =GTRANSTMP(:,:)
      NROTK=NROTK-INOT
      NROT=INROT

#ifdef debug
      IF (IU6>0) THEN
      WRITE(*,'(A)') 'NEW ISYMOP'
      DO IS=1,NROTK
         WRITE(*,'(A,I5)') 'SYM. OP',IS
         DO I=1,3
            WRITE(*,'(7I5)') (ISYMOP(I,J,IS),J=1,3),(IGRPOP(I,J,IS),J=1,3),INVMAP(IS)
         ENDDO
         WRITE(*,*) 'NEW GTRANS'
         WRITE(*,'(3F8.5)') (GTRANS(J,IS),J=1,3)
      ENDDO
      ENDIF
#endif

!> ----------------------------------------------------------------------------
!> Now try to find out which symmetry operations we have here:
!> ----------------------------------------------------------------------------
      CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
      IF (IU6>0) &
      WRITE(IU6,'(/4A)') '(MLR_SYM_REDUCE) The static configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
      IF (NROTK>NROT) THEN
         CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') '(MLR_SYM_REDUCE) The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
      END IF
      DO I=1,T_INFO%NIONS
         TAU(I,1)=T_INFO%POSION(1,I)
         TAU(I,2)=T_INFO%POSION(2,I)
         TAU(I,3)=T_INFO%POSION(3,I)
      ENDDO
      IF (ASSOCIATED(SYMM%ROTMAP)) THEN
         DEALLOCATE(SYMM%ROTMAP)
      ENDIF
      ALLOCATE(SYMM%ROTMAP(T_INFO%NIONS,NROTK,NPCELL))
      CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,SYMM%PTRANS,NPCELL,T_INFO%NIONS,1, &
     &            T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,SYMM%ROTMAP(1,1,1),TAUROT)

!> ----------------------------------------------------------------------------
!> map of the inverse elements of each (reciprocal space) group element
!> ----------------------------------------------------------------------------
      CALL INVGRP(IGRPOP,NROTK,INVMAP)

      SYMM%NROT=NROTK

      RETURN
      END SUBROUTINE MLR_SYM_REDUCE

!> ************************ SUBROUTINE NMR_TENS_WRITE *************************
!>
!> Simple routine to write NMR tensors in appropriate format
!>
!> ****************************************************************************

      SUBROUTINE NMR_WRITE_TENS(NIONS, NIOND, TENS, IO)
      use base

      INTEGER NIONS, NIOND
      REAL(q) TENS(3,3,NIOND)
      TYPE (in_struct) IO

      WRITE (IO%IU6,*) '----------------------------------------------------------------'
      WRITE (IO%IU6,*) ' ion  BDIR              X              Y              Z'
      WRITE (IO%IU6,*) '----------------------------------------------------------------'

      DO NI=1, NIONS
         IF (NI.NE.1) THEN
            WRITE (IO%IU6,*)
         ENDIF
                                          !ion bdir  x y z
         WRITE (IO%IU6,'(I5,2X,I4,3F15.6)') NI,   1, TENS(:,1,NI)
         WRITE (IO%IU6,'(5X,2X,I4,3F15.6)')       2, TENS(:,2,NI)
         WRITE (IO%IU6,'(5X,2X,I4,3F15.6)')       3, TENS(:,3,NI)
      ENDDO

      WRITE (IO%IU6,*) '----------------------------------------------------------------'

      RETURN
      END SUBROUTINE
