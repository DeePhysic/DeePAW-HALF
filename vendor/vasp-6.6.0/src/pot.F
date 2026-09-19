#include "symbol.inc"
      MODULE pot
      USE prec
      USE base, ONLY: write_pot
      USE charge
      USE lattice, ONLY: latt
      USE pot_struct_def
      USE pot_electrostat
      USE coulomb_cutoff

      COMPLEX(q), ALLOCATABLE, PRIVATE :: EXTERNAL_POTENTIAL(:)
      PRIVATE ADD_EXTERNAL_POTENTIAL

      INTERFACE str
         MODULE PROCEDURE WRT_POTENTIAL_TO_STRING
      END INTERFACE str

      CONTAINS
!************************ SUBROUTINE POTLOK ****************************
!
!> This subroutine calculates the total local potential CVTOT
!> which is the sum of the hartree potential, the exchange-correlation
!> potential and the ionic local potential. The routine also calculates
!> the total local potential SV on the small grid.
!>
!> On entry: \n
!>  CHTOT(:,1)    density \n
!>  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization \n
!> On return (LNONCOLLINEAR=.FALSE.): \n
!>  CVTOT(:,1)    potential for up \n
!>  CVTOT(:,2)    potential for down \n
!>
!> On return (LNONCOLLINEAR=.TRUE.): \n
!>  CVTOT(:,1:4)  spinor representation of potential
!>
!> @details @ref openacc :
!> On entering this subroutine OpenACC execution is switched on.
!> On exit the OpenACC execution mode reverts to its original status,
!> and SV, CHTOT, CVTOT are copied back to the host.
!>
!
! RCS:  $Id: pot.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!***********************************************************************
    SUBROUTINE POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,COMM_INTER,WDES,  &
               INFO,XC,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF, &
               COULOMB_POT)

__GPU USE moffload

      USE prec
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY: fexcf, fexcp
      USE tau_mu, ONLY : tau_handle
      USE hamil_struct_def, ONLY : ham_handle
      USE base
      USE wave
      USE mdipol
      USE Constrained_M_modular
      USE main_mpi, ONLY: COMM
      USE scpc, ONLY: SCPC_APPLY
      USE plugins, ONLY: PLUGINS_LOCAL_POTENTIAL
! solvation__
      USE solvation
! solvation__
! bexternal__
      USE bexternal
! bexternal__
! embedding__
      USE mextpot
! embedding__
      IMPLICIT NONE

      TYPE (ham_handle)  HAMILTONIAN
      TYPE (tau_handle)  KINEDEN
      TYPE (xc_info)     XC
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (wavedes)     WDES
      TYPE (transit)     SOFT_TO_C
      TYPE (info_struct) INFO
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (energy)      E
      TYPE (latt)        LATT_CUR
      TYPE (communic)    COMM_INTER
      TYPE (coulomb_potential), OPTIONAL :: COULOMB_POT

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP), CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      RGRID, INTENT(OUT) :: SV(DIMREAL(GRID%MPLWV), WDES%NCDIJ)
      COMPLEX(q), INTENT(OUT) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      COMPLEX(q), INTENT(IN), TARGET :: CHTOT(GRIDC%MPLWV, WDES%NCDIJ)
      RGRID      DENCOR(GRIDC%RL%NP)
      REAL(q)    XCSIF(3,3),TMPSIF(3,3),XDMTMP(1)
! work arrays (allocated after call to FEXCG)
      COMPLEX(q), ALLOCATABLE::  CWORK1(:),CWORK(:,:)
      REAL(q) ELECTROSTATIC,EXC,EXCG,RINPL,XCENCG,FACTM
      LOGICAL, EXTERNAL :: L_NO_LSDA_GLOBAL
      INTEGER I,ISP,MWORK1,N1,N2,N3,NC,NG, ASYNC_PRE
      COMPLEX(q) CVZERG

      RGRID , POINTER :: MU_(:,:)
      COMPLEX(q), POINTER :: MUTOT_(:,:)
      TYPE(COULOMB_POTENTIAL):: COULOMB_POT_PERIODIC
      COMPLEX(q),ALLOCATABLE :: CVCOL(:,:)
#ifdef PLUGINS
! The local potential interface sends in these two potentials to Python
! As well as the total charge density in real space
      COMPLEX(q), ALLOCATABLE :: HARTREE_POTENTIAL(:)
      COMPLEX(q), ALLOCATABLE :: ION_POTENTIAL(:)
      COMPLEX(q), ALLOCATABLE :: TOTAL_CHARGE_DENSITY(:,:)
      COMPLEX(q), ALLOCATABLE :: PLUGINS_POTENTIAL(:)
      REAL(q), ALLOCATABLE    :: DIPOLE_MOMENT(:)
#endif

#ifdef libbeef
      LOGICAL LBEEFCALCBASIS, LBEEFBAS
      COMMON /BEEFENS/ LBEEFCALCBASIS,LBEEFBAS

      REAL(q), ALLOCATABLE, SAVE :: BEEFXC(:), BEENERGIES(:)
      LOGICAL SAVELUSE_VDW
      INTEGER BEEFCOUNTER
#endif

      COMPLEX(q), POINTER :: CTMP(:,:)

      PROFILING_START('potlok')

      HAMILTONIAN%XC = XC

#ifdef OFFLOADING
      CALL OFFLOAD_ON_PUSH(.TRUE.)

      ! set ACC_ASYNC_Q to ACC_ASYNC_SYNC to avoid race conditions
      ! not all code below is checked, and store the original value
GPACC ASYNC_PRE=ACC_ASYNC_Q
GPACC ACC_ASYNC_Q=ACC_ASYNC_SYNC

      CALL COPYIN_TYPED_VAR(GRIDC)
      CALL COPYIN_TYPED_VAR(GRID_SOFT)
      CALL COPYIN_TYPED_VAR(SOFT_TO_C)
      CALL COPYIN_TYPED_VAR(P)
      CALL COPYIN_TYPED_VAR(XC)
      IF (.NOT.PRESENT(COULOMB_POT)) CALL COPYIN_TYPED_VAR(COULOMB_POT_PERIODIC)
GPACC !$ACC ENTER DATA COPYIN(CHTOT,LATT_CUR,DENCOR,CSTRF) CREATE(CVTOT,SV) __IF_ASYNC__
GPACC !$ACC ENTER DATA COPYIN(KINEDEN,KINEDEN%TAUC) IF(ASSOCIATED(KINEDEN%TAUC).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC ENTER DATA COPYIN(KINEDEN,KINEDEN%TAU) IF(ASSOCIATED(KINEDEN%TAU).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET ENTER DATA MAP(TO:CHTOT,LATT_CUR,DENCOR,CSTRF) MAP(ALLOC:CVTOT,SV) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:KINEDEN,KINEDEN%TAUC) IF(ASSOCIATED(KINEDEN%TAUC).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:KINEDEN,KINEDEN%TAU) IF(ASSOCIATED(KINEDEN%TAU).AND.OFFLOAD_ON)
#endif

      IF (XC%LNO_AUG_XC) THEN
         ALLOCATE(CTMP(GRIDC%MPLWV,WDES%NCDIJ))
GPACC !$ACC ENTER DATA COPYIN(CHDEN) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:CHDEN) IF(OFFLOAD_ON)
      ELSE
         CTMP => CHTOT
      ENDIF
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: CTMP) IF(OFFLOAD_ON)

      IF (XC%LMU) THEN
         IF (.NOT.ASSOCIATED(HAMILTONIAN%MU) .OR. .NOT.ASSOCIATED(HAMILTONIAN%MUTOT)) &
            CALL vtutor%bug("POTLOK: MU and/or MUTOT not associated.", __FILE__, __LINE__)
         MU_ => HAMILTONIAN%MU ; MUTOT_ => HAMILTONIAN%MUTOT
GPACC !$ACC ENTER DATA CREATE(MU_,MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:MU_,MUTOT_) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: MU_,MUTOT_) IF(OFFLOAD_ON)
! here the full size of MU and MUTOT has to be given, the cce 19
! compiler thinks the size is 1,1
#ifdef __DCRAYFTN_VERSION
GPOMP !$OMP TARGET MAP(PRESENT,TO:MU_,MUTOT_,HAMILTONIAN%MU,HAMILTONIAN%MUTOT) IF(OFFLOAD_ON)
         MU_ => HAMILTONIAN%MU(1:GRID_SOFT%MPLWV,1:WDES%NCDIJ)
         MUTOT_ => HAMILTONIAN%MUTOT(1:GRIDC%MPLWV,1:WDES%NCDIJ)
GPOMP !$OMP END TARGET
#endif

GPACC !$ACC KERNELS PRESENT(MU_,MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         MU_=0 ; MUTOT_=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ELSE
         ALLOCATE(MU_(1,WDES%NCDIJ),MUTOT_(1,WDES%NCDIJ))
GPACC !$ACC ENTER DATA CREATE(MU_,MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:MU_,MUTOT_) IF(OFFLOAD_ON)
      ENDIF

      MWORK1=MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV)
      ALLOCATE(CWORK1(MWORK1),CWORK(GRIDC%MPLWV,WDES%NCDIJ))
GPACC !$ACC ENTER DATA CREATE(CWORK,CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK,CWORK1) IF(OFFLOAD_ON)

!-----------------------------------------------------------------------
!
!  calculate the exchange correlation potential and the dc. correction
!
!-----------------------------------------------------------------------
      EXC     =0
      E%XCENC =0
      E%EXCG  =0
      E%CVZERO=0
      XCSIF   =0

GPACC !$ACC KERNELS PRESENT(CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVTOT) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CVTOT   =0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

#ifdef libbeef
      IF(LBEEFCALCBASIS) THEN
         BEEFCOUNTER = 1
         SAVELUSE_VDW = XC%LUSE_VDW
!        for calculation of nscf xc energies for error estimates,
!        vdW interactions don't need to be calculated
         XC%LUSE_VDW = .FALSE.
         IF(.NOT.ALLOCATED(BEEFXC)) ALLOCATE(BEEFXC(32))
         IF(.NOT.ALLOCATED(BEENERGIES)) ALLOCATE(BEENERGIES(2000))
      ENDIF
#endif

      IF (.NOT.XC%LNOXC) THEN
         ! transform the charge density to real space
         EXCG  =0
         XCENCG=0
         CVZERG=0
         TMPSIF=0

         IF (XC%LNO_AUG_XC) THEN
            ! We will use the soft charge density instead of the augmented one
GPACC !$ACC KERNELS PRESENT(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CTMP) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            CTMP=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
            DO ISP=1,WDES%NCDIJ
               CALL ADD_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CHDEN(1,ISP),CTMP(1,ISP))
               CALL SETUNB_COMPAT(CTMP(1,ISP),GRIDC)
               CALL FFT3D(CTMP(1,ISP),GRIDC,1)
            ENDDO
         ENDIF
         ! Bring the augmented charge density to real space
         DO ISP=1,WDES%NCDIJ
            CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
         ENDDO

         ! Bring the kinetic energy density to real space if needed for meta-GGA functionals
         IF (XC%LTAU) THEN
            DO ISP=1,WDES%NCDIJ
               CALL FFT3D(KINEDEN%TAU(1,ISP),GRIDC,1)
            ENDDO
#ifndef realmode
GPACC !$ACC KERNELS PRESENT(KINEDEN%TAU) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:KINEDEN%TAU) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            KINEDEN%TAU=REAL(KINEDEN%TAU,KIND=q)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#endif
         ENDIF

#ifdef libbeef
 4096    IF (LBEEFCALCBASIS) THEN
            !for calculation of nscf xc energies for error estimates,
            !cycle through beef legendre polynomial basis
            IF(BEEFCOUNTER.LT.31) THEN
               CALL BEEFSETMODE(BEEFCOUNTER-1)
            ELSE
               CALL BEEFSETMODE(BEEFCOUNTER-34)
            ENDIF
         ENDIF
#endif

         IF (WDES%ISPIN==2) THEN

            ! get the charge and the total magnetization
            CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)
! do LDA+U instead of LSDA+U
            IF (L_NO_LSDA_GLOBAL()) THEN
GPACC !$ACC KERNELS PRESENT(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
               CWORK(:,2)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
            ENDIF

            IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
               ! gradient corrections to LDA
               ! GGA potential
               CALL FEXCG(XC,2,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                    KINEDEN,CWORK,DENCOR,CVTOT,MUTOT_,XDMTMP,0,1,SIZE(MUTOT_,1))
            ENDIF

            ! add LDA part of potential
            CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
               CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
               E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
            ! add Coulomb hole
            CALL COHSM1_RGRID(XC, 2, CWORK(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
            ! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

            ! get the proper direction vx = v0 + hat m delta v
            CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)

         ELSEIF (WDES%LNONCOLLINEAR) THEN

            IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
               CALL FEXCG(XC,4,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                    KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDMTMP,0,1,SIZE(MUTOT_,1))
            ENDIF

            IF (XC%LMU) THEN
               ! rearranged from (tot,mag) to (uu,ud,du,dd) representation
               CALL POT_FLIP_RL(MUTOT_, GRIDC, WDES%NCDIJ)
            ENDIF

            ! FEXCF requires (tot,mag) density instead of (uu,ud,du,dd)
            CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)
! quick hack to do LDA+U instead of LSDA+U
            IF (L_NO_LSDA_GLOBAL()) THEN
GPACC !$ACC KERNELS PRESENT(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
               CWORK(:,2)=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
            ENDIF
! end of hack
            ALLOCATE(CVCOL(GRIDC%MPLWV, 2))
GPACC !$ACC ENTER DATA CREATE(CVCOL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CVCOL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CVCOL) __IF_ASYNC__
GPOMP !$OMP TARGET IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
            CVCOL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

            ! add LDA part of potential
            CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
                CWORK(1,1), CWORK(1,2), DENCOR, CVCOL(1,1), CVCOL(1,2), &
                  E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
              ! add Coulomb hole
            CALL COHSM1_RGRID(XC, 2, CWORK(1,1), CVCOL(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)
            ! we have now the potential for up and down stored in CVCOL(:,1) and CVCOL(:,2)
            ! get the proper direction vx = v0 + hat m delta v

            ! Use the ADD routine to add to two component term to CVTOT
            CALL MAG_DIRECTION_ADD(CTMP(1,1), CVTOT(1,1), CVCOL(1,1), GRIDC, WDES%NCDIJ)

GPACC !$ACC EXIT DATA DELETE(CVCOL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CVCOL) IF(OFFLOAD_ON)
            DEALLOCATE(CVCOL)

         ELSE
            IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
               ! gradient corrections to LDA
               CALL FEXCG(XC,1,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                    KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDMTMP,0,1,SIZE(MUTOT_,1))
            ENDIF
            ! LDA part of potential
            CALL FEXCP(GRIDC,XC,LATT_CUR%OMEGA, &
                 CTMP,DENCOR,CVTOT,CWORK,E%CVZERO,EXC,E%XCENC,XCSIF,.TRUE.)
            ! add Coulomb hole
            CALL COHSM1_RGRID(XC, 1, CTMP(1,1), CVTOT(1,1), DENCOR, GRIDC, LATT_CUR%OMEGA, .TRUE.)

         ENDIF

#ifdef libbeef
         IF(LBEEFCALCBASIS) THEN
            IF(BEEFCOUNTER.GT.0) THEN
               IF(BEEFCOUNTER.LE.31) THEN
                  BEEFXC(BEEFCOUNTER) = EXCG
                  BEEFCOUNTER = BEEFCOUNTER + 1
               ELSE
                 BEEFXC(32) = EXCG + BEEFXC(31)
#if defined(MPI) || defined(MPI_CHAIN)
                 IF (COMM%NODE_ME==COMM%IONODE) THEN
#endif
                     IF(.NOT.LBEEFBAS) THEN
                        CALL BEEFRANDINITDEF
                        CALL BEEFENSEMBLE(BEEFXC, BEENERGIES)
                        WRITE(8,*) "BEEFens 2000 ensemble energies"
                        DO BEEFCOUNTER = 1,2000
                           WRITE(8,"(E35.15)") BEENERGIES(BEEFCOUNTER)
                        ENDDO
                     ENDIF

                     WRITE(8,*) "BEEF xc energy contributions"
                     DO BEEFCOUNTER = 1,32
                        WRITE(8,*) BEEFCOUNTER, ": ", BEEFXC(BEEFCOUNTER)
                     ENDDO
#if defined(MPI) || defined(MPI_CHAIN)
                  ENDIF
#endif
                  LBEEFCALCBASIS = .FALSE.
                  BEEFCOUNTER = 0
                  XC%LUSE_VDW = SAVELUSE_VDW
                  CALL BEEFSETMODE(-1)
               ENDIF
               EXC     =0
               E%XCENC =0
               E%EXCG  =0
               E%CVZERO=0
               XCSIF   =0
               CVTOT   =0
               EXCG    =0
               XCENCG  =0
               CVZERG  =0
               TMPSIF  =0

               GOTO 4096
            ENDIF
         ENDIF
#endif

         XCSIF=XCSIF+TMPSIF
         E%EXCG=EXC+EXCG
         E%XCENC=E%XCENC+XCENCG
         E%CVZERO=E%CVZERO+CVZERG

      ELSE
         DO ISP=1,WDES%NCDIJ
            CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
         ENDDO
      ENDIF

      IF (XC%LNO_AUG_XC) THEN
         ! Store CVTOT in CTMP
GPACC !$ACC KERNELS PRESENT(CTMP,CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CTMP,CVTOT) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         CTMP=CVTOT
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ENDIF

      IF (XC%LSFBXC) CALL REMOVE_SOURCES_FROM_BXC(E, CVTOT, CHTOT, GRIDC, WDES, LATT_CUR)

!-MM- changes to accomodate constrained moments
!-----------------------------------------------------------------------
! add constraining potential
!-----------------------------------------------------------------------
#ifndef NGXhalf
#ifndef NGZhalf
      IF (M_CONSTRAINED()) THEN
      ! NB. at this point both CHTOT and CVTOT must be given
      ! in (charge,magnetization) convention in real space
GPACC !$ACC UPDATE SELF(CHTOT,CVTOT) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE FROM(CHTOT,CVTOT) IF(OFFLOAD_ON)
         CALL M_INT(CHTOT,GRIDC,WDES)
         CALL ADD_CONSTRAINING_POT(CVTOT,GRIDC,WDES)
GPACC !$ACC UPDATE DEVICE(CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(CVTOT) IF(OFFLOAD_ON)
      ENDIF
#endif
#endif
!-MM- end of addition

!-----------------------------------------------------------------------
! calculate the total potential
!-----------------------------------------------------------------------
! add external electrostatic potential
      DIP%ECORR=0
      DIP%E_ION_EXTERN=0

      IF (DIP%LCOR_DIP) THEN
         ! get the total charge and store it in CWORK
         IF  ( WDES%NCDIJ > 1) THEN
            CALL MAG_DENSITY(CHTOT,CWORK, GRIDC, WDES%NCDIJ)
         ELSE
            CALL RL_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
         ENDIF

GPACC !$ACC UPDATE SELF(CWORK,CVTOT) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE FROM(CWORK,CVTOT) IF(OFFLOAD_ON)
           CALL CDIPOL(GRIDC, LATT_CUR,P,T_INFO, &
             CWORK,CSTRF,CVTOT(1,1), WDES%NCDIJ, INFO%NELECT )
GPACC !$ACC UPDATE DEVICE(CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(CVTOT) IF(OFFLOAD_ON)

         CALL EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT(1,1))
      ELSE
         CALL EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT(1,1))
      ENDIF
! embedding__
      IF (EXTPT_LEXTPOT()) CALL EXTPT_EXTERNAL_POT_ADD(GRIDC, LATT_CUR, CVTOT(1,1))
! embedding__

      ! Take the charge density back to reciprocal space
      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO

      ! Take the kinetic energy density back to reciprocal space if necessary
      IF (XC%LTAU) THEN
         DO ISP=1,WDES%NCDIJ
            CALL FFT_RC_SCALE(KINEDEN%TAU(1,ISP),KINEDEN%TAU(1,ISP),GRIDC)
            CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, KINEDEN%TAU(1,ISP))
         ENDDO
      ENDIF

      IF (M_CONSTRAINED()) THEN
GPACC !$ACC UPDATE SELF(CHTOT) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE FROM(CHTOT) IF(OFFLOAD_ON)
      ENDIF

!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D(CVTOT(1,ISP),GRIDC,-1)
      ENDDO

      IF (XC%LMU) THEN
!-----------------------------------------------------------------------
! The same procedure to calculate MU from MUTOT
!-----------------------------------------------------------------------
         RINPL=1._q/GRIDC%NPLWV
         DO  ISP=1,WDES%NCDIJ
            CALL RL_ADD(MUTOT_(1,ISP),RINPL,MUTOT_(1,ISP),0._q,MUTOT_(1,ISP),GRIDC)
            CALL FFT3D(MUTOT_(1,ISP),GRIDC,-1)
         ENDDO
!-----------------------------------------------------------------------
! copy MUTOT_ to MU_ and set contribution of unbalanced lattice-vectors
! to zero, then FFT of MU_ and MUTOT_ to real space
!-----------------------------------------------------------------------
         DO ISP=1,WDES%NCDIJ
            CALL SETUNB_COMPAT(MUTOT_(1,ISP),GRIDC)
            CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,MUTOT_(1,ISP),CWORK1)
            CALL SETUNB(CWORK1,GRID_SOFT)
            CALL FFT3D(CWORK1,GRID_SOFT, 1)
            CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,MU_(1,ISP),GRID_SOFT)

       ! final result is only correct for first in-band-group
       ! (i.e. proc with nodeid 1 in COMM_INTER)
       ! copy to other in-band-groups using COMM_INTER
       ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
#ifdef realmode
            CALLMPI( M_bcast_d(COMM_INTER, MU_(1,ISP), GRID%RL%NP))
#else
            CALLMPI( M_bcast_z(COMM_INTER, MU_(1,ISP), GRID%RL%NP))
#endif
            CALL FFT3D(MUTOT_(1,ISP),GRIDC,1)
         ENDDO
      ENDIF
!-----------------------------------------------------------------------
! add the hartree potential and the double counting corrections
!-----------------------------------------------------------------------
      IF (PRESENT(COULOMB_POT)) THEN
         CALL POTHAR(GRIDC, LATT_CUR, CHTOT, CWORK,E%DENC, COULOMB_POT)
GPACC !$ACC UPDATE SELF(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE FROM(CWORK) IF(OFFLOAD_ON)
         CALL GET_RHO0(GRIDC,CWORK,E%POTHARZERO)
      ELSE
         CALL POTHAR(GRIDC, LATT_CUR, CHTOT, CWORK,E%DENC)
      ENDIF
#ifdef PLUGINS
#ifndef NO_LOCAL_POTENTIAL
      IF (INFO%PLUGIN%LOCAL_POTENTIAL) THEN
         ALLOCATE(HARTREE_POTENTIAL(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(HARTREE_POTENTIAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:HARTREE_POTENTIAL) IF(OFFLOAD_ON)
         CALL RC_ADD(CWORK, 1._q, CWORK, 0._q, HARTREE_POTENTIAL, GRIDC)
         CALL FFT3D(HARTREE_POTENTIAL, GRIDC, 1)
      ENDIF
#endif
#endif
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT,CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK,CVTOT) IF(OFFLOAD_ON)
      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
!-----------------------------------------------------------------------
! add external potential in reciprocal space
!-----------------------------------------------------------------------
      CALL ADD_EXTERNAL_POTENTIAL(GRIDC, CVTOT)
! solvation__
!-----------------------------------------------------------------------
! add the dielectric corrections to CVTOT and the energy
!-----------------------------------------------------------------------
      CALL SOL_Vcorrection(INFO,T_INFO,LATT_CUR,P,WDES,GRIDC,CHTOT,CVTOT)
! solvation__
!-----------------------------------------------------------------------
!  add local pseudopotential potential
!-----------------------------------------------------------------------
      IF(INFO%TURBO==0)THEN
         IF (PRESENT(COULOMB_POT)) THEN
            CALL POTION(GRIDC,P,LATT_CUR,T_INFO,CWORK,CWORK1,CSTRF,E%PSCENC,COULOMB_POT)
            CALL GET_RHO0(GRIDC,CWORK,E%POTIONZERO)
         ELSE
            CALL POTION(GRIDC,P,LATT_CUR,T_INFO,CWORK,CWORK1,CSTRF,E%PSCENC)
         ENDIF
      ELSE
__GPU    CALL OFFLOAD_ON_PUSH(.FALSE.)
         IF (PRESENT(COULOMB_POT)) THEN
            CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK,E%PSCENC,E%TEWEN,COULOMB_POT=COULOMB_POT)
         ELSE
            CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK,E%PSCENC,E%TEWEN)
         ENDIF
__GPU    CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(CWORK) IF(OFFLOAD_ON)
      ENDIF

#ifdef PLUGINS
#ifndef NO_LOCAL_POTENTIAL
      IF (INFO%PLUGIN%LOCAL_POTENTIAL) THEN
         ALLOCATE(ION_POTENTIAL(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(ION_POTENTIAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:ION_POTENTIAL) IF(OFFLOAD_ON)
         CALL RC_ADD(CWORK, 1._q, CWORK, 0._q, ION_POTENTIAL, GRIDC)
         CALL FFT3D(ION_POTENTIAL, GRIDC, 1)
      ENDIF
#endif
#endif

      ELECTROSTATIC=0
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRIDC,CWORK,CHTOT) &
GPACC !$ACC& REDUCTION(+:ELECTROSTATIC) PRIVATE(N2,N3,NG,FACTM) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:GRIDC,CWORK,CHTOT) &
GPOMP !$OMP REDUCTION(+:ELECTROSTATIC) PRIVATE(N2,N3,NG,FACTM) IF(OFFLOAD_ON)

      col: DO NC=1,GRIDC%RC%NCOL
__CPU N2= GRIDC%RC%I2(NC)
__CPU N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
__GPU   N2= GRIDC%RC%I2(NC)
__GPU   N3= GRIDC%RC%I3(NC)
        NG = N1 + (NC-1)*GRIDC%RC%NROW
        SETFACT1
        SETFACT

        ELECTROSTATIC=ELECTROSTATIC+ MULFACT CWORK(NG,1)*CONJG(CHTOT(NG,1))
      ENDDO row
      ENDDO col
      ELECTROSTATIC=ELECTROSTATIC+E%PSCENC-E%DENC+E%TEWEN

      E%PSCENC=E%PSCENC + DIP%ECORR + DIP%E_ION_EXTERN

      ! apply self consistent potential correction
      CALL SCPC_APPLY(GRIDC, LATT_CUR, COULOMB_POT_PERIODIC, CHTOT, P, T_INFO, CSTRF, CVTOT, E%ESCPC)

GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT,CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK,CVTOT) IF(OFFLOAD_ON)
      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
! bexternal__
      IF (LBEXTERNAL()) CALL BEXT_ADDV(CVTOT,GRIDC,SIZE(CVTOT,2))
! bexternal__

#ifdef PLUGINS
#ifndef NO_LOCAL_POTENTIAL
      IF (INFO%PLUGIN%LOCAL_POTENTIAL) THEN
         ALLOCATE(TOTAL_CHARGE_DENSITY(GRIDC%MPLWV, WDES%NCDIJ))
         ALLOCATE(PLUGINS_POTENTIAL(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(TOTAL_CHARGE_DENSITY,PLUGINS_POTENTIAL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:TOTAL_CHARGE_DENSITY,PLUGINS_POTENTIAL) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(TOTAL_CHARGE_DENSITY,PLUGINS_POTENTIAL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:TOTAL_CHARGE_DENSITY,PLUGINS_POTENTIAL) IF(OFFLOAD_ON)
         PLUGINS_POTENTIAL = 0.0_q
         TOTAL_CHARGE_DENSITY = 0.0_q
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

         DO ISP = 1, WDES%NCDIJ
            CALL RC_ADD(CHTOT(:,ISP), 1._q, CHTOT(:,ISP), 0._q, TOTAL_CHARGE_DENSITY(:,ISP), GRIDC)
            CALL FFT3D(TOTAL_CHARGE_DENSITY(:,ISP), GRIDC, 1)
         END DO

         IF  ( WDES%NCDIJ > 1) THEN
            CALL MAG_DENSITY(TOTAL_CHARGE_DENSITY, CWORK, GRIDC, WDES%NCDIJ)
         ELSE
            CALL RL_ADD(TOTAL_CHARGE_DENSITY, 1.0_q, TOTAL_CHARGE_DENSITY, 0.0_q, CWORK, GRIDC)
         ENDIF

         IF (DIP%LCOR_DIP) THEN
            ALLOCATE(DIPOLE_MOMENT(3))
            DIPOLE_MOMENT = DIP%DIPOLC
         ELSE
            ALLOCATE(DIPOLE_MOMENT(0))
         END IF

GPACC !$ACC UPDATE SELF(PLUGINS_POTENTIAL,HARTREE_POTENTIAL,ION_POTENTIAL,CWORK) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE FROM(PLUGINS_POTENTIAL,HARTREE_POTENTIAL,ION_POTENTIAL,CWORK) IF(OFFLOAD_ON)
__GPU    CALL OFFLOAD_ON_PUSH(.FALSE.)

         CALL PLUGINS_LOCAL_POTENTIAL(INFO, T_INFO, LATT_CUR, P, GRIDC, CWORK(:,1), &
             PLUGINS_POTENTIAL, HARTREE_POTENTIAL, ION_POTENTIAL, DIPOLE_MOMENT, E%EPLUGINS)

__GPU    CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(PLUGINS_POTENTIAL) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE TO(PLUGINS_POTENTIAL) IF(OFFLOAD_ON)

         CALL FFT_RC_SCALE(PLUGINS_POTENTIAL, PLUGINS_POTENTIAL, GRIDC)
         CALL RC_ADD(CVTOT(:,1), 1._q, PLUGINS_POTENTIAL, 1._q, CVTOT(:,1), GRIDC)

         ! clean up
GPACC !$ACC EXIT DATA DELETE(HARTREE_POTENTIAL,ION_POTENTIAL,PLUGINS_POTENTIAL,TOTAL_CHARGE_DENSITY) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:HARTREE_POTENTIAL,ION_POTENTIAL,PLUGINS_POTENTIAL,TOTAL_CHARGE_DENSITY) IF(OFFLOAD_ON)
         DEALLOCATE(HARTREE_POTENTIAL,ION_POTENTIAL,PLUGINS_POTENTIAL,TOTAL_CHARGE_DENSITY,DIPOLE_MOMENT)
      ENDIF
#endif
#endif

!=======================================================================
! if overlap is used :
! copy CVTOT to SV and set contribution of unbalanced lattice-vectors
! to zero,  then  FFT of SV and CVTOT to real space
!=======================================================================

      CALL POT_FLIP(CVTOT, GRIDC,WDES%NCDIJ )

      DO ISP=1,WDES%NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
         CALL FFT3D(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

         !  final result is only correct for first in-band-group
         ! (i.e. proc with nodeid 1 in COMM_INTER)
         !  copy to other in-band-groups using COMM_INTER
         ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
#ifdef realmode
         CALLMPI( M_bcast_d(COMM_INTER, SV(1,ISP), GRID%RL%NP))
#else
         CALLMPI( M_bcast_z(COMM_INTER, SV(1,ISP), GRID%RL%NP))
#endif
         CALL FFT3D(CVTOT(1,ISP),GRIDC,1)
      ENDDO

      IF (XC%LNO_AUG_XC) THEN
         CALL POT_FLIP(CTMP, GRIDC,WDES%NCDIJ )
         ! Subtract CTMP from CVTOT
         DO ISP=1,WDES%NCDIJ
            CALL RL_ADD(CVTOT(1,ISP),1.0_q,CTMP(1,ISP),-1.0_q,CVTOT(1,ISP),GRIDC)
         ENDDO
      ENDIF

GPACC !$ACC EXIT DATA DELETE(CTMP,CWORK1,CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP,CWORK1,CWORK) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      IF (XC%LNO_AUG_XC) THEN
GPACC !$ACC EXIT DATA DELETE(CHDEN) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CHDEN) IF(OFFLOAD_ON)
         DEALLOCATE(CTMP)
      ENDIF

! Ensure that MU and MUTOT are up-to-date on the host
GPACC !$ACC UPDATE SELF(MU_,MUTOT_) IF(OFFLOAD_ON.AND.XC%LMU) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET UPDATE FROM(MU_,MUTOT_) IF(OFFLOAD_ON.AND.XC%LMU)

GPACC !$ACC EXIT DATA DELETE(MU_,MUTOT_) IF(OFFLOAD_ON.AND.XC%LMU) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:MU_,MUTOT_) IF(OFFLOAD_ON.AND.XC%LMU)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:MU_,MUTOT_) IF(OFFLOAD_ON)
      IF (.NOT.XC%LMU) THEN
GPACC !$ACC EXIT DATA DELETE(MU_,MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:MU_,MUTOT_) IF(OFFLOAD_ON)
         DEALLOCATE(MU_,MUTOT_)
      ENDIF
      NULLIFY(MU_,MUTOT_)

      DEALLOCATE(CWORK1,CWORK)

#ifdef OFFLOADING

GPACC !$ACC EXIT DATA DELETE(KINEDEN%TAU,KINEDEN) IF(ASSOCIATED(KINEDEN%TAU).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC EXIT DATA DELETE(KINEDEN%TAUC,KINEDEN) IF(ASSOCIATED(KINEDEN%TAUC).AND.OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPACC !$ACC EXIT DATA DELETE(CHTOT,LATT_CUR,DENCOR,CSTRF) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:KINEDEN,KINEDEN%TAU) IF(ASSOCIATED(KINEDEN%TAU).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:KINEDEN,KINEDEN%TAUC) IF(ASSOCIATED(KINEDEN%TAUC).AND.OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CHTOT,LATT_CUR,DENCOR,CSTRF) IF(OFFLOAD_ON)
! Ensure that SV and CVTOT are up-to-date on the host
GPACC !$ACC UPDATE SELF(SV,CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET UPDATE FROM(SV,CVTOT) IF(OFFLOAD_ON)

GPACC !$ACC EXIT DATA DELETE(SV,CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:SV,CVTOT) IF(OFFLOAD_ON)
      IF (.NOT.PRESENT(COULOMB_POT)) CALL DELETE_TYPED_VAR(COULOMB_POT_PERIODIC)
      CALL DELETE_TYPED_VAR(GRIDC)
      CALL DELETE_TYPED_VAR(GRID_SOFT)
      CALL DELETE_TYPED_VAR(SOFT_TO_C)
      CALL DELETE_TYPED_VAR(P)
      CALL DELETE_TYPED_VAR(XC)
      CALL OFFLOAD_ON_POP
      ! set ACC_ASYNC_Q back to the original value
GPACC ACC_ASYNC_Q=ASYNC_PRE
#endif

      PROFILING_STOP('potlok')

      RETURN
    END SUBROUTINE POTLOK


!***********************************************************************
!
!> Small helper routine to set the local potential on the coarse
!> plane wave grid from the full dense (augmentation) grid.
!>
!> CVTOT must be supplied in reciprocal space (not usually the case)
!> and SV is returned in real space.
!
!***********************************************************************

    SUBROUTINE SET_SV( GRID, GRIDC, GRID_SOFT, COMM_INTER, SOFT_TO_C, NCDIJ, SV, CVTOT)

      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT NONE

      INTEGER NCDIJ
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (transit)     SOFT_TO_C

      RGRID   SV(DIMREAL(GRID%MPLWV), NCDIJ)
      COMPLEX(q) CVTOT(GRIDC%MPLWV, NCDIJ)
      TYPE (communic)    COMM_INTER
! work arrays
      COMPLEX(q) ::  CWORK1(GRID_SOFT%MPLWV)
      INTEGER ISP


      DO ISP=1,NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
         ! transform to real-space representation
         CALL FFT3D(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

    !  final result is only correct for first in-band-group
    ! (i.e. proc with nodeid 1 in COMM_INTER)
    !  copy to other in-band-groups using COMM_INTER
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
#ifdef realmode
         CALLMPI( M_bcast_d(COMM_INTER, SV(1,ISP), GRID%RL%NP))
#else
         CALLMPI( M_bcast_z(COMM_INTER, SV(1,ISP), GRID%RL%NP))
#endif
         CALL FFT3D(CVTOT(1,ISP),GRIDC,1)
      ENDDO
    END SUBROUTINE SET_SV

    SUBROUTINE ADD_EXTERNAL_POTENTIAL(GRIDC, CVTOT)
      ! adds the external potential stored as module variable
      USE mgrid, ONLY: grid_3d
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(grid_3d), INTENT(IN) :: GRIDC
      COMPLEX(q), INTENT(INOUT) :: CVTOT(:,:)
      !
      IF (.NOT.ALLOCATED(EXTERNAL_POTENTIAL)) RETURN
      !
      ! sanity checks
      IF (SIZE(EXTERNAL_POTENTIAL, 1) < GRIDC%RC%NP) &
         CALL vtutor%bug("Size of external potential incompatible with grid.", __FILE__, __LINE__)
      IF (SIZE(CVTOT, 2) > 1) &
         CALL vtutor%bug("Implementation only tested without spin polarization.", __FILE__, __LINE__)
      !
      CALL ZAXPY(GRIDC%RC%NP, zone, EXTERNAL_POTENTIAL, 1, CVTOT, 1)
      !
    END SUBROUTINE ADD_EXTERNAL_POTENTIAL

    SUBROUTINE ADD_TO_EXTERNAL_POTENTIAL(POTENTIAL)
      ! add another contribution to existing external potential
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      COMPLEX(q), INTENT(IN) :: POTENTIAL(:)
      !
      IF (.NOT.ALLOCATED(EXTERNAL_POTENTIAL)) THEN
         ALLOCATE(EXTERNAL_POTENTIAL(SIZE(POTENTIAL)))
         EXTERNAL_POTENTIAL = POTENTIAL
      ELSE
         IF (SIZE(POTENTIAL) /= SIZE(EXTERNAL_POTENTIAL)) &
            CALL vtutor%bug("Size of potentials inconsistent.", __FILE__, __LINE__)
         CALL ZAXPY(SIZE(POTENTIAL), zone, POTENTIAL, 1, EXTERNAL_POTENTIAL, 1)
      END IF
      !
    END SUBROUTINE ADD_TO_EXTERNAL_POTENTIAL

    SUBROUTINE RESET_EXTERNAL_POTENTIAL()
      ! clear the external potential
      IF (ALLOCATED(EXTERNAL_POTENTIAL)) DEALLOCATE(EXTERNAL_POTENTIAL)
    END SUBROUTINE RESET_EXTERNAL_POTENTIAL


!************************ SUBROUTINE POTXC  ****************************
!
!> This subroutine to calculate the XC-potential including gradient
!> corrections.
!>
!> This routine is required to calculate the partial core
!> corrections to the forces.
!>
!> On entry: \n
!>  CHTOT(:,1)    density \n
!>  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization \n
!> On return (LNONCOLLINEAR=.FALSE.): \n
!>  CVTOT(:,1)    average potential \n
!>  CVTOT(:,2:4)  magnetic field
!
!***********************************************************************

      SUBROUTINE POTXC(HAMILTONIAN,KINEDEN,GRIDC,GRID_SOFT,SOFT_TO_C,INFO,XC,WDES, &
                 LATT_CUR,CHDEN,CHTOT,DENCOR,CVTOT,IVDW,XDM,NPDDSC)

__GPU USE moffload_struct_def

      USE prec
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY: fexcf, fexcp
      USE tau_mu, ONLY : tau_handle
      USE hamil_struct_def, ONLY : ham_handle
      USE mpimy
      USE mgrid
      USE lattice
      USE base
      USE wave

      IMPLICIT NONE

      TYPE (ham_handle)  HAMILTONIAN
      TYPE (tau_handle)  KINEDEN
      TYPE (xc_info)     XC
      TYPE (grid_3d)     GRIDC,GRID_SOFT
      TYPE (wavedes)     WDES
      TYPE (transit)     SOFT_TO_C
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR
      INTEGER IVDW,NPDDSC,ISP
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      COMPLEX(q), TARGET :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ)
      RGRID      DENCOR(GRIDC%RL%NP)
! work arrays
      REAL(q)    XCSIF(3,3),XDM(NPDDSC),EXC,EXCG,RINPL,XCENC,XCENCG
      COMPLEX(q), ALLOCATABLE:: CWORK(:,:)
      COMPLEX(q) CVZERG,CVZERO
      COMPLEX(q), POINTER :: MUTOT_(:,:)

      COMPLEX(q), POINTER :: CTMP(:,:)
      COMPLEX(q), ALLOCATABLE ::  CVCOL(:,:)

      PROFILING_START('potxc')

      HAMILTONIAN%XC = XC

      ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ))

GPACC !$ACC ENTER DATA COPYIN(CHTOT) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:CHTOT) IF(OFFLOAD_ON)

      IF (XC%LNO_AUG_XC) THEN
         ALLOCATE(CTMP(GRIDC%MPLWV,WDES%NCDIJ))
GPACC !$ACC ENTER DATA COPYIN(CHDEN) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:CHDEN) IF(OFFLOAD_ON)
      ELSE
         CTMP => CHTOT
      ENDIF

GPACC !$ACC ENTER DATA CREATE(CTMP,CWORK,CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP,CWORK,CVTOT) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: CTMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVTOT) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CVTOT = 0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      IF (XC%LMU) THEN
         IF (.NOT.ASSOCIATED(HAMILTONIAN%MUTOT)) &
            CALL vtutor%bug("POTXC: MUTOT not associated.", __FILE__, __LINE__)
         MUTOT_ => HAMILTONIAN%MUTOT
GPACC !$ACC ENTER DATA CREATE(MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:MUTOT_) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET ENTER DATA MAP(ALWAYS, ALLOC: MUTOT_) IF(OFFLOAD_ON)
! here the full size of MU and MUTOT has to be given, the cce 19
! compiler thinks the size is 1,1
#ifdef __DCRAYFTN_VERSION
GPOMP !$OMP TARGET MAP(PRESENT,TO:MUTOT_,HAMILTONIAN%MUTOT) IF(OFFLOAD_ON)
         MUTOT_ => HAMILTONIAN%MUTOT(1:GRIDC%MPLWV,1:WDES%NCDIJ)
GPOMP !$OMP END TARGET
#endif

GPACC !$ACC KERNELS PRESENT(MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:MUTOT_) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         MUTOT_=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
      ELSE
         ALLOCATE(MUTOT_(1,WDES%NCDIJ))
GPACC !$ACC ENTER DATA CREATE(MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:MUTOT_) IF(OFFLOAD_ON)
      ENDIF

      IF (XC%LNO_AUG_XC) THEN
         ! We will use the soft charge density instead of the augmented one
GPACC !$ACC KERNELS PRESENT(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CTMP) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         CTMP=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
         DO ISP=1,WDES%NCDIJ
            CALL ADD_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CHDEN(1,ISP),CTMP(1,ISP))
            CALL SETUNB_COMPAT(CTMP(1,ISP),GRIDC)
            CALL FFT3D(CTMP(1,ISP),GRIDC,1)
         ENDDO
      ENDIF

      ! Bring the charge density to real space
      DO ISP=1,WDES%NCDIJ
         CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
      ENDDO

      ! Bring the kinetic energy density to real space if needed for meta-GGA functionals
      IF (XC%LTAU) THEN
         DO ISP=1,WDES%NCDIJ
            CALL FFT3D(KINEDEN%TAU(1,ISP),GRIDC,1)
         ENDDO
#ifndef realmode
GPACC !$ACC KERNELS PRESENT(KINEDEN%TAU) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:KINEDEN%TAU) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         KINEDEN%TAU=REAL(KINEDEN%TAU,KIND=q)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
#endif
      ENDIF

      IF (WDES%ISPIN==2) THEN

         ! get the charge and the total magnetization
         CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)

         IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
            ! gradient corrections to LDA
            CALL FEXCG(XC,2,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                 KINEDEN,CWORK,DENCOR,CVTOT,MUTOT_,XDM,IVDW,NPDDSC,SIZE(MUTOT_,1))
         ENDIF

         ! add LDA part of potential
         CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
            CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
            CVZERO,EXC,XCENC,XCSIF, .TRUE.)
         ! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

         ! get the proper direction vx = v0 + hat m delta v
         CALL MAG_DIRECTION(CTMP(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)

      ELSEIF (WDES%LNONCOLLINEAR) THEN

!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
         IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
            CALL FEXCG(XC,4,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                 KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDM,IVDW,NPDDSC,SIZE(MUTOT_,1))
         ENDIF

         ! FEXCF requires (tot,mag) density instead of (uu,ud,du,dd)
         CALL MAG_DENSITY(CTMP, CWORK, GRIDC, WDES%NCDIJ)

         ALLOCATE(CVCOL(GRIDC%MPLWV, 2))
GPACC !$ACC ENTER DATA CREATE(CVCOL) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CVCOL) IF(OFFLOAD_ON)
GPACC !$ACC KERNELS PRESENT(CVCOL) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CVCOL) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
         CVCOL=0
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC!$ACC END KERNELS
         ! add LDA part of potential
         CALL FEXCF(GRIDC,XC,LATT_CUR%OMEGA, &
            CWORK(1,1), CWORK(1,2), DENCOR, CVCOL(1,1), CVCOL(1,2), &
            CVZERO,EXC,XCENC,XCSIF, .TRUE.)

         ! we now have the LDA potential for up and down stored in CVCOL(:,1) and CVCOL(:,2)
         ! get the proper direction vx = v0 + hat m delta v
         CALL MAG_DIRECTION_ADD(CTMP(1,1), CVTOT(1,1), CVCOL(1,1), GRIDC, WDES%NCDIJ)

GPACC !$ACC EXIT DATA DELETE(CVCOL) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CVCOL) IF(OFFLOAD_ON)
         DEALLOCATE(CVCOL)

      ELSE

         IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
            ! gradient corrections to LDA
            CALL FEXCG(XC,1,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                 KINEDEN,CTMP,DENCOR,CVTOT,MUTOT_,XDM,IVDW,NPDDSC,SIZE(MUTOT_,1))
         ENDIF

         ! LDA part of potential
         CALL FEXCP(GRIDC,XC,LATT_CUR%OMEGA, &
              CTMP,DENCOR,CVTOT,CWORK,CVZERO,EXC,XCENC,XCSIF,.TRUE.)

      ENDIF

      ! Take the charge density back to reciprocal space
      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO

      ! Take the kinetic energy density back to reciprocal space if necessary
      IF (XC%LTAU) THEN
         DO ISP=1,WDES%NCDIJ
            CALL FFT_RC_SCALE(KINEDEN%TAU(1,ISP),KINEDEN%TAU(1,ISP),GRIDC)
            CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, KINEDEN%TAU(1,ISP))
         ENDDO
      ENDIF

!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3D(CVTOT(1,ISP),GRIDC,-1)
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
      ENDDO
GPACC !$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK) IF(OFFLOAD_ON)
      DEALLOCATE(CWORK)

      IF (XC%LMU) THEN
!-----------------------------------------------------------------------
! The same procedure with MUTOT
!-----------------------------------------------------------------------
         DO  ISP=1,WDES%NCDIJ
            CALL RL_ADD(MUTOT_(1,ISP),RINPL,MUTOT_(1,ISP),0._q,MUTOT_(1,ISP),GRIDC)
            CALL FFT3D(MUTOT_(1,ISP),GRIDC,-1)
            CALL SETUNB_COMPAT(MUTOT_(1,ISP),GRIDC)
         ENDDO
GPACC !$ACC EXIT DATA COPYOUT(MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:MUTOT_) IF(OFFLOAD_ON)
IFXW2 GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:MUTOT_) IF(OFFLOAD_ON)
      ELSE
GPACC !$ACC EXIT DATA DELETE(MUTOT_) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:MUTOT_) IF(OFFLOAD_ON)
         DEALLOCATE(MUTOT_)
      ENDIF
      NULLIFY(MUTOT_)

GPACC !$ACC EXIT DATA COPYOUT(CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:CVTOT) IF(OFFLOAD_ON)

GPACC !$ACC EXIT DATA DELETE(CHTOT,CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
IFXW2 !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CHTOT) IF(OFFLOAD_ON)
      IF (XC%LNO_AUG_XC) THEN
GPACC !$ACC EXIT DATA DELETE(CHDEN) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CHDEN) IF(OFFLOAD_ON)
         DEALLOCATE(CTMP)
      ENDIF

      PROFILING_STOP('potxc')

      RETURN
      END SUBROUTINE POTXC

      !> compute individual contributions to the potential and
      !> write potential information to OUTCAR, LOCPOT, and HDF5 file
     SUBROUTINE WRITE_POTENTIAL(IO, KINEDEN, GRID, GRIDC, GRID_SOFT, WDES, &
            INFO, XC_DATA, P, T_INFO, COULOMB_POT, E, LATT_CUR, CHDEN, CHTOT, CSTRF, &
            CVTOT, DENCOR, SV, HAMILTONIAN, SOFT_TO_C, XCSIF, DYN, &
            LMDIM, N_MIX_PAW, CDIJ, CQIJ)

         USE base, ONLY: in_struct, info_struct, energy
         USE fileio, ONLY: OUTPOT
         USE hamil_struct_def, ONLY: ham_handle
         USE mdipol, ONLY: DIP, WRITE_VACUUM_LEVEL
         USE mgrid, ONLY: grid_3d, transit
         USE poscar, ONLY: latt, dynamics, type_info, OUTPOS
         USE pseudo_struct_def, ONLY: potcar
         USE setexm, ONLY: xc_info, XC_DATA_NOXC, PUSH_XC_TYPE
         USE tau_mu, ONLY: tau_handle
         USE wave_struct_def, ONLY: wavedes
         USE tutor, ONLY: vtutor
         USE pot_struct_def
#ifdef VASP_HDF5
         USE vhdf5_struct_def, ONLY: GRP_RESULTS, GRP_POTENTIAL, IH5OUTFILEID
         USE vhdf5, ONLY: VH5_WRITE_POTENTIAL, VH5_WRITE_PAW_DIJ_QIJ
#endif

         IMPLICIT NONE

         TYPE(in_struct), INTENT(IN)  :: IO        !< input output unit handle
         ! begin required for POTLOK
         TYPE(tau_handle), INTENT(IN) :: KINEDEN   !< kinetic energy density handle
         TYPE(grid_3d), INTENT(IN)    :: GRID      !< grid for wavefunctions
         TYPE(grid_3d), INTENT(IN)    :: GRIDC     !< grid for potentials/charge
         TYPE(grid_3d), INTENT(IN)    :: GRID_SOFT !< soft grid for potentials/charge
         TYPE(wavedes), INTENT(IN)    :: WDES      !< wavefunction descriptor
         TYPE(info_struct), INTENT(IN):: INFO      !< info about algo etc.
         TYPE(xc_info), INTENT(INOUT) :: XC_DATA   !< xc settings
         TYPE(potcar), INTENT(IN)     :: P(:)      !< structure to support radial grid
         TYPE(type_info), INTENT(IN)  :: T_INFO    !< ionic positions in superlattice
         TYPE(coulomb_potential), INTENT(IN) :: COULOMB_POT
         TYPE(energy), INTENT(INOUT)  :: E         !< energy handle for double counting, band energies, etc.
         TYPE(latt), INTENT(IN)       :: LATT_CUR  !< lattice parameter
         COMPLEX(q), INTENT(IN)       :: CHDEN(:,:)!< pseudo charge density on soft grid
         COMPLEX(q), INTENT(IN)       :: CHTOT(:,:)!< total charge density on fine grid
         COMPLEX(q), INTENT(IN)       :: CSTRF(:,:)!< structure factor
         COMPLEX(q), INTENT(INOUT)    :: CVTOT(:,:)!< local potential
         RGRID, INTENT(INOUT)         :: DENCOR(:) !< partial core charge
         RGRID, INTENT(INOUT)         :: SV(:,:)   !< total local potential on the soft grid
         TYPE(ham_handle), INTENT(IN) :: HAMILTONIAN !< handle for MGGA potential MU and MUTOT
         TYPE(transit), INTENT(IN)    :: SOFT_TO_C !< index table between GRID_SOFT and GRIDC
         REAL(q), INTENT(INOUT)       :: XCSIF(3,3)!< stress tensor from XC
         OVERLAP, INTENT(IN)          :: CQIJ(:,:,:,:)    ! overlap of PP
         ! end required for POTLOK
         TYPE(dynamics), INTENT(IN)   :: DYN       !< contains position of ions
         INTEGER, INTENT(IN)          :: LMDIM     !< no. of augmentation channels
         INTEGER, INTENT(IN)          :: N_MIX_PAW !< number of elements of the augmentation occupancies mixed on the local node
         OVERLAP, INTENT(INOUT)       :: CDIJ(:,:,:,:) !< strength parameters
         ! local
         REAL(q) DUMMY(1)

         IF (WDES%COMM_KINTER%NODE_ME /= 1) RETURN

         IF (SHOULD_WRITE_POT(IO%LH5, IO%WRT_POTENTIAL)) THEN
            CALL OUTPOT_TO_POT("POT", .TRUE., CVTOT)
         ELSE IF (SHOULD_VASPWAVE_HAVE_POT(IO%LH5, IO%WRT_POTENTIAL)) THEN
            CALL OUTPOT_TO_VASPWAVE
         END IF

         CALL OUTPOT_TO_VASPOUT

         IF (SHOULD_WRITE_LOCPOT(IO%WRT_POTENTIAL)) &
            CALL OUTPOT_TO_LOCPOT

         IF (DIP%IDIPCO > 0 .AND. DIP%LCOR_DIP) &
            CALL OUTPUT_VACUUM_LEVEL
#ifdef VASP_HDF5
         IF (IO%WRT_POTENTIAL%PAW) THEN
             CALL VH5_WRITE_PAW_DIJ_QIJ( IH5OUTFILEID, WDES, T_INFO, P, LMDIM, CDIJ, CQIJ )
         ENDIF
#endif

      CONTAINS


      !> compute the local potential write it to LOCPOT
      SUBROUTINE OUTPOT_TO_LOCPOT
         IF (IO%WRT_POTENTIAL%LVHAR) CALL PUSH_XC_TYPE(XC_DATA,XC_DATA_NOXC)
         CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES, &
              INFO,XC_DATA,P,T_INFO,E,LATT_CUR,CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF,COULOMB_POT)
         CALL OUTPOT_TO_POT("LOCPOT", .FALSE., CVTOT)
      END SUBROUTINE OUTPOT_TO_LOCPOT


      SUBROUTINE OUTPUT_VACUUM_LEVEL
         IF (DIP%VACPOTAV .and. IO%WRT_POTENTIAL%LVTOT) THEN
            CALL vtutor%alert("You have requested that vacuum potentials are determined &
            &using the total potential (LVTOT = .True.), as opposed to just the &
            &hartree and ionic contributions. This selection requires large vacuum in your &
            &cell and a tightly converged calculations to work. Consider using &
            &LVHAR = .True. for a more rapid determination of the vacuum potentials.")
         ENDIF
         ! call the dipol routine without changing the potential
         ! It only makes sense to write out the vacuum potential
         ! if corrections to the potential are requested
         DIP%LCOR_DIP = .FALSE.
         CALL CDIPOL_CHTOT_REC(GRIDC, LATT_CUR,P,T_INFO, &
            CHTOT,CSTRF,CVTOT, WDES%NCDIJ, INFO%NELECT, E%PSCENC )
         IF (DIP%IDIPCO < 4) THEN
            CALL WRITE_VACUUM_LEVEL(IO%IU6)
         END IF
         DIP%LCOR_DIP = .TRUE.
      END SUBROUTINE OUTPUT_VACUUM_LEVEL


      !> write the information to restart from a potential to vaspwave
      SUBROUTINE OUTPOT_TO_VASPWAVE
#ifdef VASP_HDF5
         USE pawm, ONLY: SET_RHO_PAW
         USE vhdf5, ONLY: IH5WAVEFILEID, VH5_WRITE_PAW_OCCUPANCIES, GRP_POTENTIAL, VH5_DELETE, VH5_ERROR
         REAL(q) DLM_EXX(N_MIX_PAW,WDES%NCDIJ)
         INTEGER SPIN
         !
         VH5_CHECK(VH5_DELETE(IH5WAVEFILEID,"potential"))
         !
         CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)
         !
         CALL VH5_WRITE_POTENTIAL(IH5WAVEFILEID, "total", GRIDC, CVTOT, GROUP="potential")
         CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
              CDIJ, DLM_EXX)
         DO SPIN = 1, WDES%NCDIJ
            CALL VH5_WRITE_PAW_OCCUPANCIES(IH5WAVEFILEID, GRP_POTENTIAL, P, T_INFO, &
               INFO%LOVERL, DLM_EXX(:,SPIN), GRIDC%COMM, SPIN)
         END DO
         !
         CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)
#endif
      END SUBROUTINE OUTPOT_TO_VASPWAVE

      !> write locpot to vaspwave 
      SUBROUTINE LOCPOT_TO_VASPWAVE
#ifdef VASP_HDF5
         USE pawm, ONLY: SET_RHO_PAW
         USE vhdf5
         IMPLICIT NONE
         INTEGER :: SPIN
         REAL(q) DLM_EXX(N_MIX_PAW,WDES%NCDIJ)

         ! delete LOCPOT group if already present
         VH5_CHECK(VH5_DELETE(IH5WAVEFILEID,GRP_LOCPOT))

         ! dump lattice information
         CALL VH5_WRITE_LATTICE_ION(IH5WAVEFILEID," ",T_INFO,LATT_CUR%SCALE,LATT_CUR%A,&
                      .FALSE.,DYN%POSION,GROUP=GRP_LOCPOT,SUBGROUP="position")

         IF (IO%WRT_POTENTIAL%LVHAR) CALL PUSH_XC_TYPE(XC_DATA,XC_DATA_NOXC)
         CALL POTLOK(HAMILTONIAN,KINEDEN,GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,XC_DATA,P,T_INFO,E,LATT_CUR, &
                  CHDEN,CHTOT,CSTRF,CVTOT,DENCOR,SV,SOFT_TO_C,XCSIF,COULOMB_POT)

         ! the potential is complex for noncollinear calculations unless represented
         ! in terms of the Pauli matrices
         IF (WDES%NCDIJ == 4) &
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)

         ! at the moment the spin up and down potential is written to the file
         CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
              CDIJ, DLM_EXX)

         ! dump potential to vaspwave.h5 in group: /locpot/total
         CALL VH5_WRITE_GRID_QUANTITY(IH5WAVEFILEID, GRP_LOCPOT, 'total', GRIDC, CVTOT)
         DO SPIN = 1, WDES%NCDIJ
            CALL VH5_WRITE_PAW_OCCUPANCIES(IH5WAVEFILEID, GRP_LOCPOT, P, T_INFO, &
               INFO%LOVERL, DLM_EXX(:,SPIN), GRIDC%COMM, SPIN)
         END DO

         IF (WDES%NCDIJ == 4) &
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)
#endif
      END SUBROUTINE LOCPOT_TO_VASPWAVE

      !> compute the local potentials and write them to the HDF5 file
      SUBROUTINE OUTPOT_TO_VASPOUT
#ifdef VASP_HDF5
         COMPLEX(q), ALLOCATABLE :: CHARGE(:), WORK(:)
         INTEGER SPIN
         !
         IF (SHOULD_VASPOUT_HAVE_TOTAL_POT(IO%WRT_POTENTIAL)) THEN
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "total", GRIDC, CVTOT)
            CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)
         END IF
         !
         IF (SHOULD_VASPOUT_HAVE_MU_POT(IO%WRT_POTENTIAL)) THEN
            IF (.NOT.XC_DATA%LTAU) THEN
               CALL vtutor%alert("Exchange-correlation potential associated to the kinetic energy density "//&
               "(xcmu) requested in vaspout.h5, but "//&
               "the XC functional used does not depend on the kinetic energy density. "//&
               "No xcmu potential will be written.")
            ELSE
               CALL RC_FLIP_POTENTIAL(HAMILTONIAN%MUTOT, GRIDC, WDES%NCDIJ, .FALSE.)
               CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "xcmu", GRIDC, HAMILTONIAN%MUTOT)
               CALL RC_FLIP_POTENTIAL(HAMILTONIAN%MUTOT, GRIDC, WDES%NCDIJ, .TRUE.)
            END IF
         END IF
         !
         IF (SHOULD_VASPOUT_HAVE_XC_POT(IO%WRT_POTENTIAL)) THEN
            CALL POTXC(HAMILTONIAN, KINEDEN, GRIDC, GRID_SOFT, SOFT_TO_C, INFO, XC_DATA, WDES, &
                 LATT_CUR, CHDEN, CHTOT, DENCOR, CVTOT, 0, DUMMY, 1)
            DO SPIN = 1, SIZE(CVTOT, 2)
               CALL SETUNB_COMPAT(CVTOT(:,SPIN), GRIDC)
               CALL FFT3D(CVTOT(:,SPIN), GRIDC, 1)
            END DO
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "xc", GRIDC, CVTOT)
         END IF
         !
         IF (SHOULD_VASPOUT_HAVE_HARTREE_POT(IO%WRT_POTENTIAL)) THEN
            ALLOCATE(CHARGE, SOURCE=CHTOT(:,1))
            CALL FFT3D(CHARGE, GRIDC, 1)
            CALL FFT_RC_SCALE(CHARGE, CHARGE, GRIDC)
            CALL SETUNB_COMPAT(CHARGE, GRIDC)
            CALL POTHAR(GRIDC, LATT_CUR, CHARGE, CVTOT, DUMMY(1), COULOMB_POT)
            CALL SETUNB_COMPAT(CVTOT(:,1), GRIDC)
            CALL FFT3D(CVTOT(:,1), GRIDC, 1)
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "hartree", GRIDC, CVTOT(:,1:1))
         END IF
         !
         IF (SHOULD_VASPOUT_HAVE_IONIC_POT(IO%WRT_POTENTIAL)) THEN
            ALLOCATE(WORK(MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV)))
            CALL POTION(GRIDC, P, LATT_CUR, T_INFO, CVTOT, WORK, CSTRF, DUMMY(1), COULOMB_POT)
            CALL SETUNB_COMPAT(CVTOT(:,1), GRIDC)
            CALL FFT3D(CVTOT(:,1), GRIDC, 1)
            CALL VH5_WRITE_POTENTIAL(IH5OUTFILEID, "ionic", GRIDC, CVTOT(:,1:1))
         END IF
#endif
      END SUBROUTINE OUTPOT_TO_VASPOUT


      !> small routine to write the potential in POT format
      SUBROUTINE OUTPOT_TO_POT(FILENAME, WRITE_PAW, CVTOT)
         USE main_mpi, ONLY: DIR_APP, DIR_LEN
         USE pawm, ONLY: SET_RHO_PAW, WRT_RHO_PAW
         USE iso_c_binding, ONLY: c_loc, c_f_pointer 
#ifdef VASP_HDF5
         USE vhdf5
         USE string
#endif
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: FILENAME
         LOGICAL, INTENT(IN) :: WRITE_PAW
         COMPLEX(q), INTENT(INOUT), TARGET :: CVTOT(:,:)!< local potential
         REAL(q) DLM(N_MIX_PAW,WDES%NCDIJ)
         ! Local
         INTEGER :: ISP, I, NODE_ME, IONODE
         RGRID, POINTER :: RVTOT(:)
         CHARACTER(LEN=40) :: DATASETNAME

         NODE_ME = WDES%COMM%NODE_ME
         IONODE = WDES%COMM%IONODE

         IF ( .NOT. IO%LH5 ) THEN
            ! the potential is complex for noncollinear calculations unless represented
            ! in terms of the Pauli matrices
            IF (WDES%NCDIJ == 4) &
               CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)

            io_begin

            IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE=DIR_APP(1:DIR_LEN) // FILENAME,STATUS='UNKNOWN')
            REWIND IO%IUVTOT
            ! write lattice parameters and positions to specified unit
            CALL OUTPOS(IO%IUVTOT,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
            io_end

            ! at the moment the spin up and down potential is written to the file
            CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                 CDIJ, DLM)

            CALL C_F_POINTER( C_LOC( CVTOT(:,1) ), RVTOT, [ GRIDC%RL%NP ] )
            CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,RVTOT)
            IF (WRITE_PAW) &
               CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM(:,1), GRIDC%COMM, IO%IUVTOT)

            DO ISP = 2, WDES%NCDIJ
               do_io WRITE( IO%IUVTOT,'(5E20.12)') (T_INFO%ATOMOM(I),I=1,T_INFO%NIONS)
               CALL C_F_POINTER( C_LOC( CVTOT(:,ISP) ), RVTOT, [ GRIDC%RL%NP ] )
               CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,RVTOT)
               IF (WRITE_PAW) &
                  CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM(:,ISP), GRIDC%COMM, IO%IUVTOT )
            END DO
            io_begin
            CLOSE(IO%IUVTOT)
            io_end

            IF (WDES%NCDIJ == 4) &
               CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)

#ifdef VASP_HDF5
         ELSEIF( IO%LH5 ) THEN
            ! delete POT group from vaspwave.h5 if already present
            VH5_CHECK(VH5_DELETE(IH5WAVEFILEID,LOWERCASE(TRIM(FILENAME))))

            ! the potential is complex for noncollinear calculations unless represented
            ! in terms of the Pauli matrices
            IF (WDES%NCDIJ == 4) &
               CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)

            ! dump lattice information
            CALL VH5_WRITE_LATTICE_ION(IH5WAVEFILEID," ",T_INFO,LATT_CUR%SCALE,LATT_CUR%A,&
                      .FALSE.,T_INFO%POSION,GROUP=lowercase( TRIM(FILENAME) ) ,SUBGROUP="position")

            ! at the moment the spin up and down potential is written to the file
            CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                 CDIJ, DLM)
            CALL VH5_WRITE_GRID_QUANTITY(IH5WAVEFILEID, LOWERCASE( TRIM(FILENAME) ), "total", GRIDC, CVTOT)

            DO ISP = 1, WDES%NCDIJ
               IF (WRITE_PAW) &
                  CALL VH5_WRITE_PAW_OCCUPANCIES(IH5WAVEFILEID, LOWERCASE( TRIM(FILENAME) ), P, T_INFO, &
                     INFO%LOVERL, DLM(:,ISP), GRIDC%COMM, ISP)
            END DO
            !
            IF (WDES%NCDIJ == 4) &
               CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)
#endif 
         ENDIF

      END SUBROUTINE OUTPOT_TO_POT

      END SUBROUTINE WRITE_POTENTIAL

      FUNCTION WRT_POTENTIAL_TO_STRING(WRT_POTENTIAL) RESULT (RES)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         CHARACTER(LEN=:), ALLOCATABLE :: RES
         !
         IF (SHOULD_VASPOUT_HAVE_TOTAL_POT(WRT_POTENTIAL)) THEN
            RES = "total"
         ELSE
            RES = ""
         END IF
         
#ifdef __DCRAYFTN_VERSION
         IF (SHOULD_VASPOUT_HAVE_HARTREE_POT(WRT_POTENTIAL)) THEN
            RES = RES // " hartree"
            RES = TRIM(RES)
         END IF
         IF (SHOULD_VASPOUT_HAVE_IONIC_POT(WRT_POTENTIAL)) THEN
            RES = RES // " ionic"
            RES = TRIM(RES)
         END IF
         IF (SHOULD_VASPOUT_HAVE_XC_POT(WRT_POTENTIAL)) THEN
            RES = RES // " xc"
            RES = TRIM(RES)
         END IF
         IF (SHOULD_VASPOUT_HAVE_MU_POT(WRT_POTENTIAL)) THEN
            RES = RES // " xcmu"
            RES = TRIM(RES)
         END IF
         IF (WRT_POTENTIAL%PAW) THEN
            RES = RES // " paw"
            RES = TRIM(RES)
         END IF
#else
         IF (SHOULD_VASPOUT_HAVE_HARTREE_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " hartree")
         IF (SHOULD_VASPOUT_HAVE_IONIC_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " ionic")
         IF (SHOULD_VASPOUT_HAVE_XC_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " xc")
         IF (SHOULD_VASPOUT_HAVE_MU_POT(WRT_POTENTIAL)) &
            RES = TRIM(RES // " xcmu")
         IF (WRT_POTENTIAL%PAW) &
            RES = TRIM(RES // " paw")
#endif
         IF (RES == "") &
            RES = "None"
      END FUNCTION WRT_POTENTIAL_TO_STRING

      LOGICAL FUNCTION SHOULD_WRITE_POT(LH5, WRT_POTENTIAL)
         USE chi_glb, ONLY: LOEP
         USE fock_glb, ONLY: EXXOEP
         LOGICAL, INTENT(IN) :: LH5
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_WRITE_POT = .NOT.LH5 .AND. WRT_POTENTIAL%LVTOT .AND. (LOEP .OR. EXXOEP > 0)
      END FUNCTION SHOULD_WRITE_POT

      LOGICAL FUNCTION SHOULD_VASPWAVE_HAVE_POT(LH5, WRT_POTENTIAL)
         USE chi_glb, ONLY: LOEP
         USE fock_glb, ONLY: EXXOEP
         LOGICAL, INTENT(IN) :: LH5
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPWAVE_HAVE_POT = LH5 .AND. WRT_POTENTIAL%LVTOT .AND. (LOEP .OR. EXXOEP > 0)
      END FUNCTION SHOULD_VASPWAVE_HAVE_POT

      LOGICAL FUNCTION SHOULD_WRITE_LOCPOT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_WRITE_LOCPOT = WRT_POTENTIAL%TOTAL .OR. WRT_POTENTIAL%LVTOT .OR. WRT_POTENTIAL%LVHAR
      END FUNCTION SHOULD_WRITE_LOCPOT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_TOTAL_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_TOTAL_POT = WRT_POTENTIAL%TOTAL .OR. WRT_POTENTIAL%LVTOT
      END FUNCTION SHOULD_VASPOUT_HAVE_TOTAL_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_HARTREE_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_HARTREE_POT = WRT_POTENTIAL%HARTREE .OR. WRT_POTENTIAL%LVHAR
      END FUNCTION SHOULD_VASPOUT_HAVE_HARTREE_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_IONIC_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_IONIC_POT = WRT_POTENTIAL%IONIC .OR. WRT_POTENTIAL%LVHAR
      END FUNCTION SHOULD_VASPOUT_HAVE_IONIC_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_XC_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_XC_POT = WRT_POTENTIAL%XC
      END FUNCTION SHOULD_VASPOUT_HAVE_XC_POT

      LOGICAL FUNCTION SHOULD_VASPOUT_HAVE_MU_POT(WRT_POTENTIAL)
         TYPE(WRITE_POT), INTENT(IN) :: WRT_POTENTIAL
         SHOULD_VASPOUT_HAVE_MU_POT = WRT_POTENTIAL%XCMU
      END FUNCTION SHOULD_VASPOUT_HAVE_MU_POT

    END MODULE


!************************ SUBROUTINE  MAG_DIRECTION  *******************
!
!> On entry CVTOT must contain the v_xc(up) and v_xc(down).
!>
!> On return CVTOT contains \n
!> Collinear case: \n
!>  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2 \n
!>  CVTOT(:,2) =  (v_xc(up) - v_xc(down))/2 \n
!> Non collinear case: \n
!>  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2 \n
!>  CVTOT(:,2) =  hat m_x (v_xc(up) - v_xc(down))/2 \n
!>  CVTOT(:,3) =  hat m_y (v_xc(up) - v_xc(down))/2 \n
!>  CVTOT(:,4) =  hat m_z (v_xc(up) - v_xc(down))/2 \n
!> where hat m is the unit vector of the local magnetization density.
!
!***********************************************************************

#undef strict_compatibility

      SUBROUTINE MAG_DIRECTION(CHTOT, CVTOT, GRID, NCDIJ)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ

      RGRID CHTOT(DIMREAL(GRID%MPLWV), NCDIJ), &
            CVTOT(DIMREAL(GRID%MPLWV), NCDIJ)
      ! local variables
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0

      IF (NCDIJ==2) THEN
#ifndef strict_compatibility
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(V0,DELTAV) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT) &
GPOMP !$OMP PRIVATE(V0,DELTAV) IF(OFFLOAD_ON)
#else
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT,CHTOT) PRIVATE(V0,DELTAV) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT,CHTOT) &
GPOMP !$OMP PRIVATE(V0,DELTAV) IF(OFFLOAD_ON)
#endif
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
            CVTOT(K,1) = V0
#ifndef strict_compatibility
            CVTOT(K,2) = DELTAV
#else
            CVTOT(K,2) = DELTAV*REAL(CHTOT(K,2),KIND=q) / MAX(ABS(CHTOT(K,2)),1.E-20_q)
#endif
         ENDDO
      ELSE IF (NCDIJ==4) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT,CHTOT) PRIVATE(V0,DELTAV,NORM2) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT,CHTOT) &
GPOMP !$OMP PRIVATE(V0,DELTAV,NORM2) IF(OFFLOAD_ON)

         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
#ifdef realmode
            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3)+ CHTOT(K,4)*CHTOT(K,4))),1.E-20_q)
#else
            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CONJG(CHTOT(K,2))+ &
           &           CHTOT(K,3)*CONJG(CHTOT(K,3)) + CHTOT(K,4)*CONJG(CHTOT(K,4)))),1.E-20_q)
#endif

            CVTOT(K,1) = V0
            CVTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            CVTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            CVTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
         ENDDO

      ELSE
         CALL vtutor%bug("internal error: MAG_DIRECTION called with NCDIJ= " // str(NCDIJ), __FILE__, __LINE__)
      ENDIF

      END SUBROUTINE MAG_DIRECTION

!
! similar to above routine, but add two components results from CHTOT to CVTOT
!

      SUBROUTINE MAG_DIRECTION_ADD(CHTOT, CVTOT, CVCOL, GRID, NCDIJ )

__GPU USE moffload_struct_def

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)      GRID
      INTEGER NCDIJ

      RGRID CHTOT(DIMREAL(GRID%MPLWV), 4), &
            CVTOT(DIMREAL(GRID%MPLWV), 4), &
            CVCOL(DIMREAL(GRID%MPLWV), 2)

      ! local variables
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0

      IF (NCDIJ /=4) &
         CALL vtutor%bug("MAG_DIRECTION_ADD called with NCDIJ /= 4 "//str(NCDIJ), __FILE__, __LINE__)

GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT,CHTOT,CVCOL) PRIVATE(V0,DELTAV,NORM2) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT,CHTOT,CVCOL) &
GPOMP !$OMP PRIVATE(V0,DELTAV,NORM2) IF(OFFLOAD_ON)

      DO K=1,GRID%RL%NP
         V0    =(CVCOL(K,1)+CVCOL(K,2))/2 ! Use CVCOL
         DELTAV=(CVCOL(K,1)-CVCOL(K,2))/2 ! Use CVCOL
#ifdef realmode
         NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3)+ CHTOT(K,4)*CHTOT(K,4))),1.E-20_q)
#else
         NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CONJG(CHTOT(K,2))+ CHTOT(K,3)*CONJG(CHTOT(K,3)) + CHTOT(K,4)*CONJG(CHTOT(K,4)))),1.E-20_q)
#endif

         CVTOT(K,1) = CVTOT(K,1) + V0
         CVTOT(K,2) = CVTOT(K,2) + DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
         CVTOT(K,3) = CVTOT(K,3) + DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
         CVTOT(K,4) = CVTOT(K,4) + DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
      ENDDO

      END SUBROUTINE MAG_DIRECTION_ADD


!************************ SUBROUTINE MAG_DENSITY ***********************
!
!> This subroutine calculates the total charge density and the
!> absolute magnitude of the magnetization density.
!>
!> On entry: \n
!>  CHTOT  rho, m_x, m_y, m_z \n
!> On exit: \n
!>  CWORK  rho, sqrt(m_x^2 + m_y^2 + m_z^2) \n
!> in the collinear case, it this means a simple copy CHTOT to CWORK.
!
!***********************************************************************

      SUBROUTINE MAG_DENSITY(CHTOT, CWORK, GRID, NCDIJ)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ

      RGRID CHTOT(DIMREAL(GRID%MPLWV), NCDIJ), &
            CWORK(DIMREAL(GRID%MPLWV), NCDIJ)
      ! local
      INTEGER K


      IF (NCDIJ==2) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(CWORK,CHTOT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK,CHTOT) IF(OFFLOAD_ON)
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)
#ifndef strict_compatibility
            CWORK(K,2)=CHTOT(K,2)
#else
            CWORK(K,2)=ABS(CHTOT(K,2))
#endif
         ENDDO
      ELSE IF (NCDIJ==4) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(CWORK,CHTOT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK,CHTOT) IF(OFFLOAD_ON)
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)
            CWORK(K,2)=SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3) + CHTOT(K,4)*CHTOT(K,4)))
         ENDDO
      ELSE
         CALL vtutor%bug("internal error: MAG_DENSITY called with NCDIJ= " // str(NCDIJ), __FILE__, __LINE__)
      ENDIF
      END SUBROUTINE MAG_DENSITY

!******************************** REMOVE_SOURCES_FROM_BXC  **********************************
!> @brief Remove sources from the XC B-field as proposed in DOI:10.1021/acs.jctc.7b01049
! CVTOT must be given in (charge,magnetization) convention in real space
! NB: Do this before the constraining potential and the external potential is applied.
!********************************************************************************************
      SUBROUTINE REMOVE_SOURCES_FROM_BXC(E, CVTOT, CHTOT, GRIDC, WDES, LATT_CUR)
        USE wave_struct_def, only: wavedes
        USE mgrid, only: grid_3d
        USE base, only: energy
        USE lattice
        IMPLICIT NONE
        TYPE (latt),   intent(in) :: LATT_CUR !< lattice parameter
        TYPE(wavedes), intent(in) :: WDES     !< wave-function descriptor in supercell
        TYPE(grid_3d), intent(in) :: GRIDC    !< real space grid
        COMPLEX(q),    intent(in) :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ) !< charge density in (rho, mag)
        COMPLEX(q),    intent(inout) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) !< xc potential in real space
        TYPE (energy), intent(inout) :: E        !< energy
        ! local
        REAL(q) :: XCENC_SPT, XCENC_BEFORE

        ! XCENC_BEFORE = -<rho*v0>-<m*Bxc>
        CALL SPT_DCC(CVTOT, CHTOT, GRIDC, WDES, LATT_CUR, XCENC_BEFORE)
        ! update CVTOT
        CALL REMOVE_SOURCES(CVTOT, GRIDC, WDES, LATT_CUR)
        ! XCENC_SPT = -<rho*v0>-<m*Bsfxc> - XCENC_BEFORE
        CALL SPT_DCC(CVTOT, CHTOT, GRIDC, WDES, LATT_CUR, XCENC_SPT)
        XCENC_SPT = XCENC_SPT-XCENC_BEFORE
        ! add double counting corrections due to spin torque to E%XCENC
        E%XCENC = E%XCENC + XCENC_SPT

      END SUBROUTINE REMOVE_SOURCES_FROM_BXC

!***************************    SPT_DCC    *****************************
! Compute double-counting corrections for a general noncollinear run.
!
!                     1  / rho+m_z  m_x-im_y \   / v_0+B_z  B_x-iBy \
! - 0.5 * Tr[n*v] = - – |                    | * |                  |
!                     2  \ m_x+im_y rho-m_z  /   \ B_x+iBy  v_0-B_z /
!
!                 = - rho*v_0 - m_x*B_x - m_y*B_y - m_z*B_z
!
! Specifically Bxc || m is not assumed, such that spin torque is accounted for.
!
!> @brief: Return double-counting corrections -v0*rho + Bxc*m for ncl case in XCENC_SPT
!***********************************************************************
      SUBROUTINE SPT_DCC(CVTOT, CHTOT, GRIDC, WDES, LATT_CUR, XCENC_SPT)
__GPU USE moffload_struct_def
      USE prec
      USE lattice
      USE mpimy
      USE mgrid, only : grid_3d
      USE wave_struct_def, only: wavedes
      USE constant

      IMPLICIT NONE
      TYPE(grid_3d), intent(in) :: GRIDC    !< real space grid
      TYPE(wavedes), intent(in) :: WDES     !< wave-function descriptor in supercell
      TYPE (latt),   intent(in) :: LATT_CUR !< lattice parameter
      COMPLEX(q),    intent(in) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) !< xc potential in real space
      COMPLEX(q),    intent(in) :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ) !< charge density in (rho, mag)
      REAL(q),    intent(inout) :: XCENC_SPT  !< -v0*rho + Bxc*m (double counting correction)
      !-----------------------------------------------------------------
      INTEGER :: NODE_ME, IONODE, IDUMP, I, NP
      REAL(q) :: CHGMIN=1E-10
      REAL(q) :: RINPL
      REAL(q) :: RHO, MX, MY, MZ
      REAL(q) :: V0, VX, VY, VZ
! set to one for error-dumps
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
      IDUMP=0
#ifdef debug
      IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
      IDUMP=0
#ifdef debug
      IDUMP=1
#endif
#endif

      XCENC_SPT = 0._q
      NP = GRIDC%RL%NP
      ! Loop over local number of grid points in real space
GPACC !$ACC PARALLEL LOOP PRESENT(CHTOT,CVTOT) PRIVATE(RHO,MX,MY,MZ,V0,VX,VY,VZ) REDUCTION(+:XCENC_SPT) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CHTOT,CVTOT) &
GPOMP !$OMP PRIVATE(RHO,MX,MY,MZ,V0,VX,VY,VZ) REDUCTION(+:XCENC_SPT) IF(OFFLOAD_ON)

      DO I=1,NP
        RHO = CHTOT(I,1)
        MX  = CHTOT(I,2)
        MY  = CHTOT(I,3)
        MZ  = CHTOT(I,4)

        V0  = CVTOT(I,1)
        VX  = CVTOT(I,2)
        VY  = CVTOT(I,3)
        VZ  = CVTOT(I,4)

        XCENC_SPT = XCENC_SPT - RHO*V0 - VX*MX - VY*MY - VZ*MZ
      ENDDO

      XCENC_SPT = XCENC_SPT/GRIDC%NPLWV
      CALLMPI(M_sum_d(GRIDC%COMM, XCENC_SPT, 1))

#ifdef debug
! Test dumps:
      IF (IDUMP/=0) THEN
         WRITE(*,'(A,F24.14)') '<rho*v0>+<m*Bxc> =', -XCENC_SPT
      ENDIF
#endif

      END SUBROUTINE SPT_DCC

!********************* SUBROUTINE REMOVE_SOURCES ***********************
!> @brief: Remove the sources from the xc B field
!   Assuming the Helmholtz decomposition reads
!
!   B(r) = - (1/4pi) grad PHI(r) + rot A(r).
!
!   Then, the magnetic potential is obtained solving the scalar
!   Poisson equation:
!
!   nabla^2 PHI(r) = 4pi div B(r)
!           PHI(G) = 4pi*imag (G dot B)/G^2
!
!   On the plane-wave grid, we find the source free B field using FFTs
!
!   Bsf(r) = B(r) + (1/4pi) grad PHI(r)
!   Bsf(G) = B(G) -          G * (G*B(G))/G^2
!
!   Finally, CVTOT is updated.
!***********************************************************************
      SUBROUTINE REMOVE_SOURCES(CVTOT, GRIDC, WDES, LATT_CUR)
__GPU USE moffload_struct_def 
      USE prec
      USE lattice
      USE mpimy
      USE mgrid, only : grid_3d
      USE tutor, only: vtutor
      USE wave_struct_def, only: wavedes
      USE constant
      USE setexm

      IMPLICIT NONE
      TYPE(grid_3d), intent(in) :: GRIDC    !< real space grid
      TYPE(wavedes), intent(in) :: WDES     !< wave-function descriptor in supercell
      TYPE (latt),   intent(in) :: LATT_CUR !< lattice parameter
      COMPLEX(q), intent(inout) :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) !< xc potential in real space
      !-----------------------------------------------------------------
      COMPLEX(q)   :: CWORK_B(GRIDC%MPLWV,WDES%NCDIJ-1) !< components of B
      COMPLEX(q)   :: CWORK1(GRIDC%MPLWV)   !< [G*B(G)]/G^2 (mag pot= 4pi*imag*CWORK1)
      INTEGER :: ISP, N1, N2, N3, NC, NI, I
      INTEGER :: NODE_ME, IONODE, IDUMP
      REAL(q) :: GX, GY, GZ, GSQU, GXYZ, SCALE

      PROFILING_START('remove_sources')
GPACC !$ACC ENTER DATA CREATE(CWORK_B,CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CWORK_B,CWORK1) IF(OFFLOAD_ON)
      !
      ! compute div(B); inspired by xcspin.F FEXCGS_NONCOL_ddsc_
      !
GPACC !$ACC KERNELS PRESENT(CWORK1) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CWORK1) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
      CWORK1 = REAL(0.0, KIND=q)
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

      ! loop over Bx, By, Bz
      DO ISP=2, WDES%NCDIJ
        ! set CWORK_B to total real potential in real space
GPACC !$ACC PARALLEL LOOP PRESENT(CWORK_B,CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK_B,CVTOT) IF(OFFLOAD_ON)
        DO I=1,GRIDC%MPLWV
          CWORK_B(I,ISP-1) = CVTOT(I,ISP)
        END DO
        ! set CWORK_B to total real potential in reciprocal space
        CALL FFT_RC_SCALE(CWORK_B(1,ISP-1),CWORK_B(1,ISP-1),GRIDC)
        ! set contribution of unbalanced lattic-vectors in potential to zero
        CALL SETUNB_COMPAT(CWORK_B,GRIDC)
        CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK_B(1,ISP-1))

        ! loop over local number of grid points of reciprocal grid
GPACC !$ACC PARALLEL LOOP PRESENT(CWORK_B,CWORK1,LATT_CUR,GRIDC) PRIVATE(N1,NC,N2,N3,GXYZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK_B,CWORK1,LATT_CUR,GRIDC) &
GPOMP !$OMP PRIVATE(N1,NC,N2,N3,GXYZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          ! Obtain indices for reciprocal vector
          N1= MOD((I-1),GRIDC%RC%NROW) +1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          ! Compute reciprocal vector GX, GY, GZ in Cartesian coordinates
          ! Note: LATT_CUR%B contains base of reciprocal lattice vectors without factor 2pi
          GXYZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(ISP-1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(ISP-1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(ISP-1,3))
          ! Compute div(B) = GX*Bx + GY*By + GZ*Bz
          ! Note: here we account for the factor TPI=2pi
          CWORK1(I) = CWORK1(I) + CWORK_B(I,ISP-1)*GXYZ*TPI
        END DO

      END DO

      ! inspired by POTHAR
      SCALE=1.0/TPI**2
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NI,GX,GY,GZ,GSQU) &
GPACC !$ACC PRESENT(GRIDC,CWORK1,LATT_CUR) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRIDC,CWORK1,LATT_CUR) &
GPOMP !$OMP PRIVATE(N2,N3,NI,GX,GY,GZ,GSQU) IF(OFFLOAD_ON)

      DO NC=1,GRIDC%RC%NCOL
__CPU   N2= GRIDC%RC%I2(NC)
__CPU   N3= GRIDC%RC%I3(NC)
        DO N1=1,GRIDC%RC%NROW
__GPU     N2= GRIDC%RC%I2(NC)
__GPU     N3= GRIDC%RC%I3(NC)
          ! compute G^2
          NI=(NC-1)*GRIDC%RC%NROW+N1

          GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
          GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
          GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
          GSQU=GX**2+GY**2+GZ**2

          ! compute CWORK1 = G*B(G)/G^2
          ! if G^2 not 0 divide by G^2
          IF ((GRIDC%LPCTX(N1)==0).AND.(GRIDC%LPCTY(N2)==0).AND.(GRIDC%LPCTZ(N3)==0)) THEN
            CWORK1(NI)=(0.0_q,0.0_q)
          ELSE
            CWORK1(NI)=CWORK1(NI)/GSQU*SCALE
          ENDIF
        END DO
      END DO

      ! loop over Bx, By, Bz
      DO ISP=2, WDES%NCDIJ

        ! loop over reciprocal grid
GPACC !$ACC PARALLEL LOOP PRESENT(CWORK_B) PRIVATE(N1,NC,N2,N3,GXYZ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK_B, CWORK1) &
GPOMP !$OMP PRIVATE(N1,NC,N2,N3,GXYZ) IF(OFFLOAD_ON)

        DO I=1,GRIDC%RC%NP
          ! Compute reciprocal vector GX, GY, GZ
          N1= MOD((I-1),GRIDC%RC%NROW) +1
          NC= (I-1)/GRIDC%RC%NROW+1
          N2= GRIDC%RC%I2(NC)
          N3= GRIDC%RC%I3(NC)
          GXYZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(ISP-1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(ISP-1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(ISP-1,3))

          ! B_sf(G) = B(G) - G * [G*B(G)]/G^2
          CWORK_B(I,ISP-1) = CWORK_B(I,ISP-1)  - TPI*GXYZ*CWORK1(I)
        END DO

        ! Set CWORK to source free potential in real space
        CALL FFT3D(CWORK_B(1,ISP-1),GRIDC,1)
        ! Update CVTOT with source free potential
GPACC !$ACC PARALLEL LOOP PRESENT(CWORK_B,CVTOT) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CWORK_B,CVTOT) IF(OFFLOAD_ON)
        DO I=1,GRIDC%MPLWV
           CVTOT(I,ISP) = CWORK_B(I,ISP-1)
        END DO
      END DO

GPACC !$ACC EXIT DATA DELETE(CWORK1,CWORK_B) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CWORK_B,CWORK1) IF(OFFLOAD_ON)
      PROFILING_STOP('remove_sources')

      END SUBROUTINE REMOVE_SOURCES

!************************ SUBROUTINE POT_FLIP **************************
!
!> Rearranges the storage mode for spin components of potentials.
!>
!> for the collinear case calculate: \n
!>  v0 1 + v_z \n
!>  v0 1 - v_z \n
!> for the non collinear case calculate \n
!>  v  = v0 1 + sigma_x v_x + simga_y v_y + sigma_z v_z
!
!***********************************************************************

      SUBROUTINE POT_FLIP(CVTOT, GRID, NCDIJ)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID
      COMPLEX(q) :: CVTOT(GRID%MPLWV, NCDIJ)

      ! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K

      IF (NCDIJ==2) THEN
         FAC=1.0_q
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT) PRIVATE(C00,CZ) IF(OFFLOAD_ON)
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CX,CY,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT) PRIVATE(C00,CX,CY,CZ) IF(OFFLOAD_ON)         
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP

    SUBROUTINE POT_FLIP_RL(CVTOT, GRID, NCDIJ)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID

      RGRID CVTOT(DIMREAL(GRID%MPLWV), NCDIJ)
      ! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K

      IF (NCDIJ==2) THEN
         FAC=1.0_q
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT) PRIVATE(C00,CZ) IF(OFFLOAD_ON)
         DO K=1,GRID%RL%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
GPACC !$ACC PARALLEL LOOP PRESENT(CVTOT) PRIVATE(C00,CX,CY,CZ) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CVTOT) PRIVATE(C00,CX,CY,CZ) IF(OFFLOAD_ON)
         DO K=1,GRID%RL%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP_RL


!************************ SUBROUTINE EXTERNAL_POT **********************
!
!> This subroutine can be used to add an external potential.
!>
!> The units of the potential are eV.
!
!***********************************************************************

      SUBROUTINE EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT)

__GPU USE moffload_struct_def

      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      RGRID      CVTOT(DIMREAL(GRIDC%MPLWV))

      RETURN

      IF (GRIDC%RL%NFAST==3) THEN
         ! mpi version: x-> N2, y-> N3, z-> N1
         N2MAX=GRIDC%NGX
         N3MAX=GRIDC%NGY
         N1MAX=GRIDC%NGZ

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NG) &
GPACC!$ACC& PRESENT(GRIDC,CVTOT) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRIDC,CVTOT) &
GPOMP !$OMP PRIVATE(N2,N3,NG) IF(OFFLOAD_ON)

         DO NC=1,GRIDC%RL%NCOL
__CPU       N2= GRIDC%RL%I2(NC)
__CPU       N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
__GPU          N2= GRIDC%RL%I2(NC)
__GPU          N3= GRIDC%RL%I3(NC)
               NG=(NC-1)*GRIDC%RL%NROW+N1
               CVTOT(NG)=CVTOT(NG)
               IF (N1<=GRIDC%NGZ/2) THEN
                  CVTOT(NG)=CVTOT(NG)+(REAL(N1,q)-1)/GRIDC%NGZ*2*0.1
               ELSE
                  CVTOT(NG)=CVTOT(NG)+(GRIDC%NGZ-REAL(N1,q))/GRIDC%NGZ*2*0.1
               ENDIF
            ENDDO
         ENDDO
      ELSE
         ! conventional version: x-> N1, y-> N2, z-> N3
         N1MAX=GRIDC%NGX
         N2MAX=GRIDC%NGY
         N3MAX=GRIDC%NGZ

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(N2,N3,NG) &
GPACC !$ACC& PRESENT(GRIDC,CVTOT) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:GRIDC,CVTOT) &
GPOMP !$OMP PRIVATE(N2,N3,NG) IF(OFFLOAD_ON)

         DO NC=1,GRIDC%RL%NCOL
__CPU       N2= GRIDC%RL%I2(NC)
__CPU       N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
__GPU          N2= GRIDC%RL%I2(NC)
__GPU          N3= GRIDC%RL%I3(NC)
               NG=(NC-1)*GRIDC%RL%NROW+N1
               CVTOT(NG)=CVTOT(NG)
               IF (N3<=GRIDC%NGZ/2) THEN
                  CVTOT(NG)=CVTOT(NG)+(REAL(N3,q)-1)/GRIDC%NGZ*2*0.1
               ELSE
                  CVTOT(NG)=CVTOT(NG)+(GRIDC%NGZ-REAL(N3,q))/GRIDC%NGZ*2*0.1
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE
