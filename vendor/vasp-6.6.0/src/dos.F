#include "symbol.inc"



!***********************SUBROUTINE DENNP *******************************
!
! if ISMEAR=0
! subroutine DENSTA calculates a continuous density of states in the
! interval (EMIN,EMAX) by applying a gaussian broadening to the discrete
! eigenvalue spectrum contained in CELEN(NBANDS,NKPTS). The width of the
! gaussians is SIGMA. The fermi energy EFERMI is calculated from the
! integrated dos
! correction to the variational energy is calculated (-SIGMA S)
! according to A.de Vita
! if ISMEAR>0 the generalized form of Methfessel and Paxton of order
!        N=ISMEAR will be used instead of Gaussians to get the dos ...
! routine is parallelized to get full speed ...
! initially it took 2 seconds to calculate occupancies (gK)
!
!***********************************************************************
  MODULE density_of_states
    USE prec
    IMPLICIT NONE
  CONTAINS
!********************DENINI**********************************************
! RCS:  $Id: dos.F,v 1.3 2002/08/14 13:59:37 kresse Exp $
!
! initialise the fermi-weights
! if the total number of bands is larger than the number of occupied
! bands:
! initialize the fermi occupation function
! FERWE  1  for the lowest fully occupied nelect/2 bands
!       .5  for the (nelect/2+1)th band if nelect odd
!        0  for the empty bands
!
! non-metallic systems set all components of FERWE equal 1
!
!************************************************************************

    SUBROUTINE  DENINI(FERWE,NBANDS,NKPTS,ELECT,LNONCOLLINEAR,NKSTART)
      USE prec
      IMPLICIT NONE
      LOGICAL LNONCOLLINEAR
      REAL(q) ELECT
      INTEGER NBANDS
      INTEGER NKPTS
      INTEGER, OPTIONAL :: NKSTART
      REAL(q) FERWE(:,:)
! local
      INTEGER NB, KP, NBE
      INTEGER MY_NKSTART
      REAL(q) :: SPINMULT
      MY_NKSTART = 1; IF (PRESENT(NKSTART)) MY_NKSTART = NKSTART
      IF (LNONCOLLINEAR) THEN
         SPINMULT=1
      ELSE
         SPINMULT=2
      ENDIF

      IF (ABS(NBANDS-ELECT/SPINMULT)>0.000001_q) THEN
         NBE=MIN(INT(ELECT/SPINMULT),NBANDS)
         DO NB=1,NBE
         DO KP=MY_NKSTART,NKPTS
           FERWE(NB,KP)=1._q
         ENDDO
         ENDDO

         IF (NBE+1<=NBANDS) THEN
            IF (ABS(ELECT/SPINMULT-NBE)>0.000001_q) THEN
               DO KP=MY_NKSTART,NKPTS
                  FERWE(NBE+1,KP)=0.5_q*(ELECT/SPINMULT-NBE)
               ENDDO
            ELSE
               DO KP=MY_NKSTART,NKPTS
                  FERWE(NBE+1,KP)=0
               ENDDO
            ENDIF
         ENDIF

         DO NB=NBE+2,NBANDS
         DO KP=MY_NKSTART,NKPTS
            FERWE(NB,KP)=0
         ENDDO
         ENDDO

      ELSE
         DO NB=1,NBANDS
         DO KP=MY_NKSTART,NKPTS
           FERWE(NB,KP)=1
         ENDDO
         ENDDO
      ENDIF
      RETURN
    END SUBROUTINE  DENINI


    SUBROUTINE DENMP(IU0, WDES, ISPIN, RSPIN, EMIN, EMAX, EFERMI_FORCE, &
         METHOD, NELECT, ECORR, EFERMI, ISMEAR, SIGMA, FERWE, CELEN , &
         NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR, JOBPAR)

      USE prec
      USE constant
      USE fermi_energy
      USE wave
      INTEGER IU0                     ! output unit
      TYPE (wavedes)        WDES
      INTEGER ISPIN                   ! number of spin channels
      REAL(q) RSPIN                   ! spin multiplicity
      REAL(q) EMIN, EMAX              ! minimum and maximum values for DOS
      REAL(q) EFERMI_FORCE            ! force Fermi-level to this value
      REAL(q) NELECT                  ! number of electrons
      REAL(q) ECORR                   ! entropy
      REAL(q) EFERMI                  ! Fermi-level
      INTEGER NEDOS                   ! number of slots in DOS
      INTEGER LDIMP                   ! number of lm quantum numbers considered in partial dos
      INTEGER NIOND                   ! number species considered
      CHARACTER(*) METHOD             ! method used to determine the Fermi energy
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)  ! DOS and integrated DOS
      REAL(q) PAR(WDES%NB_TOT,WDES%NKPTS,LDIMP,NIOND,ISPIN)
      REAL(q) DOSPAR(NEDOS,LDIMP,NIOND,ISPIN)
      INTEGER ISMEAR                  ! type of smearing
      REAL(q) SIGMA                   ! width of broadening
      REAL(q) FERWE(:,:,:)            ! occupancies of one electron states
      COMPLEX(q) CELEN(:,:,:)         ! eigenvalues
      INTEGER :: JOBPAR
! local variables
      TYPE(default_smearing) SMEARING
      LOGICAL LOWB,HIGHB
      INTEGER ISP
      INTEGER K, N, NELOW, NEHIG, I, NI, LP, NB_GLOBAL, NITER
      REAL(q) SIGMA_, DELTAE, EPS, WEIGHT, SFUN_DONE, E, DFUN, SFUN, EPSDOS, DOSTOT
      REAL(q) EF1, EF2, ELECT, X1

      SIGMA_=ABS(SIGMA)
      IF (SIGMA_==0) RETURN

      CALL SMEARING%INITIALIZE(REAL(CELEN), RSPIN, WDES)
      SMEARING%ISMEAR = ISMEAR
      SMEARING%SIGMA = SIGMA_

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
!=======================================================================
! initialize arrays for dos and integr. dos
!=======================================================================
      IF (JOBPAR==1) THEN
         DOSPAR=0
      ENDIF

      DOS =0
      DOSI=0
!=======================================================================
! accumulate dos and integrated dos
!=======================================================================
      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      DO N=1,WDES%NBANDS
        EPS=CELEN(N,K,ISP)
        WEIGHT= RSPIN*WDES%WTKPT(K)

        NELOW=(EPS-8._q*SIGMA_-EMIN)/DELTAE+1
        NEHIG=(EPS+8._q*SIGMA_-EMIN)/DELTAE+1
        IF (NELOW<1)     NELOW=1
        IF (NELOW>NEDOS) NELOW=NEDOS
        IF (NEHIG<1)     NEHIG=1
        IF (NEHIG>NEDOS) NEHIG=NEDOS

        SFUN_DONE=0
        DO I=NELOW,NEHIG
          E=EMIN+DELTAE*(I-1)-EPS
          CALL DELSTP(ISMEAR,(E/SIGMA_),DFUN,SFUN)
          EPSDOS=DFUN/SIGMA_
!gK fix the DOS so that the integrated DOS yields accurate results
          EPSDOS=(SFUN-SFUN_DONE)/DELTAE
          SFUN_DONE=SFUN

          DOS(I,ISP) =DOS(I,ISP) +(WEIGHT*EPSDOS)
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT*SFUN
          IF (JOBPAR==1) THEN
             DO NI=1,NIOND
             DO LP=1,LDIMP
#ifdef MPI
                NB_GLOBAL=(N-1)*WDES%COMM_INTER%NCPU+WDES%COMM_INTER%NODE_ME
#else
                NB_GLOBAL=N
#endif
                ! TODO UINIT: DOSPAR and PAR are at least partially used
                ! uninitialized here (see linear_optics.F:1543)!
                ! Are the loop limits too large (bulk_CoO_wurzite_SOC_G0W0_sym_RPR)?
                DOSPAR(I,LP,NI,ISP)=DOSPAR(I,LP,NI,ISP)+ &
                         (WEIGHT*EPSDOS)*PAR(NB_GLOBAL,K,LP,NI,ISP)
             ENDDO
             ENDDO
          ENDIF
        ENDDO
        DO I=NEHIG+1,NEDOS
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALLMPI (M_sum_d(WDES%COMM_INTER, DOS(1,1),  NEDOS*ISPIN))
      CALLMPI (M_sum_d(WDES%COMM_KINTER, DOS(1,1),  NEDOS*ISPIN))

      CALLMPI (M_sum_d(WDES%COMM_INTER, DOSI(1,1), NEDOS*ISPIN))
      CALLMPI (M_sum_d(WDES%COMM_KINTER, DOSI(1,1), NEDOS*ISPIN))
      IF (JOBPAR==1) THEN
        CALLMPI (M_sum_d(WDES%COMM_INTER, DOSPAR(1,1,1,1), NEDOS*LDIMP*NIOND*ISPIN))
        CALLMPI (M_sum_d(WDES%COMM_KINTER, DOSPAR(1,1,1,1), NEDOS*LDIMP*NIOND*ISPIN))
      ENDIF

      IF (ISMEAR == -2) THEN
         ! FIXME: this recovers the old behaviour: the Fermi energy is wrong though.
         EFERMI = CALCULATE_FERMI_ENERGY_LEGACY(SMEARING, NELECT, EMIN, DELTAE, &
            DOSI, -SIGMA_, EFERMI_FORCE, 'DENMP')
         ! Return here, otherwise the Fermi weights will be recomputed.
         RETURN
      ENDIF

      IF (METHOD == 'LEGACY' .OR. METHOD == 'FORCE') THEN
         EFERMI = CALCULATE_FERMI_ENERGY_LEGACY(SMEARING, NELECT, EMIN, DELTAE, &
            DOSI, SIGMA, EFERMI_FORCE, 'DENMP')
      ELSE
         EFERMI = CALCULATE_FERMI_ENERGY(SMEARING, NELECT)
      END IF
      CALL CALCULATE_FERMI_WEIGHTS(SMEARING, EFERMI, FERWE)

!=======================================================================
! calculate entropy -SIGMA_ * S
!=======================================================================
      ECORR=0
      DO ISP=1,ISPIN
         DO K=1,WDES%NKPTS
#ifdef MPI
            IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
            DO N=1,WDES%NBANDS
               X1=(EFERMI-CELEN(N,K,ISP))/SIGMA_
               CALL DELSTP(ISMEAR,X1,DFUN,SFUN)
               ECORR=ECORR+0.5_q*DFUN*WDES%WTKPT(K)
               IF (ISMEAR>0) THEN
                  CALL DELSTP((ISMEAR-1),X1,DFUN,SFUN)
                  ECORR=ECORR-0.5_q*DFUN*WDES%WTKPT(K)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ECORR=-ECORR*SIGMA_*RSPIN
      CALLMPI( M_sum_d( WDES%COMM_INTER, ECORR, 1))
      CALLMPI( M_sum_d( WDES%COMM_KINTER, ECORR, 1))

      RETURN
    END SUBROUTINE DENMP

!***********************SUBROUTINE DENFER******************************
!
! subroutine DENFER calculates a continuous density of states in the
! interval (EMIN,EMAX) by applying a fermi broadening to the discrete
! eigenvalue spectrum contained in CELEN(NBANDS,NKPTS). The width of the
! broadening is SIGMA_ (=1/BETA). EFERMI is calculated so that
! sum over occupation-numbers is equal to number of electrons
! in addition the correction to the variational energy is calculated
! as proposed by M. Weinert J.W. Davenport, Phys Rev B 45, 13709 (1992)
!
!**********************************************************************

     SUBROUTINE DENFER(IU0, WDES, ISPIN, RSPIN, EMIN, EMAX, EFERMI_FORCE, &
              METHOD, NELECT, ECORR, EFERMI, SIGMA, FERWE, CELEN , &
              NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR, JOBPAR)

      USE prec
      USE fermi_energy
      USE wave
      INTEGER IU0                     ! output unit
      TYPE (wavedes)        WDES
      INTEGER ISPIN                   ! number of spin channels
      REAL(q) RSPIN                   ! spin multiplicity
      REAL(q) EMIN, EMAX              ! minimum and maximum values for DOS
      REAL(q) EFERMI_FORCE            ! force Fermi-level to this value
      REAL(q) NELECT                  ! number of electrons
      REAL(q) ECORR                   ! entropy
      REAL(q) EFERMI                  ! Fermi-level
      INTEGER NEDOS                   ! number of slots in DOS
      INTEGER LDIMP                   ! number of lm quantum numbers considered in partial dos
      INTEGER NIOND                   ! number species considered
      CHARACTER(*) METHOD             ! method used to determine the Fermi energy
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)  ! DOS and integrated DOS
      REAL(q) PAR(WDES%NB_TOT,WDES%NKPTS,LDIMP,NIOND,ISPIN)
      REAL(q) DOSPAR(NEDOS,LDIMP,NIOND,ISPIN)
      INTEGER, PARAMETER :: ISMEAR = -1 ! type of smearing
      REAL(q) SIGMA                   ! width of broadening
      REAL(q) FERWE(:,:,:)            ! occupancies of one electron states
      COMPLEX(q) CELEN(:,:,:)         ! eigenvalues
      INTEGER :: JOBPAR
! local variables
      TYPE(default_smearing) SMEARING
      LOGICAL LOWB,HIGHB
      INTEGER ISP
      INTEGER K, N, NELOW, NEHIG, I, NI, LP, NB_GLOBAL, NITER
      REAL(q) SIGMA_, DELTAE, EPS, WEIGHT, SFUN_DONE, E, DFUN, SFUN, EPSDOS, DOSTOT
      REAL(q) EF1, EF2, ELECT, X1, XX

      PROFILING_START( 'denfer' )
! fermi function and its derivative

      SIGMA_=ABS(SIGMA)
      IF (SIGMA_==0) THEN
         PROFILING_STOP( 'denfer' )
         RETURN
      ENDIF

      CALL SMEARING%INITIALIZE(REAL(CELEN), RSPIN, WDES)
      SMEARING%ISMEAR = ISMEAR
      SMEARING%SIGMA = SIGMA_

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
!=======================================================================
! initialize arrays for dos and integr. dos
!=======================================================================
      IF (JOBPAR==1) THEN
         DOSPAR=0
      ENDIF

      DOS =0
      DOSI=0

!=======================================================================
! accumulate dos and integrated dos
!=======================================================================

      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      DO N=1,WDES%NBANDS
        EPS=CELEN(N,K,ISP)
        WEIGHT= RSPIN*WDES%WTKPT(K)

        NELOW=(EPS-8._q*SIGMA_-EMIN)/DELTAE+1
        NEHIG=(EPS+8._q*SIGMA_-EMIN)/DELTAE+1
        IF (NELOW<1)     NELOW=1
        IF (NELOW>NEDOS) NELOW=NEDOS
        IF (NEHIG<1)     NEHIG=1
        IF (NEHIG>NEDOS) NEHIG=NEDOS

        SFUN_DONE=0
        DO I=NELOW,NEHIG
          E=EMIN+DELTAE*(I-1)-EPS
          SFUN=F(-E,SIGMA_)
          EPSDOS=G(E,SIGMA_)
!gK fix the DOS so that the integrated DOS yields accurate results
          EPSDOS=(SFUN-SFUN_DONE)/DELTAE
          SFUN_DONE=SFUN

          DOS(I,ISP) =DOS(I,ISP) +(WEIGHT*EPSDOS)
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT*SFUN
          IF (JOBPAR==1) THEN
             DO NI=1,NIOND
             DO LP=1,LDIMP
#ifdef MPI
                NB_GLOBAL=(N-1)*WDES%COMM_INTER%NCPU+WDES%COMM_INTER%NODE_ME
#else
                NB_GLOBAL=N
#endif
                DOSPAR(I,LP,NI,ISP)=DOSPAR(I,LP,NI,ISP)+ &
                          (WEIGHT*EPSDOS)*PAR(NB_GLOBAL,K,LP,NI,ISP)
             ENDDO
             ENDDO
          ENDIF
        ENDDO
        DO I=NEHIG+1,NEDOS
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT
        ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALLMPI (M_sum_d(WDES%COMM_INTER, DOS(1,1), NEDOS*ISPIN))
      CALLMPI (M_sum_d(WDES%COMM_KINTER, DOS(1,1), NEDOS*ISPIN))

      CALLMPI (M_sum_d(WDES%COMM_INTER, DOSI(1,1), NEDOS*ISPIN))
      CALLMPI (M_sum_d(WDES%COMM_KINTER, DOSI(1,1), NEDOS*ISPIN))

      IF (JOBPAR==1) THEN
        CALLMPI (M_sum_d(WDES%COMM_INTER, DOSPAR(1,1,1,1), NEDOS*LDIMP*NIOND*ISPIN))
        CALLMPI (M_sum_d(WDES%COMM_KINTER, DOSPAR(1,1,1,1), NEDOS*LDIMP*NIOND*ISPIN))
      ENDIF

      IF (METHOD == 'LEGACY' .OR. METHOD == 'FORCE') THEN
         EFERMI = CALCULATE_FERMI_ENERGY_LEGACY(SMEARING, NELECT, EMIN, DELTAE, &
            DOSI, SIGMA_, EFERMI_FORCE, 'DENFER')
      ELSE
         EFERMI = CALCULATE_FERMI_ENERGY(SMEARING, NELECT)
      END IF
      CALL CALCULATE_FERMI_WEIGHTS(SMEARING, EFERMI, FERWE)

!=======================================================================
! calculate entropy
!=======================================================================
      ECORR= 0

      DO ISP=1,ISPIN
         DO K=1,WDES%NKPTS
#ifdef MPI
            IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
            DO N=1,WDES%NBANDS
               IF (FERWE(N,K,ISP) /=0 .AND. FERWE(N,K,ISP) /=1) &
                    &  ECORR=ECORR+FERWE(N,K,ISP)    *LOG(FERWE(N,K,ISP))*WDES%WTKPT(K)+ &
                    &            (1-FERWE(N,K,ISP))*LOG(1-FERWE(N,K,ISP))*WDES%WTKPT(K)
            ENDDO
         ENDDO
      ENDDO

      ECORR=ECORR*SIGMA_*RSPIN
      CALLMPI( M_sum_d( WDES%COMM_INTER, ECORR, 1))
      CALLMPI( M_sum_d( WDES%COMM_KINTER, ECORR, 1))

      PROFILING_STOP( 'denfer' )
      RETURN
    END SUBROUTINE DENFER




    FUNCTION F(E,SIG)
      REAL(q) F,E,SIG
DOESI IF (E / SIG > 709.78_q) THEN; F = 0.0_q; RETURN; END IF ! Breaks FHEG_333_L (gnu), overflow in exp possible.
      F=  1/(1 + EXP(E /SIG))
    END FUNCTION F

    FUNCTION G(E,SIG)
      REAL(q) G,E,SIG
DOESI IF (E / SIG < -709.78_q) THEN; G = 0.0_q; RETURN; END IF ! Breaks FHEG_333_L (gnu), overflow in exp possible.
      G=  EXP(-E/SIG)/(1 + EXP(-E/SIG))**2/SIG
    END FUNCTION G

!***********************SUBROUTINE DESTA_DER **************************
!
! subroutine to set W%AUXTOT to the derivative of the one-electron
! occupancies
!
! it is assumed that DENSTA has been called beforehand so that
! EFERMI is properly set
! for the tetrahedron metho W%AUXTOT is simply set to zero, implying
! that this routine does not support the tetrahedron method
!
! the routine might not yet support all the "crazy" special
! settings DENSTA takes care off, however, since two Fermi-levels
! are passed this should not be an issue
!
!**********************************************************************

    SUBROUTINE DENSTA_DER( W, KPOINTS, FERDER)

GPACC USE moffload

      USE prec
      USE wave
      USE mkpoints
      USE fermi_energy

      IMPLICIT NONE

      TYPE (wavespin)       W
      TYPE (kpoints_struct) KPOINTS
      REAL(q) :: FERDER(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN)
   ! local
      INTEGER :: ISP, NK, N
      REAL(q) :: X, DFUN, SFUN, WEIGHT
!      REAL(q) :: EFERMI_SHIFT            ! derivative with respect to Fermi-weight (unused)
!      REAL(q) :: WSUM                    ! density of states at Fermi-level (unused)
            ! these can potentially be used to shift the Fermi-level (see rot.F)
      FERDER=0

      IF (KPOINTS%ISMEAR < -1) RETURN

!      WSUM=0.0_q
!      EFERMI_SHIFT=0.0_q
      DO ISP=1, W%WDES%ISPIN
         DO NK=1, W%WDES%NKPTS
            WEIGHT =W%WDES%WTKPT(NK)*W%WDES%RSPIN
            DO N =1,W%WDES%NB_TOTK(NK,ISP)
               X=  (W%EFERMI(ISP)- REAL( W%CELTOT(N,NK,ISP) ,KIND=q) )/KPOINTS%SIGMA
               CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)
               FERDER(N, NK, ISP)    = -DFUN/KPOINTS%SIGMA
!               WSUM        = WSUM        +(DFUN/KPOINTS%SIGMA)*WEIGHT
!               EFERMI_SHIFT = EFERMI_SHIFT+(DFUN/KPOINTS%SIGMA)*(REAL(W%CELTOT(N,NK,ISP),q)-REAL(W%CELTOT(N,NK,ISP),q))*WEIGHT
            ENDDO
         ENDDO
      ENDDO

!      CALLMPI(M_sum_d(W%WDES%COMM_KINTER,EFERMI_SHIFT,1))
!      CALLMPI(M_sum_d(W%WDES%COMM_KINTER,WSUM,1))

      RETURN
    END SUBROUTINE DENSTA_DER

  END MODULE density_of_states


!***********************SUBROUTINE DENTET******************************
!
!> subroutine DENTET calculates a continuous density of states in the
!> interval (EMIN,EMAX) applying the tetrahedron method to the discrete
!> eigenvalue spectrum in CELEN(NBANDS,NKPTS). EFERMI is calculated so
!> that sum over occupation-numbers is equal to number of electrons
!> executed on all nodes
!**********************************************************************

    SUBROUTINE DENTET(EFERMI_FORCE,METHOD,CELEN,WTKPT,NBANDS,NKPTS,DOS,DOSI,NEDOS, &
         ISPIN,RSPIN,EMIN,EMAX,IDTET,NTET,VOLWGT,NELECT,EFERMI,FERWE,ECORR,JOB, &
         IU6,PAR,DOSPAR,NKDIM,LDIMP,NIOND,JOBPAR,SIGMA,USE_DELTAS,EFERMI_NEDOS,WDES)

      USE prec
      USE fermi_energy, ONLY: calculate_fermi_energy, calculate_fermi_energy_legacy, too_small_sigma
      USE string, ONLY: str
      USE tet, ONLY: tetrahedron_calculator, new_calculator, fermi_weights_state
      USE tutor, ONLY: vtutor
      USE wave, ONLY: wavedes
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: METHOD !< method to evaluate the Fermi energy
      REAL(q), INTENT(IN) :: EFERMI_FORCE    !< force the Fermi energy to a certain value
      TYPE(wavedes), INTENT(IN) :: WDES      !< parallelize over the communicators in this descriptor
      INTEGER IU6                     !< output unit
      INTEGER NBANDS                  !< number of bands
      INTEGER NKPTS                   !< number of k-points
      INTEGER NKDIM                   !< dimension of array for k-points
      INTEGER NTET                    !< number of tetrahedrons
      INTEGER ISPIN                   !< number of spin channels
      REAL(q) RSPIN                   !< spin multiplicity
      REAL(q) EMIN, EMAX              !< minimum and maximum values for DOS
      INTEGER IDET                    !< number of tetrahedrons
      REAL(q) VOLWGT                  !< volume of each tetrahedron
      REAL(q) EFERMI                  !< Fermi-level
      REAL(q) NELECT                  !< number of electrons
      REAL(q) ECORR                   !< entropy
      INTEGER JOB                     !< type of calculation
      INTEGER NEDOS                   !< number of slots in DOS
      INTEGER LDIMP                   !< number of lm quantum numbers considered in partial dos
      INTEGER NIOND                   !< number species considered
      REAL(q) SIGMA                   !< width of the smearing
      LOGICAL USE_DELTAS              !< use exact derivatives to determine weights
      INTEGER EFERMI_NEDOS            !< number of points to compute the fermi energy using gauss-legendre integration
      COMPLEX(q) CELEN(NBANDS,NKDIM,ISPIN)
      INTEGER JOBPAR
      REAL(q) WTKPT(NKPTS)
      REAL(q) FERWE(NBANDS,NKDIM,ISPIN)
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN)
      INTEGER IDTET(0:4,NTET)
      REAL(q) DOSPAR(NEDOS,LDIMP*NIOND,ISPIN)
      REAL(q) PAR(NBANDS,NKDIM,LDIMP*NIOND,ISPIN)
      REAL(q), POINTER :: PARVIEW(:,:,:,:), DOSVIEW(:,:,:)
      TYPE(tetrahedron_calculator) CALCULATOR
      TYPE(fermi_weights_state) STATE
      REAL(q) ENERGY_STEP, MINIMAL_SIGMA
      REAL(q), ALLOCATABLE :: ENERGY_GRID(:)
      INTEGER IDOS

      ENERGY_STEP = (EMAX - EMIN) / REAL(NEDOS - 1, KIND=q)
      ALLOCATE(ENERGY_GRID(NEDOS))
      ENERGY_GRID = [(EMIN + ENERGY_STEP * (IDOS - 1), IDOS = 1, NEDOS)]

      CALCULATOR = NEW_CALCULATOR(REAL(CELEN,KIND=q), RSPIN, VOLWGT * IDTET(0,:), &
         IDTET(1:4,:), WDES%COMM, NELECT, SIGMA, USE_DELTAS, EFERMI_NEDOS)

!=======================================================================
! calculate dos and integrated dos
!=======================================================================

      CALL CALCULATOR%DOS(ENERGY_GRID, DOS, DOSI)
      IF (JOBPAR==1) THEN
          CALL CALCULATOR%DOS_WITH_WEIGHT(ENERGY_GRID, PAR, DOSPAR)
      END IF

!=======================================================================
! calculate fermi energy and weights
!=======================================================================
      IF (METHOD == 'LEGACY' .OR. METHOD == 'FORCE') THEN
         ! The Fermi energy determination has an early exit when the smearing is
         ! too small. Because the tetrahedron method may potentially use 0
         ! smearing we need to increase that value over the threshold to avoid
         ! skipping the rest of the routine. Note that the smearing is not used
         ! for the calculation of the number of electrons, only the one passed
         ! into the calculator is.
         MINIMAL_SIGMA = MAX(SIGMA, TOO_SMALL_SIGMA + 1e-10_q)
         EFERMI = CALCULATE_FERMI_ENERGY_LEGACY(CALCULATOR, NELECT, &
            EMIN, ENERGY_STEP, DOSI, MINIMAL_SIGMA, EFERMI_FORCE, 'DENTET')
      ELSE
         EFERMI = CALCULATE_FERMI_ENERGY(CALCULATOR, NELECT)
      END IF
      STATE%FERMI_ENERGY = EFERMI
      STATE%USE_BLOECHL = JOB > 0
      ALLOCATE(STATE%WEIGHT_KPOINTS, SOURCE=WTKPT)
      STATE%UNIT_OUT = IU6
      IF (JOB==0) RETURN
      CALL CALCULATOR%FERMI_WEIGHTS(STATE, FERWE)

!=======================================================================
! How to calculate the correction term to the total Energy???? Is there
! a correction term at all ('no smearing', analytical interpolation and
! then integration with 'delta-function sampling'!!!) ?????????????????
! People say generally: NO! THERE IS NO ENTROPY! --> believe it or not!
!=======================================================================
      ECORR= 0

      END SUBROUTINE DENTET


!***********************SUBROUTINE DENSTA*******************************
! DENSTA is the dispatcher to the routines calculating the
! partial occupancies
! depending on ISMEAR the required algorithm is chosen
! ISMEAR=0     gaussian smearing (DENMP)
! ISMEAR>0     the generalized form of Methfessel and Paxton of order
!              N=ISMEAR will be used instead of Gaussians to get the dos ...
! ISMEAR=-1    fermi-smearing DENSTA calls DENFER
! ISMEAR=-4,-5 the tetrahedron method will be used
! ISMEAR=-14,-15 correspond to -4 and -5, but allow for a smearing
!
! The width of the smearing is given by SIGMA
! corrections to the total energy are calculated and returned in ENTROPY
! if these corrections are added to the total energy the free
! variational energy is obtained
!
! in addition the routine also calculates the DOS and the integrated
! and if required also the partial DOS (NIOND and LDIMP must be set)
!
! this is the worst spagethi code in VASP, but still it works
!
!***********************************************************************

    SUBROUTINE DENSTA( IU0, IU, WDES, W, KPOINTS, NELECT, &
         NUP_DOWN, ENTROPY, EFERMI, SIGMA, LNOAUTO,  &
         NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR)

__GPU USE moffload

      USE prec
      USE wave
      USE mkpoints
      USE density_of_states
      IMPLICIT NONE
      TYPE (wavedes)        WDES
      TYPE (wavespin)       W
      TYPE (kpoints_struct) KPOINTS
      INTEGER NEDOS, LDIMP, NIOND, IU, IU0, JOBPAR
      REAL(q) NELECT,NUP_DOWN,ENTROPY,SIGMA,EFERMI,ENTROPY_
      LOGICAL LNOAUTO   ! determines whether KPOINTS%EMAX and KPOINTS%EMIN may be overwritten
      REAL(q) DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
! PAR and DOSPAR must not be zero sized if LDIMP or NIOND is 0
      REAL(q) PAR(WDES%NB_TOT, WDES%NKPTS, MAX(LDIMP,1), MAX(NIOND,1), WDES%NCDIJ)
      REAL(q) DOSPAR(NEDOS,MAX(LDIMP,1),MAX(NIOND,1),WDES%NCDIJ)
! local variables
      REAL(q) EADD,EPS,NEL,NEL_SHIFT
      REAL(q) RSPIN, SIGMA_TET
      LOGICAL,SAVE  :: LAUTO=.FALSE.
      INTEGER,SAVE  :: ICALLS=0
      INTEGER  :: ISP,ISP1,ISP2,K,N,ISPIN,ISPIN_MAX
      REAl(q) TMPMAX(2)
      LOGICAL USE_DELTAS

      PROFILING_START('densta')

__GPU CALL OFFLOAD_ON_PUSH(.FALSE.)

      ENTROPY=0
      ENTROPY_=0

      SELECT CASE(KPOINTS%ISMEAR)
      CASE (-15:-14)
         SIGMA_TET = SIGMA
         USE_DELTAS = .TRUE.
      CASE DEFAULT
         SIGMA_TET = dzero
         USE_DELTAS = .FALSE.
      END SELECT

      IF (NIOND == 0 .OR. LDIMP==0) THEN
        JOBPAR=0
      ELSE
        JOBPAR=1
      ENDIF

! set KPOINTS%EMAX and KPOINTS%EMIN if required
      ! first call initialize LAUTO
      IF (ICALLS==0.OR.(KPOINTS%EMAX<=KPOINTS%EMIN)) THEN
         LAUTO=(KPOINTS%EMAX<=KPOINTS%EMIN)
         ICALLS=1
      ENDIF

      IF (LAUTO .AND. .NOT. LNOAUTO) THEN
        KPOINTS%EMIN=1.E30_q
        KPOINTS%EMAX=-KPOINTS%EMIN
        DO ISP=1,WDES%ISPIN
        DO K=1,WDES%NKPTS
#ifdef MPI
        IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
        DO N=1,WDES%NB_TOTK(K,ISP)
          EPS=W%CELTOT(N,K,ISP)
          KPOINTS%EMAX=MAX(KPOINTS%EMAX,EPS)
          KPOINTS%EMIN=MIN(KPOINTS%EMIN,EPS)
        ENDDO
        ENDDO
        ENDDO

#ifdef MPI
        IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
           TMPMAX(1)=KPOINTS%EMAX
           TMPMAX(2)=-KPOINTS%EMIN
           CALLMPI ( M_max_d(WDES%COMM_KINTER,TMPMAX,2) )
           KPOINTS%EMAX=TMPMAX(1)
           KPOINTS%EMIN=-TMPMAX(2)
        ENDIF
#endif
        EADD=(KPOINTS%EMAX-KPOINTS%EMIN)*0.05_q
        EADD=MAX(EADD,10 *ABS(SIGMA))
        KPOINTS%EMIN=KPOINTS%EMIN-EADD
        KPOINTS%EMAX=KPOINTS%EMAX+EADD
      ENDIF
#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         IF (KPOINTS%ISMEAR==-4.OR.KPOINTS%ISMEAR==-5.OR.KPOINTS%LTET) THEN
!PK To avoid modifying tet.F, which is not mpi aware, synchronize eigenvalues here
            CALL KPAR_SYNC_CELTOT(WDES,W)
         END IF
      END IF
#endif
      ISPIN_MAX=1

      ISPIN=WDES%ISPIN
      RSPIN=WDES%RSPIN

      NEL_SHIFT=0
    ! calculation for a specific number of electrons in up and down component
    ! this actually requires me to do the worst fiddling ;(
      IF (NUP_DOWN >= 0 .AND. WDES%ISPIN>1 ) THEN
         ISPIN_MAX=2
         ISPIN    =1
         RSPIN    =2
         NEL_SHIFT=NUP_DOWN
    !
    ! non collinear calculation and partial DOS required
    ! also fiddle a little bit (call the DOS routines 4 times each time
    ! with a different PAR and DOSPAR)
      ELSE IF (WDES%LNONCOLLINEAR .AND. JOBPAR==1) THEN
         ISPIN_MAX=4
      ENDIF
    !
    ! now set the occupancies and calculate the DOS
    !
    spin: DO ISP=1,ISPIN_MAX
      ISP1=ISP
      IF (WDES%LNONCOLLINEAR) ISP1=1
      ISP2=ISP1+ISPIN-1

      IF (ISP==1) THEN
         NEL=NELECT+NEL_SHIFT
      ELSE
         NEL=NELECT-NEL_SHIFT
      ENDIF

      IF (KPOINTS%ISMEAR==-1) THEN
        CALL DENFER(IU0, WDES, ISPIN, RSPIN, KPOINTS%EMIN,KPOINTS%EMAX, &
           KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           NEL, ENTROPY_, EFERMI, SIGMA, W%FERWE(:,:,ISP1:ISP2), W%CELEN(:,:,ISP1:ISP2),&
           NEDOS, LDIMP, NIOND, DOS(1,ISP1), DOSI(1,ISP1), &
           PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP), JOBPAR)
        CALL MRG_FERWE(WDES,W)
      ELSE IF ((KPOINTS%ISMEAR==-4).OR.(KPOINTS%ISMEAR==-14)) THEN
        CALL DENTET(KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET, &
           KPOINTS%VOLWGT,NEL,EFERMI,W%FERTOT(1,1,ISP1), &
           ENTROPY_,-2,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKPTS,LDIMP,NIOND, &
           JOBPAR,SIGMA_TET,USE_DELTAS,KPOINTS%EFERMI_NEDOS,WDES)
      ELSE IF ((KPOINTS%ISMEAR==-5).OR.(KPOINTS%ISMEAR==-15)) THEN
        CALL DENTET(KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET,  &
           KPOINTS%VOLWGT,NEL,EFERMI,W%FERTOT(1,1,ISP1), &
           ENTROPY_,2,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKPTS,LDIMP,NIOND, &
           JOBPAR,SIGMA_TET,USE_DELTAS,KPOINTS%EFERMI_NEDOS,WDES)
      ELSE
        CALL DENMP(IU0, WDES, ISPIN, RSPIN, KPOINTS%EMIN,KPOINTS%EMAX, &
           KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           NEL, ENTROPY_, EFERMI, KPOINTS%ISMEAR, SIGMA, &
           W%FERWE(:,:,ISP1:ISP2) , W%CELEN(:,:,ISP1:ISP2),&
           NEDOS, LDIMP, NIOND,  DOS(1,ISP1), DOSI(1,ISP1), &
           PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP), JOBPAR)

        CALL MRG_FERWE(WDES,W)
      ENDIF
    ! for ISMEAR>=30, Methfessel-Paxton smearing of order ISMEAR-30 is used to
    ! calculate the  Fermi weights, but the tetrahedron method is used to get
    ! the density of states (one of the secret flag settings in VASP ...)
      ENTROPY=ENTROPY+ENTROPY_

      IF (KPOINTS%LTET .AND. KPOINTS%ISMEAR /=-4 .AND. KPOINTS%ISMEAR /= -5 .AND.&
          KPOINTS%ISMEAR/=-14 .AND. KPOINTS%ISMEAR/=-15) THEN
        CALL DENTET(KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET, &
           KPOINTS%VOLWGT,NEL,EFERMI,W%FERTOT(1,1,ISP1:ISP2), &
           ENTROPY_,0,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKPTS,LDIMP,NIOND, &
           JOBPAR,SIGMA_TET,USE_DELTAS,KPOINTS%EFERMI_NEDOS,WDES)
      ENDIF
      ! store fermi energy also in wavedes, used in chi and chi_base
      W%EFERMI(ISP1) = EFERMI
      ENDDO spin
      IF (NUP_DOWN < 0 .AND. WDES%ISPIN>1 ) THEN
         W%EFERMI(2) = W%EFERMI(1)
      ENDIF

    ! due to all this fiddling the entropy term is not right
    ! we need to correct this now
      IF (NUP_DOWN >= 0 .AND. WDES%ISPIN>1 ) THEN
         ENTROPY=ENTROPY/2
      ELSE IF (WDES%LNONCOLLINEAR .AND. JOBPAR==1) THEN
         ENTROPY=ENTROPY/4
      ENDIF

__GPU CALL OFFLOAD_ON_POP
GPACC !$ACC UPDATE DEVICE(W%FERTOT) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W,1))
GPOMP !$OMP TARGET UPDATE TO(W%FERTOT) IF(OFFLOAD_ON.AND.IS_ON_DEVICE(W))

      PROFILING_STOP('densta')

      RETURN
    END SUBROUTINE DENSTA


!***********************SUBROUTINE SYMPDOS *****************************
!
! The subroutine SYMPDOS symmetrizes the l(m)+spin+site resolved
! partial density of states
!
!***********************************************************************

    SUBROUTINE SYMPDOS(WDES,T_INFO,SYMM,LATT_CUR,LDIMP,LMDIMP,NEDOS,LMDIM,NIOND,DOSPAR)
      USE prec
      USE base
      USE wave
      USE poscar
      USE lattice

      TYPE(wavedes)   :: WDES
      TYPE(type_info) :: T_INFO
      TYPE(symmetry)  :: SYMM
      TYPE(latt)      :: LATT_CUR

      INTEGER :: NEDOS,LMDIM,NIOND,LDIMP,LMDIMP

      REAL(q) :: DOSPAR(NEDOS,LMDIM,NIOND,WDES%NCDIJ)

      REAL(q), ALLOCATABLE :: VEC(:,:,:)
      INTEGER :: I,ISP

      ! early exit
      IF ((LMDIM/=LDIMP.AND.LMDIM/=LMDIMP).OR.(SYMM%ISYM<0)) RETURN

      ALLOCATE(VEC(LMDIM,T_INFO%NIONS,WDES%NCDIJ))

      DO I=1,NEDOS
#ifdef MPI
         IF (MOD(I-1,WDES%COMM%NCPU).NE.WDES%COMM%NODE_ME-1) THEN
            DOSPAR(I,:,:,:)=0 ; CYCLE
         ENDIF
#endif
         ! copy a slot of DOSPAR into VEC
         VEC(:,:,:)=DOSPAR(I,1:LMDIM,1:T_INFO%NIONS,1:WDES%NCDIJ)
         ! bring VEC to density/magnetization storage mode (for ISPIN=2 only)
         CALL R_FLIP(VEC,LMDIM*T_INFO%NIONS,LMDIM*T_INFO%NIONS,WDES%NCDIJ,.FALSE.)

         IF (LMDIM==LDIMP) THEN
            DO ISP=1,WDES%ISPIN
               CALL PARSYL( LDIMP-1,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                           VEC(1,1,ISP),SYMM%ROTMAP(1,1,1),SYMM%MAGROT(1,1),ISP)
            ENDDO
            IF (WDES%LNONCOLLINEAR) &
               CALL MAGSYL( LDIMP-1,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                           VEC(1,1,2),SYMM%ROTMAP(1,1,1),SYMM%MAGROT(1,1),WDES%SAXIS,LATT_CUR%A,LATT_CUR%B)
         ELSE
            DO ISP=1,WDES%ISPIN
               CALL PARSYLM(LDIMP-1,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                           VEC(1,1,ISP),SYMM%ROTMAP(1,1,1),SYMM%MAGROT(1,1),LATT_CUR%A,LATT_CUR%B,ISP)
            ENDDO
            IF (WDES%LNONCOLLINEAR) &
               CALL MAGSYLM(LDIMP-1,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                           VEC(1,1,2),SYMM%ROTMAP(1,1,1),SYMM%MAGROT(1,1),WDES%SAXIS,LATT_CUR%A,LATT_CUR%B)
         ENDIF

         ! back to spin-up/down storage (for ISPIN=2 only)
         CALL R_FLIP(VEC,LMDIM*T_INFO%NIONS,LMDIM*T_INFO%NIONS,WDES%NCDIJ,.TRUE.)
         ! store back
         DOSPAR(I,1:LMDIM,1:T_INFO%NIONS,1:WDES%NCDIJ)=VEC(:,:,:)

      ENDDO

      CALLMPI(M_sum_d(WDES%COMM,DOSPAR(1,1,1,1),NEDOS*LMDIM*T_INFO%NIONS*WDES%NCDIJ))

      DEALLOCATE(VEC)

    END SUBROUTINE SYMPDOS


!***********************SUBROUTINE DENSTA_SPIN *************************
!
! this version returns two Fermi-levels, one for spin up and
! one for spin down in EFERMI(ISPIN)
!
!***********************************************************************

    SUBROUTINE DENSTA_SPIN( IU0, IU, WDES, W, KPOINTS, NELECT, &
         NUP_DOWN, ENTROPY, EFERMI, SIGMA, LNOAUTO,  &
         NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR)

      USE prec
      USE wave
      USE mkpoints
      USE density_of_states

      TYPE (wavedes)        WDES
      TYPE (wavespin)       W
      TYPE (kpoints_struct) KPOINTS
      INTEGER NEDOS, LDIMP, NIOND, IU, IU0, JOBPAR
      REAL(q) EFERMI(WDES%ISPIN)
      REAL(q) NELECT,NUP_DOWN,ENTROPY,SIGMA,ENTROPY_
      LOGICAL LNOAUTO   ! determines whether KPOINTS%EMAX and KPOINTS%EMIN may be overwritten
      REAL(q) DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
      REAL(q) PAR(WDES%NB_TOT, WDES%NKPTS, MAX(LDIMP,1), MAX(NIOND,1), WDES%NCDIJ)
      REAL(q) DOSPAR(NEDOS,MAX(LDIMP,1),MAX(NIOND,1),WDES%NCDIJ)
! local variables
      REAL(q) EADD,EPS,NEL,NEL_SHIFT
      REAL(q) RSPIN, SIGMA_TET
      LOGICAL,SAVE  :: LAUTO=.FALSE.
      INTEGER,SAVE  :: ICALLS=0
      INTEGER  :: ISP,ISP2,K,N,ISPIN,ISPIN_MAX
      REAl(q) TMPMAX(2)
      LOGICAL USE_DELTAS

      ENTROPY=0
      ENTROPY_=0

      SELECT CASE(KPOINTS%ISMEAR)
      CASE (-15:-14)
         SIGMA_TET = SIGMA
         USE_DELTAS = .TRUE.
      CASE DEFAULT
         SIGMA_TET = dzero
         USE_DELTAS = .FALSE.
      END SELECT

      IF (NIOND == 0 .OR. LDIMP==0) THEN
        JOBPAR=0
      ELSE
        JOBPAR=1
      ENDIF

! set KPOINTS%EMAX and KPOINTS%EMIN if required
      ! first call initialize LAUTO
      IF (ICALLS==0.OR.(KPOINTS%EMAX<=KPOINTS%EMIN)) THEN
         LAUTO=(KPOINTS%EMAX<=KPOINTS%EMIN)
         ICALLS=1
      ENDIF

      IF (LAUTO .AND. .NOT. LNOAUTO) THEN
        KPOINTS%EMIN=1.E30_q
        KPOINTS%EMAX=-KPOINTS%EMIN
        DO ISP=1,WDES%ISPIN
        DO K=1,WDES%NKPTS
#ifdef MPI
        IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
        DO N=1,WDES%NB_TOTK(K,ISP)
          EPS=W%CELTOT(N,K,ISP)
          KPOINTS%EMAX=MAX(KPOINTS%EMAX,EPS)
          KPOINTS%EMIN=MIN(KPOINTS%EMIN,EPS)
        ENDDO
        ENDDO
        ENDDO

#ifdef MPI
        IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
           TMPMAX(1)=KPOINTS%EMAX
           TMPMAX(2)=-KPOINTS%EMIN
           CALLMPI ( M_max_d(WDES%COMM_KINTER,TMPMAX,2) )
           KPOINTS%EMAX=TMPMAX(1)
           KPOINTS%EMIN=-TMPMAX(2)
        ENDIF
#endif
        EADD=(KPOINTS%EMAX-KPOINTS%EMIN)*0.05_q
        EADD=MAX(EADD,10 *ABS(SIGMA))
        KPOINTS%EMIN=KPOINTS%EMIN-EADD
        KPOINTS%EMAX=KPOINTS%EMAX+EADD
      ENDIF
#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         IF (KPOINTS%ISMEAR==-4.OR.KPOINTS%ISMEAR==-5.OR.KPOINTS%LTET) THEN
!PK To avoid modifying tet.F, which is not mpi aware, synchronize eigenvalues here
            CALL KPAR_SYNC_CELTOT(WDES,W)
         END IF
      END IF
#endif
      ISPIN_MAX=1

      ISPIN=WDES%ISPIN
      RSPIN=WDES%RSPIN

      NEL_SHIFT=0
    ! calculation for a specific number of electrons in up and down component
    ! this actually requires me to do the worst fiddling ;(
      IF (NUP_DOWN >= 0 .AND. WDES%ISPIN>1 ) THEN
         ISPIN_MAX=2
         ISPIN    =1
         RSPIN    =2
         NEL_SHIFT=NUP_DOWN
    !
    ! non collinear calculation and partial DOS required
    ! also fiddle a little bit (call the DOS routines 4 time each time
    ! with a different PAR and DOSPAR)
      ELSE IF (WDES%LNONCOLLINEAR .AND. JOBPAR==1) THEN
         ISPIN_MAX=4
      ENDIF
    !
    ! now set the occupancies and calculate the DOS
    !
    spin: DO ISP=1,ISPIN_MAX
      ISP1=ISP
      IF (WDES%LNONCOLLINEAR) ISP1=1
      ISP2=ISP1+ISPIN-1

      IF (ISP==1) THEN
         NEL=NELECT+NEL_SHIFT
      ELSE
         NEL=NELECT-NEL_SHIFT
      ENDIF

      IF (KPOINTS%ISMEAR==-1) THEN
        CALL DENFER(IU0, WDES, ISPIN, RSPIN, KPOINTS%EMIN,KPOINTS%EMAX, &
           KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           NEL, ENTROPY_, EFERMI(ISP), SIGMA, W%FERWE(:,:,ISP1:ISP2), W%CELEN(:,:,ISP1:ISP2),&
           NEDOS, LDIMP, NIOND, DOS(1,ISP1), DOSI(1,ISP1), &
           PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP), JOBPAR)
        CALL MRG_FERWE(WDES,W)
      ELSE IF ((KPOINTS%ISMEAR==-4).OR.(KPOINTS%ISMEAR==-14)) THEN
        CALL DENTET(KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET, &
           KPOINTS%VOLWGT,NEL,EFERMI(ISP),W%FERTOT(1,1,ISP1), &
           ENTROPY_,-2,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKPTS,LDIMP,NIOND, &
           JOBPAR,SIGMA_TET,USE_DELTAS,KPOINTS%EFERMI_NEDOS,WDES)
      ELSE IF ((KPOINTS%ISMEAR==-5).OR.(KPOINTS%ISMEAR==-15)) THEN
        CALL DENTET(KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET,  &
           KPOINTS%VOLWGT,NEL,EFERMI(ISP),W%FERTOT(1,1,ISP1), &
           ENTROPY_,2,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKPTS,LDIMP,NIOND, &
           JOBPAR,SIGMA_TET,USE_DELTAS,KPOINTS%EFERMI_NEDOS,WDES)
      ELSE
        CALL DENMP(IU0, WDES, ISPIN, RSPIN, KPOINTS%EMIN,KPOINTS%EMAX, &
           KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           NEL, ENTROPY_, EFERMI(ISP), KPOINTS%ISMEAR, SIGMA, &
           W%FERWE(:,:,ISP1:ISP2) , W%CELEN(:,:,ISP1:ISP2),&
           NEDOS, LDIMP, NIOND,  DOS(1,ISP1), DOSI(1,ISP1), &
           PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP), JOBPAR)

        CALL MRG_FERWE(WDES,W)
      ENDIF
    ! for ISMEAR>=30, Methfessel-Paxton smearing of order ISMEAR-30 is used to
    ! calculate the  Fermi weights, but the tetrahedron method is used to get
    ! the density of states (one of the secret flag settings in VASP ...)
      ENTROPY=ENTROPY+ENTROPY_

      IF (KPOINTS%LTET .AND. KPOINTS%ISMEAR /=-4 .AND. KPOINTS%ISMEAR /= -5 .AND. &
          KPOINTS%ISMEAR /= -14 .AND. KPOINTS%ISMEAR /= -15) THEN
        CALL DENTET(KPOINTS%EFERMI, KPOINTS%FERMI_METHOD, &
           W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET, &
           KPOINTS%VOLWGT,NEL,EFERMI(ISP),W%FERTOT(1,1,ISP1:ISP2), &
           ENTROPY_,0,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKPTS,LDIMP,NIOND, &
           JOBPAR,SIGMA_TET,USE_DELTAS,KPOINTS%EFERMI_NEDOS,WDES)
      ENDIF
      ! store fermi energy also in wavespin, used in chi and chi_base
      W%EFERMI(ISP1) = EFERMI(ISP)
      W%EFERMI(ISP2) = EFERMI(ISP)

      ENDDO spin

    ! due to all this fiddling the entropy term is not right
    ! we need to correct this now
      IF (NUP_DOWN >= 0 .AND. WDES%ISPIN>1 ) THEN
         ENTROPY=ENTROPY/2
      ELSE IF (WDES%LNONCOLLINEAR .AND. JOBPAR==1) THEN
         ENTROPY=ENTROPY/4
      ENDIF

    END SUBROUTINE DENSTA_SPIN



!***********************SUBROUTINE DENSTA_SPIN *************************
!
! DENSTA has a large amount of totally unnecessary stuff
! this is a simpler interface that allocates all required arrays
! temporarily and then calls DENSTA
! this version can be used if only the one-electron occupancies
! need to be recalculated
!
!***********************************************************************


    SUBROUTINE DENSTA_SIMPLE(W, KPOINTS, NELECT, NUP_DOWN )
      USE prec
      USE wave
      USE mkpoints

      IMPLICIT NONE
      TYPE (wavespin)       W
      TYPE (kpoints_struct) KPOINTS
      REAL(q) NELECT, NUP_DOWN

    ! few local variables required to call DENSTA
      REAL(q) :: EENTROPY, EFERMI
      REAL(q) ::  PAR(1,1,1,1,W%WDES%NCDIJ),DOSPAR(1,1,1,W%WDES%NCDIJ)
      INTEGER, PARAMETER :: NEDOS=100
      REAL(q) ::  DOS(NEDOS,W%WDES%ISPIN),DOSI(NEDOS,W%WDES%ISPIN)

      CALL DENSTA( -1, -1, W%WDES, W, KPOINTS, NELECT, &
           NUP_DOWN, EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
           NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

    END SUBROUTINE DENSTA_SIMPLE

