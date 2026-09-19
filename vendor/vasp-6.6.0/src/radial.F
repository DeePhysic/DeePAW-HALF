#include "symbol.inc"
!#define vector
!*******************************************************************
! RCS:  $Id: radial.F,v 1.11 2003/06/27 13:22:22 kresse Exp kresse $
!
!>  MODULE which supports operations on radial grid. 
!>  All routines written by gK with the exception of the
!>  routines required for metaGGA, which were written by
!>  Robin Hirsch, Dec 2000.
!
!*******************************************************************
  MODULE radial
    USE prec
    USE radial_struct_def

    ! auxiliary functions needed by RAD_KINEDEN and RAD_PROJ_METAGGA,
    ! allocated and constructed by RAD_AUXILIARY_FUNCTIONS_METAGGA
    REAL(q), ALLOCATABLE, PRIVATE, SAVE :: X_YLM_X_YLM(:,:,:)
    REAL(q), ALLOCATABLE, PRIVATE, SAVE :: YLMD_YLMD(:,:,:)
    REAL(q), ALLOCATABLE, PRIVATE, SAVE :: X_YLM_YLMD(:,:,:)

#ifndef noPAW
  CONTAINS

!*******************************************************************
!
! RAD_ALIGN
!
!> Aligns R%RMAX to the radial grid.
!>
!> This should be used to avoid that R%RMAX is slightly off
!> from one of the grid points (due to rounding errors on reading).
!
!*******************************************************************

    SUBROUTINE RAD_ALIGN(R,RMAX)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (rgrid) R
      REAL(q) RMAX

      RMAX=RMAX-1E-5_q

      DO N=1,R%NMAX
         IF (R%R(N).GT. RMAX) THEN
            RMAX=R%R(N)
            EXIT
         ENDIF
      ENDDO

      IF (ASSOCIATED(R%SI)) DEALLOCATE(R%SI)
      CALL SET_SIMP(R)
    END SUBROUTINE RAD_ALIGN

!*******************************************************************
!
!  RAD_CHECK_QPAW
!
!>  Checks the consistency of the array QPAW
!>  with the all electron and pseustored wavefunctions WAE and WPS
!>    Q(PAW,ll') = \int (phi^AE(l,n,r)^2  -  phi^PS(l',n,r)^2) dr
!>
!>  If there are large errors the program stops and reports
!>  the error.
!>
!>  In all cases QPAW is corrected, so that it is exactly equal to
!>  the integral.
!>
!>  Remember the wavefunction psi(r) can be obtained from the stored
!>  coefficients using:
!>      psi(r) = \sum_lmn Y_lm(r) w_ln(r) / r
!
!*******************************************************************

    SUBROUTINE RAD_CHECK_QPAW( R, CHANNELS, WAE, WPS , QPAW, QTOT, L, RDEP )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (rgrid) R
      INTEGER CHANNELS
      REAL(q) :: WAE(:,:),WPS(:,:)   !< AE and soft wavefunctions
      REAL(q) :: QPAW(:,:,0:)        !< moments of compensation charge
      REAL(q) :: QTOT(:,:)           !< L=0 moments of AE charge
      INTEGER :: L(:)
      REAL(q) :: RDEP                !< outermost radius of augmentation charge 
! local variables
      REAL(q) :: RHOT(R%NMAX)
      INTEGER CH1,CH2,I
      REAL(q) :: RES
      INTEGER :: LL,LLP,LMIN,LMAX,LMAIN
      INTEGER :: IRMAX,NMAX_STORE

      DO IRMAX=1,R%NMAX-1
         ! because of rounding errors we need to be somewhat sloppy for the comparison
         IF (RDEP>0 .AND. R%R(IRMAX)-RDEP > -5E-3) EXIT
      ENDDO

      ! set the simpson weights in accordance with IRMAX
      NMAX_STORE=R%NMAX; R%NMAX=IRMAX
!      WRITE(*,*) NMAX_STORE, IRMAX
      IF (ASSOCIATED(R%SI)) DEALLOCATE(R%SI)
      CALL SET_SIMP(R)

      DO CH1=1,CHANNELS
      DO CH2=1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
         IF (LL==LLP) THEN
            RHOT=0
            DO I=1,IRMAX
               RHOT(I)=WAE(I,CH1)*WAE(I,CH2)
            ENDDO
            CALL SIMPI(R, RHOT , RES)
            QTOT(CH1,CH2)=RES
         ENDIF
      ENDDO
      ENDDO

      ! reset the simpson weights
      DEALLOCATE(R%SI)
      R%NMAX=NMAX_STORE
      CALL SET_SIMP(R)

      DO CH1=1,CHANNELS
      DO CH2=1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         DO LMAIN=LMIN,LMAX,2
            DO I=1,R%NMAX
               RHOT(I)=(WAE(I,CH1)*WAE(I,CH2)-WPS(I,CH1)*WPS(I,CH2))*R%R(I)**LMAIN
            ENDDO
            CALL SIMPI(R, RHOT , RES)
         ! screwed if we do not have the correct integrated charge
            IF ( ABS(RES-QPAW(CH1,CH2,0)) > 1E-4 .AND. LMAIN==0 ) THEN
               CALL vtutor%bug("internal error RAD_CHECK_QPAW: QPAW is incorrect channels " // &
                  str(CH1) // " " // str(CH2) // " QPAW= " // str(RES) // " int= " // &
                  str(QPAW(CH1,CH2,0)) // " " // str(RHOT(R%NMAX-4:R%NMAX)), __FILE__, __LINE__)
            ENDIF
            QPAW(CH1,CH2,LMAIN)=RES
         ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE RAD_CHECK_QPAW


!*******************************************************************
!
!  RAD_AUG_CHARGE
!
!>  Add the compensation charge density on the radial grid
!>  to RHO
!>
!>  RHO(L,M,r) = RHO(L,M,r) +  Q(r,L) Q(PAW,ll' L) RHOLM(ll',LM)
!>
!>  Q(r,L) are the L-dependent 1-normalized compensation charges.
!
!*******************************************************************

    SUBROUTINE RAD_AUG_CHARGE( RHO, R, RHOLM, CHANNELS, L, &
         LYMAX, AUG, QPAW )
      USE constant
      IMPLICIT NONE

      REAL(q) :: RHO(:,:)    !< charge on radial grid
      TYPE (rgrid) R
      REAL(q) :: RHOLM(:)    !< occupancy of each llLM channel
      REAL(q) :: AUG(:,0:)   !< 1-normalized L-dep compensation charge
      REAL(q) :: QPAW(:,:,0:)!< moments of compensation charge Q(PAW,ll L)
      INTEGER :: LYMAX       !< maximum L
      INTEGER CHANNELS, L(:)
   ! local variables
      REAL(q) :: RHOLMT((LYMAX+1)*(LYMAX+1))
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I,LMAX
      INTEGER IBASE,JBASE,LMIN,LMAIN,MMAIN,LMMAIN

      PROFILING_START('rad_aug_charge')

!-----------------------------------------------------------------------
! first contract to L and M dependent Q(LM)
!-----------------------------------------------------------------------
     IBASE=0
      RHOLMT=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependent charge
         DO LMAIN=LMIN,MIN(LMAX,LYMAX),2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            RHOLMT(LMMAIN)=RHOLMT(LMMAIN)+QPAW(CH1,CH2,LMAIN)*RHOLM(LMMAIN+JBASE)
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
!-----------------------------------------------------------------------
! then add to charge on radial grid
!-----------------------------------------------------------------------
      ! could be replaced by matrix vecor DGEMV
      DO LMAIN=0,LYMAX
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DO I=1,R%NMAX
               RHO(I,LMMAIN)=RHO(I,LMMAIN)+AUG(I,LMAIN)*RHOLMT(LMMAIN)
            ENDDO
         ENDDO
      ENDDO
#ifdef debug
      WRITE(0,*) 'RAD_AUG_CHARGE: L compensation charges are'
      DO LL=0,LYMAX
         WRITE(0,'(I3,10F10.5)') LL,(RHOLMT(LL**2+I)*2*SQRT(PI),I=1,LL*2+1)
      ENDDO
#endif

      PROFILING_STOP('rad_aug_charge')

    END SUBROUTINE RAD_AUG_CHARGE


!*******************************************************************
!
!  RAD_AUG_CHARGE_FOCK
!
!>  Add the compensation charge density on the radial grid
!>  to RHO
!>
!>  RHO(L,M,r) = RHO(L,M,r) +  Q(r,L) Q(PAW,ll' L) RHOLM(ll',LM)
!>
!>  This version adds the fock contributions which
!>  restores the shape of the AE charge more accurately
!>
!>  Q(r,L) are the L-dependent 1-normalized compensation charges.
!
!*******************************************************************

    SUBROUTINE RAD_AUG_CHARGE_FOCK( RHO, R, RHOLM, CHANNELS, L, &
         LYMAX, AUG_FOCK, QPAW_FOCK )
      USE constant
      IMPLICIT NONE

      REAL(q) :: RHO(:,:)    !< charge on radial grid
      TYPE (rgrid) R
      REAL(q) :: RHOLM(:)    !< occupancy of each llLM channel
      REAL(q),POINTER :: AUG_FOCK(:,:,:)   !< 1-normalized L-dep compensation charge
      REAL(q),POINTER :: QPAW_FOCK(:,:,:,:)!< moments of compensation charge Q(PAW,ll L)
      INTEGER :: LYMAX       !< maximum L
      INTEGER CHANNELS, L(:)
   ! local variables
      REAL(q) :: RHOLMT((LYMAX+1)*(LYMAX+1))
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I,LMAX,N
      INTEGER IBASE,JBASE,LMIN,LMAIN,MMAIN,LMMAIN
      INTEGER :: NMAX_FOCKAE, LMAX_FOCKAE

      IF (ASSOCIATED(QPAW_FOCK)) THEN
         NMAX_FOCKAE=SIZE(QPAW_FOCK,4)
         LMAX_FOCKAE=SIZE(QPAW_FOCK,3)-1
      ELSE
         RETURN
      ENDIF

      DO N=1,NMAX_FOCKAE
!-----------------------------------------------------------------------
! first contract to L and M dependent Q(LM)
!-----------------------------------------------------------------------
      IBASE=0
      RHOLMT=0
      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependent charge
         DO LMAIN=LMIN,MIN(LMAX,LMAX_FOCKAE,LYMAX),2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            RHOLMT(LMMAIN)=RHOLMT(LMMAIN)+QPAW_FOCK(CH1,CH2,LMAIN,N)*RHOLM(LMMAIN+JBASE)
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
!-----------------------------------------------------------------------
! then add to charge on radial grid
!-----------------------------------------------------------------------
      ! could be replaced by matrix vecor DGEMV
      DO LMAIN=0,MIN(LMAX_FOCKAE,LYMAX)
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DO I=1,R%NMAX
               RHO(I,LMMAIN)=RHO(I,LMMAIN)+AUG_FOCK(I,LMAIN,N)*RHOLMT(LMMAIN)
            ENDDO
         ENDDO
      ENDDO
      ENDDO
#ifdef debug
      WRITE(0,*) 'RAD_AUG_CHARGE: L compensation charges are'
      DO LL=0,LYMAX
         WRITE(0,'(I3,10F10.5)') LL,(RHOLMT(LL**2+I)*2*SQRT(PI),I=1,LL*2+1)
      ENDDO
#endif
    END SUBROUTINE RAD_AUG_CHARGE_FOCK

!*******************************************************************
!
!  RAD_AUG_PROJ
!
!>  Calculate the integral
!>  D(ll'LM) =  \int V(r,L,M) Q(r,L) Q(PAW,ll' L) dr
!>  on a radial grid.
!
!>  Q(r,L) are the L-dependent 1-normalized compensation charges
!>  the potential is given by
!>       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!>  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
!*******************************************************************

    SUBROUTINE RAD_AUG_PROJ( POT, R, DLM, CHANNELS, L, &
         LYMAX, AUG, QPAW )
      USE constant
      IMPLICIT NONE

      REAL(q) :: POT(:,:)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q) :: AUG(:,0:)    !< 1-normalized L-dep compensation charge
      REAL(q) :: QPAW(:,:,0:) !< moments of compensation charges Q(PAW,ll L)
      INTEGER :: LYMAX        !< maximum L
      INTEGER CHANNELS,L(:)
   ! local variables
      REAL(q) :: RHOLMT((LYMAX+1)*(LYMAX+1)),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I,LMAX
      INTEGER IBASE,JBASE,LMIN,LMAIN,MMAIN,LMMAIN

      PROFILING_START('rad_aug_proj')

!-----------------------------------------------------------------------
! first calculate \int V(L,M) Q(L,M)
!-----------------------------------------------------------------------
      RHOLMT=0
      ! could be replaced by matrix vecor DGEMV
      DO LMAIN=0,LYMAX
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            SUM=0
            DO I=1,R%NMAX
               SUM=SUM+POT(I,LMMAIN)*AUG(I,LMAIN)
            ENDDO
            RHOLMT(LMMAIN)=SUM
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
! than multiply with QPAW(llp, L) and add the DLM
!-----------------------------------------------------------------------
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependet charge
         DO LMAIN=LMIN,MIN(LYMAX,LMAX),2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)-RHOLMT(LMMAIN)*QPAW(CH1,CH2,LMAIN)
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
#ifdef debug
      WRITE(0,*) 'RAD_AUG_PROJ: int V Q(LM) is'
      DO LL=0,LYMAX
         WRITE(0,'(I3,10F10.5)') LL,(RHOLMT(LL**2+I)*2*SQRT(PI),I=1,LL*2+1)
      ENDDO
#endif

      PROFILING_STOP('rad_aug_proj')

    END SUBROUTINE RAD_AUG_PROJ

!*******************************************************************
!
!  RAD_AUG_PROJ_FOCK
!
!>  Calculate the integral
!>   D(ll'LM) =  \int V(r,L,M) Q_FOCK(r,L) Q_FOCK(PAW,ll' L) dr
!>  on a radial grid.
!
!>  Q(r,L) are the L-dependent 0-normalized compensation charges
!>  the potential is given by
!>       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!>  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
!*******************************************************************

    SUBROUTINE RAD_AUG_PROJ_FOCK( POT, R, DLM, CHANNELS, L, &
         LYMAX, AUG_FOCK, QPAW_FOCK )
      USE constant
      IMPLICIT NONE

      REAL(q) :: POT(:,:)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q),POINTER :: AUG_FOCK(:,:,:)    !< 1-normalized L-dep compensation charge
      REAL(q),POINTER :: QPAW_FOCK(:,:,:,:) !< moments of compensation charges Q(PAW,ll L)
      INTEGER :: LYMAX        ! maximum L
      INTEGER CHANNELS,L(:)
   ! local variables
      REAL(q) :: RHOLMT((LYMAX+1)*(LYMAX+1)),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I,LMAX,NAE
      INTEGER IBASE,JBASE,LMIN,LMAIN,MMAIN,LMMAIN
      INTEGER :: NMAX_FOCKAE, LMAX_FOCKAE

      IF (ASSOCIATED(QPAW_FOCK)) THEN
         NMAX_FOCKAE=SIZE(QPAW_FOCK,4)
         LMAX_FOCKAE=SIZE(QPAW_FOCK,3)-1
      ELSE
         RETURN
      ENDIF

      DO NAE=1,NMAX_FOCKAE
!-----------------------------------------------------------------------
! first calculate \int V(L,M) Q(L,M)
!-----------------------------------------------------------------------
      RHOLMT=0
      ! could be replaced by matrix vecor DGEMV
      DO LMAIN=0,MIN(LYMAX,LMAX_FOCKAE)
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            SUM=0
            DO I=1,R%NMAX
               SUM=SUM+POT(I,LMMAIN)*AUG_FOCK(I,LMAIN,NAE)
            ENDDO
            RHOLMT(LMMAIN)=SUM
         ENDDO
      ENDDO
!      WRITE(*,*) '2',RHOLMT(1) ! can be compared to DLM in SETDIJ in us.F
!-----------------------------------------------------------------------
! then multiply with QPAW(llp, L) and add the DLM
!-----------------------------------------------------------------------
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependet charge
         DO LMAIN=LMIN,MIN(LMAX,LMAX_FOCKAE,LYMAX),2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)-RHOLMT(LMMAIN)*QPAW_FOCK(CH1,CH2,LMAIN,NAE)
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
#ifdef debug
      WRITE(0,*) 'RAD_AUG_PROJ: int V Q(LM) is'
      DO LL=0,LYMAX
         WRITE(0,'(I3,10F10.5)') LL,(RHOLMT(LL**2+I)*2*SQRT(PI),I=1,LL*2+1)
      ENDDO
#endif
      ENDDO
    END SUBROUTINE RAD_AUG_PROJ_FOCK


!*******************************************************************
!
!  RAD_CHARGE
!
!>  Calculate the soft pseudo/AE charge density on the radial grid
!>  from a set of frozen wavefunctions and the ll LM dependent
!>  occupancies.
!>
!>  The normalised wavefunction psi(r) can be obtained from the stored 
!>  coefficients w(l,n,r) = W(:,:) in the following way
!>      psi(r) = \sum_lmn Y_lm(r) w_ln(r) / r \n
!>  here 
!>      rho_LM(r)= \sum_n w_ln(r) w_l'n'(r)  RHOLM(ll',LM)
!>  is calculated and stored in RHO(2l+1+m,..) l=0,..,LMAX, m=0,..,2*l \n
!>  RHOLM is the occupancy of each channel (see TRANS_RHOLM)
!>  thus the charge density can be obtained from the stored
!>  coefficients rho_lm(r) in the following way
!>     rho(r) =  \sum_lm rho_lm(r) * Y_lm(r)  / r^2
!
!*******************************************************************

    SUBROUTINE RAD_CHARGE( RHOAE, R, RHOLM, CHANNELS, L, W )
      IMPLICIT NONE

      REAL(q) :: RHOAE(:,:)
      TYPE (rgrid) R
      INTEGER CHANNELS, L(:)
      REAL(q) :: RHOLM(:)
      REAL(q) :: W(:,:)
! local variables
      REAL(q) :: RHOT(R%NMAX)
      INTEGER CH1,CH2,LL,LLP,I
      INTEGER IBASE,JBASE,LMIN,LMAX,LMAIN,MMAIN,LMMAIN
   ! loop over all channels (l,epsilon)
      IBASE=0

      DO CH1=1,CHANNELS
      DO CH2=CH1,CHANNELS
         DO I=1,R%NMAX
            RHOT(I)=W(I,CH1)*W(I,CH2)
         ENDDO
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependent charge
         ! loop could be replaced by matrix vector DGEMV
         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DO I=1,R%NMAX
               RHOAE(I,LMMAIN)=RHOAE(I,LMMAIN)+RHOT(I)*RHOLM(LMMAIN+JBASE)
            ENDDO
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      ENDDO
      ENDDO

    END SUBROUTINE RAD_CHARGE


!*******************************************************************
!
!  RAD_PROJ
!
!>  Calculate
!>  D(ll',LM)= D(ll',LM)+ A *\int dr phi_l'(r) pot_lm(r) phi_l(r) dr
!>  on the radial grid.
!>
!>  The potential is given by
!>       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!>  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
!*******************************************************************

    SUBROUTINE RAD_PROJ( POT, R, A, DLM, CHANNELS, L, W )
      IMPLICIT NONE

      REAL(q) :: POT(:,:)     !< radial potential V(r,L,M)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q) :: W(:,:)       !< wavefunctions phi(r,l)
      REAL(q) :: A            !< scaling factor
      INTEGER CHANNELS, L(:)
! local variables
      REAL(q) :: RHOT(R%NMAX),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I
      INTEGER IBASE,JBASE,LMIN,LMAX,LMAIN,MMAIN,LMMAIN

      PROFILING_START('rad_proj')

   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
         DO I=1,R%NMAX
            RHOT(I)=W(I,CH1)*W(I,CH2)
         ENDDO
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN

         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            SUM=0

            ! integrate RHOT POT(L,M) (the potentials are already weighted)
            DO I=1,R%NMAX
               SUM=SUM+RHOT(I)*POT(I,LMMAIN)
            ENDDO
            DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)+SUM*A
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

      PROFILING_STOP('rad_proj')

    END SUBROUTINE RAD_PROJ

!*******************************************************************
!
!  RAD_PROJ_KINPOT
!
!>  Calculate the integral
!>   D(ll',LM) =  D(ll',LM) + A * < nabla phi_l' | V(L,M) | nabla phi_l >
!>  on the radial grid.
!>
!>  Only spherical contributions are used.
!
!*******************************************************************

    SUBROUTINE RAD_PROJ_KINPOT( POT, R, A, DLM, CHANNELS, L, W )
      IMPLICIT NONE

      REAL(q) :: POT(:)       !< radial kinetic energy potential V(r,L,M)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q) :: W(:,:)       !< wavefunctions phi(r,l)
      REAL(q) :: A            !< scaling factor
      INTEGER CHANNELS, L(:)
! local variables
      REAL(q) :: K(R%NMAX),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I
      INTEGER IBASE,JBASE,LMIN,LMAX,LMAIN,MMAIN,LMMAIN
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
         
         IF (LL == LLP) THEN

            ! here one has to calculate the spherical contribution
            ! to the kinetic energy density

            DO I=1,R%NMAX
               K(I)=W(I,CH1)*W(I,CH2)
            ENDDO

            ! Lmin and Lmax
            LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
            JBASE=IBASE-LMIN*LMIN

            LMAIN=0
            MMAIN=1

              ! currently LMMAIN is allways 1, and JBASE=IBASE
              LMMAIN=LMAIN*LMAIN+MMAIN
              ! integrate K *  POT (the potential is already weighted)
              SUM=0
              DO I=1,R%NMAX
                 SUM=SUM+K(I)*POT(I)
              ENDDO
              DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)+SUM*A
         ENDIF

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE RAD_PROJ_KINPOT

!*******************************************************************
!
!  RAD_INT
!
!>  Integrate the soft (including compensation charge)
!>  and AE chargedensity and calculate.
!>
!>  The moments (if they are not the same we have a problem).
!
!*******************************************************************

    SUBROUTINE RAD_INT(  R,  LYMAX, RHO, RHOAE )
      USE constant
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      REAL(q) :: RHOAE(:,:),RHO(:,:)
      TYPE (rgrid) R
    ! local
      INTEGER LL,I,II
      REAL(q)  RHOL(LYMAX*2+1), RHOLAE(LYMAX*2+1)
      REAL(q)  RHOT(R%NMAX)

      DO LL=0,LYMAX
         DO I=1,LL*2+1
            DO K=1,R%NMAX
               RHOT(K)=RHO(K,LL*LL+I)*R%R(K)**LL
            ENDDO
            CALL SIMPI(R, RHOT, RHOL(I))
            DO K=1,R%NMAX
               RHOT(K)=RHOAE(K,LL*LL+I)*R%R(K)**LL
            ENDDO
            CALL SIMPI(R, RHOT, RHOLAE(I))
         ENDDO
         DO II=1,LL*2+1
         ! we are screwed
         ! the moments are not correctly conserved
         IF ( ABS(RHOL(II)-RHOLAE(II)) >1E-4_q) THEN
            CALL vtutor%bug(str(LL) // " " // str(RHOL(1:LL*2+1)*2*SQRT(PI)) // "\n " // str(LL) // &
               " " // str(RHOLAE(1:LL*2+1)*2*SQRT(PI)) // "\n internal error in RAD_INT: RHOPS /= RHOAE", __FILE__, __LINE__)
         ENDIF
         ENDDO
      ENDDO

    END SUBROUTINE RAD_INT

!*******************************************************************
!
!>  Print out the charge deficit for each L,M channel.
!
!*******************************************************************

    SUBROUTINE RAD_DEF(  R,  LYMAX, RHO, RHOAE )
      USE constant
      REAL(q) :: RHOAE(:,:),RHO(:,:)
      TYPE (rgrid) R
    ! local
      INTEGER LL,I,II
      REAL(q)  RHOL(LYMAX*2+1), RHOLAE(LYMAX*2+1)
      REAL(q)  RHOT(R%NMAX)

      DO LL=0,LYMAX
         DO I=1,LL*2+1
            DO K=1,R%NMAX
               RHOT(K)=RHO(K,LL*LL+I)*R%R(K)**LL
            ENDDO
            CALL SIMPI(R, RHOT, RHOL(I))
            DO K=1,R%NMAX
               RHOT(K)=RHOAE(K,LL*LL+I)*R%R(K)**LL
            ENDDO
            CALL SIMPI(R, RHOT, RHOLAE(I))
         ENDDO
         WRITE(0,'("RADDEF",I2,14F10.6)') LL,(RHOL(II)-RHOLAE(II),II=1,LL*2+1)
      ENDDO

    END SUBROUTINE RAD_DEF

!*******************************************************************
!
!  RAD_POT
!
!>  Calculate the radial potential from the radial chargedensity
!>  the charge density rho(r) is given by
!>     rho(r) =  \sum_lm rho_lm(r) * Y_lm(r)  / r^2
!>
!>  the potential is given by
!>       V(r) =  \sum_lm pot_lm(r) * Y_lm(r) \n
!>  where rho_lm(r) is stored in RHO(2l+1+m,..) l=0,..,LMAX, m=0,..,2*l
!>  and   pot_lm(r) is stored in POT(2l+1+m,..)
!>
!> In many places we use a scaling factor 2 sqrt(pi)
!> the "real" charge density for L=0 angular quantum number is
!>   n_0= rho_00 Y_00 = RHO(r,0) / (2 sqrt(pi)) \n
!> For other channels it is
!>   n_lm= rho_lm Y_lm
!> in comments we will always distinct between n_L and rho_L
!>
!> The charge density is supplied as 
!>  total=RHO(:,:,1), magnetization=RHO(:,:,2)
!> The potential however is returned for spin up and spin down 
!>  spin up=POT(:,:,1), spin down=POT(:,:,2)
!
!*******************************************************************


    SUBROUTINE RAD_POT( R, XC, ISPIN, LMAX, LMAX_CALC, LASPH, &
         RHO, RHOC, POTC, POT, DOUBLEC, EXCG)
      USE constant
      USE setexm_struct_def, ONLY : xc_info

      IMPLICIT NONE
      INTEGER LMAX, ISPIN, LMAX_CALC
      LOGICAL :: LASPH          !< non spherical contributions
      REAL(q) :: RHOC(:)        !< core charge for exchange correlation
      REAL(q) :: POTC(:)        !< froze core potential
      REAL(q) :: RHO(:,:,:)     !< charge distribution
      ! RHO(:,:,1) contains total charge distribution
      ! RHO(:,:,2) magnetization charge distribution
      REAL(q) :: POT(:,:,:)     !< potential
      ! POT(:,:,1) up   component
      ! POT(:,:,2) down component
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      REAL(q) :: EXCG           !< exchange energy only
      REAL(q) :: DOUBLEC        !< double counting corrections
    ! local variables
      REAL(q) RHOT(R%NMAX,ISPIN)
      INTEGER K,N,I,L,M,LM
      REAL(q) SCALE,SUM
      
      LOGICAL,PARAMETER :: TREL=.TRUE. !< use relativistic corrections to exchange
      LOGICAL,PARAMETER :: TLDA=.TRUE. !< calculate LDA contribution seperately
      ! TLDA=.FALSE. works only for Perdew Burke Ernzerhof
      ! in this case non spherical contributions are missing
      REAL(q) :: DHARTREE,DEXC_LDA,DVXC_LDA,DEXC_GGA,DVXC_GGA
      REAL(q) :: TMP((LMAX+1)*(LMAX+1),ISPIN)
!
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DEXC1, DVXC1
      REAL(q) T1(R%NMAX),T2(R%NMAX),V1(R%NMAX)

      PROFILING_START('rad_pot')

      POT=0

      SCALE=2*SQRT(PI)
      N=R%NMAX

      DHARTREE=0
      DO L=0,LMAX_CALC
      DO M=0,2*L
         LM=L*L+M+1
         CALL RAD_POT_HAR(L,R,POT(:,LM,1),RHO(:,LM,1),SUM)
         IF (ISPIN==2) POT(:,LM,2)=POT(:,LM,1)
         DHARTREE=DHARTREE+SUM
         TMP(LM,1)=SUM
      ENDDO
      ! WRITE(0,'(I2,10F12.7)') L, (TMP(LM,1),LM=L*L+1,(L+1)*(L+1))
      ENDDO
      DO K=1,N
         POT(K,1,1)=POT(K,1,1)+POTC(K)*SCALE
      ENDDO
      IF (ISPIN==2) POT(:,1,2)=POT(:,1,1)
!========================================================================
! exchange correlation energy, potential
! and double counting corrections
!========================================================================
      DEXC_LDA=0
      DVXC_LDA=0
      DEXC_GGA=0
      DVXC_GGA=0

      IF (XC%LDOLDA.OR.XC%LDOGGA) THEN

      IF (ISPIN==1) THEN
        DO K=1,N
          RHOT(K,1)=(RHO(K,1,1)+RHOC(K))/ (SCALE*R%R(K)*R%R(K)) ! charge density rho_0 Y(0)
        ENDDO

      ELSE
        DO K=1,N
         RHOT(K,1)=(RHO(K,1,1)+RHOC(K))/ (SCALE*R%R(K)*R%R(K)) ! charge density n_0=rho_0 Y(0)
         RHOT(K,2)= RHO(K,1,2)/ (SCALE*R%R(K)*R%R(K))          ! magnetization
        ENDDO

      ENDIF
      ! WRITE(*,'(9F14.7)') POT(R%NMAX-10,:,:)
!========================================================================
! LDA+ GGA if required
!========================================================================
      IF (.NOT. LASPH) THEN
         ilda: IF (TLDA) THEN
            
            IF (ISPIN==1) THEN
               CALL RAD_LDA_XC( R, XC, TREL, LMAX_CALC, RHOT(:,1), RHO(:,:,1), POT(:,:,1), DEXC_LDA, DVXC_LDA, .TRUE.)
            ELSE
               CALL RAD_LDA_XC_SPIN( R, XC, TREL, LMAX_CALC, RHOT, RHO, POT, DEXC_LDA, DVXC_LDA, .TRUE.)
            ENDIF
            
         ENDIF ilda
         
         gga: IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
            
            IF (ISPIN==1) THEN
               CALL RAD_GGA_XC( R, XC, TLDA, RHOT(:,1), RHO(:,1,1), POT(:,1,1), DEXC_GGA, DVXC_GGA)
            ELSE
               DO K=1,N
                  RHOT(K,1)=(RHO(K,1,1)+RHOC(K)+RHO(K,1,2))/(2*SCALE*R%R(K)*R%R(K)) ! up
                  RHOT(K,1)=MAX(RHOT(K,1), 1E-7_q)
                  RHOT(K,2)=(RHO(K,1,1)+RHOC(K)-RHO(K,1,2))/(2*SCALE*R%R(K)*R%R(K)) ! down
                  RHOT(K,2)=MAX(RHOT(K,2), 1E-7_q)
               ENDDO
               
               CALL RAD_GGA_XC_SPIN( R, XC, TLDA, RHOT, RHO(:,1,:), POT(:,1,:), DEXC_GGA, DVXC_GGA)
            ENDIF
         ENDIF gga
      ELSE
         CALL RAD_GGA_ASPH( R, XC, ISPIN, LMAX_CALC, &
              RHO, RHOC, POTC, DEXC_LDA, DEXC_GGA, DVXC_LDA, DVXC_GGA, POT )
      ENDIF

      ENDIF

!========================================================================
! thats it
!========================================================================
      DOUBLEC= -DHARTREE/2+DEXC_LDA-DVXC_LDA+DEXC_GGA-DVXC_GGA
      EXCG= DEXC_LDA+DEXC_GGA
#ifdef debug
      WRITE(*,1)  N,-DHARTREE/2, -DVXC_LDA-DVXC_GGA, DEXC_LDA+DEXC_GGA
1     FORMAT(' -Hartree, -vxc, exc:',I4,3F14.7)
#endif

      PROFILING_STOP('rad_pot')

    END SUBROUTINE RAD_POT


!*******************************************************************
!
!  RAD_POT_HAR_ONLY
!
!>  Calculate the radial potential from the radial chargedensity.
!>
!>  Hartree contribution only.
!>  See above for details.
!
!*******************************************************************


    SUBROUTINE RAD_POT_HAR_ONLY( R, ISPIN, LMAX, LMAX_CALC, &
         RHO, POT, DOUBLEC )
      USE constant

      IMPLICIT NONE
      INTEGER LMAX, ISPIN, LMAX_CALC
      REAL(q) :: RHO(:,:,:)     !< charge distribution see above
      ! RHO(:,:,1) contains total charge distribution
      ! RHO(:,:,2) magnetization charge distribution
      REAL(q) :: POT(:,:,:)     !< potential
      ! POT(:,:,1) up   component
      ! POT(:,:,2) down component
      TYPE (rgrid) :: R
      REAL(q) :: DOUBLEC
    ! local variables
      
      REAL(q) :: DHARTREE, SCALE, SUM
      INTEGER L,M,LM,N

      POT=0

      SCALE=2*SQRT(PI)
      N=R%NMAX

      DHARTREE=0
      DO L=0,LMAX_CALC
      DO M=0,2*L
         LM=L*L+M+1
         CALL RAD_POT_HAR(L,R,POT(:,LM,1),RHO(:,LM,1),SUM)
         IF (ISPIN==2) POT(:,LM,2)=POT(:,LM,1)
         DHARTREE=DHARTREE+SUM
      ENDDO
      ENDDO
      IF (ISPIN==2) POT(:,1,2)=POT(:,1,1)

      DOUBLEC= -DHARTREE/2

    END SUBROUTINE RAD_POT_HAR_ONLY


    SUBROUTINE RAD_POT_XC_ONLY( R, XC, ISPIN, LMAX, LMAX_CALC, LASPH, &
         RHO, RHOC, POTC, POT, DOUBLEC, EXCG, POTH)
      USE constant
      USE setexm_struct_def, ONLY : xc_info

      IMPLICIT NONE
      INTEGER LMAX, ISPIN, LMAX_CALC
      LOGICAL :: LASPH          !< non spherical contributions
      REAL(q) :: RHOC(:)        !< core charge for exchange correlation
      REAL(q) :: POTC(:)        !< froze core potential
      REAL(q) :: RHO(:,:,:)     !< charge distribution
      ! RHO(:,:,1) contains total charge distribution
      ! RHO(:,:,2) magnetization charge distribution
      REAL(q) :: POT(:,:,:)     !< potential
      ! POT(:,:,1) up   component
      ! POT(:,:,2) down component
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      REAL(q) :: EXCG           !< exchange energy only
      REAL(q) :: DOUBLEC        !< double counting corrections
      REAL(q), OPTIONAL :: POTH(:,:,:) !< optional output needed when LNO_AUG_XC = .TRUE.
    ! local variables
      REAL(q) RHOT(R%NMAX,ISPIN)
      INTEGER K,N
      REAL(q) SCALE

      LOGICAL,PARAMETER :: TREL=.TRUE. !< use relativistic corrections to exchange
      LOGICAL,PARAMETER :: TLDA=.TRUE. !< calculate LDA contribution seperately
      ! TLDA=.FALSE. works only for Perdew Burke Ernzerhof
      ! in this case non spherical contributions are missing
      REAL(q) :: DEXC_LDA,DVXC_LDA,DEXC_GGA,DVXC_GGA
!
      PROFILING_START('rad_pot_xc_only')

      SCALE=2*SQRT(PI)
      N=R%NMAX

      DO K=1,N
         POT(K,1,1)=POT(K,1,1)+POTC(K)*SCALE
      ENDDO
      IF (ISPIN==2) POT(:,1,2)=POT(:,1,1)

      IF (PRESENT(POTH)) POTH = POT

!========================================================================
! exchange correlation energy, potential
! and double counting corrections
!========================================================================
      DEXC_LDA=0
      DVXC_LDA=0
      DEXC_GGA=0
      DVXC_GGA=0

      IF (XC%LDOLDA.OR.XC%LDOGGA) THEN

      IF (ISPIN==1) THEN
        DO K=1,N
          RHOT(K,1)=(RHO(K,1,1)+RHOC(K))/ (SCALE*R%R(K)*R%R(K)) ! charge density rho_0 Y(0)
        ENDDO

      ELSE
        DO K=1,N
         RHOT(K,1)=(RHO(K,1,1)+RHOC(K))/ (SCALE*R%R(K)*R%R(K)) ! charge density n_0=rho_0 Y(0)
         RHOT(K,2)= RHO(K,1,2)/ (SCALE*R%R(K)*R%R(K))          ! magnetization
        ENDDO

      ENDIF
      ! WRITE(*,'(9F14.7)') POT(R%NMAX-10,:,:)
!========================================================================
! LDA+ GGA if required
!========================================================================
      IF (.NOT. LASPH) THEN
         ilda: IF (TLDA) THEN

            IF (ISPIN==1) THEN
               CALL RAD_LDA_XC( R, XC, TREL, LMAX_CALC, RHOT(:,1), RHO(:,:,1), POT(:,:,1), DEXC_LDA, DVXC_LDA, .TRUE.)
            ELSE
               CALL RAD_LDA_XC_SPIN( R, XC, TREL, LMAX_CALC, RHOT, RHO, POT, DEXC_LDA, DVXC_LDA, .TRUE.)
            ENDIF

         ENDIF ilda

         gga: IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN

            IF (ISPIN==1) THEN
               CALL RAD_GGA_XC( R, XC, TLDA, RHOT(:,1), RHO(:,1,1), POT(:,1,1), DEXC_GGA, DVXC_GGA)
            ELSE
               DO K=1,N
                  RHOT(K,1)=(RHO(K,1,1)+RHOC(K)+RHO(K,1,2))/(2*SCALE*R%R(K)*R%R(K)) ! up
                  RHOT(K,1)=MAX(RHOT(K,1), 1E-7_q)
                  RHOT(K,2)=(RHO(K,1,1)+RHOC(K)-RHO(K,1,2))/(2*SCALE*R%R(K)*R%R(K)) ! down
                  RHOT(K,2)=MAX(RHOT(K,2), 1E-7_q)
               ENDDO

               CALL RAD_GGA_XC_SPIN( R, XC, TLDA, RHOT, RHO(:,1,:), POT(:,1,:), DEXC_GGA, DVXC_GGA)
            ENDIF
         ENDIF gga
      ELSE
         CALL RAD_GGA_ASPH( R, XC, ISPIN, LMAX_CALC, &
              RHO, RHOC, POTC, DEXC_LDA, DEXC_GGA, DVXC_LDA, DVXC_GGA, POT )
      ENDIF

      ENDIF

!========================================================================
! thats it
!========================================================================
      DOUBLEC= DOUBLEC + DEXC_LDA-DVXC_LDA+DEXC_GGA-DVXC_GGA
      EXCG= DEXC_LDA+DEXC_GGA
#ifdef debug
      WRITE(*,1)  N, -DVXC_LDA-DVXC_GGA, DEXC_LDA+DEXC_GGA
1     FORMAT(' -vxc, exc:',I4,2F14.7)
#endif

      PROFILING_STOP('rad_pot_xc_only')

    END SUBROUTINE RAD_POT_XC_ONLY


!*******************************************************************
!
!  RAD_POT_DARWIN
!
!>  Calculate darwin shift to adjust the gaps and the dispersion.
!>
!>  The potential is given by
!>       V(r) =  \sum_lm pot_lm(r) * Y_lm(r) \n
!>  and   pot_lm(r) is stored in POT(2l+1+m,..)
!>
!>  This is done by adding to the effective potentials
!>  an additional external potential V_w(r).
!>
!>  The external potential is given by
!>       V_w(r) = V_0 {r_0/r} EXP [ -(r/r_0)^2 ]
!
!*******************************************************************


    SUBROUTINE RAD_POT_DARWIN( R, ISPIN, POT, DARWIN_V, DARWIN_R )
      USE constant

      IMPLICIT NONE
      INTEGER ISPIN
      REAL(q) :: POT(:,:,:)     !< potential
      ! POT(:,:,1) up   component
      ! POT(:,:,2) down component
      REAL(q) :: DARWIN_V, DARWIN_R !< Parameters of the external potential
      TYPE (rgrid) :: R
      ! local variables
      
      REAL(q) :: SCALE, SUM
      INTEGER :: L,M,LM,N

      SCALE=2*SQRT(PI)
      N=R%NMAX

      L=0
      M=0
      LM=L*L+M+1
      DO N=1,R%NMAX
         POT(:,LM,1)=DARWIN_V*(DARWIN_R/R%R(N))*EXP(-(R%R(N)/DARWIN_R)**2_q)*SCALE
         IF (ISPIN==2) POT(:,LM,2)=DARWIN_V*(DARWIN_R/R%R(N))*EXP(-(R%R(N)/DARWIN_R)**2_q)*SCALE
      ENDDO

    END SUBROUTINE RAD_POT_DARWIN

!*******************************************************************
!
!  RAD_GGA_ASPH
!
!> Calculate the one center exchange energy and the potential
!> by performing a angular integration for each radial grid
!> point.
!>
!> Robin Hirschl   20010118 \n
!> potential added 20050225 gK
!
!*******************************************************************
    
    SUBROUTINE RAD_GGA_ASPH( R, XC, ISPIN, LMAX_CALC, &
         RHO, RHOC, POTC, EXC_LDA, EXC_GGA, DVXC_LDA, DVXC_GGA, POT )
      USE constant
      USE asa
      USE setexm_struct_def, ONLY : xc_info, PHPTS_FACT, THPTS_FACT
      USE xc_driver, ONLY : GGA_DRIVER, GGASPIN_DRIVER

      IMPLICIT NONE
      INTEGER ISPIN, LMAX_CALC
      REAL(q) :: RHOC(:)        !< core charge for exchange correlation
      REAL(q) :: POTC(:)        !< froze core potential
      REAL(q) :: RHO(:,:,:)     !< charge distribution see above
      ! RHO(:,:,1) contains total charge distribution
      ! RHO(:,:,2) magnetization charge distribution
      REAL(q), OPTIONAL :: POT(:,:,:)     !< potential
      ! POT(:,:,1) up   component
      ! POT(:,:,2) down component
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      REAL(q) :: EXC_LDA, EXC_GGA, DVXC_LDA, DVXC_GGA
    ! set potential
      LOGICAL :: LPOT           !< set potential or not
    ! local variables
      REAL(q) RHOT(R%NMAX,ISPIN),RHO_ANG(R%NMAX,ISPIN,3)
      REAL(q) POTLDA(R%NMAX,ISPIN)
      INTEGER ISP,K,I,J,LM,LMP,NP,IFAIL,LLMAX,LMMAX,LMMAX_CALC,PHPTS,THPTS,NPTS
      REAL(q) SCALE,SUM,DELTAPHI

      LOGICAL,PARAMETER :: TREL=.TRUE. !< use relativistic corrections to exchange
      LOGICAL,PARAMETER :: TLDA=.TRUE. !< calculate LDA contribution seperately
      ! TLDA=.FALSE. works only for Perdew Burke Ernzerhof

      REAL(q) :: DVINI
!
      REAL(q) SIM_FAKT
      REAL(q) TEMP1
      REAL(q), ALLOCATABLE :: RADPTS(:,:), XYZPTS(:,:),YLM(:,:),YLMD(:,:,:)
      REAL(q), ALLOCATABLE :: YLM_NABLA_YLM(:,:,:),YLM_X_YLM(:,:,:)
      REAL(q), ALLOCATABLE :: MU(:,:,:,:)
      REAL(q), ALLOCATABLE :: WEIGHT(:),ABSCIS(:)

      REAL(q) T1(R%NMAX,ISPIN)
      REAL(q) XU,YU,ZU,XD,YD,ZD,TG,EXT,DEXC,DEXC1,DEXC2,DVXC,DVXC1,DVXC2,DVC
      REAL(q) TMP(ISPIN,3)
      EXTERNAL GAUSSI2

      PROFILING_START('rad_gga_asph')

      IF (PRESENT(POT)) THEN
         LPOT=.TRUE.
      ELSE
         LPOT=.FALSE.
      ENDIF
      ! initial double counting corrections
      DVINI=0
      IF (LPOT) THEN
         CALL RAD_POT_IN_RHO( R, ISPIN, LMAX_CALC, RHO, POT, DVINI)
      ENDIF
      
      SCALE=2*SQRT(PI)                 ! 1/Y00

      LMMAX_CALC=(LMAX_CALC+1)**2
      ! restriction for the angular development
      LLMAX=MIN(6,LMAX_CALC)
      ! in GGA this needs to be increased by 1 (nabla-comment)
      IF (XC%LDOGGA.OR.XC%LDOMETAGGA) LLMAX=LLMAX+1
      LMMAX=(LLMAX+1)**2

!========================================================================
! number of theta and phi pivot points to perform angular integration
! since Exc=f(a*Yllmax,m) we need more pivot points than theoretically 
! needed to integrate Yllmax,m. 
! the factor 2 is the minium, 3 is more accurate
!========================================================================
! update 23.02.2011 by Georg Kresse
! the following lines cause sizeable errors in the one-electron energies e.g. for
! LDA+U calculations on NiO
! all versions up to vasp.5.2.11
!      PHPTS=2*(LLMAX+1)
!      THPTS=2*FLOOR(REAL(LLMAX/2+1,KIND=q))
! this version is much much better, yielding the correct degeneracies in the
! meV regime
      PHPTS=PHPTS_FACT*(LLMAX+1)
      THPTS=THPTS_FACT*FLOOR(REAL(LLMAX/2+1,KIND=q))
      NPTS=PHPTS*THPTS
      DELTAPHI=REAL(2_q*PI/PHPTS,KIND=q)
! allocate arrays
      ALLOCATE(YLM(NPTS,LMMAX),YLMD(NPTS,LMMAX,3),XYZPTS(NPTS,3),RADPTS(NPTS,2), & 
           YLM_NABLA_YLM(LMMAX,LMMAX,0:3),YLM_X_YLM(LMMAX,LMMAX,0:3), &
           WEIGHT(THPTS),ABSCIS(THPTS))
      

      IF (LPOT) THEN
         ALLOCATE (MU(R%NMAX,LMMAX,0:3,ISPIN))
         MU=0
         CALL SETYLM_NABLA_YLM(LLMAX,YLM_NABLA_YLM,YLM_X_YLM)
      ENDIF

      RADPTS=0; WEIGHT=0; ABSCIS=0
      ! set phi positions, equally spaces
      DO I=1,PHPTS
         DO J=1,THPTS
            RADPTS((J-1)*PHPTS+I,2)=(I-1)*DELTAPHI
         ENDDO
      ENDDO
     ! get theta positions (actually get cos(theta)) (Gauss integration)
      CALL GAUSSI(GAUSSI2,-1._q,1._q,0,THPTS,WEIGHT,ABSCIS,IFAIL)
      DO I=1,THPTS
         RADPTS((I-1)*PHPTS+1:I*PHPTS,1)=ABSCIS(I)
      ENDDO
      ! convert radial to cartesian coordinates
      DO I=1,NPTS
         XYZPTS(I,1)=COS(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! x
         XYZPTS(I,2)=SIN(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! y
         XYZPTS(I,3)=RADPTS(I,1)                                 ! z
      ENDDO
  ! get |r| Y_lm on a unit sphere and its derivatives
      YLM=0 ; YLMD=0

      CALL SETYLM_GRAD2(LLMAX,NPTS,YLM,YLMD,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))

!========================================================================
! main loop over all grid points
! prepare the charge density array RHOT
!========================================================================
      EXC_LDA=0
      EXC_GGA=0

! loop over all points in the angular grid
      points: DO NP=1,NPTS

         ! weight of this points
         SIM_FAKT=DELTAPHI*WEIGHT((INT((NP-1)/PHPTS)+1))

         RHOT=0
         RHO_ANG=0

         ! calculate the total charge RHOT = RHO/r^2 (up and down)
         ! and the gradient of the charge RHO_ANG
         DO K=1,R%NMAX
            DO LM=1,LMMAX_CALC
               IF (ISPIN==1) THEN
                  ! rho(r) Y_L(r)
                  RHOT(K,1)     =RHOT(K,1)     +YLM(NP,LM)   *RHO(K,LM,1)
                  ! rho(r) nabla Y_L(r)
                  RHO_ANG(K,1,:)=RHO_ANG(K,1,:)+YLMD(NP,LM,:)*RHO(K,LM,1)
               ELSE
                  RHOT(K,1)=RHOT(K,1)+YLM(NP,LM)*RHO(K,LM,1)  ! total charge
                  RHOT(K,2)=RHOT(K,2)+YLM(NP,LM)*RHO(K,LM,2)  ! magnetization
                  
                  RHO_ANG(K,1,:)=RHO_ANG(K,1,:)+YLMD(NP,LM,:)*RHO(K,LM,1)  ! total charge
                  RHO_ANG(K,2,:)=RHO_ANG(K,2,:)+YLMD(NP,LM,:)*RHO(K,LM,2)  ! magnetization
               ENDIF
            ENDDO
            ! add core charge (spherical) and divide by 1/r^2
            RHOT(K,1)=(RHOT(K,1)+RHOC(K)*YLM(NP,1))/(R%R(K)*R%R(K))
            IF (ISPIN==2) RHOT(K,2)=RHOT(K,2)/(R%R(K)*R%R(K))

            ! divide nabla Y_L through r^3 (1/r is missing in YLMD)
            RHO_ANG(K,:,:)=RHO_ANG(K,:,:)/(R%R(K)*R%R(K))/R%R(K)
         ENDDO

!========================================================================
! LDA part
!========================================================================
         ilda: IF (TLDA) THEN
            DEXC=0
            POTLDA=0
            IF (ISPIN==1) THEN
               CALL RAD_EXCOR_PARA(R, XC, TREL, RHOT(:,1), POTLDA(:,1), DEXC)
            ELSE
               CALL RAD_EXCOR_SPIN(R, XC, TREL, RHOT, POTLDA, DEXC)
            ENDIF

            EXC_LDA= EXC_LDA+DEXC*SIM_FAKT

            IF (LPOT) THEN
               DO LM=1,LMMAX_CALC
                  DO K=1,R%NMAX
                     ! POTLDA must be multiplied by Y_LM to get V_xc(r)
                     POT(K,LM,1:ISPIN)=POT(K,LM,1:ISPIN)+YLM(NP,LM)*SIM_FAKT*POTLDA(K,1:ISPIN)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF ilda
!========================================================================
! GGA if required
!========================================================================
         gga: IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
            
            IF (ISPIN==1) THEN
               CALL GRAD(R,RHOT,T1)

               PROFILING_START('rad_gga_asph|gga_driver')

               DO K=1,R%NMAX
                  ! norm of gradient for spin up
                  XU  =T1(K,1)*XYZPTS(NP,1)+RHO_ANG(K,1,1)
                  YU  =T1(K,1)*XYZPTS(NP,2)+RHO_ANG(K,1,2)
                  ZU  =T1(K,1)*XYZPTS(NP,3)+RHO_ANG(K,1,3)
                  T1(K,1)=SQRT(XU*XU+YU*YU+ZU*ZU)

                  CALL GGA_DRIVER(XC,RHOT(K,1)*AUTOA3,T1(K,1)*AUTOA4,EXT,DEXC1,DVXC1,.NOT.TLDA)
                  EXC_GGA= EXC_GGA+(EXT*RYTOEV)*RHOT(K,1)*(R%R(K)*R%R(K))*R%SI(K)*SIM_FAKT

                  ! develop
                  !             d    f_xc     grad rho
                  !             ------------  --------  = vec mu(r)
                  !             d |grad rho| |grad rho|
                  ! into spherical harmonics
                  ! afterwards vec mu can be reconstructed using
                  ! mu_i(r) =  \sum_L mu_Li(r) Y_L(r)

                  IF (LPOT) THEN
                     DEXC1=DEXC1*RYTOEV
                     DVXC1=DVXC1*RYTOEV*AUTOA

                     DVXC1=DVXC1/ MAX(T1(K,1),1.E-10_q)

                     DO LM=1,LMMAX
                     ! add to local potential to POTXC
                        MU(K,LM,0,1) = MU(K,LM,0,1)+DEXC1*SIM_FAKT*YLM(NP,LM)

                        MU(K,LM,1,1) = MU(K,LM,1,1)+(XU* DVXC1)*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,2,1) = MU(K,LM,2,1)+(YU* DVXC1)*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,3,1) = MU(K,LM,3,1)+(ZU* DVXC1)*SIM_FAKT*YLM(NP,LM)
                     ENDDO
                  ENDIF
               ENDDO

               PROFILING_STOP('rad_gga_asph|gga_driver')

            ELSE
               DO K=1,R%NMAX
                  TEMP1=RHOT(K,1)
                  RHOT(K,1)=(RHOT(K,1)+RHOT(K,2))/2    ! spin up
                  RHOT(K,2)=(TEMP1-RHOT(K,2))/2        ! spin down
                  
                  RHOT(K,1)=MAX(RHOT(K,1), 1E-7_q)
                  RHOT(K,2)=MAX(RHOT(K,2), 1E-7_q)
                ENDDO

               CALL GRAD(R,RHOT,T1)
               CALL GRAD(R,RHOT(1:R%NMAX,2),T1(1:R%NMAX,2))

               PROFILING_START('rad_gga_asph|ggaspin_driver')

               DO K=1,R%NMAX
                  ! norm of gradient for spin up
                  XU  =T1(K,1)*XYZPTS(NP,1)+(RHO_ANG(K,1,1)+RHO_ANG(K,2,1))/2
                  YU  =T1(K,1)*XYZPTS(NP,2)+(RHO_ANG(K,1,2)+RHO_ANG(K,2,2))/2
                  ZU  =T1(K,1)*XYZPTS(NP,3)+(RHO_ANG(K,1,3)+RHO_ANG(K,2,3))/2
                  T1(K,1)=SQRT(XU*XU+YU*YU+ZU*ZU)
                  ! norm of gradient for spin down
                  XD  =T1(K,2)*XYZPTS(NP,1)+(RHO_ANG(K,1,1)-RHO_ANG(K,2,1))/2
                  YD  =T1(K,2)*XYZPTS(NP,2)+(RHO_ANG(K,1,2)-RHO_ANG(K,2,2))/2
                  ZD  =T1(K,2)*XYZPTS(NP,3)+(RHO_ANG(K,1,3)-RHO_ANG(K,2,3))/2

                  T1(K,2)=SQRT(XD*XD+YD*YD+ZD*ZD)
!                  TG=T1(K,1)+T1(K,2) incorrect definition
                  TG=SQRT((XU+XD)**2+(YU+YD)**2+(ZU+ZD)**2)

                  CALL GGASPIN_DRIVER(XC,RHOT(K,1)*AUTOA3,RHOT(K,2)*AUTOA3, &
                       T1(K,1)*AUTOA4,T1(K,2)*AUTOA4,TG*AUTOA4, &
                       EXT,DEXC1,DEXC2,DVXC1,DVXC2,DVC,.NOT.TLDA)
                  EXC_GGA= EXC_GGA+(EXT*RYTOEV)*(RHOT(K,1)+RHOT(K,2))*R%R(K)*R%R(K)*R%SI(K)*SIM_FAKT

                  ! develop
                  !             d    f_xc     grad rho
                  !             ------------  --------  = vec mu(r)
                  !             d |grad rho| |grad rho|
                  ! into spherical harmonics
                  ! afterwards vec mu can be reconstructed using
                  ! mu_i(r) =  \sum_L mu_Li(r) Y_L(r)

                  IF (LPOT) THEN
                     DEXC1=DEXC1*RYTOEV
                     DEXC2=DEXC2*RYTOEV
                     DVXC1=DVXC1*RYTOEV*AUTOA
                     DVXC2=DVXC2*RYTOEV*AUTOA
                     DVC  =DVC  *RYTOEV*AUTOA

                     DVXC1=DVXC1/ MAX(T1(K,1),1.E-10_q)
                     DVXC2=DVXC2/ MAX(T1(K,2),1.E-10_q)
                     DVC  =DVC  / MAX(TG,1.E-10_q)
#ifdef debug
                     IF (K==R%NMAX-10) WRITE(0,'(I4,6F14.7)') NP, & 
                           (XU* DVXC1 + (XU+XD) * DVC ), &
                           (YU* DVXC1 + (YU+YD) * DVC ), &
                           (ZU* DVXC1 + (ZU+ZD) * DVC ), &                  
                           (XD* DVXC2 + (XU+XD) * DVC ), &
                           (YD* DVXC2 + (YU+YD) * DVC ), &
                           (ZD* DVXC2 + (ZU+ZD) * DVC )
#endif
                     DO LM=1,LMMAX
                     ! add to local potential to POTXC
                        MU(K,LM,0,1) = MU(K,LM,0,1)+DEXC1*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,0,2) = MU(K,LM,0,2)+DEXC2*SIM_FAKT*YLM(NP,LM)

                        MU(K,LM,1,1) = MU(K,LM,1,1)+(XU* DVXC1 + (XU+XD) * DVC)*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,2,1) = MU(K,LM,2,1)+(YU* DVXC1 + (YU+YD) * DVC)*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,3,1) = MU(K,LM,3,1)+(ZU* DVXC1 + (ZU+ZD) * DVC)*SIM_FAKT*YLM(NP,LM)
                               
                        MU(K,LM,1,2) = MU(K,LM,1,2)+(XD* DVXC2 + (XU+XD) * DVC)*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,2,2) = MU(K,LM,2,2)+(YD* DVXC2 + (YU+YD) * DVC)*SIM_FAKT*YLM(NP,LM)
                        MU(K,LM,3,2) = MU(K,LM,3,2)+(ZD* DVXC2 + (ZU+ZD) * DVC)*SIM_FAKT*YLM(NP,LM)
                     ENDDO
                  ENDIF
               ENDDO

               PROFILING_STOP('rad_gga_asph|ggaspin_driver')

            ENDIF
         ENDIF gga

      ENDDO points


      ! double counting corrections
      DVXC_LDA=0
      IF (LPOT) THEN
         CALL RAD_POT_IN_RHO( R, ISPIN, LMAX_CALC, RHO, POT, DVXC_LDA)
         DVXC_LDA=DVXC_LDA-DVINI
      ENDIF
!========================================================================
!
! GGA part related to d e_XC / grad n
!
! add -\sum_i Y_L  d (mu_L'i(x) Y_L'(x))/ d x_i to the potential V_L(r)
!
!                          x_i d mu_L'i(|x|)
! -\sum_L'i Y_L(x) Y_L'(x)  -  -----------  + mu_L'i(x)  Y_L nabla_i Y_L'(x)
!                          |x|   d |x|     
!
!========================================================================

      IF (LPOT .AND. (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA)) THEN
#ifdef debug
         WRITE(0,*) 'result of the potential fit'
         DO NP=1,NPTS
         WRITE(0,'(I4)',ADVANCE='NO') NP
         DO ISP=1,ISPIN
            DO I=1,3
               K=R%NMAX-10
               XU=0
               DO LM=1,LMMAX
                  XU=XU+YLM(NP,LM)*MU(K,LM,I,ISP)
               ENDDO
               WRITE(0,'(F14.7)',ADVANCE='NO') XU
            ENDDO
         ENDDO
         WRITE(0,*)
         ENDDO
#endif
         DO ISP=1,ISPIN
            DO LM=1,LMMAX_CALC
               DO K=1,R%NMAX
                  POT(K,LM,ISP)=POT(K,LM,ISP)+MU(K,LM,0,ISP)
               ENDDO
            ENDDO
         ENDDO
         ! nabla comment: 
         ! the matrix elements < Y_L| nabla | Y_L'> are non zero
         ! for L=L'+-1
         DO ISP=1,ISPIN
            DO I=1,3
               DO LMP=1,LMMAX
                  CALL GRAD(R,MU(1:R%NMAX,LMP,I,ISP),T1(1:R%NMAX,1))
                  DO LM=1,LMMAX_CALC
                     IF (ABS(YLM_NABLA_YLM(LM,LMP,I))>1E-4 .OR. ABS(YLM_X_YLM(LM,LMP,I))>1E-4) THEN
                        DO K=1,R%NMAX
                           POT(K,LM,ISP)=POT(K,LM,ISP)-( &
                                MU(K,LMP,I,ISP)*YLM_NABLA_YLM(LM,LMP,I)/ R%R(K)+ & 
                                T1(K,1)*YLM_X_YLM(LM,LMP,I) &
                              )
                        ENDDO
                  ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(MU)
      ENDIF

      ! finally calculate the double counting corrections again
      DVXC_GGA=0
      IF (LPOT) THEN
         CALL RAD_POT_IN_RHO( R, ISPIN, LMAX_CALC, RHO, POT, DVXC_GGA)
         DVXC_GGA=DVXC_GGA-DVXC_LDA-DVINI
      ENDIF

      DEALLOCATE(YLM,YLMD,XYZPTS,RADPTS,YLM_NABLA_YLM,YLM_X_YLM,WEIGHT,ABSCIS)

      PROFILING_STOP('rad_gga_asph')

    END SUBROUTINE RAD_GGA_ASPH


!*******************************************************************
!
!  RAD_POT_IN_RHO
!
!>  Small helper routine to calculate the integral double counting
!>  corrections.
!>~~~
!>    \int V(r) rho(r) r^2 dr d Omega_r
!>    \sum_LM \int pot_LM(r) (rho_LM(r) r^2) Y_LM(r) Y_LM(r) dr d Omega_r
!>    \sum_LM \int pot_LM(r) (rho_LM(r) r^2) dr 
!>    \sum_LM \sum_l POT(r,LM) (RHO(r,LM) SI
!>~~~
!*******************************************************************

    SUBROUTINE RAD_POT_IN_RHO( R, ISPIN, LMAX_CALC, RHO, POT, DVXC)
      IMPLICIT NONE
      TYPE (rgrid) :: R
      INTEGER ISPIN
      INTEGER LMAX_CALC
      REAL(q) :: RHO(:,:,:)     !< charge distribution see above
      REAL(q) :: POT(:,:,:)     !< potential
      REAL(q) :: DVXC
   ! local
      INTEGER LMMAX, LM, K

      PROFILING_START('rad_pot_in_rho')

      LMMAX=(LMAX_CALC+1)**2

      DVXC=0
      IF (ISPIN==1) THEN
         ! non spin polarized case
         DO LM=1,LMMAX
            DO K=1,R%NMAX
               DVXC=DVXC+POT(K,LM,1)*RHO(K,LM,1)*R%SI(K)
            ENDDO
         ENDDO
      ELSE
         ! non spin polarized case
         ! charge is stored as total magnetisation 
         DO LM=1,LMMAX
            DO K=1,R%NMAX
               DVXC=DVXC+POT(K,LM,1)*(RHO(K,LM,1)+RHO(K,LM,2))*R%SI(K)
            ENDDO
         ENDDO
         DO LM=1,LMMAX
            DO K=1,R%NMAX
               DVXC=DVXC+POT(K,LM,2)*(RHO(K,LM,1)-RHO(K,LM,2))*R%SI(K)
            ENDDO
         ENDDO
         DVXC=DVXC/2
      ENDIF

      PROFILING_STOP('rad_pot_in_rho')

    END SUBROUTINE RAD_POT_IN_RHO


!*******************************************************************
!
!  RAD_POT_WEIGHT
!
!>  Because POT will be required only for integration we
!>  multiply now with the weights.
!
!*******************************************************************

    SUBROUTINE RAD_POT_WEIGHT( R, ISPIN, LMAX, POT)
      IMPLICIT NONE
      REAL(q) :: POT(:,:,:)     !< radial potential
      TYPE (rgrid) :: R
      INTEGER LMAX,ISPIN,I,L,M,LM,K

      PROFILING_START('rad_pot_weight')

      DO I=1,ISPIN
      DO L=0,LMAX
      DO M=0,2*L
         LM=L*L+M+1
         DO K=1,R%NMAX
           POT(K,LM,I)=POT(K,LM,I)*R%SI(K)
         ENDDO      
      ENDDO
      ENDDO
      ENDDO

      PROFILING_STOP('rad_pot_weight')

    END SUBROUTINE RAD_POT_WEIGHT

!*******************************************************************
!
!> Coloumb potential
!
!> We use the equation
!>~~~
!> POT_lm(R) = 4 pi/(2l+1)[ 1/R^(l+1) \int_0^R RHO_lm(r) r^(l) dr
!>                       + R^l \int_R^Infinity RHO_lm(r) r^(-l-1) dr ]
!>           = 4 pi \int_0^Infinity dr RHO_lm(r) (r,R)<^l / (r,R)>^(l+1)
!>
!>  rho(r) = \sum_lm  RHO_lm(r) Y_lm(r) / r^2
!>    V(r) = \sum_lm  POT_lm(r) Y_lm(r)
!>~~~
!> which can be obtained by partial integration of
!>~~~
!> V  = \int dR 1/R^2  \int  rho(r) dr
!>~~~
!*******************************************************************


      SUBROUTINE RAD_POT_HAR(LL,R,POT,RHO,DHARTREE)
      USE constant
      IMPLICIT NONE

      TYPE (rgrid) :: R
      INTEGER LL
      REAL(q) :: RHO(:),POT(:)
      REAL(q) T1(R%NMAX),T2(R%NMAX),V1(R%NMAX),V2(R%NMAX),RL(R%NMAX)
      REAL(q) DHARTREE,H3,EXT
      INTEGER N,I,K,L

      N=R%NMAX

      I=0
    ! additional factor r stems from logarithmic grid
      DO K=N,1,-1
         RL(K)=R%R(K)**LL
         I=I+1
         T2(I)=RHO(K)/RL(K)
         T1(I)=RL(K)*R%R(K)*RHO(K)
      ENDDO
      H3=R%H/ 3.0_q
    ! integrate inward (assuming zero moment for grid point NMAX)
    ! V1 = \int_R^Inf RHO(r) r^l dr
      V1(1)=0
      DO L=3,N,2
        V1(L)  =V1(L-2)+H3*(T1(L-2)+4.0_q*T1(L-1)+T1(L))
        V1(L-1)=V1(L-2)+H3*(1.25_q*T1(L-2)+2.0_q*T1(L-1)-0.25_q*T1(L))
      ENDDO
      IF (MOD(N,2)==0) V1(N)=V1(N-2)+H3*(T1(N-2)+4.0_q*T1(N-1)+T1(N))
    ! V2 = \int_R^Inf RHO(r) r^(-l-1) dr
      V2(1)=0
      DO L=3,N,2
         V2(L)  =V2(L-2)+H3*(T2(L-2)+4.0_q*T2(L-1)+T2(L))
         V2(L-1)=V2(L-2)+H3*(1.25_q*T2(L-2)+2.0_q*T2(L-1)-0.25_q*T2(L))
      ENDDO
      IF (MOD(N,2)==0) V2(N)=V2(N-2)+H3*(T2(N-2)+4.0_q*T2(N-1)+T2(N))

     ! EXT total moment (i.e. charge)
     ! EXT-V1(K) = \int_0^R RHO(r) r^l dr 
      EXT=V1(N)
      I=0
      DO K=N,1,-1
         I=I+1
         POT(I)=(V2(K)*RL(I)+(EXT-V1(K))/(R%R(I)*RL(I)))*(FELECT*4*PI/(2*LL+1))
      ENDDO
    ! double counting corrections are simply obtained by summing over grid
    ! since int Y_lm(Omega_r) Y_lm(Omega_r) d_Omega_r= 1
      DHARTREE=0
      DO K=1,N
         DHARTREE=DHARTREE+POT(K)*(RHO(K)*R%SI(K))
      ENDDO
      END SUBROUTINE RAD_POT_HAR




!*******************************************************************
!
!> Calculate the LDA contribution to the exchange correlation
!> energy and to the potentials.
!>
!> Non spin polarised case and below the spinpolarised case.
!>
!> Following Bloechl we use a Taylor expansion of the
!> exchange correlation energy around the spherical density
!>~~~
!>  E_xc = \int d Omega r^2 dr
!>        eps_xc( n_0(r)) + 1/2 v_xc'(n_0) Y_L rho_L(r) Y_L rho_L(r)
!>~~~
!> (tick denotes derivative w.r.t. n_0).
!> The potential is defined as
!>~~~
!>  v_L(r) = var E_xc / var n_L(r) / Y_L \n
!>~~~
!> (since  the factor  Y_L is applied later in the program (see above)
!>  we have to divide by Y_L)
!
!*******************************************************************

   SUBROUTINE RAD_LDA_XC( R, XC, TREL, LMAX_CALC, RHOT, RHO, POT, DEXC, DVXC, LASPH)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY : EXCOR_DER_PARA
      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      LOGICAL TREL      !< relativistic corrections to exchange
      LOGICAL LASPH     !< wether to calculate aspherical corrections
      INTEGER LMAX_CALC !< maximum L
      REAL(q) RHOT(:)   !< spherical charge + core charge density
      REAL(q) RHO(:,:)  !< charge distribution see above
      REAL(q) POT(:,:)  !< potential
      REAL(q) DEXC      !< exchange correlation energy
      REAL(q) DVXC      !< V_xc rho
! local
      REAL(q) EXCA(4)
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DVXC1, RHOPA, SCALE
      INTEGER K,LM

      PROFILING_START('rad_lda_xc')

      SCALE=2*SQRT(PI)
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &                                                                                   
!$OMP REDUCTION(+:DEXC,DVXC) &
!$OMP PRIVATE(K,EXCA,SIM_FAKT,RHOP,EXT,VXT,LM,RHOPA,DVXC1) &
!$OMP SHARED(R,RHOT,LMAX_CALC,TREL,RHO,POT,SCALE,LASPH,XC)
      DO K=1,R%NMAX
! MM
!        IF (RHOT(K) <= 0 ) CYCLE
         IF (RHOT(K) <= 1E-99_q ) CYCLE
! MM
         CALL EXCOR_DER_PARA(XC,RHOT(K),LMAX_CALC,EXCA,TREL)

         SIM_FAKT= R%SI(K)*R%R(K)*R%R(K)
         RHOP    = RHO(K,1)/(R%R(K)*R%R(K))

         ! store v_xc(r) / Y_0
         POT(K,1)=POT(K,1)+EXCA(2)*SCALE
         !  \int v_xc(r)  rho_0(r) Y_0 4 Pi r^2 dr
         DVXC    =DVXC    +EXCA(2) *  RHOP*SCALE*SIM_FAKT
         !  \int eps_xc(r) 4 pi r^2 dr
         DEXC    =DEXC    +EXCA(1) *4*PI*SIM_FAKT
         EXT=0
         VXT=0
         IF (LASPH) THEN
         DO LM=2,(LMAX_CALC+1)*(LMAX_CALC+1)
            RHOPA=RHO(K,LM)/ (R%R(K)*R%R(K)) ! rho_L
            ! corrections to energy
            ! 1/2 \int v_xc p(n_0) Y_L rho_L(r) Y_L rho_L(r) d Omega r^2 dr
            !  = 1/2 \int v_xc p(n_0) rho_L(r) rho_L(r) r^2 dr
            EXT =EXT +(EXCA(3)*RHOPA*RHOPA)/2
            ! correction to L=0 potential
            ! 1/2  \int v_xc''(n_0) Y_L rho_L(r) Y_L rho_L(r) d Omega =
            ! \int Y_0 (v_xc''(n_0) rho_L(r) rho_L(r) Y_0/2) d Omega
            DVXC1 = (EXCA(4)*RHOPA*RHOPA)/2/SCALE
            ! \int 1/2 v_xc''(n_0) Y(L) rho_L(r) Y(L) rho_L(r) Y(0) rho_0 d Omega r^2 dr
            ! =\int  Y_0/2 v_xc''(n_0) rho_L(r) rho_L(r)  r^2 dr
            VXT=VXT + DVXC1*RHOP
            POT(K,1)=POT(K,1)+DVXC1

            ! L/=0 potentials
            ! \int Y_L (v_xc p(n_0) \rho_L)
            POT(K,LM)=POT(K,LM)+ (EXCA(3)*RHOPA)
            ! \int v_xc p(n_0) \rho_L Y_L \rho_L Y_L d Omega r^2 dr
            ! =  \int v_xc p(n_0)  \rho_L \rho_L r^2 dr
            VXT=VXT+EXCA(3)*RHOPA*RHOPA
         ENDDO
         ENDIF
         DEXC=DEXC + EXT*SIM_FAKT
         DVXC=DVXC + VXT*SIM_FAKT
      ENDDO
!$OMP END PARALLEL DO

      PROFILING_STOP('rad_lda_xc')

   END SUBROUTINE RAD_LDA_XC


   SUBROUTINE RAD_EXCOR_PARA( R, XC, TREL, RHOT, POT, DEXC)
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY : EXCOR_PARA
      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      LOGICAL TREL    !< relativistic corrections to exchange
      REAL(q) RHOT(:) !< spherical charge + core charge density
      REAL(q) POT(:)  !< potential
      REAL(q) DEXC    !< exchange correlation energy
! local
      REAL(q) EXC, EXCD, SIM_FAKT
      INTEGER K

      PROFILING_START('rad_excor_para')

      DO K=1,R%NMAX
         IF (RHOT(K) <= 0) CYCLE

         CALL EXCOR_PARA(XC, RHOT(K), EXC, EXCD, TREL)

         ! store v_xc(r)
         POT(K)=POT(K)+EXCD              ! up potential

         SIM_FAKT= R%SI(K)*R%R(K)*R%R(K) ! r^2 dr

         !  \int eps_xc(r) r^2 dr
         DEXC =DEXC  +EXC * SIM_FAKT
      ENDDO

      PROFILING_STOP('rad_excor_para')

    END SUBROUTINE RAD_EXCOR_PARA


   SUBROUTINE RAD_LDA_XC_SPIN( R, XC, TREL, LMAX_CALC, RHOT, RHO, POT, DEXC, DVXC, LASPH)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY : EXCOR_DER
      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      LOGICAL TREL      !< relativistic corrections to exchange
      LOGICAL LASPH     !< wether to calculate aspherical corrections
      INTEGER LMAX_CALC !< maximum L
      REAL(q) RHOT(:,:) !< spherical charge + core charge density
      REAL(q) RHO(:,:,:)!< charge distribution see above
      REAL(q) POT(:,:,:)!< potential
      REAL(q) DEXC      !< exchange correlation energy
      REAL(q) DVXC      !< V_xc rho
! local
      REAL(q) EXC,EXCD(2),EXCDD(2,2),EXCDDD(2,2,2)
      REAL(q) SIM_FAKT, RHOP, RHOM, EXT, VXT, DVXC1, DVXC2, RHOPA, RHOMA, SCALE, RHOMA_, RHOPA_
      INTEGER K,LM

      PROFILING_START('rad_lda_xc_spin')

      SCALE=2*SQRT(PI)

      DO K=1,R%NMAX
         IF (RHOT(K,1) <= 0) CYCLE

         CALL EXCOR_DER(XC,0.5*(RHOT(K,1)+RHOT(K,2)),0.5*(RHOT(K,1)-RHOT(K,2)), &
              LMAX_CALC,EXC,EXCD,EXCDD,EXCDDD,TREL)

         ! store v_xc(r) / Y_0
         POT(K,1,1)=POT(K,1,1)+(EXCD(1)*SCALE)                 ! up potential
         POT(K,1,2)=POT(K,1,2)+(EXCD(2)*SCALE)                 ! down potential pot_L(r)

         RHOP =0.5_q*(RHO(K,1,1)+RHO(K,1,2))/ (R%R(K)*R%R(K))  ! up valence density   rho_0
         RHOM =0.5_q*(RHO(K,1,1)-RHO(K,1,2))/ (R%R(K)*R%R(K))  ! down valence density

         SIM_FAKT= R%SI(K)*R%R(K)*R%R(K)                       ! r^2 dr

         !  \int eps_xc(r) 4 pi r^2 dr
         DEXC =DEXC  +EXC * SIM_FAKT *4*PI
         !  \int v_xc(r)  rho_0(r) Y_0 4 Pi r^2 dr
         DVXC =DVXC  +(EXCD(1)*RHOP + EXCD(2)*RHOM)* SIM_FAKT*SCALE

         ! L/=0 terms
         ! mind that the factors Y_L are not stored in POT (see above)
         EXT=0
         VXT=0

         IF (LASPH) THEN
         DO LM=2,(LMAX_CALC+1)*(LMAX_CALC+1)
            RHOPA=0.5_q*(RHO(K,LM,1)+RHO(K,LM,2))/ (R%R(K)*R%R(K)) ! rho_L
            RHOMA=0.5_q*(RHO(K,LM,1)-RHO(K,LM,2))/ (R%R(K)*R%R(K))
            ! corrections to energy
            ! 1/2 \int v_xc p(n_0) Y_L rho_L(r) Y_L rho_L(r) d Omega r^2 dr
            !  = 1/2 \int v_xc p(n_0) rho_L(r) rho_L(r) r^2 dr
            EXT =EXT +(EXCDD(1,1)*RHOPA*RHOPA + EXCDD(2,2)*RHOMA*RHOMA + &
                     2*EXCDD(1,2)*RHOPA*RHOMA)/2

            ! correction to L=0 potential
            ! 1/2 v_xc''(n_0) Y_L rho_L(r) Y_L rho_L(r)
            ! spherically averaged:  Y_0 (v_xc''(n_0) rho_L(r) rho_L(r) Y_0/2)
             DVXC1 = (EXCDDD(1,1,1)*RHOPA*RHOPA + EXCDDD(1,2,2)*RHOMA*RHOMA + &
                        2*EXCDDD(1,1,2)*RHOPA*RHOMA)/2/SCALE
             DVXC2 = (EXCDDD(2,1,1)*RHOPA*RHOPA + EXCDDD(2,2,2)*RHOMA*RHOMA + &
                        2*EXCDDD(2,1,2)*RHOPA*RHOMA)/2/SCALE
            ! \int 1/2 v_xc''(n_0) Y(L) rho_L(r) Y(L) rho_L(r) Y(0) rho_0 d Omega r^2 dr
            ! =\int  Y_0/2 v_xc''(n_0) rho_L(r) rho_L(r)  r^2 dr
             POT(K,1,1)=POT(K,1,1)+DVXC1
             POT(K,1,2)=POT(K,1,2)+DVXC2
            ! double counting correction
             VXT=VXT+(DVXC1*RHOP) + (DVXC2*RHOM)

            ! L/=0 potentials
            ! Y_L (v_xc p(n_0) \rho_L)
             POT(K,LM,1)=POT(K,LM,1)+(EXCDD(1,1)*RHOPA+EXCDD(1,2)*RHOMA)
             POT(K,LM,2)=POT(K,LM,2)+(EXCDD(1,2)*RHOPA+EXCDD(2,2)*RHOMA)

            ! double counting correction
            ! \int v_xc p(n_0) \rho_L Y_L \rho_L Y_L d Omega r^2 dr
            ! =  \int v_xc p(n_0)  \rho_L \rho_L r^2 dr
             VXT=VXT+(EXCDD(1,1)*RHOPA*RHOPA +EXCDD(2,2)*RHOMA*RHOMA + &
                    2*EXCDD(1,2)*RHOPA*RHOMA)
         ENDDO
         ENDIF

         DEXC=DEXC + EXT*SIM_FAKT
         DVXC=DVXC + VXT*SIM_FAKT
      ENDDO

      PROFILING_STOP('rad_lda_xc_spin')

    END SUBROUTINE RAD_LDA_XC_SPIN


    SUBROUTINE RAD_EXCOR_SPIN( R, XC, TREL, RHOT, POT, DEXC)
      USE setexm_struct_def, ONLY : xc_info
      USE ldalib, ONLY : EXCOR
      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      LOGICAL TREL      !< relativistic corrections to exchange
      REAL(q) RHOT(:,:) !< spherical charge + core charge density
      REAL(q) POT(:,:)  !< potential
      REAL(q) DEXC      !< exchange correlation energy
! local
      REAL(q) TMP(2), EXC, EXCD(2), SIM_FAKT
      INTEGER K

      PROFILING_START('rad_excor_spin')

      DO K=1,R%NMAX
         IF (RHOT(K,1) <= 0) CYCLE

         TMP(1) = 0.5*(RHOT(K,1)+RHOT(K,2))
         TMP(2) = 0.5*(RHOT(K,1)-RHOT(K,2))

         CALL EXCOR(XC, TMP, EXC, EXCD, TREL)

         ! store v_xc(r)
         POT(K,1)=POT(K,1)+EXCD(1)       ! up potential
         POT(K,2)=POT(K,2)+EXCD(2)       ! down potential

         SIM_FAKT= R%SI(K)*R%R(K)*R%R(K) ! r^2 dr

         !  \int eps_xc(r) r^2 dr
         DEXC =DEXC  +EXC * SIM_FAKT
      ENDDO

      PROFILING_STOP('rad_excor_spin')

    END SUBROUTINE RAD_EXCOR_SPIN


!*******************************************************************
!
!> Calculate the GGA contribution to the exchange correlation
!> energy and to the potentials.
!> Only spherical contributions are accounted for.
!
!*******************************************************************

    SUBROUTINE RAD_GGA_XC( R, XC, TLDA, RHOT, RHO, POT, DEXC_GGA, DVXC_GGA)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      USE xc_driver, ONLY : GGA_DRIVER
      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      LOGICAL TLDA      !< include LDA contributions  (usually .FALSE.)
      REAL(q) RHOT(:)   !< spherical charge + core charge density
      REAL(q) RHO(:)    !< spherical charge
      REAL(q) POT(:)    !< potential
      REAL(q) DEXC_GGA  !< exchange correlation energy
      REAL(q) DVXC_GGA  !< V_xc rho
! local
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DEXC1, DVXC1, SCALE
      REAL(q) T1(R%NMAX),T2(R%NMAX),V1(R%NMAX)
      INTEGER K

      PROFILING_START('rad_gga_xc')

      SCALE=2*SQRT(PI)

      CALL GRAD(R,RHOT,T1)
      DO K=1,R%NMAX
#ifdef vector
         CALL GGA91_WB(XC,RHOT(K)*AUTOA3,T1(K)*AUTOA4,EXT,DEXC1,DVXC1)
#else
         CALL GGA_DRIVER(XC,RHOT(K)*AUTOA3,T1(K)*AUTOA4,EXT,DEXC1,DVXC1,.NOT.TLDA)
#endif
         SIM_FAKT=R%SI(K)*SCALE
         DEXC_GGA=DEXC_GGA+(EXT*RYTOEV)*RHOT(K)*(SCALE*R%R(K)*R%R(K))*SIM_FAKT
         
         !  store d f/ d (d rho )   in T2
         T2(K) = DVXC1*RYTOEV*AUTOA
         !  store d f/ d rho  in T1
         T1(K) = DEXC1*RYTOEV
      ENDDO
      CALL GRAD(R,T2,V1)
      
      DO K=1,R%NMAX
         VXT     = T1(K) - V1(K) - 2*T2(K)/ R%R(K)
         SIM_FAKT=R%SI(K)*SCALE
         DVXC_GGA=DVXC_GGA+VXT*RHO(K)*SIM_FAKT
         POT(K)=POT(K)+VXT*SCALE
      ENDDO

      PROFILING_STOP('rad_gga_xc')

    END SUBROUTINE RAD_GGA_XC


    SUBROUTINE RAD_GGA_XC_SPIN( R, XC, TLDA, RHOT, RHO, POT, DEXC_GGA, DVXC_GGA)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      USE xc_driver, ONLY : GGASPIN_DRIVER
      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      LOGICAL TLDA      !< include LDA contributions  (usually .FALSE.)
      REAL(q) RHOT(:,:) !< up and down charge density (including core)
      REAL(q) RHO(:,:)  !< total charge and magnetisation
      REAL(q) POT(:,:)  !< potential
      REAL(q) DEXC_GGA  !< exchange correlation energy
      REAL(q) DVXC_GGA  !< V_xc rho
! local
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DEXC1, DEXC2, DVXC1, DVXC2, DVC, SCALE
      REAL(q) T1(R%NMAX,2),T2(R%NMAX,2),V1(R%NMAX,2)
      INTEGER K

      PROFILING_START('rad_gga_xc_spin')

      SCALE=2*SQRT(PI)

      CALL GRAD(R,RHOT,T1)
      CALL GRAD(R,RHOT(1:R%NMAX,2),T1(1:R%NMAX,2))

      DO K=1,R%NMAX
         CALL GGASPIN_DRIVER(XC,RHOT(K,1)*AUTOA3,RHOT(K,2)*AUTOA3, &
     &              T1(K,1)*AUTOA4,T1(K,2)*AUTOA4,(T1(K,1)+T1(K,2))*AUTOA4, &
     &                EXT,DEXC1,DEXC2,DVXC1,DVXC2,DVC,.NOT.TLDA)
         DVXC1=DVXC1+DVC
         DVXC2=DVXC2+DVC

         SIM_FAKT=R%SI(K)*SCALE
         DEXC_GGA=DEXC_GGA+(EXT*RYTOEV)*(RHOT(K,1)+RHOT(K,2))*SCALE*R%R(K)*R%R(K)*SIM_FAKT

         !   store d f/ d ( d rho )  in T2
         T2(K,1)  = DVXC1*RYTOEV*AUTOA
         T2(K,2)  = DVXC2*RYTOEV*AUTOA

         !   store d f/ d rho  in T1
         T1(K,1) = DEXC1*RYTOEV
         T1(K,2) = DEXC2*RYTOEV
      ENDDO

      CALL GRAD(R,T2,V1)
      CALL GRAD(R,T2(1,2),V1(1,2))

      DO K=1,R%NMAX
         SIM_FAKT=R%SI(K)*SCALE
         
         VXT     = T1(K,1)- V1(K,1) - 2*T2(K,1)/ R%R(K)
         DVXC_GGA=DVXC_GGA+VXT* (RHO(K,1)+RHO(K,2))*0.5_q*SIM_FAKT
         POT(K,1)=POT(K,1)+VXT*SCALE
         
         VXT     = T1(K,2)- V1(K,2) - 2*T2(K,2)/ R%R(K)
         DVXC_GGA=DVXC_GGA+VXT* (RHO(K,1)-RHO(K,2))*0.5_q*SIM_FAKT
         POT(K,2)=POT(K,2)+VXT*SCALE
      ENDDO

      PROFILING_STOP('rad_gga_xc_spin')

    END SUBROUTINE RAD_GGA_XC_SPIN


!*******************************************************************
!
!> This subroutine calculates the absolute magnitude of the 
!> magnetization density for non collinear calculations.
!
!*******************************************************************

      SUBROUTINE RAD_MAG_DENSITY( RHO, RHOCOL, LMAX, R)
      USE constant
      IMPLICIT NONE
      REAL(q) :: RHO(:,:,:)     !< density
      REAL(q) :: RHOCOL(:,:,:)  !< total density
      TYPE (rgrid) :: R
      INTEGER LMAX,ISPIN,I,L,M,LM,K
      REAL(q) AVMAG(3),ABS_AVMAG
      !
      ! L=0 just copy rho, and calculate the absolute
      ! magnitude of m
      LM=1
      AVMAG=0
      DO K=1,R%NMAX
         RHOCOL(K,LM,1)=RHO(K,LM,1)
         RHOCOL(K,LM,2)=SQRT(RHO(K,LM,2)*RHO(K,LM,2)+RHO(K,LM,3)*RHO(K,LM,3)+ &
              RHO(K,LM,4)*RHO(K,LM,4))
         AVMAG(1)=AVMAG(1)+RHO(K,LM,2)*R%SI(K)
         AVMAG(2)=AVMAG(2)+RHO(K,LM,3)*R%SI(K)
         AVMAG(3)=AVMAG(3)+RHO(K,LM,4)*R%SI(K)
      ENDDO
      ABS_AVMAG=SQRT(AVMAG(1)**2+AVMAG(2)**2+AVMAG(3)**2)

      !
      ! L/=0 calculate  m_0 / |m_0| m_lm
      !
      IF (.NOT. USE_AVERAGE_MAGNETISATION) THEN
         DO L=1,LMAX
         DO M=0,2*L
            LM=L*L+M+1
            DO K=1,R%NMAX
               RHOCOL(K,LM,1)=RHO(K,LM,1)
               RHOCOL(K,LM,2)=(RHO(K,1,2)*RHO(K,LM,2)+RHO(K,1,3)*RHO(K,LM,3)+ &
               RHO(K,1,4)*RHO(K,LM,4))/MAX(RHOCOL(K,1,2),MAGMIN)
            ENDDO
         ENDDO
         ENDDO
      ELSE
         DO L=1,LMAX
         DO M=0,2*L
            LM=L*L+M+1
            DO K=1,R%NMAX
               RHOCOL(K,LM,1)=RHO(K,LM,1)
               RHOCOL(K,LM,2)=(AVMAG(1)*RHO(K,LM,2)+AVMAG(2)*RHO(K,LM,3)+ &
               AVMAG(3)*RHO(K,LM,4))/MAX(ABS_AVMAG,MAGMIN)
            ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE RAD_MAG_DENSITY

      SUBROUTINE RAD_MAG_DENSITY_KINDENS(RHO, RHOCOL, KINDENS, KINDENSCOL, LMAX, LMAX_TAU, R, RHOPS)
      USE constant
      IMPLICIT NONE
      REAL(q) :: RHO(:,:,:)        !< density
      REAL(q) :: RHOCOL(:,:,:)     !< total density
      REAL(q) :: KINDENS(:,:,:)    !< kinetic energy density (total,magnetization)
      REAL(q) :: KINDENSCOL(:,:,:) !< return: collinear kinetic energy density (up,down)
      TYPE (rgrid) :: R
      REAL(q), OPTIONAL :: RHOPS(:,:,:) !< the non-augmented pseudo density (treated like RHO) 
      INTEGER LMAX
      INTEGER LMAX_TAU
      ! local variables
      INTEGER ISPIN,I,L,M,LM,K
      REAL(q) AVMAG(3),ABS_AVMAG,A,B

      RHOCOL=0; KINDENSCOL=0

      AVMAG=0
      DO K=1,R%NMAX
         AVMAG(1)=AVMAG(1)+RHO(K,1,2)*R%SI(K)
         AVMAG(2)=AVMAG(2)+RHO(K,1,3)*R%SI(K)
         AVMAG(3)=AVMAG(3)+RHO(K,1,4)*R%SI(K)
      ENDDO
      ABS_AVMAG=SQRT(AVMAG(1)**2+AVMAG(2)**2+AVMAG(3)**2)

      !
      ! L=0 component of the kinetic energy density is
      ! computed seperately, so that at least this component
      ! is set in case USE_AVERAGE_MAGNETISATION=.FALSE.
      DO K=1,R%NMAX
         A=KINDENS(K,1,1)
         B=(AVMAG(1)*KINDENS(K,1,2)+AVMAG(2)*KINDENS(K,1,3)+ &
        &    AVMAG(3)*KINDENS(K,1,4))/MAX(ABS_AVMAG,MAGMIN)
         KINDENSCOL(K,1,1)=(A+B)/2
         KINDENSCOL(K,1,2)=(A-B)/2
      ENDDO

      IF (.NOT. USE_AVERAGE_MAGNETISATION) THEN
         !
         ! L=0 just copy rho, and calculate the absolute
         ! magnitude of m
         DO K=1,R%NMAX
            RHOCOL(K,1,1)=RHO(K,1,1)
            RHOCOL(K,1,2)=SQRT(RHO(K,1,2)*RHO(K,1,2)+RHO(K,1,3)*RHO(K,1,3)+ &
                 RHO(K,1,4)*RHO(K,1,4))
         ENDDO
         !
         ! L/=0 calculate  m_0 / |m_0| m_lm
         !
         DO L=1,LMAX
         DO M=0,2*L
            LM=L*L+M+1
            DO K=1,R%NMAX
               RHOCOL(K,LM,1)=RHO(K,LM,1)
               RHOCOL(K,LM,2)=(RHO(K,1,2)*RHO(K,LM,2)+RHO(K,1,3)*RHO(K,LM,3)+ &
              &   RHO(K,1,4)*RHO(K,LM,4))/MAX(RHOCOL(K,1,2),MAGMIN)
            ENDDO
         ENDDO
         ENDDO
         ! possibly we have a non-augmented pseudo charge density as well
         IF (PRESENT(RHOPS)) THEN
            DO L=0,LMAX
            DO M=0,2*L
               LM=L*L+M+1
               DO K=1,R%NMAX
                  RHOPS(K,LM,1)=RHOPS(K,LM,1)
                  RHOPS(K,LM,2)=(RHO(K,1,2)*RHOPS(K,LM,2)+RHO(K,1,3)*RHOPS(K,LM,3)+ &
                 &   RHO(K,1,4)*RHOPS(K,LM,4))/MAX(RHOCOL(K,1,2),MAGMIN)
               ENDDO
            ENDDO
            ENDDO
         ENDIF
      ELSE
         ! charge density
         DO L=0,LMAX
         DO M=0,2*L
            LM=L*L+M+1
            DO K=1,R%NMAX
               RHOCOL(K,LM,1)=RHO(K,LM,1)
               RHOCOL(K,LM,2)=(AVMAG(1)*RHO(K,LM,2)+AVMAG(2)*RHO(K,LM,3)+ &
              &    AVMAG(3)*RHO(K,LM,4))/MAX(ABS_AVMAG,MAGMIN)
            ENDDO
         ENDDO
         ENDDO
         ! kinetic energy density
         DO L=1,LMAX_TAU
         DO M=0,2*L
            LM=L*L+M+1
            DO K=1,R%NMAX
               A=KINDENS(K,LM,1)
               B=(AVMAG(1)*KINDENS(K,LM,2)+AVMAG(2)*KINDENS(K,LM,3)+ &
              &    AVMAG(3)*KINDENS(K,LM,4))/MAX(ABS_AVMAG,MAGMIN)
               KINDENSCOL(K,LM,1)=(A+B)/2
               KINDENSCOL(K,LM,2)=(A-B)/2
            ENDDO
         ENDDO
         ENDDO
         ! possibly we have a non-augmented pseudo charge density as well
         IF (PRESENT(RHOPS)) THEN
            DO L=0,LMAX
            DO M=0,2*L
               LM=L*L+M+1
               DO K=1,R%NMAX
                  RHOPS(K,LM,1)=RHOPS(K,LM,1)
                  RHOPS(K,LM,2)=(AVMAG(1)*RHOPS(K,LM,2)+AVMAG(2)*RHOPS(K,LM,3)+ &
                 &    AVMAG(3)*RHOPS(K,LM,4))/MAX(ABS_AVMAG,MAGMIN)
               ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDIF

      END SUBROUTINE RAD_MAG_DENSITY_KINDENS


!************************ SUBROUTINE  MAG_DIRECTION  *******************
!
!> On entry POT must contain the v_xc(up) and v_xc(down) \n
!> on return \n 
!>  POT(:,1) =  (v_xc(up) + v_xc(down))/2 \n
!>  POT(:,2) =  m_x (v_xc(up) - v_xc(down))/2 \n
!>  POT(:,3) =  m_y (v_xc(up) - v_xc(down))/2 \n
!>  POT(:,4) =  m_z (v_xc(up) - v_xc(down))/2 \n
!> where m is the unit vector of the magnetization
!
!***********************************************************************


      SUBROUTINE RAD_MAG_DIRECTION( RHO, RHOCOL, POT, LMAX, R)
      IMPLICIT NONE
      REAL(q) :: RHO(:,:,:)     !< density
      REAL(q) :: POT(:,:,:)     !< potential
      REAL(q) :: RHOCOL(:,:,:)  !< total density
      TYPE (rgrid) :: R
      INTEGER LMAX,ISPIN,I,L,M,LM,K
      REAL(q) :: NORM,DELTAV,V0,MLMPROJ
      REAL(q) :: AVMAG(3),ABS_AVMAG

      LM=1
      AVMAG=0
      DO K=1,R%NMAX
         V0       =(POT(K,LM,1)+POT(K,LM,2))/2
         DELTAV   =(POT(K,LM,1)-POT(K,LM,2))/2
!        NORM     = MAX(RHOCOL(K,LM,2),MAGMIN)
         NORM     = MAX(ABS(RHOCOL(K,LM,2)),MAGMIN)*SIGN(1._q,RHOCOL(K,LM,2))
         POT(K,LM,1) = V0
         POT(K,LM,2) = DELTAV * RHO(K,LM,2) / NORM
         POT(K,LM,3) = DELTAV * RHO(K,LM,3) / NORM
         POT(K,LM,4) = DELTAV * RHO(K,LM,4) / NORM
         AVMAG(1)=AVMAG(1)+RHO(K,LM,2)*R%SI(K)
         AVMAG(2)=AVMAG(2)+RHO(K,LM,3)*R%SI(K)
         AVMAG(3)=AVMAG(3)+RHO(K,LM,4)*R%SI(K)
      ENDDO
      ABS_AVMAG=SQRT(AVMAG(1)**2+AVMAG(2)**2+AVMAG(3)**2)

      IF (.NOT. USE_AVERAGE_MAGNETISATION) THEN
      !
      ! for aspherical elements the subroutine has returned
      !           e_xc              
      ! d ------------------  this yields two contributions
      !    m_0 . m_lm  / |m_0|   
      !   e_xc             e_xc          m_lm    (m_0 m_lm) m_0    
      ! d ---- = d ------------------ (  ----  - ---------  ---  )
      !   m_0       m_0 m_lm  / |m_0|    |m_0|   |m_0|^2   |m_0|
      !   e_xc             e_xc            m_0
      ! d ---- = d ------------------     ----
      !   m_lm       m_0 m_lm  / |m_0|    |m_0|

         DO L=1,LMAX
         DO M=0,2*L
         LM=L*L+M+1
         DO K=1,R%NMAX

            V0    =(POT(K,LM,1)+POT(K,LM,2))/2
            DELTAV=(POT(K,LM,1)-POT(K,LM,2))/2

            NORM     = MAX(RHOCOL(K,1,2),MAGMIN)
            POT(K,LM,1) = V0
            POT(K,LM,2) = DELTAV * RHO(K,1,2) / NORM
            POT(K,LM,3) = DELTAV * RHO(K,1,3) / NORM
            POT(K,LM,4) = DELTAV * RHO(K,1,4) / NORM

            MLMPROJ = (RHO(K,LM,2) *RHO(K,1,2) +RHO(K,LM,3) *RHO(K,1,3) +RHO(K,LM,4) *RHO(K,1,4)) /NORM/NORM
            IF (RHOCOL(K,1,2)<MAGMIN) THEN
               MLMPROJ =0
            ENDIF
            POT(K,1,2)  = POT(K,1,2)+DELTAV*(RHO(K,LM,2) - MLMPROJ*RHO(K,1,2) )/ NORM
            POT(K,1,3)  = POT(K,1,3)+DELTAV*(RHO(K,LM,3) - MLMPROJ*RHO(K,1,3) )/ NORM
            POT(K,1,4)  = POT(K,1,4)+DELTAV*(RHO(K,LM,4) - MLMPROJ*RHO(K,1,4) )/ NORM
         ENDDO
         ENDDO
         ENDDO

      ELSE
         DO L=1,LMAX
         DO M=0,2*L
         LM=L*L+M+1
         DO K=1,R%NMAX

            V0    =(POT(K,LM,1)+POT(K,LM,2))/2
            DELTAV=(POT(K,LM,1)-POT(K,LM,2))/2

            NORM     = MAX(ABS_AVMAG,MAGMIN)
            POT(K,LM,1) = V0
            POT(K,LM,2) = DELTAV * AVMAG(1) / NORM
            POT(K,LM,3) = DELTAV * AVMAG(2) / NORM
            POT(K,LM,4) = DELTAV * AVMAG(3) / NORM

            MLMPROJ = (RHO(K,LM,2) *AVMAG(1) +RHO(K,LM,3) *AVMAG(2) +RHO(K,LM,4) *AVMAG(3)) /NORM/NORM
            IF (ABS_AVMAG<MAGMIN) THEN
               MLMPROJ =0
            ENDIF
            POT(K,1,2)  = POT(K,1,2)+DELTAV*(RHO(K,LM,2) - MLMPROJ*AVMAG(1) )/ NORM*R%SI(K)
            POT(K,1,3)  = POT(K,1,3)+DELTAV*(RHO(K,LM,3) - MLMPROJ*AVMAG(2) )/ NORM*R%SI(K)
            POT(K,1,4)  = POT(K,1,4)+DELTAV*(RHO(K,LM,4) - MLMPROJ*AVMAG(3) )/ NORM*R%SI(K)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE RAD_MAG_DIRECTION


      SUBROUTINE RAD_MAG_DIRECTION_KINDENS( RHO, KINDENS, LMAX_TAU, POT, MU, R)
      USE constant
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      REAL(q) :: RHO(:,:,:)        !< density
      REAL(q) :: KINDENS(:,:,:)    !< kinetic energy density (total,magnetization)
      REAL(q) :: POT(:,:,:)        !< potential (total,magnetization)
      REAL(q) :: MU(:,:,:)         !< \mu_xc= d E_xc / d \tau, on entry (up,down), on return: (total,magnetization)
      TYPE (rgrid) :: R
      INTEGER LMAX_TAU
      ! local variables
      INTEGER ISPIN,I,L,M,LM,K
      REAL(q) NORM,DELTAV,V0,MLMPROJ
      REAL(q) AVMAG(3),ABS_AVMAG

      AVMAG=0
      DO K=1,R%NMAX
         AVMAG(1)=AVMAG(1)+RHO(K,1,2)*R%SI(K)
         AVMAG(2)=AVMAG(2)+RHO(K,1,3)*R%SI(K)
         AVMAG(3)=AVMAG(3)+RHO(K,1,4)*R%SI(K)
      ENDDO
      ABS_AVMAG=SQRT(AVMAG(1)**2+AVMAG(2)**2+AVMAG(3)**2)

      IF (.NOT. USE_AVERAGE_MAGNETISATION) THEN
         CALL vtutor%error("ERROR: RAD_MAG_DIRECTION_KINDENS: USE_AVERAGE_MAGNETISATION=.FALSE. is not &
            &implemented")
      ELSE
         DO L=0,LMAX_TAU
         DO M=0,2*L
         LM=L*L+M+1
         DO K=1,R%NMAX

            V0    =(MU(K,LM,1)+MU(K,LM,2))/2
            DELTAV=(MU(K,LM,1)-MU(K,LM,2))/2

            NORM=MAX(ABS_AVMAG,MAGMIN)

            MU(K,LM,1) = V0
            MU(K,LM,2) = DELTAV * AVMAG(1) / NORM
            MU(K,LM,3) = DELTAV * AVMAG(2) / NORM
            MU(K,LM,4) = DELTAV * AVMAG(3) / NORM

            MLMPROJ = (KINDENS(K,LM,2) *AVMAG(1) +KINDENS(K,LM,3) *AVMAG(2) +KINDENS(K,LM,4) *AVMAG(3)) /NORM/NORM
!           IF (ABS_AVMAG<MAGMIN) MLMPROJ=0
            IF (ABS_AVMAG<MAGMIN) THEN
               POT(K,1,2) = POT(K,1,2)+DELTAV*(KINDENS(K,LM,2)-MLMPROJ*AVMAG(1) )/NORM/HSQDTM*R%SI(K)!*2._q
               POT(K,1,3) = POT(K,1,3)+DELTAV*(KINDENS(K,LM,3)-MLMPROJ*AVMAG(2) )/NORM/HSQDTM*R%SI(K)!*2._q
               POT(K,1,4) = POT(K,1,4)+DELTAV*(KINDENS(K,LM,4)-MLMPROJ*AVMAG(3) )/NORM/HSQDTM*R%SI(K)!*2._q
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE RAD_MAG_DIRECTION_KINDENS


!*******************************************************************
!
!  RAD_CORE_XC
!
!>  Calculate the exchange correlation energy of the core charge
!>  density on the radial grid.
!>
!>  This energy will be subtracted (definining the zero of energy).
!>  This has the advantage that small differences in the potentials
!>  in the PP-generation program and VASP will not affect the
!>  energy.
!
!*******************************************************************

    SUBROUTINE RAD_CORE_XC( R, XC, RHOC, DEXC)
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      USE xc_driver, ONLY : GGA_DRIVER
      USE ldalib, ONLY : EXCOR_PARA
      IMPLICIT NONE

      REAL(q) :: RHOC(:)        !< core charge for exchange correlation
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
    ! local variables
      REAL(q) DRH(R%NMAX),ADRH(R%NMAX),DADRH(R%NMAX),DDRH(R%NMAX)
      REAL(q) T1(R%NMAX)
      REAL(q) RHOT(R%NMAX)
      INTEGER K,N
      REAL(q) SCALE, RS, DEXC, DEXC_GGA, VXT, EXT
      REAL(q) DEXC1,DVXC1,EXC,DEXC_
      LOGICAL,PARAMETER :: TREL=.TRUE.
      REAL(q) :: SIM_FAKT

      SCALE=2*SQRT(PI)
      N=R%NMAX

      DEXC    =0

      DO K=1,N
         RHOT(K) =RHOC(K)/ (SCALE*R%R(K)*R%R(K))  ! charge density
      ENDDO
     
      DO K=1,N
         IF (RHOT(K) /=0 ) THEN
            CALL  EXCOR_PARA(XC,RHOT(K),EXC,DEXC_,TREL)
            EXT = EXC/ RHOT(K)
         ELSE
            EXT=0
         ENDIF

         SIM_FAKT=R%SI(K)*SCALE
         DEXC    =DEXC    +EXT  *RHOC(K)*SIM_FAKT
      ENDDO

! GGA part
      DEXC_GGA=0

      IF (XC%LUSE_LIBXC.OR.XC%LDOGGA.OR.XC%LDOMETAGGA) THEN
         DEXC_GGA=0
         CALL GRAD(R,RHOT,T1)

         DO K=1,N
         IF (RHOT(K) > 1E-10_q) THEN
#ifdef vector
            CALL GGA91_WB(XC,RHOT(K)*AUTOA3,T1(K)*AUTOA4,EXT,DEXC1,DVXC1)
#else
            CALL GGA_DRIVER(XC,RHOT(K)*AUTOA3,T1(K)*AUTOA4,EXT,DEXC1,DVXC1,.FALSE.)
#endif
            SIM_FAKT=R%SI(K)*SCALE
            DEXC_GGA=DEXC_GGA+EXT*RHOC(K)*RYTOEV*SIM_FAKT
         ENDIF
         ENDDO
      ENDIF
#ifdef debug
      WRITE(0,1) DEXC/RYTOEV, DEXC_GGA/RYTOEV,DEXC/RYTOEV+DEXC_GGA/RYTOEV
 1    FORMAT(' core xc energy',3F14.7)
#endif
#undef debug
      DEXC=DEXC+DEXC_GGA
    END SUBROUTINE RAD_CORE_XC


!*******************************************************************
!
!  RAD_METRIC
!
!>  Evaluate the integral
!>~~~
!>   (phi^AE_i * phi^AE_j)^2 - ( phi^PS_i * phi^PS_j + Q_ij^L )^2
!>~~~
!>  This corresponds to the metric
!>~~~
!>  int_rad  | rho_AE(r)-rho_PS(r) | ^ 2 dr
!>~~~
!*******************************************************************

   SUBROUTINE RAD_METRIC( R, AUG, QPAW, WAE1, WAE2, WPS1, WPS2, AMETRIC)
      USE constant
      IMPLICIT NONE

      TYPE (rgrid) R
      REAL(q) :: AUG(:)          !< 1-normalized compensation charge
      REAL(q) :: WAE1(:),WAE2(:) !< all electron wavefunction
      REAL(q) :: WPS1(:),WPS2(:) !< pseudo wavefunctions
      REAL(q) :: QPAW            !< moments of compensation charge
      REAL(q) AMETRIC
! local variables
      REAL(q) :: RHOT(R%NMAX)
      INTEGER I


      DO I=1,R%NMAX
         RHOT(I)=(WAE1(I)*WAE2(I))**2-(WPS1(I)*WPS2(I)+AUG(I)*QPAW)**2
      ENDDO
      CALL SIMPI(R, RHOT, AMETRIC)
      AMETRIC=AMETRIC*2*SQRT(PI)

    END SUBROUTINE RAD_METRIC


!*******************************************************************
!
!> Calculate the gradient of a function on the radial grid
!> nabla rho(r) = d/dr rho(r)
!>
!> Adopted from some routine of a program called atoms.
!
!*******************************************************************

      SUBROUTINE GRAD_(R,RH,DRH)
      IMPLICIT NONE
      TYPE (rgrid) R      !< defines the radial grid
      REAL(q) RH(*)       !< charge density
      REAL(q) DRH(*)      !< gradient of the function
      REAL(q) H
      INTEGER NR,NM2,I

      H=R%H
      NR=R%NMAX
      NM2=NR-2
! For the first and second point the forward differences are used
      DRH(1)=((6._q*RH(2)+20._q/3._q*RH(4)+1.2_q*RH(6)) &
     &      -(2.45_q*RH(1)+7.5_q*RH(3)+3.75_q*RH(5)+1._q/6._q*RH(7)))/H
      DRH(2)=((6._q*RH(3)+20._q/3._q*RH(5)+1.2_q*RH(7)) &
     &        -(2.45_q*RH(2)+7.5_q*RH(4)+3.75_q*RH(6)+1._q/6._q*RH(8)))/H
! Five points formula
      DO  I=3,NM2
         DRH(I)=((RH(I-2)+8._q*RH(I+1))-(8._q*RH(I-1)+RH(I+2)))/12._q/H
      ENDDO
! Five points formula for the last two points ('backward differences')
      DRH(NR-1)=(-1._q/12._q*RH(NR-4)+0.5_q*RH(NR-3)-1.5_q*RH(NR-2) &
     &           +5._q/6._q*RH(NR-1)+0.25_q*RH(NR))/H
      DRH(NR)=  (0.25_q*RH(NR-4)-4._q/3._q*RH(NR-3)+3._q*RH(NR-2) &
     &           -4._q*RH(NR-1)+25._q/12._q*RH(NR))/H
! account for logarithmic mesh
      DO  I=1,NR
         DRH(I)=DRH(I)/R%R(I)
      ENDDO

      RETURN
      END SUBROUTINE GRAD_

!*******************************************************************
!
!> Calculate the gradient of a function on the radial grid
!> nabla rho(r) = d/dr rho(r)
!> using different orders of the finite difference stencil
!> choosen from the INCAR
!
!*******************************************************************
      SUBROUTINE GRAD(R,RH,DRH)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE (rgrid) R      !< defines the radial grid
      REAL(q) RH(*)       !< charge density
      REAL(q) DRH(*)      !< gradient of the function
      INTEGER,PARAMETER :: ORDER=4
      ! local variables
      REAL(q) H
      INTEGER NR,I

      H=R%H
      NR=R%NMAX

      SELECT CASE(R%FDORDER)
         CASE(2)
            ! Calculate the gradient of a function on the radial grid
            ! using 2.order differentiation (this is good enough
            ! and less suseptible to noise than the previous routine -> here implemented below (option 3)
            DRH(1) = RH(2)-RH(1)
            DO I = 2, NR - 1
               DRH(I) = (RH(I+1)-RH(I-1))/2.0_q
            END DO
            DRH(NR) = RH(NR)-RH(NR - 1)
         CASE(3)
            ! Adopted from some routine of a program called atoms.
            ! Strange choice but we keep it for completeness:
            !  - 6th order accuracy for forward difference in the first 2 points
            !  - 4th order accuracy for the midpoints and backwards finite difference on the endpoints
            ! Perhaps the first points are more sensitive, hence the higher order

            ! For the first and second point the forward differences are used
            DRH(1)=((6._q*RH(2)+20._q/3._q*RH(4)+1.2_q*RH(6)) &
           &      -(2.45_q*RH(1)+7.5_q*RH(3)+3.75_q*RH(5)+1._q/6._q*RH(7)))
            DRH(2)=((6._q*RH(3)+20._q/3._q*RH(5)+1.2_q*RH(7)) &
           &        -(2.45_q*RH(2)+7.5_q*RH(4)+3.75_q*RH(6)+1._q/6._q*RH(8)))
            ! Five points formula
            DO  I=3,NR - 2
               DRH(I)=((RH(I-2)+8._q*RH(I+1))-(8._q*RH(I-1)+RH(I+2)))/12._q
            ENDDO
            ! Five points formula for the last two points ('backward differences')
            DRH(NR-1)=(-1._q/12._q*RH(NR-4)+0.5_q*RH(NR-3)-1.5_q*RH(NR-2) &
           &           +5._q/6._q*RH(NR-1)+0.25_q*RH(NR))
            DRH(NR)=  (0.25_q*RH(NR-4)-4._q/3._q*RH(NR-3)+3._q*RH(NR-2) &
           &           -4._q*RH(NR-1)+25._q/12._q*RH(NR))
         CASE(4)
            DRH(1) = -25.0d0/12.0d0*RH(1) + 4*RH(2) - 3*RH(3) + (4.0d0/3.0d0)*RH(4) - 1.0d0/4.0d0*RH(5)
            DRH(2) = -25.0d0/12.0d0*RH(2) + 4*RH(3) - 3*RH(4) + (4.0d0/3.0d0)*RH(5) - 1.0d0/4.0d0*RH(6)
            DO I = 3, NR - 2
               DRH(I) = (2.0d0/3.0d0)*RH(I + 1) - 1.0d0/12.0d0*RH(I + 2) - 2.0d0/3.0d0*RH(I - 1) + (1.0d0/12.0d0)*RH(I - 2)
            END DO
            DRH(NR)     = (1.0d0/3.0d0)*RH(NR + 1) - RH(NR - 1) + (1.0d0/6.0d0)*RH(NR - 2) + (1.0d0/2.0d0)*RH(NR)
            DRH(NR - 1) = (1.0d0/2.0d0)*RH(NR - 1) - RH(NR - 2) + (1.0d0/6.0d0)*RH( NR - 3) + (1.0d0/3.0d0)*RH(NR)
         CASE(6)
            DRH(1) = -49.0_q/20.0_q*RH(1) + 6*RH(2) - 15.0_q/2.0_q*RH(3) + (20.0_q/3.0_q)*RH(4) - 15.0_q/4.0_q*RH(5) + (6.0_q/5.0_q)*RH(6) - 1.0_q/6.0_q*RH(7)
            DRH(2) = -49.0_q/20.0_q*RH(2) + 6*RH(3) - 15.0_q/2.0_q*RH(4) + (20.0_q/3.0_q)*RH(5) - 15.0_q/4.0_q*RH(6) + (6.0_q/5.0_q)*RH(7) - 1.0_q/6.0_q*RH(8)
            DRH(3) = -49.0_q/20.0_q*RH(3) + 6*RH(4) - 15.0_q/2.0_q*RH(5) + (20.0_q/3.0_q)*RH(6) - 15.0_q/4.0_q*RH(7) + (6.0_q/5.0_q)*RH(8) - 1.0_q/6.0_q*RH(9)
            DO I = 4, NR - 3
               DRH(I) = (3.0_q/4.0_q)*RH(I + 1) - 3.0_q/20.0_q*RH(I + 2) + (1.0_q/60.0_q)*RH(I + 3) - 3.0_q/4.0_q*RH(I - 1) + (3.0_q/20.0_q)*RH(I - 2) - 1.0_q/60.0_q*RH(I - 3)
            END DO
            DRH(NR) = (1.0_q/5.0_q)*RH(NR + 1) - 2*RH(NR - 1) + RH(NR - 2) - 1.0_q/3.0_q*RH(NR - 3) + (1.0_q/20.0_q)*RH(NR - 4) + (13.0_q/12.0_q)*RH(NR)
            DRH(NR - 1) = (13.0_q/12.0_q)*RH(NR - 1) - 2*RH(NR - 2) + RH(NR - 3) - 1.0_q/3.0_q*RH(NR - 4) + (1.0_q/20.0_q)*RH(NR - 5) + (1.0_q/5.0_q)*RH(NR)
            DRH(NR - 2) = (1.0_q/5.0_q)*RH(NR - 1) + (13.0_q/12.0_q)*RH(NR - 2) - 2*RH(NR - 3) + RH(NR - 4) - 1.0_q/3.0_q*RH(NR - 5) + (1.0_q/20.0_q)*RH(NR - 6)
         CASE DEFAULT
            CALL vtutor%error('Invalid order for the finite differences stencil RADIAL_FDORDER. Must be 2, 4 or 6, found '//str(R%FDORDER))
      END SELECT

! account for logarithmic mesh
      DO  I=1,NR
         DRH(I)=DRH(I)/R%R(I)/H
      ENDDO

      END SUBROUTINE GRAD


!*******************************************************************
!
!> Calculate the gradient of a function on the radial grid
!> using 2.order differentiation (this is good enough
!> and less suseptible to noise than the previous routine)
!> nabla rho(r) = d/dr rho(r)
!
!*******************************************************************

      SUBROUTINE GRAD_OLD(R,RH,DRH)
      IMPLICIT NONE
      TYPE (rgrid) R      !< defines the radial grid
      REAL(q) RH(*)       !< charge density
      REAL(q) DRH(*)      !< gradient of the function
      REAL(q) H
      INTEGER NR,NM1,I

      H=R%H
      NR=R%NMAX
      NM1=NR-1
! 1. point use first order differantion
      DRH(1)=(RH(2)-RH(1))/H
! three point formula
      DO  I=2,NM1
         DRH(I)=(RH(I+1)-RH(I-1))/2/H
      ENDDO
! last point
      DRH(NR)=(RH(NR)-RH(NR-1))/H

      DO  I=1,NR
         DRH(I)=DRH(I)/R%R(I)
      ENDDO

      RETURN
      END SUBROUTINE GRAD_OLD

!*******************************************************************
!
!  AUG_SETE
!
!>  Calculate coefficients for generalized gaussian
!>   exp ( r^2 / alpha).
!>  Not used.
!*******************************************************************


    SUBROUTINE AUG_SETE(R, ALPHA, A, TH)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (rgrid) R

      REAL(q)  ALPHA, A, TH
      REAL(q)  TMP(R%NMAX)

      RC= R%RMAX
      ALPHA= LOG( TH) / (RC*RC)

      DO N=1,R%NMAX
         TMP(N)=0
         IF (R%R(N) <= RC) THEN
            VAL=EXP(ALPHA*R%R(N)*R%R(N))
            TMP(N)=R%R(N)*R%R(N)*VAL
         ENDIF
      ENDDO
      CALL SIMPI(R,TMP,SUM)
      A = 1/SUM

    END SUBROUTINE AUG_SETE

!*******************************************************************
!
!  AUG_SETQ
!
!>  Find a set of 2 q_i and coefficients A_i such that
!>~~~
!>     sum_i  A_i j_l(q_i,Rc)  =0
!>     sum_i  A_i j_l(q_i,Rc)' =0'
!>     sum_i  A_i j_l(q_i,Rc)''=0
!>~~~
!>  and
!>~~~
!>     sum_i int_0^Rc A_i j_l(q_i r) r^(l+2) dr = 1
!>~~~
!>  the augmentation charge is then represented by
!>~~~
!>     sum_i  A_i j_l(q_i, r)
!>~~~
!  PS: I confess I have an obsession with spherical Besselfunctions
!*******************************************************************


    SUBROUTINE AUG_SETQ(L,R,RMAX,QQ,A,LCOMPAT)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (rgrid) R

      PARAMETER (NQ=2)
      REAL(q)  QQ(NQ)
      REAL(q)  A(NQ)
      LOGICAL LCOMPAT
      REAL(q)  AMAT(NQ,NQ)
      REAL(q)  B(3)
      INTEGER  IPIV(NQ)

      REAL(q)  TMP(R%NMAX)
!-----gaussian integration
      PARAMETER (M=32)
      REAL(q)  WR(M),RA(M)
      INTEGER IFAIL
      EXTERNAL GAUSSI2
!-----------------------------------------------------------------------
! search for q so that j(q Rc)=0
!-----------------------------------------------------------------------
      CALL AUG_BEZERO(QQ,L,NQ)
!-----------------------------------------------------------------------
! second set the matrix
!-----------------------------------------------------------------------
      ITYPE=0
      CALL GAUSSI(GAUSSI2,0.0_q,RMAX, ITYPE,M,WR,RA,IFAIL)

      DO I=1,NQ
         QQ(I)=QQ(I)/RMAX
         QR  =QQ(I)*RMAX
         CALL SBESSE3( QR, BJ, BJP, BJPP, L)
         BJP =BJP *QQ(I)
         BJPP=BJPP*QQ(I)*QQ(I)
         AMAT(1,I) = BJP
!  I could not find an analytical expression for \int j_l(qr) r^(2+l) dr
!  I will check again probably it is simple enough
         IF (LCOMPAT) THEN
           DO N=1,R%NMAX
              IF (R%R(N) < RMAX) THEN
                 QR=QQ(I)*R%R(N)
                 CALL SBESSEL(QR,BJ,L)
                 TMP(N)=BJ*R%R(N)**(2+L)
              ELSE
                 TMP(N)=0
              ENDIF
           ENDDO
           CALL SIMPI(R,TMP,SUM)
           AMAT(2,I) = SUM
         ELSE

! Gauss integration, more accurate !
           SUM2=0
           DO N=1,M
              QR=QQ(I)*RA(N)
              CALL SBESSEL(QR,BJ,L)
              SUM2=SUM2+BJ*RA(N)**(2+L)*WR(N)
           ENDDO

           AMAT(2,I) = SUM2
         ENDIF
      ENDDO
!-----------------------------------------------------------------------
!  solve the linear equations  B_i = AMAT_ii' A_i'
!-----------------------------------------------------------------------
      IFAIL=0
      B(1)=0
      B(2)=1

      A(1)=0
      A(2)=0

      CALL DGETRF( NQ, NQ, AMAT, NQ, IPIV, IFAIL )
      CALL DGETRS('N', NQ, 1, AMAT, NQ, IPIV, B, NQ, IFAIL)
      A=B(1:2)

      B(1)=0
      B(2)=0
      B(3)=0

      DO I=1,NQ
         SUM=0
         QR  =QQ(I)*RMAX
         CALL SBESSE3( QR, BJ, BJP, BJPP, L)
         B(1)=B(1)+BJ*A(I)
         B(2)=B(2)+BJP*A(I)*QQ(I)
         B(3)=B(3)+BJPP*A(I)*QQ(I)*QQ(I)
      ENDDO
    END SUBROUTINE AUG_SETQ

!*******************************************************************
!  SUBROUTINE BEZERO
!
!>  Searches for NQ zeros j(qr).
!
!>  i/o: \n
!>         XNULL(NQ) result \n
!>         L           quantum number l
!  great-full spaghetti code (written by gK)
!********************************************************************

    SUBROUTINE AUG_BEZERO(XNULL,L,NQ)
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (STEP=.1_q, BREAK= 1E-10_q)
      DIMENSION XNULL(NQ)
! initialization
      X=STEP
      N=0
! entry point for next q_n
  30  CALL SBESSE2(X, BJ1, DUMMY,  L)
! coarse search
  10  X=X+STEP
      CALL SBESSE2(X, BJ2, DUMMY,  L)
! found one point
      IF(BJ1*BJ2 < 0) THEN
        ETA=0.0_q
! interval bisectioning
        SSTEP=STEP
        XX   =X
  20    SSTEP=SSTEP/2
        IF (BJ1*BJ2 < 0) THEN
          XX=XX-SSTEP
        ELSE
          XX=XX+SSTEP
        ENDIF
        CALL SBESSE2(XX, BJ2, DUMMY,  L)
        IF (SSTEP > BREAK) GOTO 20

        N=N+1
        XNULL(N)=XX
        IF (N == NQ) RETURN
        GOTO 30
      ENDIF
      GOTO 10

    END SUBROUTINE AUG_BEZERO

!********************************************************************
!
!  SUBROUTINE SIMPI
!
!>  Integrate a function on the logarithmic grid.
!>  Uses the previously setup weights.
!
!********************************************************************

    SUBROUTINE SIMPI(R,F,FI)
      IMPLICIT NONE
      TYPE (rgrid) R
      REAL(q)  F(:),FI,SUM
      INTEGER   K

      SUM=0
!OCL SCALAR
      DO K=1,R%NMAX
         SUM=SUM+F(K)*R%SI(K)
      ENDDO

      FI=SUM

    END SUBROUTINE SIMPI

!********************************************************************
!
!  SUBROUTINE SET_SIMP
!
!>  Setup weights for simpson integration on radial grid.
!>
!>  Any radial integral can then be evaluated by just summing all
!>  radial grid points with the weights SI
!>~~~
!>  int dr = sum_i si(i) * f(i)
!>~~~
!>  The factors  R%R(K)  *R%H stem from the logarithmic grid.
!********************************************************************

    SUBROUTINE SET_SIMP(R)
      IMPLICIT NONE
      TYPE (rgrid) R
      INTEGER   K

      ALLOCATE(R%SI(R%NMAX))
      R%SI=0
      DO K=R%NMAX,3,-2
         R%SI(K)=    R%R(K)  *R%H/3._q+R%SI(K)
         R%SI(K-1)=4*R%R(K-1)*R%H/3._q
         R%SI(K-2)=  R%R(K-2)*R%H/3._q
      ENDDO
    END SUBROUTINE SET_SIMP

!******************************************************************** 
!
!  SPHERICAL_BESSEL_TRANSFORM_R2Q
!
!********************************************************************
      SUBROUTINE SPHERICAL_BESSEL_TRANSFORM_R2Q( &
     & L,R,FR,G,FG &
     & )
      USE constant
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(rgrid) R
      INTEGER L
      REAL(q) FR(:),G(:),FG(:)
      ! local variables
      INTEGER IG,IR
      REAL(q) QR,BJ,SUM,SCALE
      REAL(q), ALLOCATABLE :: TMP(:),SI(:)
      
      ! quick consistency checks
      IF (SIZE(FR)<R%NMAX.OR.SIZE(FG)<SIZE(G)) THEN
         CALL vtutor%error("SPHERICAL_BESSEL_TRANSFORM: ERROR: array(s) not sufficiently large: " // &
            str(SIZE(FR)) // " " // str(R%NMAX) // " " // str(SIZE(FG)) // " " // str(SIZE(G)))
      ENDIF

!     SCALE=SQRT(2._q/PI)
      
      ALLOCATE(TMP(R%NMAX))

      DO IG=1,SIZE(G)
         TMP=0
         DO IR=1,R%NMAX
            QR=G(IG)*R%R(IR)
            CALL SBESSEL(QR,BJ,L)
!           TMP(IR)=FR(IR)*QR*BJ
            TMP(IR)=FR(IR)*BJ*R%R(IR)*R%R(IR)
         ENDDO
         CALL SIMPI(R,TMP,SUM)
!        FG(IG)=SCALE*SUM
         FG(IG)=(4._q*PI)*SUM
      ENDDO  

      RETURN
      END SUBROUTINE SPHERICAL_BESSEL_TRANSFORM_R2Q

!********************************************************************
!
!  BESSEL_TRANSFORM_RADIAL_FUNCTION
!
!********************************************************************
      SUBROUTINE BESSEL_TRANSFORM_RADIAL_FUNCTION( &
     & L,R,FR,DG,FG &
     &)
      IMPLICIT NONE
      TYPE(rgrid) R
      INTEGER L
      REAL(q) FR(:)
      REAL(q) DG
      REAL(q) FG(:,:)
      ! local variables
      INTEGER IG

      DO IG=1,SIZE(FG,1)
         FG(IG,1)=(IG-1)*DG 
      ENDDO
      CALL SPHERICAL_BESSEL_TRANSFORM_R2Q(L,R,FR(:),FG(:,1),FG(:,2))

      CALL SPLCOF(FG,SIZE(FG,1),SIZE(FG,1),0._q) 

      RETURN
      END SUBROUTINE BESSEL_TRANSFORM_RADIAL_FUNCTION


!******************** SUBROUTINE SET_PAW_METAGGA ***********************
!
!***********************************************************************
      SUBROUTINE SET_PAW_METAGGA(P)
      USE prec
      USE pseudo
      USE constant
      IMPLICIT NONE
      TYPE (potcar), POINTER :: P
      ! local variables
      INTEGER I,NMIN,NMAX
      REAL(q) SCALE
      REAL(q) RHO(P%R%NMAX),DRHO(P%R%NMAX)

      SCALE=1/(2*SQRT(PI)) ! Y00

      NMAX=P%R%NMAX

!     NMIN=NMAX
!     DO I=1,P%LMAX
!        NMIN=MIN(NMIN,IR_MATCH(P,I))
!     ENDDO
      NMIN=IR_PSMAX(P)

      RHO(:)=SCALE*P%RHOAE(:)/(P%R%R(:)*P%R%R(:))

      CALL GRAD(P%R,RHO,DRHO)
      DRHO=ABS(DRHO)
      DO I=1,NMAX
         IF ((DRHO(I)*DRHO(I)/RHO(I)/8._q*2._q*HSQDTM) > (P%TAUAE(I)*SCALE/(P%R%R(I)*P%R%R(I))) ) &
        &   P%TAUAE(I)=(DRHO(I)*DRHO(I)/RHO(I)/8._q*2._q*HSQDTM)/SCALE*(P%R%R(I)*P%R%R(I)) !*1.00001
      ENDDO

!     P%TAUAE(NMIN:NMAX)=0
      P%RHOAE(NMIN:NMAX)=P%RHOPS(NMIN:NMAX)

      RETURN
      END SUBROUTINE SET_PAW_METAGGA


!******************** SUBROUTINE RAD_KINEDEN ***************************
!
!***********************************************************************
      SUBROUTINE RAD_KINEDEN(&
     &   P,XC,W,LMAX_TAU,NCDIJ,COCC,TAU)
      USE prec
      USE asa
      USE constant
      USE pseudo
      USE setexm_struct_def, ONLY : LMAXTAU,xc_info
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE (potcar) P
      TYPE (xc_info) XC
      REAL(q) W(:,:)
      INTEGER LMAX_TAU,NCDIJ
      OVERLAP COCC(:,:,:)
      REAL(q) TAU(:,:,:)
      ! local variables
      INTEGER ISP,I,Ch1,Ch2
      INTEGER LLMAX_TAU,LMMAX_TAU,LM,LMP,LL,LLP,M,MP,MPLOW,INDYLM,INDYLMP

      REAL(q) FAKT
      REAL(q) W1(P%R%NMAX),W2(P%R%NMAX),DW1(P%R%NMAX),DW2(P%R%NMAX)

      TAU=0

      ! quick return if possible
      IF (.NOT.XC%LTAU) RETURN

      PROFILING_START("rad_kineden")

      ! check consistency
      LMMAX_TAU=(LMAX_TAU+1)**2
      IF (SIZE(TAU,2)<LMMAX_TAU) THEN
         CALL vtutor%error("RAD_KINEDEN: ERROR, size(TAU,2)<LMMAX_TAU " // str(SIZE(TAU,2)) // " " // &
            str(LMMAX_TAU))
      ENDIF

      ! restriction for the angular development of \tau and \mu
      LLMAX_TAU=MIN(LMAXTAU,LMAX_TAU); LMMAX_TAU=(LLMAX_TAU+1)**2

      spin: DO ISP=1,NCDIJ

         ! loops over channels (l,epsilon)
         LM=1; DO Ch1=1,P%LMAX
         ! d(w1/r)/dr
         W1(:)=W(:,Ch1)/P%R%R(:); CALL GRAD(P%R,W1,DW1)

         LMP=LM; DO Ch2=Ch1,P%LMAX
         ! d(w2/r)/dr
         W2(:)=W(:,Ch2)/P%R%R(:); CALL GRAD(P%R,W2,DW2)

            LL=P%LPS(Ch1); LLP=P%LPS(Ch2)
            DO M=1,2*LL+1
            MPLOW=1; IF(Ch1==Ch2) MPLOW=M
            DO MP=MPLOW,2*LLP+1
               INDYLM= LL**2 +M
               INDYLMP=LLP**2+MP

               FAKT=REAL(COCC(LM+M-1,LMP+MP-1,ISP)+COCC(LMP+MP-1,LM+M-1,ISP),q)
               IF ((LM+M-1)==(LMP+MP-1)) THEN
                  FAKT=REAL(COCC(LM+M-1,LMP+MP-1,ISP),q)
               ENDIF

               DO I=1,LMMAX_TAU
                  TAU(:P%R%NMAX,I,ISP)=TAU(:P%R%NMAX,I,ISP)+FAKT*&
                 &   (DW1(:)*X_YLM_X_YLM(INDYLM,INDYLMP,I)*DW2(:)*P%R%R(:)*P%R%R(:) + &
                 &    DW1(:)*X_YLM_YLMD(INDYLM,INDYLMP,I)*W2(:)*P%R%R(:) + &
                 &    W1(:)*X_YLM_YLMD(INDYLMP,INDYLM,I)*DW2(:)*P%R%R(:) + &
                 &    W1(:)*YLMD_YLMD(INDYLM,INDYLMP,I)*W2(:))
               ENDDO

            ENDDO
            ENDDO
         LMP=LMP+2*LLP+1
         ENDDO
         LM =LM +2*LL +1
         ENDDO

      ENDDO spin

      TAU=TAU*HSQDTM

      PROFILING_STOP("rad_kineden")

      RETURN
      END SUBROUTINE RAD_KINEDEN


!********************* SUBROUTINE RAD_POT_METAGGA **********************
!
!***********************************************************************
      SUBROUTINE RAD_POT_METAGGA(&
     &   R,XC,ISPIN,LMAX_CALC,LMAX_TAU,LASPH,CVMBJ,RHO,RHOC,POTC,TAU,POT,DOUBLEC,EXCG,GRDR,MU,NI,AUGFACT,TAUC)
      USE prec
      USE constant
      USE setexm_struct_def, ONLY : xc_info
      IMPLICIT NONE
      TYPE (rgrid) R
      TYPE (xc_info) XC
      INTEGER ISPIN
      INTEGER LMAX_CALC,LMAX_TAU
      REAL(q) RHO(:,:,:)     ! charge distribution (charge,magnetization)
      REAL(q) RHOC(:)        ! core charge distribution
      REAL(q) POTC(:)        ! frozen core potential (Hartree)
      REAL(q) TAU(:,:,:)     ! kinetic energy density (up,down)
      REAL(q) POT(:,:,:)     ! Total potential (up,down)
      REAL(q) CVMBJ          ! MBJ C-parameter
      REAL(q) GRDR           ! integral of |nabla rho(r)|/ rho(r) over the PAW sphere
      REAL(q) EXCG           ! exchange energy only
      REAL(q) DOUBLEC        ! double counting corrections
      REAL(q) MU(:,:,:)      ! d e_xc/d tau (up,down)
      REAL(q), INTENT(IN) :: AUGFACT
      REAL(q), OPTIONAL :: TAUC(:) ! kinetic energy density of the core electrons
      LOGICAL LASPH
      INTEGER NI
      ! local variables
      INTEGER N,L,M,LM,K
      REAL(q) SCALE
      REAL(q) SUM,DHARTREE
      REAL(q) DEXC,DVXC

      MU=0

      ! quick return if possible
      IF (.NOT.XC%LDOMETAGGA) RETURN

      SCALE=2*SQRT(PI)
      N=R%NMAX

      ! Exchange-correlation
      DEXC=0
      DVXC=0

      IF (PRESENT(TAUC)) THEN
         CALL RAD_METAGGA_XC(R,XC,ISPIN,LMAX_CALC,LMAX_TAU,LASPH,CVMBJ,RHO,RHOC,POTC,TAU,DEXC,DVXC,GRDR,POT,MU,NI,AUGFACT,TAUC)
      ELSE
         CALL RAD_METAGGA_XC(R,XC,ISPIN,LMAX_CALC,LMAX_TAU,LASPH,CVMBJ,RHO,RHOC,POTC,TAU,DEXC,DVXC,GRDR,POT,MU,NI,AUGFACT)
      ENDIF

      EXCG=EXCG+DEXC
      DOUBLEC=DOUBLEC+DEXC-DVXC
#ifdef debug
      WRITE(*,1) DHARTREE*0.5_q+DEXC, DOUBLEC
1     FORMAT('  Hartree + exc, -Hartree + exc -vxc',6F14.7)
#endif

      RETURN
      END SUBROUTINE RAD_POT_METAGGA


!********************* SUBROUTINE RAD_METAGGA_XC ***********************
!
!***********************************************************************
      SUBROUTINE RAD_METAGGA_XC(&
     &   R,XC,ISPIN,LMAX_CALC,LMAX_TAU,LASPH,CVMBJ,RHO,RHOC,POTC,TAU,DEXC,DVXC,GRDR,POT,MU,NI,AUGFACT,TAUC)

GPACC USE moffload_tutor, ONLY : VTUTOR_OFFLOAD_CHECK

      USE prec
      USE asa
      USE constant
      USE setexm_struct_def, ONLY : xc_info, PHPTS_FACT, THPTS_FACT, LMAXTAU
      USE xc_driver, ONLY : METAGGA_DRIVER
      USE mbj, ONLY : LscMBJ, GRHO_OVER_RHO_AUG
      USE xc_name

      IMPLICIT NONE
      TYPE (rgrid) :: R
      TYPE (xc_info) :: XC
      INTEGER ISPIN
      INTEGER LMAX_CALC,LMAX_TAU
      INTEGER ELL,LFAC,M
      REAL(q) RHOC(:)        ! core charge distribution
      REAL(q) POTC(:)        ! frozen core potential (Hartree)
      REAL(q) RHO(:,:,:)     ! charge distribution (charge,magnetization)
      REAL(q) TAU(:,:,:)     ! kinetic energy density (up,down)
      REAL(q) CVMBJ          ! MBJ C-parameter
      REAL(q) GRDR           ! integral of |nabla rho(r)|/ rho(r) over the PAW sphere
      REAL(q) DEXC, DVXC
      REAL(q) POT(:,:,:)     ! potential (up,down)
      REAL(q) MU(:,:,:)      ! d e_xc/d tau (up,down)
      REAL(q), INTENT(IN) :: AUGFACT
      REAL(q), OPTIONAL :: TAUC(:) ! kinetic energy density of the core electrons
      LOGICAL LASPH
      INTEGER NI
      ! local variables
      REAL(q) RHOT(R%NMAX,2),RHO_ANG(R%NMAX,2,3),T1(R%NMAX,2),T2(R%NMAX,2)
      REAL(q) TAUU(R%NMAX),TAUD(R%NMAX)
      INTEGER ISP,K,I,J,IFAIL
      INTEGER LM,LMP,LLMAX,LMMAX,LMMAX_alloc,LMMAX_CALC,LLMAX_TAU,LMMAX_TAU,LMAX_TMP
      INTEGER NP,PHPTS,THPTS,NPTS
      REAL(q) SCALE,SUM,DELTAPHI,EVTOH

      REAL(q) DVINI,DCMU

      REAL(q) SIM_FAKT
      REAL(q) TEMP1
      REAL(q), ALLOCATABLE :: RADPTS(:,:), XYZPTS(:,:)
      REAL(q), ALLOCATABLE :: YLM(:,:),YLMD(:,:,:),YLMDD(:,:,:)
      REAL(q), ALLOCATABLE :: YLM_NABLA_YLM(:,:,:),YLM_X_YLM(:,:,:)
      REAL(q), ALLOCATABLE :: dFXC(:,:,:,:),dFXCdT(:,:,:),dFXCdLap(:,:,:)
      REAL(q), ALLOCATABLE :: WEIGHT(:),ABSCIS(:)
      REAL(q), ALLOCATABLE :: RHOLM(:,:,:), RHOLMD(:,:,:)

      REAL(q) XU,YU,ZU,XD,YD,ZD
      REAL(q) TG,LAPLUP,LAPLDW,MUUP,MUDW,DLXC1,DLXC2
      REAL(q) EXT,DEXC1,DEXC2,DVXC1,DVXC2,DVC
      REAL(q) TMP(ISPIN,3)
      EXTERNAL GAUSSI2

      REAL(q) :: CHGMIN=1E-10_q, NABMIN=1E-10_q, TAUMIN=1E-10_q

!     TAUMIN=0 for MBJ is ok since there is no division by TAU
      IF ((XC%ID(1)==ID_XC_MBJ).OR.(XC%ID(1)==ID_XC_LMBJ)) TAUMIN=0._q
      ! initial double counting corrections
      DVINI=0
      CALL RAD_POT_IN_RHO(R,ISPIN,LMAX_CALC,RHO,POT,DVINI)

      SCALE=2*SQRT(PI) ! 1/Y00

      EVTOH=1._q/(2.*HSQDTM)*AUTOA5

      LMMAX_CALC=(LMAX_CALC+1)**2
      ! restriction for the angular development of \rho
      LLMAX=MIN(6,LMAX_CALC)
      ! in GGA this needs to be increased by 1 (nabla-comment)
      LLMAX=LLMAX+1; LMMAX=(LLMAX+1)**2

      ! restriction for the angular development of \tau and \mu
      LLMAX_TAU=MIN(LMAXTAU,LMAX_TAU); LMMAX_TAU=(LLMAX_TAU+1)**2

      ! restrict everything to spherical contributions only (LASPH=.FALSE.)
      IF (.NOT.LASPH) THEN
         LLMAX=1; LMMAX=4; LLMAX_TAU=0; LMMAX_TAU=1; LMMAX_CALC=1
      ENDIF

      LMMAX_alloc=(MAX(LLMAX,LLMAX_TAU)+1)**2

!========================================================================
! number of theta and phi pivot points to perform angular integration
! since Exc=f(a*Yllmax,m) we need more pivot points than theoretically
! needed to integrate Yllmax,m.
! the factor 2 is the minium, 3 is more accurate
!========================================================================
      PHPTS=PHPTS_FACT*(MAX(LLMAX,LLMAX_TAU)+1)
      THPTS=THPTS_FACT*FLOOR(REAL(MAX(LLMAX,LLMAX_TAU)/2+1,KIND=q))
      NPTS=PHPTS*THPTS
      DELTAPHI=REAL(2_q*PI/PHPTS,KIND=q)
! allocate arrays
      ALLOCATE(XYZPTS(NPTS,3),RADPTS(NPTS,2),WEIGHT(THPTS),ABSCIS(THPTS), &
     &     YLM(NPTS,LMMAX_alloc),YLMD(NPTS,LMMAX_alloc,3),YLMDD(NPTS,LMMAX_alloc,6), &
     &     YLM_NABLA_YLM(LMMAX,LMMAX,0:3),YLM_X_YLM(LMMAX,LMMAX,0:3), &
     &     RHOLM(R%NMAX,LMMAX,2),RHOLMD(R%NMAX,LMMAX,2))

      ALLOCATE(dFXC(R%NMAX,LMMAX,0:3,2)); dFXC=0
      CALL SETYLM_NABLA_YLM(LLMAX,YLM_NABLA_YLM,YLM_X_YLM)
      IF (XC%LDLAP) THEN
         ALLOCATE(dFXCdLap(R%NMAX,LMMAX,2)); dFXCdLap=0
      ENDIF
      IF (XC%LMU) THEN
         ALLOCATE(dFXCdT(R%NMAX,LMMAX_TAU,2)); dFXCdT=0
      ENDIF

      ! set phi positions, equally spaced
      RADPTS=0; WEIGHT=0; ABSCIS=0
      DO I=1,PHPTS
         DO J=1,THPTS
            RADPTS((J-1)*PHPTS+I,2)=(I-1)*DELTAPHI
         ENDDO
      ENDDO
      ! get theta positions (actually get cos(theta)) (Gauss integration)
      CALL GAUSSI(GAUSSI2,-1._q,1._q,0,THPTS,WEIGHT,ABSCIS,IFAIL)
      DO I=1,THPTS
         RADPTS((I-1)*PHPTS+1:I*PHPTS,1)=ABSCIS(I)
      ENDDO
      ! convert radial to cartesian coordinates
      DO I=1,NPTS
         XYZPTS(I,1)=COS(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! x
         XYZPTS(I,2)=SIN(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! y
         XYZPTS(I,3)=RADPTS(I,1)                                 ! z
      ENDDO

      ! get |r| Y_lm on a unit sphere and its derivatives
      YLM=0; YLMD=0; YLMDD=0
!     CALL SETYLM_GRAD2(LLMAX,NPTS,YLM,YLMD,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))
      CALL SETYLM_GRAD3(MAX(LLMAX,LLMAX_TAU),NPTS,YLM,YLMD,YLMDD,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))
!========================================================================
! main loop over all grid points
! prepare the charge density array RHOT
!========================================================================
      DEXC=0; GRDR=0

      ! loop over all points in the angular grid
      points: DO NP=1,NPTS
         ! weight of this points
         SIM_FAKT=DELTAPHI*WEIGHT((INT((NP-1)/PHPTS)+1))

         RHOT=0
         RHO_ANG=0
         TAUU=0; TAUD=0

         ! calculate the total charge RHOT = RHO/r^2 (up and down)
         ! and the gradient of the charge RHO_ANG
         DO K=1,R%NMAX
            DO LM=1,LMMAX_CALC
               ! total charge
               RHOT(K,1)=RHOT(K,1)+YLM(NP,LM)*RHO(K,LM,1) ! rho(r) Y_L(r)
               RHO_ANG(K,1,:)=RHO_ANG(K,1,:)+YLMD(NP,LM,:)*RHO(K,LM,1) ! rho(r) nabla Y_L(r)
               IF (ISPIN==1) THEN
               ! LM resolved charge density (up,down)
                  RHOLM(K,LM,1)=0.5_q*RHO(K,LM,1)/(R%R(K)*R%R(K))
                  RHOLM(K,LM,2)=0.5_q*RHO(K,LM,1)/(R%R(K)*R%R(K))
               ! magnetization
                  RHOT(K,2)=0
                  RHO_ANG(K,2,:)=0
               ELSE
               ! LM resolved charge density (up,down)
                  RHOLM(K,LM,1)=0.5_q*(RHO(K,LM,1)+RHO(K,LM,2))/(R%R(K)*R%R(K))
                  RHOLM(K,LM,2)=0.5_q*(RHO(K,LM,1)-RHO(K,LM,2))/(R%R(K)*R%R(K))
               ! magnetization
                  RHOT(K,2)=RHOT(K,2)+YLM(NP,LM)*RHO(K,LM,2)
                  RHO_ANG(K,2,:)=RHO_ANG(K,2,:)+YLMD(NP,LM,:)*RHO(K,LM,2)
               ENDIF
            ENDDO
            ! add core charge (spherical) and divide by 1/r^2
            RHOT(K,1)=(RHOT(K,1)+RHOC(K)*YLM(NP,1))/(R%R(K)*R%R(K))
!           RHOT(K,1)=(RHOT(K,1))/(R%R(K)*R%R(K))
            RHOT(K,2)= RHOT(K,2)/(R%R(K)*R%R(K))
            ! divide nabla Y_L through r^3 (1/r is missing in YLMD)
            RHO_ANG(K,:,:)=RHO_ANG(K,:,:)/(R%R(K)*R%R(K))/R%R(K)
            ! add core charge to RHOLM
            RHOLM(K,1,1)=RHOLM(K,1,1)+0.5_q*RHOC(K)/(R%R(K)*R%R(K))
            RHOLM(K,1,2)=RHOLM(K,1,2)+0.5_q*RHOC(K)/(R%R(K)*R%R(K))

            IF (XC%LTAU) THEN
               DO LM=1,LMMAX_TAU
                  IF (ISPIN==1) THEN
                  ! kinetic energy density
                     TAUU(K)=TAUU(K)+0.5_q*TAU(K,LM,1)*YLM(NP,LM)/(R%R(K)*R%R(K))
                     TAUD(K)=TAUU(K)
                  ELSE
                  ! kinetic energy density
                     TAUU(K)=TAUU(K)+TAU(K,LM,1)*YLM(NP,LM)/(R%R(K)*R%R(K))
                     TAUD(K)=TAUD(K)+TAU(K,LM,2)*YLM(NP,LM)/(R%R(K)*R%R(K))
                  ENDIF
               ENDDO
               ! kinetic energy density of the core electrons
               IF (PRESENT(TAUC)) THEN
                  TAUU(K)=TAUU(K)+0.5_q*TAUC(K)*YLM(NP,1)/(R%R(K)*R%R(K))
                  TAUD(K)=TAUD(K)+0.5_q*TAUC(K)*YLM(NP,1)/(R%R(K)*R%R(K))
               ENDIF
               TAUU(K)=MAX(TAUU(K),TAUMIN)
               TAUD(K)=MAX(TAUD(K),TAUMIN)
            ENDIF

         ENDDO

         ! d RHOLM / dr
         DO LM=1,LMMAX_CALC
            CALL GRAD(R,RHOLM(1:R%NMAX,LM,1),RHOLMD(1:R%NMAX,LM,1))
            CALL GRAD(R,RHOLM(1:R%NMAX,LM,2),RHOLMD(1:R%NMAX,LM,2))
         ENDDO

         DO K=1,R%NMAX
            TEMP1=RHOT(K,1)
            RHOT(K,1)=(RHOT(K,1)+RHOT(K,2))*0.5_q    ! spin up
            RHOT(K,2)=(TEMP1-RHOT(K,2))*0.5_q        ! spin down
         ENDDO

         ! d rhot / dr (up,down)
         CALL GRAD(R,RHOT,T1)
         CALL GRAD(R,RHOT(1:R%NMAX,2),T1(1:R%NMAX,2))
         ! d^2 rhot / dr^2 (up,down)
         CALL GRAD(R,T1,T2)
         CALL GRAD(R,T1(1:R%NMAX,2),T2(1:R%NMAX,2))

         DO K=1,R%NMAX

            ! calculate the laplacian of the density
            LAPLUP=0; LAPLDW=0
            IF (XC%LLAP) THEN
               DO LM=1,LMMAX_CALC
                  LAPLUP=LAPLUP+ &
                 &    YLMDD(NP,LM,1)*RHOLM(K,LM,1)/(R%R(K)*R%R(K)) &
                 &   +(YLMD(NP,LM,1)*XYZPTS(NP,1)+YLMD(NP,LM,1)*XYZPTS(NP,1))/R%R(K)*RHOLMD(K,LM,1) &
                 &   +YLMDD(NP,LM,4)*RHOLM(K,LM,1)/(R%R(K)*R%R(K)) &
                 &   +(YLMD(NP,LM,2)*XYZPTS(NP,2)+YLMD(NP,LM,2)*XYZPTS(NP,2))/R%R(K)*RHOLMD(K,LM,1) &
                 &   +YLMDD(NP,LM,6)*RHOLM(K,LM,1)/(R%R(K)*R%R(K)) &
                 &   +(YLMD(NP,LM,3)*XYZPTS(NP,3)+YLMD(NP,LM,3)*XYZPTS(NP,3))/R%R(K)*RHOLMD(K,LM,1)
                  LAPLDW=LAPLDW+ &
                 &    YLMDD(NP,LM,1)*RHOLM(K,LM,2)/(R%R(K)*R%R(K)) &
                 &   +(YLMD(NP,LM,1)*XYZPTS(NP,1)+YLMD(NP,LM,1)*XYZPTS(NP,1))/R%R(K)*RHOLMD(K,LM,2) &
                 &   +YLMDD(NP,LM,4)*RHOLM(K,LM,2)/(R%R(K)*R%R(K)) &
                 &   +(YLMD(NP,LM,2)*XYZPTS(NP,2)+YLMD(NP,LM,2)*XYZPTS(NP,2))/R%R(K)*RHOLMD(K,LM,2) &
                 &   +YLMDD(NP,LM,6)*RHOLM(K,LM,2)/(R%R(K)*R%R(K)) &
                 &   +(YLMD(NP,LM,3)*XYZPTS(NP,3)+YLMD(NP,LM,3)*XYZPTS(NP,3))/R%R(K)*RHOLMD(K,LM,2)
               ENDDO

               LAPLUP=LAPLUP- &
              &    XYZPTS(NP,1)*XYZPTS(NP,1)/(R%R(K))*T1(K,1)+XYZPTS(NP,1)*XYZPTS(NP,1)*T2(K,1) &
              &   -XYZPTS(NP,2)*XYZPTS(NP,2)/(R%R(K))*T1(K,1)+XYZPTS(NP,2)*XYZPTS(NP,2)*T2(K,1) &
              &   -XYZPTS(NP,3)*XYZPTS(NP,3)/(R%R(K))*T1(K,1)+XYZPTS(NP,3)*XYZPTS(NP,3)*T2(K,1) &
              &   +3*T1(K,1)/R%R(K)

               LAPLDW=LAPLDW- &
              &    XYZPTS(NP,1)*XYZPTS(NP,1)/(R%R(K))*T1(K,2)+XYZPTS(NP,1)*XYZPTS(NP,1)*T2(K,2) &
              &   -XYZPTS(NP,2)*XYZPTS(NP,2)/(R%R(K))*T1(K,2)+XYZPTS(NP,2)*XYZPTS(NP,2)*T2(K,2) &
              &   -XYZPTS(NP,3)*XYZPTS(NP,3)/(R%R(K))*T1(K,2)+XYZPTS(NP,3)*XYZPTS(NP,3)*T2(K,2) &
              &   +3*T1(K,2)/R%R(K)
            ENDIF

            ! norm of gradient for spin up
            XU  =T1(K,1)*XYZPTS(NP,1)+(RHO_ANG(K,1,1)+RHO_ANG(K,2,1))*0.5_q
            YU  =T1(K,1)*XYZPTS(NP,2)+(RHO_ANG(K,1,2)+RHO_ANG(K,2,2))*0.5_q
            ZU  =T1(K,1)*XYZPTS(NP,3)+(RHO_ANG(K,1,3)+RHO_ANG(K,2,3))*0.5_q
            T1(K,1)=SQRT(XU*XU+YU*YU+ZU*ZU)
            ! norm of gradient for spin down
            XD  =T1(K,2)*XYZPTS(NP,1)+(RHO_ANG(K,1,1)-RHO_ANG(K,2,1))*0.5_q
            YD  =T1(K,2)*XYZPTS(NP,2)+(RHO_ANG(K,1,2)-RHO_ANG(K,2,2))*0.5_q
            ZD  =T1(K,2)*XYZPTS(NP,3)+(RHO_ANG(K,1,3)-RHO_ANG(K,2,3))*0.5_q

            T1(K,2)=SQRT(XD*XD+YD*YD+ZD*ZD)
!           TG=T1(K,1)+T1(K,2) incorrect definition
            TG=SQRT((XU+XD)**2+(YU+YD)**2+(ZU+ZD)**2)

!           not needed for MBJ
            IF ((XC%ID(1)/=ID_XC_MBJ) .AND. (XC%ID(1)/=ID_XC_LMBJ)) THEN
               T1(K,1)=MAX(T1(K,1),NABMIN)
               T1(K,2)=MAX(T1(K,2),NABMIN)
               TG=MAX(TG,2*NABMIN)
            ENDIF

            RHOT(K,1)=MAX(RHOT(K,1), CHGMIN)
            RHOT(K,2)=MAX(RHOT(K,2), CHGMIN)

            CALL METAGGA_DRIVER(&
           &   XC,RHOT(K,1)*AUTOA3,RHOT(K,2)*AUTOA3,T1(K,1)*AUTOA4,T1(K,2)*AUTOA4,TG*AUTOA4, &
           &   LAPLUP*AUTOA5,LAPLDW*AUTOA5,TAUU(K)*EVTOH,TAUD(K)*EVTOH, &
           &   EXT,DEXC1,DEXC2,DVXC1,DVXC2,DVC,MUUP,MUDW,DLXC1,DLXC2,CVMBJ)

            IF (LscMBJ()) GRDR=GRDR+(TG*R%R(K)**2*AUTOA*R%SI(K)*SIM_FAKT)/(RHOT(K,1)+RHOT(K,2))

            IF (XC%ID(1)==ID_XC_LMBJ) THEN
            ! TR: Spin-up
              ! factor 1/2 because we handle a single spin channel only
              GRHO_OVER_RHO_AUG(NI)=GRHO_OVER_RHO_AUG(NI)+0.5_q*AUGFACT &
                      *T1(K,1)*R%R(K)**2*AUTOA*R%SI(K)*SIM_FAKT/RHOT(K,1)
            ! TR: Spin-down
              ! factor 1/2 because we handle a single spin channel only
              GRHO_OVER_RHO_AUG(NI)=GRHO_OVER_RHO_AUG(NI)+0.5_q*AUGFACT &
                      *T1(K,2)*R%R(K)**2*AUTOA*R%SI(K)*SIM_FAKT/RHOT(K,2)
            ENDIF

            DEXC=DEXC+(EXT*RYTOEV)*(RHOT(K,1)+RHOT(K,2))*R%R(K)*R%R(K)*R%SI(K)*SIM_FAKT

            ! develop
            !             d    f_xc     grad rho
            !             ------------  --------  = vec dFXC(r)
            !             d |grad rho| |grad rho|
            ! into spherical harmonics
            ! afterwards vec dFXC can be reconstructed using
            ! dFXC_i(r) =  \sum_L dFXC_Li(r) Y_L(r)

            DEXC1=DEXC1*RYTOEV
            DEXC2=DEXC2*RYTOEV
            DVXC1=DVXC1*RYTOEV*AUTOA
            DVXC2=DVXC2*RYTOEV*AUTOA
            DVC  =DVC  *RYTOEV*AUTOA

            DVXC1=DVXC1/ MAX(T1(K,1),NABMIN)
            DVXC2=DVXC2/ MAX(T1(K,2),NABMIN)
            DVC  =DVC  / MAX(TG,2*NABMIN)

            DO LM=1,LMMAX
               dFXC(K,LM,0,1) = dFXC(K,LM,0,1)+DEXC1*SIM_FAKT*YLM(NP,LM)
               dFXC(K,LM,0,2) = dFXC(K,LM,0,2)+DEXC2*SIM_FAKT*YLM(NP,LM)

               dFXC(K,LM,1,1) = dFXC(K,LM,1,1)+(XU* DVXC1 + (XU+XD) * DVC)*SIM_FAKT*YLM(NP,LM)
               dFXC(K,LM,2,1) = dFXC(K,LM,2,1)+(YU* DVXC1 + (YU+YD) * DVC)*SIM_FAKT*YLM(NP,LM)
               dFXC(K,LM,3,1) = dFXC(K,LM,3,1)+(ZU* DVXC1 + (ZU+ZD) * DVC)*SIM_FAKT*YLM(NP,LM)

               dFXC(K,LM,1,2) = dFXC(K,LM,1,2)+(XD* DVXC2 + (XU+XD) * DVC)*SIM_FAKT*YLM(NP,LM)
               dFXC(K,LM,2,2) = dFXC(K,LM,2,2)+(YD* DVXC2 + (YU+YD) * DVC)*SIM_FAKT*YLM(NP,LM)
               dFXC(K,LM,3,2) = dFXC(K,LM,3,2)+(ZD* DVXC2 + (ZU+ZD) * DVC)*SIM_FAKT*YLM(NP,LM)
            ENDDO

            ! same thing, if laplacian-dependent functional is called,
            ! resolve d e_xc / d (lap n_sigma) into spherical harmonics
            IF (XC%LDLAP) THEN
               DLXC1=DLXC1*RYTOEV*AUTOA2
               DLXC2=DLXC2*RYTOEV*AUTOA2
               DO LM=1,LMMAX
                  dFXCdLap(K,LM,1) = dFXCdLap(K,LM,1)+DLXC1*SIM_FAKT*YLM(NP,LM)
                  dFXCdLap(K,LM,2) = dFXCdLap(K,LM,2)+DLXC2*SIM_FAKT*YLM(NP,LM)
               ENDDO
            ENDIF

            IF (XC%LMU) THEN
               MUUP=0.5_q*MUUP*RYTOEV*AUTOA2
               MUDW=0.5_q*MUDW*RYTOEV*AUTOA2
               DO LM=1,LMMAX_TAU
                  dFXCdT(K,LM,1) = dFXCdT(K,LM,1)+MUUP*SIM_FAKT*YLM(NP,LM)
                  dFXCdT(K,LM,2) = dFXCdT(K,LM,2)+MUDW*SIM_FAKT*YLM(NP,LM)
               ENDDO
            ENDIF
         ENDDO
      ENDDO points

GPACC CALL VTUTOR_OFFLOAD_CHECK()

!========================================================================
!
! GGA part related to d e_XC / grad n
!
! add -\sum_i Y_L  d (dFXC_L'i(x) Y_L'(x))/ d x_i to the potential V_L(r)
!
!                          x_i d dFXC_L'i(|x|)
! -\sum_L'i Y_L(x) Y_L'(x)  -  --------------- + dFXC_L'i(x) Y_L nabla_i Y_L'(x)
!                          |x|   d |x|
!
!========================================================================
      DO ISP=1,ISPIN
         DO LM=1,LMMAX_CALC
            DO K=1,R%NMAX
               POT(K,LM,ISP)=POT(K,LM,ISP)+dFXC(K,LM,0,ISP)
            ENDDO
         ENDDO
      ENDDO
      ! nabla comment:
      ! the matrix elements < Y_L| nabla | Y_L'> are non zero
      ! for L=L'+-1
      DO ISP=1,ISPIN
         DO I=1,3
            DO LMP=1,LMMAX
               CALL GRAD(R,dFXC(1:R%NMAX,LMP,I,ISP),T1(1:R%NMAX,1))
               DO LM=1,LMMAX_CALC
                  IF (ABS(YLM_NABLA_YLM(LM,LMP,I))>1E-4 .OR. ABS(YLM_X_YLM(LM,LMP,I))>1E-4) THEN
                     DO K=1,R%NMAX
                        POT(K,LM,ISP)=POT(K,LM,ISP)-( &
                             dFXC(K,LMP,I,ISP)*YLM_NABLA_YLM(LM,LMP,I)/ R%R(K)+ &
                             T1(K,1)*YLM_X_YLM(LM,LMP,I) &
                           )
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DEALLOCATE(dFXC)

      ! v_xc_sigma = [d/dn - div d/d(grad n) + lap d/d/(lap n)] e_xc_sigma
      ! this adds the laplacian contribution
      ! lap R(r)Y_lm = R''(r)Y_lm + 2/r*R'(r) Y_lm - l(l+1)/r**2 R(r)Y_lm
      ! T1(:,1) = R'(r)
      ! T2(:,1) = R''(r), just to make the distinction clearer

      IF (XC%LDLAP) THEN
         LMAX_TMP = LMAX_CALC
         IF (.NOT.LASPH) LMAX_TMP = 0
         DO ISP=1,ISPIN
           DO ell=0,LMAX_TMP
             DO M=1,2*ell+1
               lm = ell*ell + m
               CALL GRAD(R,dFXCdLap(1:R%NMAX,LM,ISP),T1(1:R%NMAX,1))
               CALL GRAD(R,T1(1:R%NMAX,1),T2(1:R%NMAX,1))
               lfac = ell*(ell+1)
               DO K=1,R%NMAX
                 POT(K,LM,ISP)=POT(K,LM,ISP) + T2(K,1) + 2*T1(K,1)/R%R(K) - lfac*dFXCdLap(K,LM,ISP)/R%R(K)**2
               ENDDO
             enddo
           ENDDO
         ENDDO
         DEALLOCATE(dFXcdLap)
      ENDIF

      DCMU=0
      IF (XC%LMU) THEN
         DO ISP=1,ISPIN
            DO LM=1,LMMAX_TAU
               DO K=1,R%NMAX
                  MU(K,LM,ISP)=dFXCdT(K,LM,ISP)
                  DCMU=DCMU+MU(K,LM,ISP)*TAU(K,LM,ISP)*R%SI(K)/HSQDTM
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(dFXcdT)
      ENDIF

      ! finally calculate the double counting corrections again
      DVXC=0
      CALL RAD_POT_IN_RHO(R,ISPIN,LMAX_CALC,RHO,POT,DVXC)
      DVXC=DVXC-DVINI

      IF (XC%LMU) THEN
         DVXC=DVXC+DCMU
!        WRITE(*,'(A,F24.14)') 'dcmu =',DCMU
      ENDIF

      DEALLOCATE(YLM,YLMD,YLMDD,XYZPTS,RADPTS,YLM_NABLA_YLM,YLM_X_YLM,WEIGHT,ABSCIS,RHOLM,RHOLMD)

      RETURN
      END SUBROUTINE RAD_METAGGA_XC


!********************* SUBROUTINE RAD_PROJ_METAGGA *********************
!
!***********************************************************************
      SUBROUTINE RAD_PROJ_METAGGA(&
     &   P,XC,W,LMAX_TAU,MU,A,DLLMM)
      USE prec
      USE asa
      USE pseudo
      USE setexm_struct_def, ONLY : LMAXTAU,xc_info
      USE constant
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE (potcar) P
      TYPE (xc_info) XC
      INTEGER LMAX_TAU
      REAL(q) A
      REAL(q) W(:,:)
      REAL(q) MU(:,:,:)
      OVERLAP DLLMM(:,:,:)
      ! local variables
      INTEGER ISP,I,K,Ch1,Ch2
      INTEGER LLMAX_TAU,LMMAX_TAU,LM,LMP,LL,LLP,M,MP,MPLOW,INDYLM,INDYLMP

      REAL(q) SUM
      REAL(q) W1(P%R%NMAX),W2(P%R%NMAX),DW1(P%R%NMAX),DW2(P%R%NMAX)
      REAL(q), ALLOCATABLE :: TAU(:,:)

      ! quick return if possible
      IF ((.NOT.XC%LDOMETAGGA).OR.(.NOT.XC%LMU)) RETURN

      ! non-collinear magnetism
!     IF (SIZE(DLLMM,3)==4) THEN
!        WRITE(*,*) 'RAD_PROJ_METAGGA: internal error; non-collinear case not supported (yet).'
!        STOP
!     ENDIF

      ! check consistency
      LMMAX_TAU=(LMAX_TAU+1)**2
      IF (SIZE(MU,2)<LMMAX_TAU) THEN
         CALL vtutor%error("RAD_PROJ_METAGGA: ERROR, size(MU,2)<LMMAX_TAU " // str(SIZE(MU,2)) // " &
            &" // str(LMMAX_TAU))
      ENDIF

      ! restriction for the angular development of \tau and \mu
      LLMAX_TAU=MIN(LMAXTAU,LMAX_TAU); LMMAX_TAU=(LLMAX_TAU+1)**2

      ALLOCATE(TAU(P%R%NMAX,LMMAX_TAU))

      ! loops over channels (l,epsilon)
      LM=1; DO Ch1=1,P%LMAX

      ! d(w1/r)/dr
      W1(:)=W(:,Ch1)/P%R%R(:); CALL GRAD(P%R,W1,DW1)

      LMP=LM; DO Ch2=Ch1,P%LMAX
      ! d(w2/r)/dr
      W2(:)=W(:,Ch2)/P%R%R(:); CALL GRAD(P%R,W2,DW2)

         LL=P%LPS(Ch1); LLP=P%LPS(Ch2)
         DO M=1,2*LL+1
         MPLOW=1; IF(Ch1==Ch2) MPLOW=M
         DO MP=MPLOW,2*LLP+1
            INDYLM= LL**2 +M
            INDYLMP=LLP**2+MP

            TAU=0

            DO I=1,LMMAX_TAU
               TAU(:,I)=TAU(:,I)+ &
              &   (DW1(:)*X_YLM_X_YLM(INDYLM,INDYLMP,I)*DW2(:)*P%R%R(:)*P%R%R(:) + &
              &    DW1(:)*X_YLM_YLMD(INDYLM,INDYLMP,I)*W2(:)*P%R%R(:) + &
              &    W1(:)*X_YLM_YLMD(INDYLMP,INDYLM,I)*DW2(:)*P%R%R(:) + &
              &    W1(:)*YLMD_YLMD(INDYLM,INDYLMP,I)*W2(:))
            ENDDO

            DO ISP=1,SIZE(DLLMM,3)
               SUM=0
               DO I=1,LMMAX_TAU
                  ! Integrate over r (Simpson integration)
                  DO K=1,P%R%NMAX
                     SUM=SUM+MU(K,I,ISP)*TAU(K,I)*P%R%SI(K)
                  ENDDO
               ENDDO
               DLLMM(LM+M-1,LMP+MP-1,ISP)=DLLMM(LM+M-1,LMP+MP-1,ISP)+A*SUM
               DLLMM(LMP+MP-1,LM+M-1,ISP)=DLLMM(LM+M-1,LMP+MP-1,ISP)
            ENDDO

         ENDDO
         ENDDO
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

!     CALL DUMP_DLLMM( "NABLA-MU-NABLA",DLLMM(:,:,1),P)

      DEALLOCATE(TAU)

      RETURN
      END SUBROUTINE RAD_PROJ_METAGGA


!********************* SUBROUTINE RAD_AUXILIARY_FUNCTIONS_METAGGA ******
!
!***********************************************************************
      SUBROUTINE RAD_AUXILIARY_FUNCTIONS_METAGGA(XC,LLMAX)
      USE prec
      USE asa
      USE constant
      USE pseudo
      USE setexm_struct_def, ONLY : xc_info, PHPTS_FACT, THPTS_FACT, LMAXTAU
      IMPLICIT NONE
      TYPE (xc_info) XC
      INTEGER LLMAX
      ! local variables
      INTEGER I,J,Ch1,Ch2
      INTEGER IFAIL,PHPTS,THPTS,NPTS,NP
      INTEGER LMMAX,LLMAX_TAU,LMMAX_TAU
      INTEGER LM,LMP,LL,LLP,M,MP,INDYLM,INDYLMP

      REAL(q), ALLOCATABLE :: RADPTS(:,:), XYZPTS(:,:)
      REAL(q), ALLOCATABLE :: YLM(:,:),YLMD(:,:,:)
      REAL(q), ALLOCATABLE :: WEIGHT(:),ABSCIS(:)

      REAL(q) DELTAPHI,SIM_FAKT
      EXTERNAL GAUSSI2

      ! quick return if possible
      IF (.NOT.XC%LDOMETAGGA) RETURN

      IF (ALLOCATED(X_YLM_X_YLM).AND.ALLOCATED(YLMD_YLMD).AND.ALLOCATED(X_YLM_YLMD)) RETURN

      LMMAX=(LLMAX+1)**2
      LLMAX_TAU=MIN(2*LLMAX+2,LMAXTAU); LMMAX_TAU=(LLMAX_TAU+1)**2

      IF (ALLOCATED(X_YLM_X_YLM)) DEALLOCATE(X_YLM_X_YLM)
      ALLOCATE(X_YLM_X_YLM(LMMAX,LMMAX,LMMAX_TAU)); X_YLM_X_YLM=0

      IF (ALLOCATED(YLMD_YLMD)) DEALLOCATE(YLMD_YLMD)
      ALLOCATE(YLMD_YLMD(LMMAX,LMMAX,LMMAX_TAU)); YLMD_YLMD=0

      IF (ALLOCATED(X_YLM_YLMD)) DEALLOCATE(X_YLM_YLMD)
      ALLOCATE(X_YLM_YLMD(LMMAX,LMMAX,LMMAX_TAU)); X_YLM_YLMD=0


!========================================================================
! number of theta and phi pivot points to perform angular integration
! since Exc=f(a*Yllmax,m) we need more pivot points than theoretically
! needed to integrate Yllmax,m.
! the factor 2 is the minium, 3 is more accurate
!========================================================================
      PHPTS=PHPTS_FACT*(MAX(LLMAX,LLMAX_TAU)+1)
      THPTS=THPTS_FACT*FLOOR(REAL(MAX(LLMAX,LLMAX_TAU)/2+1,KIND=q))
      NPTS=PHPTS*THPTS
      DELTAPHI=REAL(2_q*PI/PHPTS,KIND=q)
! allocate arrays
      ALLOCATE(XYZPTS(NPTS,3),RADPTS(NPTS,2),WEIGHT(THPTS),ABSCIS(THPTS), &
     &     YLM(NPTS,MAX(LMMAX,LMMAX_TAU)),YLMD(NPTS,MAX(LMMAX,LMMAX_TAU),3))

      ! set phi positions, equally spaced
      RADPTS=0; WEIGHT=0; ABSCIS=0
      DO I=1,PHPTS
         DO J=1,THPTS
            RADPTS((J-1)*PHPTS+I,2)=(I-1)*DELTAPHI
         ENDDO
      ENDDO
     ! get theta positions (actually get cos(theta)) (Gauss integration)
      CALL GAUSSI(GAUSSI2,-1._q,1._q,0,THPTS,WEIGHT,ABSCIS,IFAIL)
      DO I=1,THPTS
         RADPTS((I-1)*PHPTS+1:I*PHPTS,1)=ABSCIS(I)
      ENDDO
      ! convert radial to cartesian coordinates
      DO I=1,NPTS
         XYZPTS(I,1)=COS(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! x
         XYZPTS(I,2)=SIN(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! y
         XYZPTS(I,3)=RADPTS(I,1)                                 ! z
      ENDDO

      ! get |r| Y_lm on a unit sphere and its derivatives
      YLM=0; YLMD=0
      CALL SETYLM_GRAD2(MAX(LLMAX,LLMAX_TAU),NPTS,YLM,YLMD,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))

      ! loop over all points in the angular grid
      points: DO NP=1,NPTS

         ! weight of this points
         SIM_FAKT=DELTAPHI*WEIGHT((INT((NP-1)/PHPTS)+1))

         ! loops over channels (l,epsilon)
         DO LL =0,LLMAX
         DO LLP=0,LLMAX

            DO M=1,2*LL+1
            DO MP=1,2*LLP+1
               INDYLM= LL**2 +M
               INDYLMP=LLP**2+MP
               DO I=1,LMMAX_TAU
                  X_YLM_X_YLM(INDYLM,INDYLMP,I)=X_YLM_X_YLM(INDYLM,INDYLMP,I)+ SIM_FAKT*YLM(NP,I)* &
                 &   (XYZPTS(NP,1)*XYZPTS(NP,1)+XYZPTS(NP,2)*XYZPTS(NP,2)+XYZPTS(NP,3)*XYZPTS(NP,3))* &
                 &   YLM(NP,INDYLM)*YLM(NP,INDYLMP)
                  X_YLM_YLMD(INDYLM,INDYLMP,I)=X_YLM_YLMD(INDYLM,INDYLMP,I)+ SIM_FAKT*YLM(NP,I)* &
                 &   (XYZPTS(NP,1)*YLMD(NP,INDYLMP,1)+XYZPTS(NP,2)*YLMD(NP,INDYLMP,2)+XYZPTS(NP,3)*YLMD(NP,INDYLMP,3))* &
                 &   YLM(NP,INDYLM)
                  YLMD_YLMD(INDYLM,INDYLMP,I)=YLMD_YLMD(INDYLM,INDYLMP,I)+ SIM_FAKT*YLM(NP,I)* &
                 &   (YLMD(NP,INDYLM,1)*YLMD(NP,INDYLMP,1)+YLMD(NP,INDYLM,2)*YLMD(NP,INDYLMP,2)+YLMD(NP,INDYLM,3)*YLMD(NP,INDYLMP,3))
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO

      ENDDO points

      DEALLOCATE(XYZPTS,RADPTS,WEIGHT,ABSCIS,YLM,YLMD)

      RETURN
      END SUBROUTINE RAD_AUXILIARY_FUNCTIONS_METAGGA


!******************** FUNCTION IR_MATCH ********************************
!
!***********************************************************************
      FUNCTION IR_MATCH(PP,CHANNEL)
      USE prec
      USE pseudo
      IMPLICIT NONE
      TYPE(potcar), POINTER :: PP
      INTEGER CHANNEL
      INTEGER IR_MATCH
      ! local variables
      INTEGER K
      IR_MATCH=PP%R%NMAX
      DO K=PP%R%NMAX,1,-1
         IF (PP%WAE(K,CHANNEL)/=PP%WPS(K,CHANNEL)) EXIT
         IR_MATCH=K
      ENDDO
      END FUNCTION IR_MATCH


!******************** FUNCTION IR_PSMAX ********************************
!
!***********************************************************************
      FUNCTION IR_PSMAX(PP)
      USE prec
      USE pseudo
      IMPLICIT NONE
      TYPE(potcar), POINTER :: PP
      INTEGER IR_PSMAX
      ! local variables
      INTEGER K
      IR_PSMAX=PP%R%NMAX
      DO K=PP%R%NMAX,2,-1
         IF ((PP%R%R(K)<PP%PSDMAX).OR.(PP%R%R(K)<PP%PSRMAX)) EXIT
         IR_PSMAX=K
      ENDDO
      END FUNCTION IR_PSMAX

#endif
  END MODULE radial
