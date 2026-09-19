#include "symbol.inc"
!#define RPAgamma
MODULE esf
  USE chi_glb
  USE chi_base
  USE fock 
  USE tutor, ONLY: vtutor
  USE mkpoints_struct_def
  USE esf_struct_def
  IMPLICIT NONE

  !< number of linear regression points (or individual cutoffs of chi* V )
  INTEGER, PARAMETER :: NE_REG=8

  REAL(q), SAVE :: DATAKEMAX         !< maximum energy used in Coulomb kernel

  REAL(q), SAVE :: ESF_CONV = 0.01_q !< convergence threashold for spline fit refinement
  INTEGER, SAVE :: ESF_NINTER = 15   !< number of interstitial grid points for spline fit

  CONTAINS 

   SUBROUTINE ESF_READER(INCAR)
      USE incar_reader, ONLY: incar_file, &
                              process_incar, &
                              count_elements
      USE tutor, ONLY: vtutor  
      USE string, ONLY: str 

      IMPLICIT NONE
      TYPE(incar_file),       INTENT(INOUT)  :: INCAR

      ESF_CONV = 0.01_q 
      CALL PROCESS_INCAR(INCAR, "ESF_CONV", ESF_CONV )

      ESF_NINTER = 15 
      CALL PROCESS_INCAR(INCAR, "ESF_NINTER", ESF_NINTER )
      IF ( ESF_NINTER < 0 ) THEN
         CALL vtutor%error( "negative number of interstitial points for spline fit given "//str( ESF_NINTER ) )
      ENDIF 
   END SUBROUTINE ESF_READER



!************************ SUBROUTINE ROTLN    **************************
!
!> this subroutine calculates Log(1-A)+A
!> old version no longer used
!> several test cases are implemented as well
!> the second order expansion is A*A and should correspond exactly
!> to the Coulomb term in the MP2 energy
!>
!> @todo: the routine should be rewritten to conserve memory
!
!***********************************************************************

  SUBROUTINE ROTLN(CUNI, NBANDS, LREAL, LAMBDA, IMP2, TRACE, IU6, HFEIG_ )
    INTEGER NBANDS
    GDEF       :: CUNI(:,:)
    LOGICAL    :: LREAL
    REAL(q)    :: LAMBDA
    INTEGER    :: IMP2   ! IMP2=0 use RPA, IMP2=1 yield second order MP2 term
    COMPLEX(q) :: TRACE
    INTEGER    :: IU6
    REAL(q),OPTIONAL :: HFEIG_(:)
    ! local
    INTEGER, PARAMETER :: LWORK=32
    COMPLEX(q) ::CWRK(LWORK*NBANDS)
    REAL(q) :: HFEIG(NBANDS),W(3*NBANDS)
    INTEGER :: IFAIL, NDIM, N1, N2

    GDEF, ALLOCATABLE ::  CTMP(:,:), CEIDB(:,:)

    PROFILING_START( 'rotln' )

    NDIM = SIZE(CUNI,1)
    ALLOCATE(CTMP(NBANDS,NBANDS),CEIDB(NBANDS,NBANDS))

!=======================================================================
! diagononalize the matrix to get eigenvalues (N=eigenv, V=eigenvectors)
!=======================================================================
    IF (LREAL) THEN
#ifndef gammareal
       CALL vtutor%bug("LREAL is set, but complex code", __FILE__, __LINE__)
#endif
       CALL DSYEV &
            ('V','U',NBANDS,CUNI(1,1),NDIM,HFEIG,CWRK,LWORK*NBANDS, &
            IFAIL)
    ELSE
#ifdef gammareal
       CALL vtutor%bug("LREAL is not set, but gamma only code", __FILE__, __LINE__)
#endif
       CALL ZHEEV &
         ('V','U',NBANDS,CUNI(1,1),NDIM,HFEIG,CWRK,LWORK*NBANDS, &
         W,  IFAIL)
    ENDIF
!=======================================================================
! diagononalize the matrix
!=======================================================================
    IF (IFAIL/=0) THEN
       CALL vtutor%bug("Call to routine ZHEEV failed! Error code was " // &
          str(IFAIL), __FILE__, __LINE__)
    ENDIF
    IF (IU6>=0) WRITE(IU6,'(10F10.6)') HFEIG
!=======================================================================
! set up the unitary transformation CEIDB for wavefunctions
!=======================================================================

    DO N2=1,NBANDS
       IF (1-HFEIG(N2)<1E-8) THEN
          CTMP(:,N2)=0
          WRITE(*,*) 'eigenvalues are negative or close to zero', CTMP
       ELSE
          IF (IMP2==0) THEN
             ! log (1-A) +A
             TRACE=TRACE+((LOG(1-HFEIG(N2)*LAMBDA)+HFEIG(N2)*LAMBDA))
             CTMP(1:NBANDS,N2)=CUNI(1:NBANDS,N2)*(LOG(1-HFEIG(N2)*LAMBDA)+HFEIG(N2)*LAMBDA)!/MAX(HFEIG(N2)*HFEIG(N2),1E-32_q)
!             CTMP(1:NBANDS,N2)=CUNI(1:NBANDS,N2)*(-0.5_q)/(1.0_q-HFEIG(N2))
          ELSE
             TRACE=TRACE+(-HFEIG(N2)*HFEIG(N2))/2
             ! A*A
             CTMP(1:NBANDS,N2)=CUNI(1:NBANDS,N2)*(-0.5_q)
          ENDIF
       ENDIF
       IF ( PRESENT( HFEIG_ ) )THEN
          IF ( N2 > SIZE( HFEIG_) ) THEN
             CALL vtutor%bug("HFEIG_ too small "//&
                str( NBANDS )//" "//str(SIZE( HFEIG_ )), __FILE__, __LINE__)
          ENDIF
          HFEIG_( N2 ) = HFEIG( N2 ) 
       ENDIF
    ENDDO
    CALL GGEMM( 'N', trans, NBANDS, NBANDS, NBANDS, one, CTMP, &
         &             NBANDS, CUNI(1,1), NDIM, zero, CEIDB, NBANDS)
   
    CUNI=0 
    CUNI(1:NBANDS,1:NBANDS)=CEIDB(1:NBANDS,1:NBANDS)

    DEALLOCATE(CTMP, CEIDB)

    PROFILING_STOP( 'rotln' )

  END SUBROUTINE ROTLN


!******************** SUBROUTINE SETUP_COULOMB_DIAGONAL ****************
!
!> this subroutine sets up the squareroot of the truncated Coulomb 
!> diagonal used in XI_LOCAL_FIELD_ACFDT
!
!***********************************************************************

  SUBROUTINE SETUP_COULOMB_DIAGONAL(LATT_CUR, ENCUT, ENCUTSOFT, WGWQ, MAXINDEX, LSR, OLD_INDEX, DATAKE )
    USE constant
    USE fock
    USE c2f_interface, ONLY : ERRF
    USE fermi_energy, ONLY : DELSTP
    IMPLICIT NONE
    TYPE (latt)                     :: LATT_CUR
    REAL(q), INTENT( IN )           :: ENCUT
    REAL(q), INTENT( IN )           :: ENCUTSOFT
    TYPE (wavedes1)                 :: WGWQ
    INTEGER, INTENT(INOUT)          :: MAXINDEX !< number of grid points using truncated Coulomb kernel
    LOGICAL, INTENT( IN )           :: LSR
    REAL(q), POINTER, INTENT(INOUT) :: DATAKE(:)!< result 
    INTEGER, POINTER, INTENT(INOUT) :: OLD_INDEX(:)
  ! local
    INTEGER    NI, NP, NI_
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE, POTFAK, E
    REAL(q) :: DFUN, SFUN

    REAL(q) :: OMEGBK,QC
    REAL(q) :: GAMMASCALE

! multipole correction
    INTEGER    NJ,NK
    COMPLEX(q) :: TSUM(10)
    COMPLEX(q) :: POTCORRECTION(10)
    REAL(q) :: Q1, Q2

    PROFILING_START( 'setup_coulomb_diagonal' )

    IF (ENCUTSOFT>=0) THEN
       Q1=SQRT(ENCUTSOFT/HSQDTM)
       Q2=SQRT(ENCUT/HSQDTM)
    ENDIF

    ! e^2/ volume 
    SCALE=EDEPS/LATT_CUR%OMEGA
!=======================================================================
! first set up the truncated Coulomb kernel
! smoothly going from 1/G^2 to zero between ENCUTGWSOFT and ENCUTGW
!=======================================================================
    ! convert wavevector q into cartesian corrdinates
    DKX=(WGWQ%VKPT(1))*LATT_CUR%B(1,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(1,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(1,3)
    DKY=(WGWQ%VKPT(1))*LATT_CUR%B(2,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(2,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(2,3)
    DKZ=(WGWQ%VKPT(1))*LATT_CUR%B(3,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(3,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(3,3)

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    IF( ASSOCIATED( DATAKE ) ) THEN
       !CALL vtutor%bug( "DATAKE associated, wrong usage of this routine",__FILE__,__LINE__)
       DEALLOCATE( DATAKE ) 
       NULLIFY( DATAKE )  
    ENDIF 
    ALLOCATE( DATAKE(NP) )

    IF( ASSOCIATED( OLD_INDEX ) ) THEN
       !CALL vtutor%bug( "OLD_INDEX associated, wrong usage of this routine",__FILE__,__LINE__)
       DEALLOCATE( OLD_INDEX ) 
       NULLIFY( OLD_INDEX)  
    ENDIF 
    ALLOCATE( OLD_INDEX(NP))

    DATAKEMAX=0.0
    MAXINDEX=0
    DO NI=1,NP
       NI_=NI
       IF (WGWQ%LGAMMA) NI_=(NI-1)/2+1
       ! calculate the cartesian components of G
       GX=(WGWQ%IGX(NI_)*LATT_CUR%B(1,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WGWQ%IGX(NI_)*LATT_CUR%B(2,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WGWQ%IGX(NI_)*LATT_CUR%B(3,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
       ! calculate the square of the length (G+q)^2
       GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2

       IF (ABS(GSQU)<G2ZERO) THEN
          ! head and wing
          POTFAK=SCALE
          ! if HFRCUT is set, a fixed spherical cutoff is used
          ! resulting in a finite values at G=0
          ! since the head stores the q^2, the final contribution is zero
          ! from the head
          IF (HFRCUT/=0) THEN
             POTFAK=0
          ! if LRHFCALC.AND.LRSCOR the Coulomb kernel is replaced by v_{LR}=v*exp(-q^2/(4*mu^2))
          ! for q=0  v_{LR}=1, this introduces an error for small mu because of the finite q-grid
          ! In order to correct this error, v_{LR}(q=0) is set to average of v_{LR} in a sphere  
          ! with volume corresponding to BZ volume/ # q-points       
          ELSE IF (LRHFCALC.AND.LRSCOR) THEN
             CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
             QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
             GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
            &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
            &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
             POTFAK=POTFAK*GAMMASCALE
          ! if LSR the Coulomb kernel is replaced by v_{SR}=v*[1-exp(-q^2/(4*mu^2))]
          ! for q=0  v_{SR}=0, this introduces an error for small mu because of the finite q-grid
          ! In order to correct this error, v_{SR}(q=0) is set to average of v_{SR} in a sphere
          ! with volume corresponding to BZ volume/ # q-points
          ELSE IF (LSR) THEN
             CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
             QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
             GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
            &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
            &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
             POTFAK=POTFAK*(1-GAMMASCALE)
          ENDIF

          ! switch off standard convergence correction
          IF(MCALPHA/=0) THEN
             POTFAK=0
          ENDIF

       ELSE
          ! the factor 1/(2 pi)^2 is required to obtain proper reciprocal
          ! lattice vector lenght
          POTFAK=SCALE/(GSQU*TPI**2)

          IF (HFRCUT/=0) THEN
             ! spherical cutoff on Coloumb kernel
             ! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
!             POTFAK=POTFAK*(1-COS(SQRT(GSQU)*TPI*HFRCUT)*EXP(-(SQRT(GSQU)*TPI*HFRCUT)**2*HFRCUT_SMOOTH))
!test
             CALL DELSTP(3,SQRT(GSQU)*TPI*HFRCUT/10,DFUN,SFUN)
             POTFAK=POTFAK*(SFUN-0.5)*2
!test
          ELSE IF (LSR) THEN
             ! use the short range part of the Coulomb kernel
             POTFAK=POTFAK*(1._q-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
          ELSE IF (LRHFCALC.AND.LRSCOR) THEN
             ! use the long range part of the Coulomb kernel
             POTFAK=POTFAK*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
          ENDIF
       ENDIF

       ! smooth cutoff function between  ENCUTSOFT and ENCUT
       E=HSQDTM*(GSQU*TPI**2)
       IF (ENCUT>=0 .AND. E>ENCUT) THEN
          POTFAK=0
       ELSE

          MAXINDEX=MAXINDEX+1
          OLD_INDEX(MAXINDEX)=NI
          
          IF (ENCUTSOFT>=0 .AND. E>ENCUTSOFT) THEN

          ! cosine window
             IF ( .NOT. LSCK) THEN
                POTFAK=POTFAK*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
          ! squeezed Coulomb kernel
             ELSE 
                POTFAK=POTFAK*SQUEEZED_COULOMB_KERNEL(SQRT(GSQU)*TPI, Q1, Q2 )
             ENDIF
          ENDIF

       ENDIF

! POTFAK in fock has an additional factor. Since the corrections were originally devised
! for the fock routine this ensures consistency
       IF(MCALPHA/=0) THEN
          POTFAK=POTFAK*(1.0_q/WGWQ%GRID%NPLWV)
       ENDIF

       DATAKE(NI)=SQRT(POTFAK)
       ! maximum kinetic energy
       IF(DATAKE(NI)>DATAKEMAX) THEN
          DATAKEMAX=DATAKE(NI)
       ENDIF
    ENDDO

    ! setup of multipole corrections
    IF (MCALPHA/=0) THEN
       CALL FOCK_MULTIPOLE_CORR_SETUP(LATT_CUR, GRIDHF)
       CALL GW_MULTIPOLE_PROJ_SETUP( WGWQ, LATT_CUR, DATAKE, DATAKEMAX )
    ENDIF

    PROFILING_STOP( 'setup_coulomb_diagonal' )
    
  END SUBROUTINE SETUP_COULOMB_DIAGONAL

  !>Perform implicit basis set extrapolation
  !>Coulomb kernel is enhanced at intermediate momenta to make up for cutoff at Q2
  FUNCTION SQUEEZED_COULOMB_KERNEL( X, Q1, Q2)
    USE prec
    IMPLICIT NONE
  
    REAL(q) :: SQUEEZED_COULOMB_KERNEL, X, Q1, Q2
  
    SQUEEZED_COULOMB_KERNEL= ((Q1 - Q2)*(X - Q2))/(Q1**2 + X*(-2*Q1 + Q2))**2 *X**2
        
  END FUNCTION SQUEEZED_COULOMB_KERNEL

!
!> electronic structure factor routines 
!

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> allocates electronic structure factor handle 
   !---------------------------------------------------------------------------
   SUBROUTINE ALLOCATE_ELSTUFAC_HANDLE( ESF, WGW, ENCUTGW, ENCUTGWSOFT, LATT_CUR, SYMM, WDES, S2E, T_INFO,&
      DYN, IO, NELECT ) 
      USE constant
      USE mathtools, ONLY: INVERT_REAL_MATRIX 
      USE full_kpoints, ONLY: CHECK_FULL_KPOINTS
      USE mkpoints , ONLY: HNFORM

      USE msymmetry 
      USE chi_glb, ONLY: LESF_SPLINES , screened_2e_handle 
      USE pead, ONLY : LPEAD_SYM_RED
      USE kpoints_change, ONLY: GENERATE_KPOINTS_TRANS
      USE full_kpoints, ONLY: IBZKPT_HF
      USE screened_2e, ONLY: screened_2e_handle
      USE fock_glb, ONLY: NKREDX, NKREDY, NKREDZ
      USE fock, ONLY: SKIP_THIS_KPOINT_IN_FOCK

      TYPE (elstufac_handle)          :: ESF !< electronic structure factor handle
      TYPE (wavedes), TARGET          :: WGW !< response function descriptor 
      REAL(q), INTENT(IN)             :: ENCUTGW !< cutoff for response function and coulomb potential 
      REAL(q), INTENT(IN)             :: ENCUTGWSOFT !< soft cutoff for coulomb potential
      TYPE (latt)                     :: LATT_CUR  !< lattice information 
      TYPE (symmetry)                 :: SYMM !< symmetry 
      TYPE (wavedes)                  :: WDES !< wavefunction descriptor 
      TYPE (screened_2e_handle), INTENT(IN) :: S2E
      TYPE (type_info)                    :: T_INFO !< type information
      TYPE (dynamics)                     :: DYN !< MD handle
      TYPE (in_struct)                    :: IO !< ouput  handle 
      REAL(q), INTENT(IN), OPTIONAL       :: NELECT !< number of electrons in unit cell 
      !local 
      INTEGER                         :: NUMBER_OF_NQ !< number of q-points
      INTEGER :: I,J,L
      INTEGER :: NKPTS_ORIG
      REAL(q) :: B(3,3)
      REAL(q) :: BNORM(3),K(3)
      INTEGER :: N,NK, NI(3)
      INTEGER, PARAMETER :: NTOT = 100
      INTEGER :: IHNF(3,3)
      TYPE (skpoints_full)      :: KPOINTS_I2F
      INTEGER :: NG
      INTEGER :: NKRED(3)
      TYPE (wavedes1)           :: WGWK

      IF ( ESF%NKPX > 0 .AND. ESF%NKPY >0 .AND. ESF%NKPZ > 0 ) THEN
         CALL VTUTOR%BUG( "ESF already associated" , __FILE__, __LINE__ ) 
      ENDIF

      NKRED(1) = NKREDX
      NKRED(2) = NKREDY
      NKRED(3) = NKREDZ
      NUMBER_OF_NQ = S2E%NUMBER_OF_NQ
      ALLOCATE( ESF%NQ(S2E%NUMBER_OF_NQ))
      DO NK = 1, S2E%NUMBER_OF_NQ 
         ESF%NQ( NK ) = S2E%NQ( NK ) 
      ENDDO

      CALL CHECK_FULL_KPOINTS

      !scale reciprocal grid such that generating k-point grid is a subset 
      NI = 1 
      DO NK = 1 , KPOINTS_FULL_ORIG%NKPTS
         K = KPOINTS_FULL_ORIG%VKPT(:,NK)*NKRED(:)
         DO I = 1, 3 
            IF ( ABS( K(I) ) > 1.E-8_q ) THEN
               ! find first integer J such that K(1) * J = 1
               N = NI(I) 
               num: DO L = 1, NTOT
                  DO J = 1, NTOT
                     IF( ABS( ABS( K(I) * J ) - L ) < 1.E-8_q ) EXIT num
                  ENDDO
               ENDDO num
               IF ( J> NTOT .AND. L> NTOT) THEN
                  CALL VTUTOR%BUG( "k-point multiple not found "//str([ NK, I] )//" "//str([ K(I) ] ) , __FILE__, __LINE__ )
               ENDIF
               NI(I)  = MAX( J, N )
            ENDIF
         ENDDO
      ENDDO
      IF (IO%NWRITE > 2 ) THEN
         IF( IO%IU0>=0) WRITE(IO%IU0,'(A,3I6)') "effective k-point sampling found:", NI
         IF( IO%IU6>=0) WRITE(IO%IU6,'(A,3I6)') "effective k-point sampling found:", NI
      ENDIF

      ! store effective k-point sampling 
      ESF%NKPX=NI(1)
      ESF%NKPY=NI(2)
      ESF%NKPZ=NI(3)
      B(:,1) = LATT_CUR%B(:,1) / ESF%NKPX
      B(:,2) = LATT_CUR%B(:,2) / ESF%NKPY
      B(:,3) = LATT_CUR%B(:,3) / ESF%NKPZ
      BNORM(1) = SQRT( DOT_PRODUCT( B(:,1),B(:,1) ))
      BNORM(2) = SQRT( DOT_PRODUCT( B(:,2),B(:,2) ))
      BNORM(3) = SQRT( DOT_PRODUCT( B(:,3),B(:,3) ))
      
      IF( ESF%NKPX < 1 .OR. ESF%NKPY < 1 .OR. ESF%NKPZ <1 ) THEN
          CALL vtutor%bug( "number of k-points inconsistent: "//str([ ESF%NKPX, ESF%NKPY, ESF%NKPZ ] ), __FILE__,__LINE__ )
      ENDIF

      ! determine symmetry operations to blow up diagonal of 
      ! response function diagonal from irz to full BZ        
      NULLIFY( ESF%WGW )
  
      ! point to wavefunction descriptor 
      ESF%WGW => WGW 
      ! set number of q-points
      ESF%NUMBER_OF_NQ = NUMBER_OF_NQ 
      ALLOCATE( ESF%NG( ESF%NUMBER_OF_NQ ) )
      ESF%NG = -1 
      ! symmetry info 
      ESF%ISYM = -1

#ifndef gammareal
      ESF%ISYM = SYMM%ISYM
      IF (SYMM%ISYM>=0 ) THEN
         CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,DYN%EFOR,T_INFO%LSFOR, &
              T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
              SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
              SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)

         ! Since only a reduced set S_q(G) = diag( S_q( G,G') ) is calculated 
         ! one requires a map that maps S_q(G) to S_k(G) for k in the full BZ
         ! this is done by the following call to GENERATE_KPOINTS_TRANS routine.
         ! However, care must be taken for the RPAgamma case, where a reduced "full" k-point grid is used;
         ! one that contains typically half the number of the full k-points grid as KPOINTS_FULL_ORIG%NKPTS
         ! this reduced number is stored in ESF%WGW%NKPTS when run via ALGO = RPAR, GWR, etc. see chi_super.F
         ! to index a map for q -> k for k in the entire BZ, overwrite this variable, but make a backup
         ! of the original value, since ESF%WGW is pointer to WGW that is used later on. 
         NKPTS_ORIG = ESF%WGW%NKPTS
         ESF%WGW%NKPTS= KPOINTS_FULL_ORIG%NKPTS
         CALL GENERATE_KPOINTS_TRANS(ESF%WGW%GRID, ESF%NUMBER_OF_NQ, ESF%WGW, KPOINTS_FULL_ORIG, &
            ESF%KPOINTS_TRANS)
         ESF%WGW%NKPTS = NKPTS_ORIG 

!TODO probably remove this call 
         ! first map from irz to half BZ 
         CALL IBZKPT_HF(LATT_CUR, KPOINTS_ORIG, &
              KPOINTS_I2F, T_INFO%NIONS, SYMM%ROTMAP, &
              SYMM%MAGROT,SYMM%ISYM, IO%IU6, IO%IU0)

      ENDIF 
#endif 

      ! allocate storage for electronic structure factor in q+G domain
      NG=0 
      DO J = 1, ESF%WGW%NKPTS
         K = KPOINTS_FULL_ORIG%VKPT(:,J)
         IF( SKIP_THIS_KPOINT_IN_FOCK( K ) ) CYCLE 
         NG = MAX( NG, ESF%WGW%NGVECTOR(J) )
      ENDDO
      ALLOCATE( ESF%S( m_ NG, ESF%NUMBER_OF_NQ ) )
      ESF%S = 0 
      ! store number of electrons, optionally
      ! required for exhange energy
      ! implicitly, this also decides whether ESF handle is set up for Fock-exchange 
      ! and drops allocation of Coulomb potential
      ESF%NELECT = -1._q
      IF( PRESENT( NELECT ) ) THEN
         ESF%NELECT = NELECT/ WDES%RSPIN  
      ELSE
         ALLOCATE( ESF%OLD_INDEX( m_ NG ) )
         ALLOCATE( ESF%V( m_ NG ) )
         ESF%V=0 
      ENDIF

      ! store cutoffs 
      ESF%ENCUT = ENCUTGW
      ESF%ENCUTSOFT = ENCUTGWSOFT

      ! determine the number of q+G vectors for each direction 
      ! this follows from E_cut = |q_i+G_i|^2 / (2 m ) , i = x, y, z
      ! there are in total ( E_cut / |b_i| 2 m  + 1 )  G_i-points  (+1 to be save )
      ! also, there are NKPI q-points, which results in following formula 
      ! for the number of grid points in each direction 
      ESF%N1 = CEILING(SQRT(ENCUTGW/RYTOEV)/(TPI*AUTOA/LATT_CUR%ANORM(1)) + 1)*ESF%NKPX
      ESF%N2 = CEILING(SQRT(ENCUTGW/RYTOEV)/(TPI*AUTOA/LATT_CUR%ANORM(2)) + 1)*ESF%NKPY
      ESF%N3 = CEILING(SQRT(ENCUTGW/RYTOEV)/(TPI*AUTOA/LATT_CUR%ANORM(3)) + 1)*ESF%NKPZ
      ! WRITE(*,*) ESF%WGW%GRID%NGPTAR, ESF%N1, ESF%N2,ESF%N3
      ! WRITE(*,*) NINT((SQRT(COR%ENCUTGW(1)/RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))+1)*KPOINTS_%NKPZ)
      ESF%OMEGA=LATT_CUR%OMEGA

      ! scaled reciprocal bravais matrix, scaled to k-point sampling 
      ESF%B(:,:)=B

      ! determine real-space Bravais lattice from reciprocal one 
      ESF%A(:,:)=ESF%B(:,:)
      CALL INVERT_REAL_MATRIX( ESF%A , -1 )

      ESF%LATT_CUR=LATT_CUR

      ESF%ESPLINED = 0
      ESF%ISTATUS = -1

! test if k-point map is set up correctly 
#ifndef gammareal
      IF (  ESF%ISYM >= 0 ) THEN
      IF ( SIZE( ESF%KPOINTS_TRANS%NK_OLD ) < KPOINTS_FULL_ORIG%NKPTS ) THEN
         CALL VTUTOR%BUG( " k-point grids are inconsistent: "//str([ SIZE(ESF%KPOINTS_TRANS%NK_OLD), KPOINTS_FULL_ORIG%NKPTS ]), __FILE__, __LINE__ )
      ENDIF
      ENDIF
      ! k-point map test
      ! loop over irreducible brillouin zone points
      DO I=1,ESF%NUMBER_OF_NQ
         ! loop over full brillouin zone points
         DO J=1,KPOINTS_FULL_ORIG%NKPTS
            K = KPOINTS_FULL_ORIG%VKPT(:,J)
            IF( SKIP_THIS_KPOINT_IN_FOCK( K ) ) CYCLE 
            CALL SETWDES(ESF%WGW, WGWK, J )
            IF( ESF%ISYM>=0 ) THEN
               ! loop over symmetry - equivalent points
               IF ( ESF%KPOINTS_TRANS%NK_OLD( J ) /= ESF%NQ(I) ) CYCLE
            ENDIF 
         ENDDO
      ENDDO
#endif
      ! check transformation of grid from carthesian to FFT and back, only once required 
      CALL GRID_TRAFO_CHECK( ESF ) 
   END SUBROUTINE ALLOCATE_ELSTUFAC_HANDLE

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> deallocates electronic structure factor handle 
   !---------------------------------------------------------------------------
   SUBROUTINE DEALLOCATE_ELSTUFAC_HANDLE( ESF )
      TYPE (elstufac_handle)          :: ESF !< electronic structure factor handle

      IF ( ESF%NKPX < 0 .OR. ESF%NKPY <0 .OR. ESF%NKPZ < 0 ) THEN
          CALL VTUTOR%BUG( "ESF is not associated" , __FILE__, __LINE__ ) 
      ENDIF

      ESF%NKPX=0
      ESF%NKPY=0
      ESF%NKPZ=0

      NULLIFY( ESF%WGW )
      ESF%NUMBER_OF_NQ = 0
      ESF%ISYM = -1
      DEALLOCATE( ESF%NG ) 

      DEALLOCATE( ESF%NQ ) 

      DEALLOCATE(ESF%S )

      ESF%ENCUT = 0
      ESF%ENCUTSOFT = 0

      ESF%N1 = 0
      ESF%N2 = 0
      ESF%N3 = 0

      ESF%OMEGA=0

      ! scaled reciprocal bravais matrix, scaled to k-point sampling 
      ESF%B(:,:)=0
      ! determine real-space Bravais lattice from reciprocal one 
      ESF%A(:,:)=0
      ESF%ESPLINED = 0
      ESF%ISTATUS = -1

      IF ( ALLOCATED( ESF%V ) ) DEALLOCATE( ESF%V )
      IF( ALLOCATED( ESF%OLD_INDEX ) ) DEALLOCATE( ESF%OLD_INDEX )

   END SUBROUTINE DEALLOCATE_ELSTUFAC_HANDLE
   
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> calculates electronic strcuture factor for a specific q-points
   !> and stores result into responsefunction structure 
   !---------------------------------------------------------------------------
   SUBROUTINE EXTRACT_ELSTUFAC( CHI_QG, NP, LREAL, WGWQ, LATT_CUR, ESF, &
      NQ, LAMB, OWEIGHT, ERPA, EDMP2, LSR, IU0 ) 
      GDEF, INTENT(INOUT)      :: CHI_QG(:,:) !< current chi at q+G
      INTEGER, INTENT(IN)      :: NP          !< max dimension used of CHI_QG
      LOGICAL,INTENT(IN)       :: LREAL       !< is gamma point 
      TYPE( wavedes1 )         :: WGWQ        !< information about q
      TYPE (latt)              :: LATT_CUR    ! lattice information 
      TYPE( elstufac_handle) :: ESF 
      INTEGER,INTENT(IN)       :: NQ          !< current q-point 
      REAL(q),INTENT(IN)       :: LAMB        !< coupling parameter
      REAL(q),INTENT(IN)       :: OWEIGHT     !< omega weight 
      COMPLEX(q),INTENT(INOUT) :: ERPA
      COMPLEX(q),INTENT(INOUT) :: EDMP2
      LOGICAL,INTENT(IN)       :: LSR         !< short range coulomb potential
      INTEGER,INTENT(IN)       :: IU0 ! output
      !local 
      REAL(q)                  :: HFEIG(NP)
      COMPLEX(q)               :: TRACE
      INTEGER                  :: N 
      INTEGER                  :: NI
      CHARACTER(LEN=16)        :: APP
      CHARACTER(LEN=32)        :: FM

      ! diagonalize V.chi.V and obtain ( ln(1-v.chi.v) - v.chi.v )
      ! where v = V^{-1/2}
      TRACE = 0 
      CALL ROTLN( CHI_QG, NP, LREAL, LAMB, 0, TRACE, -1, HFEIG )

#ifdef debug
      !dump diagonal for each value of omega
      WRITE(98,'(2I4,3F20.10)')  NQ, NP, OWEIGHT, REAL(ERPA,q),REAL(EDMP2,q)
      DO NI = 1, NP
         WRITE(98,'(I4,4F20.10)') NI, HFEIG( NI ), CHI_QG(NI,NI)
      ENDDO
#endif

      DO NI=1,NP
         ! divide by 
         IF ( ABS(ESF%V( ESF%OLD_INDEX( NI) )) > 1E-8_q ) THEN
             CHI_QG(1:NP,NI)=CHI_QG(1:NP,NI)/ ESF%V(ESF%OLD_INDEX( NI ))
             CHI_QG(NI,1:NP)=CHI_QG(NI,1:NP)/ ESF%V(ESF%OLD_INDEX( NI ))
         ELSE
             CHI_QG(1:NP,NI)=CHI_QG(1:NP,NI)/MAX(1E-8_q,ESF%V(ESF%OLD_INDEX( NI )))
             CHI_QG(NI,1:NP)=CHI_QG(NI,1:NP)/MAX(1E-8_q,ESF%V(ESF%OLD_INDEX( NI )))

         ENDIF
      ENDDO

!#define unitest
#ifdef unitest
      !unit test
      TRACE=0
      DO NI=1,NP
            TRACE = TRACE + CHI_QG( NI, NI) * ESF%V( ESF%OLD_INDEX( NI ) ) ** 2
      ENDDO
      ERPA = ERPA+ TRACE
!WRITE(80,'(I4,3F20.10)') NQ, OWEIGHT, TRACE
#else
      ! obtain contributions to RPA and MP2 energies
      CALL TR_EIGENVALUES( ERPA, EDMP2, HFEIG, NP )
#endif

      ! potentially apply multipole corrections
      !CALL APPLY_MULTIPOLE_CORRECTIONS( CHI_QG, NP, WGWQ%GRID%NPLWV  )
      !DO NI=1,NP
      !   ! get v^ 1/2
      !   CHI_QG(1:NP,NI)=CHI_QG(1:NP,COR%ESF(I)%OLD_INDEX(NI))
      !ENDDO

      !DO NI=1,NP
      !   ! get v^ 1/2
      !   CHI_QG(NI,1:NP)=CHI_QG(COR%ESF(I)%OLD_INDEX(NI),1:NP)
      !ENDDO

      ! store diagonal of electronic structure factor for highest cutoff 
!WRITE(*,*) "S ",NQ, OWEIGHT, ESF%NQ(NQ) , SIZE( ESF%S, 2), NP
      DO  N = 1, NP 
         ESF%S( N, NQ ) = ESF%S( N, NQ ) + CHI_QG( N, N )*OWEIGHT
!WRITE(*,'(I3,3X,40F10.4)') N, (REAL(CHI_QG(N,NI),q), NI=1, NP)
      ENDDO
      ! store number of G vectors for current q-point 
      ESF%NG( NQ )  = NP 

      CONTAINS 

      !--------------------------------------------------------------------------
      !
      ! computes direct MP2 and RPA contribution for current q-point
      !
      !--------------------------------------------------------------------------
      SUBROUTINE TR_EIGENVALUES( ERPA, EDMP2, HFEIG, NP )
         USE chi_glb , ONLY: LLTDMP2
         COMPLEX(q)               :: ERPA
         COMPLEX(q)               :: EDMP2
         REAL(q)                  :: HFEIG(:)
         INTEGER                  :: NP 
         ! local 
         INTEGER                  :: N
         
         !=======================================================================
         ! sum over eigenvalues
         !=======================================================================
         DO N=1,NP
            IF (1-HFEIG(N)<1E-8 .AND. .NOT. LLTDMP2 ) THEN
               WRITE(*,*) 'TR_EIGENVALUES reports '//&
                  'negative or close to zero eigenvalue', HFEIG(N), N
            ELSE
               ! skip RPA correlation energy for laplace transformed direct MP2
               IF ( .NOT. LLTDMP2 ) &
               ERPA=ERPA+((LOG(1-HFEIG(N)*LAMBDA)+HFEIG(N)*LAMBDA))
               IF ( LDMP1 ) THEN
                  EDMP2=EDMP2-HFEIG(N)
               ELSE
                  EDMP2=EDMP2+(-HFEIG(N)*HFEIG(N))/2
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE TR_EIGENVALUES

      SUBROUTINE DUMP_STUFAK( VKPT, VQPT, S, V, NP_NEW, CPHASE, NINDPW, LINV, LSHIFT ,INDX, IU) 
         REAL(q),INTENT(IN)     :: VKPT( 3 ) , VQPT(3)
         GDEF, INTENT(IN)       :: S(:,:)
         REAL(q), INTENT(IN)    :: V(:)
         INTEGER, INTENT(IN)    :: NP_NEW
         INTEGER, INTENT(IN)    :: NINDPW(:)
         INTEGER, INTENT(IN)    :: INDX(:)
         COMPLEX(q), INTENT(IN) :: CPHASE(:)
         LOGICAL, INTENT(IN)    :: LINV, LSHIFT
         INTEGER                :: IU

         CHARACTER( LEN=32 )    :: APP, FM
         INTEGER :: N, M, N_

         WRITE(APP,*) NP_NEW
         FM='(I6, 2F12.6, "|", '//TRIM(ADJUSTL(APP))//"F12.6)"
         WRITE(IU,'(A,I4,3F8.3,"|",3F8.3)') " REAL " , NP_NEW, VKPT(:), VQPT(:)
!         WRITE(IU,'(A,13X,22I12)') " NINDPW" , NINDPW(:)
         WRITE(IU,'(A,13X,'//TRIM(ADJUSTL(APP))//'I12)') " INDX  " , INDX(1:NP_NEW)
         WRITE(IU,'(A,13X,'//TRIM(ADJUSTL(APP))//'F12.6)') " V  " , V(1:NP_NEW)
!         WRITE(IU,'(A,13X,22F12.6)') "ESFV" , COR%ESF(I)%V(:)
!         IF( ASSOCIATED( COR%ROTMAP ) ) THEN
!            WRITE(IU,'(A,13X,22I12)') " NIND I" , ( NINDPW(INDX(N)),N=1,NP_NEW)
!            WRITE(IU,'(A,13X,22I12)') " I NIND" , ( INDX(NINDPW(N)),N=1,NP_NEW)
!         ENDIF
         DO N = 1, NP_NEW
            !WRITE(IU,TRIM(ADJUSTL(FM)))N, V(N), V(INDX(N)), REAL( S(N, 1:NP_NEW ) , q)
            N_=N
IF( WGWQ%LGAMMA) N_ =(N-1)/2+1
            WRITE(IU,TRIM(ADJUSTL(FM))) N_, V(N_), V(INDX(N)),(  REAL( S( INDX(N), INDX(M) ) , q),M=1,NP_NEW )
         ENDDO 

         WRITE(IU,*) " UNSORTED " 
         DO N = 1, NP_NEW
            !WRITE(IU,TRIM(ADJUSTL(FM)))N, V(N), V(INDX(N)), REAL( S(N, 1:NP_NEW ) , q)
            N_=N
IF( WGWQ%LGAMMA) N_ =(N-1)/2+1
            WRITE(IU,TRIM(ADJUSTL(FM))) N_, V(N_), V(INDX(N)),(  REAL( S( N, M ) , q),M=1,NP_NEW )
         ENDDO 
      END SUBROUTINE DUMP_STUFAK

   END SUBROUTINE EXTRACT_ELSTUFAC

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> gathers electronic strcuture factor for a specific q-points
   !> and stores result into responsefunction structure, used in low-scaling 
   !> RPA algorithm, called in COL2BC_AND_RPA (non-GPU version)
   !> @param[in]  NP_NEW   dimension of X.V for current cutoff NE
   !> @param[in]  EV       eigenvectors of V^1/2 .X. V^1/2 stored in distributed fashion according RDES2D
   !> @param[in]  E        eigenvalues of V^1/2 .X. V^1/2 stored in distributed fashion according RDES2D
   !> @param[in]  RDES2D   distributed response descriptor for block-cyclic 2D processor grid 
   !> @param[in]  NQ       currently treated q-point
   !> @param[in]  OWEIGHT  frequency integration weight
   !> @param[out] ERPA   RPA energy contribution of current q-point that is added to total RPA correlation energy
   !> @param[out] EMP2   MP2 energy contribution of current q-point that is added to total RPA correlation energy
   !> @param[inout] ESF  electronic structure factor handle where S_G(q) is accumulated on exit
   !---------------------------------------------------------------------------
#ifdef scaLAPACK
   SUBROUTINE EXTRACT_STUFAK_DISTRI( NP_NEW, EV, E, RDES2D, NQ, OWEIGHT,&
       ERPA, EMP2, ESF )
       USE scala 
       USE ini 
       INTEGER, INTENT(IN)          :: NP_NEW 
       GDEF                         :: EV(:)
       REAL(q), INTENT(IN)          :: E(:)
       TYPE(greens_mat_des)         :: RDES2D 
       INTEGER, INTENT(IN)          :: NQ
       REAL(q), INTENT(IN)          :: OWEIGHT
       REAL(q), INTENT(INOUT)       :: ERPA
       REAL(q), INTENT(INOUT)       :: EMP2
       TYPE(elstufac_handle)        :: ESF 
       !> local variables
       INTEGER    :: I 
       REAL(q)    :: LOGE( NP_NEW ) 
       !> require 2 more auxillary matrices 
       GDEF, ALLOCATABLE   :: GWORK1(:)
       GDEF, ALLOCATABLE   :: GWORK2(:)
       INTEGER             :: ISTAT 
       INTEGER             :: IU0 = -1

       IF( RDES2D%COMM%NODE_ME==1 ) IU0=6
       
       PROFILING_START( 'extract_stufak_distri' )

       EMP2 = 0 
       ERPA = 0 

       ! store number of G vectors for current q-point 
       ESF%NG( NQ )  = NP_NEW

       ! allocate auxillary storage 
       ALLOCATE(GWORK2(RDES2D%MY_NROWS*RDES2D%MY_NCOLS), STAT = ISTAT) 
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "COL2BC_AND_RPA is not able to allocate "//&
            str( ndata*8._q*RDES2D%MY_NROWS*RDES2D%MY_NCOLS) //&
            " kB of data on MPI rank"//str( RDES2D%COMM%NODE_ME ) )
       ENDIF
       CALL REGISTER_ALLOCATE(ndata*8._q*SIZE(GWORK2,KIND=qi8), "esf")

       ! allocate auxillary storage 
       ALLOCATE(GWORK1(RDES2D%MY_NROWS*RDES2D%MY_NCOLS), STAT = ISTAT) 
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "COL2BC_AND_RPA is not able to allocate "//&
            str( ndata*8._q*RDES2D%MY_NROWS*RDES2D%MY_NCOLS) //&
            " kB of data on MPI rank"//str( RDES2D%COMM%NODE_ME ) )
       ENDIF
       CALL REGISTER_ALLOCATE(ndata*8._q*SIZE(GWORK1,KIND=qi8), "esf")

       ! determine MP2 energy 
       DO I = 1, NP_NEW
          IF( LDMP1 ) THEN
             EMP2 = EMP2 + E(I)*OWEIGHT
          ELSE
             EMP2 = EMP2 + E(I)*E(I)/2*OWEIGHT
          ENDIF
       ENDDO
   
       ! remainer of this routine calcualtes the electronic structure factor S_G(q)
   
       ! obtain RPA formula 
       DO I = 1, NP_NEW
           IF( 1 + E( I ) > 0 ) THEN 
              LOGE(I) = ( LOG( 1 + E( I ) ) - E( I ) ) 
           ELSE IF ( ABS( OWEIGHT ) > 1.E-8_q .AND. 1+E(I) <0 ) THEN
              WRITE(*,*)' eigenvalue of CHI.V is not positive ', NQ,I, E(I)
           ENDIF
       ENDDO
!IF( IU0>= 0 ) WRITE(*,*) "LOG (1+E)-E"
!IF( IU0>= 0 ) WRITE(*,'(20F10.4)') LOGE

!#ifdef oldversion
       GWORK1 = 0 
       CALL ADD_TO_DIAGONALE_GDEF( NP_NEW, GWORK1, LOGE, RDES2D%DESC )
       ! GWORK1 <-(  U * (ln( 1-e)-e) ) U^+
       CALL PGGEMM( 'N', 'N', NP_NEW, NP_NEW, NP_NEW, one,  &
            EV(1), 1, 1, RDES2D%DESC, &
            GWORK1(1), 1, 1, RDES2D%DESC, &
            zero,  &
            GWORK2(1), 1, 1, RDES2D%DESC )
   !CALL DUMP_MAT( "U*_ln  ", GWORK2, NP_NEW, RDES2D%DESC, RDES2D%COMM_INTER, IO%IU0 )
!#else
!       GWORK2 = EV
!       CALL RIGHT_MULTIPLY_BY_R( NP_NEW, GWORK2, LOGE, RDES2D%DESC )
!#endif
       CALL PGGEMM( 'N', trans, NP_NEW, NP_NEW, NP_NEW, one,  &
            GWORK2(1), 1, 1, RDES2D%DESC, &
            EV(1), 1, 1, RDES2D%DESC, &
            zero,  &
            GWORK1(1), 1, 1, RDES2D%DESC )
!   CALL DUMP_MAT( "U*_ln_U", GWORK1, NP_NEW, RDES2D%DESC, RDES2D%COMM_INTER, IU0 )
   
#ifdef unit_test
       ! extract diagonal 
       LOGE=0
       CALL DETERMINE_DIAGONALE_GDEF_REAL ( NP_NEW, GWORK1, LOGE, RDES2D%DESC ) 
 WRITE(300+RDES2D%COMM%NODE_ME,'(A,3I4,F10.4)') "RPA", RDES2D%MY_NROWS,RDES2D%MY_NCOLS,NP_NEW, ERPA
 WRITE(300+RDES2D%COMM%NODE_ME,'(20F10.4)') LOGE
       CALLMPI( M_sum_d( RDES2D%COMM_INTRA, LOGE(1), NP_NEW ) )
       ERPA = 0 
       DO I = 1, NP_NEW
          ERPA=ERPA + LOGE(I)
       ENDDO
 WRITE(400+RDES2D%COMM%NODE_ME,'(A,3I4,F10.4)') "RPA", RDES2D%MY_NROWS,RDES2D%MY_NCOLS,NP_NEW, ERPA
 WRITE(400+RDES2D%COMM%NODE_ME,'(20F10.4)') LOGE
STOP
   PROFILING_STOP( 'extract_stufak_distri' )
   RETURN
#endif

       ! GWORK1 contains (  U * (ln( 1-e)-e) ) U^+
       ! multiply with Coulomb potential from left and right 
       LOGE = 0 
       ! setup coulomb potential 
       DO I = 1, NP_NEW
          IF ( ABS(ESF%V( ESF%OLD_INDEX(I) )) > 1E-8_q ) THEN
             LOGE(I) = 1._q/ESF%V( ESF%OLD_INDEX(I) )
          ELSE
             LOGE(I) = 1E8_q
          !   WRITE(*,*) "WARNING", NQ, I, ESF%V( ESF%OLD_INDEX(I) ) 
          ENDIF
       ENDDO

#ifdef oldversion
       EV = 0 
       CALL ADD_TO_DIAGONALE_GDEF( NP_NEW, EV, LOGE, RDES2D%DESC )
       CALL PGGEMM( 'N', 'N', NP_NEW, NP_NEW, NP_NEW, one,  &
            EV(1), 1, 1, RDES2D%DESC, &
            GWORK1(1), 1, 1, RDES2D%DESC, &
            zero,  &
            GWORK2(1), 1, 1, RDES2D%DESC )
       CALL PGGEMM( 'N', 'N', NP_NEW, NP_NEW, NP_NEW, one,  &
            GWORK2(1), 1, 1, RDES2D%DESC, &
            EV(1), 1, 1, RDES2D%DESC, &
            zero,  &
            GWORK1(1), 1, 1, RDES2D%DESC )
   !CALL DUMP_MAT( "S ", GWORK1, NP_NEW, RDES2D%DESC, RDES2D%COMM_INTER, IO%IU0 )
#else
       CALL MULTIPLY_A_BY_CELTOTINV_GDEF( NP_NEW, GWORK1, LOGE, RDES2D%DESC ) 
#endif

       ! extract electronic structure factor
       LOGE = 0 
       CALL DETERMINE_DIAGONALE_GDEF_REAL ( NP_NEW, GWORK1, LOGE, RDES2D%DESC ) 
       CALLMPI( M_sum_d( RDES2D%COMM_INTRA, LOGE(1), NP_NEW ) )
   
       ! contract with coulomb potential, this should give the RPA correlation energy
!IF( IU0>=0 ) &
!WRITE(81,'(A,I4, 2F12.6)') "#", NP_NEW, ESF%ENCUT, OWEIGHT
       DO I = 1, NP_NEW
!IF( IU0>=0 ) &
!WRITE(81,'(I4,3F20.10)') I, E(I), LOGE(I), ESF%V( ESF%OLD_INDEX(I) )
          ERPA=ERPA + ( ESF%V( ESF%OLD_INDEX(I) )**2)  * LOGE(I) * OWEIGHT
          ! accumulate electronic structure factor for this q-point 
          ESF%S( I, NQ ) = ESF%S(I, NQ ) + LOGE(I)*OWEIGHT/2
!IF( IU0>=0 ) &
!WRITE(81,'(2I4,2F20.10)') NQ, I, ESF%V( ESF%OLD_INDEX(I) ) 
       ENDDO

       ! free auxillary storage
       CALL DEREGISTER_ALLOCATE(ndata*8._q*SIZE(GWORK1,KIND=qi8), "esf")
       DEALLOCATE(GWORK1)
       CALL DEREGISTER_ALLOCATE(ndata*8._q*SIZE(GWORK2,KIND=qi8), "esf")
       DEALLOCATE(GWORK2)

       PROFILING_STOP( 'extract_stufak_distri' )
   
       CONTAINS
   
       ! helper routine that dumps a distributed matrix to stdout or file 
       SUBROUTINE DUMP_MAT( STRING, CHAM_DISTRI, NB_TOT, DESC, COMM_INTRA, IU_,IUNIT)
         USE wave
         USE scala
         IMPLICIT NONE
         CHARACTER (LEN=*) :: STRING            ! string to dump
         GDEF              :: CHAM_DISTRI(*)    ! distributed matrix
         INTEGER           :: DESC(*)
         TYPE(communic)    :: COMM_INTRA 
         INTEGER           :: NB_TOT
         INTEGER           :: IU_               ! unit to write to (not dump for IU<0)
         INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
       ! local
         INTEGER :: NDUMP=40
         GDEF,ALLOCATABLE   :: CHAM(:,:)
         INTEGER IU
         INTEGER N1, N2
         INTEGER COLUMN_HIGH, COLUMN_LOW
         LOGICAL :: LGAMMA_PRINT = .FALSE.
     
         IF ( PRESENT( IUNIT ) ) THEN 
            IU = IUNIT
            NDUMP =  NB_TOT
         ELSE
            IU = IU_
         ENDIF
         COLUMN_LOW=1
     
         COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)
     
         ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))
     
         CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  DESC, COLUMN_LOW, COLUMN_HIGH)
     
         CALLMPI( M_sum_g(COMM_INTRA, CHAM(1,1), SIZE(CHAM)))
     
     
         IF (IU>=0) THEN
            WRITE(IU,*) STRING
            DO N1=1,NDUMP
               IF ( PRESENT( IUNIT ) ) THEN 
                  DO N2 = 1, NDUMP
                     WRITE(IU,10)N1+COLUMN_LOW-1, N2+COLUMN_LOW-1, CHAM(N1+COLUMN_LOW-1,N2)
                  ENDDO
               ELSE
                  WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1+COLUMN_LOW-1,N2) ,KIND=q) ,N2=1,NDUMP)
               ENDIF
            ENDDO
         ENDIF
     
         DEALLOCATE(CHAM)
1        FORMAT(1I3,3X,40F10.4)
2        FORMAT(1I3,3X,40F10.4)
10       FORMAT(2I4,2F20.10)
     
       END SUBROUTINE DUMP_MAT
   
   END SUBROUTINE EXTRACT_STUFAK_DISTRI
#endif

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> main workhorse, tricubic spline fit of electronic structure factor
   !> first S(G,q) is used to set up STUFAK(X,Y,Z), finally the spline 
   !> interpolation is performed
   !---------------------------------------------------------------------------
    SUBROUTINE ESF_SPLINE_FIT( ESF, LCUT, IO )
       USE constant
       USE ini
       USE mpimy, ONLY :LOCAL_INDEX_CYCLIC 
       USE lattice, ONLY: LATTIC
       TYPE( elstufac_handle ) :: ESF
       LOGICAL                 :: LCUT 
       TYPE( in_struct )       :: IO 
       !> local variables
       REAL(q)              :: INTCORR(ESF_NINTER)
       INTEGER              :: L, I, J, K
       REAL(q)              :: EDIFF
       INTEGER              :: SCALING(3)
       INTEGER              :: M1, M2, M3
       INTEGER              :: I_, J_, K_

       INTEGER              :: ISTAT
       REAL(q), ALLOCATABLE :: INTSTUFAK(:,:,:)

       REAL(q), ALLOCATABLE :: P(:,:)
       REAL(q), ALLOCATABLE :: SXYZ(:,:,:)
       REAL(q)              :: FINT, FITDER, Z, Y, X
       INTEGER              :: IU0 = -1
       ! coulomb potential
       REAL(q)              :: POTFAK
       REAL(q)              :: DK(3)
       REAL(q)              :: B(3,3)
       LOGICAL              :: LSR = .FALSE. 
       INTEGER              :: NODE_ME=1
       INTEGER              :: IONODE=1
       INTEGER              :: NCPU=1
       INTEGER              :: IND
       TYPE( latt )         :: LATT_EWALD
       REAL(q)              :: FSG_INT(ESF_NINTER )
! use this for slines (first derivative at q=0)
       REAL(q), PARAMETER   :: Y1P = 0.0_q
! use this for natural splines (first derivative at q /=0)
!       REAL(q), PARAMETER   :: Y1P = 10E10_q
!#ifdef debug
       INTEGER, SAVE :: ICOUNTER = -1
       ICOUNTER = ICOUNTER + 1
!#endif
       PROFILING_START( "esf_spline_fit" )

#ifdef MPI
       NODE_ME= ESF%WGW%COMM%NODE_ME
       IONODE = ESF%WGW%COMM%IONODE
       NCPU   = ESF%WGW%COMM%NCPU
#endif

       ESF%ISTATUS = 0

       B = ESF%B
       !> it might be better to interpolate in index space,
       !> (i,j,k) are the carthesian components of x_1, x_2, x_3
       !> with x_i being the scaled unit-vector e_i/L
#define CarthesianInterpolation
#ifdef CarthesianInterpolation
       ! setup scaled euclidian metric 
       B=0 
       B(1,1) = 1._q
       B(2,2) = 1._q
       B(3,3) = 1._q
       B(:,1) = B(:,1) / ESF%NKPX
       B(:,2) = B(:,2) / ESF%NKPY
       B(:,3) = B(:,3) / ESF%NKPZ
#endif

       !> output 
       IF( IO%NWRITE > 2 ) THEN
          IU0 = IO%IU6
!       ELSE
!          IF( IO%IU0>=0 ) WRITE(*,'(A,2F12.3)') "electronic structure interpolation for E_cut:  ",ESF%ENCUT, ESF%ENCUTSOFT
       ENDIF

       IF( IU0 >=0 ) THEN
          WRITE(IU0,'(A,3I4)') "structure factor k-points:  ",ESF%NKPX, ESF%NKPY,ESF%NKPZ
          WRITE(IU0,'(A,3I4)') "structure factor lattice:", ESF%N1, ESF%N2,ESF%N3
          WRITE(IU0 ,'(A,I20,A)') "max. RAM required : ", (2*ESF%N1+1)* (2*ESF%N2+1)* (2*ESF%N3+1)*8/1024, " kB"

          WRITE(*  ,'(A,3I4)') "structure factor k-points:  ",ESF%NKPX, ESF%NKPY,ESF%NKPZ
          WRITE(*  ,'(A,3I4)') "structure factor lattice:", ESF%N1, ESF%N2,ESF%N3
          WRITE(*  ,'(A,I20,A)') "max. RAM required : ", (2*ESF%N1+1)* (2*ESF%N2+1)* (2*ESF%N3+1)*8/1024, " kB"
       ENDIF
       ! allocate storage for setup of structure factor in carthesian coordinates 
       ALLOCATE( SXYZ( 2*ESF%N1+1, 2*ESF%N2+1, 2*ESF%N3+1 ) )
       SXYZ = 0 
       CALL REGISTER_ALLOCATE(8._q*SIZE(SXYZ,KIND=8), "esf")

       ! setup ESF in carthesian coordinates
       CALL SETUP_SXYZ( ESF, SXYZ, LCUT, IO )
#ifdef debug
IF( IU0>=0 ) THEN
IF( ESF%NELECT >0 ) THEN
 ICOUNTER = 1000
ENDIF
DO K = 1, ESF%N3* 2 + 1 
   K_ = K -1 - ESF%N3
   DO J = 1, ESF%N2* 2 + 1 
      J_ = J -1 - ESF%N2
      DO I = 1, ESF%N1* 2 + 1 
         I_ = I -1 - ESF%N1

         X = I_*    B(1,1) + J_*     B(1,2) + K_*    B(1,3)
         Y = I_*    B(2,1) + J_*     B(2,2) + K_*    B(2,3)
         Z = I_*    B(3,1) + J_*     B(3,2) + K_*    B(3,3)
         IF( I_ == 0 ) WRITE(700+ICOUNTER,'(9F12.6)') Y,Z, SXYZ( I, J, K )
         IF( J_ == 0 ) WRITE(800+ICOUNTER,'(9F12.6)') X,Z, SXYZ( I, J, K )
         IF( K_ == 0 ) WRITE(900+ICOUNTER,'(9F12.6)') X,Y, SXYZ( I, J, K )
                      WRITE(2000+ICOUNTER,'(9F12.6)') X,Y,Z,  SXYZ( I, J, K )
         WRITE(1000+ICOUNTER,'(3I5, 9F12.6)') I_, J_, K_, X,Y,Z, SXYZ( I, J, K )
      ENDDO
   ENDDO
ENDDO
ENDIF
IF( ESF%NELECT >0 ) THEN
STOP
ENDIF
#endif       

       INTCORR=0
       EDIFF = 0
       FSG_INT=0
       refinement: DO L=1, ESF_NINTER
         SCALING(1)=L
         SCALING(2)=L
         SCALING(3)=L
         M1=SCALING(1)*ESF%N1
         M2=SCALING(2)*ESF%N2
         M3=SCALING(3)*ESF%N3

         IF ( ESF%NELECT > 0 ) THEN
            ! estimate ewald energy for 
            LATT_EWALD%A(:,1)=(ESF%LATT_CUR%A(:,1)*L*ESF%NKPX)/NKREDX
            LATT_EWALD%A(:,2)=(ESF%LATT_CUR%A(:,2)*L*ESF%NKPY)/NKREDY
            LATT_EWALD%A(:,3)=(ESF%LATT_CUR%A(:,3)*L*ESF%NKPZ)/NKREDZ
            CALL LATTIC(LATT_EWALD)
            CALL EWALD_MONO(ESF%WGW%GRID,FSG_INT(L),1.0_q,LATT_EWALD)
            FSG_INT(L)=FSG_INT(L)*2 * ESF%NELECT 
         ENDIF
       
         ! allocate storage for interpolated structure factor 
         ALLOCATE(INTSTUFAK( 2*M3+1, 2*M2+1, 2*ESF%N1+1), STAT = ISTAT)
         IF (ISTAT/=0) THEN
            CALL vtutor%error( " cannot allocate "//str( SIZE( INTSTUFAK, KIND=8 )*8 ) &
               // " bytes of memory for electronic structure interpolation of L="//str([ L, ESF%N1, ESF%N2, ESF%N3]) )
         ENDIF
         INTSTUFAK = 0 
         CALL REGISTER_ALLOCATE(8._q*SIZE(INTSTUFAK,KIND=8), "esf" )
         ! spline interpolation into z-direction  
         DO I = -ESF%N1, ESF%N1
            ! determine index position in interpolated structure factor
            ! position of I in [-N1:N1] in refined array [1:2M1+1] is 
            ! 1 + ( I+ESF%N1 + 1 - 1 )*SCALING(1)
            ! because I__ in [ 1:2N1+1] is at position
            !           I__  = I + ESF%N1 + 1  
            ! and mapping I__ from [ 1:2N1+1] -> [ -M1:M1 ] is essentially
            !  1 + ( I__ - 1 ) * SCALING(1) 
            I_ = 1 + I + ESF%N1 
#ifdef MPI
            !> block this loop 
            IND = LOCAL_INDEX_CYCLIC( NCPU, NODE_ME, I_ ) 
            IF ( IND < 1 ) CYCLE 
#endif
            DO J = -ESF%N2, ESF%N2
               ! determine index position in interpolated structure factor
               J_ = 1 + (J+ESF%N2)*SCALING(2)
               ! tricubic splines: 
               ! z-direction first:
               ALLOCATE( P( 2*ESF%N3 +1 , 5 ) )
               DO K = -ESF%N3, ESF%N3
                  K_ = 1 + K+ ESF%N3 !< position in [1:2N3+1] is 
                  ! z-coordinate
                  P( K + ESF%N3 + 1, 1 ) = I* B(3,1) + J *B(3,2) + K * B(3,3)
                  ! function value S(x,0,0) 
                  P( K + ESF%N3 + 1, 2 ) =  SXYZ( 1+I+ESF%N1, 1+J+ESF%N2, 1 + K+ ESF%N3 )
               ENDDO
               ! obtain natural-spline coefficients  (those with f"(boundary) = 0 )
               CALL  SPLCOF( P, 2*ESF%N3+1 , 2*ESF%N3+1, Y1P )
               ! now interpolate spline fit and store 
               DO K=-M3, M3 
                  K_ = K + M3 + 1
                  Z = I*    B(3,1) + J*    B(3,2) + K*(     B(3,3)/SCALING(3) )
                  CALL SPLVAL( Z , FINT, FITDER, P(1,1), 2*ESF%N3+1, 2*ESF%N3+1 )
                  ! store result in proper slot of interpolated structure factor 
                  INTSTUFAK( K_, J_, I_ ) = FINT 
               ENDDO
               !
               DEALLOCATE( P ) 
            ENDDO
         ENDDO
         ! now interpolate y-direction
         ! setup spline coefficients in z-direction  
         DO I = -ESF%N1, ESF%N1
            ! map index to position in interpolated function
            I_ = 1 + I + ESF%N1 
#ifdef MPI
            !> block this loop 
            IND = LOCAL_INDEX_CYCLIC( NCPU, NODE_ME, I_ ) 
            IF ( IND < 1 ) CYCLE 
#endif
            DO K = -M3, M3
               K_ = K + M3 + 1
               !
               ! spline in y-direction:
               !        
               ALLOCATE( P( 2*ESF%N2 +1 , 5 ) )
               DO J = -ESF%N2, ESF%N2
                  ! index in interpolated function 
                  J_ = 1 + (J +ESF%N2)*SCALING(2)
                  ! setup { y, S(x,y,z) }  table
                  P( J + ESF%N2 + 1, 1 ) = I*    B(2,1) + J*    B(2,2) + K * (     B(2,3)/SCALING(3) )
                  P( J + ESF%N2 + 1, 2 ) = INTSTUFAK( K_ , J_ , I_ )
               ENDDO
               ! obtain natural-spline coefficients  (those with f"(boundary) = 0 )
               CALL  SPLCOF( P, 2*ESF%N2+1 , 2*ESF%N2+1, Y1P )
               ! now interpolate spline fit and store 
               DO J=-M2, M2 
                  J_ = J + M2 + 1 
                  Y = I*    B(2,1) + J*(     B(2,2)/SCALING(2) ) + K*(     B(2,3)/SCALING(3) )
                  CALL SPLVAL( Y , FINT, FITDER, P(1,1), 2*ESF%N2+1, 2*ESF%N2+1 )
                  ! store result in proper slot of interpolated structure factor 
                  INTSTUFAK( K_, J_, I_ ) = FINT 
               ENDDO
               !
               DEALLOCATE( P ) 
! put this in a separate subroutine 
            ENDDO
         ENDDO

         CALLMPI( M_sum_d( ESF%WGW%COMM, INTSTUFAK(1,1,1), SIZE(INTSTUFAK) ) )

         ! now interpolate x-direction
         DO K = -M3, M3
           ! map index to position in interpolated function
           K_ = K + M3 + 1
#ifdef MPI
           !> block this loop 
           IND = LOCAL_INDEX_CYCLIC( NCPU, NODE_ME, K_ ) 
           IF ( IND < 1 ) CYCLE 
#endif
           DO J = -M2, M2
              J_ = J + M2 + 1 
              ! spline in x-direction:
              ALLOCATE( P( 2*ESF%N1 +1 , 5 ) )
              DO I = -ESF%N1, ESF%N1
                 ! index in interpolated function 
                 I_ = 1 + I + ESF%N1 
                 ! x-coordinate
                 P( I + ESF%N1 + 1, 1 ) = I*    B(1,1) + J*(    B(1,2)/SCALING(2)) + K * (     B(1,3)/SCALING(3) )
                 P( I + ESF%N1 + 1, 2 ) =INTSTUFAK( K_ , J_ , I_ )
              ENDDO
              ! obtain natural-spline coefficients  (those with f"(boundary) = 0 )
              CALL  SPLCOF( P, 2*ESF%N1+1 , 2*ESF%N1+1, Y1P )
              ! now interpolate spline fit and store 
              DO I=-M1, M1 
                 X = I*(    B(1,1)/SCALING(1)) + J*(     B(1,2)/SCALING(2) ) + K*(     B(1,3)/SCALING(3) )
                 CALL SPLVAL( X , FINT, FITDER, P(1,1), 2*ESF%N1+1, 2*ESF%N1+1 )
                 ! FINT stores result of interpolated structure factor 
                 I_ = I + M1 + 1 
                 ! set kpoint coordinates for refined Coulomb potential
                 DK = I * ESF%B(:,1)/SCALING(1) &
                    + J * ESF%B(:,2)/SCALING(2) &
                    + K * ESF%B(:,3)/SCALING(3)

                 POTFAK = 0.5_q
                 IF ( (.NOT. LDMP1) .AND. LCUT ) THEN
                    ! setup refined Coulomb potential 
                    CALL REFINED_POTFAK(POTFAK, ESF%ENCUT, ESF%ENCUTSOFT, ESF%OMEGA, ESF%WGW%GRID%NPLWV, DOT_PRODUCT(DK,DK), .FALSE.)
                 ENDIF

#ifdef debug
IF ( IU0>=0 ) THEN
Y = I*(    B(2,1)/SCALING(1)) + J*(     B(2,2)/SCALING(2)) + K*(     B(2,3)/SCALING(3) )
Z = I*(    B(3,1)/SCALING(1)) + J*(     B(3,2)/SCALING(2)) + K*(     B(3,3)/SCALING(3) )
IF( I == 0 ) WRITE(100+ICOUNTER*10+L, '(3F20.10)' ) Y, Z, FINT
IF( J == 0 ) WRITE(200+ICOUNTER*10+L, '(3F20.10)' ) Z, X, FINT
IF( K == 0 ) WRITE(300+ICOUNTER*10+L, '(3F20.10)' ) X, Y, FINT
WRITE(400+ICOUNTER*10+L,'(3I5, 9F12.6)') I_, J_, K_, X,Y,Z, FINT
ENDIF
#endif
                 ! contract refined Coulomb potential with structure factor 
                 INTCORR(L)=INTCORR(L)+FINT*POTFAK/&
                   (SCALING(1)*SCALING(2)*SCALING(3))
              ENDDO
              DEALLOCATE( P ) 
            ENDDO
         ENDDO

         CALLMPI( M_sum_d( ESF%WGW%COMM, INTCORR(L), 1 ) )

!#ifdef verbose
         !> check for breaking condition 
         IF (IU0 >=0 ) THEN
           WRITE(IU0,'(A,I2, F12.6,4F20.10)') "INTCORR=", L, ESF%ENCUT, INTCORR(L), ABS(INTCORR(L) - EDIFF ), FSG_INT(L)
           WRITE(*,'(A,I2, F12.6,4F20.10)') "INTCORR=", L, ESF%ENCUT, INTCORR(L), ABS(INTCORR(L) - EDIFF ), FSG_INT(L)
         ENDIF
!#endif

         !> check for breaking condition 
         IF (  ABS(INTCORR(L) - EDIFF ) < ESF_CONV ) EXIT refinement

         !back up interpolation result for this refinement parameter
         EDIFF = INTCORR( L ) 

         CALL DEREGISTER_ALLOCATE(8._q*SIZE(INTSTUFAK,KIND=8), "esf" )
         DEALLOCATE(INTSTUFAK)

       ENDDO refinement

       !> check for breaking condition 
       IF (  L > ESF_NINTER ) THEN
          IF (IO%IU0 >=0 ) THEN
             WRITE(*  , '(A,E13.6)' ) "WARNING: convergence not reached for spline interpolation: ", ABS(INTCORR(L-2) - EDIFF ) 
             WRITE(IO%IU6, '(A,E13.6)' ) "WARNING: convergence not reached for spline interpolation: ", ABS(INTCORR(L-2) - EDIFF )
             WRITE(*  , '(A,I2, F12.6,3F20.10)') "INTCORR=", L-1, ESF%ENCUT, INTCORR(L-2), ABS(INTCORR(L-2) - EDIFF )
             WRITE(IO%IU6, '(A,I2, F12.6,3F20.10)') "INTCORR=", L-1, ESF%ENCUT, INTCORR(L-2), ABS(INTCORR(L-2) - EDIFF )
          ENDIF
          ! mark this result with an error code 
          ESF%ISTATUS = -1
          ESF%ESPLINED = INTCORR(L-1) + FSG_INT(MAX(1,L-1))

       !> successful exit 
       ELSE
          IF (IU0 >=0 ) THEN
             WRITE(IU0,'(A,I2, F12.6,3F20.10)') "INTCORR=", L, ESF%ENCUT, INTCORR(L), ABS(INTCORR(L) - EDIFF )
             WRITE(  *,'(A,I2, F12.6,3F20.10)') "INTCORR=", L, ESF%ENCUT, INTCORR(L), ABS(INTCORR(L) - EDIFF )
          ENDIF
          ! store result
          ESF%ESPLINED = INTCORR(L) + FSG_INT(MAX(1,L))
          ! mark this result with a successful code
          ESF%ISTATUS =  0
       ENDIF

       CALL DEREGISTER_ALLOCATE(8._q*SIZE(SXYZ,KIND=8), "esf")
       DEALLOCATE( SXYZ ) 
       PROFILING_STOP( "esf_spline_fit" )

       CONTAINS
   
       !----------------------------------------------------------------------
       !
       ! sets the Coulomb potential in carthesian coordinates
       ! 
       !----------------------------------------------------------------------
       SUBROUTINE REFINED_POTFAK(POTFAK, ENCUT, ENCUTSOFT, OMEGA, NPLWV, GSQU, LSR)
          USE constant
          USE fock
          USE c2f_interface, ONLY : ERRF
          USE fermi_energy, ONLY : DELSTP
          IMPLICIT NONE
          REAL(q), INTENT(INOUT)  :: POTFAK    !< refined Coulomb potential
          REAL(q), INTENT(IN)     :: ENCUT     !< energy cutoff
          REAL(q), INTENT(IN)     :: ENCUTSOFT !< soft cutoff
          REAL(q), INTENT(IN)     :: OMEGA !< reciprocal cell volume
          INTEGER, INTENT(IN)     :: NPLWV !< number of splane waves
          REAL(q), INTENT(IN)     :: GSQU  !< 1/|q+G|^2
          LOGICAL, INTENT(IN)     :: LSR   !< short range ?
       !  local
          REAL(q) :: E
          REAL(q) :: DFUN, SFUN
          REAL(q) :: OMEGBK,QC
          REAL(q) :: GAMMASCALE
          REAL(q) :: Q1, Q2

          IF (ESF%ENCUTSOFT>=0) THEN
             Q1=SQRT(ESF%ENCUTSOFT/HSQDTM)
             Q2=SQRT(ESF%ENCUT/HSQDTM)
          ENDIF
         
          IF (ABS(GSQU)<G2ZERO) THEN
             ! head and wing
             POTFAK=EDEPS/OMEGA
             ! if HFRCUT is set, a fixed spherical cutoff is used
             ! resulting in a finite values at G=0
             ! since the head stores the q^2, the final contribution is zero
             ! from the head
             IF (HFRCUT/=0) THEN
                POTFAK=0
#ifdef HF 
             ! if LRHFCALC.AND.LRSCOR the Coulomb kernel is replaced by v_{LR}=v*exp(-q^2/(4*mu^2))
             ! for q=0  v_{LR}=1, this introduces an error for small mu because of the finite q-grid
             ! In order to correct this error, v_{LR}(q=0) is set to average of v_{LR} in a sphere  
             ! with volume corresponding to BZ volume/ # q-points       
             ELSE IF (LRHFCALC.AND.LRSCOR) THEN
                CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
                QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
                GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
               &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
               &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
                POTFAK=POTFAK*GAMMASCALE
             ! if LSR the Coulomb kernel is replaced by v_{SR}=v*[1-exp(-q^2/(4*mu^2))]
             ! for q=0  v_{SR}=0, this introduces an error for small mu because of the finite q-grid
             ! In order to correct this error, v_{SR}(q=0) is set to average of v_{SR} in a sphere
             ! with volume corresponding to BZ volume/ # q-points
             ELSE IF (LSR) THEN
                CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
                QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
                GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
               &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
               &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
                POTFAK=POTFAK*(1-GAMMASCALE)
#endif
             ENDIF
     
             ! switch off standard convergence correction
             IF(MCALPHA/=0) THEN
                POTFAK=0
             ENDIF
     
          ELSE
             ! the factor 1/(2 pi)^2 is required to obtain proper reciprocal
             ! lattice vector lenght
             POTFAK=EDEPS/OMEGA/(GSQU*TPI**2)
             IF (HFRCUT/=0) THEN
                ! spherical cutoff on Coloumb kernel
                ! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
                CALL DELSTP(3,SQRT(GSQU)*TPI*HFRCUT/10,DFUN,SFUN)
                POTFAK=POTFAK*(SFUN-0.5)*2
             ELSE IF (LSR) THEN
                ! use the short range part of the Coulomb kernel
                POTFAK=POTFAK*(1._q-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
             ELSE IF (LRHFCALC.AND.LRSCOR) THEN
                ! use the long range part of the Coulomb kernel
                POTFAK=POTFAK*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
             ENDIF
          ENDIF
     
          ! smooth cutoff function between  ENCUTSOFT and ENCUT
          E=HSQDTM*(GSQU*TPI**2)
          IF (ENCUT>=0 .AND. E>ENCUT) THEN
             POTFAK=0
          ELSE
             IF (ENCUTSOFT>=0 .AND. E>ENCUTSOFT) THEN
                ! cosine window
                IF ( .NOT. LSCK) THEN
                   POTFAK=POTFAK*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
                ! squeezed Coulomb kernel
                ELSE 
                   POTFAK=POTFAK*SQUEEZED_COULOMB_KERNEL(SQRT(GSQU)*TPI, Q1, Q2 )
                ENDIF
             ENDIF
          ENDIF
          ! POTFAK in fock has an additional factor. Since the corrections were originally devised
          ! for the fock routine this ensures consistency
          IF(MCALPHA/=0) THEN
             POTFAK=POTFAK*(1.0_q/NPLWV)
          ENDIF

       END SUBROUTINE REFINED_POTFAK

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> S_G(q) is used to set up STUFAK(X,Y,Z)
   !> this is done, by looping over all k-points in the full BZ and using 
   !> S_G(q) for q in IBZ to rotate the contribution to k in full BZ
   !---------------------------------------------------------------------------
       SUBROUTINE SETUP_SXYZ( ESF, SXYZ, LCUTOFF, IO )
          USE full_kpoints, ONLY : KPOINTS_FULL, LIDENTICAL_KPOINT 
          TYPE( elstufac_handle ) :: ESF
          REAL(q), INTENT(INOUT)  :: SXYZ(:,:,:)
          LOGICAL                 :: LCUTOFF ! smooth cutoff is applied, so reordering of indices is required
          TYPE(in_struct)         :: IO
          !> local 
          INTEGER :: NK, NQ
          INTEGER :: NI, NI_
          REAL(q) :: DK(3), QG(3), G(3)
          INTEGER :: NGVECTOR
          INTEGER :: I1,I2,I3
          INTEGER :: I1_,I2_,I3_
          GDEF, ALLOCATABLE :: STUFAK_DIAG(:)
          TYPE( wavedes1 )       :: WGWK
          TYPE( wavedes1 )       :: WGWQ
          INTEGER                :: NP_NEW
          INTEGER, POINTER       :: OLD_INDEX(:) =>NULL()
          REAL(q), POINTER       :: DATAKE(:) =>NULL()

          PROFILING_START( "setup_sxyz" )
!TODO: k-point lattice of FOCK_ACC differs than the one of ACFDT 

          ! loop over irreducible brillouin zone points
          qpoints: DO NQ=1,ESF%NUMBER_OF_NQ
             !WRITE(*,*) 'doing ', NQ, ESF%WGW%NKPTS , ESF%NG(NQ)
             IF ( ESF%NG( NQ ) <0 )THEN
                CALL VTUTOR%BUG("NG not set for interpolation",__FILE__,__LINE__ )
             ENDIF
             CALL SETWDES(ESF%WGW, WGWQ, ESF%NQ(NQ) )

             ! reshuffling of indices is requred, setup proper potential with proper cutoff 
             IF ( LCUTOFF ) THEN
                ESF%OLD_INDEX = 0 
                CALL SETUP_COULOMB_DIAGONAL( ESF%LATT_CUR, ESF%ENCUT, ESF%ENCUTSOFT, WGWQ, NP_NEW, .FALSE., OLD_INDEX, DATAKE )
             ELSE
                DEALLOCATE( OLD_INDEX ) 
                NULLIFY( OLD_INDEX )
                
                ALLOCATE( OLD_INDEX( ESF%NG(NQ) ) )
                DO NI = 1, ESF%NG(NQ)
                   OLD_INDEX( NI ) = NI 
                ENDDO
             ENDIF
           
             kpoints: DO NK=1,KPOINTS_FULL_ORIG%NKPTS
                IF( SKIP_THIS_KPOINT_IN_FOCK( KPOINTS_FULL_ORIG%VKPT(:,NK) ) ) CYCLE 
                CALL SETWDES(ESF%WGW, WGWK, NK )
                
                IF( ESF%ISYM>=0 ) THEN
                   ! loop over symmetry - equivalent points
                   IF ( ESF%KPOINTS_TRANS%NK_OLD( NK ) /= ESF%NQ(NQ) ) CYCLE kpoints
                ELSE
                   ! loop over equivalent points
                   IF ( .NOT. LIDENTICAL_KPOINT( WGWK%VKPT, WGWQ%VKPT ) ) CYCLE kpoints
                ENDIF 
!IF( IO%IU0>=0)&
!WRITE(*,'(3I4,3F12.6, I4)') NQ, ESF%NUMBER_OF_NQ,NK, KPOINTS_FULL_ORIG%VKPT(:,NQ), ESF%KPOINTS_TRANS%NK_OLD( NK )
             ! loop over full brillouin zone points

                ! number of plane waves used for this q-point 
                NGVECTOR = WGWK%NGVECTOR
#ifdef gammareal
                NGVECTOR = NGVECTOR * 2 
#endif
                !set kpoint coordinates
                !DK(:)=WGWK%VKPT(1)*LATT_CUR%B(:,1)+WGWK%VKPT(2)*LATT_CUR%B(:,2)+WGWK%VKPT(3)*LATT_CUR%B(:,3)
                DK(:)=WGWK%VKPT(1)*ESF%B(:,1)*ESF%NKPX&
                     +WGWK%VKPT(2)*ESF%B(:,2)*ESF%NKPY&
                     +WGWK%VKPT(3)*ESF%B(:,3)*ESF%NKPZ

                !store diagonal
                ALLOCATE( STUFAK_DIAG( NGVECTOR ) )
                STUFAK_DIAG = 0 
                DO NI= 1, ESF%NG(NQ) 
                   STUFAK_DIAG( OLD_INDEX( NI ) ) = ESF%S( NI, NQ )
                ENDDO
                ! rotate S(q+G) -> S(k+G) 
                !IF( ESF%ISYM>0  .AND.  .NOT. LIDENTICAL_KPOINT( KPOINTS_FULL%VKPT(:,NK), WGWQ%VKPT(:) ) ) THEN
                IF( ESF%ISYM>=0  .AND.  .NOT.  LIDENTICAL_KPOINT( WGWK%VKPT, WGWQ%VKPT )  ) THEN
                    ! now rotate structure factor to NK 
                    CALL ROTATE_WPOT_DIAG_INPLACE(STUFAK_DIAG, ESF%WGW%NGVECTOR(NK), &
                            ESF%KPOINTS_TRANS%CPHASE(:,NK), ESF%KPOINTS_TRANS%NINDPW(:,NK),  &
                            ESF%KPOINTS_TRANS%LINV(NK), ESF%KPOINTS_TRANS%LSHIFT(NK) )
                ENDIF 

!dump esf here 
#ifdef debug
IF ( IO%IU0 >=0 ) THEN
WRITE( 100+NQ*10+NK, '( A,2I4,F20.10,3F12.6)') "#", NQ, NK, ESF%ENCUT, WGWK%VKPT
DO I1 = 1, ESF%NG(NQ)
WRITE( 100+NQ*10+NK, '( 2I4,2F20.10)') I1, OLD_INDEX(I1), REAL(STUFAK_DIAG( I1 )), ESF%S( I1, NQ )
ENDDO
ENDIF
STOP
#endif 

                gpoints: DO NI=1,NGVECTOR
                   ! skip zero contributions 
                   IF ( ABS( STUFAK_DIAG( NI ) ) < 1.E-6_q ) CYCLE gpoints

                   NI_ = NI 
                   IF (WGWK%LGAMMA) NI_=(NI-1)/2+1
                   G(:) = WGWK%IGX( NI_ )*ESF%B(:,1)*ESF%NKPX&  
                        + WGWK%IGY( NI_ )*ESF%B(:,2)*ESF%NKPY&
                        + WGWK%IGZ( NI_ )*ESF%B(:,3)*ESF%NKPZ  

                   QG=DK+G

                   !obtain cartesian coordinates
                   I1 = NINT( DOT_PRODUCT( ESF%A( 1,:), QG(:) ) )
                   I2 = NINT( DOT_PRODUCT( ESF%A( 2,:), QG(:) ) )
                   I3 = NINT( DOT_PRODUCT( ESF%A( 3,:), QG(:) ) )

                   IF (I1< -ESF%N1 .OR. I1 > ESF%N1 .OR. &
                       I2< -ESF%N2 .OR. I2 > ESF%N2 .OR. &
                       I3< -ESF%N3 .OR. I3 > ESF%N3 ) THEN
                       CALL vtutor%bug( "q,G= "//str(WGWK%VKPT)//", "//str(G)//" out of bounds "//str( [I1,I2,I3,ESF%N1,ESF%N2,ESF%N3] )//" S= "//str( STUFAK_DIAG(NI) ), __FILE__, __LINE__)
                   ENDIF 

                   ! shift to 
                   I1_ = I1 + ESF%N1 + 1 
                   I2_ = I2 + ESF%N2 + 1 
                   I3_ = I3 + ESF%N3 + 1 
                   ! store structure factor and coulomb potential 
                   SXYZ(I1_,I2_,I3_)=SXYZ(I1_,I2_,I3_)+REAL(STUFAK_DIAG( NI ),q)/KPOINTS_FULL_ORIG%NKPTS&
                     *NKREDX* NKREDY* NKREDZ
                ENDDO gpoints 
                DEALLOCATE( STUFAK_DIAG )
             ENDDO kpoints
          ENDDO qpoints
          PROFILING_STOP( "setup_sxyz" )
       END SUBROUTINE SETUP_SXYZ

    END SUBROUTINE ESF_SPLINE_FIT

    !> this routine checks integrity from FFT indices to carthesian coordinates 
    SUBROUTINE GRID_TRAFO_CHECK( ESF ) 
       TYPE( elstufac_handle ) :: ESF 
       !> local variables
       INTEGER :: I, J, K, I_, J_, K_
       INTEGER :: I1, I2, I3
       REAL(q) :: X,Y,Z
     
       DO K = 1, ESF%N3* 2 + 1 
          K_ = K -1 - ESF%N3
          DO J = 1, ESF%N2* 2 + 1 
             J_ = J -1 - ESF%N2
             DO I = 1, ESF%N1* 2 + 1 
                I_ = I -1 - ESF%N1
                X = I_*ESF%B(1,1) + J_* ESF%B(1,2) + K_*ESF%B(1,3)
                Y = I_*ESF%B(2,1) + J_* ESF%B(2,2) + K_*ESF%B(2,3)
                Z = I_*ESF%B(3,1) + J_* ESF%B(3,2) + K_*ESF%B(3,3)
                ! use X_, Y_, Z_, to map back to I_,J_,K_
                I1 = NINT( DOT_PRODUCT( ESF%A(1,:) , [X,Y,Z] ) )
                I2 = NINT( DOT_PRODUCT( ESF%A(2,:) , [X,Y,Z] ) )
                I3 = NINT( DOT_PRODUCT( ESF%A(3,:) , [X,Y,Z] ) )

                IF(  I_ /= I1 .OR. J_ /= I2 .OR. K_ /= I3 ) THEN
                   CALL VTUTOR%BUG( " grid transform inconsistent "//str( [I_, J_, K_, I1, I2, I3 ] ), __FILE__ , __LINE__ )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    END SUBROUTINE GRID_TRAFO_CHECK 

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> rotates a vector S_G(q) for q in the irreducible BZ to a k-point in the 
   !> full BZ
   !> @param[inout] S      vector S_G(q) 
   !> @param[in]    NP     number of plane wave components G 
   !> @param[in]    CHPASE vector of phases for each G component for point k 
   !> @param[in]    NINDPW index array for (q+G) -> (k+G) 
   !> @param[in]    LINV   k-point is connect with q-point via mirror symmery ?
   !> @param[in]    LSHIFT k-point requires multiplication with phase factor CPHASE 
   !> 
   !---------------------------------------------------------------------------
   SUBROUTINE ROTATE_WPOT_DIAG_INPLACE( S, NP, CPHASE, NINDPW, LINV, LSHIFT)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT NONE

      GDEF,INTENT(INOUT)       :: S(:)
      INTEGER,   INTENT(IN)    :: NP
      COMPLEX(q),INTENT(IN)    :: CPHASE(:)
      INTEGER, INTENT(IN)      :: NINDPW(:)
      LOGICAL, INTENT(IN)      :: LINV
      LOGICAL, INTENT(IN)      :: LSHIFT
      ! local
      INTEGER N
      GDEF,ALLOCATABLE         :: S_ROT(:)

      ALLOCATE(S_ROT(NP))
      !first build the complex conjugated if required
      S_ROT = 0 
      IF (LINV) THEN
         DO N=1,NP
            S_ROT( N ) = GCONJG( S( N ) )
         ENDDO
      ELSE
         DO N=1,NP
            S_ROT( N ) =  S( N )
         ENDDO
      ENDIF

      !second apply phase factors
      IF (LSHIFT) THEN
         DO N=1,NP
            S(N)=S_ROT(N)*CONJG(CPHASE( N ))
         ENDDO
         DO N=1,NP
            S_ROT(N)=CPHASE( N )* S(N)
         ENDDO
      ENDIF
      !third re-index rows and columns
      DO N=1,NP
         S(NINDPW(N))=S_ROT(N)
      ENDDO
      DEALLOCATE(S_ROT)

   END SUBROUTINE ROTATE_WPOT_DIAG_INPLACE

END MODULE esf
