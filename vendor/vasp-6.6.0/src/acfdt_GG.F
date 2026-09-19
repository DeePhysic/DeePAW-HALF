#include "symbol.inc"

MODULE acfdt_gg
  USE chi_base
  USE GG_base 
  USE fock 
  USE acfdt
  USE esf
  USE minimax_struct, ONLY : imag_grid_handle
  USE fermi_energy, ONLY : DELSTP
  IMPLICIT NONE

  !< number of linear regression points (or individual cutoffs of chi* V )
  INTEGER, PRIVATE :: NE=8
  !< number spacing between regression points 
  REAL(q), PARAMETER, PRIVATE :: STEP_BETWEEN_FREQUENCIES=1.05
! alternative spanning the same frequency range but with only 4 points
!  INTEGER, PARAMETER, PRIVATE :: NE=4
!  REAL(q) :: STEP_BETWEEN_FREQUENCIES=1.120576983388_q

!> Fermi-wave vector for HEG kernel
!> set some default to avoid crashing 
!
  REAL(q), PRIVATE :: KFERMI=2.0_q

#ifdef scaLAPACK

CONTAINS

!********************************* XI_ACFDT_SETUP   *************************************
!
!> allocates correaltion handle used in this module, i.e. 
!> sets up the energy cutoffs at which the correlation energy is 
!> calculated 
! 
!****************************************************************************************

  SUBROUTINE XI_ACFDT_SETUP_GG( COR, ENCUTGW, ENCUTGWSOFT)
    USE constant
    TYPE (correlation), POINTER :: COR
    REAL(q) :: ENCUTGW, ENCUTGWSOFT
    INTEGER :: I

    ! use single shot method with SCK for RPA 
    IF ( LSCK ) THEN
       NE = 1 
    ELSE
       NE = NE_REG
    ENDIF 

    ALLOCATE(COR)
    COR%NE=NE
    ALLOCATE(COR%ENCUTGW(NE))
    ALLOCATE(COR%ENCUTGWSOFT(NE))
    ALLOCATE(COR%CORRELATION(NE))
    ALLOCATE(COR%CORRSOSEX(NE))
    ALLOCATE(COR%CORRELATION_K(NE))
    ALLOCATE(COR%CORRSOSEX_K(NE))
    ALLOCATE(COR%CORRMP2DIR(NE))
    ALLOCATE(COR%CORRMP2EX(NE))
    ALLOCATE(COR%CORRMP2DIR_K(NE))
    ALLOCATE(COR%CORRMP2EX_K(NE))
    ALLOCATE(COR%CORRELATION_LAMBDA(NE,0:NLAMBDA))

    COR%ENCUTGW(1)    =ENCUTGW
    COR%ENCUTGWSOFT(1)=ENCUTGWSOFT
    DO I=2,NE
       COR%ENCUTGW(I)    =COR%ENCUTGW(I-1)/STEP_BETWEEN_FREQUENCIES
       COR%ENCUTGWSOFT(I)=COR%ENCUTGWSOFT(I-1)/STEP_BETWEEN_FREQUENCIES
       IF (COR%ENCUTGWSOFT(I)<=0) COR%ENCUTGWSOFT(I)=-1
    ENDDO

    COR%CORRELATION_LAMBDA=0
    COR%CORRELATION=0
    COR%CORRSOSEX=0
    COR%CORRELATION_K=0
    COR%CORRSOSEX_K=0
    COR%CORRMP2DIR=0
    COR%CORRMP2EX=0
    COR%CORRMP2DIR_K=0
    COR%CORRMP2EX_K=0
    
  END SUBROUTINE XI_ACFDT_SETUP_GG


!********************************* XI_ACFDT_SETUP   *************************************
!
!> deallocates correlation handle 
! 
!****************************************************************************************
  SUBROUTINE XI_ACFDT_DEALLOCATE_GG( COR )
    TYPE (correlation), POINTER :: COR
    ! local 
    INTEGER :: I 
    DEALLOCATE(COR%ENCUTGW)
    DEALLOCATE(COR%CORRELATION)
    DEALLOCATE(COR%CORRSOSEX)
    DEALLOCATE(COR%CORRELATION_K)
    DEALLOCATE(COR%CORRSOSEX_K)
    DEALLOCATE(COR%CORRMP2DIR)
    DEALLOCATE(COR%CORRMP2EX)
    DEALLOCATE(COR%CORRMP2DIR_K)
    DEALLOCATE(COR%CORRMP2EX_K)
    DEALLOCATE(COR%CORRELATION_LAMBDA)

    ! in case electronic structure factor has been collected
    ! deallocate electronic structure factor
    IF ( LESF_SPLINES ) THEN
       DO I = 1, COR%NE 
          CALL DEALLOCATE_ELSTUFAC_HANDLE( COR%ESF(I) )
       ENDDO
       DEALLOCATE( COR%ESF )
       NULLIFY( COR%ESF )
    ENDIF
    DEALLOCATE(COR)
    NULLIFY(COR)

  END SUBROUTINE XI_ACFDT_DEALLOCATE_GG

!********************************* CALCULATE_RPA_CORRELATION_ENERGY *********************
!
!> determine the RPA correlation energy
!> this subroutine integrates the entire calculation in one
!> routine to conserve memory and speed things up
! 
!****************************************************************************************
SUBROUTINE CALCULATE_RPA_CORRELATION_ENERGY( WDES, CHI, GDES, WGWQ, &
    IMAG_GRIDS, LATT_CUR, COR, NQ , IDIR_MAX, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE minimax, ONLY : LOCAL_INDEX_TO_GLOBAL 
    IMPLICIT NONE
    TYPE (wavedes)          :: WDES                 !< wave function descriptor
    TYPE (responsefunction) :: CHI                  !< response function
    TYPE (wavedes1)         :: WGWQ                 !< response function descriptor for current q-point
    TYPE (imag_grid_handle) :: IMAG_GRIDS           !< imaginary grids 
    TYPE (greensfdes)       :: GDES                 !> greensfunction descriptor
    TYPE (latt)             :: LATT_CUR             !< lattice structure
    TYPE (correlation)      :: COR                  !< correlation energy
    INTEGER                 :: NQ                   !< current q-point
    INTEGER                 :: IDIR_MAX             !< maximum number of directions
    TYPE (in_struct)        :: IO                 
  ! local
    GDEF, POINTER,CONTIGUOUS:: CHI_WORK(:,:) => NULL()

    TYPE (communic)         :: COMM_INOMEGA           !communicator in tau
    INTEGER                 :: NOMEGA               !number of frequency points

    INTEGER                 :: NROWS_NEW, NCOLS_NEW !new size of redistributed chi 
    INTEGER                 :: NBLOCK                 !block size
    INTEGER ::  IDIR, NP_NEW, I,IP, ILAMBDA, NK, NCURR
    INTEGER :: I_, J
    COMPLEX(q) :: SUM, SUMMP2, SUMSOSEX, SUMMP2EX, SUM_sr, SUMMP2_sr
    REAL(q) :: LAMBDA
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group

    ! total number of points
    NOMEGA = IMAG_GRIDS%NOMEGA

    COMM_INOMEGA = IMAG_GRIDS%B%COMM_IN_GROUP

    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("internal error in CALCULATE_RPA_CORRELATION_ENERGY: \n or the Gamma point " &
          // "version CHI%LREALSTORE must be set", __FILE__, __LINE__)
    ENDIF

    PROFILING_START('calculate_rpa_correlation_energy')

    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = IMAG_GRIDS%B%NPOINTS_IN_GROUP
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    ! Merzuk's version changed the sign of RESPONSEFUN here
    ! I moved this over to the place where the RESPONSEFUN is first used

#ifdef ACC_OFFLOAD
    CALL OFFLOAD_ON_PUSH(.TRUE.)
!$ACC ENTER DATA COPYIN(CHI) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(CHI%HEAD) IF(OFFLOAD_ON .AND. CHI%LGAMMA) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(CHI%RESPONSER) IF(OFFLOAD_ON .AND. CHI%LREALSTORE) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(CHI%CWINGR,CHI%WINGR) IF(OFFLOAD_ON .AND. CHI%LREALSTORE .AND. CHI%LGAMMA) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(CHI%RESPONSEFUN) IF(OFFLOAD_ON .AND. .NOT.CHI%LREALSTORE) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(CHI%CWING,CHI%WING) IF(OFFLOAD_ON .AND. .NOT.CHI%LREALSTORE .AND. CHI%LGAMMA) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(COMM_INOMEGA) __IF_ASYNC__
#endif

    !--------------------------------------------------------------------------------------------
    lamb: DO ILAMBDA=0,NLAMBDA                       !optionally loop over lambda parameter
    LAMBDA=1.0_q*(ILAMBDA+1)/MAX(1,NLAMBDA)          !set current lambda value
       !-----------------------------------------------------------------------------------------
       cutoff: DO I=1, NE                          !loop over different energy cutoffs
       !-----------------------------------------------------------------------------------------
          !initialize correlation energy for the current energy cutoff # I:
          COR%CORRELATION_K(I)=0
          COR%CORRMP2DIR_K(I)=0
          omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
          !--------------------------------------------------------------------------------------
             ! set local omega_frequency loop             
             CALL SETUP_LOOPDES_INDICES( NCURR, IMAG_GRIDS%B )
             NOMEGA_GLOBAL = IMAG_GRIDS%B%NPOINT_GLOBAL
             !-----------------------------------------------------------------------------------
             ! directions only need to be considered for the head 
             IF ( CHI%LGAMMA ) THEN
             directions: DO IDIR = 1, IDIR_MAX          !loop over the tree directions  
             !-----------------------------------------------------------------------------------
                !for the long-wave limit we copy head and wings to RESPONSEFUN
                CALL BODY_FROM_WING_GG( CHI, IDIR, IMAG_GRIDS%B%NPOINTSC, GDES, COMM_INOMEGA)

                !CHI is distributed in groups to all processors
                !each group contains CHI at one frequency point. 
                !multiplication with the Coulomb potential needs to be done carefully,
                !CHI will be used for GW calculations, therefore
                !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
                !multiply with v^1/2 from left and right
                !furthermore set the correct energy cutoff
                CALL XI_LOCAL_FIELD_ACFDT_GG( CHI_WORK, CHI, GDES, &
                   COMM_INOMEGA, LATT_CUR, COR, I, IMAG_GRIDS%B%NPOINTSC, WGWQ, NP_NEW, .FALSE. )

                !use scaLAPACK to compute eigenvalues of CHI_WORK
                CALL COL2BC_AND_RPA(CHI_WORK, GDES, IMAG_GRIDS%B, NP_NEW, COR, NQ,&
                   IMAG_GRIDS%BOS_RE_WEIGHT(NOMEGA_GLOBAL)/IDIR_MAX, I, IO )

             !-----------------------------------------------------------------------------------
             ENDDO directions  
             ELSE 
             !-----------------------------------------------------------------------------------
                !CHI is distributed in groups to all processors
                !each group contains CHI at one frequency point. 
                !multiplication with the Coulomb potential needs to be done carefully,
                !CHI will be used for GW calculations, therefore
                !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
                !multiply with v^1/2 from left and right
                !furthermore set the correct energy cutoff
                CALL XI_LOCAL_FIELD_ACFDT_GG( CHI_WORK, CHI, GDES, &
                   COMM_INOMEGA, LATT_CUR, COR, I, IMAG_GRIDS%B%NPOINTSC, WGWQ, NP_NEW, .FALSE. )

                !use scaLAPACK to compute eigenvalues of CHI_WORK
                CALL COL2BC_AND_RPA(CHI_WORK, GDES, IMAG_GRIDS%B, NP_NEW, COR, NQ,&
                   IMAG_GRIDS%BOS_RE_WEIGHT(NOMEGA_GLOBAL), I, IO )
             !-----------------------------------------------------------------------------------
             ENDIF
          !--------------------------------------------------------------------------------------
          ENDDO omega_sim
       !-----------------------------------------------------------------------------------------
       ENDDO cutoff
       !-----------------------------------------------------------------------------------------
   
       !finally communicate the result between nodes (note that COR%CORR... is complex)
       CALLMPI( M_sum_d(IMAG_GRIDS%B%COMM_BETWEEN_GROUPS, COR%CORRELATION_K, 2*NE ))  !RPA
       CALLMPI( M_sum_d(IMAG_GRIDS%B%COMM_BETWEEN_GROUPS, COR%CORRMP2DIR_K, 2*NE ))   !MP2

       !scale appropriately  
       COR%CORRELATION_K=COR%CORRELATION_K/2*KPOINTS_ORIG%WTKPT(NQ) 
       COR%CORRMP2DIR_K =COR%CORRMP2DIR_K/2*KPOINTS_ORIG%WTKPT(NQ) 
       !save the result ot CORELATION and MP2DIR      
       COR%CORRELATION=COR%CORRELATION  +COR%CORRELATION_K
       COR%CORRMP2DIR =COR%CORRMP2DIR   +COR%CORRMP2DIR_K
                 
    !--------------------------------------------------------------------------------------------
    ENDDO lamb
    !---------------------------------------------------------------------------------------------

    ! accumulate electronic structure factor 
    IF( LESF_SPLINES ) THEN
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       IF ( NQ == COR%ESF(1)%NUMBER_OF_NQ ) THEN
#ifdef MPI
          DO I = 1,NE 
             ! require a collective sum for the structure factor 
             CALLMPI( M_sum_d(IMAG_GRIDS%B%COMM_BETWEEN_GROUPS, COR%ESF(I)%S(1,1), SIZE( COR%ESF(I)%S) ))
          ENDDO
#endif   
          CALL EVALUATE_TR( COR , IMAG_GRIDS%B%COMM_BETWEEN_GROUPS ) 
       ENDIF
    ENDIF


#ifdef ACC_OFFLOAD
!note(sm): CHI%REPSONSEFUN and CHI%RESPONSER are modified on the GPU in the CHI%LGAMMA case and might need
!          to be copied out instead. However, it looks like the data is not used anymore as the results are ok
!$ACC EXIT DATA DELETE(CHI%RESPONSER) IF(OFFLOAD_ON .AND. CHI%LREALSTORE) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(CHI%CWINGR,CHI%WINGR) IF(OFFLOAD_ON .AND. CHI%LREALSTORE .AND. CHI%LGAMMA) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(CHI%RESPONSEFUN) IF(OFFLOAD_ON .AND. .NOT.CHI%LREALSTORE) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(CHI%CWING,CHI%WING) IF(OFFLOAD_ON .AND. .NOT.CHI%LREALSTORE .AND. CHI%LGAMMA) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(CHI%HEAD) IF(OFFLOAD_ON .AND. CHI%LGAMMA) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(CHI) __IF_ASYNC__
!$ACC EXIT DATA DELETE(COMM_INOMEGA) __IF_ASYNC__
    CALL OFFLOAD_ON_POP
#endif

    PROFILING_STOP('calculate_rpa_correlation_energy')

    CONTAINS

    SUBROUTINE EVALUATE_TR( COR , COMM_BETWEENOMEGA) 
       TYPE( correlation) :: COR
       TYPE( communic ) :: COMM_BETWEENOMEGA
       ! local
       REAL(q) :: CORR
       INTEGER :: NQ, NI 
        
#ifdef debug          
       DO I = 1,NE 
          DO NQ = 1, COR%ESF(I)%NUMBER_OF_NQ 
IF( IO%IU0 >= 0) &
             WRITE(200+I, '(2I5, 3F20.10)') NQ
             DO NI = 1, COR%ESF(I)%NG( NQ ) 
IF( IO%IU0 >= 0 ) &
                 WRITE(200+I, '(I5, F20.10)') NI, COR%ESF(I)%S(NI, NQ)
             ENDDO
          ENDDO
       ENDDO
#endif 
    END SUBROUTINE EVALUATE_TR
END SUBROUTINE CALCULATE_RPA_CORRELATION_ENERGY

!******************************* XI_LOCAL_FIELD_ACFDT_GG ********************************
!
!> this subroutine truncates the response function by multiplying
!> with a truncated  Coulomb kernel smoothly going from 1/G^2 to zero
!> between ENCUTGWSOFT and ENCUTGW
!
!****************************************************************************************

  SUBROUTINE XI_LOCAL_FIELD_ACFDT_GG( CHI_WORK, CHI, GDES, &
    COMM_INOMEGA, LATT_CUR, COR, I, NOMEGA_LOCAL, WGWQ , MAXINDEX, LSR)
    USE constant
    USE fock
    USE mpimy
    USE c2f_interface, ONLY : ERRF
    IMPLICIT NONE
    GDEF,POINTER,CONTIGUOUS :: CHI_WORK(:,:)   !< current response function with proper cutoff
    TYPE (responsefunction) :: CHI             !< unchanged response function
    TYPE (greensfdes)       :: GDES            !< greensfunction descriptor
    TYPE (latt) LATT_CUR                       !< lattice structure 
    TYPE(communic)          :: COMM_INOMEGA    !< communicator inside frequency group
    TYPE(correlation)       :: COR
    INTEGER, INTENT(IN)     :: I 
    INTEGER, INTENT(IN)     :: NOMEGA_LOCAL  !< local freuquency point
    TYPE (wavedes1)         :: WGWQ          !< response function descriptor for current q-point
    INTEGER                 :: MAXINDEX      !< number of grid points using truncated Coulomb kernel
    LOGICAL, INTENT(IN)     :: LSR           !< range separated RPA is calculated
  ! local
    REAL(q) :: ENCUT                 !< energy cutoff for response function
    REAL(q) :: ENCUTSOFT             !< lower cutoff for response function 
    INTEGER    NI, NP, NI_
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE, POTFAK, E
    INTEGER :: NI_LOC, NOLD_LOC, NODE_RECV, IERROR
    REAL(q), ALLOCATABLE :: DATAKE(:)
    INTEGER, ALLOCATABLE :: OLD_INDEX(:)
    INTEGER :: ISEND(4) , IRECV(4)
    GDEF, ALLOCATABLE    :: CHI_RECV(:)
    GDEF, ALLOCATABLE    :: CHI_SEND(:)
    REAL(q) :: DFUN, SFUN
    ! neu
    REAL(q) :: OMEGBK,QC
    REAL(q) :: GAMMASCALE

    ! multipole correction
    INTEGER    NJ,NK,NO
    COMPLEX(q) :: TSUM(10)
    COMPLEX(q) :: POTCORRECTION(10)
    !for communication
    INTEGER    :: NI_GLOB, NSEND
    INTEGER    :: NODE_SEND
    INTEGER    :: NROWS                 !< number of rows 
    INTEGER    :: NCOLS                 !< number of columns 

    REAL(q) :: Q1, Q2

    PROFILING_START('xi_local_field_acfdt_gg')

    NROWS = GDES%RES_NRPLWV_ROW_DATA_POINTS
    NCOLS = GDES%RES_NRPLWV_COL_DATA_POINTS

    NULLIFY( CHI_WORK )
    ALLOCATE(CHI_WORK(NROWS,NCOLS)) 
!$ACC ENTER DATA CREATE(CHI_WORK) __IF_ASYNC__
!$ACC KERNELS PRESENT(CHI_WORK) __IF_ASYNC__
    CHI_WORK = 0
!$ACC END KERNELS


    ENCUT = COR%ENCUTGW(I)
    ENCUTSOFT = COR%ENCUTGWSOFT(I)

    IF (ENCUTSOFT>=0) THEN
       Q1=SQRT(ENCUTSOFT/HSQDTM)
       Q2=SQRT(ENCUT/HSQDTM)
    ENDIF
    ! e^2/ volume
    ! mK scales by EDEPS/LATT_CUR%OMEGA here
    SCALE=EDEPS/LATT_CUR%OMEGA
!    SCALE=1.0_q
!=======================================================================
! first set up the truncated Coulomb kernel
! smoothly going from 1/G^2 to zero between ENCUTGWSOFT and ENCUTGW
!=======================================================================
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

    ALLOCATE( DATAKE(NROWS), OLD_INDEX(NROWS))
    DO NI = 1, NROWS
     OLD_INDEX(NI) = NI
    ENDDO

    PROFILING_START('loop1|xi_local_field_acfdt_gg')
    DATAKEMAX=0.0
    MAXINDEX=0
    DO NI=1,NP
       NI_=NI
       IF (WGWQ%LGAMMA) NI_=(NI-1)/2+1
       
       GX=(WGWQ%IGX(NI_)*LATT_CUR%B(1,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WGWQ%IGX(NI_)*LATT_CUR%B(2,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WGWQ%IGX(NI_)*LATT_CUR%B(3,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
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
          POTFAK=POTFAK*(1.0_Q/WGWQ%GRID%NPLWV)
       ENDIF

       DATAKE(NI)=SQRT(POTFAK)
       ! maximum kinetic energy
       IF(DATAKE(NI)>DATAKEMAX) THEN
          DATAKEMAX=DATAKE(NI)
       ENDIF

    ENDDO

    IF ( LESF_SPLINES ) THEN
       !store index order
       COR%ESF(I)%OLD_INDEX = 0 
       COR%ESF(I)%OLD_INDEX( 1:NP ) = OLD_INDEX(1:NP )
       !store attenuated coulomb potential 
       COR%ESF(I)%V = 0 
       COR%ESF(I)%V( 1:NP ) = DATAKE(1:NP)
    ENDIF
!$ACC ENTER DATA COPYIN(DATAKE,OLD_INDEX) __IF_ASYNC__

    PROFILING_STOP ('loop1|xi_local_field_acfdt_gg')
    ! setup of multipole corrections
    IF (MCALPHA/=0) THEN
!note(sm): this is code that would fail below anyway, not porting to OpenACC
       CALL FOCK_MULTIPOLE_CORR_SETUP(LATT_CUR, GRIDHF)
       CALL GW_MULTIPOLE_PROJ_SETUP( WGWQ, LATT_CUR, DATAKE,DATAKEMAX )
    ENDIF

!=======================================================================
! now multiply the response function with the kernel
! from left and right hand side
!=======================================================================

    PROFILING_START('loop2|xi_local_field_acfdt_gg')
    !save chi to chi_work
    ! change sign of RESPONSEFUN  here 
    !CHI_WORK=0._q
    ! for some reason NEC compiler wants this loop to be spelled out 
    IF (CHI%LREALSTORE) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CHI_WORK,CHI) __IF_ASYNC__
       DO NI_= 1, NCOLS
          DO NI = 1, NROWS
             CHI_WORK(NI,NI_) = -CHI%RESPONSER( NI, NI_ , NOMEGA_LOCAL )
          ENDDO
       ENDDO
    ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CHI_WORK,CHI) __IF_ASYNC__
       DO NI_= 1, NCOLS
          DO NI = 1, NROWS
             CHI_WORK(NI,NI_) = -CHI%RESPONSEFUN( NI, NI_ , NOMEGA_LOCAL )
          ENDDO
       ENDDO
    ENDIF
    !kill elements larger than NP
    IF ( NP < NROWS ) THEN
!$ACC PARALLEL LOOP PRESENT(DATAKE) __IF_ASYNC__
       DO NI = NP + 1, NROWS
          DATAKE(NI) = 0
       ENDDO
    ENDIF

    !DATAKE(:) has the coulomb potenial, this needs to be multiplied with CHI_WORK
    !and symmetrized
!$ACC PARALLEL LOOP PRESENT(DATAKE,CHI_WORK) PRIVATE(POTFAK) __IF_ASYNC__
    DO NI=1,NROWS  !from the left
       POTFAK=DATAKE( NI )
       CHI_WORK(NI,1:NCOLS)=CHI_WORK(NI,1:NCOLS)*POTFAK
    ENDDO

    !from right
!$ACC PARALLEL LOOP PRESENT(DATAKE,CHI_WORK,COMM_INOMEGA) PRIVATE(POTFAK) __IF_ASYNC__
    DO NI=1,NCOLS
       POTFAK=DATAKE( NI + (COMM_INOMEGA%NODE_ME-1)*NCOLS )
       !IF ( NI + (COMM_INOMEGA%NODE_ME-1)*NCOLS > NP ) POTFAK = 0 
       CHI_WORK(1:NROWS,NI)=CHI_WORK(1:NROWS,NI)*POTFAK
    ENDDO


    ! multipole corrections:
    IF(MCALPHA/=0) THEN 

      CALL vtutor%error("Sorry MCALPHA not implemented")

!note(sm): The following is dead code. Won't port to OpenACC
       DO NI=1,NP
          
          TSUM=0.0
          DO NK=1,10
             DO NJ=1,NP
                TSUM(NK)=TSUM(NK)+CHI_WORK(NI,NJ)*CONJG(P_PROJ2(NJ,NK))                
             ENDDO
          ENDDO

          POTCORRECTION=(MATMUL(AAH_MAT,TSUM))

          DO NK=1,10
             DO NJ=1,NP
                CHI_WORK(NI,NJ)=CHI_WORK(NI,NJ)+P_PROJ2(NJ,NK)*POTCORRECTION(NK)
             ENDDO
          ENDDO
          
       ENDDO
       
       
       DO NI=1,NP
          
          TSUM=0.0
          DO NK=1,10
             DO NJ=1,NP
                TSUM(NK)=TSUM(NK)+CHI_WORK(NJ,NI)*P_PROJ2(NJ,NK)
             ENDDO
          ENDDO

          POTCORRECTION=(MATMUL(AAH_MAT,TSUM))

          DO NK=1,10
             DO NJ=1,NP
                CHI_WORK(NJ,NI)=CHI_WORK(NJ,NI)+CONJG(P_PROJ2(NJ,NK))*POTCORRECTION(NK)
             ENDDO
          ENDDO
          
       ENDDO
! When using MCALPHA POTFAK had an additional normalisation factor to ensure consistency with fock.F
! acfdt does not need this extra factor so we have to take it out again:
       CHI_WORK=CHI_WORK*(WGWQ%GRID%NPLWV)
    ENDIF

    !smooth energy cutoff done by replacing the NI^th coulmn by the OLD_INDEX(NI)^th column 
    !and analogously for the row. 

    !The replacement of the rows is straight forward 
    !and is be done directly by the following loop
!$ACC PARALLEL LOOP SEQ PRIVATE(NO) PRESENT(OLD_INDEX,CHI_WORK) __IF_ASYNC__ VECTOR_LENGTH(32)
    DO NI=1,MAXINDEX
      NO = OLD_INDEX(NI)
!$ACC LOOP GANG VECTOR
      DO NK=1,NCOLS
        CHI_WORK(NI,NK)=CHI_WORK(NO,NK)
        IF (NI.ne.NO) CHI_WORK(NO,NK) =0._q
      ENDDO
    ENDDO
    !for the columns this is not straight forward
    !since the whole matrix is distributed in stripes of NCOLS
    !among the nodes in each frequency group.

    PROFILING_STOP ('loop2|xi_local_field_acfdt_gg')
    PROFILING_START('comms|xi_local_field_acfdt_gg')

    !send columns between nodes
    ALLOCATE(CHI_SEND( NROWS ) )  
    ALLOCATE(CHI_RECV( NROWS ) ) 
!$ACC ENTER DATA CREATE(CHI_SEND,CHI_RECV) __IF_ASYNC__
!$ACC KERNELS PRESENT(CHI_RECV,CHI_SEND) __IF_ASYNC__
    CHI_RECV = 0 
    CHI_SEND = 0 
!$ACC END KERNELS

    replace_cols: IF( COMM_INOMEGA%NCPU == 1) THEN

!$ACC PARALLEL LOOP SEQ PRIVATE(NO) PRESENT(OLD_INDEX,CHI_WORK) __IF_ASYNC__ VECTOR_LENGTH(32)
       DO NI=1,MAXINDEX
          NO=OLD_INDEX(NI)
!$ACC LOOP GANG VECTOR
          DO NK=1,MAXINDEX
             CHI_WORK(NK,NI)=CHI_WORK(NK,NO)
             IF (NI/=NO) CHI_WORK(NK,NO)=0._q
          ENDDO
       ENDDO

    ELSE replace_cols

       CALLMPI( M_barrier(COMM_INOMEGA)) !each node should start at the same time from here
       !==============================================================================================
       !we loop over the cpus in one group and let other node send the appropriate columns
       DO NODE_RECV = 1, COMM_INOMEGA%NCPU
       !==============================================================================================
           loc_col : DO NI = 1 , NCOLS                      !search column which should be replaced
           !------------------------------------------------------------------------------------------

              NI_GLOB = NI + (NODE_RECV-1)*NCOLS   !global index for current node NODE_RECV

              !find node NODE_SEND containing OLD_INDEX(NI_GLOB)
              DO NODE_SEND = 1, COMM_INOMEGA%NCPU
                IF ( 1 + ( NODE_SEND-1)*NCOLS <= OLD_INDEX(NI_GLOB) .AND. &
                     OLD_INDEX(NI_GLOB) <= NODE_SEND*NCOLS) EXIT
              ENDDO

              NSEND = OLD_INDEX(NI_GLOB) - (NODE_SEND-1)*NCOLS   !local column index on sending node

              IF ( NSEND < 1 ) THEN
                 CALL vtutor%error("Cannot send column with negative index  "&
                    // str(NODE_RECV) // " " // str(NODE_SEND) // " " // str(NI_GLOB) // " " // &
                    str(OLD_INDEX(NI_GLOB)) // " " // str(NSEND) // " " // str(NCOLS))
              ENDIF

              !if old column is on the same node simply replace it
              IF ( NODE_SEND == COMM_INOMEGA%NODE_ME .AND. NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN
                 IF (NSEND.ne.NI) THEN
!$ACC KERNELS PRESENT(CHI_WORK) __IF_ASYNC__
                    CHI_WORK(1:NROWS, NI) = CHI_WORK(1:NROWS, NSEND)
                    !jK added
                    CHI_WORK(:,NSEND) = 0._q
!$ACC END KERNELS
                 ENDIF

              !otherwise tell NODE_SEND to send the appropriate column NSEND to NODE_RECV
              ELSEIF ( NODE_RECV < NODE_SEND ) THEN
                 IF ( NODE_SEND == COMM_INOMEGA%NODE_ME ) THEN
!$ACC KERNELS PRESENT(CHI_SEND,CHI_WORK) __IF_ASYNC__
                    CHI_SEND(1:NROWS) = CHI_WORK(1:NROWS, NSEND)
                    CHI_WORK(:,NSEND) = 0._q
!$ACC END KERNELS

                    !send this column to node NODE_RECV
                    IF (CHI%LREALSTORE) THEN
                       CALLMPI ( M_send_d( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                    ELSE
                       CALLMPI ( M_send_z( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                    ENDIF

                 ELSEIF ( NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN
                    IF (CHI%LREALSTORE) THEN
                       CALLMPI ( M_recv_d( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                    ELSE
                       CALLMPI ( M_recv_z( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                    ENDIF

                    !and replace CHI_RECV
!$ACC KERNELS PRESENT(CHI_WORK,CHI_RECV) __IF_ASYNC__
                    CHI_WORK(1:NROWS, NI ) = CHI_RECV(1:NROWS)
!$ACC END KERNELS
                 ENDIF
              ENDIF

              !mpi barrier here
              CALLMPI( M_barrier(COMM_INOMEGA) )
           !------------------------------------------------------------------------------------------
           ENDDO loc_col                           !local columns
       !==============================================================================================
       ENDDO                                    !cpus in group
    ENDIF replace_cols
    !==============================================================================================
!$ACC EXIT DATA DELETE(CHI_SEND,CHI_RECV) __IF_ASYNC__
    DEALLOCATE( CHI_SEND )
    DEALLOCATE( CHI_RECV )
    PROFILING_STOP ('comms|xi_local_field_acfdt_gg')

!$ACC EXIT DATA DELETE(DATAKE,OLD_INDEX) __IF_ASYNC__
    DEALLOCATE(DATAKE, OLD_INDEX)

    PROFILING_STOP('xi_local_field_acfdt_gg')

 END SUBROUTINE XI_LOCAL_FIELD_ACFDT_GG


SUBROUTINE COL2BC_AND_RPA(CHI_WORK, GDES, B, NP_NEW, COR,NQ, OMEGAWEIGHT, NCUT, IO)

GPACC USE moffload_struct_def

    USE base
    USE constant
    USE mpimy 
    USE ini
    USE scala
    USE chi_glb, ONLY: LLTDMP2 
    IMPLICIT NONE
    GDEF,POINTER,CONTIGUOUS :: CHI_WORK(:,:)        !< matrix that is redistributed
    TYPE (greensfdes)       :: GDES                 !< greensfunction descriptor
    TYPE (loop_des)         :: B                    !< imaginary grids for bosons
    INTEGER, INTENT(IN)     :: NP_NEW               !< dimension of submatrix
    TYPE (correlation)      :: COR                  !< correlation energy
    INTEGER, INTENT(IN)     :: NQ                   !< current q-point
    REAL(q), INTENT(IN)     :: OMEGAWEIGHT          !< integeration weight
    INTEGER, INTENT(IN)     :: NCUT                 !< current cut off 
    TYPE (in_struct)        :: IO                   !< IO 
    !local                
    GDEF, ALLOCATABLE       :: CHI_TMP(:)
    INTEGER                 :: I,J                  !some loop variables 
    REAL(q), ALLOCATABLE    :: E(:)                 !eigenvectors and global eigenvalues
    REAL(q)                 :: RTMP                 !auxilairy 
    REAL(q)                 :: SUM2                 !stores sum of squares of CHI_WORK elements 
    REAL(q)                 :: EXX                  !stores exact exchange part 
                                                    !to obtain dMP2 contribution
    LOGICAL                 :: LRPADIAG             !use diagonalisation or Cholesky
    INTEGER                 :: ISTAT 
    TYPE(greens_mat_des), POINTER :: RDES1D=>NULL() ! distributed response col.  major distributed
    TYPE(greens_mat_des), POINTER :: RDES2D=>NULL() ! distributed response on block-cyclic 2D grid

    PROFILING_START( 'col2bc_and_rpa' )

#ifdef ACC_OFFLOAD
    IF (OFFLOAD_ON) THEN
       CALL CUSOLVER_RPA(CHI_WORK, GDES, B, NP_NEW, COR, NQ,&
          OMEGAWEIGHT, NCUT, IO )
       PROFILING_STOP( 'col2bc_and_rpa' )
       RETURN
    ENDIF
#endif

    ! if linear term is calculated use diagonalization routine
    ! usually linear term is not calculated, so Cholesky decomposition is used
    ! by default
    ! diagonalization is also required if electronic structure factor is interpolated 
    !LRPADIAG = ( LDMP1 .OR. LESF_SPLINES )
    LRPADIAG = ( LESF_SPLINES .AND. .NOT. LDMP1 )
    !LRPADIAG =  LDMP1 
    IF( LESF_SPLINES ) THEN
       COR%ESF(NCUT)%NG(NQ) = NP_NEW 
    ENDIF

    !LRPADIAG=.FALSE.               ! get energy from Cholesky for RPA and from \sum |A_ij|^2 for dMP2 
    !LRPADIAG=.TRUE.               ! get energy from diagonalisation 

    
    !! initialize ScaLAPACK descriptor ( for CHI_WORK, i.e. column major distribution)
    CALL RESPONSE_DES_INIT_COLMAJ( RDES1D, B, GDES )
    !! initialize ScaLAPACK descriptor ( for CHI_TMP that is optimal for scaLAPACK, i.e. block cyclic)
    CALL RESPONSE_DES_INIT_BC( RDES2D, B, GDES )

    !---------------------------------------------------------------------------------------
    !initialize 1D column-major processor grid
    !the following is only neccessary if we are not in dummy mode
    IF( B%LDO_POINT_LOCAL ) THEN

       PROFILING_START('redis|col2bc_and_rpa')

       !allocate storage for redistributed matrix
       ALLOCATE(CHI_TMP(RDES2D%MY_NROWS*RDES2D%MY_NCOLS), STAT = ISTAT) 
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "COL2BC_AND_RPA is not able to allocate "//&
            str( ndata*8._q*RDES2D%MY_NROWS*RDES2D%MY_NCOLS) //&
            " kB of data on MPI rank"//str( RDES2D%COMM%NODE_ME ) )
       ENDIF
       CALL REGISTER_ALLOCATE(ndata*8._q*SIZE(CHI_TMP,KIND=qi8), "ACFDT_chi")
       CHI_TMP = 0 

       !now call the BLACS redistribution routine 
       CALL PGGEMR2D( RDES1D%NROWS, RDES1D%NROWS, CHI_WORK, 1, 1, RDES1D%DESC, & 
                      CHI_TMP, 1, 1, RDES2D%DESC, RDES2D%ICTXT )

       PROFILING_STOP('redis|col2bc_and_rpa')

       PROFILING_START('diag|col2bc_and_rpa')

       ALLOCATE( E( RDES2D%NROWS) )
       E=0
       IF (LRPADIAG) THEN
          ! exchange part needs no diagonalization, only diagonal of X.V
          IF ( ABS(OMEGAWEIGHT) < 1.E-8_q .OR. LDMP1 ) THEN
             E=0
             CALL DETERMINE_DIAGONALE_GDEF_REAL ( NP_NEW, CHI_TMP, E, RDES2D%DESC ) 
             CALLMPI( M_sum_d( RDES2D%COMM_INTRA, E(1), NP_NEW ) )
             RTMP=1
             IF( LDMP1 ) THEN
                RTMP = OMEGAWEIGHT 
             ELSE IF ( NQ == 1 )  THEN
                RTMP=RTMP/IDIR_MAX
             ENDIF 
             !this is the trace (note that E has wrong sign)    
             DO I=1, NP_NEW
                COR%CORRMP2DIR_K( NCUT ) = COR%CORRMP2DIR_K( NCUT ) + E(I)*RTMP
             ENDDO
             ! calculate electronic structure factor for exchange 
             ! S_G(q) =[ X.V ]
             IF( LESF_SPLINES ) THEN
                DO I=1, NP_NEW
                   COR%ESF(NCUT)%S( I, NQ ) = COR%ESF(NCUT)%S(I, NQ ) + E(I)*RTMP
                ENDDO
             ENDIF
          ELSE
             ! call scaLAPACK diagonalizer for hermitian and symmetric matrices respectively
             CALL PDSYEV_ZHEEV_DESC( CHI_TMP, E, NP_NEW, RDES2D%DESC, RDES2D%COMM_INTRA )
             ! calculate electronic structure factor  
             ! S_G(q) = V^{-1/2}.U^+.[ ln( 1-X.V ) - X.V ].U.V^{-1/2}
             IF( LESF_SPLINES ) THEN
                   CALL EXTRACT_STUFAK_DISTRI( NP_NEW, CHI_TMP, E, RDES2D, NQ, OMEGAWEIGHT,&
                       RTMP, SUM2, COR%ESF(NCUT) ) 
                   ! accumulate results 
                   COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + RTMP
                   COR%CORRMP2DIR_K( NCUT )  = COR%CORRMP2DIR_K( NCUT ) - SUM2
             ELSE
                !this is the trace (note that E has wrong sign)    
                DO I=1, NP_NEW
                   IF ( 1 + E(I) > 1E-8_q .AND. ABS(OMEGAWEIGHT) > 1.E-8_q  ) THEN
                      !RPA (Pines & Noizeres respectively Gell-Mann & Brueckner)
                      RTMP = ( LOG( 1 + E( I ) ) - E( I ) ) * OMEGAWEIGHT
                      COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + RTMP
   
                      !direct Moeller Plessett term of second order
                      RTMP = ( E( I )*E( I ) / 2 ) * OMEGAWEIGHT
                      COR%CORRMP2DIR_K( NCUT )  = COR%CORRMP2DIR_K( NCUT ) - RTMP
                   ELSE
                      WRITE(*,'(" ERROR in COL2BC_AND_RPA, eigenvalue #",I4," is negative '//&
                      'for group",I3,":",F12.5)') I, RDES2D%COMM_INTRA%NODE_ME, E( I )
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
       ELSE
          ! obtain the dMP1 energy here, this should be the EXHF part 
          IF ( ABS(OMEGAWEIGHT) < 1.E-8_q .OR. LDMP1 ) THEN
             E=0
             CALL DETERMINE_DIAGONALE_GDEF_REAL ( NP_NEW, CHI_TMP, E, RDES2D%DESC ) 
             CALLMPI( M_sum_d( RDES2D%COMM_INTRA, E(1), NP_NEW ) )

             RTMP=1
             IF( LDMP1 ) THEN
                RTMP = OMEGAWEIGHT 
             ELSE IF ( NQ == 1 )  THEN
                RTMP=RTMP/IDIR_MAX
             ENDIF 
             !this is the trace (note that E has wrong sign)    
             SUM2 = 0 
             DO I=1, NP_NEW
                SUM2= SUM2 + E(I)*RTMP
             ENDDO
             IF( LDMP1 ) THEN
                COR%CORRMP2DIR_K( NCUT ) = COR%CORRMP2DIR_K( NCUT ) + SUM2
             ENDIF
             ! calculate electronic structure factor for exchange 
             ! S_G(q) =[ X.V ]
             IF( LESF_SPLINES ) THEN
                DO I=1, NP_NEW
                   COR%ESF(NCUT)%S( I, NQ ) = COR%ESF(NCUT)%S(I, NQ ) + E(I)*RTMP
                ENDDO
             ENDIF
          ELSE
             ! obtain the dMP2 energy here
             ! this is essentially the squared Frobenius norm of CHI_TMP
             CALL MATRIX_NORM_DESC( CHI_TMP, NP_NEW, RDES2D%DESC, SUM2 )
             SUM2=0.5_q*SUM2* OMEGAWEIGHT
             COR%CORRMP2DIR_K( NCUT ) = COR%CORRMP2DIR_K( NCUT ) - SUM2
          ENDIF
          ! Use the well-known identity of the logarithm of the 
          ! determinant |A| of a positive definite matrix A 
          !            Ln | A | = Tr Ln A  
          ! in combination with the Cholesky decomposition 
          !    A = L^T . L  -> | A | = |L^T . L | = |L^T| |L| = |L|^2 
          ! to  rewrite  
          !    Tr Ln( 1 + A ) =  Ln |1+A| = 2 * Ln | L' | = 2 * Ln \Prod_i^N L'_ii
          !                   = 2 Tr Ln diag( L' )
          ! where L' is the Cholesky decomposition of 
          !    L'^T . L' = 1 + A 
          
          ! obtain trace of CHI_TMP for subtraction later on
          CALL DETERMINE_DIAGONALE_GDEF_REAL ( NP_NEW, CHI_TMP, E, RDES2D%DESC ) 
          CALLMPI( M_sum_d( RDES2D%COMM_INTRA, E(1), NP_NEW ) )
          EXX = SUM( E(1:NP_NEW ) )
          ! set up  1 + A 
          E(MIN( NP_NEW, RDES2D%NROWS):) = 0 
          E(1:NP_NEW) = 1 
          CALL ADD_TO_DIAGONALE_GDEF( NP_NEW, CHI_TMP, E, RDES2D%DESC )
          ! obtain Cholesky decompisiton
          CALL PDPOTRF_PZPOTRF_DESC( CHI_TMP, NP_NEW, RDES2D%DESC ) 
          !now, determine Tr Ln 1+A as  2* Tr Ln diag( L )
          CALL DETERMINE_DIAGONALE_GDEF_REAL ( NP_NEW, CHI_TMP, E, RDES2D%DESC ) 
          CALLMPI( M_sum_d( RDES2D%COMM_INTRA, E(1), NP_NEW ) )
          RTMP = 0 
          DO I = 1, NP_NEW  
             ! dump warning if Cholesky eigenvalue is negative
             IF ( E(I) > 0 ) THEN
                RTMP= RTMP + 2 * LOG( E(I) ) 
             ELSE IF ( ABS(OMEGAWEIGHT) > 1.E-8 .AND. E(I) < 0 ) THEN
                WRITE(*,*)"WARNING: eigenvalue of 1+X.V is negative", NCUT, I, E(I)
             ENDIF
          ENDDO
          !integrate over frequencies 
          RTMP=(RTMP-EXX)* OMEGAWEIGHT
          COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + RTMP 
       ENDIF

       DEALLOCATE(E)

       CALL DEREGISTER_ALLOCATE(ndata*8._q*SIZE(CHI_TMP,KIND=qi8), "ACFDT_chi")
       DEALLOCATE(CHI_TMP)

       PROFILING_STOP('diag|col2bc_and_rpa')
    ENDIF

    !exit blacs  
    CALL DESTROY_GDES_MAT( RDES2D )
    CALL DESTROY_GDES_MAT( RDES1D )

    ! in case laplace transformed direct MP2 was done, the RPA energy does not make sense
    IF( LLTDMP2 ) COR%CORRELATION_K( NCUT ) = 0 

    DEALLOCATE ( CHI_WORK )
    NULLIFY(CHI_WORK )
    PROFILING_STOP( 'col2bc_and_rpa' )

END SUBROUTINE COL2BC_AND_RPA


#ifdef ACC_OFFLOAD
!****************************************************************************************
!
!> redistributes column wise distributed CHI_WORK  2D block-cyclic wise 
!> also the RPA and MP2 correlation energy is determined. (openacc version)
!
!****************************************************************************************
SUBROUTINE CUSOLVER_RPA(CHI_WORK, GDES, B, NP_NEW, COR,NQ, OMEGAWEIGHT, NCUT, IO)

GPACC USE moffload

    USE base
    USE constant
    USE mpimy
    USE ini
    USE scala
    USE chi_glb, ONLY: LLTDMP2
    IMPLICIT NONE
    GDEF,POINTER,CONTIGUOUS :: CHI_WORK(:,:)        !< matrix that is redistributed
    TYPE (greensfdes)       :: GDES                 !< greensfunction descriptor
    TYPE (loop_des)         :: B                    !< loop descriptor for bosons
    INTEGER, INTENT(IN)     :: NP_NEW               !< dimension of submatrix
    TYPE (correlation)      :: COR                  !< correlation energy
    INTEGER, INTENT(IN)     :: NQ                   !< current q-point
    REAL(q), INTENT(IN)     :: OMEGAWEIGHT          !< integeration weight
    INTEGER, INTENT(IN)     :: NCUT                 !< current cut off 
    TYPE (in_struct)        :: IO                   !< IO 
    !local
    GDEF, POINTER, DIMENSION(:,:)  :: CHI_TMP(:,:)
    GDEF, POINTER, DIMENSION(:,:)  :: CHI_WORK_ORIG(:,:)
    INTEGER                 :: I,J                  !some loop variables
    REAL(q), ALLOCATABLE    :: E(:)                 !eigenvectors and global eigenvalues
    INTEGER                 :: PINFO                !scaLAPACK routine info variable
    INTEGER                 :: LWORK, LRWORK        !for scaLAPACK routine
    GDEF, ALLOCATABLE       :: WORK(:)              !working array
    COMPLEX(q), ALLOCATABLE :: RWORK(:)             !needed for complex computation
    REAL(q)                 :: RTMP                 !auxilairy
    REAL(q)                 :: SUM2                 !stores sum of squares of CHI_WORK elements
                                                    !to obtain dMP2 contribution
    COMPLEX(q)              :: CHI_TRACE            !trace of CHI_WORK to get RPA energy
    LOGICAL                 :: LRPADIAG             !use diagonalisation or Cholesky
    COMPLEX(q)              :: CTMPREL, CTMPRMP2

    INTEGER :: NCOLS, NROWS
    GDEF, ALLOCATABLE :: CTMP(:,:)
    GDEF, ALLOCATABLE :: CEIDB(:,:)
    
    REAL(q), ALLOCATABLE :: V(:) 
    REAL(q), POINTER :: V_ORIG(:) => NULL()
    REAL(q), POINTER, CONTIGUOUS :: S(:) => NULL()
    INTEGER, POINTER, CONTIGUOUS :: OLD_INDEX(:)  => NULL()

    PROFILING_START( 'cusolver_rpa' )

    ! if linear term is calculated use diagonalization routine
    ! usually linear term is not calculated, so Cholesky decomposition is used
    ! by default
    !LRPADIAG = LDMP1
    LRPADIAG = ( LDMP1 .OR. LESF_SPLINES )

    !LRPADIAG=.FALSE.               ! get energy from Cholesky for RPA and from \sum |A_ij|^2 for dMP2
    !LRPADIAG=.TRUE.               ! get energy from diagonalisation

    IF( B%LDO_POINT_LOCAL ) THEN

    PROFILING_START('diag|cusolver_rpa')

    NROWS = GDES%RES_NRPLWV_ROW_DATA_POINTS
    NCOLS = GDES%RES_NRPLWV_COL_DATA_POINTS

    IF ( LESF_SPLINES ) THEN
       COR%ESF(NCUT)%NG( NQ )  = NP_NEW
       ALLOCATE( V( NP_NEW ) ) 
!$ACC ENTER DATA CREATE(V) __IF_ASYNC__
       S => COR%ESF(NCUT)%S(:,NQ)
       OLD_INDEX => COR%ESF(NCUT)%OLD_INDEX
       V_ORIG => COR%ESF(NCUT)%V
!$ACC ENTER DATA COPYIN(S, OLD_INDEX,V_ORIG) __IF_ASYNC__
    ENDIF

    NCOLS = NCOLS*B%COMM_IN_GROUP%NCPU
    ! gather complete matrix of frequency group to rank 1 
    IF ( B%COMM_IN_GROUP%NCPU > 1 ) THEN
       ALLOCATE(CHI_TMP(NROWS,NCOLS))
!$ACC ENTER DATA CREATE(CHI_TMP) __IF_ASYNC__
       CALL M_gathero_z(B%COMM_IN_GROUP, size(CHI_WORK,1)*size(CHI_WORK,2), CHI_WORK, CHI_TMP)
       CHI_WORK_ORIG=>CHI_WORK
       CHI_WORK=>CHI_TMP
    ENDIF

    IF (LRPADIAG) THEN
       IF ( B%COMM_IN_GROUP%NODE_ME==1 ) THEN
       !allocate eigenvalue and eigenvector arrays
       ALLOCATE( E( NROWS) ) ;  E = 0
!$ACC ENTER DATA COPYIN(OMEGAWEIGHT) CREATE(E) __IF_ASYNC__

       IF (.NOT. OFFLOAD_ON ) THEN
          ALLOCATE(WORK(1))
       !first get the optimal size of working array
#ifdef gammareal
          !for real values only
          CALL DSYEV( 'V', 'U', NROWS, CHI_WORK(1,1), NROWS, E(1),&
               WORK(1), -1, PINFO)
#else
          !for a hermitian matrix
          ALLOCATE(RWORK(1))
          CALL ZHEEV( 'V', 'U', NROWS, CHI_WORK(1,1), NROWS, E(1), &
               WORK(1), -1, RWORK(1), PINFO )
#endif

          IF ( PINFO /=0 ) THEN
             CALL vtutor%error("CUSOLVER_RPA: for NCUT " // str(NCUT) // " 1st call of &
                &DSYEV returns " // str(PINFO))
          ENDIF

          !LWORK is now the optimal size of WORK
          LWORK=WORK(1)
          IF ( WORK(1)  /= 0 ) THEN   !if it's not zero reallocate WORK array
             DEALLOCATE(WORK)
             ALLOCATE(WORK(LWORK))
          ENDIF
#ifndef gammareal
          !in the complex version do the same with the RWORK array
          LRWORK=MAX(1,3*NROWS-2)
          IF ( LRWORK  /= 0 ) THEN   !if it's not zero reallocate WORK array
            DEALLOCATE(RWORK)
            ALLOCATE(RWORK(LRWORK))
          ENDIF
#endif 
       ENDIF

       !diagonalize
#ifdef gammareal
       !real version
       CALL WDSYEV( 'V', 'U', NROWS, CHI_WORK(1,1), NROWS, E(1),&
            WORK(1), LWORK , PINFO)
#else
       !hermitian version
       CALL WZHEEV( 'V', 'U', NROWS, CHI_WORK(1,1), NROWS, E(1),&
            WORK(1), LWORK, RWORK(1), PINFO )
#endif
       IF ( PINFO /= 0 ) THEN
          CALL vtutor%error("CUSOLVER_RPA: for NCUT " // str(NCUT) // " 2nd call of &
             &DSYEV/ZHEEV returns " // str(PINFO))
       ENDIF
      ENDIF ! master rank


    ELSE
       !!!!!!!
       ! obtain the dMP2 energy here
       !!!!!!!

       SUM2=0._q
       IF ( B%COMM_IN_GROUP%NODE_ME==1 ) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) REDUCTION(+:SUM2) PRESENT(CHI_WORK) __IF_ASYNC__
          DO J=1,NCOLS
             DO I=1,NROWS
                SUM2=SUM2+CHI_WORK(I,J)*GCONJG(CHI_WORK(I,J))
             ENDDO
          ENDDO
       ENDIF

       ! use Cholesky decomposition
       !calculate trace
       !add one to the diagonal
       CHI_TRACE=0
       RTMP=0
       IF ( B%COMM_IN_GROUP%NODE_ME==1 ) THEN
!$ACC PARALLEL LOOP PRESENT(CHI_WORK) REDUCTION(+:CHI_TRACE) __IF_ASYNC__
          DO I=1, NROWS
             CHI_TRACE = CHI_TRACE + CHI_WORK(I,I)
             CHI_WORK(I,I) = CHI_WORK(I,I) + 1
          ENDDO
#ifdef gammareal
          !real version
          !call Cholesky
          CALL WDPOTRF( 'U', NROWS, CHI_WORK(1,1), NROWS, PINFO)
#else
          !complex version
          !call Cholesky
          CALL WZPOTRF( 'U', NROWS, CHI_WORK(1,1), NROWS, PINFO)
#endif
       !do log of diagonal
!$ACC PARALLEL LOOP PRESENT(CHI_WORK) REDUCTION(+:RTMP) __IF_ASYNC__
          DO I=1, NROWS
             RTMP = RTMP + LOG(REAL(CHI_WORK(I,I),q))
          ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          RTMP=RTMP*2

          !energy obtained as 2\sum_i ln(A)_{ii}
          RTMP=(RTMP-CHI_TRACE)* OMEGAWEIGHT
          COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + DBLE(RTMP)

          SUM2=0.5_q*SUM2* OMEGAWEIGHT
          COR%CORRMP2DIR_K( NCUT ) = COR%CORRMP2DIR_K( NCUT ) - SUM2
       ENDIF
    ENDIF

    IF (LRPADIAG) THEN
       IF ( B%COMM_IN_GROUP%NODE_ME==1 ) THEN
          CTMPREL  = 0
          CTMPRMP2 = 0
!$ACC PARALLEL LOOP REDUCTION(+:CTMPREL) REDUCTION(-:CTMPRMP2) PRESENT(E,OMEGAWEIGHT,CTMP,CHI_WORK) PRIVATE(RTMP) __IF_ASYNC__
          DO I=1, NROWS
             IF ( 1 + E(I) > 1E-8_q  ) THEN
                !RPA (Pines & Noizeres respectively Gell-Mann & Brueckner)
                RTMP = ( LOG( 1 + E( I ) ) - E( I ) ) * OMEGAWEIGHT
                CTMPREL = CTMPREL + RTMP

                ! linear term
                IF ( LDMP1 ) THEN
                   RTMP = ( E( I ) ) * OMEGAWEIGHT
                !direct Moeller Plessett term of second order
                ELSE
                   RTMP = ( E( I )*E( I ) / 2 ) * OMEGAWEIGHT
                ENDIF
                CTMPRMP2 = CTMPRMP2 - RTMP
             ELSE
                WRITE(*,*) 'Error in CUSOLVER_RPA, eigenvalue #', I, 'is negative for group', &
                           B%COMM_BETWEEN_GROUPS%NODE_ME, ': ', E(I)
             ENDIF
          ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
          COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + CTMPREL
          COR%CORRMP2DIR_K( NCUT )  = COR%CORRMP2DIR_K( NCUT ) + CTMPRMP2
          
! optionally, collect electronic structure factor 
          IF ( LESF_SPLINES ) THEN
!$ACC PARALLEL LOOP PRIVATE( J ) PRESENT(V,V_ORIG,OLD_INDEX) __IF_ASYNC__
             DO I = 1, NP_NEW 
                J=OLD_INDEX(I)
                !V(I) = COR%ESF(NCUT)%V( COR%ESF(NCUT)%OLD_INDEX(I) )
                V(I) = V_ORIG( J )
             ENDDO

             ALLOCATE(CTMP(NROWS, NCOLS ))
             ALLOCATE(CEIDB(NROWS, NCOLS ))
!$ACC ENTER DATA CREATE(CTMP,CEIDB) __IF_ASYNC__
!$ACC PARALLEL LOOP PRESENT(E,OMEGAWEIGHT,CTMP,CHI_WORK) __IF_ASYNC__
             DO I=1, NROWS
                IF ( 1 + E(I) > 1E-8_q  ) THEN
                   CTMP(1:NROWS,I)=CHI_WORK(1:NROWS, I)*(LOG(1+E(I))-E(I))*OMEGAWEIGHT/2
                ELSE
                   CTMP(:,I)=0
                ENDIF
             ENDDO
             CALL WGGEMM( 'N', trans, NROWS, NROWS, NROWS, one, CTMP, &
                  &             NROWS, CHI_WORK(1,1), NROWS, zero, CEIDB, NROWS)
!$ACC KERNELS PRESENT(CEIDB,V,S) __IF_ASYNC__
             DO  I = 1, NP_NEW 
                S( I ) =S( I ) + CEIDB( I, I )/V(I)**2
             ENDDO
!$ACC END KERNELS

!$ACC EXIT DATA DELETE(CTMP,CEIDB) __IF_ASYNC__
             DEALLOCATE(CTMP, CEIDB)
!$ACC EXIT DATA DELETE(V_ORIG, OLD_INDEX,V) __IF_ASYNC__
!$ACC EXIT DATA COPYOUT(S) __IF_ASYNC__
             DEALLOCATE( V ) 
          ENDIF

!$ACC EXIT DATA DELETE(E,OMEGAWEIGHT) __IF_ASYNC__
          DEALLOCATE(E)
       ENDIF
    ENDIF

    IF ( B%COMM_IN_GROUP%NCPU > 1 ) THEN
       CHI_WORK=>CHI_WORK_ORIG
!$ACC EXIT DATA DELETE(CHI_TMP) __IF_ASYNC__
       DEALLOCATE(CHI_TMP)
    ENDIF

    PROFILING_STOP('diag|cusolver_rpa')

    ENDIF

    ! in case laplace transformed direct MP2 was done, the RPA energy does not make sense
    IF( LLTDMP2 ) COR%CORRELATION_K( NCUT ) = 0

!$ACC EXIT DATA DELETE(CHI_WORK) __IF_ASYNC__
    DEALLOCATE ( CHI_WORK)     

    PROFILING_STOP( 'cusolver_rpa' )
END SUBROUTINE CUSOLVER_RPA
#endif

!******************************* LIN_REG_GG**********************************************
!
!> perform linear regression of correlation energy versus
!> 1/energy^(3/2)
!
!***************************************************************************************

  SUBROUTINE LIN_REG_GG(COR, E, IO)
  USE base 
  TYPE (correlation), POINTER :: COR
  TYPE( energy ) :: E
  TYPE( in_struct ) :: IO 
 
  REAL(q) :: SX, SY, SXY, SX2, SXM, SYM, SXYM, SX2M
  REAL(q) :: AREG, BREG, AREGM, BREGM
  INTEGER :: I, ILAMBDA
  REAL(q) :: LAMBDA
  !> for interpolation of electronic structure factor
  REAL(q) :: SXI, SYI, SXYI
  REAL(q) :: AREGI, BREGI
  CHARACTER(LEN=1) :: CN
  INTEGER :: INOUT 
  INOUT = IO%IU6
  
    SX = 0.0   
    SY = 0.0
    SXY = 0.0
    SX2 = 0.0
    SXM = 0.0   
    SYM = 0.0
    SXYM = 0.0
    SX2M = 0.0
  
    AREG = 0.0
    BREG = 0.0
    AREGM = 0.0
    BREGM = 0.0

    BREGI = 0
    AREGI = 0
    SYI = 0
    SXYI = 0

    CN="2"
    IF ( LDMP1 ) CN="1"

    IF ( LESF_SPLINES ) THEN
       IF (INOUT>=0) WRITE(INOUT,'(A)')"      cutoff energy     smooth cutoff   RPA   correlation   Hartree contr. to MP"//CN&
          //"  RPA spline-interp."
       IF (INOUT>=0) WRITE(INOUT,'("-----------------------------------------------------------------------------------------------------")')       
    ELSE
       IF (INOUT>=0) WRITE(INOUT,'(A)')"      cutoff energy     smooth cutoff   RPA   correlation   Hartree contr. to MP"//CN
       IF (INOUT>=0) WRITE(INOUT,'("---------------------------------------------------------------------------------")')       
    ENDIF


    ! set correlation energy with highest cutoff 
    IF ( L2ORDER ) THEN
       E%ERPA_CUT=REAL(COR%CORRMP2DIR(1),q)
    ELSE
       E%ERPA_CUT=REAL(COR%CORRELATION(1),q)
    ENDIF

    DO I=1,COR%NE
      IF ( LESF_SPLINES ) THEN
         CALL ESF_SPLINE_FIT( COR%ESF(I), .TRUE., IO ) 

         !> dump result 
         IF (INOUT>=0) WRITE(INOUT,'(" ",2F18.3,3F20.10)') COR%ENCUTGW(I),COR%ENCUTGWSOFT(I), &
                    REAL(COR%CORRELATION(I),q), REAL(COR%CORRMP2DIR(I),q) , REAL(  COR%ESF(I)%ESPLINED, q )
      ELSE
         !> dump result 
         IF (INOUT>=0) WRITE(INOUT,'(" ",2F18.3,3F20.10)') COR%ENCUTGW(I),COR%ENCUTGWSOFT(I), &
                    REAL(COR%CORRELATION(I),q), REAL(COR%CORRMP2DIR(I),q) 
      ENDIF
      IF ( COR%NE==1 ) CYCLE   

      SX = SX+ 1/COR%ENCUTGW(I)**1.5
      SY = SY+ REAL(COR%CORRELATION(I),q)
      SXY = SXY+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%CORRELATION(I),q)

      SX2 = SX2 + (1/COR%ENCUTGW(I)**1.5)**2
      SXM = SXM+ 1/COR%ENCUTGW(I)**1.5
      SYM = SYM+ REAL(COR%CORRMP2DIR(I),q)
      SXYM = SXYM+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%CORRMP2DIR(I),q)
      SX2M = SX2M + (1/COR%ENCUTGW(I)**1.5)**2      

      ! interpolate result with splines 
      IF( LESF_SPLINES ) THEN
         SYI  = SYI + REAL( COR%ESF(I)%ESPLINED, q)
         SXYI = SXYI+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%ESF(I)%ESPLINED,q)
      ENDIF
    ENDDO

    IF ( .NOT. COR%NE==1 ) THEN
       BREG =  (SY*SX - COR%NE*SXY)/(SX*SX - COR%NE*SX2) 
       AREG =  1/(COR%NE+0.0)*(SY - BREG*SX)
       BREGM = (SYM*SXM - COR%NE*SXYM)/(SXM*SXM - COR%NE*SX2M) 
       AREGM =  1/(COR%NE+0.0)*(SYM - BREGM*SXM) 

       IF (INOUT >=0) WRITE(INOUT,'("  linear regression    ")')  
       IF( .NOT. LESF_SPLINES ) THEN
          IF (INOUT >=0) WRITE(INOUT,'("  converged value                    ", 3F20.10)') AREG, AREGM
          E%ERPA_INF = AREG
       ELSE
          BREGI =  (SYI*SX - COR%NE*SXYI)/(SX*SX - COR%NE*SX2) 
          AREGI =  1/(COR%NE+0.0)*(SYI - BREGI*SX)
          IF (INOUT >=0) WRITE(INOUT,'("  converged value                    ", 3F20.10)') AREG, AREGM, AREGI
          E%ERPA_INF = AREGI
       ENDIF
    ELSE 
        E%ERPA_INF = REAL(COR%CORRELATION(1),q)
    ENDIF

    IF (NLAMBDA>=1) THEN
       IF (INOUT >=0) WRITE(INOUT,'("  converged value lambda= ",F7.4,"     ", 2F20.10)') 0.0_q,0.0_q
       DO ILAMBDA=0,NLAMBDA
       LAMBDA=1.0_q*(ILAMBDA+1)/MAX(1,NLAMBDA)
       SX = 0.0   
       SY = 0.0
       SXY = 0.0
       SX2 = 0.0

       DO I=1,COR%NE
          SX = SX+ 1/COR%ENCUTGW(I)**1.5
          SY = SY+ REAL(COR%CORRELATION_LAMBDA(I,ILAMBDA),q)
          SXY = SXY+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%CORRELATION_LAMBDA(I,ILAMBDA),q)
          SX2 = SX2 + (1/COR%ENCUTGW(I)**1.5)**2
       ENDDO

       BREG =  (SY*SX - COR%NE*SXY)/(SX*SX - COR%NE*SX2) 
       AREG =  1/(COR%NE+0.0)*(SY - BREG*SX)

       IF (INOUT >=0) WRITE(INOUT,'("  converged value lambda= ",F7.4,"     ", 2F20.10)') LAMBDA,AREG
       ENDDO
    ENDIF


  END SUBROUTINE LIN_REG_GG

#ifdef deprecated
!************************************ STORE_VXV_TO_CHI **********************************
!
!> stores the second order screened potential W^(2) = V.X.V to CHI 
!> @deprecated
!
!****************************************************************************************
SUBROUTINE STORE_VXV_TO_CHI( WDES, CHI, GDES, WGWQ, &
    COMM_BETWEENOMEGA, COMM_INOMEGA, LATT_CUR, &
    NOMEGA_START,NOMEGA_SIMULTANEOUS, COR, NQ , IDIR_MAX, NCUTOFF, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  !<response function
    TYPE (greensfdes)       :: GDES                 !> greensfunction descriptor
    TYPE (wavedes1)         :: WGWQ
    TYPE (communic)         :: COMM_BETWEENOMEGA    !<communicator between tau
    TYPE (communic)         :: COMM_INOMEGA         !<communicator in tau
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NOMEGA_START         !<number of frequency points
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !<# of freq. points on one core
    TYPE (correlation)      :: COR                  !<correlation energy
    INTEGER                 :: NQ                   !<current q-point
    INTEGER                 :: IDIR_MAX             !<maximum number of directions
    INTEGER                 :: NCUTOFF              !<cutoff index for stoarge
    TYPE (in_struct)        :: IO                 
  ! local
    GDEF, POINTER, DIMENSION(:,:) :: CHI_WORK
    INTEGER                 :: NROWS_NEW, NCOLS_NEW !new size of redistributed chi 
    INTEGER                 :: NBLOCK                 !block size
    INTEGER ::  IDIR, NP_NEW, I,IP, NK, NCURR
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    INTEGER                 :: NROWS, NCOLS         !<number of local rows and cols of CHI
    
    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("internal error in STORE_XVX_TO_CHI: \n or the Gamma point version " &
          // "CHI%LREALSTORE must be set", __FILE__, __LINE__)
    ENDIF

    IF (NCUTOFF > NE .OR. NCUTOFF <= 0 ) THEN
       CALL vtutor%error("cutoff index negative or too large! " // str(NE) // " " // str(NCUTOFF))
    ENDIF
    NROWS = GDES%RES_NRPLWV_ROW_DATA_POINTS
    NCOLS = GDES%RES_NRPLWV_COL_DATA_POINTS

    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    !scaLAPACK block size and processor grid is set here
 
    omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
    !--------------------------------------------------------------------------------------
       IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
          !this is the global omega point treated currently in each group 
          NOMEGA_GLOBAL = NOMEGA_START + NCURR - 1 
          NOMEGA_CURRENT = NCURR
          LDMODE = .FALSE.
       ELSE
          !group having less points than the root group use their largest local frequency
          NOMEGA_GLOBAL = NOMEGA_START + NOMEGA_SIMULTANEOUS - 1 
          NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
          LDMODE = .TRUE. 
       ENDIF

       !-----------------------------------------------------------------------------------
       ! directions only need to be considered for the head 
       IF ( CHI%LGAMMA ) THEN
       directions: DO IDIR = 1, IDIR_MAX          !loop over the tree directions  
       !-----------------------------------------------------------------------------------
          !for the long-wave limit we copy head and wings to RESPONSEFUN
          CALL BODY_FROM_WING_GG( CHI, IDIR, NOMEGA_CURRENT, GDES, COMM_INOMEGA)

          !CHI is distributed in groups to all processors
          !each group contains CHI at one frequency point. 
          !multiplication with the Coulomb potential needs to be done carefully,
          !CHI will be used for GW calculations, therefore
          !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
          !using CHI_WORK 

          !first we allocate CHI_WORK in the same way as CHI
          ALLOCATE(CHI_WORK(NROWS,NCOLS)) 

!test: compute V.X.V and store the result to CHI
          CHI_WORK = 0
          CALL DETERMINE_VXV( CHI_WORK, CHI, NROWS, NCOLS,                 &
             COMM_INOMEGA, LATT_CUR, COR%ENCUTGW(NCUTOFF), &
             COR%ENCUTGWSOFT(NCUTOFF),  &
             NOMEGA_CURRENT, WGWQ , NP_NEW,     &
             .FALSE. )
!test 
          DEALLOCATE ( CHI_WORK)     
       !-----------------------------------------------------------------------------------
       ENDDO directions  
       ELSE 
       !-----------------------------------------------------------------------------------
          !CHI is distributed in groups to all processors
          !each group contains CHI at one frequency point. 
          !multiplication with the Coulomb potential needs to be done carefully,
          !CHI will be used for GW calculations, therefore
          !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
          !using CHI_WORK 
          
          !first we allocate CHI_WORK in the same way as CHI
          ALLOCATE(CHI_WORK(NROWS,NCOLS)) 
!test: compute V.X.V
          CHI_WORK = 0
          CALL DETERMINE_VXV( CHI_WORK, CHI, NROWS, NCOLS,                 &
             COMM_INOMEGA, LATT_CUR, COR%ENCUTGW(NCUTOFF), &
             COR%ENCUTGWSOFT(NCUTOFF),  &
             NOMEGA_CURRENT, WGWQ , NP_NEW,   &
             .FALSE. )
!test
          DEALLOCATE ( CHI_WORK)     
       !-----------------------------------------------------------------------------------
       ENDIF
    !--------------------------------------------------------------------------------------
    ENDDO omega_sim

 RETURN

   CONTAINS

!************************************ DETERMINE_VXV ************************************
!
!> this subroutine truncates the response function by multiplying from left and right
!> with a truncated Coulomb kernel smoothly going from 1/G^2 to zero
!> between ENCUTGWSOFT and ENCUTGW
!
!****************************************************************************************

  SUBROUTINE DETERMINE_VXV( CHI_WORK, CHI, NROWS, NCOLS, &
    COMM_INOMEGA, LATT_CUR, ENCUT, ENCUTSOFT, &
    NCURR, WGWQ , MAXINDEX, &
    LSR)
    USE constant
    USE fock
    USE mpimy
    USE c2f_interface, ONLY : ERRF
    IMPLICIT NONE
    GDEF,POINTER  :: CHI_WORK(:,:)   !< response function with proper cutoff 
    TYPE (responsefunction) :: CHI   !< unchanged response function
    TYPE (latt) LATT_CUR             !< lattice structure
    INTEGER :: NCOLS, NROWS
    TYPE(communic) :: COMM_INOMEGA   !< communicator in one group
    REAL(q) :: ENCUT                 !< energy cutoff for response function
    REAL(q) :: ENCUTSOFT             !< lower cutoff for response function 
    INTEGER :: NCURR                 !< current freuquency point
    TYPE (wavedes1) :: WGWQ          !< response function descriptor for current q-point
    INTEGER :: MAXINDEX              !< number of grid points using truncated Coulomb kernel
    LOGICAL :: LSR
  ! local
    INTEGER    NI, NP, NI_, N1, I, J, I_
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE, POTFAK, E
    INTEGER :: NI_LOC, NOLD_LOC, NODE_RECV, IERROR
    REAL(q), ALLOCATABLE :: DATAKE(:)
    INTEGER, ALLOCATABLE :: OLD_INDEX(:)
    INTEGER :: ISEND(4) , IRECV(4)
    GDEF, ALLOCATABLE    :: CHI_RECV(:)
    GDEF, ALLOCATABLE    :: CHI_SEND(:)
    REAL(q) :: DFUN, SFUN
    ! neu
    REAL(q) :: OMEGBK,QC
    REAL(q) :: GAMMASCALE

    ! multipole correction
    INTEGER    NJ,NK
    COMPLEX(q) :: TSUM(10)
    COMPLEX(q) :: POTCORRECTION(10)
    !for communication
    INTEGER    :: NI_GLOB, NSEND
    INTEGER    :: NODE_SEND

    ! e^2/ volume 
    SCALE=EDEPS/LATT_CUR%OMEGA
!=======================================================================
! first set up the truncated Coulomb kernel
! smoothly going from 1/G^2 to zero between ENCUTGWSOFT and ENCUTGW
!=======================================================================
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

    ALLOCATE( DATAKE(NROWS), OLD_INDEX(NROWS))
    DO NI = 1, NROWS
     OLD_INDEX(NI) = NI
    ENDDO

    !treated frequency point
    N1=NCURR

    DATAKEMAX=0.0
    MAXINDEX=0
    DO NI=1,NP
       NI_=NI
       IF (WGWQ%LGAMMA) NI_=(NI-1)/2+1
       
       GX=(WGWQ%IGX(NI_)*LATT_CUR%B(1,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WGWQ%IGX(NI_)*LATT_CUR%B(2,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WGWQ%IGX(NI_)*LATT_CUR%B(3,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
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
             POTFAK=POTFAK*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
          ENDIF
       ENDIF

! POTFAK in fock has an additional factor. Since the corrections were originally devised
! for the fock routine this ensures consistency
       IF(MCALPHA/=0) THEN
          POTFAK=POTFAK*(1.0_Q/WGWQ%GRID%NPLWV)
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
       CALL GW_MULTIPOLE_PROJ_SETUP( WGWQ, LATT_CUR, DATAKE,DATAKEMAX )
    ENDIF

!=======================================================================
! now multiply the response function with the kernel
! from left and right hand side
!=======================================================================

    !save chi to chi_work
    IF (CHI%LREALSTORE) THEN
          CHI_WORK(1:NROWS,1:NCOLS) = CHI%RESPONSER(   1 : NROWS , 1 : NCOLS , N1 )
    ELSE
          CHI_WORK(1:NROWS,1:NCOLS) = CHI%RESPONSEFUN( 1 : NROWS , 1 : NCOLS , N1 )
    ENDIF

    !kill elements larger than NP
    IF ( NP < NROWS ) THEN
       DO NI = NP + 1, NROWS
          DATAKE(NI) = 0
       ENDDO
    ENDIF

    !DATAKE(:) has the coulomb potenial, this needs to be multiplied with CHI_WORK
    !and symmetrized
    !from the left
    DO NI=1,NROWS  
       !DATAKE contains the square root of V!
       POTFAK=DATAKE( NI )*DATAKE( NI )
       CHI_WORK(NI,1:NCOLS)=CHI_WORK(NI,1:NCOLS)*POTFAK
    ENDDO
    !from right
    DO NI=1,NCOLS
       !DATAKE contains the square root of V!
       POTFAK=DATAKE( NI + (COMM_INOMEGA%NODE_ME-1)*NCOLS )**2
       CHI_WORK(1:NROWS,NI)=CHI_WORK(1:NROWS,NI)*POTFAK
    ENDDO


    ! multipole corrections:
    IF(MCALPHA/=0) THEN 
      CALL vtutor%error("Sorry MCALPHA not implemented")
    ENDIF
   
    !smooth energy cutoff done by replacing the NI^th coulmn by the OLD_INDEX(NI)^th column 
    !and analogously for the row. 

    !The replacement of the rows is straight forward 
    DO NI=1,MAXINDEX
       CHI_WORK(NI,1:NCOLS)=CHI_WORK(OLD_INDEX(NI),1:NCOLS)
    ENDDO

    !for the columns this is not straight forward
    !since the whole matrix is distributed in stripes of NCOLS
    !among the nodes in each frequency group.

    !send columns between nodes
    ALLOCATE(CHI_SEND( NROWS ) )  
    ALLOCATE(CHI_RECV( NROWS ) ) 
    CHI_RECV = 0 
    CHI_SEND = 0 

    CALL MPI_BARRIER(COMM_INOMEGA%MPI_COMM,IERROR)   !each node should start at the same time from here
    !==============================================================================================
    !we loop over the cpus in one group and let other node send the appropriate columns
    DO NODE_RECV = 1, COMM_INOMEGA%NCPU
    !==============================================================================================
        loc_col : DO NI = 1 , NCOLS                      !search column which should be replaced
        !------------------------------------------------------------------------------------------
 
           NI_GLOB = NI + (NODE_RECV-1)*NCOLS   !global index for current node NODE_RECV
  
 
           !find node NODE_SEND containing OLD_INDEX(NI_GLOB)
           DO NODE_SEND = 1, COMM_INOMEGA%NCPU
             IF ( 1 + ( NODE_SEND-1)*NCOLS <= OLD_INDEX(NI_GLOB) .AND. &
                  OLD_INDEX(NI_GLOB) <= NODE_SEND*NCOLS) EXIT
           ENDDO
 
           NSEND = OLD_INDEX(NI_GLOB) - (NODE_SEND-1)*NCOLS   !local column index on sending node
 
           IF ( NSEND < 1 ) THEN 
              CALL vtutor%error("Cannot send column with negative index " &
                 &// str(NODE_RECV) // " " // str(NODE_SEND) // " " // str(NI_GLOB) // " " // &
                 str(OLD_INDEX(NI_GLOB)) // " " // str(NSEND) // " " // str(NCOLS))
           ENDIF
 
           
           !if old column is on the same node simply replace it 
           IF ( NODE_SEND == COMM_INOMEGA%NODE_ME .AND. NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN 
#if debug  
              WRITE(100+WGWQ%COMM%NODE_ME,&
                  '( " COL ON SAME CPU #",I4," COL #",2I4, "REPLACED BY COL #", I4,10F10.5)')&
                  NODE_RECV, NI,NI_GLOB, NSEND,CHI_WORK(1:10,NSEND)
#endif
              CHI_WORK(1:NROWS, NI) = CHI_WORK(1:NROWS, NSEND)

           !otherwise tell NODE_SEND to send the appropriate column NSEND to NODE_RECV
           ELSEIF ( NODE_RECV < NODE_SEND ) THEN
              IF ( NODE_SEND == COMM_INOMEGA%NODE_ME ) THEN
#if debug  
                 WRITE(100+WGWQ%COMM%NODE_ME,'( " COL #",I4, "SENT TO COL", I4 , " ON CPU #", I4)')&
                 NSEND,NI,NODE_RECV
#endif
                 CHI_SEND(1:NROWS) = CHI_WORK(1:NROWS, NSEND)

                 !send this column to node NODE_RECV              
                 IF (CHI%LREALSTORE) THEN
                    CALLMPI ( M_send_d( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                 ELSE
                    CALLMPI ( M_send_z( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                 ENDIF

              ELSEIF ( NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN 
#if debug  
                 WRITE(100+WGWQ%COMM%NODE_ME,'( " RECV COL #", I4 ," FROM CPU #"'//&
                    ',I4, " FOR COL", I4)') NSEND, NODE_SEND, NI
#endif
                 IF (CHI%LREALSTORE) THEN
                    CALLMPI ( M_recv_d( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                 ELSE
                    CALLMPI ( M_recv_z( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                 ENDIF
                     
                 !and replace CHI_RECV 
                 CHI_WORK(1:NROWS, NI ) = CHI_RECV(1:NROWS)
              ENDIF
           ENDIF

        !mpi barrier here 
        CALL MPI_BARRIER(COMM_INOMEGA%MPI_COMM,IERROR)   
        !------------------------------------------------------------------------------------------
        ENDDO loc_col                                !local columns  
    !==============================================================================================
    ENDDO                                    !cpus in group
    !==============================================================================================
    DEALLOCATE( CHI_SEND )
    DEALLOCATE( CHI_RECV )

    DEALLOCATE(DATAKE, OLD_INDEX)

    !save the result V.X.V to CHI
    !save chi to chi_work
    IF (CHI%LREALSTORE) THEN
       CHI%RESPONSER(   1 : NROWS , 1 : NCOLS , N1 ) = CHI_WORK(1:NROWS,1:NCOLS)
    ELSE
       CHI%RESPONSEFUN( 1 : NROWS , 1 : NCOLS , N1 ) = CHI_WORK(1:NROWS,1:NCOLS) 
    ENDIF

 END SUBROUTINE DETERMINE_VXV


END SUBROUTINE STORE_VXV_TO_CHI
#endif 
#endif

END MODULE acfdt_gg
