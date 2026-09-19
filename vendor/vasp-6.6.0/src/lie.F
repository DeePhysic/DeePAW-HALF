!#define debug
#include "symbol.inc"

!*********************************************************************
!
!  This module contains routines for the constrained optimization of a 
!  cost function on the Lie manifold of all unitary matrices U(n).
!  The main idea is to exploit the geometry of U(n), which is quite
!  well known. For instance the spin manifold SU(2) is isomorph to 
!  SO(3).
!  The main algorithm FIND_ON_LIE_MANIFOLD is from Traian Abrudan,
!  who proposed the algorithm in 
!
!  [1] Signal Processing, vol. 89, 2009, 1704-1714.
!  [2] Transactions on Signal Processing, vol. 56, 2008, 1134-1147.
!
!  The routines COST_FUNCTION and EUCLIDIAN_GRADIENT are such that
!  off-diagonal terms (defined by NTARGET) are minimized. This has 
!  the advantage, that we can find arbitrary Jordan decompositions
!  of S (useful for CRPA).
!  By changing these functions, one can minimize arbitrary real valued
!  scalar cost functions on U(n).
!
!                                                 Merzuk Kaltak (2014)
!
!*********************************************************************

MODULE lie
   USE prec
   USE base
   USE constant
   IMPLICIT NONE
   
   PRIVATE 
   PUBLIC :: FIND_ON_LIE_MANIFOLD,JOINT_SEPARATION_IN_LIEGROUP
    
   !global variables 
   COMPLEX(q), SAVE, POINTER :: S(:,:) => NULL()     !matrix to be optimized
   COMPLEX(q), SAVE, POINTER :: SK(:,:,:) => NULL()  !matrices to be optimized
   REAL(q), SAVE             :: CRIT_THRESH          !critical threshold
   INTEGER, SAVE             :: NDIM                 !NDIM**2=dimension of matrix
   INTEGER, SAVE             :: QDIM=4               !order of matrix in cost function
   INTEGER, SAVE, POINTER    :: NTARGET(:) => NULL() !target matrix elements 
   ! some tuning parameter   
   INTEGER, SAVE             :: MITERCG=100          !max number of CG iterations
   ! some tiny number  
   REAL(q), SAVE             :: TINY_U=1.E-5_q       !threshold for unitarity
   LOGICAL, SAVE             :: LMAXTR=.FALSE.
   INTEGER, SAVE             :: NMAT=-1
CONTAINS

!*******************JOINT_SEPARATION_IN_LIEGROUP**********************
!
! A set of hermitian matrices is used to build the cost function
!
!*********************************************************************

   SUBROUTINE JOINT_SEPARATION_IN_LIEGROUP(SOPT,WOPT,N,NKTOT,NT,CT,MITER,IO)
      IMPLICIT NONE
      COMPLEX(q)              :: SOPT(:,:,:) !input matrices to be optimized
      COMPLEX(q)              :: WOPT(:,:) !input matrices to be optimized
      INTEGER                 :: N         !dimension
      INTEGER                 :: NKTOT     !# of k-points
      INTEGER                 :: NT(:)     !sets the target matrix elements  
      REAL(q)                 :: CT        !convergence threshold
      INTEGER                 :: MITER     !maximum number of iterations
      TYPE(in_struct)         :: IO        !IO handle for writing to CAR files 
      !local variables
      CHARACTER(LEN=3)        :: OPT       !minimum or maximum is decided by OPT
      INTEGER                 :: I,J,K    !some loop variables
      COMPLEX(q), ALLOCATABLE :: W0(:,:)  !starting guess for solution of problem             
      COMPLEX(q), ALLOCATABLE :: U0(:,:)  !starting guess for solution of problem             
      REAL(q)                 :: JSTART   !starting cost function value 
      REAL(q)                 :: JCURRENT !current cost function value 
      INTEGER                 :: ITER     
      INTEGER                 :: CGSTEPS    
      INTEGER                 :: SDSTEPS    
      LOGICAL                 :: LSUCC    !information if success
      INTEGER                 :: ITRY 
      !for diagonal attempts
      INTEGER,PARAMETER       :: NDIAG_TRIES=50 
      REAL(q),ALLOCATABLE     :: R1(:,:)
      REAL(q),ALLOCATABLE     :: R2(:,:)

      CGSTEPS=MITER
      SDSTEPS=MITER
      OPT='MIN'
      LSUCC=.FALSE.
   
      !--------------------------------------------------------------------------------------
      ! first of all, set global variables 
      ! dimension
      IF ( N<1 ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,I3)')'ERROR in JOINT_SEPARATION_IN_LIEGROUP:'//&
            ' N<1',N
         STOP
      ELSE
         NDIM=N
      ENDIF
      ! # of matrices to be optimized
      IF ( NKTOT<1 ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,I3)')'ERROR in JOINT_SEPARATION_IN_LIEGROUP:'//&
            'NTOT<1',NKTOT
         STOP
      ELSE
         NMAT=NKTOT
         !NMAT=2
      ENDIF
      ! target matrix elements  
      IF ( SIZE(NT)<NDIM ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,2I4)')'ERROR in JOINT_SEPARATION_IN_LIEGROUP: '//&
            'NT does not contain enough entries',SIZE(NT),NDIM
      ELSE
         IF ( ASSOCIATED(NTARGET) ) NULLIFY(NTARGET)
         ALLOCATE(NTARGET(NDIM))
         NTARGET(1:NDIM)=NT(1:NDIM) 
      ENDIF 
      ! convergence precision 
      IF ( CT<0 ) THEN
         IF(IO%IU0>=0) WRITE(IO%IU0,'(A,E10.3)')'ERROR in JOINT_SEPARATION_IN_LIEGROUP:'//&
            ' CT<0',CT
         STOP
      ELSE
         IF ( CT < 1.E-12_q ) THEN
         IF(IO%IU0>=0) WRITE(IO%IU0,'(A,E10.3)')'ERROR in JOINT_SEPARATION_IN_LIEGROUP:'//&
               ' CT maybe too small',CT
         ENDIF 
         CRIT_THRESH=CT
      ENDIF 
      ! global variables set
      !--------------------------------------------------------------------------------------

      !allocate matrix, which should be extremized w.r.t. cost function
      IF ( ASSOCIATED(SK)) NULLIFY(SK)
      ALLOCATE(SK(NDIM,NDIM,NMAT)) ; SK=0_q
      DO K = 1, NMAT
         DO I = 1, NDIM
            DO J = 1, NDIM
               SK(I,J,K)=SOPT(I,J,K) 
            ENDDO
         ENDDO
      ENDDO
#ifdef debug
WRITE(*,*)'one of the matrices to be extremized'
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')REAL(SK(I,1:NDIM,NMAT),q)
ENDDO
WRITE(*,*)'  ' 
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')IMAG(SK(I,1:NDIM,NMAT))
ENDDO
#endif
      !allocate starting guess matrix W0 
      ALLOCATE( W0(NDIM, NDIM)) 

      !intitialize W0 as the identity matrix (for the begining)
      W0=zero
      DO I = 1, NDIM
         W0(I,I)=one
      ENDDO 
      !evaluate cost function with W0
      JSTART=COST_FUNCTION(W0,IO%IU0)
      IF (IO%IU0>=0) WRITE(IO%IU0,'(A,E12.4)')' Cost function value to be minimized:',JSTART
      !optimize matrix W0 using conjugated gradient method first and SD second 
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO,LSUCC)     !CG PR A
      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO,LSUCC)     !Poly CG
      IF ( .NOT. LSUCC ) THEN
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)     !FFT SD
      ENDIF 
      !if failed start from begining with random diagonal unitary matrix 
      !(diagonal, because we still want to start from the Lie Algbra
      IF ( .NOT. LSUCC  ) THEN 
         !we will try this NDIAG_TRIES times 
         diag: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            !evaluate cost function with W0
            JSTART=COST_FUNCTION(W0,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO,LSUCC) !Poly CG
            IF ( .NOT. LSUCC ) THEN
               CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)!FFTSD
            ENDIF 
            IF ( LSUCC ) THEN
               EXIT diag
            ENDIF 
         ENDDO diag
      ENDIF 
 
      !if failed start from a random point on Lie group,
      !not necessarily an element of the Lie algebra
      IF ( .NOT. LSUCC  ) THEN 
         !we will try this NDIAG_TRIES times 
         gen: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT(NDIM,IO%IU0)
            !evaluate cost function with W0
            JSTART=COST_FUNCTION(W0,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO,LSUCC) !Poly CG
            IF ( .NOT. LSUCC ) THEN
               CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)!FFTSD
            ENDIF 
            IF ( LSUCC ) THEN
               EXIT gen
            ENDIF 
         ENDDO gen
      ENDIF 

      !maybe the solution is still not found yet. If so, 
      !we try starting with conjugate gradient FFT method 
      !but start with the zero vector on the Lie algebra
      IF ( .NOT. LSUCC ) THEN
         W0=zero
         DO I = 1, NDIM
            W0(I,I)=one
         ENDDO 
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,CGSTEPS,IO)        !FFT CG
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,SDSTEPS,IO,LSUCC)  !POLY SD
         !if still not found try random diagonal matrices         
         diag_rev: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,CGSTEPS,IO)       !FFT CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,SDSTEPS,IO,LSUCC) !POLY SD
            IF ( LSUCC ) THEN
               EXIT diag_rev
            ENDIF 
         ENDDO diag_rev
      ENDIF   

      !maybe the solution is still not found yet. If so, 
      !we try the Armijo method for the conjugate gradient
      !but start with the zero vector on the Lie algebra
      IF ( .NOT. LSUCC ) THEN
         W0=zero
         DO I = 1, NDIM
            W0(I,I)=one
         ENDDO 
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO)        !ARM CG
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)  !FFT SD
         !if still not found try random diagonal matrices         
         diag_armijo: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO)        !ARM CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)  !FFT SD
            IF ( LSUCC ) THEN
               EXIT diag_armijo
            ENDIF 
         ENDDO diag_armijo
      ENDIF   

      !if still no sucess, give up here 
      IF ( .NOT. LSUCC ) THEN
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A,I5)')'ERROR in JOINT_SEPARATION_IN_LIEGROUP:'//&
            ' Jordan decomposition was NOT found'
         STOP
      ELSE
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')'Jordan decomposition found'
         WOPT=zero
         DO I = 1, NDIM
         DO J = 1, NDIM
            WOPT(I,J)=W0(I,J)
         ENDDO
         ENDDO
!#ifdef debug
WRITE(*,*)'one of the matrices to be extremized'
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')REAL(SK(I,1:NDIM,NMAT),q)
ENDDO
WRITE(*,*)'  ' 
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')IMAG(SK(I,1:NDIM,NMAT))
ENDDO
WRITE(*,*)'Jordan decomposition is'
        ALLOCATE(U0(NDIM,NDIM))
        U0=MATMUL(CONJG(TRANSPOSE(W0)),MATMUL(SK(:,:,NMAT),W0))
DO I = 1,NDIM
   WRITE(*,'(10F14.6)')REAL(U0(I,1:NDIM),q)
ENDDO
WRITE(*,*)'  ' 
DO I = 1,NDIM
   WRITE(*,'(10F14.6)')IMAG(U0(I,1:NDIM))
ENDDO
!#endif
      ENDIF 

      DEALLOCATE(W0)
      IF (ASSOCIATED(SK)) NULLIFY(SK)
      IF (ASSOCIATED(NTARGET)) NULLIFY(NTARGET)
      RETURN
   END SUBROUTINE JOINT_SEPARATION_IN_LIEGROUP

!*************************FIND_ON_LIE_MANIFOLD************************
!
! This is the main routine of this module. It minimizes the 
! cost function COST_FUNCTION(S) on the Lie manifold costructed by
! the set of all unitary matrices U. 
! S     ... matrix to be optimized 
! N     ... NxN is the dimension of the Lie manifold/matrix
! OPT   ... either 'min' or 'max' of COST_FUNCTION(S) is found
! NT    ... target states of S to be extremized
! CT    ... convergence threshold 
! MITER ... maximum number of iterations
!
!*********************************************************************

   SUBROUTINE FIND_ON_LIE_MANIFOLD(SOPT,N,OPT,NT,CT,MITER,NK,IO,LMAXIMIZE,LTEST)
      IMPLICIT NONE
      COMPLEX(q)              :: SOPT(:,:) !input matrix to be optimized
      INTEGER                 :: N         !dimension
      CHARACTER(LEN=3)        :: OPT       !minimum or maximum is decided by OPT
      INTEGER                 :: NT(:)     !sets the target matrix elements  
      REAL(q)                 :: CT        !convergence threshold
      INTEGER                 :: MITER     !maximum number of iterations
      INTEGER                 :: NK        !k-point 
      TYPE(in_struct)         :: IO        !IO handle for writing to CAR files 
      LOGICAL                 :: LMAXIMIZE !switch to other cost function
      LOGICAL, OPTIONAL       :: LTEST     !test switch
      !local variables
      INTEGER                 :: I,J,K    !some loop variables
      COMPLEX(q), ALLOCATABLE :: W0(:,:)  !starting guess for solution of problem             
      REAL(q)                 :: JSTART   !starting cost function value 
      REAL(q)                 :: JCURRENT !current cost function value 
      INTEGER                 :: ITER     
      INTEGER                 :: CGSTEPS    
      INTEGER                 :: SDSTEPS    
      LOGICAL                 :: LSUCC    !information if success
      INTEGER                 :: ITRY 
      !for diagonal attempts
      INTEGER,PARAMETER       :: NDIAG_TRIES=50 

      REAL(q),ALLOCATABLE     :: R1(:,:)
      REAL(q),ALLOCATABLE     :: R2(:,:)
      COMPLEX(q), ALLOCATABLE :: TMP(:,:)  !starting guess for solution of problem             

      CGSTEPS=MITER
      SDSTEPS=MITER
   
      !--------------------------------------------------------------------------------------
      ! first of all, set global variables 
      ! dimension
      IF ( N<1 ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,I3)')'ERROR in FIND_ON_LIE_MANIFOLD: N<1',N
         STOP
      ELSE
         NDIM=N
      ENDIF
      ! target matrix elements  
      IF ( SIZE(NT)<NDIM ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,2I4)')'ERROR in FIND_ON_LIE_MANIFOLD: '//&
            'NT does not contain enough entries',SIZE(NT),NDIM

      ELSE
         IF ( ASSOCIATED(NTARGET) ) NULLIFY(NTARGET)
         ALLOCATE(NTARGET(NDIM))
         NTARGET(1:NDIM)=NT(1:NDIM) 
      ENDIF 
      !overwrite this if maximization of Tr is needed
      IF ( LMAXIMIZE ) THEN
         NTARGET(1:NDIM)=1
         LMAXTR=.TRUE.
      ELSE
         LMAXTR=.FALSE. ! just to be sure it's not set
      ENDIF 

      ! convergence precision 
      IF ( CT<0 ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,E10.3)')'ERROR in FIND_ON_LIE_MANIFOLD: CT<0',CT
         STOP
      ELSE
         IF ( CT < 1.E-12_q ) THEN
            IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,E10.3)')'WARNING from FIND_ON_LIE_MANIFOLD:'//&
               ' CT maybe too small',CT
         ENDIF 
         CRIT_THRESH=CT
      ENDIF 
      ! global variables set
      !--------------------------------------------------------------------------------------

      !allocate matrix, which should be extremized w.r.t. cost function
      IF ( ASSOCIATED(S)) NULLIFY(S)
      ALLOCATE(S(NDIM,NDIM)) ; S=0_q
      DO I = 1, NDIM
         DO J = 1, NDIM
            S(I,J)=SOPT(I,J) 
         ENDDO
      ENDDO

      !use a random hermitian matrix
      IF (PRESENT(LTEST).AND.LTEST) THEN
          S=RHERM(NDIM,IO%IU0)
      ENDIF 

!#ifdef debug
WRITE(*,*)'Matrix to be extremized'
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')REAL(S(I,1:NDIM),q)
ENDDO
WRITE(*,*)'  ' 
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')IMAG(S(I,1:NDIM))
ENDDO
!#endif

      !allocate starting guess matrix W0 
      ALLOCATE( W0(NDIM, NDIM)) 

      !intitialize W0 as the identity matrix (for the begining)
      W0=zero
      DO I = 1, NDIM
         W0(I,I)=one
      ENDDO 
      !evaluate cost function with W0
      JSTART=COST_FUNCTION(W0,IO%IU0)
      !optimize matrix W0 using conjugated gradient method first and SD second 
      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO)        !Poly CG
      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)        !FFT SD
      !if failed start from begining with random diagonal unitary matrix 
      !(diagonal, because we still want to start from the Lie Algbra
      IF ( .NOT. LSUCC  ) THEN 
         !we will try this NDIAG_TRIES times 
         diag: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            !evaluate cost function with W0
            JSTART=COST_FUNCTION(W0,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO)       !Poly CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC) !FFT SD
            IF ( LSUCC ) THEN
               EXIT diag
            ENDIF 
         ENDDO diag
      ENDIF 
 
      !if failed start from a random point on Lie group,
      !not necessarily an element of the Lie algebra
      IF ( .NOT. LSUCC  ) THEN 
         !we will try this NDIAG_TRIES times 
         gen: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT(NDIM,IO%IU0)
            !evaluate cost function with W0
            JSTART=COST_FUNCTION(W0,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO)       !Poly CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC) !FFT SD
            IF ( LSUCC ) THEN
               EXIT gen
            ENDIF 
         ENDDO gen
      ENDIF 

      !maybe the solution is still not found yet. If so, 
      !we try starting with conjugate gradient FFT method 
      !but start with the zero vector on the Lie algebra
      IF ( .NOT. LSUCC ) THEN
         W0=zero
         DO I = 1, NDIM
            W0(I,I)=one
         ENDDO 
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,CGSTEPS,IO)        !FFT CG
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,SDSTEPS,IO,LSUCC)  !POLY SD
         !if still not found try random diagonal matrices         
         diag_rev: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,CGSTEPS,IO)        !FFT CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,SDSTEPS,IO,LSUCC)  !POLY SD
            IF ( LSUCC ) THEN
               EXIT diag_rev
            ENDIF 
         ENDDO diag_rev
      ENDIF   

      !maybe the solution is still not found yet. If so, 
      !we try the Armijo method for the conjugate gradient
      !but start with the zero vector on the Lie algebra
      IF ( .NOT. LSUCC ) THEN
         W0=zero
         DO I = 1, NDIM
            W0(I,I)=one
         ENDDO 
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO)        !ARM CG
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)  !FFT SD
         !if still not found try random diagonal matrices         
         diag_armijo: DO ITRY = 1, NDIAG_TRIES
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO)        !ARM CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)  !FFT SD
            IF ( LSUCC ) THEN
               EXIT diag_armijo
            ENDIF 
         ENDDO diag_armijo
      ENDIF   

      !if still no sucess, give up here 
      IF ( .NOT. LSUCC ) THEN
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A,I5)')'ERROR in FIND_ON_LIE_MANIFOLD:'//&
            ' no Jordan decomposition was found for NK=',NK
         STOP
      ELSE
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A,I5)')'Jordan decomposition found for NK=',NK
!          SOPT=W0
         !otherwise overwrith Hamiltonian with decomposition 
         ALLOCATE(TMP(NDIM,NDIM))
         TMP=zero
         !TODO: use LAPACK's ZGEMM at this stage 
         TMP=MATMUL(MATMUL(CONJG(TRANSPOSE(W0)),S),W0)
         IF ( .NOT.PRESENT(LTEST)) THEN
            SOPT=TMP
         ENDIF 
!#ifdef debug
WRITE(*,*)'Extremized matrix'
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')REAL(TMP(I,1:NDIM),q)
ENDDO
WRITE(*,*)'  ' 
DO I = 1,NDIM
   WRITE(*,'(16F14.6)')IMAG(TMP(I,1:NDIM))
ENDDO
!#endif

         DEALLOCATE(TMP)
      ENDIF 

      DEALLOCATE(W0)
      IF (ASSOCIATED(S)) NULLIFY(S)
      IF (ASSOCIATED(NTARGET)) NULLIFY(NTARGET)
      RETURN
   END SUBROUTINE FIND_ON_LIE_MANIFOLD

!testing routine 
   SUBROUTINE FIND_ON_LIE_MANIFOLD_(SOPT,N,OPT,NT,CT,MITER,NK,IO)
      IMPLICIT NONE
      COMPLEX(q)              :: SOPT(:,:)!input matrix to be optimized
      INTEGER                 :: N        !dimension
      CHARACTER(LEN=3)        :: OPT      !minimum or maximum is decided by OPT
      INTEGER                 :: NT(:)    !sets the target matrix elements  
      REAL(q)                 :: CT       !convergence threshold
      INTEGER                 :: MITER    !maximum number of iterations
      INTEGER                 :: NK       !k-point 
      TYPE(in_struct)         :: IO       !IO handle for writing to CAR files 
      !local variables
      INTEGER                 :: I,J,K    !some loop variables
      COMPLEX(q), ALLOCATABLE :: W0(:,:)  !starting guess for solution of problem             
      REAL(q)                 :: JSTART   !starting cost function value 
      REAL(q)                 :: JCURRENT !current cost function value 
      INTEGER                 :: ITER     
      INTEGER,PARAMETER       :: CGSTEPS=75
      INTEGER,PARAMETER       :: SDSTEPS=100
      LOGICAL                 :: LSUCC    !information if success
      INTEGER                 :: ITRY 
      !for diagonal attempts
      INTEGER,PARAMETER       :: NDIAG_TRIES=50 

      REAL(q),ALLOCATABLE     :: R1(:,:)
      REAL(q),ALLOCATABLE     :: R2(:,:)
      COMPLEX(q), ALLOCATABLE :: TMP(:,:)  !starting guess for solution of problem             

   
      !--------------------------------------------------------------------------------------
      ! first of all, set global variables 
      ! dimension
      IF ( N<1 ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,I3)')'ERROR in FIND_ON_LIE_MANIFOLD: N<1',N
         STOP
      ELSE
         NDIM=N
      ENDIF
      ! target matrix elements  
      IF ( SIZE(NT)<NDIM ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,2I4)')'ERROR in FIND_ON_LIE_MANIFOLD: '//&
            'NT does not contain enough entries',SIZE(NT),NDIM

      ELSE
         IF ( ASSOCIATED(NTARGET) ) NULLIFY(NTARGET)
         ALLOCATE(NTARGET(NDIM))
         NTARGET(1:NDIM)=NT(1:NDIM) 
      ENDIF 
      ! convergence precision 
      IF ( CT<0 ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,E10.3)')'ERROR in FIND_ON_LIE_MANIFOLD: CT<0',CT
         STOP
      ELSE
         IF ( CT < 1.E-12_q ) THEN
            IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,E10.3)')'WARNING from FIND_ON_LIE_MANIFOLD:'//&
               ' CT maybe too small',CT
         ENDIF 
         CRIT_THRESH=CT
      ENDIF 
      ! global variables set
      !--------------------------------------------------------------------------------------

      !allocate matrix, which should be extremized w.r.t. cost function
      IF ( ASSOCIATED(S)) NULLIFY(S)
      ALLOCATE(S(NDIM,NDIM)) ; S=0_q
      DO I = 1, NDIM
         DO J = 1, NDIM
            S(I,J)=SOPT(I,J) 
         ENDDO
      ENDDO

test: DO J = 1, 200

WRITE(*,*)J
S=RHERM(NDIM,IO%IU0)

!#ifdef debug
WRITE(100+J,*)' real part of input matrix:' 
DO I = 1,NDIM
   WRITE(100+J,'(16F30.16)')REAL(S(I,1:NDIM))
ENDDO
WRITE(100+J,*)' imag part of input matrix:' 
DO I = 1,NDIM
   WRITE(100+J,'(16F30.16)')IMAG(S(I,1:NDIM))
ENDDO
!#endif
      !allocate starting guess matrix W0 
      ALLOCATE( W0(NDIM, NDIM)) 

      !intitialize W0 as the identity matrix (for the begining)
      W0=zero
      DO I = 1, NDIM
         W0(I,I)=one
      ENDDO 
      !evaluate cost function with W0
      JSTART=COST_FUNCTION(W0,IO%IU0)

      !optimize matrix W0 using conjugated gradient method first and SD second 
      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO)        !Poly CG
!WRITE(100+J,'(A,2E20.10)')'After CGPR P ',JSTART,JCURRENT
      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)        !FFT SD
!WRITE(100+J,'(I4,A,2E20.10,L4)')J,' SDSA F ',JSTART,JCURRENT,LSUCC

      !if failed start from begining with random diagonal unitary matrix 
      !(diagonal, because we still want to start from the Lie Algbra
      IF ( .NOT. LSUCC  ) THEN 
         !we will try this NDIAG_TRIES times 
         diag: DO ITRY = 1, NDIAG_TRIES
!WRITE(*,'(I4,A,I4)')J,' trying diagonal starting matrix ',ITRY
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            !evaluate cost function with W0
            JSTART=COST_FUNCTION(W0,IO%IU0)
            WRITE(100+J,'(A,I4)')'continuing with random diagonal unitary matrix ',ITRY
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,CGSTEPS,IO)       !Poly CG
!WRITE(100+J,'(A,2E20.10)')'After CGPR P ',JSTART,JCURRENT
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC) !FFT SD
WRITE(100+J,'(I4,A,2E20.10,L4)')J,' SDSA F ',JSTART,JCURRENT,LSUCC
           
            IF ( LSUCC ) THEN
               EXIT diag
            ENDIF 
         ENDDO diag
      ENDIF 

      !maybe the solution is still not found yet. If so, 
      !we try starting with conjugate gradient FFT method 
      !but start with the zero vector on the Lie algebra
      IF ( .NOT. LSUCC ) THEN
         WRITE(*,'(I4,A)')J,' trying reverse order '
         W0=zero
         DO I = 1, NDIM
            W0(I,I)=one
         ENDDO 
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,CGSTEPS,IO)        !FFT CG
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,SDSTEPS,IO,LSUCC)  !POLY SD

         !if still not found try random diagonal matrices         
         diag_rev: DO ITRY = 1, NDIAG_TRIES
            WRITE(*,'(I4,A,I4)')J,' trying diagonal starting matrix ',ITRY
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,CGSTEPS,IO)        !FFT CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,SDSTEPS,IO,LSUCC)  !POLY SD
            IF ( LSUCC ) THEN
               EXIT diag_rev
            ENDIF 
         ENDDO diag_rev
      ENDIF   

      !maybe the solution is still not found yet. If so, 
      !we try the Armijo method for the conjugate gradient
      !but start with the zero vector on the Lie algebra
      IF ( .NOT. LSUCC ) THEN
         WRITE(*,'(I4,A)')J,' Armijo method '
         W0=zero
         DO I = 1, NDIM
            W0(I,I)=one
         ENDDO 
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO)        !ARM CG
         CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)  !FFT SD

         !if still not found try random diagonal matrices         
         diag_armijo: DO ITRY = 1, NDIAG_TRIES
            WRITE(*,'(I4,A,I4)')J,' trying diagonal starting matrix ',ITRY
            W0=RUNIT_DIAG(NDIM,IO%IU0)
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,CGSTEPS,IO)        !ARM CG
            CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,SDSTEPS,IO,LSUCC)  !FFT SD
            IF ( LSUCC ) THEN
               EXIT diag_armijo
            ENDIF 
         ENDDO diag_armijo
      ENDIF   
  


IF (LSUCC) THEN 
WRITE(100+J,'(2E20.10,A)')JSTART,JCURRENT,'  SUCCESS!!!'
!print the result 
ALLOCATE(TMP(NDIM,NDIM))
 TMP=zero
 TMP=MATMUL(MATMUL(CONJG(TRANSPOSE(W0)),S),W0)
WRITE(100+J,*)' real part of solution:' 
DO I = 1,NDIM
   WRITE(100+J,'(16F15.8)')REAL(TMP(I,1:NDIM),q)
ENDDO
WRITE(100+J,*)' imag part of solution:' 
DO I = 1,NDIM
   WRITE(100+J,'(16F15.8)')IMAG(TMP(I,1:NDIM))
ENDDO
DEALLOCATE(TMP)
ELSE
WRITE(100+J,'(2E20.10)')JSTART,JCURRENT
EXIT test 
ENDIF 
!#endif 

DEALLOCATE(W0)
ENDDO  test

STOP

!      All possible calls are listed below 
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','P',OPT,10,IO)        !Poly CG
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','P',OPT,10,IO)        !Poly SD
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','A',OPT,10,IO)        !Armijo SD
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','A',OPT,10,IO)        !Armijo SD
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGFR','P',OPT,10,IO)        !Poly CGFR
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGFR','A',OPT,10,IO)        !Armijo CGFR
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'SDSA','F',OPT,10,IO)        !FFT SD
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGPR','F',OPT,10,IO)        !Poly CG
!      CALL EXTREMIZE_COST_FUNCTION(W0,JCURRENT,'CGFR','F',OPT,10,IO)        !Poly CG

      RETURN
   END SUBROUTINE FIND_ON_LIE_MANIFOLD_

!************************EXTREMIZE_COST_FUNCTION**********************
!
! main calculation routine:
! this routine finds the minimum on the U(n) Lie manifold, starting
! at the point W0.
!
!*********************************************************************

   SUBROUTINE EXTREMIZE_COST_FUNCTION(W0,JF,METHOD,ALGO,OPT,MITER,IO,LSUCC)
      IMPLICIT NONE
      COMPLEX(q)             :: W0(:,:)  !starting guess for solution
      REAL(q)                :: JF       !final cost function value
      CHARACTER(LEN=*)       :: METHOD   !which method is chosen CG or SD
      CHARACTER(LEN=*)       :: ALGO     !which line search algorithm are we using 
      CHARACTER(LEN=*)       :: OPT      !either 'min' or 'max'
      INTEGER                :: MITER    !maximum number of iterations 
      TYPE(in_struct)        :: IO       !IO handle for writing to CAR files 
      LOGICAL,OPTIONAL       :: LSUCC    !information if success
      !local 
      INTEGER                :: N         !that's the matrix size 
      REAL(q)                :: R         !norm of matrix 
      INTEGER                :: I         !loop variable
      LOGICAL                :: LUNITARY  !check string 
      REAL(q)                :: SIGN_MU   !moving direction of gradient  
      REAL(q)                :: MU        !magnitude of direction gradient  
      COMPLEX(q),ALLOCATABLE :: WK(:,:)   !auxillary iteration matrix
      INTEGER                :: KITER     !iteration variable 
      REAL(q)                :: J         !cost function value 
      REAL(q)                :: JINI      !cost function value 
      COMPLEX(q),ALLOCATABLE :: DEK(:,:)  !Euclidian gradient of cost function  
      COMPLEX(q),ALLOCATABLE :: DK(:,:)   !Riemanian gradient of cost function  
      COMPLEX(q),ALLOCATABLE :: HK(:,:)   !line search gradient of cost function  
      COMPLEX(q),ALLOCATABLE :: DKP1(:,:) !auxillary gradients for iteration
      COMPLEX(q),ALLOCATABLE :: DEKP1(:,:)!auxillary gradients for iteration
      COMPLEX(q),ALLOCATABLE :: HKP1(:,:) !auxillary gradients for iteration
      COMPLEX(q),ALLOCATABLE :: WKP1(:,:) !auxillary gradients for iteration
      COMPLEX(q),ALLOCATABLE :: R1(:,:)   !matrix exponential
      COMPLEX(q),ALLOCATABLE :: AUX(:,:)  !auxilary matrix
      INTEGER                :: NPOLY     !order of polynomial-based line search method
      REAL(q)                :: GAMMA_K   !factor for direction
      INTEGER                :: COUNTER   !counter for displaying warnings 
      LOGICAL                :: LBREAK    !breaking switch
 
      !initialize information 
      IF(PRESENT(LSUCC)) THEN
         LSUCC=.FALSE.    
      ENDIF 

      ! a few consistency check for starting guess:
      ! get the size of the matrix
      N=SIZE(W0,1)
      IF ( N /= SIZE(W0,2) ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,2I3)')'ERROR in EXTREMIZE_COST_FUNCTION:'//&
            ' W0 is not a square matrix:',N,SIZE(W0,2)
         STOP
      ENDIF 
      !also the matrix should have the size NDIM
      IF ( N /= NDIM ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,2I3)')'ERROR in EXTREMIZE_COST_FUNCTION:'//&
            ' W0 does not have proper size :',N,NDIM
         STOP
      ENDIF 

      !is W0 a unitary matrix ?       
      LUNITARY=LUNITCHECK(W0,IO%IU0)
      IF ( .NOT.LUNITARY ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A,L3)')'ERROR in EXTREMIZE_COST_FUNCTION:'//&
            ' W0 is not unitary:',LUNITARY
         STOP
      ENDIF 

      !are we computing the minimum of the maximum of the cost function 
      IF (OPT=='min'.OR.OPT=='MIN') THEN
         SIGN_MU=-1._q  !negative gradient
      ELSEIF(OPT=='max'.OR.OPT=='MAX') THEN 
         SIGN_MU=1._q   !positve gradient
      ELSE
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A)')'ERROR in EXTREMIZE_COST_FUNCTION: '//&
            'OPT not set properly '//OPT//'.'
         STOP
      ENDIF 

      !we need a few matrices in the following  
      !allocate memory for Euclidian and Riemanian gradient of cost function
      ALLOCATE(DEK(NDIM,NDIM),DK(NDIM,NDIM))
      !and initialize them
      DEK=zero
      DK=zero
      !also a matrix for the search direction will be needed 
      ALLOCATE(HK(NDIM,NDIM))
      HK=zero
      !updated gradients and search direction
      ALLOCATE(DKP1(NDIM,NDIM),DEKP1(NDIM,NDIM),HKP1(NDIM,NDIM))
      DKP1=zero
      DEKP1=zero
      HKP1=zero
      !updated solution
      ALLOCATE(WKP1(NDIM,NDIM))
      WKP1=zero
      !updated solution
      ALLOCATE(R1(NDIM,NDIM))
      R1=zero
      !some auxilary matrix
      ALLOCATE(AUX(NDIM,NDIM))
      AUX=zero
      !allocate memory for itaration matrix 
      ALLOCATE(WK(NDIM,NDIM))

      !JINI is always the cost function of the identity
      WK=zero
      DO I = 1, NDIM
        WK(I,I)=one
      ENDDO
      JINI=COST_FUNCTION(WK,IO%IU0)

      !initialize WK and cost function value J
      WK(:,:)=W0(:,:) 
      !cost function
      J=0._q
#ifdef debug
      DO I=1,NDIM
         WRITE(*,'(16F20.10)')REAL(S(I,1:NDIM),q)
      ENDDO 
      DO I=1,NDIM
         WRITE(*,'(16F20.10)')IMAG(S(I,1:NDIM))
      ENDDO 
#endif 
        
      !-------------------------------------------------------------------------------
      search: DO KITER=1,MITER 
      !-------------------------------------------------------------------------------
          !determine cost function value for current matrix WK
          J=COST_FUNCTION(WK,IO%IU0)
          !determine euclidian gradient, Riemannian gradient and search direction
          !note: index value KITER=1 corresponds to index =0 in Ref. [1] of Abrudan.
          !we can reset the search direction periodically, thanks to the topology
          !of the vector bundle of the Lie manifold.
          IF ( MOD(KITER-1,NDIM**2) == 0 ) THEN
              !Euclidian gradient of J at WK is saved to DEK
              DEK=EUCLIDIAN_GRADIENT(WK,IO%IU0)   
              !Riemanian gradient of J at Wk is saved to DK
              DK=RIEMANNIAN_GRADIENT(DEK,WK,IO%IU0)
              !periodically reset the search direction to the gradient DK 
              HK=DK
          ELSE 
              !otherwise DK = old DK + 1
              DK=DKP1   
              !otherwise HK = old HK + 1
              HK=HKP1
          ENDIF 

          !decide which geodesic line search algorithm is used 
          IF (ALGO=='a'.OR.ALGO=='A') THEN
             !Armijo line search method
!             IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')'Armijo line search method chosen'
             !find step size using formula of Armijo (see Ref. [2])
             MU=ARMIJO_STEP_SIZE(WK,DK,HK,OPT,IO%IU0)
             AUX=SIGN_MU*MU*SKEW(HK,IO%IU0)  
             !the matrix exponential of a skew hermitian matrix gives a nice unitary matrix R1
             R1=MATRIX_EXPONENTIAL(AUX,IO%IU0)  
          ELSEIF(ALGO=='p'.OR.ALGO=='P') THEN 
             !polynomial approximation-based line search method
!             IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')' polynomial search method is chosen '
             NPOLY=5  !fifth order is hard coded (maybe change this in future applications)
             MU=POLY_STEP_SIZE(WK,HK,NPOLY,OPT,IO%IU0)  
             AUX=SIGN_MU*MU*SKEW(HK,IO%IU0)   
             !the matrix exponential of a skew hermitian matrix gives a nice unitary matrix R1 
             R1=MATRIX_EXPONENTIAL(AUX,IO%IU0)
          ELSEIF(ALGO=='f'.OR.ALGO=='F') THEN
             !discrete Fourier transformation approximation-based line search method
!             IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')' FFT search method is chosen '
             !we use 5 almost periods for the FFT-based line search 
             NPOLY=5  ! abuse NPOLY as number of considered FFT periods
             MU=FFT_STEP_SIZE(WK,HK,NPOLY,OPT,IO%IU0)
             AUX=SIGN_MU*MU*SKEW(HK,IO%IU0)   
             !the matrix exponential of a skew hermitian matrix gives a nice unitary matrix R1 
             R1=MATRIX_EXPONENTIAL(AUX,IO%IU0)
          ELSE
             IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')'ERROR in EXTREMIZE_COST_FUNCTION:'//&
                ' bad line search method chosen:'//ALGO//'.'
             STOP
          ENDIF 

          !update solution
          WKP1(:,:)=zero
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,&
             zero,WKP1(1,1),NDIM)

          !determine the new Euclidian gradient at WKP1
          DEKP1=EUCLIDIAN_GRADIENT(WKP1,IO%IU0)
          !determine the new Riemannian gradient at WKP1
          DKP1=RIEMANNIAN_GRADIENT(DEKP1,WKP1,IO%IU0)

          !now determine the factor gamma_k, this depends on the choice of METHOD
          IF(METHOD=='sdsa'.OR.METHOD=='SDSA') THEN
             GAMMA_K=0._q
          ELSEIF(METHOD=='cgpr'.OR.METHOD=='CGPR') THEN
             GAMMA_K=INNER_PRODUCT(DKP1-DK,DKP1,IO%IU0)/INNER_PRODUCT(DK,DK,IO%IU0)
          ELSEIF(METHOD=='cgfr'.OR.METHOD=='CGFR') THEN
             GAMMA_K=INNER_PRODUCT(DKP1,DKP1,IO%IU0)/INNER_PRODUCT(DK,DK,IO%IU0)
          ELSE
             IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')'ERROR in EXTREMIZE_COST_FUNCTION:'//&
                ' bad method chosen:'//METHOD//'.'
             STOP
          ENDIF 

          !update search direction HKP1
          HKP1=DKP1+GAMMA_K*HK
 
          !check search direction
          R=INNER_PRODUCT(DKP1,HKP1,IO%IU0)
          !if wrong -> reset to steepest descent  
          IF ( R < 0 ) THEN
             COUNTER=COUNTER+1
             IF (COUNTER < 2 ) THEN
                IF ( IO%IU0>=0) WRITE(IO%IU0,'(A)')'EXREMIZE_COST_FUNCTION:'//&
                   ' HK set to SD/SA direction ... surpressing further messages'
             ENDIF 
             HKP1=DKP1
          ENDIF 

          !breaking condition
          !depends if SK is associated we use the change in WK as criterion
          IF ( ASSOCIATED(SK) ) THEN
             J=COST_FUNCTION(WKP1,IO%IU0)
             GAMMA_K=INNER_PRODUCT(DK,DK,IO%IU0)
             R=ABS(GAMMA_K)
             !GAMMA_K=NORM(DKP1-DK,0,IO%IU0)
             !WK should not change
             IF (ABS(GAMMA_K) < CRIT_THRESH) THEN
                !we oppose a constraint on the cost function as well 
                IF (OPT=='min'.OR.OPT=='MIN') THEN
                   IF ( J < JINI ) THEN
                      LBREAK = .TRUE.
                   ELSE
                      LBREAK = .FALSE.
                   ENDIF
                ELSEIF(OPT=='max'.OR.OPT=='MAX') THEN 
                   IF ( J > JINI ) THEN
                      LBREAK = .TRUE.
                   ELSE
                      LBREAK = .FALSE.
                   ENDIF
                ENDIF
             ELSE
                LBREAK = .FALSE.
             ENDIF 
          ELSE
             R=ABS(COST_FUNCTION(WKP1,IO%IU0))
             IF ( R<CRIT_THRESH) THEN
                LBREAK = .TRUE.
             ELSE
                LBREAK = .FALSE.
             ENDIF 
          ENDIF

          !update WK
          WK=WKP1

          IF ( LBREAK ) THEN
              EXIT search
          ENDIF 

      !-------------------------------------------------------------------------------
      ENDDO search
      !-------------------------------------------------------------------------------
      
      !overwrite W0 and use R as final cost function value
      W0=WK
      JF=COST_FUNCTION(W0,IO%IU0)

      !now it's time to get rid of all allocated arrays
      DEALLOCATE(WK)
      DEALLOCATE(DEK,DK,HK)
      DEALLOCATE(DKP1,DEKP1,HKP1,WKP1)
      DEALLOCATE(R1,AUX)

      IF ( KITER > MITER ) THEN  !Failure
!          IF ( IO%IU0>=0) WRITE(IO%IU0,1000)METHOD//' '//ALGO,KITER,JF
          IF (PRESENT(LSUCC)) THEN
             LSUCC=.FALSE.
          ENDIF 
      ELSE                       !Success
          IF ( IO%IU0>=0) WRITE(IO%IU0,1001)METHOD//' '//ALGO,KITER,R,JF
          IF (PRESENT(LSUCC)) THEN
             LSUCC=.TRUE.
          ENDIF 
      ENDIF 

      1001 FORMAT(A," converged after",I4," iterations with Error ",E12.4," and J=",E12.4)
      1000 FORMAT("WARNING: ",A," NOT converged after",I4," iterations, CF=",E12.4)
    END SUBROUTINE EXTREMIZE_COST_FUNCTION

!*********************************************************************
!
!                           F U N C T I O N S
!
!*********************************************************************
 
!***********************COST_FUNCTION*********************************
!
! this is the cost function, which is minimized on the Lie manifold
!
!*********************************************************************

   FUNCTION COST_FUNCTION(W,IU0)
      IMPLICIT NONE
      REAL(q)                :: COST_FUNCTION !value of cost function for W
      COMPLEX(q)             :: W(:,:)        !unitary matrix W
      INTEGER                :: IU0           !writing unit
      !local 
      INTEGER                :: L, M, K       !loop variables 
      COMPLEX(q),ALLOCATABLE :: AUX(:,:)      !auxilary matrix  W^+.S.W 
      COMPLEX(q),ALLOCATABLE :: CAUX(:,:)     !auxilary matrix  (W^+.S.W)*
      COMPLEX(q),ALLOCATABLE :: SAUX(:,:)     !auxilary matrix  (W^+.S.W)*
      COMPLEX(q),ALLOCATABLE :: A(:,:)        !auxilary matrix 
      REAL(q)                :: R             !function value 
 
      !for the following we will need NTARGET to be associated
      IF ( .NOT. ASSOCIATED(NTARGET) ) THEN
         IF ( IU0>=0) WRITE(IU0,'(A)') 'ERROR in COST_FUNCTION: NTARGET is not associated'
         STOP
      ENDIF    
      !also NTARGET should be at least >= NDIM
      IF ( SIZE(NTARGET) < NDIM ) THEN
         IF ( IU0>=0) WRITE(IU0,'(A)') 'ERROR in COST_FUNCTION: NTARGET is too small',&
            SIZE(NTARGET),NDIM
         STOP
      ENDIF 
 
      !allocate some memory for AUX and CAUX
      ALLOCATE(AUX(NDIM,NDIM), CAUX(NDIM,NDIM))
      ALLOCATE(A(NDIM,NDIM))
 
      ! decide if a set of matrices is optimized 
      IF ( ASSOCIATED(SK) ) THEN
         ALLOCATE(SAUX(NDIM,NDIM))
         SAUX=zero

         R=0._q
         DO K=1, NMAT
            A(:,:)=zero       
            AUX(:,:)=zero       
            SAUX(1:NDIM,1:NDIM)=SK(1:NDIM,1:NDIM,K)

            !compute A=Sk.W
            CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,SAUX(1,1),NDIM,W(1,1),NDIM,&
               zero,A(1,1),NDIM)
            CAUX=W 
            !compute AUX = W^+.Sk.W and CAUX = (W^+.Sk.W)*
            CALL ZGEMM('C','N',NDIM,NDIM,NDIM,one,CAUX(1,1),NDIM,A(1,1),NDIM,&
               zero,AUX(1,1),NDIM)
            CAUX=CONJG(AUX)

            !now we can evaluate the cost function
            !we define it as Jk=\sum_lm ( |(W^+.Sk.W )_lm|^2) 
            !where m is in ntarget and l is not in ntarget
            !we support only real function here!
            DO L=1, NDIM
               DO M=1, NDIM
                  IF ( NTARGET(M)==1 .AND. NTARGET(L) == 0) THEN            
                  !IF ( M>L) THEN   !Joint off-diagonal minimization        
                      R=R+REAL(CAUX(L,M)*AUX(L,M)+CAUX(M,L)*AUX(M,L),q)
                  ENDIF 
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(SAUX)
      ELSE
         A(:,:)=zero       
         AUX(:,:)=zero       
         !compute A=S.W
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,S(1,1),NDIM,W(1,1),NDIM,zero,A(1,1),NDIM)
         !CAUX=TRANSPOSE(W) 
         CAUX=W 
         !compute AUX = W^+.S.W and CAUX = (W^+.S.W)*
         CALL ZGEMM('C','N',NDIM,NDIM,NDIM,one,CAUX(1,1),NDIM,A(1,1),NDIM,zero,AUX(1,1),NDIM)
         CAUX=CONJG(AUX)

         !if this switch is used we maximize the trace of A'*A
         IF ( LMAXTR ) THEN
            !now we can evaluate the cost function
            !we define it as J=\sum_lm ( |(W^+.S.W )_lm|^2) 
            !where m is l 
            !we support only real functions here!
            R=0._q
            DO L=1, NDIM
               DO M=1, NDIM
                  IF ( M==L) THEN            
                      R=R+REAL(CAUX(L,M)*AUX(L,M)+CAUX(M,L)*AUX(M,L),q)
                 ENDIF 
               ENDDO
            ENDDO
         ELSE
            !now we can evaluate the cost function
            !we define it as J=\sum_lm ( |(W^+.S.W )_lm|^2) 
            !where m is in ntarget and l is not in ntarget
            !we support only real function here!
            R=0._q
            DO L=1, NDIM
               DO M=1, NDIM
                  IF ( NTARGET(M)==1 .AND. NTARGET(L) == 0) THEN            
                      R=R+REAL(CAUX(L,M)*AUX(L,M)+CAUX(M,L)*AUX(M,L),q)
                 ENDIF 
               ENDDO
            ENDDO
         ENDIF  
      ENDIF 
      !------------------------------------------------------------------------------------

      !get rid of auxilary matrices 
      DEALLOCATE(A,AUX,CAUX)

      !save value of cost function and get out 
      COST_FUNCTION=R

      RETURN
   ENDFUNCTION COST_FUNCTION

!***********************EUCLIDIAN_GRADIENT****************************
!
! this is the Euclidian gradient of the cost function at W 
!
!*********************************************************************

   FUNCTION EUCLIDIAN_GRADIENT(W,IU0)
      IMPLICIT NONE
      COMPLEX(q)  :: EUCLIDIAN_GRADIENT(NDIM,NDIM) !value of cost function for W
      COMPLEX(q)  :: W(:,:)                        !unitary matrix W
      INTEGER     :: IU0                           !writing unit
      !local 
      INTEGER                :: NK                 !loop variables 
      INTEGER                :: I,J,K,L            !loop variables 
      COMPLEX(q),ALLOCATABLE :: A(:,:)             !auxilary matrix W.S^*
      COMPLEX(q),ALLOCATABLE :: AC(:,:)            !auxilary matrix A^*
      COMPLEX(q),ALLOCATABLE :: B(:,:)             !auxilary matrix W^+.S.W
      COMPLEX(q),ALLOCATABLE :: BC(:,:)            !auxilary matrix B^*
      COMPLEX(q),ALLOCATABLE :: C(:,:)             !auxilary matrix S.W
      COMPLEX(q),ALLOCATABLE :: R(:,:)             !Function value 
 
      !for the following we will need NTARGET to be associated
      IF ( .NOT. ASSOCIATED(NTARGET) ) THEN
         IF ( IU0>=0) WRITE(IU0,'(A)') 'ERROR in EUCLIDIAN_GRADIENT: NTARGET is not associated'
         STOP
      ENDIF    
      !also NTARGET should be at least >= NDIM
      IF ( SIZE(NTARGET) < NDIM ) THEN
         IF ( IU0>=0) WRITE(IU0,'(A)') 'ERROR in EUCLIDIAN_GRADIENT: NTARGET is too small',&
            SIZE(NTARGET),NDIM
         STOP
      ENDIF 
 
      !allocate some memory for auxilary matrices A, AC, B, BC and C
      ALLOCATE(A(NDIM,NDIM),AC(NDIM,NDIM))
      ALLOCATE(B(NDIM,NDIM),BC(NDIM,NDIM))
      ALLOCATE(C(NDIM,NDIM), R(NDIM,NDIM))
      
      IF (ASSOCIATED(SK)) THEN
         !the gradient of cost function is  
         R=zero 
         DO NK=1, NMAT
            !compute A=W.Sk^* 
            A=zero       
            !compute AC=Sk*
            AC(:,:)=CONJG(SK(:,:,NK))
            !compute A=W.Sk^*
            CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,W(1,1),NDIM,AC(1,1),NDIM,zero,A(1,1),NDIM)
            !compute B=W^+.SkW
            B=W
            R=zero
            CALL ZGEMM('C','N',NDIM,NDIM,NDIM,one,B(1,1),NDIM,SK(1,1,NK),NDIM,&
               zero,R(1,1),NDIM)
            B=zero
            CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R(1,1),NDIM,W(1,1),NDIM,zero,B(1,1),NDIM)
            BC=CONJG(B)
            C=zero
            CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,SK(1,1,NK),NDIM,W(1,1),NDIM,&
               zero,C(1,1),NDIM)

            DO K=1,NDIM
               DO L=1,NDIM
                  DO I=1,NDIM
                      DO J=1,NDIM
                         IF ( NTARGET(I)==0 .AND. NTARGET(J)==1 ) THEN
                         !IF ( I>J ) THEN !Joint diagonalization
                         !IF ( I==J ) THEN !Joint diagonalization
                            R(K,L)=R(K,L)+A(I,L)*B(I,J)*DELTA(J,K)
                            R(K,L)=R(K,L)+DELTA(L,I)*BC(I,J)*C(K,J)
                            R(K,L)=R(K,L)+A(J,L)*B(J,I)*DELTA(I,K)
                            R(K,L)=R(K,L)+DELTA(J,L)*BC(J,I)*C(K,I)
                         ENDIF 
                      ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSE
         !compute A=W.S^* 
         A=zero       
         !compute AC=S*
         AC=CONJG(S)
         !compute A=W.S^*
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,W(1,1),NDIM,AC(1,1),NDIM,zero,A(1,1),NDIM)
         !compute B=W^+.S.W
         B=W
         R=zero
         CALL ZGEMM('C','N',NDIM,NDIM,NDIM,one,B(1,1),NDIM,S(1,1),NDIM,zero,R(1,1),NDIM)
         B=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R(1,1),NDIM,W(1,1),NDIM,zero,B(1,1),NDIM)
         BC=CONJG(B)
         C=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,S(1,1),NDIM,W(1,1),NDIM,zero,C(1,1),NDIM)
         
         !see FUNCTION COST_FUNCTION for explanation of this if clause
         IF ( LMAXTR ) THEN
            !the gradient of cost function is  
            R=zero 
            DO K=1,NDIM
               DO L=1,NDIM
                  DO I=1,NDIM
                      DO J=1,NDIM
                         IF ( I==J) THEN
                            R(K,L)=R(K,L)+A(I,L)*B(I,J)*DELTA(J,K)
                            R(K,L)=R(K,L)+DELTA(L,I)*BC(I,J)*C(K,J)
                            R(K,L)=R(K,L)+A(J,L)*B(J,I)*DELTA(I,K)
                            R(K,L)=R(K,L)+DELTA(J,L)*BC(J,I)*C(K,I)
                         ENDIF 
                      ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            !the gradient of cost function is  
            R=zero 
            DO K=1,NDIM
               DO L=1,NDIM
                  DO I=1,NDIM
                      DO J=1,NDIM
                         IF ( NTARGET(I)==0 .AND. NTARGET(J)==1 ) THEN
                            R(K,L)=R(K,L)+A(I,L)*B(I,J)*DELTA(J,K)
                            R(K,L)=R(K,L)+DELTA(L,I)*BC(I,J)*C(K,J)
                            R(K,L)=R(K,L)+A(J,L)*B(J,I)*DELTA(I,K)
                            R(K,L)=R(K,L)+DELTA(J,L)*BC(J,I)*C(K,I)
                         ENDIF 
                      ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

#ifdef debug 
WRITE(*,*)'GEk' 
DO I = 1,NDIM
   WRITE(*,'(16E16.6)')REAL(R(I,1:NDIM))
ENDDO
WRITE(*,*)' ' 
DO I = 1,NDIM
   WRITE(*,'(16E16.6)')IMAG(R(I,1:NDIM))
ENDDO
#endif

      EUCLIDIAN_GRADIENT=R
       
      !get rid of auxilary matrices 
      DEALLOCATE(A,AC,B,BC,C,R) 
      RETURN
   ENDFUNCTION EUCLIDIAN_GRADIENT

!*********************************************************************
!
!              L I N E  S E A R C H  A L G O R I T H M S
!
!*********************************************************************

!*************************ARMIJO_STEP_SIZE****************************
!
! This function finds the step size for the search algorithm using 
! Armijo's formula (see Ref. [2])
!
!*********************************************************************

   FUNCTION ARMIJO_STEP_SIZE(WK,DK,HK,OPT,IU0)
      IMPLICIT NONE
      REAL(q)          :: ARMIJO_STEP_SIZE !result 
      COMPLEX(q)       :: WK(:,:)          !unitary matrix
      COMPLEX(q)       :: HK(:,:)          !search direction matrix
      COMPLEX(q)       :: DK(:,:)          !Riemannian gradient
      CHARACTER(LEN=3) :: OPT              !minimization or maximization 
      INTEGER          :: IU0              !for writing 
      !local        
      INTEGER                :: N, M, I, J   !loop and auxillary
      REAL(q)                :: MU           !function value
      COMPLEX(q),ALLOCATABLE :: R1(:,:)      !matrix exponent
      COMPLEX(q),ALLOCATABLE :: R2(:,:)      !matrix exponent
      COMPLEX(q),ALLOCATABLE :: A(:,:)       !auxilary matrix
      COMPLEX(q),ALLOCATABLE :: B(:,:)       !auxilary matrix
      REAL(q)                :: JF           !cost function value
      REAL(q)                :: JF1          !cost function value
      REAL(q)                :: JF2          !cost function value
      REAL(q)                :: DK_HK        !<HK,DK>
      REAL(q)                :: ACC=1.E-12_q !threshold for values being equal
      INTEGER,PARAMETER      :: MAXIT1=100
      INTEGER,PARAMETER      :: MAXIT2=100
 
      !consistency check for WKa
      N=SIZE(WK,1)
      M=SIZE(WK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in ARMIJO_STEP_SIZE:'//&
            ' Size of WK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
      !consistency check for DK
      N=SIZE(DK,1)
      M=SIZE(DK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in ARMIJO_STEP_SIZE:'//&
            ' Size of DK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
      !consistency check for HK
      N=SIZE(HK,1)
      M=SIZE(HK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in ARMIJO_STEP_SIZE:'//&
            ' Size of HK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
        
      !allocate auxilary quantities
      ALLOCATE(R1(NDIM,NDIM),R2(NDIM,NDIM))
      ALLOCATE(A(NDIM,NDIM),B(NDIM,NDIM))

      !initialize step size to 1
      MU=1._q

      !decide wether we are minimizing or maximizing the cost function
      IF (OPT=='min'.OR.OPT=='MIN') THEN
         !towards the descent direction    
         A=zero
         A=-MU*SKEW(HK,IU0)
         R1=MATRIX_EXPONENTIAL(A,IU0)  !matrix exponential of strict skew-hermitian HK
         !R2 is it's square 
         R2=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,R1(1,1),NDIM,zero,R2(1,1),NDIM)
         !we need JF=J(WK),JF1=J(A=R1*WK),JF2=J(B=R2*WK)
         !A=R1*WK 
         A=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,zero,A(1,1),NDIM)
         !B=R2*WK 
         B=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R2(1,1),NDIM,WK(1,1),NDIM,zero,B(1,1),NDIM)
         JF =COST_FUNCTION(WK,IU0) 
         JF1=COST_FUNCTION(A,IU0)
         JF2=COST_FUNCTION(B,IU0)
         !check the result in order to save calculation time
         DK_HK=INNER_PRODUCT(DK,HK,IU0)
         IF ( (JF - JF2 >= MU*DK_HK) .AND. (JF - JF1 < 0.5_q*MU*DK_HK) ) THEN
            !in this case, basically nothing has to be done
            MU=MU
         ELSE
            !---------------------------------------------------- 
            !first condition: step size needs to be doubled?
            min_cond1: DO I = 1, MAXIT1          
               IF ( JF - JF2 >= MU*DK_HK) THEN   !that's the upper bound for step size
                  MU=2._q*MU
                  !souble step size by squaring the matrix
                  R1=R2      
                  R2=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,R1(1,1),NDIM,&
                     zero,R2(1,1),NDIM)
                  !in this case JF1 and JF2 need to be updated
                  !A=R1*WK 
                  A=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,&
                     zero,A(1,1),NDIM)
                  !B=R2*WK 
                  B=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R2(1,1),NDIM,WK(1,1),NDIM,&
                     zero,B(1,1),NDIM)
                  JF1=COST_FUNCTION(A,IU0)
                  JF2=COST_FUNCTION(B,IU0)
               ELSE  !otherwise get out
                  EXIT min_cond1
               ENDIF 
            ENDDO min_cond1
            !maybe we have not found upper bound for the step size  
            IF ( I > MAXIT1 ) THEN 
               IF (IU0>=0) WRITE(IU0,'(A,2I5,E13.4)')'ARMIJO_STEP_SIZE '//&
                  'reports: maximum number in min_cond1 exceeded ',I,MAXIT2,MU
            ENDIF 
            !---------------------------------------------------- 
             
            
            !---------------------------------------------------- 
            !second condition: step size is halved 
            min_cond2: DO I = 1, MAXIT2          
               IF ( JF - JF1 < 0.5_q*MU*DK_HK ) THEN  !that's the lower bound for step size
                  MU=0.5_q*MU
                  !again: need the matrix exponent  
                  A=zero
                  A=-MU*SKEW(HK,IU0)
                  R1=MATRIX_EXPONENTIAL(A,IU0) !matrix exponential of strict skew-hermitian HK
                  !and JF1=J(A=R1*WK) needs to be updated
                  !A=R1*WK 
                  A=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,&
                     zero,A(1,1),NDIM)
                  JF1=COST_FUNCTION(A,IU0)
               ELSE  !otherwise get out
                  EXIT min_cond2
               ENDIF 
               !also get out if MU is too small 
               IF ( ABS(MU) < ACC) THEN
                  EXIT min_cond2
               ENDIF 
 
            ENDDO min_cond2
            !maybe we have not found lower bound for the step size  
            IF ( I > MAXIT2 ) THEN 
               IF (IU0>=0) WRITE(IU0,'(A,2I5,E13.4)')'ARMIJO_STEP_SIZE '//&
                  'reports: maximum number in min_cond2 exceeded ',I,MAXIT2,MU
            ENDIF 
            !---------------------------------------------------- 
         ENDIF 
      ELSEIF(OPT=='max'.OR.OPT=='MAX') THEN 
         !towards the descent direction    
         A=zero
         A=MU*SKEW(HK,IU0)
         R1=MATRIX_EXPONENTIAL(A,IU0)  !matrix exponential of strict skew-hermitian HK
         !R2 is it's square 
         R2=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,R1(1,1),NDIM,zero,R2(1,1),NDIM)
         !we need JF=J(WK),JF1=J(A=R1*WK),JF2=J(B=R2*WK)
         !A=R1*WK 
         A=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,zero,A(1,1),NDIM)
         !B=R2*WK 
         B=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R2(1,1),NDIM,WK(1,1),NDIM,zero,B(1,1),NDIM)
         JF =COST_FUNCTION(WK,IU0) 
         JF1=COST_FUNCTION(A,IU0)
         JF2=COST_FUNCTION(B,IU0)
         !check the result in order to save calculation time
         DK_HK=INNER_PRODUCT(DK,HK,IU0)
         IF ( (JF - JF2 <= -MU*DK_HK) .AND. (JF - JF1 > 0.5_q*MU*DK_HK) ) THEN
            !in this case, basically nothing has to be done
            MU=MU
         ELSE
            !---------------------------------------------------- 
            !check first condition: step size needs to be doubled?
            max_cond1: DO I = 1, MAXIT1          
               IF ( JF - JF2 <= -MU*DK_HK) THEN   !that's the upper bound for step size
                  MU=2._q*MU
                  !souble step size by squaring the matrix
                  R1=R2      
                  R2=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,R1(1,1),NDIM,&
                     zero,R2(1,1),NDIM)
                  !in this case JF1 and JF2 need to be updated
                  !A=R1*WK 
                  A=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,&
                     zero,A(1,1),NDIM)
                  !B=R2*WK 
                  B=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R2(1,1),NDIM,WK(1,1),NDIM,&
                     zero,B(1,1),NDIM)
                  JF1=COST_FUNCTION(A,IU0)
                  JF2=COST_FUNCTION(B,IU0)
               ELSE  !otherwise get out
                  EXIT max_cond1
               ENDIF 
            ENDDO max_cond1
            !maybe we have not found upper bound for the step size  
            IF ( I > MAXIT1 ) THEN 
               IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ARMIJO_STEP_SIZE '//&
                  'reports: maximum number in max_cond1 exceeded ',I,MAXIT2
            ENDIF 
            !---------------------------------------------------- 
             
            
            !---------------------------------------------------- 
            !check second condition: step size is halved 
            max_cond2: DO I = 1, MAXIT2          
               IF ( JF - JF1 > -0.5_q*MU*DK_HK ) THEN  !that's the lower bound for step size
                  MU=0.5_q*MU
                  !again: need the matrix exponent  
                  A=zero
                  A=MU*SKEW(HK,IU0)
                  R1=MATRIX_EXPONENTIAL(A,IU0)  !matrix exponential of strict skew-hermitian HK
                  !and JF1=J(A=R1*WK) needs to be updated
                  !A=R1*WK 
                  A=zero
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R1(1,1),NDIM,WK(1,1),NDIM,&
                     zero,A(1,1),NDIM)
                  JF1=COST_FUNCTION(A,IU0)
               ELSE  !otherwise get out
                  EXIT max_cond2
               ENDIF 
               IF ( ABS(MU) < ACC) THEN
                  EXIT max_cond2
               ENDIF 
            ENDDO max_cond2
            !maybe we have not found lower bound for the step size  
            IF ( I > MAXIT2 ) THEN 
               IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ARMIJO_STEP_SIZE '//&
                  'reports: maximum number in max_cond2 exceeded ',I,MAXIT2
            ENDIF 
            !---------------------------------------------------- 
         ENDIF 
      ELSE
         IF ( IU0>=0) WRITE(IU0,'(A)')'ERROR in ARMIJO_STEP_SIZE: '//&
            'OPT not set properly '//OPT//'.'
         STOP
      ENDIF 

      !final step size is MU at this stage
      ARMIJO_STEP_SIZE=MU  
      RETURN
   ENDFUNCTION ARMIJO_STEP_SIZE

!***************************POLY_STEP_SIZE****************************
!
! This function finds the step size for the search algorithm using 
! a polynomial approximation of 
!   
!               J(W(mu))=J(E^(mu H)W)  
!
! differentiating w.r.t. mu gives then
! J'(W(mu))=-2 Tr(DK.W(mu).W'*W(mu)'*H)
!
! Choose this option if you are close to an extremum.
! see Eq. 14 of Ref. [2] for more details
!
!*********************************************************************

   FUNCTION POLY_STEP_SIZE(WK,HK,NPOLY,OPT,IU0)
      IMPLICIT NONE
      REAL(q)          :: POLY_STEP_SIZE !result 
      COMPLEX(q)       :: WK(:,:)        !unitary matrix
      COMPLEX(q)       :: HK(:,:)        !search direction matrix
      INTEGER          :: NPOLY          !polynomial order 
      CHARACTER(LEN=3) :: OPT            !minimization or maximization 
      INTEGER          :: IU0            !for writing 
      !local        
      REAL(q),PARAMETER   :: THRESH=1.E-10_q !threshold for derivative 
      LOGICAL                :: LUNITARY     !consistency check
      INTEGER                :: N, M, I, J   !loop and auxillary
      REAL(q)                :: MU           !function value
      REAL(q)                :: SIGN_MU      !sign of mu
      REAL(q)                :: OMEGA_MAX    !largest eigenvalue
      REAL(q)                :: T_MU         !period 
      REAL(q)                :: MU_STEP      !uniform step size 
      COMPLEX(q),ALLOCATABLE :: RP(:,:)      !matrix exponent of polynom
      COMPLEX(q),ALLOCATABLE :: RM(:,:)      !matrix exponent of stepsize
      COMPLEX(q),ALLOCATABLE :: A(:,:)       !auxilary matrix
      REAL(q)                :: JF           !cost function value
      INTEGER                :: NITER        !poly loop
      COMPLEX(q),ALLOCATABLE :: DK(:,:)      !Euclidian gradient
      COMPLEX(q),ALLOCATABLE :: B(:,:)       !auxilary matrix
      REAL(q),ALLOCATABLE    :: D1POLY(:)    !gradient at different points
      REAL(q),ALLOCATABLE    :: D1POLYRED(:) !reduced gradient at different points 
      REAL(q),ALLOCATABLE    :: MU_I(:)      !points where derivative DK is evaluated 
      REAL(q),ALLOCATABLE    :: C(:,:)       !Vandermonde matrix 
      REAL(q),ALLOCATABLE    :: COEFF(:)     !coefficient vector of polynom
      REAL(q),ALLOCATABLE    :: REALROOTS(:) !real roots 
      INTEGER                :: NROOTS       !number of real roots  
      INTEGER                :: NPR          !number of positive roots  
      REAL(q),ALLOCATABLE    :: PR(:)        !positive roots 

      !consistency check for WK
      N=SIZE(WK,1)
      M=SIZE(WK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in POLY_STEP_SIZE:'//&
            ' Size of WK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
      !consistency check for HK
      N=SIZE(HK,1)
      M=SIZE(HK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in POLY_STEP_SIZE:'//&
            ' Size of HK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
        
      !also check if WK is unitary  and HK is skew hermitian
      LUNITARY=LUNITCHECK(WK,IU0)
      IF ( .NOT. LUNITARY ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,L4)')'ERROR in POLY_STEP_SIZE:'//&
            ' WK is not unitary ',LUNITARY
         STOP
      ENDIF 

      !skew hermiticity 
      ALLOCATE(A(NDIM,NDIM))
      A=HK+CONJG(TRANSPOSE(HK))
      JF=NORM(A,2,IU0)
      IF ( JF > TINY_U ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,E10.3)')'ERROR in POLY_STEP_SIZE:'//&
            ' HK is not skew hermitian:',JF
         STOP
      ENDIF 

      !also NPOLY needs to be checked.
      !we allow NPOLY to be an integer between 3 and 5
      IF ( NPOLY < 3 .OR. NPOLY > 5 ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,I4)')'ERROR in POLY_STEP_SIZE:'//&
            ' NPOLY needs to be an integer between 3 and 5:',NPOLY
      ENDIF 

      !are we computing the minimum of the maximum of the cost function 
      IF (OPT=='min'.OR.OPT=='MIN') THEN
         SIGN_MU=-1._q  !negative gradient
      ELSEIF(OPT=='max'.OR.OPT=='MAX') THEN 
         SIGN_MU=1._q   !positve gradient
      ELSE
         IF ( IU0>=0) WRITE(IU0,'(A)')'ERROR in POLY_STEP_SIZE: '//&
            'OPT not set properly '//OPT//'.'
         STOP
      ENDIF 

      !we need the largest eigenvalue 
      OMEGA_MAX=NORM(HK,2,IU0)
   
      OMEGA_MAX=MAX_EIG(HK,IU0)

      !determine the (almost) period T_MU
      T_MU=TPI/(QDIM*OMEGA_MAX)
   
      !determine uniform step size
      MU_STEP=T_MU/REAL(NPOLY,q)
  
      !determine roation matrix RP 
      ALLOCATE(RP(NDIM,NDIM))
      A=CMPLX(SIGN_MU*MU_STEP,0._q,q)*HK
      RP=MATRIX_EXPONENTIAL(A,IU0)

      !rotation is inialized to unity  
      ALLOCATE(RM(NDIM,NDIM))
      RM=zero
      DO I=1,NDIM
         RM(I,I)=one
      ENDDO

      !we need some memory for the Euclidian gradient
      ALLOCATE(DK(NDIM,NDIM))
      ALLOCATE(B(NDIM,NDIM))
      ALLOCATE(D1POLY(NPOLY+1))

      DO NITER = 1, NPOLY+1
         !determine the Euclidian gradient at RM*WK
         A=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,RM(1,1),NDIM,WK(1,1),NDIM,&
            zero,A(1,1),NDIM)
         DK = EUCLIDIAN_GRADIENT(A,IU0)
         !get the coefficients as the derivative of DK*Wk'*RM'*HK'
         A=zero
         B=WK
         CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,DK(1,1),NDIM,B(1,1),NDIM,&
            zero,A(1,1),NDIM)
         B=RM
         DK=zero
         CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,A(1,1),NDIM,B(1,1),NDIM,&
            zero,DK(1,1),NDIM)
         B=HK
         A=zero
         CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,DK(1,1),NDIM,B(1,1),NDIM,&
            zero,A(1,1),NDIM)
         D1POLY(NITER)=-2._q*SIGN_MU*TRACE_OF(A,IU0)
         !obtain RM=RM*RP     
         A=RM
         B=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,A(1,1),NDIM,RP(1,1),NDIM,&
            one,B(1,1),NDIM)
         RM=B
#ifdef debug
IF(NITER==2) THEN
WRITE(*,'(A,F30.16)')'RM',D1POLY(NITER)
DO I=1,NDIM
  WRITE(*,'(16F30.16)')REAL(RM(I,1:NDIM),q)
ENDDO
WRITE(*,*)' ' 
DO I=1,NDIM
  WRITE(*,'(16F30.16)')IMAG(RM(I,1:NDIM))
ENDDO
ENDIF
#endif 
      ENDDO 
      !---------------------------------------------------------------
      !if the derivative is too small MU = 0  
      IF ( ABS(D1POLY(1)) < THRESH ) THEN
      !---------------------------------------------------------------
         MU=0._q
      !---------------------------------------------------------------
      ELSE   !otherwise continue and find smallest positive root of polynom
      !---------------------------------------------------------------
         ! determine the NPOLY+1 uniform evaluation points MU_I
         ! for the Euclidian gradient DK
         ALLOCATE(MU_I(NPOLY+1))
         DO I=1,NPOLY+1
            MU_I(I) = 0._q + (I-1)*T_MU/REAL(NPOLY,q)
         ENDDO      

         !allocate memory for the Vandermonde matrix 
         ALLOCATE(C(NPOLY,NPOLY))
         !reduced coefficient vector 
         ALLOCATE(D1POLYRED(NPOLY))
    
         !determine Vandermonde matrix (without 1 column)
         DO I = 2, NPOLY+1
            DO J = 1, NPOLY
               C(I-1,J)=MU_I(I)**J
            ENDDO 
         ENDDO 
#ifdef debug
WRITE(*,*)' VANDERMONDE' 
DO I=1,NPOLY
  WRITE(*,'(16E12.4)')C(I,1:NPOLY)
ENDDO
#endif 
         !set up reduced coefficient vector
         DO J = 1, NPOLY
            D1POLYRED(J) = D1POLY(J+1)-D1POLY(1)
         ENDDO 
        
         !solve the real linear system of equations
         CALL R_LINEAR_SOLVE(C,D1POLYRED,IU0)

         !save coefficient in reverse order to COEFF
         ALLOCATE(COEFF(NPOLY+1))
         DO I=1,NPOLY
            COEFF(I)=D1POLYRED(NPOLY-I+1)
         ENDDO
         !attach the first coefficient to the end
         COEFF(NPOLY+1)=D1POLY(1)
   
         !find all real valued roots of polynomial,
         !which is given by the coefficient vector COEFF (NPOLY in total)
         ALLOCATE(REALROOTS(NPOLY))
         CALL R_POLY_ROOTS(COEFF,REALROOTS,NROOTS,IU0)

         !if real roots found than 
         IF( NROOTS > 0 ) THEN
            !find positive roots 
            ALLOCATE(PR(NROOTS))
            NPR=0
            DO I = 1, NROOTS 
               IF ( REALROOTS(I) > 0 ) THEN
                  NPR=NPR+1
                  PR(NPR)=REALROOTS(I)
               ENDIF 
            ENDDO 
            !if no positive root found set MU to zero 
            IF ( NPR == 0 ) THEN
               MU = 0._q
            ELSE  !otherwise find smallest positive root 
               MU=MINVAL(PR(1:NPR)) 
            ENDIF 
            !get rid of PR
            DEALLOCATE(PR)
         ELSE  !if no real root found than MU = 0 
            MU = 0._q
         ENDIF
          
         !kill arrays
         DEALLOCATE(REALROOTS)     
         DEALLOCATE(COEFF)  
         DEALLOCATE(D1POLYRED)
         DEALLOCATE(C)
         DEALLOCATE(MU_I)
      !---------------------------------------------------------------
      ENDIF 
      !---------------------------------------------------------------

      !get rid of arrays 
      DEALLOCATE(D1POLY)
      DEALLOCATE(B)
      DEALLOCATE(DK)
      DEALLOCATE(RM)
      DEALLOCATE(RP)
      DEALLOCATE(A)

      !final step size is MU at this stage
      POLY_STEP_SIZE=MU  
      RETURN
   ENDFUNCTION POLY_STEP_SIZE

!****************************FFT_STEP_SIZE****************************
!
! This function finds the step size for the search algorithm using 
! the FFT method from Ref. [2].
! Choose this option if you are far away from an extremum.
!
!*********************************************************************

   FUNCTION FFT_STEP_SIZE(WK,HK,NT,OPT,IU0)
      IMPLICIT NONE
      REAL(q)          :: FFT_STEP_SIZE !result 
      COMPLEX(q)       :: WK(:,:)        !unitary matrix
      COMPLEX(q)       :: HK(:,:)        !search direction matrix
      INTEGER          :: NT             !number of considered periods 
      CHARACTER(LEN=3) :: OPT            !minimization or maximization 
      INTEGER          :: IU0            !for writing 
      !local        
      INTEGER,PARAMETER:: KSAMPLE=3          !# of samples (use 3 for now)
      REAL(q),PARAMETER:: MAX_DIST=0.01_q    !maximum distance to unit circle     
      REAL(q),PARAMETER:: DJ_THRESH=1.E-12_q !threshold for the derivative    
      INTEGER                :: N, M, I      !loop and auxillary
      LOGICAL                :: LUNITARY     !consistency check
      COMPLEX(q),ALLOCATABLE :: A(:,:)       !auxilary matrix
      REAL(q)                :: JMIN         !cost function value minimum
      REAL(q)                :: MU           !function value
      REAL(q)                :: SIGN_MU      !sign of mu
      REAL(q)                :: OMEGA_MAX    !largest eigenvalue
      REAL(q)                :: T_MU         !period 
      REAL(q)                :: T_FFT        !length of FFT interval 
      INTEGER                :: N_FFT        !FFT size 
      COMPLEX(q),ALLOCATABLE :: R_FFT(:,:)   !rotation matrix 
      COMPLEX(q),ALLOCATABLE :: R_MU(:,:)    !rotation matrix at different MU
      REAL(q),ALLOCATABLE    :: J(:)         !cost function at equidistant points 
      REAL(q),ALLOCATABLE    :: D1J(:)       !derivative (w.r.t. MU) of J 
      COMPLEX(q),ALLOCATABLE :: D1J_FFT(:)   !FFT of derivative (w.r.t. MU) of J 
      COMPLEX(q),ALLOCATABLE :: DK(:,:)      !Euclidian gradient
      COMPLEX(q),ALLOCATABLE :: B(:,:)       !auxilary matrix
      REAL(q),ALLOCATABLE    :: H_WINDOW(:)  !Hanning window steps 
      COMPLEX(q),ALLOCATABLE :: CROOTS(:)    !roots 
      COMPLEX(q),ALLOCATABLE :: CIRCROOTS(:) !roots on unit circle
      REAL(q),ALLOCATABLE    :: PHASEROOTS(:)!phase of roots on unit circle
      REAL(q),ALLOCATABLE    :: MUROOTS(:)   !roots w.r.t. mu on unit circle
      REAL(q),ALLOCATABLE    :: MU_FFT(:)    !equispaced evaluation points 
      INTEGER                :: IMIN(1)      !position of minimum in J_FFT
      REAL(q)                :: MU_FFT_MIN   !Minimum of J w.e.t. FFT

      !consistency check for WK
      N=SIZE(WK,1)
      M=SIZE(WK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in FFT_STEP_SIZE:'//&
            ' Size of WK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
      !consistency check for HK
      N=SIZE(HK,1)
      M=SIZE(HK,2)
      IF ( N/=M .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in FFT_STEP_SIZE:'//&
            ' Size of HK inconsistent:',N,M,NDIM
         STOP
      ENDIF 
      !number of considered periods has to be a positive integer
      IF ( NT <= 0 ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,I5)')'ERROR in FFT_STEP_SIZE:'//&
            ' # of considered period has to be a positive number',NT
         STOP
      ENDIF 
      !are we computing the minimum of the maximum of the cost function 
      IF (OPT=='min'.OR.OPT=='MIN') THEN
         SIGN_MU=-1._q  !negative gradient
      ELSEIF(OPT=='max'.OR.OPT=='MAX') THEN 
         SIGN_MU=1._q   !positve gradient
      ELSE
         IF ( IU0>=0) WRITE(IU0,'(A)')'ERROR in POLY_STEP_SIZE: '//&
            'OPT not set properly '//OPT//'.'
         STOP
      ENDIF 

      !determine the almost period T_MU
      !for this we need the largest eigenvalue 
      OMEGA_MAX=NORM(HK,2,IU0)
      !the (almost) period T_MU is then 
      T_MU=TPI/(QDIM*OMEGA_MAX)
       
      !determine the length of the FFT interval 
      T_FFT=NT*T_MU

      !FFT size, ensure odd value
      N_FFT=2*FLOOR(KSAMPLE*NT/2._q)+1

      !we need to evaluate the rotation matrix on a uniformly space grid 
      !this is done by evaluating the step matrix R_FFT once 
      !as an matrix exponential, and using matrix multiplication for R_MU
      ALLOCATE(A(NDIM,NDIM))
      ALLOCATE(R_FFT(NDIM,NDIM),R_MU(NDIM,NDIM)) 
      A=SIGN_MU*(T_FFT/REAL(N_FFT,q))*HK
      R_FFT=MATRIX_EXPONENTIAL(A,IU0)
      !initialize R_MU to identity (MU=0)
      R_MU=zero
      DO I = 1, NDIM
         R_MU(I,I)=one 
      ENDDO
#ifdef debug
WRITE(*,*)'R_FFT'
DO I = 1, NDIM 
   WRITE(*,'(16E12.3)')REAL(R_FFT(I,1:NDIM),q)
ENDDO
WRITE(*,*)''
DO I = 1, NDIM 
   WRITE(*,'(16E12.3)')IMAG(R_FFT(I,1:NDIM))
ENDDO
#endif    
      !evaluate the cost function at the points -> J     
      !evaluate the derivative dJ/dMu at the points -> D1J     
      ALLOCATE(J(N_FFT),D1J(N_FFT))
      J=0._q
      D1J=0._q

      !we need some memory for the Euclidian gradient
      ALLOCATE(DK(NDIM,NDIM))
      ALLOCATE(B(NDIM,NDIM))
      DO I = 1, N_FFT
         !determine the Euclidian gradient at RM*WK
         A=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,R_MU(1,1),NDIM,WK(1,1),NDIM,&
            zero,A(1,1),NDIM)
         DK = EUCLIDIAN_GRADIENT(A,IU0)
         !get the coefficients as the derivative of DK*Wk'*R_MU'*HK'
         A=zero
         B=WK
         CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,DK(1,1),NDIM,B(1,1),NDIM,&
            zero,A(1,1),NDIM)
         B=R_MU
         DK=zero
         CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,A(1,1),NDIM,B(1,1),NDIM,&
            zero,DK(1,1),NDIM)
         B=HK
         CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,DK(1,1),NDIM,B(1,1),NDIM,&
            zero,A(1,1),NDIM)
         !derivative of J at R_MU*WK
         D1J(I)=-2*SIGN_MU*TRACE_OF(A,IU0)
         A=WK
         DK=R_MU
         B=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,DK(1,1),NDIM,A(1,1),NDIM,&
            zero,B(1,1),NDIM)
         !J at R_MU*WK
         J(I)=COST_FUNCTION(B,IU0)
         !R_MU = R_MU*R_FFT
         A=R_MU
         B=zero
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,one,A(1,1),NDIM,R_FFT(1,1),NDIM,&
            one,B(1,1),NDIM)
         R_MU=B
      ENDDO 

#ifdef debug
WRITE(*,*)'J and D1J'
DO I = 1, N_FFT
   WRITE(*,'(2E20.10)')J(I),D1J(I)
ENDDO
#endif    
 
      !compute N_FFT steps of a Hanning window
      ALLOCATE(H_WINDOW(N_FFT))
      H_WINDOW=HANNING(N_FFT,IU0)
      !and multiplicate derivative with Hanning window
      DO I = 1, N_FFT
         D1J(I)=D1J(I)*H_WINDOW(I)
      ENDDO

      !Fourier transform the vector D1J
      !we therefore need a complex version of D1J
      ALLOCATE(D1J_FFT(N_FFT))
      CALL R_FFTSHIFT(D1J,D1J_FFT,IU0)

#ifdef debug
WRITE(*,*)'FFT'
DO I =1 , N_FFT
   WRITE(*,'(I4,3E20.10,"i")')I,D1J(I),D1J_FFT(I)
ENDDO
#endif 
      !find complex roots of polynomial 
      ALLOCATE(CROOTS(N_FFT-1))
      !stupidly, this routine needs a reverse order of the array 
      CALL C_POLY_ROOTS(D1J_FFT(N_FFT:1:-1),CROOTS,IU0)        

      !find roots closest to unit circle 
      ALLOCATE(CIRCROOTS(N_FFT-1))
      ALLOCATE(PHASEROOTS(N_FFT-1))
      CIRCROOTS=zero
      PHASEROOTS=0._q
      M=0
      DO I = 1, N_FFT-1
         IF ( ABS(ABS(CROOTS(I))-1._q)<MAX_DIST ) THEN
            M = M+1
            CIRCROOTS(M) = CROOTS(I)
            !find argument as well 
            PHASEROOTS(M) = IMAG(LOG(CROOTS(I)))/(2*PI)
         ENDIF 
      ENDDO

      !we continue only if root on circle are found 
      IF ( M > 0 ) THEN 
         !find zero-crossing values corresponding to the roots on the circle
         ALLOCATE(MUROOTS(M))
         MUROOTS=0._q

         !shift roots to [0,T_FFT]
         MUROOTS=MODULO(PHASEROOTS*T_FFT,T_FFT)
         !sort them 
         CALL SORT(M,MUROOTS,IU0) 
         !and use only odd entries 
         PHASEROOTS=0._q
         PHASEROOTS=MUROOTS
         MUROOTS=0._q
         N=0
         DO I = 1, M, 2
            N=N+1
            MUROOTS(N)=PHASEROOTS(I) 
         ENDDO
      ELSE      
         !if no roots on circle are found their total number is N=0 
         N=0
      ENDIF  
      IF ( (N <= 0).OR.(ABS(D1J(1))< DJ_THRESH) ) THEN
         ! in case no root found or derivative is too small, step size is 0
         MU=0._q
      ELSE 
         !generate array of equally spaced points in interval [0,NT*T_MU]
         ALLOCATE(MU_FFT(N_FFT))
         MU_FFT=0._q
         DO I=1,N_FFT
            MU_FFT(I) = (I-1)*NT*T_MU/REAL(N_FFT,q)
         ENDDO
          
         !now decide, what kind of extremum is needed
         IF ( SIGN_MU < 0 ) THEN   !in case of a minimum of J is needed
            !global minimum is 
            JMIN=MINVAL(J)       
            IMIN=MINLOC(J) !corresponding index of 'global' minimum
            !in case something went wrong (should not be the case :)
            IF ( IMIN(1) > N_FFT) THEN
               IF (IU0>=0) WRITE(IU0,'(A,I6)')'ERROR in FFT_STEP_SIZE:'//&
                  'IMIN(1) (of min) beyond bounds',IMIN(1),N_FFT
               STOP
            ENDIF
            !this is the MU value corresponding to global minimum
            MU_FFT_MIN=MU_FFT(IMIN(1))
            !which index of equispaced grid is closest to MUROOTS -> optimum MU
            IMIN=MINLOC(ABS(MU_FFT_MIN-MUROOTS(1:N))) 
            !in case something went wrong (should not be the case :)
            IF ( IMIN(1) > N) THEN
               IF (IU0>=0) WRITE(IU0,'(A,I6)')'ERROR in FFT_STEP_SIZE:'//&
                  'IMIN(1) (of min) beyond bounds (2)',IMIN(1),N
               STOP
            ENDIF
            MU=MUROOTS(IMIN(1))
         ELSE                      !in case of a maximum of J is needed 
            !global maximum is 
            JMIN=MAXVAL(J)       
            IMIN=MAXLOC(J) !corresponding index of 'global' maximum
            !in case something went wrong (should not be the case :)
            IF ( IMIN(1) > N_FFT) THEN
               IF (IU0>=0) WRITE(IU0,'(A,I6)')'ERROR in FFT_STEP_SIZE:'//&
                  'IMIN(1) (of max) beyond bounds',IMIN(1),N_FFT
               STOP
            ENDIF
            !this is the MU value corresponding to global minimum
            MU_FFT_MIN=MU_FFT(IMIN(1))
            !which index of equispaced grid is closest to MUROOTS -> optimum MU
            IMIN=MINLOC(ABS(MU_FFT_MIN-MUROOTS(1:N))) 
            !in case something went wrong (should not be the case :)
            IF ( IMIN(1) > N) THEN
               IF (IU0>=0) WRITE(IU0,'(A,I6)')'ERROR in FFT_STEP_SIZE:'//&
                  'IMIN(1) (of max) beyond bounds (2)',IMIN(1),N
               STOP
            ENDIF
            MU=MUROOTS(IMIN(1))
         ENDIF 
      ENDIF 

      !final step size is MU at this stage
      FFT_STEP_SIZE=MU  

      !get rid of all allocated arrays
      IF(ALLOCATED(A))DEALLOCATE(A)
      IF(ALLOCATED(R_FFT)) DEALLOCATE(R_FFT)
      IF(ALLOCATED(R_MU))DEALLOCATE(R_MU)
      IF(ALLOCATED(J)) DEALLOCATE(J)
      IF(ALLOCATED(D1J)) DEALLOCATE(D1J)
      IF(ALLOCATED(DK)) DEALLOCATE(DK)
      IF(ALLOCATED(B)) DEALLOCATE(B)
      IF(ALLOCATED(H_WINDOW)) DEALLOCATE(H_WINDOW)
      IF(ALLOCATED(D1J_FFT)) DEALLOCATE(D1J_FFT)
      IF(ALLOCATED(CROOTS)) DEALLOCATE(CROOTS)
      IF(ALLOCATED(CIRCROOTS)) DEALLOCATE(CIRCROOTS)
      IF(ALLOCATED(PHASEROOTS)) DEALLOCATE(PHASEROOTS)
      IF(ALLOCATED(MUROOTS)) DEALLOCATE(MUROOTS)
      IF(ALLOCATED(MU_FFT)) DEALLOCATE(MU_FFT)

      RETURN
   ENDFUNCTION FFT_STEP_SIZE

!*********************************************************************
!
!                A U X I L A R Y  F U N C T I O N S 
!
!*********************************************************************

!***********************RIEMANNIAN_GRADIENT***************************
!
! this is the Riemannian gradient of the cost function at W 
! (needs the Euclidian gradient DE)
!
!*********************************************************************

   FUNCTION RIEMANNIAN_GRADIENT(DE,W,IU0)
      IMPLICIT NONE
      COMPLEX(q)  :: RIEMANNIAN_GRADIENT(NDIM,NDIM)!value of cost function for W
      COMPLEX(q)  :: DE(:,:)                       !Euclidian gradient
      COMPLEX(q)  :: W(:,:)                        !point on Lie manifold
      INTEGER     :: IU0                           !writing unit
      !local 
      COMPLEX(q),ALLOCATABLE :: A(:,:)             !auxilary matrix
      COMPLEX(q),ALLOCATABLE :: R(:,:)             !Function value 
      COMPLEX(q),ALLOCATABLE :: L(:,:)             !Function value 
      INTEGER :: I 
 
      !allocate some memory for auxilary matrices A, AC, B, BC and C
      ALLOCATE(R(NDIM,NDIM), A(NDIM,NDIM), L(NDIM, NDIM))
      
      R=zero       
!      A=TRANSPOSE(W)
      A=W
      !compute R=DE.W^+
      CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,DE(1,1),NDIM,A(1,1),NDIM,zero,R(1,1),NDIM)
      L=zero       
!      A=TRANSPOSE(DE)
      A=DE
      !compute R=DE.W^+ - W.DE^+
      CALL ZGEMM('N','C',NDIM,NDIM,NDIM,one,W(1,1),NDIM,A(1,1),NDIM,zero,L(1,1),NDIM)

#ifdef debug
WRITE(*,*)' DE.W^+ ' 
      DO I=1,NDIM
         WRITE(*,'(16F20.10)')REAL(R(I,1:NDIM),q)
      ENDDO 
      DO I=1,NDIM
         WRITE(*,'(16F20.10)')IMAG(R(I,1:NDIM))
      ENDDO 
WRITE(*,*)' W.DE^+ ' 
      DO I=1,NDIM
         WRITE(*,'(16F20.10)')REAL(L(I,1:NDIM),q)
      ENDDO 
      DO I=1,NDIM
         WRITE(*,'(16F20.10)')IMAG(L(I,1:NDIM))
      ENDDO 
#endif 
      
      RIEMANNIAN_GRADIENT=R-L
       
      !get rid of auxilary matrices 
      DEALLOCATE(A,R,L) 

      RETURN

   ENDFUNCTION RIEMANNIAN_GRADIENT

!*****************************TRACE_OF********************************
!
! computes the (real) race of a matrix W
!
!*********************************************************************

   FUNCTION TRACE_OF(W,IU0)
      IMPLICIT NONE
      REAL(q)     :: TRACE_OF!result 
      COMPLEX(q)  :: W(:,:)               !unitary matrix W
      INTEGER     :: IU0                  !writing unit
      !local 
      INTEGER                :: I,J,K,L   !loop variables 
      COMPLEX(q),ALLOCATABLE :: A(:,:)    !auxilary matrix 
      REAL(q)                :: R         !temporary result
      INTEGER                :: N

      !some consistency check 
      N=SIZE(W,1)
      IF ( N/=SIZE(W,2) .OR. N /= NDIM ) THEN 
         IF ( IU0>=0 ) WRITE(IU0,'(A,3I4)')'ERROR in TRACE_OF:'//&
            ' Size of W inconsistent:',N,SIZE(W,2),NDIM
         STOP
      ENDIF 
 
      !allocate some memory for auxilary matrix A
      ALLOCATE(A(NDIM,NDIM))
      A=W      

      !the gradient of cost function is  
      R=0._q
      DO I=1,NDIM
         R=R+REAL(A(I,I),q)
      ENDDO

      TRACE_OF=R
       
      !get rid of auxilary matrices 
      DEALLOCATE(A) 

      RETURN
   ENDFUNCTION TRACE_OF 

!*********************************SKEW********************************
!
! this function makes a strict skew hermitian matrix of A
! recall: skew hermitian matrices are a basis for the Lie algebra 
!
!*********************************************************************

   FUNCTION SKEW(A,IU0)
      IMPLICIT NONE
      COMPLEX(q)    :: SKEW(NDIM,NDIM)
      COMPLEX(q)    :: A(:,:)
      INTEGER       :: IU0
      !local        
      INTEGER                :: N, M, I, J
      COMPLEX(q),ALLOCATABLE :: R(:,:)
     
      N=SIZE(A,1)
      M=SIZE(A,2)
      IF( N/=M ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I4)')'ERROR in SKEW:'//&
            'A is not a square matrix',N,M
         STOP
      ENDIF 
      
      IF( N>NDIM ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I4)')'WARNING from SKEW:'//&
            'input matrix A is not of size NDIM',N,NDIM
      ENDIF 
 
      ALLOCATE(R(NDIM,NDIM))

      DO I=1,NDIM
         DO J=1,NDIM
            R(I,J)=CONJG(A(J,I))
         ENDDO
      ENDDO
       
      R=(A-R)/2._q
  
      SKEW=R

      DEALLOCATE(R)
      RETURN
   ENDFUNCTION SKEW

!**************************MATRIX_EXPONENTIAL*************************
!
! this function gives the matrix exponential of a matrix A.
! use this function only for skew-hermitian matrices A, because
! only then E^A coincides with the exponential map on U(n).
! E^A is computed using a scaling and squaring algorithm with a 
! Pade approximation. 
!
!*********************************************************************

   FUNCTION MATRIX_EXPONENTIAL(A,IU0)
      COMPLEX(q)    :: MATRIX_EXPONENTIAL(NDIM,NDIM)
      COMPLEX(q)    :: A(:,:)
      INTEGER       :: IU0
      !local        
      INTEGER                :: N, M, I, J
      COMPLEX(q),ALLOCATABLE :: R(:,:)
      !for ZGPADM 
      INTEGER,PARAMETER      :: PADEDEG=6
      REAL(q),PARAMETER      :: T=1._q
      COMPLEX(q),ALLOCATABLE :: WSP(:)
      INTEGER                :: LWSP
      INTEGER,ALLOCATABLE    :: IPIV(:)
      INTEGER                :: IEXPA
      INTEGER                :: NS
      INTEGER                :: IERROR
     
      N=SIZE(A,1)
      M=SIZE(A,2)
      IF( N/=M ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I4)')'ERROR in MATRIX_EXPONENTIAL:'//&
            'A is not a square matrix',N,M
         STOP
      ENDIF 
      
      IF( N>NDIM ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I4)')'WARNING from MATRIX_EXPONENTIAL:'//&
            'input matrix A is not of size NDIM',N,NDIM
      ENDIF 
  
      ALLOCATE(R(NDIM,NDIM))
      CALL C8MAT_EXPM1(NDIM,A,R,IU0)
      MATRIX_EXPONENTIAL=R

      DEALLOCATE(R)
      RETURN
   ENDFUNCTION MATRIX_EXPONENTIAL

!***************************INNER_PRODUCT*****************************
!
! computes the bilinear form Re[Tr(X'*Y)] for two points X, Y on U(n). 
!
!*********************************************************************

   FUNCTION INNER_PRODUCT(X,Y,IU0)
      IMPLICIT NONE
      REAL(q)       :: INNER_PRODUCT
      COMPLEX(q)    :: X(:,:)
      COMPLEX(q)    :: Y(:,:)
      INTEGER       :: IU0
      !local 
      REAL(q)                :: R
      COMPLEX(q),ALLOCATABLE :: A(:,:)
      COMPLEX(q),ALLOCATABLE :: B(:,:)
      INTEGER                :: I, J, M, N

      I=SIZE(X,1)
      J=SIZE(X,2)
      M=SIZE(Y,1)
      N=SIZE(Y,2)

      IF ( I /= J ) THEN
         IF ( IU0>=0) WRITE(IU0,'(A,2I4)')'ERROR in INNER_PRODUCT:'//&
            'X is not a square matrix',I,J
         STOP
      ENDIF 
      IF ( M /= N ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I4)')'ERROR in INNER_PRODUCT:'//&
            'Y is not a square matrix',M,N
         STOP
      ENDIF 
      IF ( I /= N ) THEN
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I4)')'ERROR in INNER_PRODUCT:'//&
            'X and Y must be of same size',I,N
         STOP
      ENDIF 

      ALLOCATE(A(N,N),B(N,N))
      A=X
      B=zero
      CALL ZGEMM('C','N',N,N,N,one,A(1,1),N,Y(1,1),N,zero,B(1,1),N )
#ifdef debug
WRITE(*,*)'X' 
    DO J=1,NDIM
        WRITE(*,'(16F20.10)')REAL(X(J,1:NDIM),q)
     ENDDO 
     DO J=1,NDIM
        WRITE(*,'(16F20.10)')IMAG(X(J,1:NDIM))
     ENDDO 
WRITE(*,*)'Y'
     DO J=1,NDIM
        WRITE(*,'(16F20.10)')REAL(Y(J,1:NDIM),q)
     ENDDO 
     DO J=1,NDIM
        WRITE(*,'(16F20.10)')IMAG(Y(J,1:NDIM))
     ENDDO 

WRITE(*,*)'B'
     DO J=1,NDIM
        WRITE(*,'(16F20.10)')REAL(B(J,1:NDIM),q)
     ENDDO 
     DO J=1,NDIM
        WRITE(*,'(16F20.10)')IMAG(B(J,1:NDIM))
     ENDDO 
#endif
      R=0._q
      DO I = 1, N
         R = R + REAL(B(I,I),q)
      ENDDO

      R=R/2._q
      INNER_PRODUCT=R

      DEALLOCATE(A,B)
      RETURN
   ENDFUNCTION INNER_PRODUCT

!**************************LUNITCHECK**********************************
!
! computes the norm of a complex matrix W as the largest singular value
! of W*.W and check's if the input matrix W is unitary 
!
!*********************************************************************

   FUNCTION LUNITCHECK(W,IU0)
      IMPLICIT NONE
      LOGICAL                :: LUNITCHECK!the function value is real
      COMPLEX(q)             :: W(:,:)    !input matrix   
      INTEGER                :: IU0       !writing unit
      !local
      INTEGER                :: N,M      !dimension of matrix
      INTEGER                :: I, J     !loop variables 
      REAL(q)                :: R        !auxilary variable 
      COMPLEX(q),ALLOCATABLE :: A(:,:)   !auxilary matrix
      COMPLEX(q),ALLOCATABLE :: B(:,:)   !auxilary matrix

      !get the size of the matrix
      N=SIZE(W,1)
      M=SIZE(W,2) 
   
      !allocate auxillary matrices 
      ALLOCATE(A(M,N),B(M,M)) 
      A=W
      B=zero
      !compute the "square of the matix", i.e. W'*W in matlab language
      CALL ZGEMM('C','N',M,M,N,one,A(1,1),M,W(1,1),N,zero,B(1,1),M )
   
      !when checking unitarity, the identity of W'*W is subtracted
      DO I=1,M
         B(I,I)=B(I,I)-CMPLX(1._q,0._q,q)
      ENDDO 
   
      !the norm is the largest singular value       
#ifdef debug
      DO I=1,N
         WRITE(*,'(16F20.10)')REAL(B(I,1:N),q)
      ENDDO 
      WRITE(*,*)
      DO I=1,N
         WRITE(*,'(16F20.10)')IMAG(B(I,1:N))
      ENDDO 
#endif 
      
      R=MAX_SINGULAR_VALUE(B,IU0) 
      IF ( ABS(R)>TINY_U ) THEN
         LUNITCHECK=.FALSE.
      ELSE
         LUNITCHECK=.TRUE.
      ENDIF  
   
      DEALLOCATE(A,B)
      RETURN
   ENDFUNCTION LUNITCHECK

!****************************DELTA************************************
!
! a simple Kronecker delta function (real valued)
!
!*********************************************************************

   FUNCTION DELTA(I,J)
      IMPLICIT NONE
      REAL(q)  :: DELTA !is either 0 or 1 
      INTEGER  :: I,J   !arguments 
      IF ( I == J ) THEN
         DELTA=1._q
      ELSE
         DELTA=0._q
      ENDIF 
      RETURN
   ENDFUNCTION DELTA

!****************************NORM*************************************
!
! gives the P norm of the matrix A
!
!*********************************************************************

   FUNCTION NORM(W,P,IU0)
      IMPLICIT NONE
      REAL(q)                :: NORM     !function value
      INTEGER                :: P        !what kind of norm 
      COMPLEX(q)             :: W(:,:)   !input matrix   
      INTEGER                :: IU0      !writing unit
      !local
      INTEGER                :: N,M      !dimension of matrix
      INTEGER                :: I, J     !loop variables 
      REAL(q)                :: R        !auxilary variable 
      REAL(q),ALLOCATABLE    :: RP(:)    !auxilary vector
      COMPLEX(q),ALLOCATABLE :: A(:,:)   !auxilary matrix
      COMPLEX(q),ALLOCATABLE :: B(:,:)   !auxilary matrix

      !get the size of the matrix
      N=SIZE(W,1)
      M=SIZE(W,2) 
      ALLOCATE(A(N,M))   
      A=W

      IF ( P == 0 ) THEN       !p==0 is chosen to be the Frobenius norm of A
         DO I=1,N 
            DO J=1,M 
               R=R+ABS(A(I,J))**2
            ENDDO
         ENDDO
         NORM=SQRT(R)
      ELSEIF( P==1 ) THEN      !the p=1 norm is the maximum of the column sums
         R=0._q
         ALLOCATE(RP(N))
         RP=0._q
         DO I=1,N
            RP(I)=SUM(ABS(A(I,1:M)))
         ENDDO 
         R=MAXVAL(RP(1:N))
         NORM=R
         DEALLOCATE(RP)
      ELSE                     !the norm is the largest singular value       
         R=MAX_SINGULAR_VALUE(A,IU0) 
         NORM=R
      ENDIF 
  
      DEALLOCATE(A)
      RETURN
   ENDFUNCTION NORM

!***********************MAX_SINGULAR_VALUE****************************
!
! computes the largest singular value a complex matrix A
! (not sure if function is working for rectangular matrices)
!
!*********************************************************************
 
   FUNCTION MAX_SINGULAR_VALUE(A,IU0)
      REAL(q)                :: MAX_SINGULAR_VALUE
      COMPLEX(q)             :: A(:,:)
      INTEGER                :: IU0
      !local
      INTEGER                :: M, N
      INTEGER                :: LDA, LDU, LDVT
      INTEGER,PARAMETER      :: LWMAX=1000
      INTEGER                :: INFO, LWORK
      COMPLEX(q),ALLOCATABLE :: AA(:,:)
      REAL(q),ALLOCATABLE    :: SIGMA(:)
      REAL(q),ALLOCATABLE    :: RWORK(:)
      COMPLEX(q),ALLOCATABLE :: U(:,:)
      COMPLEX(q),ALLOCATABLE :: VT(:,:)
      COMPLEX(q),ALLOCATABLE :: WORK(:)

      !get the size of the matrix and the LDs
      M=SIZE(A,2)
      N=SIZE(A,1)
      LDA = M
      LDU = M
      LDVT = N

      !we work in AA
      ALLOCATE(AA(N,M))
      AA=A

      !RWORK dimension should be at least MAX( 1, 5*MIN(M,N) )
      ALLOCATE(SIGMA(MIN(M,N)),RWORK(MAX( 1, 5*MIN(M,N))))

      ALLOCATE(U(LDU,M),VT(LDVT,N),WORK(LWMAX))

      LWORK = -1
      CALL ZGESVD( 'N', 'N', M, N, AA, LDA, SIGMA, U, LDU, VT, LDVT,&
         WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      !Compute only singular values 
      CALL ZGESVD( 'N', 'N', M, N, AA, LDA, SIGMA, U, LDU, VT, LDVT,&
         WORK, LWORK, RWORK, INFO )

      !Check for convergence.
      IF( INFO>0 ) THEN
         IF(IU0>=0) &
            WRITE(IU0,*)'ERROR in MAX_SINGULAR_VALUE: ZGESVD not converged',INFO
         STOP
      ENDIF

      !the function value is the first singular value of A
      MAX_SINGULAR_VALUE=SIGMA(1)
      DEALLOCATE(SIGMA,RWORK,U,VT,WORK)
      DEALLOCATE(AA)
      RETURN
   ENDFUNCTION MAX_SINGULAR_VALUE

!***************************MAX_EIG***********************************
!
! computes the largest absolute eigenvalue a complex matrix C
! (not sure if function is working for rectangular matrices)
!
!*********************************************************************
 
   FUNCTION MAX_EIG(C,IU0)
      REAL(q)                :: MAX_EIG
      COMPLEX(q)             :: C(:,:)
      INTEGER                :: IU0
      !local
      INTEGER                :: M, N
      INTEGER                :: LDC
      INTEGER                :: LDVL
      INTEGER                :: LDVR
      INTEGER,PARAMETER      :: LWMAX=10000
      INTEGER                :: INFO, LWORK
      REAL(q),ALLOCATABLE    :: RWORK(:) 
      COMPLEX(q),ALLOCATABLE :: CC(:,:)            
      COMPLEX(q),ALLOCATABLE :: VL(:,:)            
      COMPLEX(q),ALLOCATABLE :: VR(:,:)            
      COMPLEX(q),ALLOCATABLE :: W(:)            
      COMPLEX(q),ALLOCATABLE :: WORK(:)            

      !get the size of the matrix and the LDs
      M=SIZE(C,2)
      N=SIZE(C,1)
      LDC = M
      LDVL = 1
      LDVR = 1

      !we work in CC
      ALLOCATE(CC(N,M))
      CC=C

      ALLOCATE(RWORK(2*LDC))
      ALLOCATE(VL(LDC,LDC),VR(LDC,LDC))
      ALLOCATE(W(LDC), WORK(LWMAX))

      ! Query the optimal workspace.
      LWORK = -1
      CALL ZGEEV( 'N', 'N', LDC, CC, LDC, W, VL, LDVL,&
         VR, LDVR, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      ! Solve eigenproblem.
      CALL ZGEEV( 'N', 'N', LDC, CC, LDC, W, VL, LDVL,&
         VR, LDVR, WORK, LWORK, RWORK, INFO )

      ! Check for convergence.
      IF( INFO > 0 ) THEN
         IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in MAX_EIG:'//&
            'ZGEEV returns',INFO
         STOP
      ENDIF

      !the function value is the first singular value of A
      MAX_EIG=MAXVAL(ABS(W))
      !get rid of auxilary arrays
      DEALLOCATE(CC)
      DEALLOCATE(RWORK,VL,VR,W,WORK)
      RETURN
   ENDFUNCTION MAX_EIG

!******************************HANNING*******************************
!
! Computes the Hanning window in a symmetric fashion
!
!*********************************************************************
 
   FUNCTION HANNING(NLENGTH,IU0)
      REAL(q)                :: HANNING(NLENGTH)
      INTEGER                :: NLENGTH
      INTEGER                :: IU0
      !local
      INTEGER                :: I
      REAL(q),ALLOCATABLE    :: R(:)

      IF ( NLENGTH <= 0 ) THEN
         IF ( IU0 >=0) WRITE(IU0,'(A,I6)')'ERROR in HANNING:'//&
            'NLENGHT is negative',NLENGTH
         STOP
      ENDIF 

      ALLOCATE(R(NLENGTH)) 
      
      DO I = 1, NLENGTH      
         R(I) = (1._q-COS(2*PI*I/REAL(NLENGTH+1,q)))/2._q
      ENDDO 

      HANNING(1:NLENGTH)=R(1:NLENGTH)
      DEALLOCATE(R)
      RETURN
   ENDFUNCTION HANNING

!*********************************************************************
!
!                 A U X I L A R Y  R O U T I N E S 
!
!*********************************************************************

!**************************R_LINEAR_SOLVE*****************************
!
! This solves a real linear system of equations 
! A ... coefficient matrix (left untouched)
! B ... r.h.s vector (contains solution on exit) 
!
!*********************************************************************
 
   SUBROUTINE R_LINEAR_SOLVE(A,B,IU0)
      REAL(q)                :: A(:,:)
      REAL(q)                :: B(:)
      INTEGER                :: IU0
      !local
      INTEGER                :: M, N
      INTEGER                :: LDA
      REAL(q),ALLOCATABLE    :: AA(:,:)
      REAL(q),ALLOCATABLE    :: BB(:)
      INTEGER                :: INFO
      INTEGER,ALLOCATABLE    :: IPIV(:)

      LDA=SIZE(A,2)
        N=SIZE(A,1)
 
      !some checks (we allow for square matrices A only!)
      IF ( SIZE(B) /= LDA ) THEN
         IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ERROR in R_LINEAR_SOLVE:'//&
            ' Sizes of B and A are inconsistent',SIZE(B),LDA
         STOP
      ENDIF 
      IF ( N /= LDA ) THEN
         IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ERROR in R_LINEAR_SOLVE:'//&
            ' A is not a square coefficient matrix',N,LDA
         STOP
      ENDIF 
 
      !allocate memory 
      ALLOCATE(AA(N,LDA),BB(LDA))  
      AA=A
      BB=B
 
      !also we need a pivot array IPIV
      ALLOCATE(IPIV(LDA))
  
      !now it's time to call LAPACK, first a LU decomposition of AA
      CALL DGETRF(LDA,LDA,AA(1,1),LDA,IPIV,INFO)
      IF (INFO /= 0 ) THEN
         IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in R_LINEAR_SOLVE:'//&
            ' DGETRF failed with output ',INFO
         STOP
      ENDIF 

      ! now solve the system with DGETRS
      CALL DGETRS( 'N', LDA, 1, AA( 1, 1 ), LDA, IPIV, BB(1), LDA, INFO)
      IF (INFO /= 0 ) THEN
         IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in R_LINEAR_SOLVE:'//&
            ' DGETRS failed with output ',INFO
         STOP
      ENDIF 
  
      !BB is the solution vector, save it to B
      B=BB
 
      !get rid of all auxilary arrays
      DEALLOCATE(AA,BB)
      DEALLOCATE(IPIV)
 
      RETURN
   END SUBROUTINE R_LINEAR_SOLVE

!**************************C_LINEAR_SOLVE*****************************
!
! This solves a complex linear system of equations 
! A ... coefficient matrix (left untouched)
! B ... r.h.s vector (contains solution on exit) 
!
!*********************************************************************
 
   SUBROUTINE C_LINEAR_SOLVE(A,B,IU0)
      COMPLEX(q)             :: A(:,:)
      COMPLEX(q)             :: B(:)
      INTEGER                :: IU0
      !local
      INTEGER                :: M, N
      INTEGER                :: LDA
      COMPLEX(q),ALLOCATABLE :: AA(:,:)
      COMPLEX(q),ALLOCATABLE :: BB(:)
      INTEGER                :: INFO
      INTEGER,ALLOCATABLE    :: IPIV(:)

      LDA=SIZE(A,2)
        N=SIZE(A,1)
 
      !some checks (we allow for square matrices A only!)
      IF ( SIZE(B) /= LDA ) THEN
         IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ERROR in C_LINEAR_SOLVE:'//&
            ' Sizes of B and A are inconsistent',SIZE(B),LDA
         STOP
      ENDIF 
      IF ( N /= LDA ) THEN
         IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ERROR in C_LINEAR_SOLVE:'//&
            ' A is not a square coefficient matrix',N,LDA
         STOP
      ENDIF 
 
      !allocate memory 
      ALLOCATE(AA(N,LDA),BB(LDA))  
      AA=A
      BB=B
 
      !also we need a pivot array IPIV
      ALLOCATE(IPIV(LDA))
  
      !now it's time to call LAPACK, first a LU decomposition of AA
      CALL ZGETRF(LDA,LDA,AA(1,1),LDA,IPIV,INFO)
      IF (INFO /= 0 ) THEN
         IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in C_LINEAR_SOLVE:'//&
            ' ZGETRF failed with output ',INFO
         STOP
      ENDIF 

      ! now solve the system with DGETRS
      CALL ZGETRS( 'N', LDA, 1, AA( 1, 1 ), LDA, IPIV, BB(1), LDA, INFO)
      IF (INFO /= 0 ) THEN
         IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in C_LINEAR_SOLVE:'//&
            ' ZGETRS failed with output ',INFO
         STOP
      ENDIF 
  
      !BB is the solution vector, save it to B
      B=BB
 
      !get rid of all auxilary arrays
      DEALLOCATE(AA,BB)
      DEALLOCATE(IPIV)
 
      RETURN
   END SUBROUTINE C_LINEAR_SOLVE

!**************************C_POLY_ROOTS*******************************
!
! This routine calculates all (complex) roots of the polynomial
! given by the coefficient vector A. The solutions are found by 
! calculating the eigenvalues of the companion matrix 
! 
! A ... coefficient vector (complex vector)
! B ... solution on exit 
! 
!*********************************************************************
 
   SUBROUTINE C_POLY_ROOTS(A,B,IU0)
      COMPLEX(q)             :: A(:)
      COMPLEX(q)             :: B(:)
      INTEGER                :: IU0
      !local
      INTEGER                :: I, J 
      INTEGER                :: LDA
      COMPLEX(q),ALLOCATABLE :: AA(:)
      COMPLEX(q),ALLOCATABLE :: BB(:)
      REAL(q),PARAMETER      :: VERYTINY=1.E-16_q
      INTEGER                :: NONZERO
      COMPLEX(q),ALLOCATABLE :: C(:,:)
      COMPLEX(q),ALLOCATABLE :: D(:)
      !for ZGEEV
      INTEGER                :: LDC
      INTEGER                :: LDVL
      INTEGER                :: LDVR
      INTEGER,PARAMETER      :: LWMAX=10000
      INTEGER                :: INFO, LWORK
      REAL(q),ALLOCATABLE    :: RWORK(:) 
      COMPLEX(q),ALLOCATABLE :: CC(:,:)            
      COMPLEX(q),ALLOCATABLE :: VL(:,:)            
      COMPLEX(q),ALLOCATABLE :: VR(:,:)            
      COMPLEX(q),ALLOCATABLE :: W(:)            
      COMPLEX(q),ALLOCATABLE :: WORK(:)            

      LDA=SIZE(A)
 
      ! there are SIZE(A)-1 complex roots for a polynomial
      ! of degree SIZE(A)-1
      IF ( SIZE(B) /= LDA-1 ) THEN
         IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ERROR in C_POLY_ROOTS:'//&
            ' Sizes of A and B are inconsistent',LDA-1,SIZE(B)
         STOP
      ENDIF 
 
      !allocate memory 
      ALLOCATE(AA(LDA),BB(LDA))  
      AA=A
      BB=B
      
      !first find all non-zero (maximally LDA) elements of AA
      NONZERO=0
      DO I = 1, LDA
         IF (ABS(AA(I)) > VERYTINY ) THEN
            NONZERO=NONZERO+1
         ENDIF 
      ENDDO
 
      !do this if AA contains non-zero entries 
      IF ( NONZERO > 0  ) THEN
         LDC=LDA-1
         !allocate memory for companion matrix 
         ALLOCATE(C(LDC,LDC))
         C=zero
         DO I = 1, LDC-1
            C(I+1,I)=one
         ENDDO
         !d-vector is -A/A(1) and is the upper most vector of C
         ALLOCATE(D(LDC))
         D=zero
         D(1:LDC) = -AA(2:LDA)/AA(1)
         DO I=1,LDC
            C(1,I)=D(I)
         ENDDO
#ifdef debug
WRITE(*,*)' companion' 
DO I=1,LDC
  WRITE(*,'(16E12.4)')C(I,1:LDC)
ENDDO
#endif 
         !call LAPACK for Eigenvalues of C
         !they can be complex!!
         ALLOCATE(CC(LDC,LDC))
         DO I = 1, LDC
            DO J = 1, LDC
               CC(I,J)=C(I,J)
            ENDDO
         ENDDO
 
         !working arrays
         LDVL=1
         LDVR=1
         ALLOCATE(RWORK(2*LDC))
         ALLOCATE(VL(LDC,LDC),VR(LDC,LDC))
         ALLOCATE(W(LDC), WORK(LWMAX))

         ! Query the optimal workspace.
         LWORK = -1
         CALL ZGEEV( 'N', 'N', LDC, CC, LDC, W, VL, LDVL,&
            VR, LDVR, WORK, LWORK, RWORK, INFO )
         LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

         ! Solve eigenproblem.
 
         CALL ZGEEV( 'N', 'N', LDC, CC, LDC, W, VL, LDVL,&
            VR, LDVR, WORK, LWORK, RWORK, INFO )

         ! Check for convergence.
         IF( INFO > 0 ) THEN
            IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in C_POLY_ROOTS:'//&
               'ZGEEV returns',INFO
            STOP
         ENDIF

#ifdef debug
WRITE(*,*)' eigenvalues ' 
DO I=1,LDC
  WRITE(*,'(2E12.4)')W(I)
ENDDO
#endif 
         !the complex eigenvalues are the roots of the problem
         BB=W
         !get rid of auxilary arrays
         DEALLOCATE(C)
         DEALLOCATE(D)
         DEALLOCATE(CC)
         DEALLOCATE(RWORK,VL,VR,W,WORK)
      ELSE   ! if all elements are zero the solution is trivial
         BB=zero 
      ENDIF

      B=BB     
 
      !get rid of auxillary arrays
      DEALLOCATE(AA,BB)
 
      RETURN
   END SUBROUTINE C_POLY_ROOTS

!**************************R_POLY_ROOTS*******************************
!
! This routine calculates all (real) roots of the polynomial
! given by the coefficient vector A. The solutions are found by 
! calculating the eigenvalues of the companion matrix 
! 
! A      ... coefficient vector (real vector)
! NROOTS ... contains the number of all real roots on exit
! B      ... real roots  
! 
!*********************************************************************
 
   SUBROUTINE R_POLY_ROOTS(A,B,NROOTS,IU0)
      REAL(q) :: A(:) 
      REAL(q) :: B(:)
      INTEGER :: NROOTS 
      INTEGER :: IU0
      !local
      INTEGER                :: I, J 
      INTEGER                :: LDA
      REAL(q),ALLOCATABLE    :: AA(:)
      REAL(q),ALLOCATABLE    :: BB(:)
!      REAL(q),PARAMETER      :: VERYTINY=2.2204E-16_q !(floating point prec)
      REAL(q),PARAMETER      :: VERYTINY=2.2204E-12_q !(larger than theoretical prec)
      INTEGER                :: NONZERO
      REAL(q),ALLOCATABLE    :: C(:,:)
      REAL(q),ALLOCATABLE    :: D(:)
      !for ZGEEV
      INTEGER                :: LDC
      INTEGER                :: LDVL
      INTEGER                :: LDVR
      INTEGER,PARAMETER      :: LWMAX=10000
      INTEGER                :: INFO, LWORK
      REAL(q),ALLOCATABLE    :: RWORK(:) 
      COMPLEX(q),ALLOCATABLE :: CC(:,:)            
      COMPLEX(q),ALLOCATABLE :: VL(:,:)            
      COMPLEX(q),ALLOCATABLE :: VR(:,:)            
      COMPLEX(q),ALLOCATABLE :: W(:)            
      COMPLEX(q),ALLOCATABLE :: WORK(:)            
      INTEGER                :: NZ

      LDA=SIZE(A)
 
      ! there are SIZE(A)-1 complex roots for a polynomial
      ! of degree SIZE(A)-1
      IF ( SIZE(B) /= LDA-1 ) THEN
         IF (IU0>=0) WRITE(IU0,'(A,2I5)')'ERROR in R_POLY_ROOTS:'//&
            ' Sizes of B and A are inconsistent',SIZE(B),LDA-1
         STOP
      ENDIF 
 
      !allocate memory 
      ALLOCATE(AA(LDA),BB(LDA))  
      AA=A
      BB=0._q
      
      !first find all non-zero (maximally LDA) elements of AA
      NONZERO=0
      DO I = 1, LDA
         IF (ABS(AA(I)) > VERYTINY ) THEN
            NONZERO=NONZERO+1
         ENDIF 
      ENDDO
 
      !do this if AA contains non-zero entries 
      IF ( NONZERO > 0  ) THEN
         LDC=LDA-1
         !allocate memory for companion matrix 
         ALLOCATE(C(LDC,LDC))
         C=0._q
         DO I = 1, LDC-1
            C(I+1,I)=one
         ENDDO
         !d-vector is -A/A(1) and is the upper most vector of C
         ALLOCATE(D(LDC))
         D=0._q
         D(1:LDC) = -AA(2:LDA)/AA(1)
         DO I=1,LDC
            C(1,I)=D(I)
         ENDDO
#ifdef debug
WRITE(*,*)' companion' 
DO I=1,LDC
  WRITE(*,'(16E12.4)')C(I,1:LDC)
ENDDO
#endif 

         !call LAPACK for Eigenvalues of C
         !they can be complex!!
         ALLOCATE(CC(LDC,LDC))
         DO I = 1, LDC
            DO J = 1, LDC
               CC(I,J)=CMPLX(C(I,J),0._q,q)
            ENDDO
         ENDDO
 
         !working arrays
         LDVL=1
         LDVR=1
         ALLOCATE(RWORK(2*LDC))
         ALLOCATE(VL(LDC,LDC),VR(LDC,LDC))
         ALLOCATE(W(LDC), WORK(LWMAX))

         ! Query the optimal workspace.
         LWORK = -1
         CALL ZGEEV( 'N', 'N', LDC, CC, LDC, W, VL, LDVL,&
            VR, LDVR, WORK, LWORK, RWORK, INFO )
         LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

         ! Solve eigenproblem.
 
         CALL ZGEEV( 'N', 'N', LDC, CC, LDC, W, VL, LDVL,&
            VR, LDVR, WORK, LWORK, RWORK, INFO )

         ! Check for convergence.
         IF( INFO > 0 ) THEN
            IF(IU0>=0) WRITE(IU0,'(A,I5)')'ERROR in R_POLY_ROOTS:'//&
               'ZGEEV returns',INFO
            STOP
         ENDIF

#ifdef debug
WRITE(*,*)' eigenvalues ' 
DO I=1,LDC
  WRITE(*,'(2E12.4)')W(I)
ENDDO
#endif 
         !get rid of auxilary arrays
         DEALLOCATE(C)
         DEALLOCATE(D)
         DEALLOCATE(CC)
         DEALLOCATE(RWORK,VL,VR,WORK)

         !we need to find all real roots 
         NZ = 0 
         DO I = 1, LDC
            IF ( ABS(IMAG(W(I))) < VERYTINY ) THEN
               NZ = NZ + 1 
               BB(NZ) = REAL(W(I),q)
            ENDIF 
         ENDDO
         DEALLOCATE(W)
      ELSE   ! if all elements are zero the solution is trivial
         BB=zero 
         !in this case we only have one trivial root at 0 
         NZ=1
      ENDIF

      B=BB   
      NROOTS=NZ  
 
      !get rid of auxillary arrays
      DEALLOCATE(AA,BB)
 
      RETURN
   END SUBROUTINE R_POLY_ROOTS

!****************************R_FFTSHIFT*******************************
!
! Performs a FFT of the real vector V and shifts the zero-frequency 
! component to the center of the Fourier spectrum.
! FV contains the (complex) transformed vector.
! (If V has length 1 the FFT is not performed)
!
!*********************************************************************

   SUBROUTINE R_FFTSHIFT(V,FV,IU0)
      IMPLICIT NONE 
      include 'fftw3.f'
      REAL(q)     :: V(:)      !vector to be Fourier transformed
      COMPLEX(q)  :: FV(:)     !transformed vector 
      INTEGER     :: IU0       !for writing messages to the user 
      !local
      INTEGER                :: N         !length of the vector 
      INTEGER(8)             :: PLAN      !FFTW plan
      COMPLEX(q),ALLOCATABLE :: DATAIN(:) !auxilary 
      COMPLEX(q),ALLOCATABLE :: DATAOUT(:)!auxilary 
      INTEGER                :: I         !loop 
      INTEGER                :: MID       !aux
      
      !determine size of V
      N=SIZE(V)
      !maybe FV array is too small 
      IF ( SIZE(FV) < N ) THEN 
         IF (IU0>=0) WRITE(IU0,'(A,2I6)')'ERROR in R_FFTSHIFT:'//&
            ' FV is too small',SIZE(FV),N
         STOP
      ENDIF 
 
      IF ( N > 1 ) THEN   
         !allocate some auxilary arrays
         ALLOCATE(DATAIN(N),DATAOUT(N))
         !save input data to DATAIN
         DO I = 1, N
            DATAIN(I)=CMPLX(V(I),0._q,q)
         ENDDO

         !create plan
         call dfftw_plan_dft_1d( PLAN , N , DATAIN , DATAOUT, &
            FFTW_FORWARD,FFTW_ESTIMATE )
 
         !a value 0 for PLAN is bad
         IF ( PLAN == 0 ) THEN
            IF ( IU0>=0) WRITE(IU0,'(A,3I6)')'ERROR in R_FFTSHIFT:'//&   
               ' Creating a FFT plan failed',PLAN
            STOP
         ENDIF 

         !perform fft 
         call dfftw_execute_dft( PLAN, DATAIN, DATAOUT )
         
         !destroy plans 
         call dfftw_destroy_plan( PLAN ) 

         !largest frequency is at DATAOUT(1)
         !we want it in the center MID=N/2
         MID=N/2+1
         DO I = 1, MID-1
            DATAIN(I)=DATAOUT(MID+I)
            DATAIN(MID+I-1)=DATAOUT(I)
         ENDDO
         DATAIN(N)=DATAOUT(MID)
         !save result to FV (divided by N)
         FV(1:N)=DATAIN(1:N)/REAL(N,q)
         DEALLOCATE(DATAIN,DATAOUT)
      ELSEIF ( N==1 ) THEN
         !a single number is not transformed here
         FV(1)=CMPLX(V(1),0._q,q)
      ELSE  !this will probably not happen
         IF ( IU0>=0 ) WRITE(IU0,'(A,2I6)')'ERROR in R_FFTSHIFT:'//&
            ' Size of V is ',N
         STOP
      ENDIF 
  

      RETURN
   END SUBROUTINE R_FFTSHIFT

!******************************SORT***********************************
!
! Sorts array C
!
!*********************************************************************

   SUBROUTINE SORT(N,C,IU0)
      REAL(q) :: A,C(N)
      REAL(q) :: ARR(N) 
      INTEGER :: IU0
      INTEGER :: I, N, J

      IF (N>100) THEN 
         IF(IU0>=0) WRITE(IU0,'(A,I6)')'ERROR in SORT:'//&
            ' Array too large',N
      ENDIF
      ARR(1:N)=C(1:N)
      DO J=2, N
            A=ARR(J)
         DO I=J-1,1,-1
            IF (ARR(I)<=A) GOTO 10
            ARR(I+1)=ARR(I)
         ENDDO
         I=0
10       ARR(I+1)=A
      ENDDO
      C(1:N)=ARR(1:N)
      RETURN
   END SUBROUTINE SORT

!*********************************************************************
!
! the following routines are for the matrix exponential
!
!*********************************************************************

!****************************c8mat_expm1******************************
!
! C8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential 
!  algorithm.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Author:
!    FORTRAN90 version by John Burkardt
!  Reference:
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!  Parameters:
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!    Input, complex ( kind = 8 ) A(N,N), the matrix.
!    Output, complex ( kind = 8 ) E(N,N), the estimate for exp(A).
!*********************************************************************

   SUBROUTINE C8MAT_EXPM1 ( N, A, E, IU0)
      IMPLICIT NONE
      INTEGER           :: N
      COMPLEX(q)        :: A(N,N)
      COMPLEX(q)        :: E(N,N)
      INTEGER           :: IU0
      !local 
      COMPLEX(q)        :: A2(N,N)
      REAL(q)           :: A_NORM
      REAL(q)           :: C
      COMPLEX(q)        :: D(N,N)
      INTEGER           :: EE
      INTEGER           :: K
      LOGICAL           :: P
      COMPLEX(q)        :: X(N,N)
      INTEGER           :: I
      INTEGER,PARAMETER :: Q = 6
      INTEGER           :: S
    !
    !  make a copy of the matrix.
    !
      A2(1:N,1:N) = A(1:N,1:N)
    !
    !  compute the l-infinity norm.
    !
      A_NORM = NORM(A2,1,IU0)
    !
    !  determine a scaling factor for the matrix.
    !
      EE = INT ( LOG_2 ( A_NORM ) ) + 1
      S = MAX ( 0, EE + 1 )
      A2(1:N,1:N) = A2(1:N,1:N) / 2.d+0**S
      X(1:N,1:N) = A2(1:N,1:N)
    
      C = 0.5d+0
    
      E = 0.d+0
      DO I = 1, N
         E(I,I)=1.0d+0 
      ENDDO
     
      E(1:N,1:N) = E(1:N,1:N) + C * A2(1:N,1:N)
    
      D = 0.d+0
      DO I = 1, N
         D(I,I)=1.0d+0 
      ENDDO
    
      D(1:N,1:N) = D(1:N,1:N) - C * A2(1:N,1:N)
    
      P = .TRUE.
    
      DO K = 2, Q
    
        C = C*REAL( Q - K + 1, kind=8 )/REAL ( K * ( 2 * Q - K + 1 ), kind=8 )
    
        X(1:N,1:N) = MATMUL( A2(1:N,1:N), X(1:N,1:N) )
    
        E(1:N,1:N) = E(1:N,1:N) + C * X(1:N,1:N)
    
        IF ( P ) THEN
          D(1:N,1:N) = D(1:N,1:N) + C * X(1:N,1:N)
        ELSE
          D(1:N,1:N) = D(1:N,1:N) - C * X(1:N,1:N)
        END IF
        P = .NOT. P
      END DO
    !
    !  E -> INVERSE(D) * E
    !
      CALL C8MAT_MINVM ( N, N, D, E, E )
    !
    !  E -> E^(2*S)
    !
      DO K = 1, S
        E(1:N,1:N) = MATMUL( E(1:N,1:N), E(1:N,1:N) )
      END DO
      RETURN
   END SUBROUTINE C8MAT_EXPM1

!***************************C8MAT_MINVM*******************************
!
! C8MAT_MINVM computes inverse(A) * B for C8MAT's.
!  Discussion:
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) N1, N2, the order of the matrices.
!    Input, complex ( kind = 8 ) A(N1,N1), B(N1,N2), the matrices.
!    Output, complex ( kind = 8 ) C(N1,N2), the result, 
!            C = inverse(A) * B.
!*********************************************************************

   SUBROUTINE C8MAT_MINVM( N1, N2, A, B, C )
     IMPLICIT NONE
     INTEGER ( KIND = 4 ) N1
     INTEGER ( KIND = 4 ) N2
     COMPLEX (q) A(N1,N1)
     COMPLEX (q) ALU(N1,N1)
     COMPLEX (q) B(N1,N2)
     COMPLEX (q) C(N1,N2)
     INTEGER ( KIND = 4 ) INFO
   
     ALU(1:N1,1:N1) = A(1:N1,1:N1)
     C(1:N1,1:N2) = B(1:N1,1:N2)
   
     CALL C8MAT_FSS ( N1, ALU, N2, C, INFO )
     IF ( INFO /= 0 ) THEN
       WRITE ( *, '(A)' ) ' '
       WRITE ( *, '(A)' ) 'C8MAT_MINVM - FATAL ERROR!'
       WRITE ( *, '(A)' ) '  THE MATRIX A WAS NUMERICALLY SINGULAR.'
       STOP
     END IF
     RETURN
   END SUBROUTINE C8MAT_MINVM

!******************************C8MAT_FSS******************************
!
!  C8MAT_FSS factors and solves a system with multiple right
!  hand sides.
!  Discussion:
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!    This routine does not save the LU factors of the matrix, 
!    and hence cannot be used to efficiently solve multiple 
!    linear systems, or even to factor A at one time, and solve 
!    a single linear system at a later time. This routine uses partial 
!    pivoting, but no pivot vector is required.
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!  Author:
!    John Burkardt
!  Parameters:

!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, complex ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input/output, complex ( kind = 8 ) B(N,NB).
!    On input, the right hand sides of the linear system.
!    On output, the solutions of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
!*********************************************************************

   SUBROUTINE C8MAT_FSS( N, A, NB, B, INFO )
      IMPLICIT NONE
      INTEGER (KIND = 4) N
      INTEGER (KIND = 4) NB
      COMPLEX (q) A(N,N)
      COMPLEX (q) B(N,NB)
      INTEGER ( KIND = 4 ) I
      INTEGER ( KIND = 4 ) INFO
      INTEGER ( KIND = 4 ) IPIV
      INTEGER ( KIND = 4 ) J
      INTEGER ( KIND = 4 ) JCOL
      REAL (q) PIV
      COMPLEX (q) ROW(N)
      COMPLEX (q) T(NB)
      COMPLEX (q) TEMP
    
      INFO = 0
    
      DO JCOL = 1, N
    !
    !  find the maximum element in column i.
    !
        PIV = ABS ( A(JCOL,JCOL) )
        IPIV = JCOL
        DO I = JCOL + 1, N
          IF ( PIV < ABS ( A(I,JCOL) ) ) THEN
            PIV = ABS ( A(I,JCOL) )
            IPIV = I
          END IF
        END DO
    
        IF ( PIV == 0._q ) THEN
          INFO = JCOL
          WRITE ( *, '(A)' ) ' '
          WRITE ( *, '(A)' ) 'C8MAT_FSS - FATAL ERROR!'
          WRITE ( *, '(A,I8)' ) '  ZERO PIVOT ON STEP ', INFO
          STOP
        END IF
    !
    !  switch rows jcol and ipiv, and b.
    !
        IF ( JCOL /= IPIV ) THEN
    
          ROW(1:N) = A(JCOL,1:N)
          A(JCOL,1:N) = A(IPIV,1:N)
          A(IPIV,1:N) = ROW(1:N)
    
          T(1:NB)      = B(JCOL,1:NB)
          B(JCOL,1:NB) = B(IPIV,1:NB)
          B(IPIV,1:NB) = T(1:NB)
    
        END IF
    !
    !  scale the pivot row.
    !
        A(JCOL,JCOL+1:N) = A(JCOL,JCOL+1:N) / A(JCOL,JCOL)
        B(JCOL,1:NB) = B(JCOL,1:NB) / A(JCOL,JCOL)
        A(JCOL,JCOL) = 1._q
    !
    !  use the pivot row to eliminate lower entries in that column.
    !
        DO I = JCOL + 1, N
          IF ( A(I,JCOL) /= 0._q ) THEN
            TEMP = - A(I,JCOL)
            A(I,JCOL) = 0._q
            A(I,JCOL+1:N) = A(I,JCOL+1:N) + TEMP * A(JCOL,JCOL+1:N)
            B(I,1:NB) = B(I,1:NB) + TEMP * B(JCOL,1:NB)
          END IF
        END DO
    
      END DO
    !
    !  back solve.
    !
      DO J = 1, NB
        DO JCOL = N, 2, -1
          B(1:JCOL-1,J) = B(1:JCOL-1,J) - A(1:JCOL-1,JCOL) * B(JCOL,J)
        END DO
      END DO
    
      RETURN
   END SUBROUTINE C8MAT_FSS

!*****************************LOG_2***********************************
!  LOG_2 returns the logarithm base 2 of an double precision real.
!  Discussion:
!    value = Log ( |X| ) / Log ( 2.0 )
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Author:
!    John Burkardt
!  Parameters:
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm 
!    is desired. X should not be 0. 
!    Output, real ( kind = 8 ) LOG_2, the logarithm base 2 
!    of the absolute value of X.  
!    It should be true that |X| = 2^LOG_2.
!*********************************************************************

  FUNCTION LOG_2 ( X )
     IMPLICIT NONE
     REAL(q) LOG_2
     REAL(q) X
   
     IF ( X == 0._q ) THEN
        LOG_2 = - HUGE ( X )
     ELSE
        LOG_2 = LOG ( ABS ( X ) ) / LOG ( 2._q )
     END IF
   
     RETURN
   ENDFUNCTION LOG_2

!*****************************RHERM***********************************
!
! Generates a random hermitian matrix of size N
!
!*********************************************************************

   FUNCTION RHERM(N,IU0)
      IMPLICIT NONE
      COMPLEX(q)  :: RHERM(N,N) !value of cost function for W
      INTEGER     :: N          !size of matrix 
      INTEGER     :: IU0        !for writing in parallel
      !local
      COMPLEX(q)  :: R(N,N) 
      INTEGER     :: I,J
      INTEGER     :: ISEED 
      REAL(q)     :: A,B
      REAL(q)     :: XR(N*N),XI(N*N)
      REAL(q)     :: ISTART,IFIN

      !the interval is hard coded to [-1,1]
      ISTART=-1._q
      IFIN=1._q
 
      !initialize seed 
      ISEED=1996
      CALL RANDOM_SEED(SIZE=ISEED) 
      !we need 2*N*N random numbers in total 
      DO I = 1, N*N
         CALL RANDOM_NUMBER(HARVEST=A)
         CALL RANDOM_NUMBER(HARVEST=B)
         !scale the numbers to interval (ISTART,IFIN)
         XR(I) = REAL(2*A-1,q)
         XI(I) = REAl(2*B-1,q)
      ENDDO  
      
      DO I =1, N
         DO J =1, N 
            IF ( I == J ) THEN  !fill diagonal with real numbers only 
               R(I,J) = CMPLX(XR(I+N*(J-1)),0._q,q)
            ELSE
               R(I,J) = CMPLX(XR(I+N*(J-1)),XI(I+N*(J-1)),q)
            ENDIF 
         ENDDO
      ENDDO 

      RHERM=(R+CONJG(TRANSPOSE(R)))/2._q
      RETURN
   ENDFUNCTION RHERM  

!*****************************RUNIT***********************************
!
! Generates a random unitary matrix from an random hermitian one.
!
!*********************************************************************

   FUNCTION RUNIT(N,IU0)
      IMPLICIT NONE
      COMPLEX(q)  :: RUNIT(N,N) !value of cost function for W
      INTEGER     :: N          !size of matrix 
      INTEGER     :: IU0        !for writing in parallel
      !local
      COMPLEX(q)  :: R(N,N) 
      INTEGER     :: I,J
      !for ZHEEV
      INTEGER,PARAMETER :: LWMAX=1000
      INTEGER           :: INFO, LWORK
      REAL(q)           :: W( N )
      REAL(q)           :: RWORK( 3*N-2 )
      COMPLEX(q)        :: WORK( LWMAX )


      !first generate a hermitian matrix 
      R=RHERM(N,IU0)    

      !we obtain the unitary matrix from matrix diagonalization
      !Query the optimal workspace.
      LWORK = -1
      CALL ZHEEV( 'Vectors', 'Lower', N, R, N, W, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      !Solve eigenproblem.
      CALL ZHEEV( 'Vectors', 'Lower', N, R, N, W, WORK, LWORK, RWORK, INFO )
      !Check for convergence.
      IF( INFO>0 ) THEN
         WRITE(IU0,'(A,I6)')'ERROR in RUNIT:'//&
            ' ZHEEV reports',INFO
         STOP
      ENDIF
      !otherwise R is the unitary matrix 
      RUNIT=R
      RETURN
   ENDFUNCTION RUNIT 

!*****************************RUNIT_DIA*******************************
!
! Generates a random diagonal unitary matrix.
!
!*********************************************************************

   FUNCTION RUNIT_DIAG(N,IU0)
      IMPLICIT NONE
      COMPLEX(q)  :: RUNIT_DIAG(N,N) !value of cost function for W
      INTEGER     :: N               !size of matrix 
      INTEGER     :: IU0             !for writing in parallel
      !local
      COMPLEX(q)  :: R(N,N) 
      INTEGER     :: I,J
      INTEGER     :: ISEED 
      REAL(q)     :: A,B
      REAL(q)     :: ISTART,IFIN
      REAL(q),PARAMETER :: VERYTINY=1.E-16

      !the interval is hard coded to [-1,1]
      ISTART=-1._q
      IFIN=1._q
      !initialize seed 
      ISEED=1996
      CALL RANDOM_SEED(SIZE=ISEED) 

      !clean R
      R=zero 
      !we need N random numbers in total 
      DO I = 1, N
         CALL RANDOM_NUMBER(HARVEST=A)
         CALL RANDOM_NUMBER(HARVEST=B)
         !scale the numbers to interval (ISTART,IFIN)
         R(I,I)=CMPLX(REAL(2*A-1,q),REAl(2*B-1,q),q)
         IF ( ABS(R(I,I)) > VERYTINY ) THEN
            R(I,I)=R(I,I)/ABS(R(I,I))
         ENDIF 
      ENDDO  

      !otherwise R is the unitary matrix 
      RUNIT_DIAG=R
      RETURN
   ENDFUNCTION RUNIT_DIAG

ENDMODULE lie
