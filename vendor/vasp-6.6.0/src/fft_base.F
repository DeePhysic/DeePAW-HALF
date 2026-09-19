#include "symbol.inc"
#define inversion_through_mid_point

!************************* SUBROUTINE FFTCHK ***************************
!
!> Returns the next correct setting for the three dimensional FFT
!
!***********************************************************************

      SUBROUTINE FFTCHK(NFFT)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION NFFT(3)
      LOGICAL FFTCH1

      DO 100 IND=1,3
  200 CONTINUE
        IF (FFTCH1(NFFT(IND))) GOTO 100
        NFFT(IND)=NFFT(IND)+1
        GOTO 200
  100 CONTINUE
      END SUBROUTINE FFTCHK


      LOGICAL FUNCTION FFTCH1(NIN)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (NFACT=4)
      DIMENSION IFACT(NFACT),NCOUNT(NFACT)
      DATA      IFACT /2,3,5,7/
      N=NIN
      DO 100 I=1,NFACT
         NCOUNT(I)=0
120      NEXT=N/IFACT(I)
         IF (NEXT*IFACT(I)==N) THEN
            N=NEXT
            NCOUNT(I)=NCOUNT(I)+1
            GOTO 120
         ENDIF
100   ENDDO
      IF (N==1 .AND. (NCOUNT(1)/=0)) &
           &  THEN
         FFTCH1=.TRUE.
      ELSE
         FFTCH1=.FALSE.
      ENDIF
      RETURN
      END FUNCTION FFTCH1


!************************* SUBROUTINE FFTGRIDPLAN **********************
!
!***********************************************************************

    SUBROUTINE FFTGRIDPLAN(GRID)
      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE (grid_3d) GRID

      ! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
      INTEGER, PARAMETER :: NADD=1024, NTIMES=4
      INTEGER :: I
#ifdef MPI
      IF (GRID%RL%NFAST/=1 .OR. GRID%RL_FFT%NFAST/=1) THEN
         CALL FFTGRIDPLAN_MPI(GRID)
         RETURN
      ENDIF
#endif
      ALLOCATE(CWORK(GRID%MPLWV+NADD))

      DO I=1,NTIMES
         CALL INIDAT(GRID%MPLWV+NADD,CWORK)
         CALL FFTMAKEPLAN(CWORK(I),GRID)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
!$       CALL FFTMAKEPLAN(CWORK(I),GRID)
!$OMP END MASTER
!$OMP END PARALLEL
      ENDDO

      DEALLOCATE(CWORK)

      RETURN
    END SUBROUTINE FFTGRIDPLAN


!************************* SUBROUTINE FFTINI ***************************
!
!>  If necessary this routine performes initialization
!>  for #FFTWAV and #FFTEXT
!>
!>  Usually this is only necessary for the Gamma-only
!>  single k-point version
!>
!>   FFTSCA(.,1) is the scaling factor for extracting the wavefunction
!>               from the FFT grid (#FFTEXT)
!>
!>   FFTSCA(.,2) is the scaling factor for puting the wavefunction on
!>               the grid
!
!***********************************************************************

    SUBROUTINE FFTINI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d)  GRID
      INTEGER :: NRPLW               !< maximum number of plane wave coefficients
      INTEGER :: NKPTS               !< number of k-points
      INTEGER :: NINDPW(NRPLW,NKPTS) !< index array from coefficients in sphere to 3D grid
      INTEGER :: NPLWKP(NKPTS)       !< actual number of plane wave coefficients for each k-points
    ! local
      INTEGER :: NK, NX, NY, NZ, NPL, N, IND, N1, N2, N3, NC
      INTEGER :: N2INV, N3INV, NUMBER_OF_CONJG, NCP
      REAL(q) :: FACTM

      IF (GRID%REAL2CPLX) THEN
#ifdef MPI
         IF (GRID%RL%NFAST/=1) THEN
            CALL FFTINI_MPI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
            RETURN
         ENDIF
#endif
         IF (NKPTS>1) THEN
            CALL vtutor%error("FFT3D: real version works only for 1 k-point")
         ENDIF

         NK=1
         NX=GRID%NGPTAR(1)
         NY=GRID%NGPTAR(2)
         NZ=GRID%NGPTAR(3)
         NPL=NPLWKP(NK)
         NULLIFY(GRID%FFTSCA)
         ALLOCATE(GRID%FFTSCA(NPL,2))
#ifdef inversion_through_mid_point
         NUMBER_OF_CONJG=0
         DO N=1,NPL
            IND=NINDPW(N,NKPTS)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)
            IF (N1==1 .AND. N/=1) THEN
               NUMBER_OF_CONJG=NUMBER_OF_CONJG+1
               ! invert second and third index
               N2INV=MOD(-GRID%LPCTY(N2)+GRID%NGY,GRID%NGY)+1
               N3INV=MOD(-GRID%LPCTZ(N3)+GRID%NGZ,GRID%NGZ)+1
               ! now determine the corresponding column index
               NCP=N2INV+(N3INV-1)*GRID%NGY
               ! check whether correct
               IF ( GRID%LPCTY(N2)+GRID%LPCTY(GRID%RC%I2(NCP)) /=0 .OR.  GRID%LPCTZ(N3)+GRID%LPCTZ(GRID%RC%I3(NCP))/=0) THEN
                  CALL vtutor%bug("internal error in FFTINI: could not determine conjugated &
                     &coefficient " // str(GRID%LPCTY(N2)) // " " // str(GRID%LPCTY(GRID%RC%I2(NCP))) &
                     // " " // str(GRID%LPCTZ(N3) + GRID%LPCTZ(GRID%RC%I3(NCP))), __FILE__, __LINE__)
               ENDIF
            ENDIF
         ENDDO
         NULLIFY(GRID%NINDPWCONJG, GRID%IND_IN_SPHERE)
         ! allocate the required index arrays
         ALLOCATE(GRID%NINDPWCONJG(NUMBER_OF_CONJG), GRID%IND_IN_SPHERE(NUMBER_OF_CONJG))
#endif
         NUMBER_OF_CONJG=0
         DO N=1,NPL
            IND=NINDPW(N,NKPTS)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)
            FACTM=SQRT(2._q)
            IF (N==1) FACTM=1
            ! in the real version the coefficients stored in the compressed mode (sphere)
            ! are multiplied by a factor sqrt(2) compared to complex version (except at Gamma)
            ! this allows to use DGEMM calls to calculate inproducts
            GRID%FFTSCA(N,1)= FACTM
            GRID%FFTSCA(N,2)= 1/FACTM
#ifdef inversion_through_mid_point
            IF (N1==1 .AND. N/=1) THEN
               NUMBER_OF_CONJG=NUMBER_OF_CONJG+1
               ! invert second and third index
               N2INV=MOD(-GRID%LPCTY(N2)+GRID%NGY,GRID%NGY)+1
               N3INV=MOD(-GRID%LPCTZ(N3)+GRID%NGZ,GRID%NGZ)+1
               ! now determine the corresponding column index
               NCP=N2INV+(N3INV-1)*GRID%NGY
               ! store original index in the compressed storage mode (sphere)
               GRID%IND_IN_SPHERE(NUMBER_OF_CONJG)=N
               ! store index to the 3D grid
               GRID%NINDPWCONJG(NUMBER_OF_CONJG)=N1+(NCP-1)*GRID%RC%NROW
            ENDIF
#else
            ! in the simple version, for the index nx=0
            ! only half of the FFT components are set on the 3D grid by FFTWAV
            ! however the complex to real FFT's (reciprocal space to real space)
            ! expects that the entire x=0 plane is set
            ! multiplying in a factor 2/sqrt(2) yields the correct result for many FFT versions
            IF (N1==1) GRID%FFTSCA(N,2)=FACTM
#endif
         ENDDO
      ENDIF

    END SUBROUTINE FFTINI


!************************* SUBROUTINE FFTWAV ***************************
!
!> Transforms a wavefunction C defined within the cutoff-sphere
!> to real space CR
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0 (wNGXhalf)
!> are multiplied by a factor sqrt(2) on the reduced plane wave grid
!> this factor has to be removed before the FFT transformation
!> (scaling with FFTSCA(M,2))!
!
!***********************************************************************

    SUBROUTINE FFTWAV(NPL,NINDPW,CR,C,GRID)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      ! local variables
      INTEGER :: M

__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftwav')

#ifdef MPI
      IF (GRID%RL%NFAST/=1) THEN
         CALL FFTWAV_MPI(NPL,NINDPW,CR,C,GRID)
         PROFILING_STOP('fftwav')
         RETURN
      ENDIF
#endif

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(CR).AND. &
__GPU    IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      ! zero all elements on the grid
GPACC !$ACC PARALLEL LOOP PRESENT(CR,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID) PRIVATE(M)
      DO M=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
         CR(M)=(0.0_q,0.0_q)
      ENDDO
__CPU !$OMP END PARALLEL DO

      ! now fill in non zero elements from
      ! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,NINDPW,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,NINDPW,GRID,C) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
__CPU !$OMP END PARALLEL DO

#ifdef inversion_through_mid_point
GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID,C) PRIVATE(M)
         DO M=1,SIZE(GRID%IND_IN_SPHERE)
            CR(GRID%NINDPWCONJG(M))=CONJG(C(GRID%IND_IN_SPHERE(M)))*GRID%FFTSCA(GRID%IND_IN_SPHERE(M),2)
         ENDDO
__CPU !$OMP END PARALLEL DO
#endif
      ELSE

GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,NINDPW) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,NINDPW,C) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)
         ENDDO
__CPU !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D(CR,GRID,1)

      PROFILING_STOP('fftwav')

      RETURN
    END SUBROUTINE FFTWAV


!************************* SUBROUTINE FFTWAV_USEINV ********************
!
!> Transforms a wavefunction C defined within the cutoff-sphere
!> to real space CR
!>
!> @note This routine differs from #FFTWAV in that it allows for the input
!> data to be only defined on half the grid,
!> and uses an inversion through the mid point to determine
!> the PW coefficients at -G-k
!
!***********************************************************************

    SUBROUTINE FFTWAV_USEINV(NPL,NINDPW,NINDPW_INV,FFTSCA,CR,C,GRID)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL),NINDPW_INV(NPL)
      INTEGER    :: NPL
      ! local variables
      INTEGER :: M

__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftwav_useinv')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(NINDPW_INV).AND.IS_ON_DEVICE(FFTSCA).AND. &
__GPU    IS_ON_DEVICE(CR).AND.IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

GPACC !$ACC PARALLEL LOOP PRESENT(CR,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID) PRIVATE(M)
      DO M=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
         CR(M)=(0.0_q,0.0_q)
      ENDDO
__CPU !$OMP END PARALLEL DO

      ! now fill in non zero elements from
      ! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP PRESENT(CR,NINDPW,C,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,NINDPW,C,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,C,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
__CPU !$OMP END PARALLEL DO

#ifdef inversion_through_mid_point
GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID,C) PRIVATE(M)
         DO M=1,SIZE(GRID%IND_IN_SPHERE)
            CR(GRID%NINDPWCONJG(M))=CONJG(C(GRID%IND_IN_SPHERE(M)))*GRID%FFTSCA(GRID%IND_IN_SPHERE(M),2)
         ENDDO
__CPU !$OMP END PARALLEL DO
#endif

      ELSE

GPACC !$ACC PARALLEL LOOP PRESENT(CR,NINDPW_INV,C,FFTSCA) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,NINDPW_INV,C,FFTSCA) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,C,NINDPW_INV,FFTSCA) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW_INV(M))=CONJG(C(M))*FFTSCA(M)
         ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,FFTSCA) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,NINDPW,FFTSCA) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,NINDPW,C,FFTSCA) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*FFTSCA(M)
         ENDDO
__CPU !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D(CR,GRID,1)

      PROFILING_STOP('fftwav_useinv')

      RETURN
    END SUBROUTINE FFTWAV_USEINV


!************************* SUBROUTINE FFTEXT ***************************
!
!> Performes a FFT to reciprocal space and extracts data from the FFT-mesh
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0
!> are multiplied by a factor sqrt(2) on the reduced grid
!> this factor has to be applied after the FFT transformation
!> (scaling with FFTSCA(M))!
!
!***********************************************************************

    SUBROUTINE FFTEXT(NPL,NINDPW,CR,C,GRID,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      LOGICAL    :: LADD
      ! local variables
      INTEGER :: M

__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftext')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(CR).AND. &
__GPU    IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      CALL FFT3D(CR,GRID,-1)

!gK LUSEINV TODO: REAL2CPLX could possibly go
      IF (LADD .AND. GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
__CPU !$OMP END PARALLEL DO

      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
__CPU !$OMP END PARALLEL DO

      ELSE IF (GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
__CPU !$OMP END PARALLEL DO

     ELSE

GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW) IF(OFF_ACTIVE)

!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
__CPU !$OMP END PARALLEL DO

      ENDIF

      PROFILING_STOP('fftext')

      RETURN
    END SUBROUTINE FFTEXT


!************************* SUBROUTINE FFTEXT_USEINV ********************
!
!> Performes a FFT to reciprocal space and extracts data from the FFT-mesh
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0
!> are multiplied by a factor sqrt(2) on the reduced grid
!> this factor has to be applied after the FFT transformation
!> (scaling with FFTSCA(M))!
!
!***********************************************************************

    SUBROUTINE FFTEXT_USEINV(NPL,NINDPW,FFTSCA,CR,C,GRID,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      LOGICAL    :: LADD
      ! local variables
      INTEGER :: M
__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftext_useinv')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(FFTSCA).AND.IS_ON_DEVICE(CR).AND. &
__GPU    IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      CALL FFT3D(CR,GRID,-1)

!gK LUSEINV TODO: REAL2CPLX could possibly go
      IF (LADD .AND. GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,FFTSCA) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,FFTSCA) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*FFTSCA(M)
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
__CPU !$OMP END PARALLEL DO
     ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,FFTSCA) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,FFTSCA) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*FFTSCA(M)
        ENDDO
__CPU !$OMP END PARALLEL DO
      ENDIF

      PROFILING_STOP('fftext_useinv')

      RETURN
    END SUBROUTINE FFTEXT_USEINV


!************************* SUBROUTINE FFT3D ****************************
!
!> 3-d fast fourier transform (possibly real to complex and vice versa)
!> for charge densities and potentials
!>
!> @param ISN = +1  q->r   vr= sum(q) vq exp(+iqr) (might be C2R)\n
!>            = -1  r->q   vq= sum(r) vr exp(-iqr) (might be R2C)
!
!***********************************************************************

    SUBROUTINE FFT3D(C,GRID,ISN)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: C(*)
      INTEGER ::ISN
      ! local variables
      INTEGER :: NX, NY, NZ
      INTEGER :: IL, NDEST, NSRC

      PROFILING_START('fft3d')
#ifdef MPI
      IF (GRID%RL%NFAST/=1 .OR. GRID%RL_FFT%NFAST/=1) THEN
         CALL FFT3D_MPI(C, GRID, ISN)
         PROFILING_STOP('fft3d')
         RETURN
      ENDIF
#endif
      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

!-------------------------------------------------------------------------------
!  complex to complex version
!-------------------------------------------------------------------------------
      IF (.NOT. GRID%REAL2CPLX .AND. .NOT. GRID%LREAL ) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 1 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         CALL WFFTBAS(C,GRID,ISN)
!-------------------------------------------------------------------------------
!  complex to complex version, but with a real array in real space
!-------------------------------------------------------------------------------
      ELSE IF (.NOT. GRID%REAL2CPLX .AND. GRID%LREAL) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 2 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         !     q->r FFT
         IF (ISN==1) THEN
            CALL WFFTBAS(C,GRID,ISN)

            !  go from complex stride 2 to 1
!DIR$ IVDEP
!OCL NOVREC
            DO IL=0,NX*NY*NZ-1
               NDEST=IL+1
               NSRC =IL*2+1
               C(NDEST)=C(NSRC)
            ENDDO
         ELSE

            !     r->q FFT
            !  go from stride 1 to stride 2
!DIR$ IVDEP
!OCL NOVREC
            DO IL=NX*NY*NZ-1,0,-1
               NSRC =IL+1
               NDEST=IL*2+1
               C(NDEST)=C(NSRC)
               C(NDEST+1)=0
            ENDDO
            CALL WFFTBAS(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
!  real to complex FFT  only half grid mode in X direction supported
!  data are stored as real array in real space
!-------------------------------------------------------------------------------
      ELSE IF (GRID%LREAL) THEN
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 3 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         !  in real space the first dimension in VASP is NGX (REAL data)
         !  but the FFT required NGX+2 (real data)
         !  therefore some data movement is required

         IF (ISN==1) THEN
            !     q->r FFT
            CALL WFFTBRC(C,GRID,ISN)
            !  concat  x-lines (go from stride NX+2 to NX)
            CALL RESTRIDE_Q2R(NY*NZ,NX,C)
         ELSE
            !     x-lines (go from stride NX to NX+2)
            CALL RESTRIDE_R2Q(NY*NZ,NX,C)
            !     r->q FFT
            CALL WFFTBRC(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
! same as above (real to complex FFT) but this time the data layout
! is complex in real space
!-------------------------------------------------------------------------------
      ELSE
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 4 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         IF (ISN==1) THEN
            !     q->r FFT
            CALL WFFTBRC(C,GRID,ISN)
            ! concat  x-lines (go from "real" stride NX+2 to complex stride NX)
            CALL RESTRIDE_Q2R_CMPLX(NY*NZ,NX,C)
         ELSE
            !     x-lines (go from complex stride NX to real stride NX+2)
            CALL RESTRIDE_R2Q_CMPLX(NY*NZ,NX,C)
            !     r->q FFT
            CALL WFFTBRC(C,GRID,ISN)
         ENDIF
      ENDIF

      PROFILING_STOP('fft3d')

      RETURN
    END SUBROUTINE FFT3D


!************************* SUBROUTINE EXTBAS ***************************
!
!> Extracts data from the FFT-mesh
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0
!> are multiplied by a factor sqrt(2) on the reduced grid
!> this factor has to be applied after the FFT transformation
!> (scaling with FFTSCA(M))!
!
!***********************************************************************

    SUBROUTINE EXTBAS(NPL,NINDPW,CR,C,GRID,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%MPLWV)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      LOGICAL    :: LADD
      ! local variables
      INTEGER :: M
__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('extbas')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(CR).AND. &
__GPU    IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

!gK LUSEINV TODO: REAL2CPLX could possibly go
      IF (LADD .AND. GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
__CPU !$OMP END PARALLEL DO
     ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW) IF(OFF_ACTIVE)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
__CPU !$OMP END PARALLEL DO
      ENDIF

      PROFILING_STOP('extbas')

      RETURN
    END SUBROUTINE EXTBAS


!************************* SUBROUTINE MULZ *****************************
!
!> Multiplies the Z!=0 components by a factor FACT
!
!***********************************************************************

    SUBROUTINE MULZ(C,NGX,NGY,NGZ,FACT)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(0:NGX-1,0:NGY-1,0:NGZ-1)

      DO N3=1,NGZ/2-1
       DO N2=0,NGY-1
        DO N1=0,NGX-1
          C(N1,N2,N3)= C(N1,N2,N3)*FACT
      ENDDO
      ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULZ


    SUBROUTINE RESTRIDE_Q2R(NCOL,N,C)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
    ! local
      INTEGER :: ICOL,NDEST,NSRC,I
#ifdef OFFLOADING
    ! to be able to collapse the loops ICOL and I below simultaneously
    ! data needs to be stored in a temporary array
      REAL(q), ALLOCATABLE :: CTMP(:)

      IF (NCOL<2) RETURN

      ALLOCATE(CTMP(NCOL*(N+2)))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CTMP,C(1:(NCOL-1)*(N+2)+N)) IF(OFFLOAD_ON)

      DO I=1,(NCOL-1)*(N+2)+N
         CTMP(I)=C(I)
      ENDDO

#endif
!DIR$ IVDEP
!OCL NOVREC

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1:NCOL*N),CTMP) &
GPOMP !$OMP PRIVATE(NSRC,NDEST) IF(OFFLOAD_ON)

      DO ICOL=1,NCOL-1
__CPU    NDEST=ICOL* N
__CPU    NSRC =ICOL*(N+2)
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
         DO I=1,N
__GPU       NDEST=ICOL* N
__GPU       NSRC =ICOL*(N+2)
__CPU       C(NDEST+I)=   C(NSRC+I)
__GPU       C(NDEST+I)=CTMP(NSRC+I)
         ENDDO
      ENDDO
#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_Q2R


    SUBROUTINE RESTRIDE_R2Q(NCOL,N,C)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
    ! local
      INTEGER :: ICOL,NDEST,NSRC,I
#ifdef OFFLOADING
    ! to be able to collapse the loops ICOL and I below simultaneously
    ! data needs to be stored in a temporary array
      REAL(q), ALLOCATABLE :: CTMP(:)

      IF (NCOL<2) RETURN

      ALLOCATE(CTMP(NCOL*N))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CTMP,C(1:NCOL*N)) IF(OFFLOAD_ON)

      DO I=1,NCOL*N
         CTMP(I)=C(I)
      ENDDO
#endif

!DIR$ IVDEP
!OCL NOVREC

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1:(NCOL-1)*(N+2)+N),CTMP) &
GPOMP !$OMP PRIVATE(NSRC,NDEST) IF(OFFLOAD_ON)

      DO ICOL=NCOL-1,1,-1
__CPU    NSRC =ICOL* N
__CPU    NDEST=ICOL*(N+2)
! ifc10.1 has trouble vectorizing this statement
!!DIR$ IVDEP
!!OCL NOVREC
         DO I=N,1,-1
__GPU       NSRC =ICOL* N
__GPU       NDEST=ICOL*(N+2)
__CPU       C(NDEST+I)=   C(NSRC+I)
__GPU       C(NDEST+I)=CTMP(NSRC+I)
         ENDDO
      ENDDO

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_R2Q


    SUBROUTINE RESTRIDE_Q2R_CMPLX(NCOL,N,C)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
    ! local
      INTEGER :: ICOL,NDEST,NSRC,I
#ifdef OFFLOADING
    ! to be able to collapse the loops ICOL and I below simultaneously
    ! data needs to be stored in a temporary array
      REAL(q), ALLOCATABLE :: CTMP(:)

      IF (NCOL<1) RETURN

      ALLOCATE(CTMP(NCOL*(N+2)))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPACC !$ACC PARALLEL LOOP PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CTMP,C(1:NCOL*(N+2))) IF(OFFLOAD_ON)
      DO I=1,(NCOL-1)*(N+2)+N
         CTMP(I)=C(I)
      ENDDO
#endif
!DIR$ IVDEP
!OCL NOVREC
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CTMP,C(1:NCOL*(N+2))) &
GPOMP !$OMP PRIVATE(NSRC,NDEST) IF(OFFLOAD_ON)
      DO ICOL=NCOL-1,0,-1
__CPU    NDEST=ICOL* N*2
__CPU    NSRC =ICOL*(N+2)
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
         DO I=N,1,-1
__GPU       NDEST=ICOL* N*2
__GPU       NSRC =ICOL*(N+2)
__CPU       C(NDEST+I*2-1)=C(NSRC+I)
__GPU       C(NDEST+I*2-1)=CTMP(NSRC+I)
            C(NDEST+I*2  )=0
         ENDDO
      ENDDO

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_Q2R_CMPLX


    SUBROUTINE RESTRIDE_R2Q_CMPLX(NCOL,N,C)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      USE tutor, ONLY: vtutor
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
    ! local
      INTEGER :: ICOL,NDEST,NSRC,I
#ifdef OFFLOADING
    ! to be able to collapse the loops ICOL and I below simultaneously
    ! data needs to be stored in a temporary array
      REAL(q), ALLOCATABLE :: CTMP(:)

      IF (NCOL<1) RETURN

      ALLOCATE(CTMP(2*NCOL*N))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPACC !$ACC PARALLEL LOOP PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CTMP,C(1:2*NCOL*N)) IF(OFFLOAD_ON)
      DO I=1,2*NCOL*N
         CTMP(I)=C(I)
      ENDDO
#else
      IF (NCOL<1) RETURN
      IF (N<2) CALL vtutor%bug("RESTRIDE_R2Q_CMPLX: can not restride array for N=1", __FILE__, __LINE__)
#endif
!DIR$ IVDEP
!OCL NOVREC
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CTMP,C(1:2*NCOL*N)) &
GPOMP !$OMP PRIVATE(NSRC,NDEST) IF(OFFLOAD_ON)
      DO ICOL=0,NCOL-1
__CPU    NSRC =ICOL* N*2
__CPU    NDEST=ICOL*(N+2)
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
         DO I=1,N
__GPU       NSRC =ICOL* N*2
__GPU       NDEST=ICOL*(N+2)
__CPU       C(NDEST+I)=   C(NSRC+I*2-1)
__GPU       C(NDEST+I)=CTMP(NSRC+I*2-1)
         ENDDO
      ENDDO
#ifdef OFFLOADING
GPACC  !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_R2Q_CMPLX


!************************* SUBROUTINE FFT1D ****************************
!
!> One dimensional FFT
!> Currently only used in stm.F (which is pretty irrelevant)
!
!***********************************************************************

    SUBROUTINE FFT1D(C,NFFT,ISIGN)
      USE prec

      IMPLICIT NONE

      COMPLEX(q) :: C(*)
      INTEGER    :: NFFT
      INTEGER    :: ISIGN

      CALL FFT1D_C2C_PLAN(C,NFFT,ISIGN)

    END SUBROUTINE FFT1D


!************************* SUBROUTINE FFTWAV_MU ************************
!
!***********************************************************************

    SUBROUTINE FFTWAV_MU(NPL,N,NINDPW,CR,LDR,C,LDC,GRID)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC

      ! local variables
      INTEGER :: I,J
__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftwav_mu')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(CR).AND. &
__GPU   IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      ! zero all elements on the grid
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR(1,1),GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,GRID,CR) PRIVATE(I,J)
      DO I=1,N
         DO J=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
            CR(J,I)=(0.0_q,0.0_q)
         ENDDO
      ENDDO
__CPU !$OMP END PARALLEL DO

      ! now fill in non zero elements from
      ! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR(1,1),C(1,1),NINDPW,GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*GRID%FFTSCA(J,2)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO
#ifdef inversion_through_mid_point
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR(1,1),C(1,1),GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,CR,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,SIZE(GRID%IND_IN_SPHERE)
               CR(GRID%NINDPWCONJG(J),I)=CONJG(C(GRID%IND_IN_SPHERE(J),I))*GRID%FFTSCA(GRID%IND_IN_SPHERE(J),2)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO
#endif

      ELSE

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR(1,1),C(1,1),NINDPW) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C) PRIVATE(I,J)
         DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D_MU(N,CR,LDR,GRID,1)

      PROFILING_STOP('fftwav_mu')

    END SUBROUTINE FFTWAV_MU


!************************* SUBROUTINE FFTWAV_USEINV_MU *****************
!
!***********************************************************************

    SUBROUTINE FFTWAV_USEINV_MU(NPL,N,NINDPW,NINDPW_INV,FFTSCA,CR,LDR,C,LDC,GRID)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL),NINDPW_INV(NPL)
      INTEGER    :: NPL, N, LDR, LDC

      ! local variables
      INTEGER :: I,J
__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftwav_useinv_mu')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(NINDPW_INV).AND.IS_ON_DEVICE(FFTSCA).AND. &
__GPU    IS_ON_DEVICE(CR).AND.IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      ! zero all elements on the grid
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CR(1,1),GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,GRID,CR) PRIVATE(I,J)
      DO I=1,N
         DO J=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
            CR(J,I)=(0.0_q,0.0_q)
         ENDDO
      ENDDO
__CPU !$OMP END PARALLEL DO

      ! now fill in non zero elements from
      ! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CR(1,1),C(1,1),NINDPW,GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,GRID) PRIVATE(I,J)
         DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*GRID%FFTSCA(J,2)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO
#ifdef inversion_through_mid_point
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,GRID) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CR(1,1),C(1,1),GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,CR,C,GRID) PRIVATE(I,J)
         DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO J=1,SIZE(GRID%IND_IN_SPHERE)
               CR(GRID%NINDPWCONJG(J),I)=CONJG(C(GRID%IND_IN_SPHERE(J),I))*GRID%FFTSCA(GRID%IND_IN_SPHERE(J),2)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO
#endif

      ELSE

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW_INV,FFTSCA) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CR(1,1),C(1,1),NINDPW_INV,FFTSCA) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW_INV,C,FFTSCA) PRIVATE(I,J)
         DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO J=1,NPL
               CR(NINDPW_INV(J),I)=CONJG(C(J,I))*FFTSCA(J)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,FFTSCA) IF(OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CR(1,1),C(1,1),NINDPW,FFTSCA) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,FFTSCA) PRIVATE(I,J)
         DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*FFTSCA(J)
            ENDDO
         ENDDO
__CPU !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D_MU(N,CR,LDR,GRID,1)

      PROFILING_STOP('fftwav_useinv_mu')

    END SUBROUTINE FFTWAV_USEINV_MU


!************************* SUBROUTINE FFTEXT_MU ************************
!
!***********************************************************************

    SUBROUTINE FFTEXT_MU(NPL,N,NINDPW,CR,LDR,C,LDC,GRID,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC
      LOGICAL    :: LADD

      ! local variables
      INTEGER :: I,J
__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftext_mu')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(CR).AND. &
__GPU    IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      CALL FFT3D_MU(N,CR,LDR,GRID,-1)

      IF (GRID%REAL2CPLX) THEN
         IF (LADD) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW,GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ELSE
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW,GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ENDIF
      ELSE
         IF (LADD) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ELSE
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ENDIF
      ENDIF

      PROFILING_STOP('fftext_mu')

    END SUBROUTINE FFTEXT_MU


!************************* SUBROUTINE FFTEXT_USEINV_MU *****************
!
!***********************************************************************

    SUBROUTINE FFTEXT_USEINV_MU(NPL,N,NINDPW,FFTSCA,CR,LDR,C,LDC,GRID,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC
      LOGICAL    :: LADD

      ! local variables
      INTEGER :: I,J
__GPU LOGICAL :: OFF_ACTIVE

      PROFILING_START('fftext_useinv_mu')

__GPU OFF_ACTIVE=IS_ON_DEVICE(NINDPW).AND.IS_ON_DEVICE(FFTSCA).AND.IS_ON_DEVICE(CR).AND. &
__GPU    IS_ON_DEVICE(C).AND.IS_ON_DEVICE(GRID,1).AND.OFFLOAD_ON

      CALL FFT3D_MU(N,CR,LDR,GRID,-1)

      IF (GRID%REAL2CPLX) THEN
         IF (LADD) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW,GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW,GRID) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ENDIF
      ELSE
         IF (LADD) THEN
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,FFTSCA) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW,FFTSCA) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,FFTSCA) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*FFTSCA(J)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (OFF_ACTIVE) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C(1,1),CR(1,1),NINDPW,FFTSCA) IF(OFF_ACTIVE)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,FFTSCA) PRIVATE(I,J)
            DO I=1,N
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*FFTSCA(J)
               ENDDO
            ENDDO
__CPU !$OMP END PARALLEL DO
         ENDIF
      ENDIF

      PROFILING_STOP('fftext_useinv_mu')

    END SUBROUTINE FFTEXT_USEINV_MU


!************************* SUBROUTINE FFT3D_MU *************************
!
!***********************************************************************

    SUBROUTINE FFT3D_MU(N,C,LDC,GRID,ISN)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(LDC,*)
      INTEGER :: N,LDC,ISN

      ! local variables
      INTEGER :: NX,NY,NZ

      PROFILING_START('fft3d_mu')

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IF (.NOT. GRID%REAL2CPLX) THEN
!-------------------------------------------------------------------------------
!  complex to complex FFTs
!-------------------------------------------------------------------------------
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("FFT3D_MU: grid dimensions not properly set (1) " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // &
               str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // str(GRID%NGZ_rd), &
               __FILE__, __LINE__)
         ENDIF

         IF  (.NOT. GRID%LREAL ) THEN

            CALL WFFTBAS_MU(N,C,LDC,GRID,ISN)

         ELSE

            CALL vtutor%bug("FFT3D_MU: batched ffts for REAL2CPLX = " // str(GRID%REAL2CPLX) // &
               " and LREAL = " // str(GRID%LREAL) //" not supported yet",__FILE__,__LINE__)

         ENDIF
      ELSE
!-------------------------------------------------------------------------------
!  real to complex FFTs
!-------------------------------------------------------------------------------
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("FFT3D_MU: grid dimensions not properly set (2)" // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // &
               str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " //  str(GRID%NGZ_rd), &
                __FILE__, __LINE__)
         ENDIF

         IF (.NOT. GRID%LREAL) THEN
!
! the data in real space is stored as an array of complex numbers (with Im(c)=0)
!
            IF (ISN==1) THEN
               ! q -> r FFT
               CALL WFFTBRC_MU(N,C,LDC,GRID,ISN)
               ! contract columns along x (go from stride NX+2 to stride NX)
               ! and store real valued transform as a array of complex numbers
               ! (with Im(c) = 0)
               CALL RESTRIDE_Q2R_CMPLX_MU(N,NY*NZ,NX,C,LDC)
            ELSE
               ! expand columns along x (go from stride NX to stride NX+2)
               ! and store array of complex numbers with Im(c)=0 as array
               ! of reals
               CALL RESTRIDE_R2Q_CMPLX_MU(N,NY*NZ,NX,C,LDC)
               ! r -> q FFT
               CALL WFFTBRC_MU(N,C,LDC,GRID,ISN)
            ENDIF

         ELSE
!
! the data in real space is stored as an array of real numbers
!
            IF (ISN==1) THEN
               ! q -> r FFT
               CALL WFFTBRC_MU(N,C,LDC,GRID,ISN)
               ! contract columns along x (go from stride NX+2 to stride NX)
               CALL RESTRIDE_Q2R_MU(N,NY*NZ,NX,C,LDC)
            ELSE
               ! expand columns along x (go from stride NX to stride NX+2)
               CALL RESTRIDE_R2Q_MU(N,NY*NZ,NX,C,LDC)
               ! q -> r FFT
               CALL WFFTBRC_MU(N,C,LDC,GRID,ISN)
            ENDIF

         ENDIF
      ENDIF

      PROFILING_STOP('fft3d_mu')

    END SUBROUTINE FFT3D_MU


!************************* SUBROUTINE RESTRIDE_Q2R_CMPLX_MU ************
!
!> The output of an inplace complex-to-real FFT on a grid of NX*NY*NZ
!> points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!> dimension is padded with zeros.
!>
!> RESTRIDE_Q2R_CMPLX_MU does two things:
!>
!> * The padding is removed, i.e., the data is re-stored as an
!>   ( NX,NY,NZ ) array.
!> * The data is stored as an array of complex numbers (with Im(c)=0).
!
!***********************************************************************

    SUBROUTINE RESTRIDE_Q2R_CMPLX_MU(N,NCOL,NX,C,LDC)

__GPU USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef OFFLOADING
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(NCOL*(NX+2),N))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CTMP,C(1,1)) IF(OFFLOAD_ON)
      DO J=1,N
         DO I=1,(NCOL-1)*(NX+2)+NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C(1,1),CTMP) PRIVATE(NDEST,NSRC) IF(OFFLOAD_ON)
__CPU !$OMP PARALLEL DO PRIVATE(J,ICOL,NDEST,NSRC,I)
      DO J=1,N
         DO ICOL=NCOL-1,0,-1
__CPU       NDEST=ICOL* NX*2
__CPU       NSRC =ICOL*(NX+2)
            DO I=NX,1,-1
__GPU          NDEST=ICOL* NX*2
__GPU          NSRC =ICOL*(NX+2)
__CPU          C(NDEST+I*2-1,J)=C   (NSRC+I,J)
__GPU          C(NDEST+I*2-1,J)=CTMP(NSRC+I,J)
               C(NDEST+I*2  ,J)=0
            ENDDO
         ENDDO
      ENDDO
__CPU !$OMP END PARALLEL DO

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_Q2R_CMPLX_MU


!************************* SUBROUTINE RESTRIDE_R2Q_CMPLX_MU ************
!
!> RESTRIDE_R2Q_CMPLX_MU adresses two peculiarities:
!>
!> * The input of an inplace real-to-complex FFT on a grid of NX*NY*NZ
!>   points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!>   dimension is padded with zeros.
!>
!> * In many cases the data (even though it is real) in real space is
!>   stored as an array of complex numbers (with Im(c)=0)
!>
!> RESTRIDE_R2Q_CMPLX_MU does two things:
!>
!> * Padding is added, i.e., the data is re-stored as an ( NX+2,NY,NZ )
!>   array.
!> * The data is stored as an array of real numbers.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_R2Q_CMPLX_MU(N,NCOL,NX,C,LDC)

__GPU USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef OFFLOADING
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(2*NCOL*NX,N))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CTMP,C(1,1)) IF(OFFLOAD_ON)
      DO J=1,N
         DO I=1,2*NCOL*NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C(1,1),CTMP) PRIVATE(NSRC,NDEST) IF(OFFLOAD_ON)
__CPU !$OMP PARALLEL DO PRIVATE(J,ICOL,NSRC,NDEST,I)
      DO J=1,N
         DO ICOL=0,NCOL-1
__CPU       NSRC =ICOL* NX*2
__CPU       NDEST=ICOL*(NX+2)
            DO I=1,NX
__GPU          NSRC =ICOL* NX*2
__GPU          NDEST=ICOL*(NX+2)
__CPU          C(NDEST+I,J)=C   (NSRC+I*2-1,J)
__GPU          C(NDEST+I,J)=CTMP(NSRC+I*2-1,J)
            ENDDO
         ENDDO
      ENDDO
__CPU !$OMP END PARALLEL DO

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_R2Q_CMPLX_MU


!************************* SUBROUTINE RESTRIDE_Q2R_MU ******************
!
!> The output of an inplace complex-to-real FFT on a grid of NX*NY*NZ
!> points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!> dimension is padded with zeros.
!>
!> RESTRIDE_Q2R_MU removes the padding, i.e., the data is re-stored as
!> an ( NX,NY,NZ ) array.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_Q2R_MU(N,NCOL,NX,C,LDC)

__GPU USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef OFFLOADING
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(NCOL*(NX+2),N))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CTMP,C(1,1)) IF(OFFLOAD_ON)
      DO J=1,N
         DO I=1,(NCOL-1)*(NX+2)+NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C(1,1),CTMP) PRIVATE(NDEST,NSRC) IF(OFFLOAD_ON)
__CPU !$OMP PARALLEL DO PRIVATE(J,ICOL,NDEST,NSRC,I)
      DO J=1,N
         DO ICOL=1,NCOL-1
__CPU       NDEST=ICOL* NX
__CPU       NSRC =ICOL*(NX+2)
            DO I=1,NX
__GPU          NDEST=ICOL* NX
__GPU          NSRC =ICOL*(NX+2)
__CPU          C(NDEST+I,J)=C   (NSRC+I,J)
__GPU          C(NDEST+I,J)=CTMP(NSRC+I,J)
            ENDDO
         ENDDO
      ENDDO
__CPU !$OMP END PARALLEL DO

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_Q2R_MU


!************************* SUBROUTINE RESTRIDE_R2Q_MU ******************
!
!> The input of an inplace real-to-complex FFT on a grid of NX*NY*NZ
!> points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!> dimension is padded with zeros.
!>
!> RESTRIDE_R2Q_MU add this padding, i.e., the data is re-stored as an
!> ( NX+2,NY,NZ ) array.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_R2Q_MU(N,NCOL,NX,C,LDC)

__GPU USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef OFFLOADING
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(NCOL*NX,N))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CTMP) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:CTMP,C(1,1)) IF(OFFLOAD_ON)
      DO J=1,N
         DO I=1,NCOL*NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C(1,1),CTMP) PRIVATE(NSRC,NDEST) IF(OFFLOAD_ON)
__CPU !$OMP PARALLEL DO PRIVATE(J,ICOL,NSRC,NDEST,I)
      DO J=1,N
         DO ICOL=NCOL-1,1,-1
__CPU       NSRC =ICOL* NX
__CPU       NDEST=ICOL*(NX+2)
            DO I=NX,1,-1
__GPU          NSRC =ICOL* NX
__GPU          NDEST=ICOL*(NX+2)
__CPU          C(NDEST+I,J)=C   (NSRC+I,J)
__GPU          C(NDEST+I,J)=CTMP(NSRC+I,J)
            ENDDO
         ENDDO
      ENDDO
__CPU !$OMP END PARALLEL DO

#ifdef OFFLOADING
GPACC !$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CTMP) IF(OFFLOAD_ON)
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_R2Q_MU


#ifdef MPI
!************************* SUBROUTINE FFTGRIDPLAN_MPI ******************
!
!***********************************************************************

    SUBROUTINE FFTGRIDPLAN_MPI(GRID)
      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE (grid_3d) GRID

      ! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
      INTEGER, PARAMETER :: NADD=1024, NTIMES=4
      INTEGER :: I

      ALLOCATE(CWORK(GRID%MPLWV+NADD))

      DO I=1,NTIMES
         CALL INIDAT(GRID%MPLWV+NADD,CWORK)
         CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
!$       CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP END MASTER
!$OMP END PARALLEL
      ENDDO

      DEALLOCATE(CWORK)

      RETURN
    END SUBROUTINE FFTGRIDPLAN_MPI


!************************* SUBROUTINE FFTINI_MPI ***********************
!
!  if necessary this routine performes initialization
!  for FFTWAV and FFTEXT
!  usually this is only necessary for the Gamma point only
!  1-kpoint version
!
!   FFTSCA(.,1) is the scaling factor for extracting the wavefunction
!               from the FFT grid (FFTEXT)
!   FFTSCA(.,2) is the scaling factor for puting the wavefunction on
!               the grid
!***********************************************************************

    SUBROUTINE FFTINI_MPI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
      USE prec
      USE mgrid_struct_def
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      INTEGER :: NKPTS, NRPLW
      INTEGER :: NPLWKP(NKPTS)
      INTEGER :: NINDPW(NRPLW,NKPTS)
      ! local variables
      REAL(q) :: FACTM
      INTEGER :: NK, NPL, N, IND, N1, NC, N2, N3

      IF (GRID%REAL2CPLX) THEN
         IF (NKPTS>1) THEN
            CALL vtutor%error("FFT3D: real version works only for 1 k-point")
         ENDIF

         NK=1
         NPL=NPLWKP(NK)
         NULLIFY(GRID%FFTSCA)
         ALLOCATE(GRID%FFTSCA(NPL,2))

         DO N=1,NPL
            IND=NINDPW(N,NK)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)

            FACTM=SQRT(2._q)
            IF (N1==1 .AND. N2==1 .AND. N3==1) FACTM=1
            GRID%FFTSCA(N,1)= FACTM
            GRID%FFTSCA(N,2)= 1/FACTM
            ! this statment is required
            ! because for z==0 only half of the FFT components are set
            ! upon calling FFTWAV
            IF (N3==1) GRID%FFTSCA(N,2)=FACTM
         ENDDO
      END IF
      RETURN
    END SUBROUTINE FFTINI_MPI


!************************* SUBROUTINE FFTWAV_MPI ***********************
!
!  this subroutine transforms a wavefunction C defined  within  the
!  cutoff-sphere to real space CR
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGZ != 0
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be removed before the FFT transformation !
! (scaling with   FFTSCA(M,2))
!
!> @details @ref openmp :
!> the loops that fill CR with the reciprocal space components of
!> an orbital are distributed over all available OpenMP threads.
!
!***********************************************************************

    SUBROUTINE FFTWAV_MPI(NPL,NINDPW,CR,C,GRID)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL), CR(GRID%MPLWV)
      INTEGER :: NPL, NINDPW(NPL)
      ! local variables
      INTEGER :: M

      PROFILING_START('fftwav_mpi')

      IF (GRID%LREAL) THEN
!DIR$ IVDEP
!OCL NOVREC
GPACC !$ACC PARALLEL LOOP PRESENT(CR) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR) IF(OFFLOAD_ON)
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(CR,GRID)
         DO M=1,GRID%RL%NCOL*GRID%NGZ/2
            CR(M)=0.0_q
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(CR,GRID)
GPACC !$ACC PARALLEL LOOP PRESENT(CR) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR) IF(OFFLOAD_ON)
         DO M=1,GRID%RL%NCOL*GRID%NGZ
            CR(M)=0.0_q
         ENDDO
__CPU !$OMP END PARALLEL DO
      ENDIF

      IF (GRID%REAL2CPLX) THEN
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW,GRID)
GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,GRID%FFTSCA) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,NINDPW,GRID%FFTSCA) IF(OFFLOAD_ON)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW)
GPACC !$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CR,C,NINDPW) IF(OFFLOAD_ON)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)
         ENDDO
__CPU !$OMP END PARALLEL DO
      ENDIF

      CALL FFT3D_MPI(CR,GRID,1)

      PROFILING_STOP('fftwav_mpi')

      RETURN
    END SUBROUTINE FFTWAV_MPI


!************************* SUBROUTINE FFTEXT_MPI ***********************
!
! this subroutine performes a FFT to reciprocal space and extracts data
! from the FFT-mesh
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGX != 0
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be applied after the FFT transformation !
!  (scaling with   FFTSCA(M))
!
!> @details @ref openmp :
!> the loops that extract the reciprocal space components of an orbital
!> from the full FFT grid are distributed over all available OpenMP
!> threads.
!
!***********************************************************************

    SUBROUTINE FFTEXT_MPI(NPL,NINDPW,CR,C,GRID,LADD)

__GPU USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%MPLWV)
      INTEGER :: NPL, NINDPW(NPL)
      LOGICAL :: LADD
      ! local variables
      INTEGER :: M

      PROFILING_START('fftext_mpi')

      CALL FFT3D_MPI(CR,GRID,-1)

      IF (LADD .AND. GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID%FFTSCA) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID%FFTSCA) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW, GRID)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(NPL,CR,NINDPW) REDUCTION(+:C)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID%FFTSCA) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW,GRID%FFTSCA) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW,GRID)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
__CPU !$OMP END PARALLEL DO
      ELSE
GPACC !$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CR,NINDPW) IF(OFFLOAD_ON)
!DIR$ IVDEP
!OCL NOVREC
__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
__CPU !$OMP END PARALLEL DO
      ENDIF

      PROFILING_STOP('fftext_mpi')

      RETURN
    END SUBROUTINE FFTEXT_MPI


!************************* SUBROUTINE FFT3D_MPI ************************
!
!    3-d fast fourier transform (possibly real to complex and vice versa)
!    for chardensities and potentials
!     +1  q->r   vr= sum(q) vq exp(+iqr)    (might be complex to real)
!     -1  r->q   vq= sum(r) vr exp(-iqr)    (might be real to complex)
!
!***********************************************************************

    SUBROUTINE FFT3D_MPI(C,GRID,ISN)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: C(*)
      INTEGER :: ISN
      ! local variables
      INTEGER :: NX, NY, NZ

      PROFILING_START('fft3d_mpi')

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IF (GRID%RL%NFAST==1) THEN
!-------------------------------------------------------------------------------
! parallel FFT with serial data layout  (GRID%RL%NFAST==1)
!-------------------------------------------------------------------------------
         ! complex to complex case
         IF (.NOT. GRID%REAL2CPLX) THEN
            IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
               CALL vtutor%bug("internal error 1 in FFT3D_MPI: something not properly set " // &
                  str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) &
                  // " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
                  str(GRID%NGZ_rd), __FILE__, __LINE__)
            ENDIF

            IF (ISN==1) THEN
               ! q->r FFT
               CALL WFFTBAS_MPI(C,GRID,ISN)
               ! bring real space data from parallel to serial data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTPAR_TO_SER(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
            ELSE
               ! bring real space data from serial to parallel data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTSER_TO_PAR(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
               ! r->q FFT
               CALL WFFTBAS_MPI(C,GRID,ISN)
            ENDIF

         ! complex to real case
         ELSE IF (GRID%LREAL) THEN
            IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
               CALL vtutor%bug("internal error 2 in FFT3D_MPI: something not properly set " // &
                  str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) &
                  // " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
                  str(GRID%NGZ_rd), __FILE__, __LINE__)
            ENDIF

            !  in real space the first dimension in VASP is NGZ (REAL data)
            !  but the FFT requires NGZ+2 (real data)
            !  therefore some data movement is required
            IF (ISN==1) THEN
               ! q->r FFT
               CALL WFFTBRC_MPI(C,GRID,ISN)
               ! x-lines (go from stride NZ+2 to NZ)
               CALL RESTRIDE_Q2R(GRID%RL_FFT%NCOL,NZ,C)
               ! bring real space data from parallel to serial data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTPAR_TO_SER_REAL(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
            ELSE
               ! bring real space data from serial to parallel data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTSER_TO_PAR_REAL(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
               ! x-lines (go from stride NZ to NZ+2)
               CALL RESTRIDE_R2Q(GRID%RL_FFT%NCOL,NZ,C)
               ! r->q FFT
               CALL WFFTBRC_MPI(C,GRID,ISN)
            ENDIF
         ELSE
            CALL vtutor%error("ERROR in FFT3D_MPI: this version does not support the required half &
               &grid mode \n " // str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // str(GRID%NGZ_rd))
         ENDIF
!-------------------------------------------------------------------------------
!  complex parallel FFT
!-------------------------------------------------------------------------------
      ELSE IF (.NOT. GRID%REAL2CPLX) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 3 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF
         CALL WFFTBAS_MPI(C,GRID,ISN)
!-------------------------------------------------------------------------------
!  real to complex parallel FFT
!-------------------------------------------------------------------------------
      ELSE IF (GRID%LREAL) THEN
         IF (.NOT.(NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 4 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         !  in real space the first dimension in VASP is NGZ (REAL data)
         !  but the FFT requires NGZ+2 (real data)
         !  therefore some data movement is required
         IF (ISN==1) THEN
            ! q->r FFT
            CALL WFFTBRC_MPI(C,GRID,ISN)
            ! concat z-lines (go from stride NZ+2 to NZ)
            CALL RESTRIDE_Q2R(GRID%RL%NCOL,NZ,C)
         ELSE
            ! concat z-lines (go from stride NZ to NZ+2)
            CALL RESTRIDE_R2Q(GRID%RL%NCOL,NZ,C)
            ! r->q FFT
            CALL WFFTBRC_MPI(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
!  real to complex parallel FFT with complex storage layout in real space
!-------------------------------------------------------------------------------
      ELSE
         IF (.NOT.(NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 5 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         IF (ISN==1) THEN
            ! q->r FFT
            CALL WFFTBRC_MPI(C,GRID,ISN)
            ! concat z-lines (go from stride NZ+2 to NZ)
            CALL RESTRIDE_Q2R_CMPLX(GRID%RL%NCOL,NZ,C)
         ELSE
            ! z-lines (go from complex stride NZ to real stride NZ+2)
            CALL RESTRIDE_R2Q_CMPLX(GRID%RL%NCOL,NZ,C)
            ! r->q FFT
            CALL WFFTBRC_MPI(C,GRID,ISN)
         ENDIF
      ENDIF

      PROFILING_STOP('fft3d_mpi')

    END SUBROUTINE FFT3D_MPI


!************************* SUBROUTINE FFTPAR_TO_SER ********************
!
! change data layout from parallel to serial data layout
! and vice versa for complex and real arrays
! operates usually in real space
!
!***********************************************************************

    SUBROUTINE FFTPAR_TO_SER(NGX, NGY, NGZ, CORIG)

__GPU USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      COMPLEX(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      COMPLEX(q) :: C(NGX*NGY*NGZ)

GPACC !$ACC ENTER DATA CREATE(C) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:C) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=0,NGX-1
         DO IY=0,NGY-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IX,IY,IZ)=CORIG(IZ,IX,IY)
               C(1+IX+NGX*(IY+NGY*IZ))=CORIG(1+IZ+NGZ*(IX+NGX*IY))
            ENDDO
         ENDDO
      ENDDO
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
GPACC !$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
GPACC !$ACC EXIT DATA DELETE(C) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:C) IF(OFFLOAD_ON)
    END SUBROUTINE FFTPAR_TO_SER


    SUBROUTINE FFTPAR_TO_SER_REAL(NGX, NGY, NGZ, CORIG)

__GPU USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      REAL(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      REAL(q) :: C(NGX*NGY*NGZ)

GPACC !$ACC ENTER DATA CREATE(C) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:C) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=0,NGX-1
         DO IY=0,NGY-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IX,IY,IZ)=CORIG(IZ,IX,IY)
               C(1+IX+NGX*(IY+NGY*IZ))=CORIG(1+IZ+NGZ*(IX+NGX*IY))
            ENDDO
         ENDDO
      ENDDO
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
GPACC !$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
GPACC !$ACC EXIT DATA DELETE(C) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:C) IF(OFFLOAD_ON)
    END SUBROUTINE FFTPAR_TO_SER_REAL


    SUBROUTINE FFTSER_TO_PAR(NGX, NGY, NGZ, CORIG)

__GPU USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      COMPLEX(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      COMPLEX(q) :: C(NGX*NGY*NGZ)

GPACC !$ACC ENTER DATA CREATE(C) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:C) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=0,NGX-1
         DO IY=0,NGY-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IZ,IX,IY)=CORIG(IX,IY,IZ)
               C(1+IZ+NGZ*(IX+NGX*IY))=CORIG(1+IX+NGX*(IY+NGY*IZ))
            ENDDO
         ENDDO
      ENDDO
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
GPACC !$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
GPACC !$ACC EXIT DATA DELETE(C) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:C) IF(OFFLOAD_ON)
    END SUBROUTINE FFTSER_TO_PAR


    SUBROUTINE FFTSER_TO_PAR_REAL(NGX, NGY, NGZ, CORIG)

__GPU USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      REAL(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      REAL(q) :: C(NGX*NGY*NGZ)

GPACC !$ACC ENTER DATA CREATE(C) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:C) IF(OFFLOAD_ON)

GPACC !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=0,NGX-1
         DO IY=0,NGY-1
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IZ,IX,IY)=CORIG(IX,IY,IZ)
               C(1+IZ+NGZ*(IX+NGX*IY))=CORIG(1+IX+NGX*(IY+NGY*IZ))
            ENDDO
         ENDDO
      ENDDO
__CPU !DIR$ IVDEP
__CPU !OCL NOVREC
GPACC !$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:C,CORIG) IF(OFFLOAD_ON)
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
GPACC !$ACC EXIT DATA DELETE(C) __IF_ASYNC__
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:C) IF(OFFLOAD_ON)
    END SUBROUTINE FFTSER_TO_PAR_REAL

#endif // MPI
