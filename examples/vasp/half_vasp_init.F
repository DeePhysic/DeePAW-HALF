#include "symbol.inc"
MODULE half_vasp_init
  USE iso_c_binding, ONLY: c_char,c_double,c_double_complex,c_int,c_int32_t,c_int64_t,c_loc,c_null_char,c_ptr,c_sizeof
  USE prec
  USE poscar_struct_def, ONLY: latt,type_info
  USE wave
  USE wave_mpi, ONLY: DIS_PW_BAND
  USE mpimy
  USE tutor, ONLY: vtutor
  IMPLICIT NONE
  INTEGER(c_int),PARAMETER::HALF_SUCCESS=0,HALF_XC_PBE=2,HALF_BACKEND_AUTO=0,HALF_SOLVER_EVD=1
  TYPE,BIND(C)::half_request_geometry_v1
    INTEGER(c_int32_t)::struct_size,flags,nions,ntypes,grid(3),reserved_i32
    REAL(c_double)::lattice(9)
    TYPE(c_ptr)::species,positions_fractional
    INTEGER(c_int64_t)::reserved(8)
  END TYPE
  INTERFACE
    INTEGER(c_int) FUNCTION half_create_from_files_backend(charge,potential,encut,xc,use_uspp,backend,solver,handle,error,error_capacity) &
        BIND(C,NAME='half_create_from_files_backend')
      IMPORT::c_char,c_double,c_int,c_int64_t
      CHARACTER(c_char),INTENT(IN)::charge(*),potential(*)
      REAL(c_double),VALUE::encut
      INTEGER(c_int),VALUE::xc,use_uspp,backend,solver,error_capacity
      INTEGER(c_int64_t),INTENT(OUT)::handle
      CHARACTER(c_char),INTENT(OUT)::error(*)
    END FUNCTION
    INTEGER(c_int) FUNCTION half_set_request_geometry(handle,geometry,error,error_capacity) &
        BIND(C,NAME='half_set_request_geometry')
      IMPORT::c_char,c_int,c_int64_t,c_ptr
      INTEGER(c_int64_t),VALUE::handle
      TYPE(c_ptr),VALUE::geometry
      CHARACTER(c_char),INTENT(OUT)::error(*)
      INTEGER(c_int),VALUE::error_capacity
    END FUNCTION
    INTEGER(c_int) FUNCTION half_solve_kpoint_mapped(handle,kpoint,nbands,host_npw,host_gvec,eigenvalues,eigenvectors,ld,smin,smax,error,error_capacity) &
        BIND(C,NAME='half_solve_kpoint_mapped')
      IMPORT::c_char,c_double,c_int,c_int32_t,c_int64_t,c_ptr
      INTEGER(c_int64_t),VALUE::handle,host_npw,ld
      INTEGER(c_int),VALUE::nbands,error_capacity
      REAL(c_double),INTENT(IN)::kpoint(3)
      INTEGER(c_int32_t),INTENT(IN)::host_gvec(*)
      REAL(c_double),INTENT(OUT)::eigenvalues(*),smin,smax
      TYPE(c_ptr),VALUE::eigenvectors
      CHARACTER(c_char),INTENT(OUT)::error(*)
    END FUNCTION
    INTEGER(c_int) FUNCTION half_destroy(handle,error,error_capacity) BIND(C,NAME='half_destroy')
      IMPORT::c_char,c_int,c_int64_t
      INTEGER(c_int64_t),VALUE::handle
      CHARACTER(c_char),INTENT(OUT)::error(*)
      INTEGER(c_int),VALUE::error_capacity
    END FUNCTION
  END INTERFACE
CONTAINS
  SUBROUTINE HALF_INITIALIZE_WAVEFUNCTIONS(WDES,W,ENCUT,CHARGE_PATH,POTENTIAL_PATH,LATT_CUR,T_INFO,NGRID,IU0)
    TYPE(wavedes),TARGET,INTENT(INOUT)::WDES
    TYPE(wavespin),INTENT(INOUT)::W
    REAL(q),INTENT(IN)::ENCUT
    CHARACTER(LEN=*),INTENT(IN)::CHARGE_PATH,POTENTIAL_PATH
    TYPE(latt),INTENT(IN)::LATT_CUR
    TYPE(type_info),INTENT(IN)::T_INFO
    INTEGER,INTENT(IN)::NGRID(3)
    INTEGER,INTENT(IN)::IU0
    TYPE(wavedes1)::WDES1
    INTEGER(c_int64_t)::HANDLE
    INTEGER(c_int)::STATUS
    INTEGER::NK,NB,I,J,NC,IND,SERIAL,NPL_TOT,ROOT
    INTEGER,ALLOCATABLE::GVEC_WORK(:)
    INTEGER(c_int32_t),ALLOCATABLE,TARGET::HOST_GVEC(:)
    INTEGER(c_int32_t),ALLOCATABLE,TARGET::REQUEST_SPECIES(:)
    REAL(c_double),ALLOCATABLE,TARGET::EIGENVALUES(:)
    REAL(c_double),ALLOCATABLE,TARGET::REQUEST_POSITIONS(:)
    COMPLEX(c_double_complex),ALLOCATABLE,TARGET::VECTORS(:,:)
    REAL(c_double)::KPOINT(3),SMIN,SMAX
    TYPE(half_request_geometry_v1),TARGET::REQUEST
    CHARACTER(c_char),ALLOCATABLE::CHARGE_C(:),POTENTIAL_C(:)
    CHARACTER(c_char)::ERROR_C(512)

    IF (WDES%ISPIN/=1 .OR. WDES%NRSPINORS/=1 .OR. WDES%LNONCOLLINEAR .OR. WDES%LSORBIT) &
      CALL vtutor%error("LHALF_INIT currently supports only collinear ISPIN=1 without SOC")
    IF (WDES%LGAMMA) CALL vtutor%error("LHALF_INIT requires the full-complex VASP build")
    IF (WDES%COMM_KINTER%NCPU/=1) CALL vtutor%error("LHALF_INIT currently requires KPAR=1")

    ROOT=WDES%COMM%IONODE
    HANDLE=0_c_int64_t
    ERROR_C=c_null_char
    CALL MAKE_C_STRING(TRIM(CHARGE_PATH),CHARGE_C)
    CALL MAKE_C_STRING(TRIM(POTENTIAL_PATH),POTENTIAL_C)
    STATUS=HALF_SUCCESS
    IF (WDES%COMM%NODE_ME==ROOT) THEN
      STATUS=half_create_from_files_backend(CHARGE_C,POTENTIAL_C,REAL(ENCUT,c_double),HALF_XC_PBE,1_c_int, &
        HALF_BACKEND_AUTO,HALF_SOLVER_EVD,HANDLE,ERROR_C,INT(SIZE(ERROR_C),c_int))
      IF (STATUS/=HALF_SUCCESS .AND. IU0>=0) WRITE(IU0,'(A)') 'HALF: '//TRIM(C_ERROR_STRING(ERROR_C))
      IF (STATUS==HALF_SUCCESS) THEN
        ALLOCATE(REQUEST_SPECIES(T_INFO%NIONS),REQUEST_POSITIONS(3*T_INFO%NIONS))
        REQUEST%struct_size=INT(c_sizeof(REQUEST),c_int32_t);REQUEST%flags=0_c_int32_t
        REQUEST%nions=INT(T_INFO%NIONS,c_int32_t);REQUEST%ntypes=INT(T_INFO%NTYP,c_int32_t)
        REQUEST%grid=INT(NGRID,c_int32_t);REQUEST%reserved_i32=0_c_int32_t
        REQUEST%lattice=0.0_c_double;REQUEST%reserved=0_c_int64_t
        DO I=1,3
          DO J=1,3;REQUEST%lattice(3*(I-1)+J)=REAL(LATT_CUR%A(J,I),c_double);ENDDO
        ENDDO
        DO I=1,T_INFO%NIONS
          REQUEST_SPECIES(I)=INT(T_INFO%ITYP(I),c_int32_t)
          DO J=1,3;REQUEST_POSITIONS(3*(I-1)+J)=REAL(T_INFO%POSION(J,I),c_double);ENDDO
        ENDDO
        REQUEST%species=c_loc(REQUEST_SPECIES);REQUEST%positions_fractional=c_loc(REQUEST_POSITIONS)
        STATUS=half_set_request_geometry(HANDLE,c_loc(REQUEST),ERROR_C,INT(SIZE(ERROR_C),c_int))
        IF (STATUS==HALF_SUCCESS .AND. IU0>=0) WRITE(IU0,'(A,I0,A,I0,A,3(I0,1X))') &
          'HALF retained VASP host system: ',T_INFO%NIONS,' ions, ',T_INFO%NTYP,' types, grid ',NGRID
        IF (STATUS/=HALF_SUCCESS .AND. IU0>=0) WRITE(IU0,'(A)') 'HALF: '//TRIM(C_ERROR_STRING(ERROR_C))
        DEALLOCATE(REQUEST_POSITIONS,REQUEST_SPECIES)
      ENDIF
    ENDIF
    CALLMPI(M_bcast_i(WDES%COMM,STATUS,1))
    IF (STATUS/=HALF_SUCCESS) CALL vtutor%error("DeePAW-HALF context creation failed on the root rank")

    ALLOCATE(EIGENVALUES(WDES%NB_TOT))
    DO NK=1,WDES%NKPTS
      CALL SETWDES(WDES,WDES1,NK)
      NPL_TOT=WDES1%NPL_TOT
      ALLOCATE(GVEC_WORK(3*NPL_TOT),HOST_GVEC(3*NPL_TOT),VECTORS(NPL_TOT,WDES%NB_TOT))
      GVEC_WORK=0
      IF (WDES1%COMM_INB%NCPU==1) THEN
        DO I=1,WDES1%NPL
          GVEC_WORK(3*I-2)=WDES1%IGX(I)
          GVEC_WORK(3*I-1)=WDES1%IGY(I)
          GVEC_WORK(3*I  )=WDES1%IGZ(I)
        ENDDO
      ELSE
        IND=0
        DO NC=1,WDES1%NCOL
          DO I=1,WDES1%PL_COL(NC)
            IND=IND+1
            SERIAL=WDES1%PL_INDEX(NC)+I
            GVEC_WORK(3*SERIAL-2)=WDES1%IGX(IND)
            GVEC_WORK(3*SERIAL-1)=WDES1%IGY(IND)
            GVEC_WORK(3*SERIAL  )=WDES1%IGZ(IND)
          ENDDO
        ENDDO
      ENDIF
      CALLMPI(M_sum_i(WDES1%COMM_INB,GVEC_WORK,SIZE(GVEC_WORK)))
      HOST_GVEC=INT(GVEC_WORK,c_int32_t)
      KPOINT=REAL(WDES%VKPT(:,NK),c_double)
      STATUS=HALF_SUCCESS
      IF (WDES%COMM%NODE_ME==ROOT) THEN
        STATUS=half_solve_kpoint_mapped(HANDLE,KPOINT,INT(WDES%NB_TOT,c_int),INT(NPL_TOT,c_int64_t),HOST_GVEC, &
          EIGENVALUES,c_loc(VECTORS),INT(NPL_TOT,c_int64_t),SMIN,SMAX,ERROR_C,INT(SIZE(ERROR_C),c_int))
        IF (STATUS/=HALF_SUCCESS .AND. IU0>=0) WRITE(IU0,'(A)') 'HALF: '//TRIM(C_ERROR_STRING(ERROR_C))
      ENDIF
      CALLMPI(M_bcast_i(WDES%COMM,STATUS,1))
      IF (STATUS/=HALF_SUCCESS) CALL vtutor%error("DeePAW-HALF k-point solve failed on the root rank")
      CALLMPI(M_bcast_d(WDES%COMM,EIGENVALUES,WDES%NB_TOT))
      W%CELTOT(:,NK,1)=CMPLX(EIGENVALUES,0.0_q,q)
      DO NB=1,WDES%NB_TOT
        CALL DIS_PW_BAND(WDES1,NB,VECTORS(:,NB),W%CPTWFP(1,1,NK,1))
      ENDDO
      IF (IU0>=0) WRITE(IU0,'(A,I0,A,I0,A,2ES13.5)') 'HALF initialized k-point ',NK, &
        ' with ',NPL_TOT,' plane waves; overlap range ',SMIN,SMAX
      DEALLOCATE(VECTORS,HOST_GVEC,GVEC_WORK)
    ENDDO

    STATUS=HALF_SUCCESS
    IF (WDES%COMM%NODE_ME==ROOT) STATUS=half_destroy(HANDLE,ERROR_C,INT(SIZE(ERROR_C),c_int))
    CALLMPI(M_bcast_i(WDES%COMM,STATUS,1))
    IF (STATUS/=HALF_SUCCESS) CALL vtutor%error("DeePAW-HALF context destruction failed")
    DEALLOCATE(EIGENVALUES,CHARGE_C,POTENTIAL_C)
  END SUBROUTINE

  SUBROUTINE MAKE_C_STRING(INPUT,OUTPUT)
    CHARACTER(LEN=*),INTENT(IN)::INPUT
    CHARACTER(c_char),ALLOCATABLE,INTENT(OUT)::OUTPUT(:)
    INTEGER::I,N
    N=LEN_TRIM(INPUT);ALLOCATE(OUTPUT(N+1))
    DO I=1,N;OUTPUT(I)=INPUT(I:I);ENDDO
    OUTPUT(N+1)=c_null_char
  END SUBROUTINE

  FUNCTION C_ERROR_STRING(INPUT) RESULT(OUTPUT)
    CHARACTER(c_char),INTENT(IN)::INPUT(:)
    CHARACTER(LEN=511)::OUTPUT
    INTEGER::I
    OUTPUT=' '
    DO I=1,MIN(SIZE(INPUT),LEN(OUTPUT))
      IF (INPUT(I)==c_null_char) EXIT
      OUTPUT(I:I)=INPUT(I)
    ENDDO
  END FUNCTION
END MODULE half_vasp_init
