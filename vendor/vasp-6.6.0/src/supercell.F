#include "symbol.inc"
#if defined(gammareal) && defined(RPAgamma)
#undef RPAgamma
#endif
#if defined(ACC_OFFLOAD) && defined(RPAgamma)
#undef RPAgamma
#endif

MODULE msupercell
  USE base
  USE chi_glb, ONLY: greensfdes
  USE poscar_struct_def, ONLY: type_info, latt
  USE wave_struct_def, ONLY: wavedes, wavedes1
  USE nonlr_struct_def, ONLY: nonlr_struct
  USE nonl_struct_def, ONLY: nonl_struct
  USE mgrid_struct_def, ONLY: grid_3d
  USE string, ONLY: str
  USE tutor, ONLY: vtutor
  USE wave
  USE mpimy
  IMPLICIT NONE
  PRIVATE


  TYPE supercell
     TYPE(type_info) :: T_INFO            !< ionic positions in superlattice
     TYPE (latt)     :: LATT_CUR          !< lattice constants in superlattice
     INTEGER, POINTER:: POS_IN_PRIM(:)=>NULL()    !< position in primitive lattice
     REAL(q), POINTER:: POSION(:,:)=>NULL()       !< position of each ion in fractional coordinates of original cell
     INTEGER         :: NKPX              !< replicate counter in x direction
     INTEGER         :: NKPY              !< replicate counter in y direction
     INTEGER         :: NKPZ              !< replicate counter in z direction
     INTEGER         :: NREP              !< numer of replications of primitive cell
     TYPE (wavedes)  :: WDES              !< wavefunction descriptor in supercell
     TYPE (wavedes1) :: WDES1             !< since we have only a single k-point, store WDES1 as well
     TYPE (wavedes)  :: WGW               !< descriptor for response functions in supercell
     TYPE (wavedes1) :: WGW1              !< since we have only a single k-point, store WDES1 as well
     TYPE (grid_3d)  :: GRID              !< FFT grid for orbitals in supercell
     TYPE (grid_3d)  :: GRID_RES          !< FFT grid for charges or response functions in supercell
              !< in real space the data distribution is strictly identical to GRID
              !< however input and output arrays (to FFT's) 
              !< must be declared as GDEF in real space (as opposed to COMPLEX for FFT's using GRID)
     TYPE (nonlr_struct) :: FAST_AUG      !< fast aug structure for supercell
     REAL(q), POINTER, CONTIGUOUS :: TRANS_MATRIX(:,:,:,:)=>NULL()
     TYPE (wavedes1) :: AUG_DES           !< descriptor for augmenting densities

     INTEGER, POINTER :: NGVECTOR(:)=>NULL()      !< number of indices per k-point (same values as WDES%NGVECTOR)
     !<, used for collecting orbital related quantitites from k-points into supercell
     INTEGER, POINTER :: INDEX(:,:)=>NULL()       !< index array for collecting g vectors for orbital related quantities
     !<, used for collecting orbital related quantitites from k-points into supercell
     INTEGER, POINTER :: INDEX_INV(:,:)=>NULL()   !< index array for collecting g vectors related by inversion symmetry
     !<, used for collecting orbital related quantitites from k-points into supercell

     !> special arrays used if k-points allow for use of inversion symmetry (see also greensfdes in GG_base.F)
     !> at some k-points inversion symmetry can be used to reduce the number of PW coefficients
     !> this is possible since G(g+k,g'+k) = G*(-g-k,-g'-k)
     !> GDES elements LUSEINV and NGVECTOR_INV are linked to entries in GDES (see GG_base.F)
     LOGICAL,POINTER :: LUSEINV(:)=>NULL()        !< specifies that inversion symmetry can be used at this k-point
     INTEGER,POINTER :: NGVECTOR_INV(:)=>NULL()   !< number of actual PW components if inversion symmetry is used
     INTEGER,POINTER :: IND_IF_INV(:,:)=>NULL()   !< index array for collecting g vectors for orbital rel. quantities, if inversion symmetry is used

     INTEGER, POINTER :: NGVECTOR_RES(:)=>NULL()   !< number of indices per k-point
     !<, related to collecting response function related quantities from k-points into supercell
     INTEGER, POINTER :: INDEX_RES(:,:)=>NULL()    !< index array for collecting g vectors for response related quantities
     !<, related to collecting response function related quantities from k-points into supercell
     INTEGER, POINTER :: INDEX_RES_INV(:,:)=>NULL()!< index array for collecting g vectors related by inversion symmetry
     !<, related to collecting response function related quantities from k-points into supercell
          !< primitive cell into supercell arrays (for each k-point in primitive cell)
          !< this one is for response function like quantities
     TYPE (nonl_struct),POINTER  :: NONL_S=>NULL()!< debugging only, required to obtain CPROJ in supercell
     INTEGER, POINTER    :: KWEIGHT(:)=>NULL()    !< weight of individual k-points (if ISYM=0 is used in GW)
  END TYPE supercell

  PUBLIC :: supercell
  PUBLIC :: INIT_SUPERCELL
  PUBLIC :: INIT_SUPERCELL_FFT_WF
  PUBLIC :: INIT_SUPERCELL_FFT
  PUBLIC :: DEALLOCATE_SUPERCELL

CONTAINS

!*********************************************************************
!
!> the subroutine INIT_SUPERCELL initializes the supercell structure 
!> at this point the simplest possible algorithm is chosen
!> A supercell is created by multiplying the original cell with
!> NKPX, Y and Z in the three directions
!> then the replicated ions are created 
!
!*********************************************************************

  SUBROUTINE INIT_SUPERCELL( S, T_INFO, LATT_CUR, NKPX, NKPY, NKPZ, B, IU6) 
    USE lattice
    TYPE (supercell), POINTER :: S !< super cell handle set by this routine
    TYPE(type_info) :: T_INFO     !< ionic positions in original cell
    TYPE (latt)     :: LATT_CUR   !< lattice constants of original cell
    INTEGER         :: NKPX       !< integer division along rec. lattice 1-vector for generating the mesh
    INTEGER         :: NKPY       !< integer division along rec. lattice 2-vector for generating the mesh
    INTEGER         :: NKPZ       !< integer division along rec. lattice 3-vector for generating the mesh
    REAL(q)         :: B(3,3)     !< alternatively generating lattice for k-point mesh
    INTEGER         :: IU6        !< IO unit for output

  ! local
    INTEGER :: N, M, N1, N2, N3, NI, NI_NEW
    INTEGER         :: IA(3,3)    ! supercell in terms of original cell in integer units
    REAL(q)         :: TMP(3,3)   ! temporary vector
    REAL(q), ALLOCATABLE :: TRANS(:,:)

    IF (ASSOCIATED(S)) &
       CALL vtutor%bug("INIT_SUPERCELL: S is allready associated", __FILE__, __LINE__)

    ALLOCATE(S)

    S%NKPX=NKPX
    S%NKPY=NKPY
    S%NKPZ=NKPZ

    IF (NKPX/=-1 .AND. NKPY/=-1 .AND. NKPZ/=-1) THEN
!
! (TODO: the supercell could be Minkovski reduced to yield near
!  orthogonal vectors)
      S%LATT_CUR%A(:,1)=LATT_CUR%A(:,1)*NKPX
      S%LATT_CUR%A(:,2)=LATT_CUR%A(:,2)*NKPY
      S%LATT_CUR%A(:,3)=LATT_CUR%A(:,3)*NKPZ

      S%NREP=NKPX*NKPY*NKPZ
      ALLOCATE(TRANS(3,S%NREP))

      CALL LATTIC(S%LATT_CUR)
      N=0
      DO N1=0,NKPX-1
         DO N2=0,NKPY-1
            DO N3=0,NKPZ-1
               N=N+1
               TRANS(1,N)=N1
               TRANS(2,N)=N2
               TRANS(3,N)=N3
            ENDDO
         ENDDO
      ENDDO
    ELSE
      ! set unit cell to generating lattice of k-point grid	
      S%LATT_CUR%A=B
      ! determine reciprocal lattice of generating lattice -> B
      CALL LATTIC(S%LATT_CUR)

      ! TMP now stores the supercell
      TMP=S%LATT_CUR%B

      CALL KARDIR(3, TMP, LATT_CUR%B)
      IA=NINT(TMP)
      IF (IU6>=0) THEN
        WRITE(IU6,*)
        WRITE(IU6,'(A)') ' Supercell in multiplies of origin unit cell:'
        WRITE(IU6,'(A)') ' ============================================'
        WRITE(IU6,'(2X,3I5)') IA
      ENDIF
      ! just in case recalculate A to reduce user rounding errors
      
      TMP=IA
      CALL DIRKAR(3, TMP, LATT_CUR%A)

      S%LATT_CUR%A=TMP
      CALL LATTIC(S%LATT_CUR)

      ! now determine generating basis in terms of original lattice vectors
      S%NREP=NINT(S%LATT_CUR%OMEGA/LATT_CUR%OMEGA)
      ALLOCATE(TRANS(3,S%NREP))
      N=0
      ! brute force search for integer Bravais points that are within the supercell 
      DO N1=0,S%NREP
         DO N2=0,S%NREP
            DO N3=0,S%NREP
               TMP(1,1)=N1
               TMP(2,1)=N2
               TMP(3,1)=N3
               ! bring to cartesian coordinates
               CALL DIRKAR(1, TMP, LATT_CUR%A)
               ! convert to direct coordinates in super cell units
               CALL KARDIR(1, TMP, S%LATT_CUR%B)
               ! bring to Wigner Seitz  super cell
               TMP(1,1)=MODULO(TMP(1,1),1.0_q)
               TMP(2,1)=MODULO(TMP(2,1),1.0_q)
               TMP(3,1)=MODULO(TMP(3,1),1.0_q)

               ! rounding issues, very close to one, make it zero 
               IF (ABS(TMP(1,1)-1.0_q)<1E-4_q) TMP(1,1)=0
               IF (ABS(TMP(2,1)-1.0_q)<1E-4_q) TMP(2,1)=0
               IF (ABS(TMP(3,1)-1.0_q)<1E-4_q) TMP(3,1)=0

               ! cartesian coordinates
               CALL DIRKAR(1, TMP, S%LATT_CUR%A)
               ! direct coordinates of original primitive cell
               CALL KARDIR(1, TMP, LATT_CUR%B)
               IF (ABS(TMP(1,1)-NINT(TMP(1,1)))>0.1  .OR. & 
                   ABS(TMP(2,1)-NINT(TMP(2,1)))>0.1  .OR. &
                   ABS(TMP(3,1)-NINT(TMP(3,1)))>0.1) THEN
                  CALL vtutor%bug("INIT_SUPERCELL: not integer multiples " // &
                     str(TMP(:,1)), __FILE__, __LINE__)
               ENDIF
               TMP(:,1)=NINT(TMP(:,1))
               
               ! compare to yet found positions
               DO M=1,N
                  IF ((TMP(1,1)-TRANS(1,M))==0 .AND. & 
                      (TMP(2,1)-TRANS(2,M))==0 .AND. &
                      (TMP(3,1)-TRANS(3,M))==0) EXIT
               ENDDO
               IF (M>N) THEN
                  N=N+1
                  IF (N>S%NREP) THEN
                     CALL vtutor%bug(str(TRANS) // "\n INIT_SUPERCELL: too many &
                        &replicated positions " // str(N) // " " // str(S%NREP) // " " // str(TMP(:,1)), &
                        __FILE__, __LINE__)
                  ENDIF
                  TRANS(1,N)=TMP(1,1)
                  TRANS(2,N)=TMP(2,1)
                  TRANS(3,N)=TMP(3,1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
    ENDIF

    IF (N/=S%NREP) THEN
       CALL vtutor%bug(str(TRANS) // "\n INIT_SUPERCELL: not all replicated " &
          // "positions found " // str(N) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF

    N=S%NREP

    S%T_INFO%NTYPD =T_INFO%NTYPD       ! dimension for types
    S%T_INFO%NTYP  =T_INFO%NTYP        ! number of types
    S%T_INFO%NTYPPD=T_INFO%NTYPD       ! dimension for types inc. empty spheres
    S%T_INFO%NTYPP =T_INFO%NTYP        ! number of types empty spheres
    S%T_INFO%NIOND =N*T_INFO%NIOND     ! dimension for ions
    S%T_INFO%NIONPD=N*T_INFO%NIOND     ! dimension for ions inc. empty spheres
    S%T_INFO%NIONS =N*T_INFO%NIONS     ! actual number of ions
    S%T_INFO%NIONP =N*T_INFO%NIONS     ! actual number of ions inc. empty spheres
    S%T_INFO%LSDYN =T_INFO%LSDYN       ! selective dynamics (yes/ no)
    S%T_INFO%LDIRCO=T_INFO%LDIRCO      ! positions in direct/recproc. lattice

#ifdef VASP_HDF5
    IF(ASSOCIATED(T_INFO%SHA256)) THEN
      ALLOCATE(S%T_INFO%SHA256(SIZE(T_INFO%SHA256)))
      S%T_INFO%SHA256 = T_INFO%SHA256
    ENDIF
#endif


    N=S%T_INFO%NTYPD

    ALLOCATE(S%T_INFO%NITYP(N), S%T_INFO%TYPE(N), S%T_INFO%VCA(N))
    S%T_INFO%NITYP(1:N)=T_INFO%NITYP(1:N)*S%NREP  ! number of  ions per type
    S%T_INFO%TYPE(1:N) =T_INFO%TYPE(1:N)   ! type information
    S%T_INFO%VCA(1:N)  =T_INFO%VCA(1:N)

    N=S%T_INFO%NIOND

    ALLOCATE(S%T_INFO%POSION(3,N), S%T_INFO%LSFOR(3,N), S%T_INFO%ITYP(N))
    ALLOCATE(S%POS_IN_PRIM(N), S%POSION(3,N))
    NI_NEW=0
    DO NI=1,T_INFO%NIONS
       ! maybe we should chose a minimum image convention here
       ! but results should be entirely independent of this choice
       DO M=1,S%NREP
          NI_NEW=NI_NEW+1
          ! position in fractional coordinates of original cell
          S%POSION(1,NI_NEW)=T_INFO%POSION(1,NI)+TRANS(1,M)
          S%POSION(2,NI_NEW)=T_INFO%POSION(2,NI)+TRANS(2,M)
          S%POSION(3,NI_NEW)=T_INFO%POSION(3,NI)+TRANS(3,M)
          TMP(:,1)=S%POSION(:,NI_NEW)
          ! cartesian coordinates and then fractional coordinates of supercell
          CALL DIRKAR(1, TMP, LATT_CUR%A)
          CALL KARDIR(1, TMP, S%LATT_CUR%B)

          ! bring to Wigner Seitz  super cell
!          TMP(1,1)=MODULO(TMP(1,1),1.0_q)
!          TMP(2,1)=MODULO(TMP(2,1),1.0_q)
!          TMP(3,1)=MODULO(TMP(3,1),1.0_q)
          
          ! position in fractional coordinates of super cell
          S%T_INFO%POSION(:,NI_NEW)=TMP(:,1)
          S%T_INFO%LSFOR(:,NI_NEW) =T_INFO%LSFOR(:,NI)
          S%T_INFO%ITYP(NI_NEW)    =T_INFO%ITYP(NI)
          S%POS_IN_PRIM(NI_NEW)=NI

          ! convert from supercell coordinates to unit cell coordinates
!          CALL DIRKAR(1, TMP, S%LATT_CUR%A)
!          CALL KARDIR(1, TMP, LATT_CUR%B)
!          S%POSION(:,NI_NEW)=TMP(:,1)
       ENDDO
    ENDDO
    ! pretty sure none of the remaining elements are important
    ! nullify so that we stop when something is missing
    NULLIFY( S%T_INFO%POMASS, S%T_INFO%RWIGS, S%T_INFO%ROPT, S%T_INFO%ATOMOM, &
             S%T_INFO%DARWIN_R, S%T_INFO%DARWIN_V, S%T_INFO%ZCT, S%T_INFO%RGAUS)

    DEALLOCATE(TRANS)

  END SUBROUTINE INIT_SUPERCELL

!*********************************************************************
!
!> initializes the FFT's for the supercell only for wavefunctions
!
!*********************************************************************

  SUBROUTINE INIT_SUPERCELL_FFT_WF( S, P, LMDIM, WDES, LATT_CUR, IU6)
    USE pseudo
    USE pseudo_struct_def
    USE fock
    
    TYPE (supercell), POINTER :: S !< super cell handle set by this routine
    TYPE (potcar),TARGET::  P(:) !< potcar handle
    INTEGER        :: LMDIM      !< leading dimension of arrays like CDIJ
    TYPE (wavedes) :: WDES       !< original orbital descriptor for primitive cell
    TYPE (wavedes1):: WDES1      !< original orbital descriptor for primitive cell
    TYPE (latt)    :: LATT_CUR   !< lattice constants of original cell
    INTEGER        :: IU6        !< IO unit for output
  ! local
    REAL(q) :: G1, G2, G3, F1, F2, F3, SPACING
    INTEGER :: M, NK, IG1, IG2, IG3, NG, M_INV
    INTEGER, ALLOCATABLE ::  INDEX_INTO_SUPERCELL(:,:,:)
    INTEGER, ALLOCATABLE ::  KFOUND(:,:,:)
    LOGICAL :: LINV
    REAL(q), PARAMETER :: G2ZERO=1E-12_q

    S%WDES=WDES
    ! set number of k-points to 1
    S%WDES%NKDIM=1
    S%WDES%NKPTS=1
    S%WDES%NKPTS_FOR_GEN_LAYOUT=1

    CALL WDES_SET_NPRO(S%WDES, S%T_INFO, P, WDES%LOVERL)
 
    ! set the single k-point to Gamma
    NULLIFY(S%WDES%VKPT, S%WDES%WTKPT)
    ALLOCATE(S%WDES%VKPT(3,1), S%WDES%WTKPT(1))
    S%WDES%VKPT =0
    S%WDES%WTKPT=1

    S%GRID=WDES%GRID
    IF (S%NKPX/=-1 .AND. S%NKPY/=-1 .AND. S%NKPZ/=-1) THEN
       S%GRID%NGPTAR(1)=WDES%GRID%NGPTAR(1)*S%NKPX
       S%GRID%NGPTAR(2)=WDES%GRID%NGPTAR(2)*S%NKPY
       S%GRID%NGPTAR(3)=WDES%GRID%NGPTAR(3)*S%NKPZ
    ELSE
       ! try to get same spacing as for conventional cell
       ! unfortunately the mapping is not always exact
       ! hence the FFT grids will be necessarily incompatible
       SPACING=(WDES%GRID%NGPTAR(1)/LATT_CUR%ANORM(1)+ & 
                WDES%GRID%NGPTAR(2)/LATT_CUR%ANORM(2)+ &
                WDES%GRID%NGPTAR(3)/LATT_CUR%ANORM(3))/3.0_q
       S%GRID%NGPTAR(1)=SPACING*S%LATT_CUR%ANORM(1)+0.99
       S%GRID%NGPTAR(2)=SPACING*S%LATT_CUR%ANORM(2)+0.99
       S%GRID%NGPTAR(3)=SPACING*S%LATT_CUR%ANORM(3)+0.99
    ENDIF

    CALL FFTCHK(S%GRID%NGPTAR)
    ! reinitialize the loop counters
    CALL INILGRD(S%GRID%NGPTAR(1),S%GRID%NGPTAR(2),S%GRID%NGPTAR(3),S%GRID)

    IF (IU6>=0) THEN
       WRITE(IU6,20) S%GRID%NGPTAR
20     FORMAT(/' FFT grid for supercell:   NGX =',I4,'; NGY =',I4,'; NGZ =',I4,/ & 
               ' =======================',/)
    ENDIF

    ! set S%WDES%GRID to S%GRID and generate the data layout for FFT's
    CALL GEN_LAYOUT(S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, .TRUE.)
    CALL GEN_INDEX (S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, -1, .TRUE.)

    CALLMPI( MAPSET(S%GRID))  ! generate the communication maps (patterns) for FFT
    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(S%WDES%NINDPW(1,1),S%WDES%NGVECTOR(1),S%WDES%NKPTS,S%WDES%NGDIM,S%GRID)
    ! create (and destroy) the FFT plans to generate "wisdom"
    CALL FFTGRIDPLAN(S%GRID)

    ! set up the structure (FAST_AUG_FOCK) to perform the fast augmentation 
    ! of the charge density

    IF (S%WDES%LOVERL .AND. LMAX_FOCK>=0 .AND. LHFCALC) THEN
       ! note FAST_AUG_FOCK, LMAX_FOCKAE, and NMAX_FOCKAE are global in fock.F
       ! imported via USE
       CALL FASTAUG_SUPER( S%T_INFO, P, S%WDES, S%GRID, S%LATT_CUR,  LMDIM, & 
            FAST_AUG_FOCK, S%FAST_AUG, S%TRANS_MATRIX, LMAX_FOCKAE, NMAX_FOCKAE )
       CALL RSPHER(S%GRID, S%FAST_AUG, S%LATT_CUR)
       S%FAST_AUG%RPROJ= S%FAST_AUG%RPROJ*SQRT(S%LATT_CUR%OMEGA)
       ! set up AUG_DES describing data layout of one center terms
       CALL SETUP_AUG_DES(S%GRID, S%WDES, S%AUG_DES, S%FAST_AUG )
    ELSE
       NULLIFY(S%TRANS_MATRIX)
       ALLOCATE(S%AUG_DES%LMMAX(SIZE(S%WDES%LMMAX)))
       ! to make live easier set LMMAX (number of channels) to zero for each type
       ! also set NPRO, NPROD etc. to zero
       S%AUG_DES%LMMAX=0
       S%AUG_DES%NPRO =0
       S%AUG_DES%NPROD=0
       S%AUG_DES%NPRO_TOT=0
    ENDIF

!=====================================================================
! set up the index array to go from original orbitals in primitive
! cell to supercell indexing
!=====================================================================
    ! allocate a "huge" array to track the reciprocal lattice vectors in the supercell
    CALL SETWDES(S%WDES, S%WDES1, 1)
    NG=S%WDES1%NGVECTOR

    ALLOCATE(INDEX_INTO_SUPERCELL( & 
         MINVAL(S%WDES1%IGX(1:NG)):MAXVAL(S%WDES1%IGX(1:NG)), &
         MINVAL(S%WDES1%IGY(1:NG)):MAXVAL(S%WDES1%IGY(1:NG)), &
         MINVAL(S%WDES1%IGZ(1:NG)):MAXVAL(S%WDES1%IGZ(1:NG))))
    ALLOCATE(KFOUND( & 
         MINVAL(S%WDES1%IGX(1:NG)):MAXVAL(S%WDES1%IGX(1:NG)), &
         MINVAL(S%WDES1%IGY(1:NG)):MAXVAL(S%WDES1%IGY(1:NG)), &
         MINVAL(S%WDES1%IGZ(1:NG)):MAXVAL(S%WDES1%IGZ(1:NG))))
    ! setup the indices in that array
    INDEX_INTO_SUPERCELL=0  ! clear all entries
    KFOUND=0
    DO M=1,S%WDES1%NGVECTOR
       ! store the index M at position G_x, G_y, G_z
       INDEX_INTO_SUPERCELL(S%WDES1%IGX(M),S%WDES1%IGY(M), S%WDES1%IGZ(M))=M
    ENDDO
    ALLOCATE(S%INDEX(WDES%NGDIM, WDES%NKPTS), S%NGVECTOR(WDES%NKPTS))
    ALLOCATE(S%INDEX_INV(WDES%NGDIM, WDES%NKPTS), S%IND_IF_INV(WDES%NGDIM, WDES%NKPTS))
    S%INDEX=0
    S%INDEX_INV=0
    S%IND_IF_INV=0

    DO NK=1,WDES%NKPTS
       CALL SETWDES(WDES, WDES1, NK)
       S%NGVECTOR(NK)=WDES1%NGVECTOR
       DO M=1,WDES1%NGVECTOR
          ! G vector in original cell
          G1=WDES1%IGX(M)+WDES1%VKPT(1)
          G2=WDES1%IGY(M)+WDES1%VKPT(2)
          G3=WDES1%IGZ(M)+WDES1%VKPT(3)
          ! cartesian coordinates
          F1=G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)
          F2=G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)
          F3=G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)
          ! supercell reciprocal coordinates
          G1=F1*S%LATT_CUR%A(1,1)+F2*S%LATT_CUR%A(2,1)+F3*S%LATT_CUR%A(3,1)
          G2=F1*S%LATT_CUR%A(1,2)+F2*S%LATT_CUR%A(2,2)+F3*S%LATT_CUR%A(3,2)
          G3=F1*S%LATT_CUR%A(1,3)+F2*S%LATT_CUR%A(2,3)+F3*S%LATT_CUR%A(3,3)

          IG1=NINT(G1)
          IG2=NINT(G2)
          IG3=NINT(G3)
          IF( ABS(G1-IG1)>1E-6 .OR. ABS(G2-IG2)>1E-6 .OR. ABS(G3-IG3)>1E-6) THEN
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT_WF: can not map g+k to supercell: &
                &" // str(G1) // " " // str(IG1) // " " // str(G2) // " " // str(IG2) // " " // &
                str(G3) // " " // str(IG3), __FILE__, __LINE__)
          ENDIF
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 .OR. KFOUND(IG1, IG2, IG3)>0 ) THEN
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT_WF: can not find g+k: " // str(IG1) &
                // " " // str(IG2) // " " // str(IG3) // " " // str(INDEX_INTO_SUPERCELL(IG1,IG2,IG3)) &
                // " " // str(KFOUND(IG1,IG2,IG3)), __FILE__, __LINE__)
          ENDIF
          S%INDEX(M, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
          ! a one to one mapping must exist, so clear the one we have already found
          KFOUND(IG1, IG2, IG3)=NK ! set k-point where it was found
       ENDDO
       ! store number of required G vectors, if inversion symmetry is used
    ENDDO

    ! check whether all G vectors in supercell have been found
    DO M=1,S%WDES1%NGVECTOR
       IF (KFOUND(S%WDES1%IGX(M),S%WDES1%IGY(M), S%WDES1%IGZ(M))==0) THEN
          CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT_WF: can not map G to primitive cell: " &
             // str(S%WDES1%IGX(M)) // " " // str(S%WDES1%IGY(M)) // " " // str(S%WDES1%IGZ(M)), __FILE__, __LINE__)
       ENDIF
    ENDDO

    DEALLOCATE(INDEX_INTO_SUPERCELL, KFOUND)
!=====================================================================
! finally set the weight for the projectors when they are
! collected from primitive to supercell
!=====================================================================
    ALLOCATE(S%KWEIGHT(WDES%NKPTS))
    DO NK=1,WDES%NKPTS
       S%KWEIGHT(NK)=WDES%WTKPT(NK)*S%NREP
    ENDDO
!=====================================================================
! for debugging it is conventient to get set NONL_S structure
! so that the CPROJ array can be calculated directly in the supercell
!=====================================================================
    NULLIFY(S%NONL_S)

  END SUBROUTINE INIT_SUPERCELL_FFT_WF

!*********************************************************************
!
!> initializes the FFT's for the supercell
!
!*********************************************************************

  SUBROUTINE INIT_SUPERCELL_FFT( S, P, LMDIM, WDES, WGW, GDES, LATT_CUR, ENCUTGW, IU6)
    USE pseudo
    USE pseudo_struct_def
    USE FOCK
    
    TYPE (supercell), POINTER :: S !< super cell handle set by this routine
    TYPE (potcar),TARGET::  P(:) !< potcar handle
    INTEGER        :: LMDIM      !< leading dimension of arrays like CDIJ
    TYPE (wavedes) :: WDES       !< original orbital descriptor for primitive cell
    TYPE (wavedes1):: WDES1      !< original orbital descriptor for primitive cell
    TYPE (wavedes) :: WGW        !< original orbital response descriptor for primitive cell
    TYPE (greensfdes) :: GDES    !< descriptor for distributed Green's function & polarization functions
    TYPE (latt)    :: LATT_CUR   !< lattice constants of original cell
    REAL(q)        :: ENCUTGW    !< cutoff energy 
    INTEGER        :: IU6        !< IO unit for output
  ! local
    REAL(q) :: G1, G2, G3, F1, F2, F3, SPACING
    INTEGER :: M, NK, IG1, IG2, IG3, NG, M_INV
    INTEGER, ALLOCATABLE ::  INDEX_INTO_SUPERCELL(:,:,:)
    INTEGER, ALLOCATABLE ::  KFOUND(:,:,:)
    LOGICAL :: LINV
    REAL(q), PARAMETER :: G2ZERO=1E-12_q

    S%WDES=WDES
    ! set number of k-points to 1
    S%WDES%NKDIM=1
    S%WDES%NKPTS=1
    S%WDES%NKPTS_FOR_GEN_LAYOUT=1

    CALL WDES_SET_NPRO(S%WDES, S%T_INFO, P, WDES%LOVERL)
 
    ! set the single k-point to Gamma
    NULLIFY(S%WDES%VKPT, S%WDES%WTKPT)
    ALLOCATE(S%WDES%VKPT(3,1), S%WDES%WTKPT(1))
    S%WDES%VKPT =0
    S%WDES%WTKPT=1

    S%GRID=WDES%GRID
    IF (S%NKPX/=-1 .AND. S%NKPY/=-1 .AND. S%NKPZ/=-1) THEN
       S%GRID%NGPTAR(1)=WDES%GRID%NGPTAR(1)*S%NKPX
       S%GRID%NGPTAR(2)=WDES%GRID%NGPTAR(2)*S%NKPY
       S%GRID%NGPTAR(3)=WDES%GRID%NGPTAR(3)*S%NKPZ
    ELSE
       ! try to get same spacing as for conventional cell
       ! unfortunately the mapping is not always exact
       ! hence the FFT grids will be necessarily incompatible
       SPACING=(WDES%GRID%NGPTAR(1)/LATT_CUR%ANORM(1)+ & 
                WDES%GRID%NGPTAR(2)/LATT_CUR%ANORM(2)+ &
                WDES%GRID%NGPTAR(3)/LATT_CUR%ANORM(3))/3.0_q
       S%GRID%NGPTAR(1)=SPACING*S%LATT_CUR%ANORM(1)+0.99
       S%GRID%NGPTAR(2)=SPACING*S%LATT_CUR%ANORM(2)+0.99
       S%GRID%NGPTAR(3)=SPACING*S%LATT_CUR%ANORM(3)+0.99
    ENDIF

    CALL FFTCHK(S%GRID%NGPTAR)
    ! reinitialize the loop counters
    CALL INILGRD(S%GRID%NGPTAR(1),S%GRID%NGPTAR(2),S%GRID%NGPTAR(3),S%GRID)

    IF (IU6>=0) THEN
       WRITE(IU6,20) S%GRID%NGPTAR
20     FORMAT(/' FFT grid for supercell:   NGX =',I4,'; NGY =',I4,'; NGZ =',I4,/ & 
               ' =======================',/)
    ENDIF

    ! set S%WDES%GRID to S%GRID and generate the data layout for FFT's
#ifdef RPAgamma
    CALL GEN_LAYOUT(S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, .TRUE., LNOGAMMA=.FALSE.)
#else
    CALL GEN_LAYOUT(S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, .TRUE.)
#endif
    CALL GEN_INDEX (S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, -1, .TRUE.)

    CALLMPI( MAPSET(S%GRID))  ! generate the communication maps (patterns) for FFT
    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(S%WDES%NINDPW(1,1),S%WDES%NGVECTOR(1),S%WDES%NKPTS,S%WDES%NGDIM,S%GRID)
    ! create (and destroy) the FFT plans to generate "wisdom"
    CALL FFTGRIDPLAN(S%GRID)

    ! set up the structure (FAST_AUG_FOCK) to perform the fast augmentation 
    ! of the charge density

    IF (S%WDES%LOVERL .AND. LMAX_FOCK>=0) THEN
       ! note FAST_AUG_FOCK, LMAX_FOCKAE, and NMAX_FOCKAE are global in fock.F
       ! imported via USE
       CALL FASTAUG_SUPER( S%T_INFO, P, S%WDES, S%GRID, S%LATT_CUR,  LMDIM, & 
            FAST_AUG_FOCK, S%FAST_AUG, S%TRANS_MATRIX, LMAX_FOCKAE, NMAX_FOCKAE )
       CALL RSPHER(S%GRID, S%FAST_AUG, S%LATT_CUR)
#ifdef shmem_rproj
       IF (S%FAST_AUG%COMM_shmem%NODE_ME==1) S%FAST_AUG%RPROJ= S%FAST_AUG%RPROJ*SQRT(S%LATT_CUR%OMEGA)
       CALLMPI( M_barrier(S%FAST_AUG%COMM_shmem))
#else
       S%FAST_AUG%RPROJ= S%FAST_AUG%RPROJ*SQRT(S%LATT_CUR%OMEGA)
#endif
       ! set up AUG_DES describing data layout of one center terms
       CALL SETUP_AUG_DES(S%GRID, S%WDES, S%AUG_DES, S%FAST_AUG )
    ELSE
       NULLIFY(S%TRANS_MATRIX)
       ALLOCATE(S%AUG_DES%LMMAX(SIZE(S%WDES%LMMAX)))
       ! to make live easier set LMMAX (number of channels) to zero for each type
       ! also set NPRO, NPROD etc. to zero
       S%AUG_DES%LMMAX=0
       S%AUG_DES%NPRO =0
       S%AUG_DES%NPROD=0
       S%AUG_DES%NPRO_TOT=0
    ENDIF
!
! WGW descriptor is used for FFT of overlap density related quantities
!
    S%WGW=S%WDES
    S%WGW%ENMAX=ENCUTGW
    S%GRID_RES=S%GRID
    IF (IU6>=0) THEN 
       WRITE(IU6,*)
       WRITE(IU6,*) 'Basis sets for responsefunction in the supercell:'
       WRITE(IU6,*) '================================================='
    ENDIF
#ifdef RPAgamma
    CALL GEN_LAYOUT(S%GRID_RES, S%WGW, S%LATT_CUR%B, S%LATT_CUR%B, IU6,.TRUE.,LNOGAMMA=.FALSE., LREAL=.TRUE.)
#else
    CALL GEN_LAYOUT(S%GRID_RES, S%WGW, S%LATT_CUR%B, S%LATT_CUR%B, IU6,.TRUE.)
#endif
    CALL GEN_INDEX (S%GRID_RES, S%WGW, S%LATT_CUR%B, S%LATT_CUR%B, IU6, -1, .TRUE.)

    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(S%WGW%NINDPW(1,1), S%WGW%NGVECTOR(1), S%WGW%NKPTS, S%WGW%NGDIM, S%GRID_RES)

#ifdef RPAgamma
    CALL GDES_ALLOC_K(GDES, WDES)
    S%LUSEINV     =>GDES%LUSEINV
    S%NGVECTOR_INV=>GDES%NGVECTOR_INV
#else
    NULLIFY(S%LUSEINV, S%NGVECTOR_INV)
#endif
!=====================================================================
! set up the index array to go from original orbitals in primitive
! cell to supercell indexing
!=====================================================================
    ! allocate a "huge" array to track the reciprocal lattice vectors in the supercell
    CALL SETWDES(S%WDES, S%WDES1, 1)
    NG=S%WDES1%NGVECTOR

    ALLOCATE(INDEX_INTO_SUPERCELL( & 
         MINVAL(S%WDES1%IGX(1:NG)):MAXVAL(S%WDES1%IGX(1:NG)), &
         MINVAL(S%WDES1%IGY(1:NG)):MAXVAL(S%WDES1%IGY(1:NG)), &
         MINVAL(S%WDES1%IGZ(1:NG)):MAXVAL(S%WDES1%IGZ(1:NG))))
    ALLOCATE(KFOUND( & 
         MINVAL(S%WDES1%IGX(1:NG)):MAXVAL(S%WDES1%IGX(1:NG)), &
         MINVAL(S%WDES1%IGY(1:NG)):MAXVAL(S%WDES1%IGY(1:NG)), &
         MINVAL(S%WDES1%IGZ(1:NG)):MAXVAL(S%WDES1%IGZ(1:NG))))
    ! setup the indices in that array
    INDEX_INTO_SUPERCELL=0  ! clear all entries
    KFOUND=0
    DO M=1,S%WDES1%NGVECTOR
       ! store the index M at position G_x, G_y, G_z
       INDEX_INTO_SUPERCELL(S%WDES1%IGX(M),S%WDES1%IGY(M), S%WDES1%IGZ(M))=M
    ENDDO
    ALLOCATE(S%INDEX(WDES%NGDIM, WDES%NKPTS), S%NGVECTOR(WDES%NKPTS))
    ALLOCATE(S%INDEX_INV(WDES%NGDIM, WDES%NKPTS), S%IND_IF_INV(WDES%NGDIM, WDES%NKPTS))
    S%INDEX=0
    S%INDEX_INV=0
    S%IND_IF_INV=0

    DO NK=1,WDES%NKPTS
#ifdef RPAgamma
       ! determine whether we can reduce the number of PW coefficients for this k-point
       ! this reduction will be used for G(g,g'), G(g,r'), G(r',g)
       IF (ABS(MOD(WDES%VKPT(1,NK)*2+100,1._q))> SQRT(G2ZERO) .OR. &
           ABS(MOD(WDES%VKPT(2,NK)*2+100,1._q))> SQRT(G2ZERO) .OR. &
           ABS(MOD(WDES%VKPT(3,NK)*2+100,1._q))> SQRT(G2ZERO)) THEN
          S%LUSEINV(NK)=.FALSE.
       ELSE
          S%LUSEINV(NK)=.TRUE.
       ENDIF
       M_INV=0    ! counter for the number of PW after use of inversion symmetry
#endif
       CALL SETWDES(WDES, WDES1, NK)
       S%NGVECTOR(NK)=WDES1%NGVECTOR
       DO M=1,WDES1%NGVECTOR
          ! G vector in original cell
          G1=WDES1%IGX(M)+WDES1%VKPT(1)
          G2=WDES1%IGY(M)+WDES1%VKPT(2)
          G3=WDES1%IGZ(M)+WDES1%VKPT(3)
          ! cartesian coordinates
          F1=G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)
          F2=G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)
          F3=G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)
          ! supercell reciprocal coordinates
          G1=F1*S%LATT_CUR%A(1,1)+F2*S%LATT_CUR%A(2,1)+F3*S%LATT_CUR%A(3,1)
          G2=F1*S%LATT_CUR%A(1,2)+F2*S%LATT_CUR%A(2,2)+F3*S%LATT_CUR%A(3,2)
          G3=F1*S%LATT_CUR%A(1,3)+F2*S%LATT_CUR%A(2,3)+F3*S%LATT_CUR%A(3,3)

          IG1=NINT(G1)
          IG2=NINT(G2)
          IG3=NINT(G3)
          IF( ABS(G1-IG1)>1E-6 .OR. ABS(G2-IG2)>1E-6 .OR. ABS(G3-IG3)>1E-6) THEN
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not map g+k to supercell: &
                &" // str(G1) // " " // str(IG1) // " " // str(G2) // " " // str(IG2) // " " // &
                str(G3) // " " // str(IG3), __FILE__, __LINE__)
          ENDIF
#ifdef RPAgamma
          ! skip wave vectors that are not considered in half grid mode
          LINV=.FALSE.
          IF (IG1<0) THEN
             LINV=.TRUE.
          ELSE IF (IG1==0) THEN
             IF (IG2<0) THEN
                LINV=.TRUE.
             ELSE IF (IG2==0) THEN
                IF (IG3<0) THEN
                   LINV=.TRUE.
                ENDIF
             ENDIF
          ENDIF
          IF (LINV) THEN
             IF (INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)<=0 ) THEN
                CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not find inv g+k: " // &
                   str(-IG1) // " " // str(-IG2) // " " // str(-IG3) // " " // &
                   str(INDEX_INTO_SUPERCELL(-IG1,-IG2,-IG3)), __FILE__, __LINE__)
             ENDIF
             S%INDEX_INV(M, NK)=INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)             
             KFOUND(-IG1, -IG2, -IG3)=NK ! set k-point where it was found
             CYCLE
          ELSE
             ! Green's function G_k(g,g'), G_k(r,g'), etc. are stored applying inversion symmetry
             IF (S%LUSEINV(NK)) THEN
                M_INV=M_INV+1
                S%IND_IF_INV(M_INV, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
                GDES%MAP_TO_FULL(M_INV,1,NK)=M
             ENDIF
          ENDIF
          ! test whether this G vector has already been hit does not work for RPAgamma
          ! since some vectors are found twice (once for the conjugated plane wave)
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 ) THEN
#else
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 .OR. KFOUND(IG1, IG2, IG3)>0 ) THEN
#endif
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not find g+k: " // str(IG1) &
                // " " // str(IG2) // " " // str(IG3) // " " // str(INDEX_INTO_SUPERCELL(IG1,IG2,IG3)) &
                // " " // str(KFOUND(IG1,IG2,IG3)), __FILE__, __LINE__)
          ENDIF
          S%INDEX(M, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
          ! a one to one mapping must exist, so clear the one we have already found
          KFOUND(IG1, IG2, IG3)=NK ! set k-point where it was found
       ENDDO
       ! store number of required G vectors, if inversion symmetry is used
#ifdef RPAgamma
       S%NGVECTOR_INV(NK)=M_INV    ! zero if inversion symmetry does not apply at this k-points
#endif
    ENDDO

    ! check whether all G vectors in supercell have been found
    DO M=1,S%WDES1%NGVECTOR
       IF (KFOUND(S%WDES1%IGX(M),S%WDES1%IGY(M), S%WDES1%IGZ(M))==0) THEN
          CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not map G to primitive cell: " &
             // str(S%WDES1%IGX(M)) // " " // str(S%WDES1%IGY(M)) // " " // str(S%WDES1%IGZ(M)), __FILE__, __LINE__)
       ENDIF
    ENDDO

#ifdef RPAgamma
    CALL GDES_SET_K(GDES, WDES)
#endif
    DEALLOCATE(INDEX_INTO_SUPERCELL, KFOUND)
!=====================================================================
! same for response function: index array to go from original orbitals in primitive
! cell to supercell indexing
!=====================================================================
    ! allocate a "huge" array to track the position in the supercell
    CALL SETWDES(S%WGW, S%WGW1, 1)
    NG=S%WGW1%NGVECTOR

    ALLOCATE(INDEX_INTO_SUPERCELL( & 
         MINVAL(S%WGW1%IGX(1:NG)):MAXVAL(S%WGW1%IGX(1:NG)), &
         MINVAL(S%WGW1%IGY(1:NG)):MAXVAL(S%WGW1%IGY(1:NG)), &
         MINVAL(S%WGW1%IGZ(1:NG)):MAXVAL(S%WGW1%IGZ(1:NG))))
    ALLOCATE(KFOUND( & 
         MINVAL(S%WGW1%IGX(1:NG)):MAXVAL(S%WGW1%IGX(1:NG)), &
         MINVAL(S%WGW1%IGY(1:NG)):MAXVAL(S%WGW1%IGY(1:NG)), &
         MINVAL(S%WGW1%IGZ(1:NG)):MAXVAL(S%WGW1%IGZ(1:NG))))
    ! setup the indices in that array
    INDEX_INTO_SUPERCELL=0  ! clear all entries
    KFOUND=0
    DO M=1,S%WGW1%NGVECTOR
       ! store the index M at position G_x, G_y, G_z
       INDEX_INTO_SUPERCELL(S%WGW1%IGX(M),S%WGW1%IGY(M), S%WGW1%IGZ(M))=M
    ENDDO

    ALLOCATE(S%INDEX_RES(WGW%NGDIM, WGW%NKPTS),S%NGVECTOR_RES(WGW%NKPTS))
    ALLOCATE(S%INDEX_RES_INV(WGW%NGDIM, WGW%NKPTS))
    S%INDEX_RES=0
    S%INDEX_RES_INV=0

    DO NK=1,WGW%NKPTS
       CALL SETWDES(WGW, WDES1, NK)
       S%NGVECTOR_RES(NK)=WDES1%NGVECTOR
       DO M=1,WDES1%NGVECTOR
          G1=WDES1%IGX(M)+WDES1%VKPT(1)
          G2=WDES1%IGY(M)+WDES1%VKPT(2)
          G3=WDES1%IGZ(M)+WDES1%VKPT(3)

          ! cartesian coordinates
          F1=G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)
          F2=G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)
          F3=G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)
          ! to supercell reciprocal coordinates
          G1=F1*S%LATT_CUR%A(1,1)+F2*S%LATT_CUR%A(2,1)+F3*S%LATT_CUR%A(3,1)
          G2=F1*S%LATT_CUR%A(1,2)+F2*S%LATT_CUR%A(2,2)+F3*S%LATT_CUR%A(3,2)
          G3=F1*S%LATT_CUR%A(1,3)+F2*S%LATT_CUR%A(2,3)+F3*S%LATT_CUR%A(3,3)
 
          IG1=NINT(G1)
          IG2=NINT(G2)
          IG3=NINT(G3)
          IF( ABS(G1-IG1)>1E-6 .OR. ABS(G2-IG2)>1E-6 .OR. ABS(G3-IG3)>1E-6) THEN
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not map g+k to supercell: &
                &" // str(G1) // " " // str(IG1) // " " // str(G2) // " " // str(IG2) // " " // &
                str(G3) // " " // str(IG3), __FILE__, __LINE__)
          ENDIF
#ifdef RPAgamma
          ! skip wave vectors that are not considered in half grid mode
          LINV=.FALSE.
          IF (IG1<0) THEN
             LINV=.TRUE.
          ELSE IF (IG1==0) THEN
             IF (IG2<0) THEN
                LINV=.TRUE.
             ELSE IF (IG2==0) THEN
                IF (IG3<0) THEN
                   LINV=.TRUE.
                ENDIF
             ENDIF
          ENDIF
          IF (LINV) THEN
             IF ( & 
                  -IG1<LBOUND(INDEX_INTO_SUPERCELL,1).OR. -IG1>UBOUND(INDEX_INTO_SUPERCELL,1) .OR. &
                  -IG2<LBOUND(INDEX_INTO_SUPERCELL,2).OR. -IG2>UBOUND(INDEX_INTO_SUPERCELL,2) .OR. &
                  -IG3<LBOUND(INDEX_INTO_SUPERCELL,3).OR. -IG3>UBOUND(INDEX_INTO_SUPERCELL,3)) THEN
                CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: index error " // str(IG1) // &
                   " " // str(IG2) // " " // str(IG3) // " " // str(LBOUND(INDEX_INTO_SUPERCELL)) // &
                   " " // str(UBOUND(INDEX_INTO_SUPERCELL)), __FILE__, __LINE__)
             ENDIF
             
             IF (INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)<=0) THEN
                CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not find inv g+k: " // &
                   str(-IG1) // " " // str(-IG2) // " " // str(-IG3) // " " // &
                   str(INDEX_INTO_SUPERCELL(-IG1,-IG2,-IG3)), __FILE__, __LINE__)
             ENDIF
             S%INDEX_RES_INV(M, NK)=INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)             
             KFOUND(-IG1, -IG2, -IG3)=NK ! set k-point where it was found
             CYCLE
          ENDIF
          ! test whether this G vector has already been hit does not work for RPAgamma
          ! since some vectors are found twice (once for the conjugated plane wave
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 ) THEN
#else
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 .OR. KFOUND(IG1, IG2, IG3)>0 ) THEN
#endif
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not find g+k: " // &
                str(IG1) // " " // str(IG2) // " " // str(IG3) // " " // str(INDEX_INTO_SUPERCELL(IG1,&
                IG2,IG3)) // " " // str(KFOUND(IG1,IG2,IG3)), __FILE__, __LINE__)
          ENDIF
          IF ( & 
              IG1<LBOUND(INDEX_INTO_SUPERCELL,1).OR. IG1>UBOUND(INDEX_INTO_SUPERCELL,1) .OR. &
              IG2<LBOUND(INDEX_INTO_SUPERCELL,2).OR. IG2>UBOUND(INDEX_INTO_SUPERCELL,2) .OR. &
              IG3<LBOUND(INDEX_INTO_SUPERCELL,3).OR. IG3>UBOUND(INDEX_INTO_SUPERCELL,3)) THEN
              CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: index error " // str(IG1) // &
                 " " // str(IG2) // " " // str(IG3) // " " // str(LBOUND(INDEX_INTO_SUPERCELL)) // " " &
                // str(UBOUND(INDEX_INTO_SUPERCELL)), __FILE__, __LINE__)
          ENDIF
          S%INDEX_RES(M, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
          ! a one to one mapping must exist, so clear the one we have already found
          KFOUND(IG1, IG2, IG3)=NK ! set k-point where it was found
       ENDDO
    ENDDO
    ! check whether all G vectors in supercell have been found
    DO M=1,S%WGW1%NGVECTOR
       IF (KFOUND(S%WGW1%IGX(M),S%WGW1%IGY(M), S%WGW1%IGZ(M))==0) THEN
          CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not map G to primitive cell: " &
             // str(S%WGW1%IGX(M)) // " " // str(S%WGW1%IGY(M)) // " " // str(S%WGW1%IGZ(M)), __FILE__, __LINE__)
       ENDIF
    ENDDO

    DEALLOCATE(INDEX_INTO_SUPERCELL, KFOUND)
!=====================================================================
! finally set the weight for the projectors when they are
! collected from primitive to supercell
!=====================================================================
    ALLOCATE(S%KWEIGHT(WDES%NKPTS))
    DO NK=1,WDES%NKPTS
       S%KWEIGHT(NK)=WDES%WTKPT(NK)*S%NREP
    ENDDO
!=====================================================================
! for debugging it is conventient to get set NONL_S structure
! so that the CPROJ array can be calculated directly in the supercell
!=====================================================================
    NULLIFY(S%NONL_S)
#ifdef debug
    ALLOCATE(S%NONL_S)
    CALL NONL_ALLOC(S%NONL_S,S%T_INFO,P,S%WDES, .FALSE.)
    ! set phase once and forever
    CALL SPHER(S%GRID, S%NONL_S, P, S%WDES, S%LATT_CUR,  1)
    CALL PHASE(S%WDES,S%NONL_S,1)
#endif

  END SUBROUTINE INIT_SUPERCELL_FFT

!*********************************************************************
!
!> deallocates all data structures
!> related to supercell
!
!*********************************************************************

  SUBROUTINE DEALLOCATE_SUPERCELL( S )
    TYPE (supercell), POINTER :: S !< super cell handle that is being deallocated 

    IF ( ASSOCIATED(S) ) THEN
       DEALLOCATE(S%T_INFO%NITYP, S%T_INFO%TYPE, S%T_INFO%VCA)
       DEALLOCATE(S%T_INFO%POSION, S%T_INFO%LSFOR, S%T_INFO%ITYP)
       DEALLOCATE(S%POS_IN_PRIM, S%POSION)
#ifdef VASP_HDF5
       IF(ASSOCIATED(S%T_INFO%SHA256)) DEALLOCATE(S%T_INFO%SHA256)
#endif

       DEALLOCATE(S%WDES%VKPT, S%WDES%WTKPT)
       DEALLOCATE(S%INDEX, S%INDEX_INV, S%NGVECTOR)

       IF (ASSOCIATED(S%INDEX_RES)) DEALLOCATE(S%INDEX_RES)
       IF (ASSOCIATED(S%INDEX_RES_INV)) DEALLOCATE(S%INDEX_RES_INV)
       IF (ASSOCIATED(S%NGVECTOR_RES)) DEALLOCATE(S%NGVECTOR_RES)
       DEALLOCATE(S%KWEIGHT)
       ! S%FAST_AUG reuses some elements of FAST_AUG
       ! deallocation is therefore not save 
       ! S%AUG_DES%LMMAX is not properly deallocated either
       ! so here we have some memory leakage
       CALL DEALLOC_FASTAUG_SUPER( S%FAST_AUG, S%TRANS_MATRIX )
       CALL DEALLOCWDES(S%WDES, .TRUE.)
       CALL DEALLOC_GRD(S%GRID)
       DEALLOCATE(S)
       NULLIFY(S)
    ENDIF


  END SUBROUTINE DEALLOCATE_SUPERCELL

!***************** SUBROUTINE DEALLOC_FASTAUG_SUPER ********************
!
!> Similar to #DEALLOCATE_FASTAUG, but compatible with FASTAUG_SUPER
!>
!> Deallocation of FAST_AUG\%BETA\%LPS is avoided, since it can point 
!> to some other array. 
!>
!> @note Should *not* be called to deallocate a structure that was
!        set up by #SETUP_FASTAUG
!
!***********************************************************************

  SUBROUTINE DEALLOC_FASTAUG_SUPER( FAST_AUG, TRANS_MATRIX )
    USE nonlr
    IMPLICIT NONE

    TYPE (nonlr_struct) :: FAST_AUG
    REAL(q), POINTER, CONTIGUOUS :: TRANS_MATRIX(:,:,:,:) 

    IF (.NOT. ASSOCIATED(TRANS_MATRIX)) RETURN

    DEALLOCATE(TRANS_MATRIX)
    NULLIFY( TRANS_MATRIX ) 
    CALL  NONLR_DEALLOC(FAST_AUG)

  END SUBROUTINE DEALLOC_FASTAUG_SUPER

!***********************************************************************
!
!> allocate inital entries in GDES that are required if inversion
!> symmetry is used
!
!***********************************************************************
    
  SUBROUTINE GDES_ALLOC_K(GDES, WDES)
    USE  wave
    USE  poscar
    USE  pseudo

    TYPE (greensfdes)  GDES
    TYPE (wavedes)     WDES          !< descriptor for orbitals
    ! nasty part, we need additional data in GDES
    ALLOCATE(GDES%LUSEINV(WDES%NKPTS))
    GDES%LUSEINV=.FALSE.
    ALLOCATE(GDES%NGVECTOR_INV(WDES%NKPTS), GDES%NRPLWV_ROW_DATA_POINTS_NK(WDES%NKPTS), & 
         GDES%NRPLWV_COL_DATA_POINTS_NK(WDES%NKPTS), GDES%NRPLWV_COL_MAX_DATA_POINTS_NK(WDES%NKPTS))
    GDES%NGVECTOR_INV=0
    GDES%NRPLWV_ROW_DATA_POINTS_NK=0
    GDES%NRPLWV_COL_DATA_POINTS_NK=0
    GDES%NRPLWV_COL_MAX_DATA_POINTS_NK=0

    ALLOCATE(GDES%MAP_TO_FULL(WDES%NGDIM, 2, WDES%NKPTS))

  END SUBROUTINE GDES_ALLOC_K


!***********************************************************************
!
!> set entries in GDES that are required if inversion symmetry is used
!>
!> input: entries in GDES:
!> output:
!> @param[in] GDES\%NGVECTOR_INV number of plane waves after applying inversion
!>                    symmetry (half grid mode)
!> @param[in] GDES\%MAP_TO_FULL(:,1,NK) index the PW G+k has in the standard mode 
!>                     without inversion symmetry
!> @param[out] GDES%\NRPLWV_ROW_DATA_POINTS_NK
!> @param[out] GDES%\NRPLWV_COL_DATA_POINTS_NK
!> @param[out] GDES%\NRPLWV_COL_MAX_DATA_POINTS_NK
!> @param[out] GDES%\MAP_TO_FULL(:,2,NK) index of the PW -G-k has has in the standard mode 
!>                     without inversion symmetry
!***********************************************************************

  SUBROUTINE GDES_SET_K(GDES, WDES)
    USE  wave
    USE  poscar
    USE  pseudo

    TYPE (greensfdes)  GDES
    TYPE (wavedes)     WDES           ! descriptor for orbitals
    INTEGER, ALLOCATABLE :: IGRIDIND(:,:,:)
    INTEGER :: NK, NGX, NGY, NGZ, NG1, NG2, NG3, NG1I, NG2I, NG3I, NI, NI_ORIG
    
#ifdef gammareal
    GDES%NRPLWV_ROW_DATA_POINTS_NK=GDES%NGVECTOR_INV*2
#else
    GDES%NRPLWV_ROW_DATA_POINTS_NK=GDES%NGVECTOR_INV
#endif

#ifdef MPI
    GDES%NRPLWV_COL_DATA_POINTS_NK=(GDES%NRPLWV_ROW_DATA_POINTS_NK+GDES%COMM%NCPU-1)/GDES%COMM%NCPU
    GDES%NRPLWV_ROW_DATA_POINTS_NK=GDES%NRPLWV_COL_DATA_POINTS_NK*GDES%COMM%NCPU
    GDES%NRPLWV_COL_MAX_DATA_POINTS_NK=GDES%NRPLWV_COL_DATA_POINTS_NK
    CALLMPI( M_max_i(GDES%COMM, GDES%NRPLWV_COL_MAX_DATA_POINTS_NK ,SIZE(GDES%NRPLWV_COL_MAX_DATA_POINTS_NK )))

#else
    GDES%NRPLWV_COL_DATA_POINTS_NK=GDES%NRPLWV_ROW_DATA_POINTS_NK
    GDES%NRPLWV_COL_MAX_DATA_POINTS_NK=GDES%NRPLWV_ROW_DATA_POINTS_NK
#endif
    
!   now we need to set up  GDES%MAP_TO_FULL(:,2,NK)
! allocate arrays and set some values to zero
    NGX=WDES%GRID%NGX
    NGY=WDES%GRID%NGY
    NGZ=WDES%GRID%NGZ

    ALLOCATE(IGRIDIND(NGX,NGY,NGZ))

! now index the wavefunctions at the new k-point
! take G-index from old k-point, rotate G-vec, enter new index
    DO NK=1,WDES%NKPTS
       IF (GDES%LUSEINV(NK)) THEN
          ! first we set up an 3d-array, that stores the index into
          ! the compact plane wave array at the new k-point
          IGRIDIND=0
          DO NI=1,GDES%NGVECTOR_INV(NK)
             ! original index  GDES%MAP_TO_FULL(NI,1,NK)
             NI_ORIG=GDES%MAP_TO_FULL(NI,1,NK)
             NG1=WDES%IGX(NI_ORIG,NK)
             NG2=WDES%IGY(NI_ORIG,NK)
             NG3=WDES%IGZ(NI_ORIG,NK)
             IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))=NI
          ENDDO
          GDES%MAP_TO_FULL(:,2,NK)=0
          ! now loop over all G vectors and search inverted ones
          ! and store their position 
          DO NI=1,WDES%NGVECTOR(NK)
             NG1I=WDES%IGX(NI,NK)
             NG2I=WDES%IGY(NI,NK)
             NG3I=WDES%IGZ(NI,NK)
             ! apply time inversion symmetry
             ! (the actual G vector is G+k -> -G-k, and then k needs to be subtracted
             NG1=-NG1I-NINT(WDES%VKPT(1,NK)*2)
             NG2=-NG2I-NINT(WDES%VKPT(2,NK)*2)
             NG3=-NG3I-NINT(WDES%VKPT(3,NK)*2)

             IF (IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))>0) THEN
                IF (GDES%MAP_TO_FULL(IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ)),2,NK)/=0) THEN
                   CALL vtutor%bug("internal error in GDES_SET_K: found a wave vector twice " // &
                      str(NG1) // " " // str(NG2) // " " // str(NG3), __FILE__, __LINE__)
                ENDIF
                GDES%MAP_TO_FULL(IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ)),2,NK)=NI
                IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))=0
             ENDIF
          ENDDO
          ! we must have cleared IGRIDIND completely, otherwise some problem occured
          IF (SUM(IGRIDIND)/=0) THEN
             CALL vtutor%bug("internal error in GDES_SET_K: the inverted grid does not cover all " &
                // "original wave vectors G " // str(SUM(IGRIDIND)), __FILE__, __LINE__)
          ENDIF
          ! values in GDES%MAP_TO_FULL(:,:,NK) must be all distinct, except for the Gamma-point
          ! which might be found in both maps
       ENDIF
    ENDDO
    DEALLOCATE(IGRIDIND)

    CONTAINS
!
! small helper function, produces positive indices
!

    FUNCTION MP(IND,MAXI)
      IMPLICIT NONE
      INTEGER MP,MAXI,IND
      IF (IND<=0) THEN
         MP=MAXI+IND
      ELSE
         MP=IND
      ENDIF
    END FUNCTION MP
    

  END SUBROUTINE GDES_SET_K

END MODULE msupercell
