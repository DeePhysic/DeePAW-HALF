#include "symbol.inc"
!************************************************************************
! RCS:  $Id: mkpoints.F,v 1.3 2003/06/27 13:22:20 kresse Exp kresse $
!
!> This module contains the functionality to read the KPOINTS file
!> and generate the k-points in the irreducuble Brillouin zone
!
!***********************************************************************
MODULE mkpoints
  USE prec
  USE mkpoints_struct_def
  IMPLICIT NONE

  PRIVATE :: determinant,swap_row,swap_column, &
             get_minmax_indices,matrix_inverse

  INTERFACE determinant
     module procedure determinant_real, determinant_integer
  END INTERFACE

  CONTAINS

!************************************************************************
!
!> Stores the KPOINT set presently used in the static structure KPOINTS_
!> (mkpoints_struct_def::KPOINTS_)
!
!************************************************************************

    SUBROUTINE SETUP_KPOINTS_STATIC(KPOINTS)
      USE prec
      TYPE (kpoints_struct) KPOINTS   !< k-points structure
      
      NULLIFY(KPOINTS_)
      ALLOCATE(KPOINTS_)
      KPOINTS_=KPOINTS

    END SUBROUTINE SETUP_KPOINTS_STATIC

!************************************************************************
!
!> Deallocates the globally defined structure KPOINTS_
!> (mkpoints_struct_def::KPOINTS_)
!
!************************************************************************

    SUBROUTINE DEALLOC_KPOINTS_STATIC
      USE prec
    
      IF (ASSOCIATED(KPOINTS_)) THEN
         DEALLOCATE(KPOINTS_%VKPT,KPOINTS_%WTKPT,KPOINTS_%IDTET)
         DEALLOCATE(KPOINTS_)
         NULLIFY(KPOINTS_)
      ENDIF

    END SUBROUTINE DEALLOC_KPOINTS_STATIC

!************************************************************************
!
!> Same as #DEALLOC_KPOINTS_STATIC but for a kpoints_struct variable
!> KPOINTS that has to be specified as an argument
!
!************************************************************************

    SUBROUTINE DEALLOC_KPOINTS(KPOINTS)
      USE prec
      TYPE (kpoints_struct), POINTER :: KPOINTS

      IF (ASSOCIATED(KPOINTS)) THEN
         DEALLOCATE(KPOINTS%VKPT,KPOINTS%WTKPT,KPOINTS%IDTET)
         DEALLOCATE(KPOINTS)
         NULLIFY(KPOINTS)
      ENDIF

    END SUBROUTINE DEALLOC_KPOINTS

!************************************************************************
!
!> Set up a KPOINTS (mkpoints_struct_def::kpoints_struct) structure
!>
!> Reads the KPOINTS file (UNIT=14):
!> number of k-points and k-points in reciprocal lattice
!> then calls the symmetry package to generate the k-points
!>
!> The call to SETUP_KPOINTS also initialises the static KPOINTS_ structure
!> (mkpoints_struct_def::KPOINTS_) which stores the KPOINT set presently used
!
!************************************************************************

    SUBROUTINE SETUP_KPOINTS(KPOINTS, LATT_CUR, LINVERSION, LNOSYM, IU6, IU0)
      USE prec
      USE lattice
      IMPLICIT NONE

      TYPE (kpoints_struct) KPOINTS   !< k-points structure
      TYPE (latt) LATT_CUR            !< lattice parameters
      LOGICAL LINVERSION              !< use time inversion symmetry
                                      !< if SO-coupling is used time inversion symmetry misses
      LOGICAL LNOSYM                  !< LNOSYM no symmetry in k-point generation
      INTEGER IU0,IU6                 !< units for output

      CALL RD_KPOINTS(KPOINTS, LATT_CUR, LINVERSION, LNOSYM, IU6, IU0)

      CALL SETUP_KPOINTS_STATIC(KPOINTS)

    END SUBROUTINE SETUP_KPOINTS


    SUBROUTINE RD_KPOINTS(KPOINTS, LATT_CUR, LINVERSION, LNOSYM, IU6, IU0, MULTIPL, FILE_NAME)
      USE prec
      USE sym_prec
      USE lattice
      USE constant
      USE main_mpi
      USE string, ONLY: str, lowercase
      USE tutor, ONLY: vtutor, isError, isAlert, newArgument, GammaK, LINTET, KPOINTS_IO
      USE mathtools, ONLY: matrix_inverse_3x3i
#ifdef VASP_HDF5
      USE vhdf5_base
#endif

      IMPLICIT NONE

      TYPE (kpoints_struct) KPOINTS   !< k-points structure
      TYPE (latt) LATT_CUR            !< lattice parameters
      LOGICAL LINVERSION              !< use time inversion symmetry
                                      !< time inversion symmetry is lacking if SO is used
      LOGICAL LNOSYM                  !< LNOSYM no symmetry in k-point generation
      INTEGER IU0,IU6                 !< units for output
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL :: MULTIPL !< multiplier
      CHARACTER(*), OPTIONAL, INTENT(IN)  :: FILE_NAME       !< read from file FILE_NAME instead of KPOINTS

! local
      CHARACTER (1) CSEL,CLINE
      REAL(q) :: BK(3,3),SHIFT(3),SUPL_SHIFT(3),BK_REC(3,3),AK(3,3),RS(3,3)
      INTEGER :: IHNF(3,3)
! required for reallocation
      REAL(q), POINTER :: VKPT(:,:),WTKPT(:),VKPT2(:,:)
      INTEGER, POINTER :: IDTET(:,:)=>NULL()
      INTEGER :: IERR,INDEX,NINTER,N
      LOGICAL :: WRITE_ALERT

      INTEGER KTH,NKPX,NKPY,NKPZ,NKP,ITET,I,J,NT,NK
      REAL(q) RKLEN,WSUM
! warnings from tutor
      INTEGER ITUT(1)
      LOGICAL DUMPXML
      CHARACTER (32) CSEL_LONG
! kpoints labels
      INTEGER, PARAMETER      :: KLABEL_SIZE = 20 !< max size of each label
      INTEGER, PARAMETER      :: KLABELS_INCREMENT = 10 !< steps in which to increase the size of the klabels array
      INTEGER                 :: NKL, NKLABELS
      INTEGER                 :: COMMENTPOS, COMMENTPOS_TMP
      CHARACTER (200) :: KLINE
      CHARACTER (KLABEL_SIZE) :: KLABEL, KLABEL_TMP
      CHARACTER (KLABEL_SIZE), ALLOCATABLE :: KLABELS(:), KLABELS_TMP(:)
      INTEGER, ALLOCATABLE  :: KLABELSPOS(:), KLABELSPOS_TMP(:)

#ifdef VASP_HDF5
      INTEGER(HID_T) :: IH5OUTKPTSGROUP_ID, IH5INKPTSGROUP_ID
      CHARACTER (200) :: MY_SUBGRP_KPOINTS
      CHARACTER(LEN=:), ALLOCATABLE :: SZNAMK, CSEL_HDF
#endif

      ! If no output to IU6 or IU0, also assume no output to XML
      DUMPXML = .TRUE. ; IF (IU6==-1.AND.IU0==-1) DUMPXML = .FALSE.

      IERR=0
      
      IF (PRESENT(FILE_NAME)) THEN
         OPEN(UNIT=14,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IERR)
      ELSE
         OPEN(UNIT=14,FILE=DIR_APP(1:DIR_LEN)//'KPOINTS',STATUS='OLD',IOSTAT=IERR)
         IF (IERR/=0) THEN
            OPEN(UNIT=14,FILE='KPOINTS',STATUS='OLD',IOSTAT=IERR)
         ENDIF
      ENDIF
      WRITE_ALERT = IERR/=0
#ifdef VASP_HDF5
      KPOINTS_FOUND = .FALSE. ; IF (IERR==0) KPOINTS_FOUND = .TRUE.

      KPOINTS_HDF5_FOUND = .FALSE.
      MY_SUBGRP_KPOINTS = SUBGRP_KPOINTS; IF (PRESENT(FILE_NAME)) MY_SUBGRP_KPOINTS = LOWERCASE(FILE_NAME)
      IF (.NOT.KPOINTS_FOUND .AND. HDF5_FOUND) THEN
        IERR = VH5_GROUP_OPEN(IH5ININPUTGROUP_ID, TRIM(MY_SUBGRP_KPOINTS), IH5INKPTSGROUP_ID)
        IF (IERR == 0) KPOINTS_HDF5_FOUND = .TRUE.
      ENDIF
      VH5_CHECK( VH5_GROUP_OPEN_OR_CREATE(IH5OUTINPUTGROUP_ID, TRIM(MY_SUBGRP_KPOINTS), IH5OUTKPTSGROUP_ID) )
      WRITE_ALERT = WRITE_ALERT.AND..NOT.KPOINTS_HDF5_FOUND
#endif
      IF (WRITE_ALERT .AND. .NOT.LKSPACING_WAS_READ) THEN
         CALL vtutor%alert('The requested file ' // trim(FILE_NAME) // &
            ' could not be found or opened for reading k-point information, &
            &and the KSPACING tag was also not set. &
            &Automatic k-point generation is used as a fallback, which may &
            &lead to unwanted results.')
      ENDIF

      KPOINTS%NTET=0
      KPOINTS%NKPTS=0
      KPOINTS%NKPTS_NON_ZERO=0
      KPOINTS%SZNAMK="read from INCAR"

      IF (KPOINTS%LGAMMA) THEN
        CSEL='G'
      ELSE
        CSEL='M'
      ENDIF

      IF (IU6>=0) WRITE(IU6,*)
#ifdef VASP_HDF5
      found1: IF (KPOINTS_FOUND .OR. KPOINTS_HDF5_FOUND) THEN
#else
      found1: IF (IERR==0) THEN
#endif

        ! no error, e.g. KPOINTS file exists, read in the header
        ITUT(1)=1
#ifdef VASP_HDF5
        IF (KPOINTS_FOUND) THEN
#endif
        READ(14,'(A40)',ERR=70111,END=70111) KPOINTS%SZNAMK
#ifdef VASP_HDF5
        ENDIF
        IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
          VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "system", SZNAMK) )
          IF (ALLOCATED(SZNAMK)) KPOINTS%SZNAMK = SZNAMK
        ENDIF
        VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "system", KPOINTS%SZNAMK, SKIP_PRESENT=.TRUE.))
#endif
        IF (IU6>=0) WRITE(IU6,*)'KPOINTS: ',KPOINTS%SZNAMK

        ! read the number of k-points
        ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
        IF (KPOINTS_FOUND) THEN
#endif
        READ(14,*,ERR=70111,END=70111) KPOINTS%NKPTS
#ifdef VASP_HDF5
        ENDIF
        IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
          VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "number_kpoints", KPOINTS%NKPTS) )
        ENDIF
        VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "number_kpoints", KPOINTS%NKPTS, SKIP_PRESENT=.TRUE.))
#endif

        ! read the mode selector
        ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
        IF (KPOINTS_FOUND) THEN
#endif
        READ(14,*,ERR=70111,END=70111) CSEL_LONG
        CALL STRIP( CSEL_LONG, I, 'A' )
        ! only first letter matters in the following
        CSEL = CSEL_LONG(1:1)
#ifdef VASP_HDF5
        ENDIF
        IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
          VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "mode", CSEL_HDF) )
          IF (ALLOCATED(CSEL_HDF)) CSEL = CSEL_HDF
        ENDIF
#endif

        ! If CSEL is starting with l for line, k-points are interpolated
        ! between the points read from KPOINTS
        IF (CSEL=='L'.OR.CSEL=='l') THEN
#ifdef VASP_HDF5
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'mode', 'l', SKIP_PRESENT=.TRUE.))
#endif
          CLINE='L'
          ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
          IF (KPOINTS_FOUND) THEN
#endif
          READ(14,'(A1)',ERR=70111,END=70111) CSEL
#ifdef VASP_HDF5
          ENDIF
          IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "coordinate_space", CSEL_HDF) )
            IF (ALLOCATED(CSEL_HDF)) CSEL = CSEL_HDF
          ENDIF
#endif
          KPOINTS%NKPTS=MAX( KPOINTS%NKPTS,2)
        ELSE
          CLINE=" "
        ENDIF

        IF (CSEL=='K'.OR.CSEL=='k'.OR. CSEL=='C'.OR.CSEL=='c') THEN
          CSEL='K'
!!          IF (IU6 >=0 .AND. KPOINTS%NKPTS>0)  WRITE(IU6,*)' k-points in cartesian coordinates'
#ifdef VASP_HDF5
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'coordinate_space', 'C', SKIP_PRESENT=.TRUE.))
#endif
        ELSE
!!          IF (IU6 >= 0 .AND. KPOINTS%NKPTS>0)  WRITE(IU6,*)' k-points in reciprocal lattice'
#ifdef VASP_HDF5
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'coordinate_space', 'R', SKIP_PRESENT=.TRUE.))
#endif
        ENDIF

      ENDIF found1


!=======================================================================
! read in a set of k-points
!=======================================================================
      auto: IF (KPOINTS%NKPTS>0) THEN
#ifndef oldsym
        CALL SET_SPINROT_WRAPPER(LATT_CUR%B(1,1),-1)
#endif
!=======================================================================
! ... and interpolate NKPTS between each
!=======================================================================
        kr: IF (CLINE=='L') THEN

          KPOINTS%MODE=LineMode

          IF (IU6>=0) WRITE(IU6,'(/A)')'Interpolating k-points between supplied coordinates.'

          IF (KPOINTS%LTET) THEN
            CALL vtutor%write(isAlert, LINTET)
            KPOINTS%ISMEAR = 1
            KPOINTS%LTET = .FALSE.
          ENDIF

          ALLOCATE(VKPT2(3,NKDIMD))
          ALLOCATE(KLABELS(KLABELS_INCREMENT))
          ALLOCATE(KLABELSPOS(KLABELS_INCREMENT))

          NINTER=KPOINTS%NKPTS
#ifdef VASP_HDF5
          IF (KPOINTS_FOUND .AND. .NOT.KPOINTS_HDF5_FOUND) THEN
#endif
          NKP=0  ! counter for the number of k-points already read in
          NKL=0  ! counter for the number of labels already read in
          DO
            NKP=NKP+1
            IF (NKP>NKDIMD) THEN
              CALL vtutor%bug("NKDIMD too small " // str(NKDIMD),__FILE__,__LINE__)
            ENDIF
            ITUT(1)=ITUT(1)+1
            KLINE=''
            IERR=0
            ! read until finding a non-empty line in the KPOINTS file
            DO WHILE (LEN(TRIM(KLINE))==0)
                READ(14,'(A)',IOSTAT=IERR) KLINE
                IF (IERR/=0) EXIT
            ENDDO
            KLABEL_TMP=''
            ! try to read with label
            READ(KLINE,*,IOSTAT=IERR) VKPT2(1,NKP),VKPT2(2,NKP),VKPT2(3,NKP),KLABEL_TMP
            ! strip comments from the label
            COMMENTPOS_TMP = SCAN(KLABEL_TMP,'!')
            COMMENTPOS = LEN(KLABEL_TMP); IF (COMMENTPOS_TMP>0) COMMENTPOS = COMMENTPOS_TMP
            KLABEL = KLABEL_TMP(1:COMMENTPOS-1)
            ! if the label was read then store it
            IF (LEN(TRIM(KLABEL))>0) THEN
                NKL=NKL+1
                ! check if we need to increase the size of the klabels array
                IF (NKL>SIZE(KLABELS)) THEN
                   NKLABELS=SIZE(KLABELS)
                   ALLOCATE(KLABELS_TMP(NKLABELS+KLABELS_INCREMENT))
                   KLABELS_TMP(:NKLABELS) = KLABELS
                   CALL MOVE_ALLOC(KLABELS_TMP,KLABELS)
                   ALLOCATE(KLABELSPOS_TMP(NKLABELS+KLABELS_INCREMENT))
                   KLABELSPOS_TMP(:NKLABELS) = KLABELSPOS
                   CALL MOVE_ALLOC(KLABELSPOS_TMP,KLABELSPOS)
                ENDIF
                KLABELSPOS(NKL) = NKP
                KLABELS(NKL) = KLABEL
            ENDIF
            ! if it does not work read only the coordinates
            IF (IERR/=0) READ(KLINE,*,IOSTAT=IERR) VKPT2(1,NKP),VKPT2(2,NKP),VKPT2(3,NKP)
            ! if this fails continue
            IF (IERR/=0) EXIT
          ENDDO
#ifdef VASP_HDF5
          ELSE IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "number_line_segments", NKP) )
            IF (NKP>NKDIMD) CALL vtutor%bug("NKDIMD too small " // str(NKDIMD),__FILE__,__LINE__)
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "coordinates_kpoints", VKPT2(:,:NKP)) )
            NKP = NKP+1
          ELSE
            CALL vtutor%bug("KPOINTS file not found neither KPOINTS section in hdf5 file.",__FILE__,__LINE__)
          ENDIF
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'number_line_segments', NKP-1, SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'coordinates_kpoints', VKPT2(:,1:NKP-1), SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'labels_kpoints', KLABELS(:NKL), SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'positions_labels_kpoints', KLABELSPOS(:NKL), SKIP_PRESENT=.TRUE.))
#endif
          KPOINTS%NKPTS=NKP-1
          IF (CSEL=='K') THEN
            VKPT2(:,1:KPOINTS%NKPTS)=VKPT2(:,1:KPOINTS%NKPTS)/LATT_CUR%SCALE
            CALL KARDIR(KPOINTS%NKPTS,VKPT2,LATT_CUR%A)
          ENDIF

          INDEX=0
          ! make NKPTS even
          KPOINTS%NKPTS=(KPOINTS%NKPTS/2)*2

          ALLOCATE(VKPT(3,(KPOINTS%NKPTS/2)*NINTER),WTKPT((KPOINTS%NKPTS/2)*NINTER))
          KPOINTS%NKDIM=(KPOINTS%NKPTS/2)*NINTER

          IF (IU6>=0) WRITE(IU6,'(/X,A)') 'Line segments in fractional coordinates (reciprocal lattice):'

          DO NKP=1,KPOINTS%NKPTS-1,2
            IF(IU6>=0) WRITE(IU6,'(2X,3F10.6,2X,"...",3F10.6)') VKPT2(:,NKP),VKPT2(:,NKP+1)
            SHIFT=(VKPT2(:,NKP+1)-VKPT2(:,NKP))/(NINTER-1)
            DO N=0,NINTER-1
              INDEX=INDEX+1
              IF (INDEX>KPOINTS%NKDIM) THEN
                CALL vtutor%bug("VKPT and WTKPT too small " // str((KPOINTS%NKPTS/2)*NINTER),__FILE__,__LINE__)
              ENDIF
              VKPT(:,INDEX)=VKPT2(:,NKP)+SHIFT*N
              WTKPT(INDEX)=1._q/(KPOINTS%NKPTS/2*NINTER)
            ENDDO
          ENDDO
          KPOINTS%NKPTS=(KPOINTS%NKPTS/2)*NINTER

          IF (IU6>=0) WRITE(IU6,'(/X,A,X,I5,X,A,X,I6,X,A/)') &
             'Generated',NINTER,'points on each line segment, for a total of',KPOINTS%NKPTS,'k-points:'

          ! write the kpoints to the vasprun.xml file
          IF (DUMPXML) CALL XML_KPOINTS_3(KPOINTS%NKPTS, VKPT, WTKPT, NINTER, NKL, KLABELS, KLABELSPOS)

          DEALLOCATE(KLABELS)
          DEALLOCATE(KLABELSPOS)
          DEALLOCATE(VKPT2)

        ELSE kr
!=======================================================================
!  Just read in a given set of arbitrary k-points:
!=======================================================================

          KPOINTS%MODE=ExplicitList

          IF (IU6>=0 ) WRITE(IU6,'(/A)') 'Read explicit list of k-points.'

#ifdef VASP_HDF5
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "mode", "e", SKIP_PRESENT=.TRUE.))
#endif

          ALLOCATE(VKPT(3,KPOINTS%NKPTS),WTKPT(KPOINTS%NKPTS))
          ALLOCATE(KLABELS(KLABELS_INCREMENT))
          ALLOCATE(KLABELSPOS(KLABELS_INCREMENT))
          KPOINTS%NKDIM=KPOINTS%NKPTS

          NKL=0 ! counter for the number of labels already read in
          DO NKP=1,KPOINTS%NKPTS
            ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
            IF (KPOINTS_FOUND) THEN
#endif
            KLINE=''
            IERR=0
            ! read until finding a non-empty line in the KPOINTS file
            DO WHILE (LEN(TRIM(KLINE))==0)
                READ(14,'(A)',IOSTAT=IERR) KLINE
                IF (IERR/=0) EXIT
            ENDDO
            KLABEL_TMP=''
            ! try to read with label
            READ(KLINE,*,IOSTAT=IERR) VKPT(1,NKP),VKPT(2,NKP),VKPT(3,NKP),WTKPT(NKP),KLABEL_TMP
            ! strip comments from the label
            COMMENTPOS_TMP = SCAN(KLABEL_TMP,'!')
            COMMENTPOS = LEN(KLABEL_TMP); IF (COMMENTPOS_TMP>0) COMMENTPOS = COMMENTPOS_TMP
            KLABEL = KLABEL_TMP(1:COMMENTPOS-1)
            ! if the label was read then store it
            IF (LEN(TRIM(KLABEL))>0) THEN
                NKL=NKL+1
                ! check if we need to increase the size of the klabels array
                IF (NKL>SIZE(KLABELS)) THEN
                   NKLABELS=SIZE(KLABELS)
                   ALLOCATE(KLABELS_TMP(NKLABELS+KLABELS_INCREMENT))
                   KLABELS_TMP(:NKLABELS) = KLABELS
                   CALL MOVE_ALLOC(KLABELS_TMP,KLABELS)
                   ALLOCATE(KLABELSPOS_TMP(NKLABELS+KLABELS_INCREMENT))
                   KLABELSPOS_TMP(:NKLABELS) = KLABELSPOS
                   CALL MOVE_ALLOC(KLABELSPOS_TMP,KLABELSPOS)
                ENDIF
                KLABELSPOS(NKL) = NKP
                KLABELS(NKL) = KLABEL
            ENDIF
            ! if it does not work read only the coordinates
            IF (IERR/=0) READ(KLINE,*,ERR=70111,END=70111) VKPT(1,NKP),VKPT(2,NKP),VKPT(3,NKP),WTKPT(NKP)
#ifdef VASP_HDF5
            ENDIF
#endif
          ENDDO
#ifdef VASP_HDF5
          IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "coordinates_kpoints", VKPT(:,1:KPOINTS%NKPTS)) )
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "weights_kpoints", WTKPT(1:KPOINTS%NKPTS)) )
          ENDIF
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'coordinates_kpoints', VKPT(:,1:KPOINTS%NKPTS), SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'weights_kpoints', WTKPT(1:KPOINTS%NKPTS), SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'labels_kpoints', KLABELS(:NKL), SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'positions_labels_kpoints', KLABELSPOS(:NKL), SKIP_PRESENT=.TRUE.))
#endif

          WSUM=0
          DO NKP=1,KPOINTS%NKPTS
            WSUM=WSUM+WTKPT(NKP)
          ENDDO

          IF (WSUM==0) CALL vtutor%error("RD_KPOINTS: sum of weights is zero")

          IF (CSEL=='K') THEN
            VKPT(:,1:KPOINTS%NKPTS)=VKPT(:,1:KPOINTS%NKPTS)/LATT_CUR%SCALE
            CALL KARDIR(KPOINTS%NKPTS,VKPT,LATT_CUR%A)
          ENDIF

          WTKPT(1:KPOINTS%NKPTS)=WTKPT(1:KPOINTS%NKPTS)/WSUM

          IF (KPOINTS%LTET) THEN
            ! Read in tetrahedra if you want to use tetrahedron method:
            ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
            IF (KPOINTS_FOUND) THEN
#endif
            READ(14,'(A)',ERR=70111,END=70111) CSEL
            IF (CSEL/='T' .AND. CSEL/='t') GOTO 70111
            ! read the number of tetrahedra
            ITUT(1)=ITUT(1)+1
            READ(14,*,ERR=70111,END=70111) KPOINTS%NTET,KPOINTS%VOLWGT

            ! allocate and read the tetrahedra ids
            ALLOCATE(IDTET(0:4,KPOINTS%NTET))
            DO ITET=1,KPOINTS%NTET
              ITUT(1)=ITUT(1)+1
              READ(14,*,ERR=70111,END=70111) (IDTET(KTH,ITET),KTH=0,4)
            ENDDO
#ifdef VASP_HDF5
            ENDIF
#endif
#ifdef VASP_HDF5
            IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
              VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "number_tetrahedras", KPOINTS%NTET) )
              VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "volume_weight_tetrahedra", KPOINTS%VOLWGT) )

              ALLOCATE(IDTET(0:4,KPOINTS%NTET))
              VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "coordinate_id_tetrahedra", IDTET(:,1:KPOINTS%NTET)) )
            ENDIF
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'number_tetrahedras', KPOINTS%NTET, SKIP_PRESENT=.TRUE.))
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'volume_weight_tetrahedra', KPOINTS%VOLWGT, SKIP_PRESENT=.TRUE.))
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'coordinate_id_tetrahedra', IDTET(:,1:KPOINTS%NTET), SKIP_PRESENT=.TRUE.))
#endif
          ELSE
            ALLOCATE(IDTET(0:4,0))
          ENDIF

          IF (IU6>=0) WRITE(IU6,'(/X,A,X,I6,X,A/)') 'Found',KPOINTS%NKPTS,'k-points:'

          ! write the kpoints and tetrahedra to the vasprun.xml file
          IF (DUMPXML) CALL XML_KPOINTS_1(KPOINTS%NKPTS, VKPT, WTKPT, KPOINTS%NTET, IDTET, KPOINTS%VOLWGT, NKL, KLABELS, KLABELSPOS)

          DEALLOCATE(KLABELS)
          DEALLOCATE(KLABELSPOS)

        ENDIF kr

        ! write the list of k-points to OUTCAR
        IF (IU6>=0) THEN
          WRITE(IU6,'(X,A)' ) 'Following reciprocal coordinates:'
          WRITE(IU6,'(X,A)' ) '           Coordinates               Weight'
          DO I=1,KPOINTS%NKPTS
             WRITE(IU6,'(3F10.6,4X,F10.6)') VKPT(1,I),VKPT(2,I),VKPT(3,I),WTKPT(I)
          ENDDO
          WRITE(IU6,'(/X,A)') 'Following cartesian coordinates:'
          WRITE(IU6,'(X,A)' ) '           Coordinates               Weight'
          DO I=1,KPOINTS%NKPTS
             WRITE(IU6,'(3F10.6,4X,F10.6)') MATMUL(LATT_CUR%B,VKPT(:,I))*LATT_CUR%SCALE,WTKPT(I)
          ENDDO
          WRITE(IU6,*)
        ENDIF

        ! These settings signify the set of k-points was specified explicitly
        KPOINTS%NKPX=-1
        KPOINTS%NKPY=-1
        KPOINTS%NKPZ=-1

        KPOINTS%B=0
        KPOINTS%SHIFT=0

      ELSE auto
!=======================================================================
! Automatic generation of a mesh if KPOINTS%NKPTS<=0:
!=======================================================================
        IF (IU6>=0 ) WRITE(IU6,'(/A)') 'Automatic generation of k-mesh.'

        ALLOCATE(VKPT(3,NKDIMD),WTKPT(NKDIMD),IDTET(0:4,NTETD))
        WTKPT = 0
        KPOINTS%NKDIM=NKDIMD

        SHIFT(1)=0._q
        SHIFT(2)=0._q
        SHIFT(3)=0._q
        SUPL_SHIFT=SHIFT

!=======================================================================
! Read in the generating vectors for the k-lattice and a shift of
! this lattice w.r.t. (0, 0, 0):
!=======================================================================
        mode: IF ((CSEL/='M').AND.(CSEL/='m').AND. &
     &            (CSEL/='G').AND.(CSEL/='g').AND. &
     &            (CSEL/='A').AND.(CSEL/='a')) THEN

          KPOINTS%MODE=GeneratingLattice

          ! k-lattice basis vectors in cartesian or reciprocal coordinates?
          IF (CSEL=='K'.OR.CSEL=='k'.OR. CSEL=='C'.OR.CSEL=='c') THEN
            CSEL='K'
!!            IF (IU6>=0 ) WRITE(IU6,*)' k-lattice basis in cartesian coordinates'
!!          ELSE
!!            IF (IU6>=0 )WRITE(IU6,*)' k-lattice basis in reciprocal lattice'
          ENDIF

          ! Read in the basis vectors for the k-lattice (unscaled!):
#ifdef VASP_HDF5
          IF (KPOINTS_FOUND) THEN
#endif
          ITUT(1)=ITUT(1)+1
          READ(14,*,ERR=70111,END=70111) BK(1,1),BK(2,1),BK(3,1)
          ITUT(1)=ITUT(1)+1
          READ(14,*,ERR=70111,END=70111) BK(1,2),BK(2,2),BK(3,2)
          ITUT(1)=ITUT(1)+1
          READ(14,*,ERR=70111,END=70111) BK(1,3),BK(2,3),BK(3,3)
#ifdef VASP_HDF5
          ENDIF
          IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "basis_vectors", BK) )
          ENDIF
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "mode", "b", SKIP_PRESENT=.TRUE.))
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "basis_vectors", BK, SKIP_PRESENT=.TRUE.))
#endif

          ! Correct scaling with LATT_CUR%SCALE ('lattice constant'):
          IF (CSEL=='K') BK=BK/LATT_CUR%SCALE
          ! Routine IBZKPT needs cartesian coordinates:
          IF (CSEL/='K') CALL DIRKAR(3,BK,LATT_CUR%B)
          KPOINTS%B=BK

          ! Read in the shift of the k-mesh (w.r.t. Gamma): these values must be given in
          ! k-lattice basis coordinates (usually 0 or 1/2 ...):
          ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
          IF (KPOINTS_FOUND) THEN
#endif
          READ(14,*,ERR=70112,END=70112) SHIFT(1),SHIFT(2),SHIFT(3)
#ifdef VASP_HDF5
          ENDIF
          IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
            IERR = VH5_READ(IH5INKPTSGROUP_ID, "shift", SHIFT)
            IF (IERR==3) GOTO 70112
          ENDIF
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'shift', SHIFT, SKIP_PRESENT=.TRUE.))
#endif
          SUPL_SHIFT=SHIFT
70112     CONTINUE

          ! Use HNFORM to determine how many points of the generating k-lattice (BK)
          ! fall into a unit cell of the reciprocal space lattice (B).
          CALL HNFORM(LATT_CUR%B,BK,IHNF,IERR)
          ! Sanity check: the generating k-lattice (BK) and the reciprocal lattice (B) must be
          ! commensurate in the sense that B = RS * BK, where the elements of RS are integer
          CALL RECIPS(1._q,BK(1,1),BK(1,2),BK(1,3),AK(1,1),AK(1,2),AK(1,3))
          RS=MATMUL(TRANSPOSE(AK),LATT_CUR%B)
          ! RS should be an integer matrix
          !IF (ANY(ABS(RS-NINT(RS))>TINY)) &
          IF (IERR>0) &
             CALL vtutor%error("RD_KPOINTS: the k-point generating vectors and reciprocal lattice are incommensurate.")
          NKPX=IHNF(1,1)
          NKPY=IHNF(2,2)
          NKPZ=IHNF(3,3)

          ! When KBLOWUP was not specifically set in the INCAR file we will skip the blow-up phase
          ! per default for k-point meshes that are specified by means of a generating lattice
          IF (.NOT.LKBLOWUP_WAS_READ.AND.KPOINTS%LKBLOWUP) THEN
             CALL vtutor%warning('RD_KPOINTS: internal override: KBLOWUP is changed to .FALSE.\n ... skipping the blow-up phase ...')
             KPOINTS%LKBLOWUP=.FALSE.
          ENDIF

          IF (IU6>=0) WRITE(IU6,'(X,A)') 'Grid dimensions derived from generating k-lattice:'


!=======================================================================
! Gamma-centered MP grids: the number of grid points along each
! of the reciprocal space lattice vectors is determined by a single
! parameter: an inverse reciproval length, RKLEN
!=======================================================================
        ELSE IF ((CSEL=='A').OR.(CSEL=='a')) THEN mode

          KPOINTS%MODE=RkLength

          ! read in RKLEN
#ifdef VASP_HDF5
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'mode', 'a', SKIP_PRESENT=.TRUE.))
#endif
          ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
          IF (KPOINTS_FOUND) THEN
#endif
          READ(14,*) RKLEN
#ifdef VASP_HDF5
          ENDIF
          IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
            VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "rk_length", RKLEN) )
          ENDIF
          VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'rk_length', RKLEN, SKIP_PRESENT=.TRUE.))
#endif
          NKPX =MAX(1._q,RKLEN*LATT_CUR%BNORM(1)+0.5_q)
          NKPY =MAX(1._q,RKLEN*LATT_CUR%BNORM(2)+0.5_q)
          NKPZ =MAX(1._q,RKLEN*LATT_CUR%BNORM(3)+0.5_q)
          ! use multiplier if supplied
          IF(PRESENT(MULTIPL)) THEN
            NKPX=NKPX*MULTIPL(1)
            NKPY=NKPY*MULTIPL(2)
            NKPZ=NKPZ*MULTIPL(3)
          ENDIF
!!          IF (IU0>=0) WRITE(IU0,99502) NKPX,NKPY,NKPZ
!!          IF (IU6>=0) WRITE(IU6,99502) NKPX,NKPY,NKPZ
!!99502     FORMAT( ' generate k-points for:',3I5)

          DO I=1,3
            BK(I,1)=LATT_CUR%B(I,1)/FLOAT(NKPX)
            BK(I,2)=LATT_CUR%B(I,2)/FLOAT(NKPY)
            BK(I,3)=LATT_CUR%B(I,3)/FLOAT(NKPZ)
          ENDDO

          IF (IU6>=0) WRITE(IU6,'(X,A)') 'Grid dimensions derived from RKLEN:'

!=======================================================================
! Monkhorst-Pack grids: the number of grid points along each of the
! reciprocal space lattice vectors is read from the KPOINTS file
!=======================================================================
        ELSE mode

          KPOINTS%MODE=MonkhorstPack

          ! the default is determined by the INCAR-tag KSPACING
          NKPX=MAX(1,CEILING(LATT_CUR%BNORM(1)*PI*2/KPOINTS%SPACING))
          NKPY=MAX(1,CEILING(LATT_CUR%BNORM(2)*PI*2/KPOINTS%SPACING))
          NKPZ=MAX(1,CEILING(LATT_CUR%BNORM(3)*PI*2/KPOINTS%SPACING))
!!          NKPX=KPOINTS%NKPX
!!          NKPY=KPOINTS%NKPY
!!          NKPZ=KPOINTS%NKPZ
          IF(PRESENT(MULTIPL)) THEN
            ! use multiplier if supplied
            NKPX = NKPX*MULTIPL(1)
            NKPY = NKPY*MULTIPL(2)
            NKPZ = NKPZ*MULTIPL(3)
          ENDIF

          ! if k-points file was successfully opened, continue reading it
#ifdef VASP_HDF5
          found2: IF (KPOINTS_FOUND .OR. KPOINTS_HDF5_FOUND) THEN
#else
          found2: IF (IERR==0) THEN
#endif
            ! read number of grid points along b1, b2, and b3
            ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
            IF (KPOINTS_FOUND) THEN
#endif
            READ(14,*,ERR=70111,END=70111) NKPX,NKPY,NKPZ
            ! use multiplier if supplied
            IF(PRESENT(MULTIPL)) THEN
              NKPX = NKPX*MULTIPL(1)
              NKPY = NKPY*MULTIPL(2)
              NKPZ = NKPZ*MULTIPL(3)
            ENDIF
#ifdef VASP_HDF5
            ENDIF
            IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
              VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "nkpx", NKPX) )
              VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "nkpy", NKPY) )
              VH5_CHECK( VH5_READ(IH5INKPTSGROUP_ID, "nkpz", NKPZ) )
              IF (IERR ==3) GOTO 70111
              ! use multiplier if supplied
              IF(PRESENT(MULTIPL)) THEN
                NKPX = NKPX*MULTIPL(1)
                NKPY = NKPY*MULTIPL(2)
                NKPZ = NKPZ*MULTIPL(3)
              ENDIF
            ENDIF
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'nkpx', NKPX, SKIP_PRESENT=.TRUE.))
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'nkpy', NKPY, SKIP_PRESENT=.TRUE.))
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'nkpz', NKPZ, SKIP_PRESENT=.TRUE.))
#endif
               
            ! read (a possible) shift (in units of the reciprocal lattice vectors!):
            ITUT(1)=ITUT(1)+1
#ifdef VASP_HDF5
            IF (KPOINTS_FOUND) THEN
#endif
            READ(14,*,ERR=70113,END=70113) SHIFT(1),SHIFT(2),SHIFT(3)
#ifdef VASP_HDF5
            ENDIF
            IF (.NOT. KPOINTS_FOUND .AND. KPOINTS_HDF5_FOUND) THEN
              IERR = VH5_READ(IH5INKPTSGROUP_ID, "shift", SHIFT)
              IF (IERR==3) GOTO 70113
            ENDIF
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, 'shift', SHIFT, SKIP_PRESENT=.TRUE.))
#endif
            SUPL_SHIFT=SHIFT
70113       CONTINUE

            IF (IU6>=0) WRITE(IU6,'(X,A)') 'Grid dimensions read from file:'

          ELSE found2
          ! use the values determined by KSPACING

              KPOINTS%MODE=kSpacing

              IF (IU6>=0) WRITE(IU6,'(X,A)') 'Grid dimensions derived from KSPACING:'
!!            IF (IU0>=0) WRITE(IU0,99502) NKPX,NKPY,NKPZ
!!            IF (IU6>=0) WRITE(IU6,99502) NKPX,NKPY,NKPZ
          ENDIF found2

          ! Internal rescaling and centering according to Monkhorst-Pack:
          IF ((CSEL=='M').OR.(CSEL=='m')) THEN
            SHIFT(1)=SHIFT(1)+0.5_q*MOD(NKPX+1,2)
            SHIFT(2)=SHIFT(2)+0.5_q*MOD(NKPY+1,2)
            SHIFT(3)=SHIFT(3)+0.5_q*MOD(NKPZ+1,2)
#ifdef VASP_HDF5
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "mode", "m", SKIP_PRESENT=.TRUE.))
          ELSE
            VH5_CHECK(VH5_WRITE(IH5OUTKPTSGROUP_ID, "mode", "g", SKIP_PRESENT=.TRUE.))
#endif
          ENDIF

          ! The generating vectors of the Monkhorst-Pack grid
          DO I=1,3
            BK(I,1)=LATT_CUR%B(I,1)/FLOAT(NKPX)
            BK(I,2)=LATT_CUR%B(I,2)/FLOAT(NKPY)
            BK(I,3)=LATT_CUR%B(I,3)/FLOAT(NKPZ)
          ENDDO

        ENDIF mode

        IF (IU6>=0) THEN
           WRITE(IU6,'(X,A,3I5)') 'generate k-points for:',NKPX,NKPY,NKPZ

           BK_REC=BK; CALL KARDIR(3,BK_REC,LATT_CUR%A)
           WRITE(IU6,'(/X,A/)') 'Generating k-lattice:'
           WRITE(IU6,'(2X,A,21X,A)') 'Cartesian coordinates','Fractional coordinates (reciprocal lattice)'
           WRITE(IU6,'(3(2(3X,3F13.9),/))') ((BK(I,J), I=1,3), (BK_REC(I,J), I=1,3), J=1,3)

           WRITE(IU6,'(2X,A)') 'Length of vectors'
           WRITE(IU6,'(3X,3F13.9/)') (SQRT(SUM(BK(:,I)*BK(:,I))), I=1,3)

           WRITE(IU6,'(2X,A)') 'Shift w.r.t. Gamma in fractional coordinates (k-lattice)'
           WRITE(IU6,'(3X,3F13.9/)') SHIFT
        ENDIF

!=======================================================================
! Find all irreducible points in the first Brillouin zone ... :
!=======================================================================
        CALL IBZKPT(LATT_CUR%B,BK,SHIFT,KPOINTS%NKPTS,VKPT,WTKPT, &
                    KPOINTS%LTET,KPOINTS%LKMAP,KPOINTS%NTET,IDTET,KPOINTS%FBZ2IBZ,KPOINTS%INVMAP,KPOINTS%VOLWGT,LATT_CUR%SCALE, &
                    LINVERSION,LNOSYM,KPOINTS%LKBLOWUP,LSHIFT_KPOINTS,IU6)

        ! store the smith normal form
        CALL SNFORM(LATT_CUR%B,BK,KPOINTS%ISNF,KPOINTS%ILSNF)
        KPOINTS%ILSNFINV = MATRIX_INVERSE_3X3I(KPOINTS%ILSNF)
        KPOINTS%NKFULL = KPOINTS%ISNF(1,1)*KPOINTS%ISNF(2,2)*KPOINTS%ISNF(3,3)
        ! store the generating vectors
        KPOINTS%B = BK
        CALL RECIPS(1._q,BK(1,1),BK(1,2),BK(1,3),AK(1,1),AK(1,2),AK(1,3))
        KPOINTS%A = AK
        ! and the shift
        KPOINTS%SHIFT = SHIFT
        ! and NKPX/Y/Z
        KPOINTS%NKPX=NKPX
        KPOINTS%NKPY=NKPY
        KPOINTS%NKPZ=NKPZ

! mM: 03062020: at the moment there are many places in the code where a Monkhorts-Pack
! grid is assumed when KPOINTS%NKPX/Y/Z /= -1. Until such time that those cases have been
! made conditional upon KPOINTS%MODE we will have to stick to KPOINTS%NKPX/Y/Z=-1 for
! generalized regular grids (GeneratingLattice).
        IF (KPOINTS%MODE==GeneratingLattice) THEN
           KPOINTS%NKPX=-1
           KPOINTS%NKPY=-1
           KPOINTS%NKPZ=-1
        ENDIF

        ! NKPTS_NON_ZERO is the number of k-points with non-zero weight
        ! presently this will differ from NKPTS only if LSHIFT_KPOINTS
        ! is used resulting in the addition of the difference vectors
        ! (differences between all k-points)
        KPOINTS%NKPTS_NON_ZERO=KPOINTS%NKPTS
        DO NK=1,KPOINTS%NKPTS
          IF (WTKPT(NK)==0) KPOINTS%NKPTS_NON_ZERO=KPOINTS%NKPTS_NON_ZERO-1
        ENDDO

        ! write output to the vasprun.xml file
        BK_REC=BK; CALL KARDIR(3,BK_REC,LATT_CUR%A)
        IF (DUMPXML) CALL XML_KPOINTS_2(KPOINTS%NKPTS, VKPT, WTKPT, KPOINTS%NTET, IDTET, KPOINTS%VOLWGT, &
                                        CSEL, RKLEN, NKPX, NKPY, NKPZ, SUPL_SHIFT, SHIFT, BK_REC )

      ENDIF auto

#ifdef gammareal
      IF (.NOT. PRESENT(FILE_NAME)) THEN
         ! In the gamma-only version we can not have more than one kpoint
         ! and it has to be the Gamma-point (k=0), otherwise exit in error ...
         IF (KPOINTS%NKPTS/=1 .OR. VKPT(1,1)**2 + VKPT(2,1)**2 + VKPT(3,1)**2 > 1E-10_q) THEN
            CALL vtutor%write(isError, GammaK)
         ENDIF
      ENDIF
#endif

! reallocate everthing with the minimal number of kpoints and tetrahedra
      NK=KPOINTS%NKPTS
      NT=MAX(KPOINTS%NTET,1)

      ALLOCATE(KPOINTS%VKPT(3,NK),KPOINTS%WTKPT(NK),KPOINTS%IDTET(0:4,NT))
      KPOINTS%VKPT = VKPT(1:3,1:NK)
      KPOINTS%WTKPT= WTKPT(1:NK)
      IF (KPOINTS%NTET==0) THEN
        KPOINTS%IDTET= 0
      ELSE
        KPOINTS%IDTET= IDTET(0:4,1:NT)
      ENDIF

      DEALLOCATE(VKPT,WTKPT) ; IF (ASSOCIATED(IDTET)) DEALLOCATE(IDTET)
      NULLIFY(VKPT,WTKPT,IDTET)

      KPOINTS%NKDIM=NK
      ! if NKPTS_NON_ZERO is not yet set, simply copy KPOINTS%NKPTS
      IF (KPOINTS%NKPTS_NON_ZERO==0) KPOINTS%NKPTS_NON_ZERO=KPOINTS%NKPTS

! ... done
      CLOSE(UNIT=14)
      RETURN

! Here we arrive only when something went wrong reading the KPOINTS file ...
! ... write error message and stop execution
70111 CONTINUE
      CALL vtutor%write(isError, KPOINTS_IO, newArgument(ival = itut))

    END SUBROUTINE RD_KPOINTS


!***********************************************************************
!
!> Search a particular k-point and return the equivalent k-point index
!> in the array KPOINTS_F, if not found return -1
!>
!> This is equal to the routine full_kpoints::KPOINT_IN_FULL_GRID,
!> except that this routine returns -1 and no error if point is
!> not found in the KPOINTS_F array
!
!***********************************************************************

  FUNCTION KPOINT_IN_FULL_GRID_KINTER(VKPT,KPOINTS_F)
    !USE sym_prec
    USE lattice
    INTEGER :: KPOINT_IN_FULL_GRID_KINTER
    REAL(q) VKPT(:)
    TYPE (kpoints_struct) :: KPOINTS_F
    ! local
    REAL(q), PARAMETER :: TINY=1E-8_q
    INTEGER NK
    DO NK=1,KPOINTS_F%NKPTS
       IF ( &
         (ABS(MOD(VKPT(1)-KPOINTS_F%VKPT(1,NK)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(VKPT(2)-KPOINTS_F%VKPT(2,NK)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(VKPT(3)-KPOINTS_F%VKPT(3,NK)+10.5_q,1._q)-0.5_q)<TINY)) EXIT
    ENDDO

    IF (NK>KPOINTS_F%NKPTS) THEN
       ! no kpoint found, set nk=-1
       NK=-1
    ENDIF

    KPOINT_IN_FULL_GRID_KINTER=NK
  END FUNCTION KPOINT_IN_FULL_GRID_KINTER


!*************************************************************************
!
!> Copies one kpoints_struct structure onto another
!
!*************************************************************************

    SUBROUTINE COPY_KPOINTS(K1,K2)
      IMPLICIT NONE
      TYPE(kpoints_struct) :: K1,K2
      
      K2%NKDIM=K1%NKDIM
      K2%NKPTS=K1%NKPTS
      ! pointers
      IF(ASSOCIATED(K1%VKPT))THEN
         ALLOCATE(K2%VKPT(SIZE(K1%VKPT,1),SIZE(K1%VKPT,2)))
         K2%VKPT(:,:)=K1%VKPT(:,:)
      END IF
      IF(ASSOCIATED(K1%WTKPT))THEN
         ALLOCATE(K2%WTKPT(SIZE(K1%WTKPT,1)))
         K2%WTKPT(:)=K1%WTKPT(:)
      END IF
      K2%VOLWGT=K1%VOLWGT
      IF(ASSOCIATED(K1%IDTET))THEN
         ALLOCATE(K2%IDTET(SIZE(K1%IDTET,1),SIZE(K1%IDTET,2)))
         K2%IDTET(:,:)=K1%IDTET(:,:)
      END IF
      ! allocatables
      IF(ALLOCATED(K1%INVMAP))THEN
         IF (ALLOCATED(K2%INVMAP)) DEALLOCATE(K2%INVMAP)
         ALLOCATE(K2%INVMAP(SIZE(K1%INVMAP)))
         K2%INVMAP(:)=K1%INVMAP(:)
      END IF
      IF(ALLOCATED(K1%FBZ2IBZ))THEN
         IF (ALLOCATED(K2%FBZ2IBZ)) DEALLOCATE(K2%FBZ2IBZ)
         ALLOCATE(K2%FBZ2IBZ(SIZE(K1%FBZ2IBZ,1),SIZE(K1%FBZ2IBZ,2)))
         K2%FBZ2IBZ(:,:)=K1%FBZ2IBZ(:,:)
      END IF
      K2%VOLWGT=K1%VOLWGT
      K2%NTET=K1%NTET
      K2%LTET=K1%LTET
      K2%LKMAP=K1%LKMAP
      K2%ISMEAR=K1%ISMEAR
      K2%SIGMA=K1%SIGMA
      K2%EMIN=K1%EMIN
      K2%EMAX=K1%EMAX
      K2%EFERMI=K1%EFERMI
      K2%FERMI_METHOD=K1%FERMI_METHOD
      K2%EFERMI_NEDOS=K1%EFERMI_NEDOS
      K2%BANDGAP=K1%BANDGAP
      K2%NKPX=K1%NKPX
      K2%NKPY=K1%NKPY
      K2%NKPZ=K1%NKPZ
      K2%A=K1%A
      K2%B=K1%B
      K2%SHIFT=K1%SHIFT
      K2%SZNAMK=K1%SZNAMK
      K2%MODE=K1%MODE
      K2%SPACING=K1%SPACING
      K2%ISNF=K1%ISNF
      K2%ILSNF=K1%ILSNF
      K2%ILSNFINV=K1%ILSNFINV
      K2%NKFULL=K1%NKFULL
      K2%LGAMMA=K1%LGAMMA
      K2%LKBLOWUP=K1%LKBLOWUP
    END SUBROUTINE COPY_KPOINTS

    SUBROUTINE KPOINTS_ASSIGN(K1,K2)
      IMPLICIT NONE
      TYPE(kpoints_struct) :: K1,K2

      K2=K1
    END SUBROUTINE KPOINTS_ASSIGN

!******************** SUBROUTINE IBZKPT ********************************
!
!> Reduces a given k-lattice in the full 1st Brillouin zone
!> to its irreducible part (setting also the correct weights) ...
!>
!> This routine needs:
!> @param B      The reciprocal lattice vectors [cartesian]
!> @param BK     The basis vectors of the k-lattice [cartesian]
!> @param SHIFT  The shift of the k-lattice [BK-coordinates]
!> @param LTET   A switch for usage of the tetrahedron method
!> @param SCALE  The 'scaling factor' for the Bs
!> @param IU6    The unit where to write informations
!>
!> It returns:
!> @param NKPT   The number of irreducible k-points
!> @param VKPT   The k-point coordinates [B-coordinates]
!> @param WTKPT  The weighting factors for each k-point
!>
!> and if LTET=T it additionally returns:
!> @param NTET   The number of irreducible tetrahedra
!> @param VOLWGT The volume weight factor of a tetrahedron
!> @param IDTET  The degeneracy of each tetrahedron [IDTET(0,...)]\n
!>               and the corners of the tetrahedra [IDTET(1-4,...)]
!
!***********************************************************************

    SUBROUTINE IBZKPT(B,BK,SHIFT,NKPT,VKPT,WTKPT,LTET,LKMAP,NTET,IDTET,FBZ2IBZ,INVMAP_STORE,VOLWGT,SCALE, &
                      LINVERSION,LNOSYM,LKBLOWUP,LSHIFT_KPOINTS,IU6)
      USE prec
      USE main_mpi
      USE string, ONLY : str
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      REAL(q), POINTER :: VKPT(:,:),WTKPT(:)
      INTEGER :: NKPT

      INTEGER, POINTER :: IDTET(:,:)
      INTEGER :: NTET

      INTEGER, ALLOCATABLE :: FBZ2IBZ(:,:) !< Mapping of the k-points from the FBZ to the IBZ
                                           !< the first dimension contains: 1 - index of the IBZ point, 2 - index of the rotation
      INTEGER, ALLOCATABLE :: INVMAP_STORE(:) !< map symemtry operations to their inverse

      REAL(q) :: B(3,3),BK(3,3),SHIFT(3)

      REAL(q) :: SCALE,VOLWGT
      LOGICAL :: LTET,LINVERSION,LNOSYM,LKBLOWUP,LSHIFT_KPOINTS
      LOGICAL :: LKMAP

      INTEGER :: IU6

      INTEGER :: ISYMOP,IGRPOP,NROT,NROTK,INVMAP,NPCELL
      REAL(q) :: GTRANS,AP
      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &              GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      INTEGER, PARAMETER, DIMENSION(9) :: INVERS=(/-1,0,0,0,-1,0,0,0,-1/)
      REAL(q), PARAMETER :: TINY=1.E-6_q

      ! local variables
      REAL(q) :: OMEGB,OMEGBK,WEISUM
      REAL(q) :: A(3,3),AK(3,3),AKT(3,3),V(3),VR(3),T(3),BKT(3,3),TMP(3),RS(3,3)
      REAL(q) :: PB1(3),PB2(3),PB3(3),P1(3),P2(3),P3(3),CDIMR(6),CDIMK(6)

      INTEGER :: IBRAVR,IBRAVK,IBRAVB,IBRAVBK
      INTEGER :: IHNF(3,3),NKX,NKY,NKZ
      INTEGER :: NPX, NPY, NPZ
      INTEGER :: ISNF(3,3),ILSNF(3,3)

      INTEGER :: IFULLG(3,3,48),IFULLS(3,3,48),ISWRK(3,3,48), &
     &           NRFULL,NDUM,NDUMK,NA(1),INDEX(1),IRFULL
      REAL(q) :: RFULLG(3,3,48),GWRK(3,3,48),TFULL(1,3,1),TRFULL(1,3,1),SWRK(1)

      REAL(q) :: TESTV(48,3),PWORK(50,3),TWORK(48,3),WORK(50,3),SHIFTKT(3)
      INTEGER :: NCELL,IBRAVKT,IBRAVBKT,IWORK(50)

      REAL(q) :: S(3,3,48),G(3,3,48)
      INTEGER :: KGRPOP(3,3,48),NROTKM,NTEST,NK_SHIFT,NKTMP
      LOGICAL :: LINV,LEXIST,LREGULAR,LMATCH
      LOGICAL :: LZERO,LGAMMA,LIBRAVRK
      LOGICAL :: SGREQL,OCCUP

      INTEGER :: I,I1,I2,I3,IK,INDK,ITEST,ITET,KTET,IKFBZ
      INTEGER :: IU

      INTEGER, ALLOCATABLE :: HASH1(:),HASH2(:)
      CHARACTER(LEN=:), ALLOCATABLE :: MSG
      REAL(q), ALLOCATABLE :: RTMP(:,:),VTEST(:,:)
      INTEGER :: NSIZE
      INTEGER :: KINT(3)
      INTEGER :: IERR

      CHARACTER(LEN=43) :: STRIBAV

! First check: The generating cell must have a non-zero volume
      CALL CELVOL(BK(1,1),BK(1,2),BK(1,3),OMEGBK)
      IF (ABS(OMEGBK)<(TINY*TINY*TINY)) &
         CALL vtutor%error('IBZKPT: the volume of the k-mesh generating cell is zero.')
! ... and it should be positive
      IF (OMEGBK<0._q) THEN
         IF (IU6>=0) WRITE(IU6,*) 'IBZKPT: WARNING: the k-mesh generating cell has negative volume ...'
         BK=-BK
      ENDIF

! Copy INVMAP to kpoints structure
      IF (ALLOCATED(INVMAP_STORE)) DEALLOCATE(INVMAP_STORE)
      ALLOCATE(INVMAP_STORE(NROTK))
      INVMAP_STORE(:) = INVMAP(:NROTK)

!=======================================================================
! Consistency check: we allow usage of other k-mesh lattice types than
! the lattice type of the reciprocal lattice (for example bcc-reciprocal
! lattice but simple cubic k-mesh ...). But two restrictions apply here:
! both types of lattices should at least belong to the same class of
! bravais lattices (cubic, tetragonal, orthorhombic, monoclinic, ...)
! or to some higher class which is a special case of the given class -
! that means we may mix sc,fcc,bcc but not bcc and monoclinic or so ...
! - and they should be commensurable (that means the reciprocal lattice
! vectors represented in the basis of the k-lattice should have INTEGER
! components!). Otherwise PRINT ERROR MESSAGE and STOP EXECUTION !!!!!
!=======================================================================
      CALL RECIPS(1._q,B(1,1),B(1,2),B(1,3),A(1,1),A(1,2),A(1,3))
      CALL RECIPS(1._q,BK(1,1),BK(1,2),BK(1,3),AK(1,1),AK(1,2),AK(1,3))

      RS=MATMUL(TRANSPOSE(AK),B)
      IF (ANY(ABS(RS-NINT(RS))>TINY)) &
         CALL vtutor%error('IBZKPT: the k-point generating vectors and reciprocal lattice are incommensurate.')

! Determine the type of the reciprocal space lattice (B)
      PB1=B(:,1); PB2=B(:,2); PB3=B(:,3)
      CALL LATTYP(PB1,PB2,PB3,IBRAVR,CDIMR,100)
      IBRAVB=IBRAVR
      ! All cubic types (matched onto simple cubic):
      IF (IBRAVR<=3) IBRAVR=1
      ! All tetragonal types (matched onto simple tetragonal):
      IF ((IBRAVR<=6).AND.(IBRAVR>=5)) IBRAVR=5
      ! Trigonal (rhomboedric) belong to the hexagonal system ... :
      IF (IBRAVR==7) IBRAVR=4
      ! All orthorhombic types (matched onto simple orthorhombic):
      IF ((IBRAVR<=11).AND.(IBRAVR>=8)) IBRAVR=8
      ! All monoclinic types (matched onto simple monoclinic):
      IF ((IBRAVR<=13).AND.(IBRAVR>=12)) IBRAVR=12

! Determine the type of the generating k-lattice (BK)
      P1=BK(:,1); P2=BK(:,2); P3=BK(:,3)
      CALL LATTYP(P1,P2,P3,IBRAVK,CDIMK,100)
      IBRAVBK=IBRAVK
      ! All cubic types (matched onto simple cubic):
      IF (IBRAVK<=3) IBRAVK=1
      ! All tetragonal types (matched onto simple tetragonal):
      IF ((IBRAVK<=6).AND.(IBRAVK>=5)) IBRAVK=5
      ! Trigonal (rhomboedric) belong to the hexagonal system ... :
      IF (IBRAVK==7) IBRAVK=4
      ! All orthorhombic types (matched onto simple orthorhombic):
      IF ((IBRAVK<=11).AND.(IBRAVK>=8)) IBRAVK=8
      ! All monoclinic types (matched onto simple monoclinic):
      IF ((IBRAVK<=13).AND.(IBRAVK>=12)) IBRAVK=12

! Now it is not yet the full truth ... : some lattice types can be very
! special cases of other lattice types (cubic is also tetragonal or
! tetragonal is also orthorhombic, 'all is triclinic' ...).
      LIBRAVRK=.TRUE.
      IF ((IBRAVR<IBRAVK).OR.(((IBRAVR==4).OR.(IBRAVR>=12)).AND.(IBRAVR/=IBRAVK))) THEN
         CALL vtutor%alert( &
            'Your reciprocal lattice and k-lattice belong to different lattice classes:\n\n' // &
            '   The reciprocal lattice is ' // TRIM(STRIBAV(IBRAVB )) // ',\n' // &
            '   whereas your k-lattice is ' // TRIM(STRIBAV(IBRAVBK)) // '.\n\n' // &
            'Results are often still useful ...' )
         LIBRAVRK=.FALSE.
      ENDIF
! Map shifts to the interval (-0.5,0.5] ... :
      SHIFT=MOD(SHIFT+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
      ! is the shift zero?
      LZERO=ALL(ABS(SHIFT)<TINY)

! The symmetry operations of the crystalline lattice should be members
! of the group of symmetry operations of the reciprocal lattice as well.
! Get the symmetry operations of the reciprocal space lattice
      TFULL=0._q; NA(1)=1
      CALL SETGRP(ISWRK,IFULLS,GWRK,NDUM,NDUMK,TFULL,IBRAVB,NRFULL, &
     &                                   TRFULL,1,1,NA,1,INDEX,SWRK,100)
      ! consistency test: check whether all symmetry operations
      ! IGRPOP (the "crystalline" lattice) are present in the above
      ! group.
      ! First transform the symmetry operations to the basis
      ! of the reciprocal lattice (B)
      CALL SGRCON(IFULLS,IFULLG,NRFULL,PB1,PB2,PB3,B(1,1),B(1,2),B(1,3))
      ! and try to match all members of IGRPOP
      LMATCH=.TRUE.
      op1: DO I1=1,NROTK
         DO I2=1,NRFULL
            IF (SGREQL(IGRPOP(1,1,I1),IFULLG(1,1,I2))) CYCLE op1
         ENDDO
         ! this point will be reached in case we find no match
         LMATCH=.FALSE.
         IF (.NOT.LMATCH) THEN
            CALL vtutor%error( &
               'IBZKPT: not all point group operations associated with the &
               &symmetry operators of the crystalline lattice are members of the &
               &point group of the reciprocal lattice.')
         ENDIF
      ENDDO op1

      LKBLOWUP=(LKBLOWUP.AND.(.NOT.LNOSYM))
! Get the symmetry operation that will be used in the "blow-up" phase
      IF (LKBLOWUP) THEN
         ! transform IFULLS to basis of the generating k-lattice
!        CALL SGRCON(IFULLS,IFULLG,NRFULL,PB1,PB2,PB3,BK(1,1),BK(1,2),BK(1,3))
         CALL SGRTRF(REAL(IFULLS,KIND=q),RFULLG,NRFULL,PB1,PB2,PB3,BK(1,1),BK(1,2),BK(1,3))
         ! Are any of the symmetry operations represented by
         ! non-integer rotation matrices? If so, skip the "blow-up" phase, after all.
         IF (ANY(ABS(RFULLG(:,:,1:NRFULL)-NINT(RFULLG(:,:,1:NRFULL)))>TINY)) THEN
            CALL vtutor%warning( &
               'IBZKPT: internal override: KBLOWUP is changed to .FALSE.\n ... skipping the blow-up phase ...')
            LKBLOWUP=.FALSE.
         ENDIF
      ENDIF
! To skip the "blow-up" phase, we define only the identity operator
      IF (.NOT.LKBLOWUP) THEN
         NRFULL=1
         RFULLG=0; RFULLG(1,1,1)=1; RFULLG(2,2,1)=1; RFULLG(3,3,1)=1
      ENDIF
#ifndef oldsym
      CALL SET_SPINROT_WRAPPER(B(1,1),-1)
#endif
!=======================================================================
! We may always use inversion (even if it is no real symmetry operation
! of the atomic lattice) - this is due to time-reversal symmetry ... :
!=======================================================================
      NROTKM=NROTK
      LINV=.FALSE.
      DO I=1,NROTK
         ! Test existence of inversion operator:
         IF (SGREQL(IGRPOP(1,1,I),INVERS)) LINV=.TRUE.
         ! ... and copy the symmetry operators to some temporary array ... :
         KGRPOP(:,:,I)=IGRPOP(:,:,I)
      ENDDO

      IF (LINVERSION) THEN
         ! If inversion is not yet included add inversion!
         IF (.NOT.LINV) THEN
            DO I=1,NROTK
            ! ... this means in full consequence multiply all operators with I ... :
               CALL SGRPRD(INVERS,KGRPOP(1,1,I),KGRPOP(1,1,I+NROTK))
            ENDDO
            ! ... so we have now a group with twice as much symmetry operations:
            NROTKM=2*NROTK
         ENDIF
      ENDIF
!=======================================================================
! Transform the symmetry operators to the representation for the
! basis generating the k-mesh (basis BK):
!=======================================================================
      S(:,:,1:NROTKM)=REAL(KGRPOP(:,:,1:NROTKM),KIND=q)
      CALL SGRTRF(S,G,NROTKM,B(1,1),B(1,2),B(1,3),BK(1,1),BK(1,2),BK(1,3))

!=======================================================================
! How many points of the k-lattice (BK) fall into a unit cell
! of the reciprocal space lattice (B).
!=======================================================================
      CALL HNFORM(B,BK,IHNF,IERR)
      IF (IERR>0) &
         CALL vtutor%error("IBZKPT: the k-point generating vectors and reciprocal lattice are incommensurate.")
      NKX=IHNF(1,1); NKY=IHNF(2,2); NKZ=IHNF(3,3)
      ! and compute the Smith Normal form matrices
      CALL SNFORM(B,BK,ISNF,ILSNF)

!=======================================================================
! "Blow-up" phase: apply the symmetry operation of the reciprocal space
! lattice to all points of the generating k-lattice that lie within the
! first unit cell of the reciprocal space lattice.
! Since the reciprocal space lattice and the k-lattice are not necessarily
! of the same type, and since the "shift" applied to the k-points might
! not be symmetry conserving, this phase may generate additional k-points
! that are not integer multiples of the k-lattice basis vectors (BK).
! In that case LREGULAR=.FALSE. at the end.
!=======================================================================

      ALLOCATE(HASH1(NKX*NKY*NKZ),HASH2(NKX*NKY*NKZ),VTEST(3,NKX*NKY*NKZ))

      LREGULAR=.TRUE.

      NTEST=0
      HASH1=0

#ifdef hoist_blowup
      blowup: DO IRFULL=1,NRFULL
#endif
      nk3: DO I3=0,NKZ-1
         nk2: DO I2=0,NKY-1
            nk1: DO I1=0,NKX-1
               ! k-point coordinates in the basis BK:
               V(1)=REAL(I1,KIND=q)+SHIFT(1)
               V(2)=REAL(I2,KIND=q)+SHIFT(2)
               V(3)=REAL(I3,KIND=q)+SHIFT(3)
               ! Take only k-points lying within the Brillouin zone (if we do not use
               ! the same lattice type for the k-lattice as for the reciprocal lattice
               ! it might happen that we have to cut some points ...):
               VR=MATMUL(BK,V)
               ! k-point coordinates in basis B:
               V=MATMUL(TRANSPOSE(A),VR)
               ! shift into 1st Brillouin zone (unit cell of basis B):
               V=MOD(V+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
               ! back to cartesian coordinates
               VR=MATMUL(B,V)
               ! k-point coordinates in basis BK:
               V=MATMUL(TRANSPOSE(AK),VR)
#ifndef hoist_blowup
               blowup: DO IRFULL=1,NRFULL
#endif
                  ! Rotate k-point and apply periodic boundary conditions to this k-point:
!                 VR=MATMUL(TRANSPOSE(REAL(IFULLG(:,:,IRFULL),KIND=q)),V)
                  VR=MATMUL(TRANSPOSE(RFULLG(:,:,IRFULL)),V)
                  ! bring it to cartesian coordinates
                  T=MATMUL(BK,VR)
                  ! express it in the basis B
                  VR=MATMUL(TRANSPOSE(A),T)
                  ! bring to the first unit cell (-0.5,0.5]
                  VR=MOD(VR+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
                  ! back to cartesian coordinates
                  T=MATMUL(B,VR)
                  ! and back to the basis BK
                  VR=MATMUL(TRANSPOSE(AK),T)
                  ! First check whether this points is an element of the
                  ! k-lattice generated by BK
                  T=VR-SHIFT
                  IF (ANY(MOD(ABS(T-NINT(T)),1._q)>TINY)) THEN
                     ! If it is not:
                     LREGULAR=.FALSE.
                     ! And we have to test whether point has already
                     ! been generated the hard way ... :
                     LEXIST=.FALSE.
                     DO ITEST=1,NTEST
                        IF (ALL(ABS(VR-VTEST(:,ITEST))<TINY)) LEXIST=.TRUE.
                     ENDDO
                  ELSE
                     ! If this point is an element of the k-lattice generated by BK
                     ! we can use a hashtable to find out whether it was reached already
                     INDK=KINDEX_SNFORM(T,ISNF,ILSNF)
                     LEXIST=.FALSE.
                     IF (HASH1(INDK)/=0) LEXIST=.TRUE.
                     HASH1(INDK)=HASH1(INDK)+1
                  ENDIF

                  IF (.NOT.LEXIST) THEN
                     NTEST=NTEST+1
                     IF (NTEST>SIZE(VTEST,2)) THEN
                        ! store VTEST in temporary array
                        ALLOCATE(RTMP(3,SIZE(VTEST,2))); RTMP=VTEST
                        ! destroy VTEST
                        DEALLOCATE(VTEST)
                        ! reallocate VTEST with larger second dimension
                        NSIZE=SIZE(RTMP,2)*(2-REAL(I3*NKX*NKY+I2*NKX+I1+1,KIND=q)/REAL(NKX*NKY*NKZ,KIND=q))+100
                        ALLOCATE(VTEST(3,NSIZE))
                        ! copy back to VTEST
                        VTEST(:,1:SIZE(RTMP,2))=RTMP
                       ! destroy temporary array
                        DEALLOCATE(RTMP)
                     ENDIF
                     IF (NTEST>SIZE(VTEST,2)) CALL vtutor%bug('NTEST>SIZE(VTEST)',__FILE__,__LINE__)
                     VTEST(:,NTEST)=VR
                     IF (ALL(MOD(ABS(T-NINT(T)),1._q)<TINY)) HASH1(INDK)=1
                  ENDIF
#ifndef hoist_blowup
               ENDDO blowup
#endif
            ENDDO nk1
         ENDDO nk2
      ENDDO nk3
#ifdef hoist_blowup
      ENDDO blowup
#endif

!=======================================================================
! Now run over all mesh points and find out the irreducible points ...
! Reduce w.r.t. the symmetry operations of the "crystalline" lattice
! (possibly inversion symmetry has been added to this set of operations)
!=======================================================================
      LGAMMA=.FALSE.

      NKPT =0
      HASH2=0

      DO ITEST=1,NTEST
         V=VTEST(:,ITEST)
         ! is gamma part of the set of k-points?
         LGAMMA=LGAMMA.OR.ALL(ABS(V)<TINY)
         ! Set test flag ...
         LEXIST=.FALSE.
         ! For all symmetry operations
         reduce: DO I=1,NROTKM
            ! Rotate k-point
            VR=MATMUL(TRANSPOSE(G(:,:,I)),V)
            ! bring it to cartesian coordinates
            T=MATMUL(BK,VR)
            ! express it in basis B
            VR=MATMUL(TRANSPOSE(A),T)
            ! bring to the first unit cell (-0.5,0.5]
            VR=MOD(VR+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
            ! back to cartesian coordinates
            T=MATMUL(B,VR)
            ! and back to the basis BK
            VR=MATMUL(TRANSPOSE(AK),T)
            ! First check whether this points is an element of the
            ! k-lattice generated by BK
            T=VR-SHIFT
            IF (ANY(MOD(ABS(T-NINT(T)),1._q)>TINY)) THEN
               ! If not, test the rotated k-point against all previously found k-points ... the hard way:
               DO IK=1,NKPT
                  IF (ALL(ABS(VR-VKPT(:,IK))<TINY)) THEN
                  ! Here we found that the rotated k-point is equivalent to some k-point
                  ! which has already been stored in our list ... ! This means that our
                  ! original k-point is also equivalent to this k-point! So set test flag
                  ! to remark existence and add up the weight for this k-point ... :
                     LEXIST=.TRUE.
                     WTKPT(IK)=WTKPT(IK)+1
                     EXIT reduce
                  ENDIF
               ENDDO
            ELSE
               ! If this point is an element of the k-lattice generated by BK
               ! we can use a hashtable to find out whether it was reached already
               INDK=KINDEX_SNFORM(T,ISNF,ILSNF)
               IF (HASH2(INDK)/=0) THEN
                  ! point was already found before: update weight and set the test flag
                  LEXIST=.TRUE.
                  WTKPT(HASH2(INDK))=WTKPT(HASH2(INDK))+1
                  EXIT reduce
               ENDIF
            ENDIF
         ENDDO reduce

         IF (.NOT.LEXIST) THEN
         ! No symmetry operation was able to reproduce any of the k-points found
         ! previously - so we must have found a new k-point! Count it, store it!
            NKPT=NKPT+1
            IF (NKPT>SIZE(WTKPT)) THEN
               ! store VKPT in temporary array
               ALLOCATE(RTMP(3,SIZE(WTKPT))); RTMP=VKPT
               ! destroy VKPT
               DEALLOCATE(VKPT); NULLIFY(VKPT)
               ! reallocate VKPT with larger second dimension
               NSIZE=SIZE(WTKPT)*(2-REAL(ITEST,KIND=q)/REAL(NTEST,KIND=q))+100
               ALLOCATE(VKPT(3,NSIZE))
               ! copy back to VKPT
               VKPT(:,1:SIZE(RTMP,2))=RTMP
               ! store WTKPT in temporary array
               RTMP(1,:)=WTKPT
               ! destroy WTKPT
               DEALLOCATE(WTKPT); NULLIFY(WTKPT)
               ! reallocate WTKPT with larger size
               ALLOCATE(WTKPT(NSIZE))
               ! copy back to WTKPT
               WTKPT(1:SIZE(RTMP,2))=RTMP(1,:)
               ! destroy temporary array
               DEALLOCATE(RTMP)
            ENDIF
            IF (NKPT>SIZE(WTKPT)) CALL vtutor%bug('NKPT>SIZE(WTKPT)',__FILE__,__LINE__)
            VKPT(:,NKPT)=V
            ! ... and initialize the 'weight counter':
            WTKPT(NKPT)=1
            ! and possibly update the hashtable
            T=V-SHIFT
            IF (ALL(MOD(ABS(T-NINT(T)),1._q)<TINY)) THEN
               INDK=KINDEX_SNFORM(T,ISNF,ILSNF)
               HASH2(INDK)=NKPT
            ENDIF

         ENDIF
      ENDDO

      DEALLOCATE(HASH1,HASH2,VTEST)

!=======================================================================
! if LSHIFT_KPOINTS is set we need to include all k-points
! that are differences of two other k-points in the full
! k-point set
!=======================================================================
      NK_SHIFT=0
      IF (LSHIFT_KPOINTS.AND.LNOSYM) THEN
         NKTMP=NKPT
         DO I1=1,NKPT
            DO I2=1,NKPT
               ! make new k-point
               VR=VKPT(:,I1)-VKPT(:,I2)
               ! bring it to cartesian coordinates
               T=MATMUL(BK,VR)
               ! express it in basis B
               VR=MATMUL(TRANSPOSE(A),T)
               ! bring to the first unit cell (-0.5,0.5]
               VR=MOD(VR+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
               ! back to cartesian coordinates
               T=MATMUL(B,VR)
               ! and back to the basis BK
               VR=MATMUL(TRANSPOSE(AK),T)
               ! check against the current list of points
               LEXIST=.FALSE.
               DO IK=1,NKTMP
                  IF (ALL(ABS(VR-VKPT(:,IK))<TINY)) THEN
                     LEXIST=.TRUE.; EXIT
                  ENDIF
               ENDDO
               ! add if it does not exist yet
               IF (.NOT.LEXIST) THEN
                  NK_SHIFT=NK_SHIFT+1
                  NKTMP=NKTMP+1
                  IF (NKTMP>SIZE(WTKPT)) THEN
                     ! store VKPT in temporary array
                     ALLOCATE(RTMP(3,SIZE(WTKPT))); RTMP=VKPT
                     ! destroy VKPT
                     DEALLOCATE(VKPT); NULLIFY(VKPT)
                     ! reallocate VKPT with larger second dimension
                     NSIZE=SIZE(WTKPT)*(2-REAL(I1,KIND=q)/REAL(NKPT,KIND=q))+100
                     ALLOCATE(VKPT(3,NSIZE))
                     ! copy back to VKPT
                     VKPT(:,1:SIZE(RTMP,2))=RTMP
                     ! store WTKPT in temporary array
                     RTMP(1,:)=WTKPT
                     ! destroy WTKPT
                     DEALLOCATE(WTKPT); NULLIFY(WTKPT)
                     ! reallocate WTKPT with larger size
                     ALLOCATE(WTKPT(NSIZE))
                     ! copy back to WTKPT
                     WTKPT(1:SIZE(RTMP,2))=RTMP(1,:)
                     ! destroy temporary array
                     DEALLOCATE(RTMP)
                  ENDIF
                  IF (NKTMP>SIZE(VKPT,2)) CALL vtutor%bug('NKPT>SIZE(VKPT,2)',__FILE__,__LINE__)
                  VKPT(:,NKTMP)=VR
                  WTKPT(NKTMP) =0
               ENDIF
            ENDDO
         ENDDO
         NKPT=NKTMP
      ENDIF

!=======================================================================
! Following some code for setting up the "connection table" for the
! k-points on the full mesh with respect to the irreducible points
! if we want to use the tetrahedron integration method ... :
!=======================================================================
      IF (LTET.OR.LKMAP) THEN
         IF (LREGULAR.OR.LNOSYM) THEN
         ! In case the grid is a regular grid generated by the
         ! k-lattice vectors (BK) and a possible non-zero shift
         ! has not violated its symmetry.
            BKT=BK; AKT=AK; SHIFTKT=SHIFT
            S(:,:,1:NROTKM)=G(:,:,1:NROTKM)
         ELSE
         ! When B has higher symmetry than BK and/or a non-zero
         ! shift has violated the symmetry of B we probably have ended
         ! up with a non-regular mesh of k-points.
         ! There is still one thing we may try, though: the star of
         ! shift vectors generated by the symmetry operations of the
         ! reciprocal lattice may constitute a supercell of a generating
         ! k-lattice (BKT) different (and smaller) than BK.
            P1=BK(:,1)
            P2=BK(:,2)
            P3=BK(:,3)
            ! build the star of shifts
            NTEST=1
            TESTV(1,:)=SHIFT
            DO IRFULL=1,NRFULL
               ! rotate shift
!              VR=MATMUL(TRANSPOSE(REAL(IFULLG(:,:,IRFULL),KIND=q)),V)
               VR=MATMUL(TRANSPOSE(RFULLG(:,:,IRFULL)),TESTV(1,:))
               ! bring to unit cell of BK: (-0.5,0.5]
               VR=MOD(VR+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
               LEXIST=.FALSE.
               DO ITEST=1,NTEST
                  IF (ALL(ABS(VR-TESTV(ITEST,:))<TINY)) LEXIST=.TRUE.
               ENDDO
               IF (.NOT.LEXIST) THEN
                  NTEST=NTEST+1; TESTV(NTEST,:)=VR
               ENDIF
            ENDDO
            ! Which primitive cell BKT builds the supercell BK?
            NA(1)=NTEST
            CALL PRICEL(IBRAVK,CDIMK,P1,P2,P3,TESTV,BKT(1,1),BKT(1,2),BKT(1,3), &
           &            PWORK,NCELL,IBRAVKT,CDIMR,1,1,NA,48,TWORK,IWORK,WORK,-1)
            IBRAVBKT=IBRAVKT
            ! All cubic types (matched onto simple cubic):
            IF (IBRAVKT<=3) IBRAVKT=1
            ! All tetragonal types (matched onto simple tetragonal):
            IF ((IBRAVKT<=6).AND.(IBRAVKT>=5)) IBRAVKT=5
            ! Trigonal (rhomboedric) belong to the hexagonal system ... :
            IF (IBRAVKT==7) IBRAVKT=4
            ! All orthorhombic types (matched onto simple orthorhombic):
            IF ((IBRAVKT<=11).AND.(IBRAVKT>=8)) IBRAVKT=8
            ! All monoclinic types (matched onto simple monoclinic):
            IF ((IBRAVKT<=13).AND.(IBRAVKT>=12)) IBRAVKT=12

            ! Our purpose is only accomplished when as many cells of BKT make up
            ! the original k-lattice cell BK as there are shifts in the star
            ! and the Bravais lattice types of both are the same:
            IF (NCELL==NTEST.AND.IBRAVK==IBRAVKT) THEN
               ! express the shift in basis BKT
               CALL RECIPS(1._q,BKT(1,1),BKT(1,2),BKT(1,3),AKT(1,1),AKT(1,2),AKT(1,3))
               VR=MATMUL(BK,SHIFT)
               SHIFTKT=MATMUL(TRANSPOSE(AKT),VR)
               ! how many points does BKT generate within a unit cell
               ! of the reciprocal lattice B
               CALL HNFORM(B,BKT,IHNF,IERR)
               IF (IERR>0) &
                  CALL vtutor%error("IBZKPT: the k-point generating vectors and reciprocal lattice are incommensurate.")
               NKX=IHNF(1,1); NKY=IHNF(2,2); NKZ=IHNF(3,3)
               ! Get the symmetry operations represented in the basis BKT
               CALL SGRTRF(G,S,NROTKM,BK(1,1),BK(1,2),BK(1,3),BKT(1,1),BKT(1,2),BKT(1,3))
               ! and recompute the Smith Normal Form matrices
               CALL SNFORM(B,BKT,ISNF,ILSNF)
            ELSE
            ! Unfortunately TETIRR can only decompose a grid of regular
            ! parallelepipeds into tetrahedra, therefore we exit in error:
               MSG = 'IBZKPT: unable to construct a generating k-lattice suitable for use with &
                  &the tetrahedron method \n... use a different smearing method'
               IF (LKBLOWUP) &
                  MSG = MSG // '\n... and/or try setting KBLOWUP=.FALSE.'
               IF (.NOT.LGAMMA) &
                  MSG = MSG // '\n... and/or try a Gamma-centered k-lattice'
               IF (.NOT.LIBRAVRK) &
                  MSG = MSG // '\n... and/or try a different generating k-lattice'
               IF (.NOT.LZERO) &
                  MSG = MSG // '\n... and/or try a different shift'
               CALL vtutor%error(MSG)
            ENDIF
         ENDIF

         ! Refill the hashtable
         ALLOCATE(HASH2(NKX*NKY*NKZ)); HASH2=0
         DO IK=1,NKPT
            V=VKPT(:,IK)
            V=MATMUL(BK,V)
            V=MATMUL(TRANSPOSE(AKT),V)
            T=V-SHIFTKT
            ! all points in the irreducible set should be members of the regular
            ! shifted grid generated by BKT and SHIFTKT
            IF (ANY(MOD(ABS(T-NINT(T)),1._q)>TINY)) THEN
               MSG = 'IBZKPT: not all points in the irreducible set are regular points &
                  &of the (shifted) generating k-lattice \n... use a different smearing method'
               IF (LKBLOWUP) &
                  MSG = MSG // '\n... and/or try setting KBLOWUP=.FALSE.'
               IF (.NOT.LGAMMA) &
                  MSG = MSG // '\n... and/or try a Gamma-centered k-lattice'
               IF (.NOT.LIBRAVRK) &
                  MSG = MSG // '\n... and/or try a different generating k-lattice'
               IF (.NOT.LZERO) &
                  MSG = MSG // '\n... and/or try a different shift'
               CALL vtutor%error(MSG)
            ENDIF
            INDK=KINDEX_SNFORM(T,ISNF,ILSNF)
            HASH2(INDK)=IK
         ENDDO

         NPX = ISNF(1,1)
         NPY = ISNF(2,2)
         NPZ = ISNF(3,3)

         IF (ALLOCATED(FBZ2IBZ)) DEALLOCATE(FBZ2IBZ)
         ALLOCATE(FBZ2IBZ(2,NPX*NPY*NPZ))
         FBZ2IBZ=-1
         ! Now go through the k-lattice and fill the connection table:
         DO I3=0,NKZ-1
            DO I2=0,NKY-1
               DO I1=0,NKX-1
                  ! k-point coordinates in the basis BKT:
                  V(1)=REAL(I1,KIND=q)+SHIFTKT(1)
                  V(2)=REAL(I2,KIND=q)+SHIFTKT(2)
                  V(3)=REAL(I3,KIND=q)+SHIFTKT(3)
                  ! Take only k-points lying within the Brillouin zone (if we do not use
                  ! the same lattice type for the k-lattice as for the reciprocal lattice
                  ! it might happen that we have to cut some points ...):
                  VR=MATMUL(BKT,V)
                  ! k-point coordinates in basis B:
                  V=MATMUL(TRANSPOSE(A),VR)
                  ! shift into 1st Brillouin zone (unit cell of basis B):
                  V=MOD(V+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
                  VR=MATMUL(B,V)
                  ! k-point coordinates in basis BKT:
                  V=MATMUL(TRANSPOSE(AKT),VR)

                  LEXIST=.FALSE.
                  match: DO I=1,NROTKM
                     ! Rotate k-point
                     VR=MATMUL(TRANSPOSE(S(:,:,I)),V)
                     T=VR-SHIFTKT
                     IF (ALL(MOD(ABS(T-NINT(T)),1._q)<TINY)) THEN
                        INDK=KINDEX_SNFORM(T,ISNF,ILSNF)
                        IF (HASH2(INDK)/=0) THEN
                           LEXIST=.TRUE.
                           IKFBZ = KINDEX_SNFORM(V-SHIFTKT,ISNF,ILSNF)
                           FBZ2IBZ(1,IKFBZ)=HASH2(INDK)
                           ! If the list of symmetry operations has been extended to include
                           ! time-reversal symmetry then NROTKM=2*NROTK. In FBZ2IBZ I encode the
                           ! the time-reversal symmetry operations using negative numbers
                           ! if index of the symmetry operation is larger than NROTK then I subtract NROTK
                           ! to recover the original index of the symmetry operation and make it negative to
                           ! know that we use time-reversal symmetry
                           IF (I<=NROTK) THEN
                              FBZ2IBZ(2,IKFBZ)=I
                           ELSE
                              FBZ2IBZ(2,IKFBZ)=-(I-NROTK)
                           ENDIF
                           EXIT match
                        ENDIF
                     ENDIF
!!                     ! bring it to cartesian coordinates
!!                     T=MATMUL(BKT,VR)
!!                     ! express it in basis B
!!                     VR=MATMUL(TRANSPOSE(A),T)
!!                     ! bring to the first unit cell (-0.5,0.5]
!!                     VR=MOD(VR+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
!!                     ! back to cartesian coordinates
!!                     T=MATMUL(B,VR)
!!                     ! and back to the basis BK
!!                     VR=MATMUL(TRANSPOSE(AK),T)
!!                     ! Test the rotated k-point against all previously found k-points ... :
!!                     DO IK=1,NKPT
!!                        IF (ALL(ABS(VR-VKPT(:,IK))<TINY)) THEN
!!                        ! Here we found that the rotated k-point is equivalent to some k-point
!!                        ! which has already been stored in our list ... ! This means that our
!!                        ! original k-point is also equivalent to this k-point! So set test flag
!!                        ! to remark existence and add up the weight for this k-point ... :
!!                           LEXIST=.TRUE.
!!                           IKPT(I1+1,I2+1,I3+1)=IK
!!                           EXIT match
!!                        ENDIF
!!                     ENDDO
                  ENDDO match
                  IF (.NOT.LEXIST) THEN
                     ! We were not able to match this k-point to any of the
                     ! NKPT points previously determined and stored in VKPT
                     MSG = 'IBZKPT: not all points of the generating k-lattice map onto &
                        &the irreducible set of k-points \n... use a different smearing method'
                     IF (LKBLOWUP) &
                        MSG = MSG // '\n... and/or try setting KBLOWUP=.FALSE.'
                     IF (.NOT.LGAMMA) &
                        MSG = MSG // '\n... and/or try a Gamma-centered k-lattice'
                     IF (.NOT.LIBRAVRK) &
                        MSG = MSG // '\n... and/or try a different generating k-lattice'
                     IF (.NOT.LZERO) &
                        MSG = MSG // '\n... and/or try a different shift'
                     CALL vtutor%error(MSG)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! We need at least 4 k-points
         IF (NKX*NKY*NKZ<4.AND..NOT.LKMAP) &
            CALL vtutor%error('IBZKPT: tetrahedron method fails for NKPT<4. NKPT = ' // str(NKX*NKY*NKZ))

         ! Set up the list of tetrahedra
         IF (LTET) CALL TETIRR(NTET,IDTET,BK,ILSNF,ISNF,NKX,NKY,NKZ,FBZ2IBZ,IU6)

         ! Get the volume weight of each tetrahedron (according to the full mesh)
         ! where the factor 6 accounts for the fact that around each k-point we
         ! have examined a microcell giving 6 tetrahedra.
         VOLWGT=1._q/REAL(6*NKX*NKY*NKZ,KIND=q)

         DEALLOCATE(HASH2)
      ENDIF !LTET.OR.LKMAP

!=======================================================================
! last task to do: transform the coordinates from basis BK to basis B
! and supply the correct scaling to the weighting factors ...
!=======================================================================
      WEISUM=0._q
      DO IK=1,NKPT
         ! express cartesian coordinates
         V=MATMUL(BK,VKPT(:,IK))
         ! k-point coordinates in basis B:
         VKPT(:,IK)=MATMUL(TRANSPOSE(A),V)
         ! Sum up weights:
         WEISUM=WEISUM+WTKPT(IK)
      ENDDO

      ! Print the result (if desired):
      IF ((IU6>=0).AND.(IU6<=99)) THEN

         WRITE(IU6,'(A)') ' '
         WRITE(IU6,'(A)') ' Subroutine IBZKPT returns following result:'
         WRITE(IU6,'(A)') ' ==========================================='
         WRITE(IU6,'(A)') ' '

         WRITE(IU6,'(A,I6,A)') ' Found ',NKPT,' irreducible k-points:'
         IF (NK_SHIFT>0) THEN
            WRITE(IU6,'(A,I6,A)') ' Added ',NK_SHIFT,' k-points corresponding to difference vectors'
         ENDIF
         WRITE(IU6,'(A)') ' '
         WRITE(IU6,'(A)') ' Following reciprocal coordinates:'
         WRITE(IU6,'(A)') '            Coordinates               Weight'
         DO IK=1,NKPT
            WRITE(IU6,'(3F10.6,4X,F10.6)') VKPT(1,IK),VKPT(2,IK),VKPT(3,IK),WTKPT(IK)
         ENDDO
         WRITE(IU6,'(A)') ' '
         WRITE(IU6,'(A)') ' Following cartesian coordinates:'
         WRITE(IU6,'(A)') '            Coordinates               Weight'
         DO IK=1,NKPT
            V=MATMUL(B,VKPT(:,IK))*SCALE
            WRITE(IU6,'(3F10.6,4X,F10.6)') V,WTKPT(IK)
         ENDDO
         WRITE(IU6,'(A)') ' '
      ENDIF

      ! Result on file IBZKPT (can be copied to KPOINTS for following runs):
      IF (IU6>=0) THEN
         IU=-1
         DO I=0,99
            IF ((I==0).OR.(I==5).OR.(I==6)) CYCLE
            INQUIRE(UNIT=I,OPENED=OCCUP)
            IF (.NOT.OCCUP) THEN
               IU=I; EXIT
            ENDIF
         ENDDO
         IF (IU<0) THEN
            WRITE(*,*) &
           &      'Sorry: No IO-unit available! Cannot save data to IBZKPT!'
         ELSE
            OPEN(UNIT=IU,FILE=DIR_APP(1:DIR_LEN)//'IBZKPT')
            WRITE(IU,'(A)') 'Automatically generated mesh'
            WRITE(IU,'(I8)') NKPT
            WRITE(IU,'(A)') 'Reciprocal lattice'
            DO IK=1,NKPT
               WRITE(IU,'(3F20.14,4X,I10)') VKPT(1,IK),VKPT(2,IK),VKPT(3,IK),NINT(WTKPT(IK))
            ENDDO
            IF (LTET) THEN
               WRITE(IU,'(A)') 'Tetrahedra'
               WRITE(IU,'(I10,F20.14)') NTET,VOLWGT
               DO ITET=1,NTET
                  WRITE(IU,'(5I10)') (IDTET(KTET,ITET),KTET=0,4)
               ENDDO
            ENDIF
            CLOSE(IU)
         ENDIF
      ENDIF

      ! Correct scaling for weights ... :
      WTKPT=WTKPT/WEISUM

      RETURN
    END SUBROUTINE IBZKPT


!******************** SUBROUTINE TETIRR ********************************
!
!> Finds inequivalent tetrahedra in an equally spaced k-mesh
!> (setting also the correct weights) ...
!
!***********************************************************************

    SUBROUTINE TETIRR(NTET,IDTET,BK,ILSNF,ISNF,NKX,NKY,NKZ,FBZ2IBZ,IU6)
      USE prec

      INTEGER,INTENT(OUT) :: NTET !< total number of generated tetrahedra
      INTEGER,POINTER,INTENT(INOUT) :: IDTET(:,:) !< weighting factors [IDTET(0,...)] and k-index of the four corners [IDTET(1-4,...)]
      REAL(q),INTENT(IN) :: BK(3,3) !< basis vectors of the generating k-lattice in cartesian coordinates
      INTEGER,INTENT(IN) :: ILSNF(3,3) !< left transform of Smith normal form decomposition
      INTEGER,INTENT(IN) :: ISNF(3,3) !< Smith normal form matrix
      INTEGER,INTENT(IN) :: NKX !< number of k-point divisions in the x direction
      INTEGER,INTENT(IN) :: NKY !< number of k-point divisions in the y direction
      INTEGER,INTENT(IN) :: NKZ !< number of k-point divisions in the z direction
      INTEGER,INTENT(IN) :: FBZ2IBZ(:,:) !< array mapping the integer points to the irreducible briullouin zone
      INTEGER,INTENT(IN) :: IU6 !< unit wehre to write informations
      ! local variables
      INTEGER :: NPX, NPY, NPZ
      INTEGER :: KINT(3)
      INTEGER KCUT(3,4,6), IQ(4)
      INTEGER IC, ITET, I, J, I1, I2, I3, K1, K2, K3
      INTEGER J1, J2, J3, II, N, NSIZE, IKFBZ

      INTEGER, ALLOCATABLE :: IDTMP(:,:)

      ! these are the multiplicities required to generate all the points in the BZ
      ! and are taken from the diagonal of the Smith normal form matrix
      NPX=ISNF(1,1)
      NPY=ISNF(2,2)
      NPZ=ISNF(3,3)
! Get integer displacements corresponding to 6 tetrahedra with the smallest diagonal
      CALL GET_TETDISP(BK,KCUT)
! Transform the kcut for non-generalized k-grids
      DO ITET=1,6
       DO IC=1,4
        KCUT(:,IC,ITET)=MATMUL(ILSNF,KCUT(:,IC,ITET))
       ENDDO
      ENDDO
! Now start searching the tetrahedra ... :
      NTET=0
! For all k-points ...
      DO I3=1,NPZ
       DO I2=1,NPY
        DO I1=1,NPX
! From this microcell we can cut out six tetrahedra:
         TET_LOOP: DO ITET=1,6
! Set the 4 corners (identifiers) of the actual tetrahedron:
            DO IC=1,4
               KINT = KIND_MOD_LOW(NPX,NPY,NPZ,[I1+KCUT(1,IC,ITET),I2+KCUT(2,IC,ITET),I3+KCUT(3,IC,ITET)])
               IKFBZ = KINT_KINDEX_LOW(NPX,NPY,NPZ,KINT)
               IQ(IC)=FBZ2IBZ(1,IKFBZ)
            ENDDO
! Order the identifiers of the corners ...
            DO J=1,3
             DO I=1,4-J
               IF (IQ(I)>IQ(I+1)) THEN
                  II=IQ(I)
                  IQ(I)=IQ(I+1)
                  IQ(I+1)=II
               END IF
             ENDDO
            ENDDO
! Now test all tetrahedra found previously:
            DO N=1,NTET
               IF (ALL(IDTET(1:4,N)==IQ(:))) THEN
! We have found the same combination previously, so increment the
! counter for this type of tetrahedron ...
                  IDTET(0,N)=IDTET(0,N)+1
! ... and go to the next tetrahedron:
                  CYCLE TET_LOOP
               END IF
            ENDDO
! New tetrahedron found if arriving here:
! Count it, ...
            NTET=NTET+1
            IF (NTET>SIZE(IDTET,2)) THEN
               ! store IDTET in temporary array
               ALLOCATE(IDTMP(0:4,SIZE(IDTET,2))); IDTMP=IDTET
               ! destroy IDTET
               DEALLOCATE(IDTET); NULLIFY(IDTET)
               ! reallocate IDTET with larger second dimension
               NSIZE = SIZE(IDTMP,2)*(2-REAL((I3-1)*NPX*NPY + (I2-1)*NPY + (I1-1),KIND=q)/REAL(NPX*NPY*NPZ,KIND=q)) + 1000
               ALLOCATE(IDTET(0:4,NSIZE))
               ! copy back to IDTET
               IDTET(:,1:SIZE(IDTMP,2))=IDTMP
               ! destroy temporary array
               DEALLOCATE(IDTMP)
            ENDIF
            IF (NTET>SIZE(IDTET,2)) CALL ERROR(' IBZKPT',' NTET>SIZE(IDTET)',NTET-SIZE(IDTET,2))
! ... store the corner coordinates (identifier) ...
            DO I=1,4
               IDTET(I,NTET)=IQ(I)
            ENDDO
! ... and initialize the counter for this new type of tetrahedron:
            IDTET(0,NTET)=1
         ENDDO TET_LOOP
        ENDDO
       ENDDO
      ENDDO
! Now tell us the result ... :
      IF (IU6>=0) THEN
        WRITE(IU6,'(1X,A,I6,A,I8)') 'TETIRR: Found ',NTET,' inequivalent tetrahedra from ',&
                                     6*NKX*NKY*NKZ
      ENDIF
    END SUBROUTINE TETIRR

    !> @brief Compute irreducible tetrahedron using the star of each IBZ points
    SUBROUTINE TETIRR_STAR(KPOINTS,SYMM,LATT_CUR,NTET,IDTET,IU6)
      use base, only: symmetry
      use tutor, only: vtutor
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(symmetry),INTENT(IN) :: SYMM
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: IU6 !< Unit to write output
      INTEGER,INTENT(OUT) :: NTET
      INTEGER,ALLOCATABLE,INTENT(OUT) :: IDTET(:,:)
      ! local variables
      INTEGER,PARAMETER :: NTET_STEP=24
      INTEGER,PARAMETER :: MAX_NTET_IBZ=144 !< maximum number of tetrahedra that one IBZ k-point can create
      INTEGER :: I1,I2,I3,ITET,IC
      INTEGER :: IX,IY,IZ,ISYM
      INTEGER :: IKIBZ, THIS_IKIBZ, IKPT, IDX
      INTEGER :: IKSTAR, NKSTAR, NTET_IBZ, NSIZE
      INTEGER :: KCUT24(3,4,24) !< integers defining 24 tetrahedra that contribute to each IBZ point
      INTEGER :: TETRA_IBZ(4), ORDER(4)
      INTEGER :: KINT(3), KINT_ROT(3)
      INTEGER :: KINT_STAR(3,48) !< k-point star of integers
      INTEGER :: TET_IBZ(4,MAX_NTET_IBZ) !< max list of tetrahedron contributing to one IBZ
      INTEGER :: TETW_IBZ(MAX_NTET_IBZ) !< weight of each tetrahedron in the ibz
      INTEGER :: IOP(3,3)
      REAL(q) :: VKPT(3), VKPT_NEW(3)
      REAL(q) :: ROP(3,3)
      INTEGER :: NKFULL, IKFBZ
      INTEGER,ALLOCATABLE :: IDTMP(:,:)
      ! initial allocation for IDTET
      ALLOCATE(IDTET(0:4,NTET_STEP*KPOINTS%NKPTS))
      NTET=0
      NKFULL=0
      ! get 24 tetrahedra contributing to any k-point
      CALL GET_TETDISP_24TET(KPOINTS%B,KCUT24)
      ! transform kcut
      DO ITET=1,24
        DO IC=1,4
          KCUT24(:,IC,ITET) = MATMUL(KPOINTS%ILSNF,KCUT24(:,IC,ITET))
        ENDDO
      ENDDO

      ! loop over all the irreducible k-points
      DO IKIBZ=1,KPOINTS%NKPTS
        VKPT = KPOINTS%VKPT(:,IKIBZ)
        !CALL KPOINT_INTEGERS(KPOINTS,B,VKPT,KINT)

        ! determine all the unique k-points in the star
        NKSTAR=0
        DO ISYM=1,SYMM%NROT
          ROP = SYMM%IGRPOP(:,:,ISYM)
          VKPT_NEW(1)=VKPT(1)*ROP(1,1)+VKPT(2)*ROP(2,1)+VKPT(3)*ROP(3,1)
          VKPT_NEW(2)=VKPT(1)*ROP(1,2)+VKPT(2)*ROP(2,2)+VKPT(3)*ROP(3,2)
          VKPT_NEW(3)=VKPT(1)*ROP(1,3)+VKPT(2)*ROP(2,3)+VKPT(3)*ROP(3,3)
          call KPOINT_INTEGERS(KPOINTS,LATT_CUR,VKPT_NEW,KINT)

          ! Find this k-point in the list of k-points in the star
          IDX = IDX_LIST(3,KINT,NKSTAR,KINT_STAR)
          ! If this point does not exist in the star list yet, then store it
          IF (IDX==-1) THEN
            NKSTAR=NKSTAR+1
            KINT_STAR(:,NKSTAR)=KINT
          ENDIF
        ENDDO !ISYM
        NKFULL=NKFULL+NKSTAR

        ! create a list that contains all the tetrahedron that contribute to this IBZ point
        NTET_IBZ = 0
        DO IKSTAR=1,NKSTAR
          ! loop over tetrahedra that can contribute to this k-point
          TET_LOOP: DO ITET=1,24
            ! loop over the corners of the tetrahedra
            DO IC=1,4
              KINT = KINT_MOD(KPOINTS,KINT_STAR(:,IKSTAR)+KCUT24(:,IC,ITET))
              IKFBZ = KINT_KINDEX(KPOINTS,KINT)
              ! store the ibz index point corresponding to this corner
              THIS_IKIBZ = KPOINTS%FBZ2IBZ(1,IKFBZ)
              ! if the IBZ point corresponding to this corner is lower than the IBZ
              ! we used to generate the 24 tetrahedra then this tetrahedron was already
              ! included in the list so we skip it
              IF (THIS_IKIBZ<IKIBZ) CYCLE TET_LOOP
              TETRA_IBZ(IC) = THIS_IKIBZ
            ENDDO ! IC
            ! Sort the indexes of corners of the tetrahedra
            CALL SORT_TETRA(TETRA_IBZ,ORDER)
            ! we only keep the tetrahedra whose lower index is equal to the current IBZ point
            IF (TETRA_IBZ(1)/=IKIBZ) CYCLE TET_LOOP
            ! check if this tetrahedra already exists in the list
            IDX = IDX_LIST(4,TETRA_IBZ,NTET_IBZ,TET_IBZ)
            ! if it does not exist then add it
            IF (IDX==-1) THEN
              NTET_IBZ=NTET_IBZ+1
              IF (NTET_IBZ>MAX_NTET_IBZ) CALL VTUTOR%BUG('Need more tetrahedra',__FILE__,__LINE__)
              TET_IBZ(:,NTET_IBZ)=TETRA_IBZ
              TETW_IBZ(NTET_IBZ)=1
            ! if it exists then increase multiplicity
            ELSE
              TETW_IBZ(IDX)=TETW_IBZ(IDX)+1
            ENDIF
          ENDDO TET_LOOP
        ENDDO ! IKSTAR

        ! add the tetrahedron contributing to this IBZ to the main list
        DO ITET=1,NTET_IBZ
          NTET = NTET+1
          ! If the size of the list of tetrahedra is not large enough, then increase it
          NSIZE =  SIZE(IDTET,2)
          IF (NSIZE<NTET) THEN
            ALLOCATE(IDTMP(0:4,NSIZE+NTET_STEP*KPOINTS%NKPTS))
            IDTMP(0:4,:NSIZE)=IDTET(0:4,:)
            CALL MOVE_ALLOC(IDTMP,IDTET)
          ENDIF
          IDTET(0,NTET)=TETW_IBZ(ITET)
          IDTET(1:4,NTET)=TET_IBZ(:,ITET)
        ENDDO
      ENDDO ! IKIBZ
      IF (IU6>=0) THEN
        WRITE(IU6,'(1X,A,I6,A,I8)') 'TETIRR_STAR: Found ',NTET,' inequivalent tetrahedra from ',&
                                     6*KPOINTS%NKPX*KPOINTS%NKPY*KPOINTS%NKPZ
      ENDIF

      CONTAINS
      !> @brief Find an array of N integers K(N) on a list with NL arrays of N elements
      INTEGER FUNCTION IDX_LIST(N,K,NL,LIST)
        INTEGER :: N !< size of the array
        INTEGER :: K(N) !< the array
        INTEGER :: NL !< size of the array of arrays
        INTEGER :: LIST(N,NL) !< an array of arrays
        ! local variables
        DO IDX_LIST=1,NL
           IF (ALL(LIST(:,IDX_LIST)==K)) RETURN
        ENDDO
        IDX_LIST = -1
      END FUNCTION IDX_LIST
    END SUBROUTINE TETIRR_STAR

    !> @brief Sort the tetrahedra according to the indexes representing each corner
    PURE SUBROUTINE SORT_TETRA(LIST,PERM)

     INTEGER, INTENT(INOUT) :: PERM(4)
     INTEGER, INTENT(INOUT) :: LIST(4)

     INTEGER :: IA,IB,IC,ID
     INTEGER :: ILOW1,ILOW2,IHIGH1,IHIGH2
     INTEGER :: ILOWEST,IHIGHEST
     INTEGER :: IMIDDLE1,IMIDDLE2
     REAL(q) :: VA,VB,VC,VD
     REAL(q) :: VLOW1,VLOW2,VHIGH1,VHIGH2
     REAL(q) :: VLOWEST,VHIGHEST
     REAL(q) :: VMIDDLE1,VMIDDLE2

     VA = LIST(1); IA = PERM(1)
     VB = LIST(2); IB = PERM(2)
     VC = LIST(3); IC = PERM(3)
     VD = LIST(4); ID = PERM(4)

     IF (VA < VB) THEN
         VLOW1 = VA; VHIGH1 = VB
         ILOW1 = IA; IHIGH1 = IB
     ELSE
         VLOW1 = VB; VHIGH1 = VA
         ILOW1 = IB; IHIGH1 = IA
     ENDIF

     IF (VC < VD) THEN
         VLOW2 = VC; VHIGH2 = VD
         ILOW2 = IC; IHIGH2 = ID
     ELSE
         VLOW2 = VD; VHIGH2 = VC
         ILOW2 = ID; IHIGH2 = IC
     ENDIF

     IF (VLOW1 < VLOW2) THEN
         VLOWEST  = VLOW1; VMIDDLE1 = VLOW2
         ILOWEST  = ILOW1; IMIDDLE1 = ILOW2
     ELSE
         VLOWEST  = VLOW2; VMIDDLE1 = VLOW1
         ILOWEST  = ILOW2; IMIDDLE1 = ILOW1
     ENDIF

     IF (VHIGH1 > VHIGH2) THEN
         VHIGHEST = VHIGH1; VMIDDLE2 = VHIGH2
         IHIGHEST = IHIGH1; IMIDDLE2 = IHIGH2
     ELSE
         VHIGHEST = VHIGH2; VMIDDLE2 = VHIGH1
         IHIGHEST = IHIGH2; IMIDDLE2 = IHIGH1
     ENDIF

     IF (VMIDDLE1 < VMIDDLE2) THEN
         LIST = [VLOWEST,VMIDDLE1,VMIDDLE2,VHIGHEST]
         PERM = [ILOWEST,IMIDDLE1,IMIDDLE2,IHIGHEST]
     ELSE
         LIST = [VLOWEST,VMIDDLE2,VMIDDLE1,VHIGHEST]
         PERM = [ILOWEST,IMIDDLE2,IMIDDLE1,IHIGHEST]
     ENDIF

    END SUBROUTINE SORT_TETRA

    !> @brief Get integer displacements corresponding to 6 tetrahedra with the smallest diagonal
    SUBROUTINE GET_TETDISP(BK,KCUT)
      REAL(q),INTENT(IN) :: BK(3,3) !< coordinates of the generating lattice
      INTEGER,INTENT(OUT) :: KCUT(3,4,6) !< integers defining 6 tetrahedra
      ! local variables
      REAL(q) :: EDMIN,EDMAX,EDGMIN,EDGMAX,XX
      REAL(q) :: P(3,4)
      INTEGER :: ITET,IC,I,J
      INTEGER :: LX, LY, LZ
      INTEGER :: LXX, LYY
      INTEGER :: IDTET

      INTEGER :: KCUT0(3,4,6)
      SAVE KCUT0
      DATA KCUT0/ 0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1, &
                  0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1, &
                  0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /

      ! Setting up the tetrahedra will be done cutting a microcell with eight
      ! corners into six tetrahedra. The edges of the tetrahedra are given by
      ! three edges, two face diagonals and one space diagonal of the cell.
      ! Giving the space diagonal, the way how to choose the rest is uniquely
      ! determined ... . But now there are four different possibilities how to
      ! choose the space diagonal! Prefer the choice which gives the shortest
      ! edges for all tetrahedra ('the most compact tetrahedra') - just to
      ! avoid very long 'interpolation distances' ... !
      LZ=0
      LXX=0
      LYY=0
      EDGMAX=1.E30_q
      EDGMIN=0._q
! For the four choices ...
      DO LX=0,1
       DO LY=0,1
! ... we set up the 'trial division' of a given cell into 6 tetrahedra:
         DO ITET=1,6
          DO IC=1,4
            KCUT(1,IC,ITET)=KCUT0(1,IC,ITET)
            KCUT(2,IC,ITET)=KCUT0(2,IC,ITET)
            KCUT(3,IC,ITET)=KCUT0(3,IC,ITET)
            IF (LX==1) KCUT(1,IC,ITET)=1-KCUT0(1,IC,ITET)
            IF (LY==1) KCUT(2,IC,ITET)=1-KCUT0(2,IC,ITET)
          ENDDO
         ENDDO
         EDMIN=1.E30_q
         EDMAX=0._q
! For this trial setting, loop over all tetrahedra ...,
         DO ITET=1,6
! ... set up the cartesian coordinates of the four corner points ...,
            DO IC=1,4
               P(1,IC)=KCUT(1,IC,ITET)*BK(1,1)+ &
     &                 KCUT(2,IC,ITET)*BK(1,2)+ &
     &                 KCUT(3,IC,ITET)*BK(1,3)
               P(2,IC)=KCUT(1,IC,ITET)*BK(2,1)+ &
     &                 KCUT(2,IC,ITET)*BK(2,2)+ &
     &                 KCUT(3,IC,ITET)*BK(2,3)
               P(3,IC)=KCUT(1,IC,ITET)*BK(3,1)+ &
     &                 KCUT(2,IC,ITET)*BK(3,2)+ &
     &                 KCUT(3,IC,ITET)*BK(3,3)
            ENDDO
! ... and get the shortest and longest distance between two points in
! each tetrahedron (minimum/maximum taken over all tetrahedra ...):
            DO I=1,3
             DO J=I+1,4
               XX=ANRM2(P(1,I)-P(1,J),P(2,I)-P(2,J),P(3,I)-P(3,J))
               EDMAX=MAX(EDMAX,XX)
               EDMIN=MIN(EDMIN,XX)
             ENDDO
            ENDDO
         ENDDO
! Now look at the global maximum: Have we found a cut with smaller
! maximum distance between two points within one tetrahedron than
! before? If yes: store it  (until we find something better ...)!
         IF (EDMAX<EDGMAX) THEN
            LXX=LX
            LYY=LY
            EDGMAX=EDMAX
            EDGMIN=EDMIN
         END IF
       ENDDO
      ENDDO
! Now set up the 'correct' cutup giving the most compact tetrahdra ... :
      DO ITET=1,6
       DO IC=1,4
         KCUT(1,IC,ITET)=KCUT0(1,IC,ITET)
         KCUT(2,IC,ITET)=KCUT0(2,IC,ITET)
         KCUT(3,IC,ITET)=KCUT0(3,IC,ITET)
         IF (LXX==1) KCUT(1,IC,ITET)=1-KCUT0(1,IC,ITET)
         IF (LYY==1) KCUT(2,IC,ITET)=1-KCUT0(2,IC,ITET)
       ENDDO
      ENDDO
      CONTAINS
      REAL(q) FUNCTION ANRM2(X,Y,Z)
          REAL(q) X, Y, Z
          ANRM2=X*X*1.00001E0_q+&
                Y*Y*1.00002E0_q+&
                Z*Z*1.00003E0_q-&
                X*0.000004E0_q-&
                Y*0.000003E0_q-&
                Z*0.000002E0_q
      END FUNCTION
    END SUBROUTINE GET_TETDISP

    !> @brief In the FBZ there are 24 tetrahedra contributing to each k-point
    !> get the integer displacements generating those 24 tetrahedra
    !> For this we displace the microzone of 6 tetrahedron around the (0,0,0)
    !> integer point and store the tetrahedron that contain it
    SUBROUTINE GET_TETDISP_24TET(BK,KCUT24)
      use tutor, only: vtutor
      REAL(q) :: BK(3,3)
      INTEGER :: KCUT24(3,4,24)
      !local variables
      LOGICAL :: STORE
      INTEGER :: IX, IY, IZ
      INTEGER :: IC,ITET
      INTEGER :: ITETCOUNT ! count tetrahedron up to 24
      INTEGER :: TETDISP(3) ! displacement to apply to the 6 tetrahedrons
      INTEGER :: KCUT(3,4,6) ! 6 tetrahedrons with the best possible cut
      INTEGER :: TET(3,4) ! integers representing a tetrahedron

      CALL GET_TETDISP(BK,KCUT)

      ITETCOUNT=0
      ! apply 8 displacements
      DO IX=0,1
        TETDISP(1)=IX
        DO IY=0,1
          TETDISP(2)=IY
          DO IZ=0,1
            TETDISP(3)=IZ
            ! loop over all the 6 tetrahedron
            DO ITET=1,6
              ! apply displacement to each corner of each tetrahedron
              STORE=.FALSE.
              DO IC=1,4
                  TET(:,IC)=KCUT(:,IC,ITET)-TETDISP
                  ! if any of the corners is (0,0,0) then we store this tetrahedron
                  IF (ALL(TET(:,IC)==0)) STORE=.TRUE.
              ENDDO
              IF (STORE) THEN
                ITETCOUNT=ITETCOUNT+1
                IF (ITETCOUNT>24) THEN
                  CALL vtutor%bug('Number of tetrahedra contibuting to one k-point > 24',__FILE__,__LINE__)
                ENDIF
                KCUT24(:,:,ITETCOUNT)=TET(:,:)
              ENDIF
            ENDDO
          ENDDO !IZ
        ENDDO !IY
      ENDDO !IX
      IF (ITETCOUNT/=24) THEN
        CALL VTUTOR%BUG('Number of tetrahedra contibuting to one k-point /= 24',__FILE__,__LINE__)
      ENDIF
    END SUBROUTINE GET_TETDISP_24TET

    !> Get integer mapping corresponding to this k-point
    !> The integers are always between 0 and NKPX-1,NKPY-1,NKPZ-1 using the periodicity
    !>
    !> KINT -> 3 integers indentifying the k-point
    !> KINDEX -> Hash value computed from kint (recommended iterator: ik,ikp,iq)
    !> VKPT -> fractional coordinates of the k-point
    !> IBZ -> index of vkpt point in the list of irreducible zone (recommended iterator: iik,iikp,iiq)
    PURE SUBROUTINE KPOINT_INTEGERS(KPOINTS,LATT_CUR,VKPT,KINT,EPS)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      REAL(q),INTENT(IN) :: VKPT(3) !< coordinates of the k-point in fractional coordinated
      INTEGER,INTENT(OUT) :: KINT(3) !< integers corresponding to the given k-point
      REAL(q),OPTIONAL,INTENT(IN) :: EPS !< error to check if k-point matches in terms of the reciprocal basis
      ! local variables
      REAL(q) :: T(3),SHIFTK(3)
      REAL(q) :: MY_EPS

      MY_EPS=1e-6_q; IF (PRESENT(EPS)) MY_EPS = EPS
      SHIFTK = MATMUL(TRANSPOSE(KPOINTS%A),MATMUL(KPOINTS%B,KPOINTS%SHIFT))
      ! to cartesian coordinates
      T = MATMUL(LATT_CUR%B,VKPT)
      ! to the bk basis
      T = MATMUL(TRANSPOSE(KPOINTS%A),T)
      ! subtract shift
      T = T-SHIFTK
      IF (ANY(ABS(T-NINT(T))>MY_EPS)) THEN
        KINT=-1
        RETURN
      ENDIF
      ! round and mod
      KINT = KINT_MOD(KPOINTS,MATMUL(KPOINTS%ILSNF,NINT(T)))
    END SUBROUTINE KPOINT_INTEGERS

    !> Bring each of the there integers corresponding to a k-point
    !> back to [0,NPX] [0,NPY] and [0,NPZ]
    PURE FUNCTION KINT_MOD(KPOINTS,KINT) RESULT(KOUT)
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      INTEGER,INTENT(IN) :: KINT(3)
      INTEGER :: KOUT(3)
      KOUT=KIND_MOD_LOW(KPOINTS%ISNF(1,1),KPOINTS%ISNF(2,2),KPOINTS%ISNF(3,3),KINT)
    END FUNCTION KINT_MOD

    PURE FUNCTION KIND_MOD_LOW(NPX,NPY,NPZ,KINT) RESULT(KOUT)
      INTEGER,INTENT(IN) :: NPX
      INTEGER,INTENT(IN) :: NPY
      INTEGER,INTENT(IN) :: NPZ
      INTEGER,INTENT(IN) :: KINT(3)
      INTEGER :: KOUT(3)
      KOUT=MODULO(KINT,[NPX,NPY,NPZ])
      !This is equivalent to the above
      !KOUT=MOD(KINT,[NPX,NPY,NPZ])
      !IF (KOUT(1)<0) KOUT(1)=KOUT(1)+NPX
      !IF (KOUT(2)<0) KOUT(2)=KOUT(2)+NPY
      !IF (KOUT(3)<0) KOUT(3)=KOUT(3)+NPZ
    END FUNCTION KIND_MOD_LOW

    !> Compute KINT from KINDEX
    !> We have KINDEX = KINT(3)*NPY*NPX+&
    !>                  KINT(2)*NPX+&
    !>                  KINT(1)+1
    PURE FUNCTION KINDEX_KINT(KPOINTS,KINDEX) RESULT(KINT)
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      INTEGER,INTENT(IN) :: KINDEX
      INTEGER :: KINT(3)
      !local variables
      INTEGER :: NPX, NPY, NPZ
      NPX = KPOINTS%ISNF(1,1)
      NPY = KPOINTS%ISNF(2,2)
      NPZ = KPOINTS%ISNF(3,3)
      KINT(1) = MODULO(KINDEX-1,NPX)
      KINT(3) = (KINDEX-1)/(NPY*NPX)
      KINT(2) = (KINDEX-1-KINT(3)*NPY*NPX)/NPX
      KINT = KINT_MOD(KPOINTS,KINT)
    END FUNCTION KINDEX_KINT

    !> Convert from three integers denoting a k-point to a single integer
    PURE FUNCTION KINT_KINDEX(KPOINTS,KINT) RESULT(KINDEX)
      use tutor, only: vtutor
      use string, only: str
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      INTEGER,INTENT(IN) :: KINT(3)
      INTEGER :: KINDEX
      ! local variables
      INTEGER :: NPX,NPY,NPZ
      integer :: KINT_COMP(3)
      NPX = KPOINTS%ISNF(1,1)
      NPY = KPOINTS%ISNF(2,2)
      NPZ = KPOINTS%ISNF(3,3)
      KINDEX = KINT_KINDEX_LOW(NPX,NPY,NPZ,KINT)

      ! check if from kindex I get original KINT
      !KINT_COMP = KINDEX_KINT(KPOINTS,KINDEX)
      !IF (ANY(KINT_COMP/=KINT)) THEN
      !  call vtutor%bug('The mapping from KINDEX to KINT '//str(KINT_COMP)//&
      !                  ' is not the same as KINT '//str(KINT)//' to KINDEX '//str(KINDEX),__FILE__,__LINE__)
      !ENDIF
    END FUNCTION KINT_KINDEX

    PURE FUNCTION KINT_KINDEX_LOW(NPX,NPY,NPZ,KINT) RESULT(KINDEX)
      INTEGER,INTENT(IN) :: NPX
      INTEGER,INTENT(IN) :: NPY
      INTEGER,INTENT(IN) :: NPZ
      INTEGER,INTENT(IN) :: KINT(3)
      INTEGER :: KINDEX
      KINDEX = KINT(3)*NPY*NPX+&
               KINT(2)*NPX+&
               KINT(1)+1
    END FUNCTION KINT_KINDEX_LOW

    !> use the indeces of q and k' to determine the index of k
    !> The relation between the coordinates is:
    !>   K - KP = Q + G
    FUNCTION KINDEX_FROM_K_AND_KP(KPOINTS,LATT_CUR,IKFBZ,IKPFBZ) RESULT(KINDEX)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: IKFBZ
      INTEGER,INTENT(IN) :: IKPFBZ
      INTEGER :: KINDEX
      ! local variables
      REAL(q) :: VK(3), VQ(3), VKP(3)
      VK = KINDEX_VKPT(KPOINTS,LATT_CUR,IKFBZ)
      VKP  = KINDEX_VKPT(KPOINTS,LATT_CUR,IKPFBZ)
      VQ = VK-VKP
      ! TODO: the G vector is removed in this routine.
      ! we could compute it and return it if needed
      KINDEX = VKPT_KINDEX(KPOINTS,LATT_CUR,VQ)
    END FUNCTION KINDEX_FROM_K_AND_KP

    !> use the indeces of q and k' to determine the index of k
    !> The relation between the coordinates is:
    !>   K = Q + KP + G
    FUNCTION KINDEX_FROM_KP_AND_Q(KPOINTS,LATT_CUR,IKPFBZ,IQFBZ) RESULT(KINDEX)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: IKPFBZ
      INTEGER,INTENT(IN) :: IQFBZ
      INTEGER :: KINDEX
      ! local variables
      REAL(q) :: VK(3), VQ(3), VKP(3)
      VKP = KINDEX_VKPT(KPOINTS,LATT_CUR,IKPFBZ)
      VQ  = KINDEX_VKPT(KPOINTS,LATT_CUR,IQFBZ)
      VK = VQ+VKP
      ! TODO: the G vector is removed in this routine.
      ! we could compute it and return it if needed
      KINDEX = VKPT_KINDEX(KPOINTS,LATT_CUR,VK)
    END FUNCTION KINDEX_FROM_KP_AND_Q

    !> Get index of k point from its fractional coordinate
    FUNCTION VKPT_KINDEX(KPOINTS,LATT_CUR,VKPT) RESULT(KINDEX)
      use tutor, only: vtutor
      use string, only: str
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      REAL(q),INTENT(IN) :: VKPT(3)
      INTEGER :: KINDEX
      ! local variables
      INTEGER :: KINT(3)
      CALL KPOINT_INTEGERS(KPOINTS,LATT_CUR,VKPT,KINT)
      IF (ANY(KINT==-1)) THEN
         call vtutor%error('Failed to map k-point: '//str(VKPT,.false.)//' to integers')
      ENDIF
      KINDEX = KINT_KINDEX(KPOINTS,KINT)
    END FUNCTION VKPT_KINDEX

    !> Get coordinates of the k-point from integer indexes
    FUNCTION KINT_VKPT(KPOINTS,LATT_CUR,KINT) RESULT(VKPT)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: KINT(3)
      REAL(q) :: VKPT(3)
      !local variables
      REAL(q),PARAMETER :: TINY=1e-6_q
      INTEGER :: KINT_(3)
      PROFILING_START('kint_vkpt')
      KINT_ = MATMUL(KPOINTS%ILSNFINV,KINT)
      VKPT(1) = (KINT_(1)+KPOINTS%SHIFT(1))
      VKPT(2) = (KINT_(2)+KPOINTS%SHIFT(2))
      VKPT(3) = (KINT_(3)+KPOINTS%SHIFT(3))
      VKPT = MATMUL(KPOINTS%B,VKPT)
      VKPT = MATMUL(TRANSPOSE(LATT_CUR%A),VKPT)
      ! TODO: This is written in this way such as to be consistent with what is done
      ! in IBZKPT. This shuold be changed however
      VKPT = MOD(VKPT+60.5_q-0.5_q*TINY,1._q)-0.5_q+0.5_q*TINY
      PROFILING_STOP('kint_vkpt')
    END FUNCTION KINT_VKPT

    !> Get kindex from the fractional coordinates of the k-point
    FUNCTION KINDEX_VKPT(KPOINTS,LATT_CUR,KINDEX) RESULT(VKPT)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: KINDEX
      REAL(q) :: VKPT(3)
      ! local variables
      INTEGER :: KINT(3)
      KINT = KINDEX_KINT(KPOINTS,KINDEX)
      VKPT = KINT_VKPT(KPOINTS,LATT_CUR,KINT)
    END FUNCTION KINDEX_VKPT

    !> Get KINDEX for the IKIBZ point in the IBZ k-points list
    FUNCTION IKIBZ_KINDEX(KPOINTS,LATT_CUR,IKIBZ) RESULT(KINDEX)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: IKIBZ
      INTEGER :: KINDEX
      KINDEX = VKPT_KINDEX(KPOINTS,LATT_CUR,KPOINTS%VKPT(:,IKIBZ))
    END FUNCTION IKIBZ_KINDEX

    !> Get KINDEX for the IKIBZ point cooresponding to a point in the FBZ (IKFBZ)
    FUNCTION IKFBZ_KINDEX(KPOINTS,LATT_CUR,IKFBZ) RESULT(KINDEX)
      use lattice, only: latt
      TYPE(kpoints_struct),INTENT(IN) :: KPOINTS !< kpoints structure
      TYPE(latt),INTENT(IN) :: LATT_CUR !< lattice structure
      INTEGER,INTENT(IN) :: IKFBZ
      INTEGER :: KINDEX
      !local variables
      INTEGER :: IKIBZ
      IKIBZ = KPOINTS%FBZ2IBZ(1,IKFBZ)
      KINDEX = IKIBZ_KINDEX(KPOINTS,LATT_CUR,IKIBZ)
    END FUNCTION IKFBZ_KINDEX

    !> This routine takes as an argument the index of a k-point in the FBZ and returns the corresponding point in the IBZ + the rotation
    !> that when applied to the IBZ k-point generates the point in the FBZ.
    SUBROUTINE KPOINTS_FBZ_TO_IBZ(KPOINTS,IKFBZ,IKIBZ,IROT)
      TYPE(KPOINTS_STRUCT),TARGET,INTENT(IN) :: KPOINTS
      INTEGER,INTENT(IN) :: IKFBZ !< index of the k-point in the FBZ
      INTEGER,INTENT(OUT) :: IKIBZ !< index of the k-point in the IBZ
      INTEGER,INTENT(OUT) :: IROT !< index of the rotation that expands the k-point from the IBZ to the FBZ
      ! local variables
      integer :: IROT_
      IKIBZ = KPOINTS%FBZ2IBZ(1,IKFBZ)
      IROT_ = KPOINTS%FBZ2IBZ(2,IKFBZ)
      IROT  = SIGN(1,IROT_)*KPOINTS%INVMAP(ABS(IROT_))
    END SUBROUTINE KPOINTS_FBZ_TO_IBZ

    !> @brief Get IBZ indexes of the corners of 24 tetrahedra
    SUBROUTINE GET_INDEXES(KPOINTS,KINT,KCUT24,IDXS)
        TYPE(KPOINTS_STRUCT),INTENT(IN) :: KPOINTS
        INTEGER,INTENT(IN) :: KINT(3)
        INTEGER,INTENT(IN) :: KCUT24(3,4,24)
        INTEGER,INTENT(OUT) :: IDXS(4,24)
        ! local variables
        INTEGER :: PINT(3)
        INTEGER :: ITET,IC,IKIBZ,IKFBZ
        DO ITET=1,24
            DO IC=1,4
                PINT = KINT_MOD(KPOINTS,KINT+KCUT24(:,IC,ITET))
                IKFBZ = KINT_KINDEX(KPOINTS,PINT)
                IKIBZ = KPOINTS%FBZ2IBZ(1,IKFBZ)
                IDXS(IC,ITET) = IKIBZ
            ENDDO
        ENDDO
    END SUBROUTINE GET_INDEXES

    FUNCTION KINDEX_SNFORM(K,IS,ILS)
      USE prec
      INTEGER :: IS(3,3),ILS(3,3),KINDEX_SNFORM
      REAL(q) :: K(3)
      ! local variables
      REAL(q) :: TMP(3)
      TMP=MOD(MATMUL(ILS,NINT(K)),(/IS(1,1),IS(2,2),IS(3,3)/))
      IF (TMP(1)<-1.E-8_q) TMP(1)=TMP(1)+IS(1,1)
      IF (TMP(2)<-1.E-8_q) TMP(2)=TMP(2)+IS(2,2)
      IF (TMP(3)<-1.E-8_q) TMP(3)=TMP(3)+IS(3,3)
      KINDEX_SNFORM=TMP(3)*IS(2,2)*IS(1,1)+TMP(2)*IS(1,1)+TMP(1)+1
    END FUNCTION KINDEX_SNFORM

    !> Initialize some arrays needed for using the tetrahedron method
    subroutine kpoints_init_tetrahedron(kpoints,kcut24,kcut24_origin)
        type(kpoints_struct),intent(in) :: kpoints
        integer,intent(out) :: kcut24(:, :, :) !< displacements of 24 tetrahedra that contribute to a k-point
        integer,intent(out) :: kcut24_origin(:) !< index of the corner of each of the 24 tetrahedra for which the k-point has 0 displacement
        ! local variables
        integer :: ic, itet

        ! get tetrahedron displacements
        call get_tetdisp_24tet(kpoints%b,kcut24)
        ! get [0,0,0] point of each kcut
        do itet=1,24
            do ic=1,4
                if( all(kcut24(:,ic,itet)==[0,0,0])) exit
            enddo
            kcut24_origin(itet)=ic
        enddo
        ! transform kcut
        do itet=1,24
          do ic=1,4
            kcut24(:,ic,itet) = matmul(kpoints%ilsnf,kcut24(:,ic,itet))
          enddo
        enddo
    end subroutine kpoints_init_tetrahedron

    !> Set the number of integer divisions along x, y and z
    !> based on a desired k-spacing and the norm of the reciprocal space vectors
    subroutine kpoints_set_nkp_from_kspacing(kpoints,bnorm,kspacing)
      use constant, only: pi
      type(kpoints_struct),intent(inout) :: kpoints
      real(q) :: bnorm(3) !< norm of the reciprocal lattice vectors
      real(q) :: kspacing !< the maximum spacing between nearby k-points
      kpoints%nkpx=max(1,ceiling(bnorm(1)*pi*2/kspacing))
      kpoints%nkpy=max(1,ceiling(bnorm(2)*pi*2/kspacing))
      kpoints%nkpz=max(1,ceiling(bnorm(3)*pi*2/kspacing))
    end subroutine kpoints_set_nkp_from_kspacing

    !> Read a KPOINTS file in line mode and find the k-points in the kpoints FBZ
    !> which fall along the different line segments
    subroutine kpoints_along_path(kpoints,latt_cur,points_dists,points_coord,points_ikibz)
        !use mkpoints_struct_def, only: kpoints_struct
        !use mkpoints, only: rd_kpoints, copy_kpoints, kint_mod, kpoint_integers
        use poscar_struct_def, only: latt
        use m_mrgrnk, only: mrgrnk
        use tutor, only: vtutor
        type(kpoints_struct) :: kpoints
        type(latt) :: latt_cur
        integer, allocatable :: points_ikibz(:) !< contains the indexes of the points in the IBZ
        real(q), allocatable :: points_dists(:) !< contains the distance between the points
        real(q), allocatable :: points_coord(:,:) !< contains the coordinates of the points
        ! local variables
        type(kpoints_struct) :: kpoints_line
        integer :: i, j, k, ik, ib, isp, ik2, ikibz, ikfbz, n
        integer :: nkpt_line, nkpt_total
        integer :: ikx,iky,ikz
        integer :: kint(3)
        real(q) :: vkpt1(3), vkpt2(3), vkpt3(3)
        real(q) :: v1(3), v2(3), v3(3)
        real(q) :: dir1(3), dir2(3)
        real(q) :: v(3),vr(3)
        real(q) :: dist1, dist2, dist3, dist_total, dist_pline
        integer,parameter :: ndisp=7
        integer,parameter :: disp(7) = [0,1,-1,2,-2,3,-3]
        integer,parameter :: LIST_INIT=20
        integer,parameter :: LIST_STEP=10
        real(q),parameter :: tiny=1e-6
        real(q) :: t(3,3)
        real(q),allocatable :: vkpt_line(:,:)
        real(q),allocatable :: points_tmp_dists(:)
        real(q),allocatable :: points_tmp_coord(:,:)
        integer,allocatable :: points_tmp_ikibz(:)
        real(q),allocatable :: points_seg_dists(:)
        real(q),allocatable :: points_seg_coord(:,:)
        integer,allocatable :: points_seg_ikibz(:)
        integer,allocatable :: points_idx(:)
        logical,allocatable :: mask(:)

        PROFILING_START('kpoints_along_path')
        ! read kpoints_line
        call copy_kpoints(kpoints,kpoints_line)
        call rd_kpoints(kpoints_line,latt_cur,.false.,.false.,-1,-1,file_name='KPOINTS_LINE')

        t = matmul(transpose(latt_cur%a),kpoints%b)

        ! detect edges
        allocate(vkpt_line(3,kpoints_line%nkpts))
        nkpt_line = 1
        vkpt_line(:,nkpt_line) = kpoints_line%vkpt(:,1)
        do ik=2,kpoints_line%nkpts-1
            vkpt1=kpoints_line%vkpt(:,ik-1)
            vkpt2=kpoints_line%vkpt(:,ik)
            vkpt3=kpoints_line%vkpt(:,ik+1)
            dir1=vkpt2-vkpt1
            dir2=vkpt3-vkpt2
            if (any(abs(dir1-dir2)>1e-6)) then
                nkpt_line = nkpt_line + 1
                vkpt_line(:,nkpt_line) = kpoints_line%vkpt(:,ik)
            endif
        enddo
        nkpt_line = nkpt_line+1
        vkpt_line(:,nkpt_line) = kpoints_line%vkpt(:,kpoints_line%nkpts)

        ! write coordinates of segments
        !do ik=1,nkpt_line
        !    write(*,*) ik, vkpt_line(:,ik)
        !enddo

        nkpt_total=0
        dist_total=0
        allocate(points_dists(nkpt_total))
        allocate(points_coord(3,nkpt_total))
        allocate(points_ikibz(nkpt_total))
        ! loop over line segments detected above
        do ik=1,nkpt_line/2
           ! get extremeties of the line segment
           vkpt1 = vkpt_line(:,(ik-1)*2+1)
           vkpt3 = vkpt_line(:,(ik-1)*2+2)

           !write(*,*) '========================'
           !write(*,*) vkpt1
           !write(*,*) vkpt3
           !write(*,*) '------------------------'

           v3 = vkpt1-vkpt3
           dist3 = norm2(matmul(latt_cur%b,v3))

           ! find all the points in the FBZ that fall along this line segment
           ik2=0
           allocate(points_seg_dists(LIST_INIT))
           allocate(points_seg_coord(3,LIST_INIT))
           allocate(points_seg_ikibz(LIST_INIT))
           do ikx=1,kpoints%isnf(1,1)
              do iky=1,kpoints%isnf(2,2)
                 do ikz=1,kpoints%isnf(3,3)
                    ! k-point coordinates in the basis BK:
                    !V(1)=REAL(ikx,KIND=q)+kpoints%SHIFT(1)
                    !V(2)=REAL(iky,KIND=q)+kpoints%SHIFT(2)
                    !V(3)=REAL(ikz,KIND=q)+kpoints%SHIFT(3)
                    v=real(matmul(kpoints%ilsnf,[ikx,iky,ikz]),q)+kpoints%shift
                    ! Take only k-points lying within the Brillouin zone (if we do not use
                    ! the same lattice type for the k-lattice as for the reciprocal lattice
                    ! it might happen that we have to cut some points ...):
                    !vr=matmul(kpoints%b,v)
                    ! k-point coordinates in basis B:
                    !vkpt2=matmul(transpose(latt_cur%a),vr)

                    ! loops over periodic images
                    do i=1,ndisp
                    do j=1,ndisp
                    do k=1,ndisp
                    vkpt2=matmul(t,v+matmul(kpoints%ilsnf,[disp(i)*kpoints%isnf(1,1),disp(j)*kpoints%isnf(2,2),disp(k)*kpoints%isnf(3,3)]))

                    ! check if collinear
                    v1=vkpt2-vkpt1
                    dist1=norm2(matmul(latt_cur%b,v1))
                    v2=vkpt2-vkpt3
                    dist2=norm2(matmul(latt_cur%b,v2))
                    if (abs(dist1+dist2-dist3)>1e-6) cycle
                    ! distance from point to line
                    dist_pline = norm2(crossprod(matmul(latt_cur%b,vkpt2-vkpt1),matmul(latt_cur%b,vkpt2-vkpt3)))/&
                                 norm2(matmul(latt_cur%b,vkpt3-vkpt1))
                    if (dist_pline>1e-8) cycle
                    ! found point so add it to the list
                    ik2 = ik2 + 1
                    call kpoint_integers(kpoints,latt_cur,vkpt2,kint)
                    ikfbz = kint_kindex(kpoints,kint)
                    ikibz = kpoints%fbz2ibz(1,ikfbz)
                    ! if needed increase the size of the list
                    n = size(points_seg_ikibz)
                    if (size(points_seg_ikibz)<ik2) then
                        allocate(points_tmp_ikibz(n+LIST_STEP))
                        allocate(points_tmp_dists(n+LIST_STEP))
                        allocate(points_tmp_coord(3,n+LIST_STEP))
                        points_tmp_ikibz(:n) = points_seg_ikibz(:n)
                        points_tmp_dists(:n) = points_seg_dists(:n)
                        points_tmp_coord(:,:n) = points_seg_coord(:,:n)
                        call move_alloc(points_tmp_ikibz,points_seg_ikibz)
                        call move_alloc(points_tmp_dists,points_seg_dists)
                        call move_alloc(points_tmp_coord,points_seg_coord)
                    endif
                    ! compute actual distance (no pbc)
                    points_seg_coord(:,ik2) = vkpt2
                    points_seg_dists(ik2) = dist1
                    points_seg_ikibz(ik2) = ikibz
                    !write(*,*) ik2, ikibz, dist1, dist2, dist1+dist2, vkpt2
                    goto 100
                    enddo
                    enddo
                    enddo
                    100 continue
                 enddo !ikz
              enddo !iky
           enddo !ikx

           ! now sort the points by distance
           allocate(points_idx(ik2))
           do i=1,ik2
              points_idx(i) = i
           enddo
           call mrgrnk(points_seg_dists(:ik2),points_idx)
           allocate(points_tmp_dists(ik2))
           allocate(points_tmp_coord(3,ik2))
           allocate(points_tmp_ikibz(ik2))
           points_tmp_dists(:ik2) = points_seg_dists(points_idx)
           points_tmp_coord(:,:ik2) = points_seg_coord(:,points_idx)
           points_tmp_ikibz(:ik2) = points_seg_ikibz(points_idx)
           deallocate(points_seg_dists)
           deallocate(points_seg_coord)
           deallocate(points_seg_ikibz)
           deallocate(points_idx)

           !write(*,*) points_tmp_dists(:ik2)
           !write(*,*) points_tmp_ikibz(:ik2)

           ! keep unique points
           allocate(mask(ik2))
           mask = .false.
           mask(1:ik2-1) = points_tmp_ikibz(1:ik2-1)==points_tmp_ikibz(2:ik2)

           n = count(.not.mask)
           allocate(points_seg_dists(n))
           allocate(points_seg_ikibz(n))
           allocate(points_seg_coord(3,n))
           n = 0
           do i=1,ik2
              if (mask(i)) cycle
              n=n+1
              points_seg_dists(n)=points_tmp_dists(i)
              points_seg_coord(:,n)=points_tmp_coord(:,i)
              points_seg_ikibz(n)=points_tmp_ikibz(i)
           enddo
           deallocate(mask)
           deallocate(points_tmp_dists)
           deallocate(points_tmp_coord)
           deallocate(points_tmp_ikibz)

           !write(*,*) points_seg_dists(:n)
           !write(*,*) points_seg_ikibz(:n)

           ! add to the final list
           call move_alloc(points_dists,points_tmp_dists)
           call move_alloc(points_coord,points_tmp_coord)
           call move_alloc(points_ikibz,points_tmp_ikibz)
           nkpt_total = nkpt_total + n
           allocate(points_ikibz(nkpt_total))
           allocate(points_dists(nkpt_total))
           allocate(points_coord(3,nkpt_total))
           points_ikibz(:nkpt_total-n) = points_tmp_ikibz
           points_dists(:nkpt_total-n) = points_tmp_dists
           points_coord(:,:nkpt_total-n) = points_tmp_coord
           deallocate(points_tmp_dists)
           deallocate(points_tmp_coord)
           deallocate(points_tmp_ikibz)
           points_ikibz(nkpt_total-n+1:) = points_seg_ikibz
           points_dists(nkpt_total-n+1:) = points_seg_dists+dist_total
           points_coord(:,nkpt_total-n+1:) = points_seg_coord
           deallocate(points_seg_dists)
           deallocate(points_seg_ikibz)
           deallocate(points_seg_coord)
           dist_total = dist_total + dist3
        enddo ! loop over segments
        PROFILING_STOP('kpoints_along_path')

        contains
        function crossprod(a,b)
            integer :: ndim
            real(q) :: a(3),b(3)
            real(q) :: x,y,z
            real(q) :: crossprod(3)
            x=a(2)*b(3)-b(2)*a(3)
            y=a(3)*b(1)-b(3)*a(1)
            z=a(1)*b(2)-b(1)*a(2)
            crossprod=[x,y,z]
        end function crossprod
    end subroutine kpoints_along_path

    !> write the k-points along path to a file
    subroutine kpoints_along_path_write(kpoints_ikibz,kpoints_dists,kpoints_coord,iu6)
       real(q),intent(in) :: kpoints_dists(:)
       integer,intent(in) :: kpoints_ikibz(:)
       real(q),intent(in) :: kpoints_coord(:,:)
       integer,intent(in) :: iu6
       ! local variables
       integer :: i
       if (iu6>0) then
           write(iu6,'(a)')
           write(iu6,'(a)') 'Indexes of the k-points in the IBZ specified in the KPOINTS_DENSE file'
           write(iu6,'(a)') 'that lie along the path specified by KPOINTS_LINE'
           write(iu6,'(a)')
           write(iu6,'(a)') 'ikibz    distance'
           do i=1,size(kpoints_ikibz)
               write(iu6,'(i5,e12.4,3e12.4)') kpoints_ikibz(i), kpoints_dists(i), kpoints_coord(:,i)
           enddo
           write(iu6,*)
       endif
    end subroutine kpoints_along_path_write

    !> Check if kpoints struct contains information on a regular grid
    logical function kpoints_is_regular(self)
        type(kpoints_struct) :: self
        kpoints_is_regular = (self%mode==MonkhorstPack.or.&
                              self%mode==kSpacing.or.&
                              self%mode==RkLength.or.&
                              self%mode==GeneratingLattice)
    end function

!******************** SUBROUTINE HNFORM ********************************
!
!> Computes the Hermite Normal Form representation of the
!> integer matrix S, where:
!> ~~~
!>   B = BK S
!> ~~~
!> where the columns of B are the reciprocal lattive vectors and the
!> columns of BK the generating vectors of the k-point grid.
!>
!> @param [in]  B   Reciprocal lattice vectors (cartesian coordinates)
!> @param [in]  BK  Generating vectors of the k-point grid (cartesian coordinates)
!> @param [out] HNF Hermite Normal Form of S
!>
!> @note Basically this subroutine is not more than a wrapper
!> around HermiteNormalForm (written by Gus Hart).
!
!***********************************************************************

    SUBROUTINE HNFORM(B,BK,HNF,IERR)
      USE sym_prec
      USE tutor, ONLY: vtutor
      ! arguments
      REAL(q), INTENT(IN), DIMENSION(3,3) :: B, BK
      INTEGER, INTENT(OUT),DIMENSION(3,3) :: HNF
      INTEGER, INTENT(OUT) :: IERR
      ! local variables
      REAL(q), DIMENSION(3,3) :: AK,RS
      INTEGER, DIMENSION(3,3) :: S,HNB,L,SNF,R
      IERR = 0
      ! Express B in the basis BK
      CALL RECIPS(1._q,BK(1,1),BK(1,2),BK(1,3),AK(1,1),AK(1,2),AK(1,3))
      RS=MATMUL(TRANSPOSE(AK),B)
      ! RS should be an integer matrix
      IF (ANY(ABS(RS-NINT(RS))>TINY)) THEN
         IERR = 1
         RETURN
      ENDIF
      S=NINT(RS)
      ! and calculate the Hermite Normal Form of S
      CALL HermiteNormalForm(S,HNF,HNB)
    END SUBROUTINE HNFORM


!******************** SUBROUTINE SNFORM ********************************
!
!> Computes the Smith Normal Form representation of the
!> integer matrix S, where:
!> ~~~
!>   B = BK S
!> ~~~
!> where the columns of B are the reciprocal lattive vectors and the
!> columns of BK the generating vectors of the k-point grid.
!>
!> @param [in]  B   Reciprocal lattice vectors (cartesian coordinates)
!> @param [in]  BK  Generating vectors of the k-point grid (cartesian coordinates)
!> @param [out] SNF Smith Normal Form of S
!> @param [out] L   Left Transform
!>
!> @note Basically this subroutine is not more than a wrapper
!> around SmithNormalForm (written by Gus Hart).
!
!***********************************************************************

    SUBROUTINE SNFORM(B,BK,SNF,L)
      USE sym_prec
      USE tutor, ONLY: vtutor
      ! arguments
      REAL(q), INTENT(IN), DIMENSION(3,3) :: B, BK
      INTEGER, INTENT(OUT),DIMENSION(3,3) :: SNF,L
      ! local variables
      REAL(q), DIMENSION(3,3) :: AK,RS
      INTEGER, DIMENSION(3,3) :: S,R
      LOGICAL :: LERR
      ! Express B in the basis BK
      CALL RECIPS(1._q,BK(1,1),BK(1,2),BK(1,3),AK(1,1),AK(1,2),AK(1,3))
      RS=MATMUL(TRANSPOSE(AK),B)
      ! RS should be an integer matrix
      IF (ANY(ABS(RS-NINT(RS))>TINY)) THEN
         CALL vtutor%error("SNFORM: k-point generating vectors and reciprocal lattice are &
            &incommensurate.")
      ENDIF
      S=NINT(RS)
      ! and calculate the Smith Normal Form of S
      ! SNF = L S R
      CALL SmithNormalForm(S,L,SNF,R)
    END SUBROUTINE SNFORM

    !> Find the Hermite normal form of a given integer matrix.
    !> Similar to the SNF finder above but a little simpler.
    !> This routine is not very elegant, just brute force.
    !> Don't be disappointed.
    !> @param [in]  S The 3x3 integer matrix describing the relationship between two commensurate lattices.
    !> @param [out] H The resulting HNF matrix.
    !> @param [out] B The transformation matrix such that H = SB.
    SUBROUTINE HermiteNormalForm(S,H,B)
      use tutor, only: vtutor
      integer, intent(in) :: S(3,3)
      integer, intent(out), dimension(3,3) :: H, B

      !!<local name="tempcol">When two columns need to be swapped,
      !!tempcol mediates the swap.</local>
      !!<local name="check">Stores the HNF as a vector; only used as
      !!failsafe for meeting HNF criteria.</local>
      !!<local name="maxidx, minidx">The indices of max/min values in a
      !!column.</local>
      !!<local name="multiple">When reducing the columns to HNF form,
      !!the multiplicative relationship between two elements of two
      !!columns; reused.</local>
      !!<local name="minm">The smallest element in row 1 of the HNF
      !!being constructed.</local>

      integer :: i, minm, maxidx, minidx, multiple, j,  check(9), tempcol(3)

      if (determinant(S) == 0) &
         call vtutor%error("Singular matrix passed to HNF routine")
      B = 0; H = S ! H starts out as S, the input matrix
      forall(i=1:3); B(i,i) = 1; end forall ! B = identity

      do ! Keep doing column operations until all elements in row 1 are
         ! zero except the one on the diagonal.
         ! Divide the column with the smallest value into the largest
         do while (count(H(1,:)/=0) > 1) ! Keep going until only zeros beyond first element
            call get_minmax_indices(H(1,:),minidx,maxidx)
            minm = H(1,minidx)
            ! Subtract a multiple of the column containing the smallest element from
            ! the row containing the largest element
            multiple = H(1,maxidx)/minm ! Factor to multiply by
            H(:,maxidx) = H(:,maxidx) - multiple*H(:,minidx)
            B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
            if (any(matmul(S,B)/=H)) &
               call vtutor%error("COLS: Transformation matrices didn't work")
         enddo ! End of step row 1
         if (H(1,1) == 0) call swap_column(H,B,1) ! swap columns if (1,1) is zero
         if (H(1,1)<0)then; H(:,1) = -H(:,1); B(:,1) = -B(:,1) ! Change sign if (1,1) is negative
         endif
         if (count(H(1,:)/=0) > 1) &
             call vtutor%error("Didn't zero out the rest of the row")
         if (any(matmul(S,B)/=H)) &
             call vtutor%error("COLSWAP: Transformation matrices didn't work")
         exit
      enddo
      ! Now work on element H(2,3)
      do
         do while (H(2,3)/=0)
            if (H(2,2) == 0) then
               tempcol = H(:,2); H(:,2) = H(:,3); H(:,3) = tempcol
               tempcol = B(:,2); B(:,2) = B(:,3); B(:,3) = tempcol
               if (H(2,3) == 0) exit
            endif
            if (abs(H(2,3))<abs(H(2,2))) then; maxidx = 2; minidx = 3
            else; maxidx = 3; minidx = 2; endif
            multiple = H(2,maxidx)/H(2,minidx)
            H(:,maxidx) = H(:,maxidx) - multiple*H(:,minidx)
            B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
            if (any(matmul(S,B)/=H)) &
               call vtutor%error("COLS: Transformation matrices didn't work")
         enddo
         if (H(2,2) == 0) then
            tempcol = H(:,2); H(:,2) = H(:,3); H(:,3) = tempcol
         endif
         if (H(2,2)<0) then ! Change signs
         H(:,2) = -H(:,2); B(:,2) = -B(:,2); endif
         if (H(2,3)/=0) &
            call vtutor%error("Didn't zero out last element")

         if (any(matmul(S,B)/=H)) &
            call vtutor%error("COLSWAP: Transformation matrices didn't work")
         exit
      enddo
      if (H(3,3)<0) then ! Change signs
         H(:,3) = -H(:,3); B(:,3) = -B(:,3);
      endif
      check = reshape(H,(/9/))
      if (any(check((/4,7,8/))/=0)) &
         call vtutor%error("Not lower triangular")
      if (any(matmul(S,B)/=H)) &
         call vtutor%error("END PART1: Transformation matrices didn't work")

      ! Now that the matrix is in lower triangular form, make sure the lower off-diagonal
      ! elements are non-negative but less than the diagonal elements
      do while (H(2,2) <= H(2,1) .or. H(2,1)<0)
         if (H(2,2) <= H(2,1)) then; multiple = 1
         else; multiple = -1; endif
         H(:,1) = H(:,1) - multiple*H(:,2)
         B(:,1) = B(:,1) - multiple*B(:,2)
      enddo
      do j = 1,2
         do while (H(3,3) <= H(3,j) .or. H(3,j)<0)
            if (H(3,3) <= H(3,j)) then; multiple = 1
            else; multiple = -1; endif
            H(:,j) = H(:,j) - multiple*H(:,3)
            B(:,j) = B(:,j) - multiple*B(:,3)
         enddo
      enddo
      if (any(matmul(S,B)/=H)) &
         call vtutor%error("END: Transformation matrices didn't work")
      check = reshape(H,(/9/))
      if (any(check((/4,7,8/))/=0)) &
         call vtutor%error("Not lower triangular")
      if (any(check((/2,3,6/))<0)) &
         call vtutor%error("Negative elements in lower triangle")
      if (check(2) > check(5) .or. check(3) > check(9) .or. check(6) > check(9)) &
         call vtutor%error("Lower triangular elements bigger than diagonal")
    END SUBROUTINE HermiteNormalForm

    !> This routine takes an integer 3x3 matrix and computes its Smith Normal Form.
    !> M = A H B
    !> @param [in]  H Input matrix.
    !> @param [out] A Left Transform.
    !> @param [out] M Smith Normal Form matrix.
    !> @param [out] B Right Transform.
    SUBROUTINE SmithNormalForm(H,A,M,B)
      use tutor, only: vtutor
      integer, intent(in) :: H(3,3)
      integer, intent(out), dimension(3,3) :: M, A, B

      integer :: i, minm, maxidx, minidx, multiple, j, nondividx(2), check(9)
      logical :: Ldiv(2,2)
      integer :: itCnt

      if(determinant(H)<1) &
         call vtutor%error("SmithNormalForm routine failed because the input matrix had a negative determinant")
      A = 0; B = 0; M = H ! M starts out as H, the input matrix
      forall(i=1:3); A(i,i) = 1; B(i,i) = 1; end forall ! A & B = identity

      j=1
      itCnt = 0

      do ! Keep doing steps (1) and (2) until all elements in j-th row and column are
         ! zero except the one on the diagonal. When j == 1, if the (1,1) element doesn't
         ! divide every other non-zero element left in the matrix at the end, then add the
         ! offending row to row 1 and try again.
         !(1) Use row operations to zero out the column
         ! Divide the row with the smallest value into the largest

         itCnt = itCnt + 1
         if (itCnt>=100) &
            call vtutor%error("ERROR bad programming in SmithNormalForm")

         do while (count(M(:,j)/=0) > 1) ! Keep going until only 1 non-zero element in column j

            call get_minmax_indices(M(:,j),minidx,maxidx)
            minm = M(minidx,j)
            ! Subtract a multiple of the row containing the smallest element from
            ! the row containing the largest element
            multiple = M(maxidx,j)/minm
            M(maxidx,:) = M(maxidx,:) - multiple*M(minidx,:)
            A(maxidx,:) = A(maxidx,:) - multiple*A(minidx,:)
            if (any(matmul(matmul(A,H),B)/=M)) &
               call vtutor%error("ROW: Transformation matrices didn't work")
         enddo ! End of step 1

         if (M(j,j) == 0) then; call swap_row(A,M,j)
         endif
         if (any(matmul(matmul(A,H),B)/=M)) &
            call vtutor%error("ROWSWAP: Transformation matrices didn't work")

         !(2) Use column operations to zero out first row
         ! Divide the colum with the smallest value into the largest
         do while (count(M(j,:)/=0) > 1) ! Keep going until only 1 non-zero element in row 1
            !call printAMB(A,M,B)

            call get_minmax_indices(M(j,:),minidx,maxidx)
            minm = M(j,minidx)
            ! Subtract a multiple of the column containing the smallest element from
            ! the row containing the largest element
            multiple = M(j,maxidx)/minm ! Factor to multiply by
            M(:,maxidx) = M(:,maxidx) - multiple*M(:,minidx)
            B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
            if (any(matmul(matmul(A,H),B)/=M)) &
               call vtutor%error("COLS: Transformation matrices didn't work")
         enddo ! End of step 2

         if (M(j,j)<0) then ! Change signs
            M(:,j) = -M(:,j); B(:,j) = -B(:,j)
         elseif (M(j,j) == 0) then
            call swap_column(M,B,j)
         endif
         if (count(M(j,:)/=0) > 1 .or. count(M(:,j)/=0) > 1) cycle

         if (any(matmul(matmul(A,H),B)/=M)) &
            call vtutor%error("COLSWAP: Transformation matrices didn't work")

         Ldiv = mod(M(2:,2:),M(1,1)) == 0
         if (j==1 .and. any(Ldiv .eqv. .false.)) then! Add the offending row to row 1
            !      nondividx = maxloc(mod(M(2:,2:),M(1,1))) ! Find one of the elements that is not
            nondividx = maxloc(abs( mod(M(2:,2:),M(1,1)) )) ! Find one of the elements that is not
            M(1,:) = M(1,:) + M(nondividx(1)+1,:)    ! divided by the diagonal element, and
            A(1,:) = A(1,:) + A(nondividx(1)+1,:)    ! add the row it's in to row 1
            cycle ! Go back to the top of the outer do loop
         endif
         if (j==2) then
            if (mod(M(3,3),M(2,2))/=0) then
               M(2,:) = M(2,:) + M(3,:)
               A(2,:) = A(2,:) + A(3,:)
               cycle
            endif
         else
            j = 2;cycle;endif ! Start row/col 2
         ! Try again if the matrix still isn't diagonal
         if (j == 2 .and. (M(3,2)/=0 .or. M(2,3)/=0)) then; cycle; endif
         exit ! We should be done if we hit this point
      enddo
      if (M(3,3)<0) then ! Change signs
         M(:,3) = -M(:,3); B(:,3) = -B(:,3);
      endif
      if (any(matmul(matmul(A,H),B)/=M)) &
         call vtutor%error("END: Transformation matrices didn't work")
      check = reshape(M,(/9/))
      if (any(check((/2,3,4,6,7,8/))/=0)) &
         call vtutor%error("Not diagonal")
      if (mod(M(2,2),M(1,1))/=0 .or. mod(M(3,3),M(2,2))/=0) &
         call vtutor%error("SNF conditions not met")
    END SUBROUTINE SmithNormalForm

    !!<summary>Support routines for the Smith and Hermite normal form finders.</summary>
    !!<parameter name="A" regular="true"></parameter>
    !!<parameter name="M" regular="true"></parameter>
    !!<parameter name="k" regular="true"></parameter>
    SUBROUTINE swap_row(A,M,k) ! Swap rows of M (and A)
      integer, intent(inout) :: M(3,3), A(3,3)
      integer, intent(in) :: k
      integer :: tmpRow(3), maxidx(1)

      maxidx = maxloc(abs(M(k:,k)))+k-1  ! find index of the non-zero element in col k
      tmpRow = A(k,:); A(k,:) = A(maxidx(1),:); A(maxidx(1),:) = tmpRow
      tmpRow = M(k,:); M(k,:) = M(maxidx(1),:); M(maxidx(1),:) = tmpRow
    END SUBROUTINE swap_row

    !!<summary>Swaps the column 'k' with whichever column has the
    !!highest value (out of the columns to the right of 'k' in row
    !!'k'). The swap is performed in both matrices 'M' and
    !!'B'.</summary>
    !!<parameter name="M" regular="true">The matrix being
    !!transformed.</parameter>
    !!<parameter name="B" regular="true">Usually a matrix to keep track
    !!of the transformation steps on 'M'.</parameter>
    !!<parameter name="k" regular="true">The column to swap, as
    !!described in summary.</parameter>
    SUBROUTINE swap_column(M,B,k)
      integer, intent(inout) :: M(3,3), B(3,3)
      integer, intent(in) :: k
      integer :: tmpCol(3), maxidx(1)

      maxidx = maxloc(abs(M(k,k:)))+k-1 ! find index of the non-zero element in row k
      tmpCol = B(:,k); B(:,k) = B(:,maxidx(1)); B(:,maxidx(1)) = tmpCol
      tmpCol = M(:,k); M(:,k) = M(:,maxidx(1)); M(:,maxidx(1)) = tmpCol
    END SUBROUTINE swap_column

    !!<summary>Finds the indices corresponding the minimum and maximum
    !!values in an integer vector.</summary>
    !!<parameter name="invec" regular="true"></parameter>
    !!<parameter name="min" regular="true"></parameter>
    !!<parameter name="max" regular="true"></parameter>
    SUBROUTINE get_minmax_indices(invec,min,max)
      integer, intent(in) :: invec(3)
      integer, intent(out) :: min, max

      integer :: tmpmin(1), tmpmax(1), vec(3)
      vec = abs(invec)
      tmpmin = minloc(vec,vec>0)
      ! Search from the right for the max so it will be different from the minimum
      ! even if the min and max are the same value
      tmpmax = 4 - maxloc(vec(3:1:-1),vec(3:1:-1)>0)
      min = tmpmin(1)
      max = tmpmax(1)
    END SUBROUTINE get_minmax_indices

    !!<summary>This routine finds the determinant of 2x2 or 3x3 matrices</summary>
    !!<parameter name="a" regular="true"></parameter>
    FUNCTION determinant_real(a)
      use tutor, only: vtutor
      real(q) :: determinant_real,a(:,:)

      integer :: n
      n=size(a,dim=1)
      if (n==2) then
         determinant_real=-a(1,2)*a(2,1)+a(1,1)*a(2,2)
      else if (n==3) then
         determinant_real=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
              a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
              a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
      else
         call vtutor%error("DETERMINANT: Matrix dimension exceeds 3.")
      end if
    END FUNCTION determinant_real

    !!<summary>This routine finds the determinant of 2x2 or 3x3 matrices</summary>
    !!<parameter name="a" regular="true"></parameter>
    FUNCTION determinant_integer(a)
      use tutor, only: vtutor
      integer :: determinant_integer,a(:,:)
      integer :: n
      n=size(a,dim=1)
      if (n==2) then
         determinant_integer=-a(1,2)*a(2,1)+a(1,1)*a(2,2)
      else if (n==3) then
         determinant_integer=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
              a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
              a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
      else
         call vtutor%error("DETERMINANT: Matrix dimension exceeds 3.")
      end if
    END FUNCTION determinant_integer

    !!<summary>Given the matrix a, finds its inverse b</summary>
    !!<parameter name="a" regular="true"></parameter>
    !!<parameter name="b" regular="true"></parameter>
    !!<parameter name="err_" regular="true"></parameter>
    !!<parameter name="eps_" regular="true"></parameter>
    SUBROUTINE matrix_inverse(a, b, err_, eps_)
      real(q), intent(in):: a(3,3)
      real(q),intent(out):: b(3,3)

      logical, optional :: err_
      real(q), optional :: eps_

      real(q) :: c,avec(9)
      real(q) :: eps

      if(present(err_)) err_ = .false.
      if(present(eps_)) then; eps=eps_; else; eps=10d-14; endif

      c=a(1,3)*(-a(2,2)*a(3,1)+a(2,1)*a(3,2))+ &
           a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+ &
           a(1,1)*(-a(2,3)*a(3,2)+a(2,2)*a(3,3))
      avec=(/-a(2,3)*a(3,2)+a(2,2)*a(3,3),a(1,3)*a(3,2)-a(1,2)*a(3,3), &
           -a(1,3)*a(2,2)+a(1,2)*a(2,3),a(2,3)*a(3,1)-a(2,1)*a(3,3), &
           -a(1,3)*a(3,1)+a(1,1)*a(3,3),a(1,3)*a(2,1)-a(1,1)*a(2,3), &
           -a(2,2)*a(3,1)+a(2,1)*a(3,2),a(1,2)*a(3,1)-a(1,1)*a(3,2), &
           -a(1,2)*a(2,1)+a(1,1)*a(2,2)/)
      if(abs(c) < eps) then; if(present(err_)) err_ = .true.
      else; b=1/c*transpose(reshape(avec,(/3,3/))); endif
      ! transpose due to column major storage
    END SUBROUTINE matrix_inverse

END MODULE MKPOINTS
