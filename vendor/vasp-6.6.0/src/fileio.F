#include "symbol.inc"
!***********************************************************************
! RCS:  $Id: fileio.F,v 1.12 2003/06/27 13:22:18 kresse Exp kresse $
!
!> collection of subroutines to perform File IO
!>
!>  not yet really supported heavily
!>  but we want to switch all IO from the main program to this module
!
!***********************************************************************
      MODULE fileio
      USE prec
      IMPLICIT NONE

!
!> unit for the density matrix GAMMA (file GAMMA)
!
      INTEGER :: IU_GAMMA=77

      INTERFACE str
         MODULE PROCEDURE WRT_DENSITY_TO_STRING
      END INTERFACE

      CONTAINS

!****************** SUBROUTINE WRITE_GRID_QUANTITY *********************
!> @brief: Write any grid quantity to text file and/or hdf5 file.
!>
!> Standard restart files are treated as special cases:
!>      CHGCAR, TAUCAR, JDNCAR
!>      vaspwave.h5
!> This is why there are many (rarely used) optional arguments.
!>
!> To write any grid quantity as
!>     text to MYQNT                          and
!>     /results/myqnt Group to vaspout.h5
!> use, e.g.,
!> CALL WRITE_GRID_QUANTITY("MYQNT",WDES,T_INFO%POSION,T_INFO,LATT_CUR,GRIDC,CQNT)
!>
!> To write only to vaspout.h5, set LTEXTFILE=.FALSE. as optional
!> argument, use the group defined in vhdf5_base as KEYNAME and protect
!> the call with VASP_HDF5:
!> #ifdef VASP_HDF5
!>    CALL WRITE_GRID_QUANTITY(GRP_MYQNT,WDES,T_INFO%POSION,T_INFO,LATT_CUR,GRIDC,&
!>                               CQNT,LTEXTFILE=.FALSE.)
!> #endif
!***********************************************************************
      SUBROUTINE WRITE_GRID_QUANTITY(KEYNAME,WDES,POSION,T_INFO,LATT_CUR,GRIDC,CQNT,&
              SZNAM,LCHARG,LCHARGH5,&
              P,RHOLM,LOVERL,&
              LLONG,LFFT,LTEXTFILE &
              )
        USE mgrid, ONLY: grid_3d
        USE poscar, ONLY: latt, type_info, OUTPOS
        USE pseudo_struct_def, ONLY: potcar
        USE pawm, ONLY: WRT_RHO_PAW
        USE string
        USE main_mpi, ONLY: DIR_APP, DIR_LEN
        USE tutor, ONLY: vtutor
        USE wave, ONLY: wavedes
#ifdef VASP_HDF5
        USE vhdf5, ONLY: IH5OUTFILEID, IH5WAVEFILEID, GRP_RESULTS, GRP_CHARGE, GRP_TAU, &
             VASPWAVE_FOUND, VH5_WRITE_LATTICE_ION, VH5_WRITE_GRID_QUANTITY, VH5_WRITE_ATOMMOM, &
             VH5_WRITE_PAW_OCCUPANCIES, VH5_DELETE, VH5_EXISTS, VH5_ERROR
#endif
        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: KEYNAME     !< filename (and keyword to set default for h5 GROUP)
        TYPE(wavedes), INTENT(IN)    :: WDES        !< descriptor of the wavefunction (need: WDES%COMM_KINTER, WDES%NCDIJ)
        REAL(q), INTENT(IN)          :: POSION(:,:) !< fractional coordinates in unit cell, usually DYN%POSION, DYN%POSIOC, or T_INFO%POSION
        TYPE(type_info), INTENT(IN)  :: T_INFO      !< ion-related info (required for OUTPOS)
        TYPE(latt), INTENT(IN)       :: LATT_CUR    !< lattice parameter, use only LATT_CUR%SCALE,LATT_CUR%A
        TYPE(grid_3d), INTENT(IN)    :: GRIDC       !< grid for quantity
        COMPLEX(q), INTENT(IN)       :: CQNT(:,:)   !< real-valued grid quantity in complex array, e.g., CHTOT total charge density on fine grid
        ! --- optional
        CHARACTER(LEN=*),INTENT(IN), OPTIONAL  :: SZNAM     !< first comment line in text file, default: " "*40, usually INFO%SZNAM1 for INCAR SYSTEM tag (or T_INFO%SZNAM2 for POSCAR header)
        LOGICAL, INTENT(IN), OPTIONAL          :: LCHARG    !< toggle to write CHGCAR, TAUCAR, and JDNCAR, usually IO%LCHARG,
        LOGICAL, INTENT(IN), OPTIONAL          :: LCHARGH5  !< toggle to write densities to vaspwave.h5, usually IO%LCHARGH5,
        TYPE(potcar), INTENT(IN), OPTIONAL     :: P(:)      !< structure to support radial grid, only required for CHGCAR
        REAL(q), INTENT(IN), OPTIONAL          :: RHOLM(:,:)!< channel occupancies, (N_MIX_PAW,WDES%NCDIJ) (can actually be complex, uses UBOUND to hack it in), only required for CHGCAR
        LOGICAL, INTENT(IN), OPTIONAL          :: LOVERL    !< true if overlap matrix used, default=T, usually INFO%LOVERL, only required for CHGCAR
        LOGICAL, INTENT(IN), OPTIONAL          :: LLONG     !< T: double precision, F: single precision
        LOGICAL, INTENT(IN), OPTIONAL          :: LFFT      !< T: CQNT is in reciprocal space and an FFT needs to be done before writing to file (default: F, except for special cases like CHGCAR)
        LOGICAL, INTENT(IN), OPTIONAL          :: LTEXTFILE !< default:T, set .FALSE. to write only vaspout.h5 and no text file, ignored for restart information (use LH5 instead)
        ! --- local
        INTEGER    :: QNTDIM    ! no. of components of grid quantity, i.e., second dimension of CQNT, usually same as WDES%NCDIJ
        INTEGER    :: ID, NI
        INTEGER    :: IU=18, NODE_ME, IONODE
        LOGICAL    :: LSDYN=.FALSE.
        COMPLEX(q), ALLOCATABLE :: WORK(:,:)
        ! --- initialize defaults
        CHARACTER(40)         :: SZNAM_
        LOGICAL               :: LCHARG_
        LOGICAL               :: LCHARGH5_
        LOGICAL               :: LOVERL_
        LOGICAL               :: LLONG_
        LOGICAL               :: LTEXTFILE_
        LOGICAL               :: LFFT_
#ifdef VASP_HDF5
        CHARACTER(LEN=40)     :: GROUP_
#endif

        ! --- update defaults
        SZNAM_=REPEAT(" ",40) ; IF (PRESENT(SZNAM)) SZNAM_=SZNAM
        LCHARG_=.FALSE.  ; IF (PRESENT(LCHARG)) LCHARG_=LCHARG
        LCHARGH5_=.FALSE.; IF (PRESENT(LCHARGH5)) LCHARGH5_=LCHARGH5
        LOVERL_=.TRUE.   ; IF (PRESENT(LOVERL)) LOVERL_=LOVERL
        LLONG_=.TRUE.    ; IF (PRESENT(LLONG)) LLONG_=LLONG
        LTEXTFILE_=.TRUE.; IF (PRESENT(LTEXTFILE)) LTEXTFILE_=LTEXTFILE
        QNTDIM=SIZE(CQNT, 2)
        ! add cases where FFT should be performed by default here
        IF (ANY((/"CHGCAR","TAUCAR"/)==KEYNAME)) THEN
          LFFT_=.TRUE.
        ELSE
          LFFT_=.FALSE.
        ENDIF
        IF (PRESENT(LFFT)) LFFT_=LFFT

        ! checks for charge density
        IF (KEYNAME=="CHGCAR".AND.(.NOT.PRESENT(P).OR..NOT.PRESENT(RHOLM).OR..NOT.PRESENT(LOVERL))) THEN
          CALL vtutor%bug("Internal error in WRITE_GRID_QUANTITY:", __FILE__, __LINE__)
        ENDIF

        ! write text file
        NODE_ME = WDES%COMM%NODE_ME
        IONODE = WDES%COMM%IONODE
#ifdef MPI
        IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
        IF ((LCHARG_.AND.ANY((/"CHGCAR","TAUCAR","JDNCAR"/)==KEYNAME))&
            .OR.(LTEXTFILE_.AND..NOT.ANY((/"CHGCAR","TAUCAR","JDNCAR"/)==KEYNAME))) THEN
          io_begin
          OPEN(UNIT=IU,FILE=DIR_APP(1:DIR_LEN)//KEYNAME,STATUS='REPLACE')
          CALL OUTPOS(IU,.FALSE.,SZNAM_,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,POSION)
          io_end
          DO ID=1,QNTDIM
            IF (KEYNAME=="CHGCAR".AND.ID>1)  THEN
              do_io WRITE(IU,'(5E20.12)') (T_INFO%ATOMOM(NI),NI=1,T_INFO%NIONS)
            ENDIF
            CALL OUTCHG(GRIDC,IU,.TRUE.,CQNT(:,ID),LFFT_)
            IF (KEYNAME=="CHGCAR") CALL WRT_RHO_PAW(P,T_INFO,LOVERL_,RHOLM(:,ID),GRIDC%COMM,IU)
          ENDDO
          do_io CLOSE(IU)
        ENDIF
#ifdef MPI
        ENDIF
#endif

#ifdef VASP_HDF5
        ! group in h5, e.g. partial_charge , default based on: KEYNAME
        GROUP_=trim(lowercase(KEYNAME))
        ! map filenames to h5 groups
        IF (KEYNAME=="CHGCAR") GROUP_=GRP_CHARGE
        IF (KEYNAME=="TAUCAR") GROUP_=GRP_TAU

        ! write vaspwave.h5 IH5WAVEFILEID
        IF (LCHARGH5_.AND.ANY((/"CHGCAR","TAUCAR","JDNCAR"/)==KEYNAME)) THEN
          ! restart quantities should all be given for the same structure
          ! so we write the structure only once
          CALL VH5_WRITE_LATTICE_ION(IH5WAVEFILEID,SZNAM_,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,&
                  .FALSE.,POSION,GROUP="structure")

          IF (KEYNAME=="CHGCAR") THEN
            ! bring charge density to real space
            ALLOCATE(WORK,SOURCE=CQNT)
            DO ID=1,QNTDIM
              CALL FFT3D(WORK(:,ID),GRIDC,1)
            ENDDO
            ! delete CHARGE group if already present
            VH5_CHECK( VH5_DELETE(IH5WAVEFILEID,GROUP_) )
            ! write to file
            CALL VH5_WRITE_GRID_QUANTITY(IH5WAVEFILEID,GROUP_,GROUP_,GRIDC,WORK)
            CALL VH5_WRITE_ATOMMOM(IH5WAVEFILEID,GROUP_,T_INFO)
            DO ID=1,QNTDIM
              CALL VH5_WRITE_PAW_OCCUPANCIES(IH5WAVEFILEID,GROUP_,P,T_INFO,LOVERL_,RHOLM(:,ID),GRIDC%COMM,ID)
            ENDDO
          ENDIF

          IF (KEYNAME=="TAUCAR") THEN
            ! bring kinetic energy density to real space
            ALLOCATE(WORK,SOURCE=CQNT)
            DO ID=1,QNTDIM
              CALL FFT3D(WORK(:,ID),GRIDC,1)
            ENDDO
            ! delete tau group if already present
            VH5_CHECK( VH5_DELETE(IH5WAVEFILEID,GROUP_) )
            CALL VH5_WRITE_GRID_QUANTITY(IH5WAVEFILEID,GROUP_,"values",GRIDC,WORK)
          ENDIF
        ENDIF

        ! write vaspout.h5 IH5OUTFILEID
        IF (.NOT.ANY((/"CHGCAR","TAUCAR","JDNCAR"/)==KEYNAME)) THEN
          GROUP_ = GRP_RESULTS//"/"//GROUP_
          ! write structure for each grid quantity
          ! could be different choice for e.g. POSION
          CALL VH5_WRITE_LATTICE_ION(IH5OUTFILEID,SZNAM_,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,&
                      .FALSE.,POSION,GROUP=GROUP_,SUBGROUP="structure/position")
          IF (LFFT_) THEN
            ALLOCATE(WORK,SOURCE=CQNT)
            ! go to real space
            DO ID=1,QNTDIM
              CALL FFT3D(WORK(:,ID),GRIDC,1)
            ENDDO
            ! write to file
            CALL VH5_WRITE_GRID_QUANTITY(IH5OUTFILEID,GROUP_,"values",GRIDC,WORK)
            DEALLOCATE(WORK)
          ELSE
            ! directly write to file
            CALL VH5_WRITE_GRID_QUANTITY(IH5OUTFILEID,GROUP_,"values",GRIDC,CQNT)
          ENDIF
        ENDIF
#endif
      END SUBROUTINE WRITE_GRID_QUANTITY

!*************************SUBROUTINE OPENWAV ***************************
!
!>   open WAVECAR file
!>
!>   optionally an three character extension can be supplied allowing the
!>   open an alternative WAVECAR file
!>
!> - the unit 12 is always connected with the WAVECAR file
!> - for the mpi chain version (repeated images) the files
!>   XXX/WAVECAR are opended
!> - the record length is read from the WAVECAR file
!
!***********************************************************************

      SUBROUTINE OPENWAV(IO, COM, EXT)

      USE tutor, ONLY: vtutor
      USE base
      USE main_mpi
      USE mpimy
#ifdef VASP_HDF5
      USE vhdf5_base
#endif
      IMPLICIT NONE

      TYPE (in_struct) IO
      TYPE (communic) COM
      LOGICAL :: JUNK
      REAL(q) :: RDUM, RISPIN
      INTEGER :: IDUM
      INTEGER(qi8) :: WAVECAR_SIZE
      CHARACTER (LEN=3), OPTIONAL :: EXT


      JUNK=.TRUE.
      IF (PRESENT(EXT)) THEN
         INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,EXIST=junk,SIZE=WAVECAR_SIZE)
      ELSE
         INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',EXIST=junk,SIZE=WAVECAR_SIZE)
      ENDIF
      CALLMPI( M_bcast_l(COM, JUNK, 1 ))
      IO%LFOUND=JUNK
#ifdef VASP_HDF5
      IF (JUNK.AND.WAVECAR_SIZE>0) THEN
        WAVECAR_FOUND=.TRUE.
      ENDIF
      ! vaspwave always present
      IF(IO%LH5) THEN
         VASPWAVE_FOUND = .TRUE.
         IO%LFOUND = .TRUE.
      ENDIF
      IF(WAVECAR_FOUND .AND. IO%LH5) THEN
         CALL vtutor%warning("Both vaspwave.h5 and WAVECAR are found. vaspwave.h5 will be used.")
         WAVECAR_FOUND=.FALSE.
      ENDIF

#endif

! first reopen with assumed (wrong) record length ICMPLX
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%ICMPLX)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%ICMPLX)
      ENDIF
! the first record contains the record length, get it ...
      RDUM=0._q
      READ(12,REC=1,ERR=17421) RDUM,RISPIN ; IDUM=NINT(RDUM)
! in the worst case IDUM could contain completely useless data and useless is
! all <=0 or all >10000000 (since on output we use RECL=IO%ICMPLX*MAX(NRPLWV,6)
! or RECL=(NB_TOT+2)*ICMPLX more than ten millions sounds not very realistic)
      IF (IDUM<=0) IDUM=IO%ICMPLX  ! -> error reading WAVECAR
      GOTO 17422

17421 CONTINUE
      IDUM=IO%ICMPLX  ! -> error reading WAVECAR
17422 CONTINUE
      CLOSE(12)
      IO%IRECLW=IDUM
! reopen with correct record length (clumsy all that, I know ...)
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%IRECLW)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%IRECLW)
      ENDIF

      END SUBROUTINE OPENWAV

!*************************SUBROUTINE CLOSEWAV **************************
!
!>   close unit 12 (usually connected to the  WAVECAR file)
!
!***********************************************************************

      SUBROUTINE CLOSEWAV
        CLOSE(12)
      END SUBROUTINE CLOSEWAV

!*************************SUBROUTINE INWAV  ****************************
!
!>   read wavefunctions header from file
!>
!>   just parse the old cell shape old cutoff etc.
!>   this is a new version that uses only real quantities and allows
!>   the WAVECAR file to be exchanged between different machines
!>   more easily (IEEE standard)
!
!***********************************************************************

      SUBROUTINE INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI, ISTART, IU0)
      USE prec
      USE wave
      USE lattice

      TYPE (wavedes)  WDES
      TYPE (latt)     LATT_INI,LATT_CUR
      INTEGER NODE_ME, IONODE, ISTART, I, J, IU0, NKPTSF, NBANDF, NALLOC
      REAL(q) RKPTSF, RBANDF, ENMAXI
      LOGICAL LDIFF
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif

      READ(12,REC=2,ERR=100) RKPTSF,RBANDF,ENMAXI, &
     &                       ((LATT_INI%A(I,J),I=1,3),J=1,3)

      IF (IU0 >= 0) WRITE(IU0,*) 'found WAVECAR, reading the header'

      NKPTSF=NINT(RKPTSF)
      NBANDF=NINT(RBANDF)

      CALL LATTIC(LATT_INI)

      IF (ISTART==2 .AND. ENMAXI /= WDES%ENMAX) THEN
        CALL vtutor%error("ERROR: ENMAX changed please set ISTART to 1")
      ENDIF

      IF (.NOT. (NBANDF==WDES%NB_TOT .OR. (WDES%NRSPINORS==2 .AND.  NBANDF*WDES%NRSPINORS==WDES%NB_TOT))) THEN
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'number of bands has changed, file:',NBANDF,' present:',WDES%NB_TOT
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'trying to continue reading WAVECAR, but it might fail'
      ENDIF
      IF (NKPTSF/=WDES%NKPTS) THEN
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'number of k-points has changed, file:',NKPTSF,' present:',WDES%NKPTS
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'trying to continue reading WAVECAR, but it might fail'
      ENDIF
      IF (ISTART ==1) THEN

      LDIFF=.FALSE.
      DO I=1,3
        DO J=1,3
          IF (ABS(LATT_INI%A(I,J)-LATT_CUR%A(I,J)) > 1E-4) LDIFF=.TRUE.
        ENDDO
      ENDDO
      IF (ENMAXI /= WDES%ENMAX) LDIFF=.TRUE.
      IF (LDIFF) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WAVECAR: different cutoff or change in lattice found'
      ENDIF

      ENDIF

      RETURN
!
! WAVECAR does not exist or can not be read at all
!
  100 CONTINUE
      LATT_INI%A=LATT_CUR%A
      CALL LATTIC(LATT_INI)

      IF (ISTART==3) THEN
        CALL vtutor%error("ERROR: can't restart (ISTART=3) with wavefunctions on file")
      ELSE
        ISTART=0
      ENDIF
      RETURN

      END SUBROUTINE


!*************************SUBROUTINE INWAV_FAST ************************
!
!>   read wavefunctions from file
!>
!>   the routine can read WAVECAR files generated for different
!>   lattice vectors and cutoff (both in the serial and parallel case)
!>
!>   it is also possible to restart a spin polarised calculation
!>   from a non spin polarised wavefunction file or
!>   to restart a non collinear calculation from a collinear calculation
!
!***********************************************************************

      SUBROUTINE INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, ISTART, EFERMI)
      USE prec
      USE base
      USE wave_high
      USE mgrid
      USE lattice
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      TYPE (wavespin) W
      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR,LATT_INI
      REAL (q) :: EFERMI

      INTEGER NODE_ME, IONODE, IU0
      INTEGER ISPINREAD, IREC, ISP, IFAIL, I, J, K, IRECL, ISTART, IS
      INTEGER NALLOC, MALLOC, NBANDF, NKPTSF, NPL, NPLREAD, NCOMM
      REAL(q) RDUM, RKPTSF, RBANDF, RISPIN, RTAG, RNPL, ENMAXF

! local work arrays
      REAL(q) VKPT(3)
      INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
      COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
      LOGICAL ::  SINGLE_PREC, SPAN_RECORDS
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()
      REAL(q), ALLOCATABLE :: INBUF(:)
      CHARACTER(LEN=:), ALLOCATABLE :: MSG

      IU0=IO%IU0

      NALLOC = MAXVAL(WDES%NPLWKP_TOT)

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif
!
! parse the header
!
        RTAG=0
        READ(12,REC=1,ERR=200) RDUM,RISPIN,RTAG
        IF (RTAG==53300) THEN
           IF (IU0>=0) WRITE(IU0,*) "VASP.5.3 WAVECAR encountered"
           SINGLE_PREC=.TRUE.
           SPAN_RECORDS=.TRUE.
        ELSE IF (RTAG==53310) THEN
           IF (IU0>=0) WRITE(IU0,*) "VASP.5.3 double precision WAVECAR encountered"
           SINGLE_PREC=.FALSE.
           SPAN_RECORDS=.TRUE.
        ELSE IF (RTAG==45200) THEN
           SINGLE_PREC=.TRUE.
           SPAN_RECORDS=.FALSE.
        ELSE
           IF (IU0>=0) WRITE(IU0,*) "double precision WAVECAR encountered, converting it"
           SINGLE_PREC=.FALSE.
           SPAN_RECORDS=.FALSE.
        ENDIF

        ISPINREAD=NINT(RISPIN)
        READ(12,REC=2,ERR=200) RKPTSF,RBANDF,ENMAXF, &
     &                         ((LATT_INI%A(I,J),I=1,3),J=1,3),EFERMI

        IREC=2
        NKPTSF=NINT(RKPTSF)
        NBANDF=NINT(RBANDF)

        CALL LATTIC(LATT_INI)
!=======================================================================
! read WAVECAR file, number of bands agree
!=======================================================================
      IF (.NOT.(WDES%NRSPINORS==2 .AND. NBANDF*WDES%NRSPINORS == WDES%NB_TOT)) THEN

        spin:    DO ISP=1, MIN(WDES%ISPIN, ISPINREAD)
        kpoints: DO K=1,WDES%NKPTS

        IF ( K> NKPTSF ) THEN
           ! continue reading last k-point 
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        CALL SETWDES(WDES,WDES1,K)
#ifdef gammareal
        CALL GEN_INDEX_GAMMA(WDES1, LATT_CUR%B, WDES%ENMAX, INDEX, LCONJG)
#endif

        io_begin
        IFAIL=0
        IF (.NOT. SPAN_RECORDS) THEN
           IREC=IREC+1
           READ(12,REC=IREC,ERR=230) RNPL,VKPT, &
                (W%CELTOT(J,K,ISP),W%FERTOT(J,K,ISP),J=1,MIN(WDES%NB_TOT,NBANDF))
        ELSE
           ALLOCATE(INBUF(4+NBANDF*3))
           CALL READ_TO_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, INBUF, IFAIL)
           RNPL   =INBUF(1)
           VKPT(1)=INBUF(2)
           VKPT(2)=INBUF(3)
           VKPT(3)=INBUF(4)
           DO J=1,MIN(WDES%NB_TOT,NBANDF)
              W%CELTOT(J,K,ISP)=CMPLX(INBUF(5+(J-1)*3),INBUF(6+(J-1)*3),q)
              W%FERTOT(J,K,ISP)=INBUF(7+(J-1)*3)
           ENDDO
           DEALLOCATE(INBUF)

        ENDIF
        NPLREAD=NINT(RNPL)
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 230

        CALLMPI( M_bcast_i( WDES%COMM, IREC, 1))
        CALLMPI( M_bcast_i( WDES%COMM, NPLREAD, 1))
        NPL=WDES%NPLWKP_TOT(K)

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(MALLOC),CW2(MALLOC),CRD(MALLOC),IND(MALLOC),INDI(MALLOC))

        io_begin
        ! create an index to allow for change of cutoff or cell size
        IF (ISTART==2) THEN
           ! gK: fix ISTART=2 in INWAV_FAST
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
              WDES%ENMAX, ENMAXF, NPL/WDES%NRSPINORS, NPLREAD/WDES%NRSPINORS, IND, INDI, MALLOC, IFAIL )
        ELSE
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
              WDES%ENMAX, ENMAXF, NPL/WDES%NRSPINORS, NPLREAD/WDES%NRSPINORS, IND, INDI, MALLOC, IFAIL )
        ENDIF
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 220

        band: DO J=1,NBANDF
          IREC=IREC+1
          IF (J>WDES%NB_TOT) CYCLE
          io_begin
            IF (SINGLE_PREC) THEN
               READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
               CW2(1:NPLREAD)=CRD(1:NPLREAD)
            ELSE
               READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
            ENDIF

            CW1=0
            DO IS=1,WDES%NRSPINORS
               ! store the wave function coefficients according to new cutoff
               CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                  CW1((IS-1)*NPL/WDES%NRSPINORS+1), CW2((IS-1)*NPLREAD/WDES%NRSPINORS+1))
            ENDDO

          io_end
#ifdef MPI
          IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
             CALLMPI( M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1)) )
          END IF
#endif
#ifdef gammareal
          ! gamma only go from seriel layout to parallel layout if required
          CALL DIS_PW_BAND_GAMMA(WDES1, J, CW1, W%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
          CALL DIS_PW_BAND(WDES1, J, CW1, W%CW(1,1,K,ISP))
#endif
        ENDDO band

#ifdef gammareal
        CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
        DEALLOCATE(CW1,CW2,CRD,IND,INDI)

        ENDDO kpoints

        IF (NKPTSF > WDES%NKPTS) THEN
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(WDES%NB_TOT+1)
        ENDIF

        ! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALLMPI( M_bcast_z(WDES%COMM, W%CELTOT(1,1,ISP),NCOMM ))
        CALLMPI( M_bcast_d(WDES%COMM, W%FERTOT(1,1,ISP),NCOMM ))

        ENDDO spin

        IF (ISPINREAD > WDES%ISPIN .AND. IU0>=0 ) THEN
           WRITE(IU0,*) 'down-spin wavefunctions not read'
        ENDIF

        IF (NBANDF<WDES%NB_TOT) THEN
           IF (IU0>=0) WRITE(IU0,*) 'random initialization beyond band ',NBANDF
           CALL WFINIT(WDES, W, 1E10_q, NBANDF+1) ! ENINI=1E10 not cutoff restriction
        ENDIF

        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read successfully'

        IF (WDES%ISPIN<=ISPINREAD) THEN
           RETURN
        ENDIF
!
!  spin down is missing
!
        IF (IU0>=0) &
        WRITE(IU0,*) 'No down-spin wavefunctions found', &
     &             ' --> setting down-spin equal up-spin ...'

        DO K=1,WDES%NKPTS
          W%CELTOT(1:WDES%NB_TOT,K,2)=W%CELTOT(1:WDES%NB_TOT,K,1)
          W%FERTOT(1:WDES%NB_TOT,K,2)=W%FERTOT(1:WDES%NB_TOT,K,1)
          NPL=WDES%NPLWKP(K)
          W%CW(1:NPL,1:WDES%NBANDS,K,2)=W%CW(1:NPL,1:WDES%NBANDS,K,1)
        ENDDO


        RETURN
!=======================================================================
! read collinear WAVECAR file for a non collinear run
!=======================================================================
      ELSE
        W%CW=0

        spin2:    DO ISP=1,ISPINREAD
        IF (IU0>=0.AND. ISP==1) &
        WRITE(IU0,*) 'reading wavefunctions of collinear run, up'
        IF (IU0>=0.AND. ISP==2) &
        WRITE(IU0,*) 'reading wavefunctions of collinear run, down'
        kpoints2: DO K=1,NKPTSF

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        CALL SETWDES(WDES,WDES1,K)

        io_begin
        IFAIL=0
        IF (.NOT. SPAN_RECORDS) THEN
           IREC=IREC+1
           READ(12,REC=IREC,ERR=230) RNPL, VKPT(1:3), &
                       (W%CELTOT(J+NBANDF*(ISP-1),K,1),W%FERTOT(J+NBANDF*(ISP-1),K,1),J=1,NBANDF)
        ELSE
           ALLOCATE(INBUF(4+NBANDF*3))

           CALL READ_TO_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, INBUF, IFAIL)
           RNPL   =INBUF(1)
           VKPT(1)=INBUF(2)
           VKPT(2)=INBUF(3)
           VKPT(3)=INBUF(4)
           DO J=1,NBANDF
              W%CELTOT(J+NBANDF*(ISP-1),K,1)=CMPLX(INBUF(5+(J-1)*3),INBUF(6+(J-1)*3),q)
              W%FERTOT(J+NBANDF*(ISP-1),K,1)=INBUF(7+(J-1)*3)
           ENDDO
           DEALLOCATE(INBUF)

        ENDIF
        NPLREAD=NINT(RNPL)
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 230

        CALLMPI( M_bcast_i( WDES%COMM, IREC, 1))
        CALLMPI( M_bcast_i( WDES%COMM, NPLREAD, 1))
        NPL=WDES%NPLWKP_TOT(K)/2

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(2*MALLOC),CW2(2*MALLOC),CRD(2*MALLOC),IND(MALLOC),INDI(MALLOC))

        io_begin
        IF (ISTART==2) THEN
           ! gK: fix ISTART=2 in INWAV_FAST
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
        ELSE
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
        ENDIF
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 220

        IF (ISP==1) THEN
           DO J=1,NBANDF
              IREC=IREC+1
              io_begin
                IF (SINGLE_PREC) THEN
                   READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                   CW2(1:NPLREAD)=CRD(1:NPLREAD)
                ELSE
                   READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
                ENDIF


                CW1=0
                CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1, CW2)
              io_end

#ifdef MPI
              IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
                 CALLMPI( M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1)) )
              END IF
#endif
              CALL DIS_PW_BAND(WDES1, J, CW1, W%CW(1,1,K,1))

              ! copy immediately to second panel
              W%CELTOT(J+NBANDF,K,1)=W%CELTOT(J,K,1)
              W%FERTOT(J+NBANDF,K,1)=W%FERTOT(J,K,1)

              CW1(NPL+1:2*NPL)=CW1(1:NPL)
              CW1(1:NPL)=0

              CALL DIS_PW_BAND(WDES1, J+NBANDF, CW1, W%CW(1,1,K,1))
           ENDDO
        ELSE
           DO J=1,NBANDF
              IREC=IREC+1
              io_begin
                IF (SINGLE_PREC) THEN
                   READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                   CW2(1:NPLREAD)=CRD(1:NPLREAD)
                ELSE
                   READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
                ENDIF
                CW1=0
                CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1(NPL+1), CW2)
              io_end
#ifdef MPI
              IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
                 CALLMPI( M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1)) )
              END IF
#endif
              CALL DIS_PW_BAND(WDES1, J+NBANDF, CW1, W%CW(1,1,K,1))
           ENDDO

        ENDIF

        DEALLOCATE(CW1,CW2,CRD,IND,INDI)

        ENDDO kpoints2

        IF (NKPTSF > WDES%NKPTS) THEN
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(NBANDF+1)
        ENDIF

        ENDDO spin2

        ! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALLMPI( M_bcast_z(WDES%COMM, W%CELTOT(1,1,1),NCOMM ))
        CALLMPI( M_bcast_d(WDES%COMM, W%FERTOT(1,1,1),NCOMM ))

        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read successfully'
        RETURN
      ENDIF

      MSG = 'ERROR: while reading WAVECAR, file is incompatible '
      IF (WDES%ENMAX /= ENMAXF) &
         MSG = MSG // 'the energy cutoff has changed (new,old) ' // str(WDES%ENMAX) // ' ' // str(ENMAXF)
      IF (NBANDF /= WDES%NB_TOT) &
         MSG = MSG // 'the number of bands has changed (new,old) ' // str(WDES%NB_TOT) // ' ' // str(NBANDF)
      IF (NKPTSF /= WDES%NKPTS) &
         MSG = MSG // 'the number of k-points has changed (new,old) ' //      str(WDES%NKPTS) // ' ' // str(NKPTSF)
      CALL vtutor%error(MSG)
!=======================================================================
! can not do anything with WAVECAR
! hard stop, pull all breaks
!=======================================================================
  200 CONTINUE

      CALL vtutor%error("ERROR: while reading WAVECAR, header is corrupt")

  220 CONTINUE
#if defined( wNGZhalf ) || defined( wNGXhalf) || defined( NGZhalf ) || defined( NGXhalf) 
      CALL vtutor%error("ERROR: while reading WAVECAR, plane wave coefficients changed " // str(NPL) &
         // " " // str(NPLREAD))
#else
      CALL vtutor%error("ERROR: while reading WAVECAR, plane wave coefficients changed.& 
              & If you intend to read the wavefunction from a vasp_std job,&
              & double NBANDS in INCAR to continue.")
#endif

  230 CONTINUE

      CALL vtutor%error("ERROR: while reading eigenvalues from WAVECAR " // str(K) // " " // str(ISP))

  240 CONTINUE

      CALL vtutor%error("ERROR: while reading plane wave coef. from WAVECAR " // str(K) // " " // &
         str(ISP) // " " // str(J))

      END SUBROUTINE

!*************************SUBROUTINE INWAV_ALTERNATIVE *****************
!
!>   read wavefunctions from alternative file
!>
!>   read alternative WAVECAR file
!
!***********************************************************************

      SUBROUTINE INWAV_ALTERNATIVE(IO, WDES, W, GRID, LATT_CUR, LREAD, EXT )

      USE base
      USE wave_high
      USE lattice
      USE tutor, ONLY: vtutor, isAlert, AuxWAVECAR

      IMPLICIT NONE

      TYPE (in_struct)  IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR
      LOGICAL         LREAD
      CHARACTER (LEN=3)  EXT
    ! local
      TYPE (latt)     LATT_INI
      INTEGER         ISTART
      REAL(q)         ENMAXI
      REAL (q)        EFERMI

      LREAD = .FALSE.

      ISTART=1
      CALL OPENWAV(IO, WDES%COMM, EXT)
      CALL INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI, ISTART, IO%IU0)

      IF (ISTART/=0) THEN
         CALL ALLOCW( WDES, W)
         CALL INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, ISTART, EFERMI)
         LREAD = .TRUE.
         CALL vtutor%write(isAlert, AuxWAVECAR)
      ENDIF

      CALL CLOSEWAV

      END SUBROUTINE INWAV_ALTERNATIVE

!*************************SUBROUTINE OUTWAV    ************************
!
!>   write wavefunctions to file
!>
!>   the first version writes a VASP 4 file
!>   the second version writes the newer VASP 5 file format
!>   which requires less storage if many bands are calculated
!
!***********************************************************************

      SUBROUTINE OUTWAV_4(IO, WDES, W, LATT_INI, EXT )
      USE base
      USE wave
      USE lattice
      USE main_mpi

      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (latt)     LATT_INI
      CHARACTER (LEN=3), OPTIONAL :: EXT
! local work arrays
      TYPE (wavedes1) WDES1
      COMPLEX(q),ALLOCATABLE :: CW(:),EIG(:)
#ifdef WAVECAR_double
      COMPLEX(q), ALLOCATABLE :: CRD(:)
#else
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
#endif
    ! local
      INTEGER NODE_ME, IONODE, IU0, IREC, IRECLW
      INTEGER NPL_TOT, NPL, I, J, K, ISP
      REAL(q) RDUM, RISPIN, RTAG, RNB_TOT, RNKPTS, RNPL
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()


      NODE_ME=0
      IONODE=0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif

!
! report to  stdout
!
      IF (IO%IU6>=0) WRITE(IO%IU6,*)'writing wavefunctions'
      IF (IO%IU0>=0) WRITE(IO%IU0,*)'writing wavefunctions'
!
! open file
!
      IU0   =IO%IU0
      NPL_TOT = MAXVAL(WDES%NPLWKP_TOT)
#ifdef WAVECAR_double
      IO%IRECLW=MAX(MAX(NPL_TOT,6),((WDES%NB_TOT*3+1)/2+2))*IO%ICMPLX
#else
      IO%IRECLW=MAX(MAX((NPL_TOT+1)/2,6),((WDES%NB_TOT*3+1)/2+2))*IO%ICMPLX
#endif
      IRECLW=IO%IRECLW

    ! open file
      io_begin
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ENDIF
      io_end

      ALLOCATE(CW(NPL_TOT),CRD(NPL_TOT),EIG(WDES%NB_TOT))

      RISPIN=WDES%ISPIN
      RDUM  =IRECLW

#ifdef WAVECAR_double
      RTAG  =45210
#else
      RTAG  =45200
#endif
      do_io WRITE(12,REC=1) RDUM,RISPIN,RTAG

      IREC=2
      io_begin
! in order to increase exchangeability of WAVECAR files across IEEE platforms
! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
      RNKPTS =WDES%NKPTS
      RNB_TOT=WDES%NB_TOT
      WRITE(12,REC=2) RNKPTS,RNB_TOT,WDES%ENMAX,((LATT_INI%A(I,J),I=1,3),J=1,3)
      io_end
!
! loop over spin, kpoints, bands
!
      spin: DO ISP=1,WDES%ISPIN
      kpoints: DO K=1,WDES%NKPTS

#ifdef MPI
       IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
          IF (WDES%COMM_KINTER%NODE_ME.NE.1) THEN
             IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                CYCLE
             ELSE
                CALLMPI( M_send_z(WDES%COMM_KINTER, 1, W%CW(1,1,K,ISP),SIZE(W%CW,1)*SIZE(W%CW,2)) )
             END IF
          ELSE
             IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                CALLMPI( M_recv_z(WDES%COMM_KINTER, MOD(K-1,WDES%COMM_KINTER%NCPU)+1, W%CW(1,1,K,ISP),SIZE(W%CW,1)*SIZE(W%CW,2)) )
             ENDIF
          ENDIF
       END IF
#endif
       CALL SETWDES(WDES,WDES1,K)
#ifdef gammareal
       CALL GEN_INDEX_GAMMA(WDES1, LATT_INI%B, WDES%ENMAX, INDEX, LCONJG)
#endif
       NPL=WDES%NPLWKP_TOT(K)
       RNPL=NPL
       IREC=IREC+1
! write number of plane waves, k-point coordinates and all eigenvalues and
! occupation numbers for current k-point K
       io_begin
! write eigenvalues in real format
         DO J=1,WDES%NB_TOT
           EIG(J)=REAL(W%CELTOT(J,K,ISP),KIND=q)
         ENDDO
         WRITE(12,REC=IREC) RNPL,WDES%VKPT(1,K),WDES%VKPT(2,K), &
                       WDES%VKPT(3,K),(EIG(J),W%FERTOT(J,K,ISP),J=1,WDES%NB_TOT)
       io_end

       DO J=1,WDES%NB_TOT
#ifdef gammareal
         CALL MRG_PW_BAND_GAMMA(WDES1, J, CW, W%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
         CALL MRG_PW_BAND(WDES1, J, CW, W%CW(1,1,K,ISP))
#endif

         IREC=IREC+1
         do_io CRD(1:NPL)=CW(1:NPL)
         do_io WRITE(12,REC=IREC) (CRD(I),I=1,NPL)
       ENDDO

#ifdef gammareal
       CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
      ENDDO kpoints
      ENDDO spin

      DEALLOCATE(CW,CRD,EIG)

      do_io CALL CLOSEWAV

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE OUTWAV    ************************
!
!>   write wavefunctions to file
!>
!>   this version uses smaller record lenght
!>   and if required eigenvalues/ occupancies will span
!>   several records
!
!***********************************************************************

      SUBROUTINE OUTWAV(IO, WDES, W, LATT_INI, EFERMI, EXT, NBANDS_DUMP)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      use tutor, only: VTUTOR

      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (latt)     LATT_INI
      CHARACTER (LEN=3), OPTIONAL :: EXT
      INTEGER,INTENT(IN),  OPTIONAL :: NBANDS_DUMP
      REAL (q) EFERMI
! local work arrays
      TYPE (wavedes1) WDES1
      COMPLEX(q),ALLOCATABLE :: CW(:)
#ifdef WAVECAR_double
      COMPLEX(q), ALLOCATABLE :: CRD(:)
#else
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
#endif
    ! local
      INTEGER IU0, NODE_ME, IONODE, IREC, IRECLW_OLD
      INTEGER NPL, NPL_TOT, I, J, K, ISP
      REAL(q) RDUM, RISPIN, RNKPTS, RNB_TOT, RTAG, RNPL
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()
      REAL(q), ALLOCATABLE :: OUTBUF(:)
      INTEGER  NBANDS_OLD 

      NODE_ME=0
      IONODE=0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif

      NBANDS_OLD = WDES%NB_TOT 
      IF ( PRESENT( NBANDS_DUMP ) ) THEN
         IF( NBANDS_DUMP > 1 ) THEN
            WDES%NB_TOT = NBANDS_DUMP
         ELSE
            CALL VTUTOR%BUG(" silly number of bands passed "//STR(NBANDS_DUMP),__FILE__, __LINE__ )
         ENDIF
      ENDIF
!
! report to  stdout
!
      IF (IO%IU6>=0) WRITE(IO%IU6,*)'writing wavefunctions'
      IF (IO%IU0>=0) WRITE(IO%IU0,*)'writing wavefunctions'
!
! open file
!
      IU0   =IO%IU0
      NPL_TOT = MAXVAL(WDES%NPLWKP_TOT)
#ifdef WAVECAR_double
      IO%IRECLW=MAX(NPL_TOT,7)*IO%ICMPLX
      IRECLW_OLD=MAX(MAX(NPL_TOT,6),((WDES%NB_TOT*3+1)/2+2))*IO%ICMPLX
#else
      IO%IRECLW=MAX((NPL_TOT+1)/2,7)*IO%ICMPLX
      IRECLW_OLD=MAX(MAX((NPL_TOT+1)/2,6),((WDES%NB_TOT*3+1)/2+2))*IO%ICMPLX
#endif

    ! open file
      io_begin
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ENDIF
      io_end

      ALLOCATE(CW(NPL_TOT),CRD(NPL_TOT))

      RISPIN=WDES%ISPIN
      RDUM  =IO%IRECLW

#ifdef WAVECAR_double
      RTAG  =53310
#else
      RTAG  =53300
#endif
      ! if the record lenght is identical to previous vasp version we can claim
      ! to write the old file format
      IF (IRECLW_OLD==IO%IRECLW) THEN
#ifdef WAVECAR_double
         RTAG  =45210
#else
         RTAG  =45200
#endif
      ENDIF

      do_io WRITE(12,REC=1) RDUM,RISPIN,RTAG

      IREC=2
      io_begin
! in order to increase exchangeability of WAVECAR files across IEEE platforms
! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
      RNKPTS =WDES%NKPTS
      RNB_TOT=WDES%NB_TOT
      WRITE(12,REC=2) RNKPTS,RNB_TOT,WDES%ENMAX,((LATT_INI%A(I,J),I=1,3),J=1,3),EFERMI
      io_end
!
! loop over spin, kpoints, bands
!
      spin: DO ISP=1,WDES%ISPIN
      kpoints: DO K=1,WDES%NKPTS

#ifdef MPI
       IF (WDES%COMM_KINTER%NCPU > 1) &
          CALL M_bcast_z8_from(WDES%COMM_KINTER, W%CW(1,1,K,ISP),SIZE(W%CW,1,KIND=qi8)*SIZE(W%CW,2,KIND=qi8),MOD(K-1,WDES%COMM_KINTER%NCPU)+1)
#endif
       CALL SETWDES(WDES,WDES1,K)
#ifdef gammareal
       CALL GEN_INDEX_GAMMA(WDES1, LATT_INI%B, WDES%ENMAX, INDEX, LCONJG)
#endif
       NPL=WDES%NPLWKP_TOT(K)
       RNPL=NPL
! write number of plane waves, k-point coordinates and all eigenvalues and
! occupation numbers for current k-point K
! write eigenvalues in real format
       io_begin
       ALLOCATE(OUTBUF(4+WDES%NB_TOT*3))

       OUTBUF(1)=RNPL
       OUTBUF(2)=WDES%VKPT(1,K)
       OUTBUF(3)=WDES%VKPT(2,K)
       OUTBUF(4)=WDES%VKPT(3,K)
       DO I=1,WDES%NB_TOT
          OUTBUF(5+(I-1)*3) =REAL(W%CELTOT(I,K,ISP),q)
          OUTBUF(6+(I-1)*3) =AIMAG(W%CELTOT(I,K,ISP))
          OUTBUF(7+(I-1)*3) =W%FERTOT(I,K,ISP)
       ENDDO
       CALL WRITE_FROM_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, OUTBUF)
       DEALLOCATE(OUTBUF)
       io_end

       DO J=1,WDES%NB_TOT
#ifdef gammareal
         CALL MRG_PW_BAND_GAMMA(WDES1, J, CW, W%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
         CALL MRG_PW_BAND(WDES1, J, CW, W%CW(1,1,K,ISP))
#endif

         IREC=IREC+1
         do_io CRD(1:NPL)=CW(1:NPL)
         do_io WRITE(12,REC=IREC) (CRD(I),I=1,NPL)
       ENDDO

#ifdef gammareal
       CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
      ENDDO kpoints
      ENDDO spin

      DEALLOCATE(CW,CRD)

      do_io CALL CLOSEWAV

      ! restore bands to original number 
      WDES%NB_TOT = NBANDS_OLD 
      RETURN
      END SUBROUTINE


      SUBROUTINE WRITE_FROM_BUF(IU, IREC, IRECL_REAL, OUTBUF)
      INTEGER :: IU         ! unit for writing
      INTEGER :: IREC       ! record to write
      INTEGER :: IRECL_REAL ! record lenght in real(q) words
      REAL(q) :: OUTBUF(:)  ! buffer to written to file

      ! local
      INTEGER :: I

      DO I=1,(SIZE(OUTBUF)+IRECL_REAL-1)/IRECL_REAL
         IREC=IREC+1
         WRITE(12,REC=IREC) OUTBUF((I-1)*IRECL_REAL+1 : MIN(I*IRECL_REAL, SIZE(OUTBUF)))
      ENDDO

      END SUBROUTINE

      SUBROUTINE READ_TO_BUF(IU, IREC, IRECL_REAL, INBUF, IFAIL)
      INTEGER :: IU         ! unit for reading
      INTEGER :: IREC       ! record to read from
      INTEGER :: IRECL_REAL ! record lenght in real(q) words
      REAL(q) :: INBUF(:)   ! buffer to read to
      INTEGER :: IFAIL      ! error condition

      ! local
      INTEGER :: I

      IFAIL=0
      DO I=1,(SIZE(INBUF)+IRECL_REAL-1)/IRECL_REAL
         IREC=IREC+1
         READ(12,REC=IREC, IOSTAT=IFAIL) INBUF((I-1)*IRECL_REAL+1 : MIN(I*IRECL_REAL, SIZE(INBUF)))
         IF (IFAIL/=0) RETURN
      ENDDO

      END SUBROUTINE

      !> compute specific densities and write them to file
      SUBROUTINE WRITE_DENSITIES(IO,XC_DATA,INFO,T_INFO,LATT_CUR,DYN,P,SYMM,&
            GRID,GRID_SOFT,GRIDC,SOFT_TO_C,WDES,W,KINEDEN,CHTOT,DENCOR,RHOLM,CHDEN)

        USE base, ONLY: in_struct, info_struct, symmetry
        USE wave_struct_def, ONLY: wavedes, wavespin
        USE poscar, ONLY: dynamics, type_info, latt
        USE mgrid, ONLY: grid_3d, transit
        USE pseudo_struct_def, ONLY: potcar
        USE tau_mu, ONLY: tau_handle, ALLOCATE_TAU, SET_KINEDEN
        USE setexm_struct_def, ONLY : xc_info
        USE tutor, ONLY: vtutor
        ! from charge.F use CALC_TOTAL_CHARGE, GRADIENT_CHG, LAPLACIAN_CHG

        IMPLICIT NONE

        TYPE(in_struct), INTENT(IN)  :: IO        !< input output unit handle
        ! begin required for WRITE_GRID_QUANTITY
        TYPE(wavedes), INTENT(IN)    :: WDES       !< wavefunction descriptor
        TYPE(dynamics), INTENT(IN)   :: DYN        !< contains position of ions
        TYPE(type_info), INTENT(IN)  :: T_INFO     !< ionic positions in superlattice
        TYPE(latt), INTENT(IN)       :: LATT_CUR   !< lattice parameter
        TYPE(grid_3d), INTENT(IN)    :: GRIDC      !< grid for potentials/charge
        TYPE(info_struct), INTENT(IN):: INFO       !< info about algo etc.
        TYPE(potcar), INTENT(IN)     :: P(:)       !< structure to support radial grid; required for CHGCAR
        REAL(q), INTENT(IN)          :: RHOLM(:,:) !< channel occupancies, (N_MIX_PAW,WDES%NCDIJ); required for CHGCAR
        ! end required for WRITE_GRID_QUANTITY
        ! begin required for ALLOCATE_TAU and SET_KINEDEN
        TYPE(tau_handle), INTENT(IN)  :: KINEDEN   !< kinetic energy density and paramagnetic current density
        TYPE(grid_3d), INTENT(IN)     :: GRID      !< grid for wavefunctions
        TYPE(grid_3d), INTENT(IN)     :: GRID_SOFT !< soft grid for potentials/charge
        TYPE(transit), INTENT(IN)     :: SOFT_TO_C !< index table between GRID_SOFT and GRIDC
        TYPE(xc_info), INTENT(IN)     :: XC_DATA   !< xc settings
        TYPE(symmetry), INTENT(IN)    :: SYMM      !< symmetry information
        TYPE(wavespin), INTENT(IN)    :: W         !< wavefunctions
        ! end required for ALLOCATE_TAU and SET_KINEDEN
        COMPLEX(q), INTENT(IN)         :: CHDEN(:,:) !< pseudo charge density on soft grid
        COMPLEX(q), TARGET, INTENT(IN) :: CHTOT(:,:) !< charge density
        RGRID, INTENT(IN)              :: DENCOR(:)  !< core charge density

        CALL OUTDEN_RESTART_INFORMATION

        CALL OUTDEN_XC_DENSITIES

        !CALL OUTDEN_TO_AECCAR

      CONTAINS

      SUBROUTINE OUTDEN_RESTART_INFORMATION
        !---- write charge density to CHGCAR and/or vaspwave.h5
        CALL WRITE_GRID_QUANTITY("CHGCAR",WDES,DYN%POSION,T_INFO,LATT_CUR,&
                  GRIDC,CHTOT,INFO%SZNAM1,IO%WRT_DENSITY%LCHARG.AND..NOT.IO%LH5,IO%WRT_DENSITY%LCHARGH5,P,RHOLM,INFO%LOVERL)
        !---- write kinetic energy density
        IF (IO%WRT_DENSITY%LTAU) THEN
          IF (IO%WRT_DENSITY%LCHARG.AND..NOT.IO%LH5 .OR. IO%WRT_DENSITY%LCHARGH5) THEN
            ! all good, LCHARG or LCHARGH5 is set
            ! if the kinetic energy density is not computed already, compute it
            IF (.NOT.XC_DATA%LTAU) THEN
               CALL ALLOCATE_TAU(IO%WRT_DENSITY%LTAU, KINEDEN%TAU,KINEDEN%TAUC,KINEDEN%TAUL,GRIDC,WDES%NCDIJ)
               CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
                 IO%WRT_DENSITY%LTAU,T_INFO%NIONS,W,WDES,KINEDEN)
            ENDIF
          ELSE
            ! problem, LCHARG and LCHARGH5 are both false
            IF (.NOT.XC_DATA%LTAU) THEN
               CALL vtutor%alert("Evaluation of the kinetic energy density requested (LTAU=T), but &
                  the output of densities to file (LCHARG/LCHARGH5) is disabled. & 
                  Please enable the density output, e.g. using LCHARG=T, in the INCAR file and &
                  rerun the calculation.")
            ENDIF
          ENDIF
          ! write kinetic energy density to TAUCAR and/or vaspwave.h5
          CALL WRITE_GRID_QUANTITY("TAUCAR",WDES,DYN%POSION,T_INFO,LATT_CUR,&
              GRIDC,KINEDEN%TAU,INFO%SZNAM1,IO%WRT_DENSITY%LCHARG.AND..NOT.IO%LH5,IO%WRT_DENSITY%LCHARGH5)
        ENDIF
      END SUBROUTINE OUTDEN_RESTART_INFORMATION

      SUBROUTINE OUTDEN_XC_DENSITIES
__GPU USE moffload
        COMPLEX(q), ALLOCATABLE,TARGET :: CWORK(:,:),CTMP(:,:)
        INTEGER :: ISP
        COMPLEX(q),ALLOCATABLE:: CGRAD(:,:,:)  ! gradient of charge density
        COMPLEX(q),ALLOCATABLE:: CLAP(:,:)     ! laplacian of charge density
        !---- write total charge (used to call the XC functional) to file
        IF (IO%WRT_DENSITY%LDENSITY.OR.IO%WRT_DENSITY%LGRADIENT.OR.IO%WRT_DENSITY%LLAPLACIAN) THEN
          ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ))
          ! total density returned in CWORK in rec space in VASP units, i.e., [A^{-3}]
          IF (IO%WRT_DENSITY%LWRT_AUGMENTED_DENSITY) THEN
            CALL CALC_TOTAL_CHARGE(LATT_CUR,GRIDC,WDES%NCDIJ,CHTOT,DENCOR,CWORK,.FALSE.,.FALSE.,.FALSE.)
          ELSE
            ALLOCATE(CTMP(GRIDC%MPLWV,WDES%NCDIJ))
GPACC !$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC: CTMP) IF(OFFLOAD_ON)

GPACC !$ACC KERNELS PRESENT(CTMP) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:CTMP) IF(OFFLOAD_ON)
            CTMP=0._q
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS
            DO ISP=1,WDES%NCDIJ
              CALL ADD_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CHDEN(1,ISP),CTMP(1,ISP))
              CALL SETUNB_COMPAT(CTMP(1,ISP),GRIDC)
            ENDDO
            CALL CALC_TOTAL_CHARGE(LATT_CUR,GRIDC,WDES%NCDIJ,CTMP,DENCOR,CWORK,.FALSE.,.FALSE.,.FALSE.)
          ENDIF
          IF (IO%WRT_DENSITY%LDENSITY) THEN
            CALL WRITE_GRID_QUANTITY("charge_density",WDES,DYN%POSION,T_INFO,LATT_CUR,&
              GRIDC,CWORK,LFFT=.True.,LTEXTFILE=.FALSE.)
          ENDIF
        ENDIF
        !---- write gradient of charge density
        IF (IO%WRT_DENSITY%LGRADIENT) THEN
          ALLOCATE(CGRAD(GRIDC%MPLWV,WDES%NCDIJ,3))
          DO ISP=1,WDES%NCDIJ
            CALL GRADIENT_CHG(LATT_CUR,GRIDC,CWORK(1,ISP),CGRAD(1,ISP,1),CGRAD(1,ISP,2),CGRAD(1,ISP,3),.TRUE.,.FALSE.)
          ENDDO
          CALL WRITE_GRID_QUANTITY("gradient",WDES,DYN%POSION,T_INFO,LATT_CUR,&
            GRIDC,&
            RESHAPE(CGRAD,[SIZE(CGRAD,1),SIZE(CGRAD,2)*SIZE(CGRAD,3)]),&
            LTEXTFILE=.FALSE.)
        ENDIF
        !---- write Laplacian of charge density
        IF (IO%WRT_DENSITY%LLAPLACIAN) THEN
          ALLOCATE(CLAP(GRIDC%MPLWV,WDES%NCDIJ))
          DO ISP=1,WDES%NCDIJ
            CALL LAPLACIAN_CHG(LATT_CUR,GRIDC,CWORK(1,ISP),CLAP(1,ISP),.TRUE.,.FALSE.)
          ENDDO
          CALL WRITE_GRID_QUANTITY("laplacian",WDES,DYN%POSION,T_INFO,LATT_CUR,&
            GRIDC,CLAP,LTEXTFILE=.FALSE.)
        ENDIF
      END SUBROUTINE OUTDEN_XC_DENSITIES

      END SUBROUTINE WRITE_DENSITIES

      FUNCTION WRT_DENSITY_TO_STRING(WRT_DENSITY) RESULT (RES)
         USE base, ONLY: write_dens
         IMPLICIT NONE
         TYPE(write_dens), INTENT(IN) :: WRT_DENSITY
         CHARACTER(LEN=:), ALLOCATABLE :: RES
         !
         IF (WRT_DENSITY%LDENSITY) THEN
            RES = "density"
         ELSE
            RES = ""
         END IF

#ifdef __DCRAYFTN_VERSION
         IF (WRT_DENSITY%LGRADIENT) THEN
            RES = RES // " gradient"
            RES = TRIM(RES)
         END IF
         IF (WRT_DENSITY%LLAPLACIAN) THEN
            RES = RES // " laplacian"
            RES = TRIM(RES)
         END IF
#else
         IF (WRT_DENSITY%LGRADIENT) &
            RES = TRIM(RES // " gradient")
         IF (WRT_DENSITY%LLAPLACIAN) &
            RES = TRIM(RES // " laplacian")
#endif
         IF (RES == "") &
            RES = "None"
      END FUNCTION WRT_DENSITY_TO_STRING

!*************************SUBROUTINE OUTCHG ****************************
!
!>   write a component of a grid quantity to a specified unit
!>
!>   HEADER is created in WRITE_GRID_QUANTITY or manually with OUTPOS
!>   grid quantity CHTOT is expected to be in reciprocal space
!>   and is written in real space, i.e., FFT is performed
!>   optionally set LFFT=.FALSE.
!
!***********************************************************************

      SUBROUTINE OUTCHG(GRIDC,IU,LLONG,CHTOT,LFFT)
      USE prec
      USE mpimy
      USE mgrid
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (grid_3d),INTENT(IN)     :: GRIDC              !< grid for quantity
      INTEGER, INTENT(IN)           :: IU                 !< unit
      COMPLEX(q), INTENT(IN)        :: CHTOT(GRIDC%RC%NP) !< real-valued quantity in complex array, e.g. charge density
      LOGICAL, INTENT(IN)           :: LLONG              !< TRUE: E17.11, FALSE: G11.5
      LOGICAL, INTENT(IN), OPTIONAL :: LFFT               !< if T, do FFT (default: T)
! work arrays
      RGRID,  ALLOCATABLE ::  CWORK(:)
      REAL(q),ALLOCATABLE ::  WORK (:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      CHARACTER (40) FORM
      INTEGER ISTAT, NODE_ME, IONODE, N, IERR
      LOGICAL LFFT_
#ifdef __NEC_TUNE__
      INTEGER NHEAD, NFLINES, NTAIL
      CHARACTER (1024) FORM2
#endif // __NEC_TUNE__

      NALLOC=GRIDC%NGX*GRIDC%NGY

      ALLOCATE(CWORK(DIMREAL(GRIDC%MPLWV)),WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the charge immediate exit

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif
      LFFT_=.TRUE.
      IF (PRESENT(LFFT)) LFFT_=LFFT

      IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
        CALL vtutor%bug("internal ERROR: OUTCHG NPLWV is not compatibel with NGX,NGY,NGZ \n " // &
           str(GRIDC%NPLWV) // " " // str(GRIDC%NGX) // " " // str(GRIDC%NGY) // " " // str(GRIDC%NGZ), __FILE__, __LINE__)
      ENDIF

      IF (LFFT_) THEN
        ! transfer to CWORK: copy RC%NP elements
        CALL RC_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
        ! transform to real space representation
        CALL FFT3D(CWORK,GRIDC,1)
      ELSE
         ! just transfer to CWORK in real space: copy RL%NP elements
         CALL RL_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
      ENDIF

      IF (LLONG) THEN
        FORM='(1X,E17.11)'
        NWRITE=5
      ELSE
        FORM='(1X,G11.5)'
        NWRITE=10
      ENDIF

      do_io WRITE(IU,'(3I5)') GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ

      NWRITTEN=0
      DO NZ=1,GRIDC%NGZ
         CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CWORK, NZ)
         io_begin
#ifndef __NEC_TUNE__
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
#else // __NEC_TUNE__
!NEC_RP: generate one format string for a full NALLOC block
         NHEAD = NWRITE - NWRITTEN
         NFLINES = (NALLOC - NHEAD) / NWRITE
         NTAIL = MOD(NALLOC - NHEAD, NWRITE)
         WRITE(FORM2, *) "("
         IF (NHEAD /= 0) THEN
            WRITE(FORM2, *) TRIM(FORM2), NHEAD, TRIM(FORM), ",/,"
         END IF
         IF (NFLINES /= 0) THEN
            WRITE(FORM2, *) TRIM(FORM2), NFLINES, "(", NWRITE, TRIM(FORM), ",/)"
         ENDIF
         IF (NTAIL /= 0) THEN
            WRITE(FORM2, *) TRIM(FORM2), NTAIL, TRIM(FORM)
         END IF
         WRITE(FORM2, *) TRIM(FORM2), ")"

!NEC_RP write all data using one single format string
         NWRITTEN=NTAIL
         WRITE(IU,FORM2,ADVANCE="NO") WORK
#endif // __NEC_TUNE__
         io_end
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(CWORK,WORK)

      RETURN
      END SUBROUTINE

!*************************SUBROUTINE OUTCHG_RGRID **********************
!
!>   write chargedensity to a specified unit
!>
!>   simple copy from OUTCHG, however
!>   no FFT to real space is done, so CHTOT must be provided on 
!>   real-space grid 
!
!***********************************************************************

      SUBROUTINE OUTCHG_RGRID(GRID, IU, LLONG,SV)
         USE poscar_struct_def
         USE prec
         USE mpimy
         USE mgrid
         USE string, ONLY: str
         USE tutor, ONLY: vtutor

         TYPE (grid_3d) GRID

         RGRID  SV(DIMREAL(GRID%MPLWV))
         LOGICAL LLONG        ! long or short format
   ! work arrays
         REAL(q),ALLOCATABLE ::  WORK (:)
         INTEGER NALLOC,NZ, NWRITE, NWRITTEN
         CHARACTER (40) FORM
         INTEGER ISTAT, NODE_ME, IONODE, IU, N, IERR

         NALLOC=GRID%NGX*GRID%NGY

         ALLOCATE(WORK(NALLOC),STAT=ISTAT)
         IF (ISTAT>0) RETURN ! can not write the charge immediate exit

         NODE_ME=0
         IONODE =0
#ifdef MPI
         NODE_ME=GRID%COMM%NODE_ME
         IONODE =GRID%COMM%IONODE
#endif

         IF (GRID%NPLWV/= GRID%NGX*GRID%NGY*GRID%NGZ) THEN
           CALL vtutor%bug("OUTCHG_RGRID NPLWV is not compatibel with NGX,NGY,NGZ \n " // str(GRID%NPLWV) // " " // str(GRID%NGX) // " " // str(GRID%NGY) // " " // str(GRID%NGZ), __FILE__, __LINE__ )
         ENDIF
         ! original subroutine stores CHTOT (complex) to SV (real)
         ! and transforms SV to real-space 
         ! note that SV has 2 times the dimension of CHTOT
         ! so that Re and Im part of CHTOT fits into SV 
         ! this requires that RC_ADD to reside outside of a module 
         ! for the fortran compiler to understand this
!        ! transfer to SV and FFT
!        CALL RC_ADD(CHTOT,1.0_q,CHTOT,0.0_q,SV,GRID)
!        ! transform to real space representation
!        CALL FFT3D(SV,GRID,1)

         IF (LLONG) THEN
           FORM='(1(1X,E17.11))'
           NWRITE=5
         ELSE
           FORM='(1(1X,G11.5))'
           NWRITE=10
         ENDIF

         do_io WRITE(IU,'(3I5)') GRID%NGX,GRID%NGY,GRID%NGZ

         NWRITTEN=0
         DO NZ=1,GRID%NGZ
            CALL MRG_GRID_RL_PLANE(GRID, WORK, SV, NZ)
            io_begin
            DO N=1,NALLOC
               NWRITTEN=NWRITTEN+1
               IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
                  WRITE(IU,FORM) WORK(N)
               ELSE
                  WRITE(IU,FORM,ADVANCE='NO') WORK(N)
               ENDIF
            ENDDO
            io_end
         ENDDO
         IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

         DEALLOCATE(WORK)

         RETURN
   END SUBROUTINE OUTCHG_RGRID


!*************************SUBROUTINE  INCHG ****************************
!>   read charge density form a specified unit
!>
!>   HEADER must have been read previously
!***********************************************************************

      SUBROUTINE INCHG(IU,CHTOT,GRID,IERR,LFFT)
      USE prec
      USE charge
      USE mpimy
      USE mgrid
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      INTEGER,INTENT(IN)            :: IU    !< file unit
      TYPE(grid_3d),INTENT(IN)      :: GRID  !< grid
      COMPLEX(q),INTENT(INOUT)      :: CHTOT(GRID%RC%NP) !< grid quantity read from file
      INTEGER, INTENT(INOUT)        :: IERR  !< error flag
      LOGICAL, INTENT(IN),OPTIONAL  :: LFFT  !< go to reciprocal space
! work array
      REAL(q),ALLOCATABLE ::  WORK(:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      INTEGER ISTAT, NODE_ME, IONODE
      INTEGER NGXFIL, NGYFIL, NGZFIL, N
      LOGICAL LFFT_

      LFFT_ = .TRUE.; IF (PRESENT(LFFT)) LFFT_ = LFFT

      NALLOC=GRID%NGX*GRID%NGY
      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif
      CHTOT=0.0_q

      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) THEN
         IERR=2
         RETURN                 ! can not read potential, immediate exit
      ENDIF

      IERR=0
      IF (GRID%NPLWV/= GRID%NGX*GRID%NGY*GRID%NGZ) THEN
        CALL vtutor%bug("internal ERROR: INCHG GRID%NPLWV is not compatible with GRID%NGX,GRID%NGY,&
           &NGZC \n " // str(GRID%NPLWV) // " " // str(GRID%NGX) // " " // str(GRID%NGY) // " " // &
           str(GRID%NGZ), __FILE__, __LINE__)
      ENDIF

      IERR=0
      NGXFIL=0; NGYFIL=0; NGZFIL=0

#ifdef MPI
      IF (NODE_ME==IONODE) THEN
#endif
        READ(IU,*,ERR=120,END=120) NGXFIL,NGYFIL,NGZFIL
 120    CONTINUE
        IF (NGXFIL==0 .AND. NGYFIL==0 .AND. NGZFIL==0) THEN
           IERR=2
        ELSE IF (NGXFIL/=GRID%NGX .OR. NGYFIL/=GRID%NGY .OR. NGZFIL/=GRID%NGZ ) THEN
          IERR=1
        ENDIF
#ifdef MPI
      ENDIF
#endif

      CALLMPI( M_bcast_i( GRID%COMM, IERR, 1))
      IF (IERR /=0 ) GOTO 200

      NWRITE=5
      NWRITTEN=0
      DO NZ=1,GRID%NGZ
         io_begin
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            ! non advancing read is not the same on all machine
            ! this is the only version that works on all tested platforms
            IF ( MOD(NWRITTEN,NWRITE)==0 .OR. (N==NALLOC .AND. NZ==GRID%NGZ)) THEN
               READ(IU,'(1(1X,E17.11))',ERR=100,END=100) WORK(N)
            ELSE
               READ(IU,'(1(1X,E17.11))',ADVANCE='NO',ERR=100,END=100) WORK(N)
            ENDIF
         ENDDO
         io_end
         CALL DIS_GRID_RL_PLANE(GRID, WORK, CHTOT, .TRUE., NZ )
      ENDDO

      IF (LFFT_) THEN
        CALL FFT_RC_SCALE(CHTOT,CHTOT,GRID)
        CALL SETUNB_COMPAT(CHTOT,GRID)
      ENDIF

  200 CONTINUE
      DEALLOCATE(WORK)
      RETURN

  100 CONTINUE
      IERR=2
      DEALLOCATE(WORK)
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE READNI ****************************
!>   reads total number of ions from a line consisting of numbers
!>   for several species
!***********************************************************************

      SUBROUTINE READNI(ZPARSE,NIONSF)
      USE prec
      USE tutor, ONLY: vtutor
      CHARACTER (80) ZPARSE
      CHARACTER (80) ZITEM
      CHARACTER (80) ZWORK
      CHARACTER (15) ZFORM
      CHARACTER (1)  ZCHAR
      INTEGER IDATA, NDATA, NIONSF, NI, LFORM
      INTEGER, EXTERNAL :: NITEMS

!  get number of data items
      NDATA=NITEMS(ZPARSE,ZWORK,.TRUE.,'I')

!  parse and read list:
      NIONSF=0

      DO 100 IDATA=1,NDATA
         CALL SUBWRD(ZPARSE,ZITEM,IDATA,1)
         CALL CHKINT(ZITEM,ZWORK,ZCHAR,ZFORM)
         IF (ZCHAR=='Y') THEN
            ZWORK='('//ZFORM//')'
            CALL STRIP(ZWORK,LFORM,'A')
            READ(ZITEM,ZWORK(1:LFORM)) NI
         ELSE
            CALL vtutor%error("Fatal error in READNI: Invalid data found ...")
         ENDIF

         NIONSF=NIONSF+NI

  100 CONTINUE

      RETURN
      END SUBROUTINE

!*************************SUBROUTINE READCH ****************************
!
!>  This routine reads in a charge density from file IU
!>
!>  When the charge density of the file is incompatible,
!>  e.g., NGX,Y,Z differ from that used in the programm, a warning is
!>  reported and ICHARG is set to 0.
!>  If the magnetization could not be read ICHARG is set to -1.
!>  On entry to the routine:
!>     CHTOT  may be uninitialised
!>     RHOLM  must be already correctly initialised according to MAGMOM
!>  On exit:
!>     CHTOT  as read from file
!>     RHOLM  as read from file, or unchanged when
!>     KEEP_SMOOTH_ON_MISSING_PAW is true and the optional PAW block is absent
!>  The routine is complicated by several special cases:
!>  o if the magnetization could not be read in
!>     ICHARG is set to -1, CHTOT(:,2:4)=0
!>  o noncollinear case: if only mz could be read
!>     the x and y components are set to 0 in CHTOT and RHOLM
!
!***********************************************************************

      SUBROUTINE READCH(GRIDC, LOVERL, T_INFO, CHTOT, RHOLM, ICHARG, ISPIN, &
           LATT_CUR, P, CSTRF, IU, IU0, KEEP_SMOOTH_ON_MISSING_PAW)
      USE prec
      USE lattice
      USE poscar
      USE mgrid
      USE pseudo
      USE charge
      USE pawm
      USE main_mpi, ONLY: DIR_APP, DIR_LEN

      TYPE (grid_3d), INTENT(IN)  :: GRIDC
      LOGICAL, INTENT(IN)         :: LOVERL
      TYPE (type_info), INTENT(IN):: T_INFO
      INTEGER, INTENT(IN)         :: ISPIN
      TYPE (latt), INTENT(IN)     :: LATT_CUR
      TYPE (potcar), INTENT(IN)   :: P (T_INFO%NTYP)
      COMPLEX(q), INTENT(IN)      :: CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      INTEGER, INTENT(IN)         :: IU
      INTEGER, INTENT(IN)         :: IU0
      LOGICAL, INTENT(IN), OPTIONAL :: KEEP_SMOOTH_ON_MISSING_PAW
      INTEGER, INTENT(INOUT)      :: ICHARG
      REAL(q), INTENT(INOUT)      :: RHOLM(:,:)
      COMPLEX(q)                  :: CHTOT(GRIDC%MPLWV,ISPIN)
! local work arrays
      TYPE (type_info)   T_INFO_OLD
      REAL(q)         :: TMP(T_INFO%NIONS)
      CHARACTER (40) TITEL
      CHARACTER (80) ZPARSE
      CHARACTER (1)  CHARAC
      COMPLEX(q), ALLOCATABLE :: CD(:),CSTRF_OLD(:,:)
      INTEGER I, NIONSF, ISP, N
      INTEGER NODE_ME, IONODE, IERR

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif
      OPEN(UNIT=IU,FILE=DIR_APP(1:DIR_LEN)//'CHGCAR',STATUS='UNKNOWN')

!=======================================================================
! read in old ionic positions and header
!=======================================================================
      T_INFO_OLD=T_INFO

      NULLIFY(T_INFO_OLD%POSION)
      ALLOCATE(T_INFO_OLD%POSION(3,T_INFO_OLD%NIONS))

      READ(IU,*,ERR=1000,END=1000) TITEL
      DO I=1,4
        READ(IU,*,ERR=1000,END=1000)
      ENDDO
      READ(IU,'(A80)') ZPARSE
      READ(ZPARSE,*) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) READ(IU,'(A80)') ZPARSE

      CALL READNI(ZPARSE,NIONSF)
      IF (NIONSF/=T_INFO%NIONS) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: number of atoms are different on CHGCAR file'
        ICHARG=0
        CLOSE(IU)
        RETURN
      ENDIF
      READ(IU,*,ERR=1000,END=1000)

      DO  I=1,NIONSF
        READ(IU,*) T_INFO_OLD%POSION(:,I)
      ENDDO
!=======================================================================
! read in the charge density
!=======================================================================
      ! read charge density in real space and FFT to reciprocal space
      CALL INCHG(IU,CHTOT(1,1),GRIDC,IERR)
      IF (IERR==1) THEN
        IF (IU0>=0) &
        WRITE(*,*) 'WARNING: dimensions on CHGCAR file are different'
        ICHARG=0
        CLOSE(IU)
        RETURN
      ENDIF
      IF (IERR==2) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: CHGCAR file is incomplete'
        ICHARG=0
        CLOSE(IU)
        RETURN
      ENDIF
!=======================================================================
! now subtract charge density according to old position and add
! that one according to new positions
!=======================================================================
      ALLOCATE(CD(GRIDC%MPLWV),CSTRF_OLD(GRIDC%MPLWV,T_INFO%NTYP))

!---- subtract atomic charge density for old positions
      CALL STUFAK(GRIDC,T_INFO_OLD,CSTRF_OLD)
      CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO_OLD,LATT_CUR%B,P,CSTRF_OLD,CD)
      DO N=1,GRIDC%RC%NP
        CHTOT(N,1)= CHTOT(N,1)- CD(N)
      ENDDO

!---- add atomic charge density for new positions
      CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CD)
      DO N=1,GRIDC%RC%NP
        CHTOT (N,1)= CHTOT(N,1)+ CD(N)
      ENDDO

      DEALLOCATE(CD,CSTRF_OLD,T_INFO_OLD%POSION)
!=======================================================================
! read augmentation occupancies
!=======================================================================
      IF (PRESENT(KEEP_SMOOTH_ON_MISSING_PAW)) THEN
        IF (KEEP_SMOOTH_ON_MISSING_PAW) THEN
          IF (IU0>=0) WRITE(IU0,*) &
            'HALF: smooth CHGCAR accepted; using initialized atomic PAW occupancies'
          CLOSE(IU)
          RETURN
        ENDIF
      ENDIF
      CALL RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM(:,1), GRIDC%COMM, IU, IERR )
      IF (IERR/=0) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: PAW occupancies are missing on CHGCAR '
        ICHARG=0
        CLOSE(IU)
        RETURN
      ENDIF
!=======================================================================
! set/read spin density
!=======================================================================
      ! just in case initialize the magnetization to 0
      DO ISP=2,ISPIN
        CALL RC_ADD(CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),GRIDC)
      ENDDO

      IF (IU0>=0) &
        WRITE(IU0,*) 'charge density read from file: CHGCAR, ',TITEL

      DO ISP=2,ISPIN
! read in the spin density
        IERR=0
        ! in vasp.4.4 the ATOMOM array was read from the CHGCAR file
        ! reading the magnetic moments from the file really makes no sense
        ! and spoiles noncollinear calculations
        do_io READ(IU,*,IOSTAT=IERR) (TMP(I),I=1,T_INFO%NIONS)
        CALLMPI( M_bcast_i( GRIDC%COMM, IERR, 1))

        IF (IERR==0) CALL INCHG(IU,CHTOT(1,ISP),GRIDC,IERR)
        IF (IERR==0) CALL RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM(:,ISP), GRIDC%COMM, IU, IERR )

        ! error occured
        IF (IERR/=0) THEN
           IF (ISP==2) THEN
           ! we could not read any entry in the magnetization density return ICHARG=-1
           ! (magnetization according to overlapping atoms)
              ICHARG=-1
              CLOSE(IU)
              RETURN
           ELSE
           ! noncollinear case we could read at least m_z(r)
           ! copy that to the right place and initialise everything else to 0
              CALL RC_ADD(CHTOT(1,2),1.0_q,CHTOT(1,2),0.0_q,CHTOT(1,4),GRIDC) ! copy
              CALL RC_ADD(CHTOT(1,2),0.0_q,CHTOT(1,2),0.0_q,CHTOT(1,2),GRIDC) ! =0
              CALL RC_ADD(CHTOT(1,3),0.0_q,CHTOT(1,3),0.0_q,CHTOT(1,3),GRIDC) ! =0

              RHOLM(:,4)=RHOLM(:,2)
              RHOLM(:,2)=0
              RHOLM(:,3)=0

              CLOSE(IU)
              RETURN
           ENDIF
        ENDIF

        IF (IU0>=0) WRITE(IU0,'(A,I2)') 'magnetization read from file: component ',ISP-1

      ENDDO
      CLOSE(IU)
      RETURN

      ICHARG=-1
      CLOSE(IU)
      RETURN

 1000 CONTINUE
      IF (IU0>=0) &
       WRITE(IU0,*) 'WARNING: CHGCAR is incomplete'
      ICHARG=0
      CLOSE(IU)
      RETURN
      END SUBROUTINE

!*********************** SUBROUTINE READTAU ****************************
!
!>  This routine reads a kinetic energy density from file 
!>      TAUCAR or vaspwave.h5
!>
!>  When the kinetic energy density of the file is incompatible,
!>  e.g., NGX,Y,Z differ from that used in the programm, a warning is
!>  reported and IO%INITAU is set to 0.
!>  If tau_m could not be read IO%INITAU is set to -1.
!> 
!>  Special cases:
!>  o if the tau_m could not be read in
!>     IO%INITAU is set to -1, KINEDEN%TAU(:,2:4)=0
!>  o noncollinear case: if only tau_mz could be read
!>     tau_mx and tau_my are set to 0
!
!***********************************************************************

      SUBROUTINE READTAU(LTAU,GRIDC,T_INFO,LATT_CUR,WDES,IO,KINEDEN)
      USE prec
      USE lattice
      USE poscar
      USE mgrid
      USE pseudo
      USE charge
      USE pawm
      USE wave, ONLY: wavedes
      USE main_mpi, ONLY: DIR_APP, DIR_LEN
      USE tau_mu,   ONLY: tau_handle
      USE base,     ONLY: in_struct
      USE tutor,    ONLY: vtutor

      LOGICAL, INTENT(IN)           :: LTAU     ! T if we should attempt reading tau
      TYPE (grid_3d), INTENT(IN)    :: GRIDC    ! grid
      TYPE (type_info), INTENT(IN)  :: T_INFO
      TYPE (latt), INTENT(IN)       :: LATT_CUR
      TYPE (wavedes), INTENT(IN)    :: WDES
      TYPE(in_struct),INTENT(INOUT) :: IO
      TYPE (tau_handle) :: KINEDEN

      ! local work arrays
      TYPE (type_info)   T_INFO_OLD
      CHARACTER (6) FNAME
      CHARACTER (40) TITEL
      CHARACTER (80) ZPARSE
      CHARACTER (1)  CHARAC
      LOGICAL LTAUCAR_FOUND
      INTEGER IUNIT
      INTEGER I, NIONSF 
      INTEGER NCDIJ ! dimension in spinor space
      INTEGER ISP
      INTEGER N
      INTEGER NODE_ME, IONODE, IERR

      FNAME='TAUCAR'
      IUNIT=18

      ! try early return
      IF (.NOT.LTAU) RETURN
      IF (.NOT.ASSOCIATED(KINEDEN%TAU)) &
       CALL vtutor%bug("Error in READTAU: KINEDEN not allocated.", __FILE__, __LINE__)
      INQUIRE(FILE=DIR_APP(1:DIR_LEN)//FNAME,EXIST=LTAUCAR_FOUND)
      IF (.NOT.LTAUCAR_FOUND) THEN
        IF (IO%IU0>=0) WRITE(IO%IU0,*) 'no '//FNAME//' file present'
        IO%INITAU=0; RETURN
      ENDIF

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif

! reading TAUCAR
      OPEN(UNIT=IUNIT,FILE=DIR_APP(1:DIR_LEN)//FNAME,STATUS='UNKNOWN')

! --- read header: old ionic positions
      T_INFO_OLD=T_INFO

      NULLIFY(T_INFO_OLD%POSION)
      ALLOCATE(T_INFO_OLD%POSION(3,T_INFO_OLD%NIONS))

      READ(IUNIT,*,ERR=1000,END=1000) TITEL
      DO I=1,4
        READ(IUNIT,*,ERR=1000,END=1000)
      ENDDO
      READ(IUNIT,'(A80)') ZPARSE
      READ(ZPARSE,*) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) READ(IUNIT,'(A80)') ZPARSE

      CALL READNI(ZPARSE,NIONSF)
      IF (NIONSF/=T_INFO%NIONS) THEN
        CALL vtutor%warning('Number of atoms are different on '//&
                         FNAME//' file.')
        IO%INITAU=0; RETURN
      ENDIF
      READ(IUNIT,*,ERR=1000,END=1000)

      DO  I=1,NIONSF
        READ(IUNIT,*) T_INFO_OLD%POSION(:,I)
      ENDDO

      IF (ANY(T_INFO_OLD%POSION/=T_INFO%POSION)) THEN
        CALL vtutor%warning('Positions on'//FNAME//' file changed. '//&
          'Currently the extrapolation of the kinetic '//&
          'energy density is not supported.')
      ENDIF

! --- read kinetic energy density
      CALL INCHG(IUNIT,KINEDEN%TAU(:,1),GRIDC,IERR)

      IF (IERR==1) THEN
        CALL vtutor%warning('Dimensions on '//&
                        FNAME//' file are different.')
        IO%INITAU=0; RETURN
      ENDIF
      IF (IERR==2) THEN
        CALL vtutor%warning(FNAME//' file is incomplete.')
        IO%INITAU=0; RETURN
      ENDIF

      IF (IO%IU0>=0) &
        WRITE(IO%IU0,*) 'kinetic energy density read from file: '//&
          FNAME//', '//TITEL

! --- set/read tau_m
      ! initialize tau_m to 0
      DO ISP=2,WDES%NCDIJ
        CALL RC_ADD(KINEDEN%TAU(:,ISP),0.0_q,KINEDEN%TAU(:,ISP),0.0_q,KINEDEN%TAU(:,ISP),GRIDC)
      ENDDO

      DO ISP=2,WDES%NCDIJ
        IERR=0
        IF (IERR==0) CALL INCHG(IUNIT,KINEDEN%TAU(:,ISP),GRIDC,IERR)
        ! error occured
        IF (IERR/=0) THEN
           IF (ISP==2) THEN
              ! we could not read any entry in tau_m
              IO%INITAU=-1; RETURN
           ELSE
           ! noncollinear case we could read at least tau_mz(r)
           ! copy that to the right place and initialize everything else to 0
              CALL RC_ADD(KINEDEN%TAU(1,2),1.0_q,KINEDEN%TAU(1,2),0.0_q,KINEDEN%TAU(1,4),GRIDC) ! copy
              CALL RC_ADD(KINEDEN%TAU(1,2),0.0_q,KINEDEN%TAU(1,2),0.0_q,KINEDEN%TAU(1,2),GRIDC) ! =0
              CALL RC_ADD(KINEDEN%TAU(1,3),0.0_q,KINEDEN%TAU(1,3),0.0_q,KINEDEN%TAU(1,3),GRIDC) ! =0
              RETURN
           ENDIF
        ENDIF

        IF (IO%IU0>=0) &
          WRITE(IO%IU0,'(A,I2)') 'tau_m read from file: component ',ISP-1

      ENDDO
      CLOSE(IUNIT)

      IO%INITAU=1; RETURN

 1000 CONTINUE
      CALL vtutor%warning(FNAME//' file is incomplete.')
      IO%INITAU=0; RETURN
      END SUBROUTINE

!*************************SUBROUTINE READPOT ***************************
!
!>  This routine reads in a chargedensity from file IU
!>
!>  When the chargedensity of the file is incompatible
!>  i.e. NGX,Y,Z differ from that used in the programm a warning is
!>  reported and ICHARG is set to 0
!>  if the magnetisation charge density could not be read in
!>  ICHARG is set to -1
!>  on entry to the routine:
!>     CHTOT  may be uninitialised
!>     RHOLM  must be already correctly initialised according to MAGMOM
!>  on exit:
!>     CHTOT  as read from file
!>     RHOLM  as read from file
!>  the routine is complicated by several special cases:
!>  o if the magnetisation charge density could not be read in
!>     ICHARG is set to -1, CHTOT(:,2:4)=0
!>  o non collinear case: if only mz could be read
!>     the x and y components are set to 0 in CHTOT and RHOLM
!
!***********************************************************************

      SUBROUTINE READPOT(GRIDC, T_INFO, CVTOT, ICHARG, WDES, &
              LATT_CUR, P, LMDIM, CDIJ, N_MIX_PAW, IU, IU0)

      USE prec
      USE lattice
      USE poscar
      USE mgrid
      USE pseudo
      USE charge
      USE pawm
      USE wave

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER LMDIM
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      INTEGER N_MIX_PAW
      INTEGER IU, IU0, ICHARG
! local work arrays
      REAL(q)         :: DLM(N_MIX_PAW,WDES%NCDIJ)
      TYPE (type_info)   T_INFO_OLD
      REAL(q)         :: TMP(T_INFO%NIONS)
      COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      CHARACTER (40) TITEL
      CHARACTER (80) ZPARSE
      CHARACTER (1)  CHARAC
      INTEGER I, NIONSF, N, ISP
      INTEGER NODE_ME, IONODE, IERR

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif

!=======================================================================
! read in old ionic-positions and header
!=======================================================================
      T_INFO_OLD=T_INFO

      NULLIFY(T_INFO_OLD%POSION)
      ALLOCATE(T_INFO_OLD%POSION(3,T_INFO_OLD%NIONS))

      READ(IU,*,ERR=1000,END=1000) TITEL
      DO I=1,4
        READ(IU,*,ERR=1000,END=1000)
      ENDDO
      READ(IU,'(A80)') ZPARSE
      READ(ZPARSE,*) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) READ(IU,'(A80)') ZPARSE

      CALL READNI(ZPARSE,NIONSF)
      IF (NIONSF/=T_INFO%NIONS) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: number of atoms are different on CHGCAR file'
        ICHARG=0
        RETURN
      ENDIF
      READ(IU,*,ERR=1000,END=1000)

      DO  I=1,NIONSF
        READ(IU,*) T_INFO_OLD%POSION(:,I)
      ENDDO
!=======================================================================
! read in the charge-density
!=======================================================================
      CALL INCHG(IU,CVTOT(1,1),GRIDC,IERR)
      CALL FFT3D(CVTOT(1,1),GRIDC,1)

      IF (IERR==1) THEN
        IF (IU0>=0) &
        WRITE(*,*) 'WARNING: dimensions on LOCPOT file are different'
        ICHARG=0
        RETURN
      ENDIF
      IF (IERR==2) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: potential file is incomplete'
        ICHARG=0
        RETURN
      ENDIF

      DEALLOCATE(T_INFO_OLD%POSION)

      ! just in case initialize the magnetization density to 0
      DO ISP=2,WDES%NCDIJ
        CALL RC_ADD(CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
      ENDDO

      DLM(:,1)=0
      CALL RD_RHO_PAW(P, T_INFO, WDES%LOVERL, DLM(:,1), GRIDC%COMM, IU, IERR )
      IF (IERR/=0) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: PAW occupancies are missing on POT file'
        ICHARG=0
        RETURN
      ENDIF
!=======================================================================
! set/read spin density
!=======================================================================
      IF (IU0>=0) &
        WRITE(IU0,*) 'local potential read from file: ',TITEL

      DO ISP=2,WDES%NCDIJ
! read in the spin-density
        IERR=0
        ! in vasp.4.4 the ATOMOM array was read from the CHGCAR file
        ! reading the magnetic moments from the file really makes no sense
        ! and spoiles non collinear calculations
        do_io READ(IU,*,IOSTAT=IERR) (TMP(I),I=1,T_INFO%NIONS)
        CALLMPI( M_bcast_i( GRIDC%COMM, IERR, 1))

        IF (IERR==0) CALL INCHG(IU,CVTOT(1,ISP),GRIDC,IERR)
        CALL FFT3D(CVTOT(1,ISP),GRIDC,1)

        IF (IERR==0) CALL RD_RHO_PAW(P, T_INFO, WDES%LOVERL, DLM(:,ISP), GRIDC%COMM, IU, IERR )

        ! error occured
        IF (IERR/=0) THEN
           IF (ISP==2) THEN
           ! we could not read any entry in the magnetisation density return ICHARG=-1
           ! (magnetisation according to overlapping atoms)
              ICHARG=-1
              RETURN
           ELSE
           ! non collinear case we could read at least m_z(r)
           ! copy that to the right place and initialise everything else to 0
              CALL RC_ADD(CVTOT(1,2),1.0_q,CVTOT(1,2),0.0_q,CVTOT(1,4),GRIDC) ! copy
              CALL RC_ADD(CVTOT(1,2),0.0_q,CVTOT(1,2),0.0_q,CVTOT(1,2),GRIDC) ! =0
              CALL RC_ADD(CVTOT(1,3),0.0_q,CVTOT(1,3),0.0_q,CVTOT(1,3),GRIDC) ! =0

              DLM(:,4)=DLM(:,2)
              DLM(:,2)=0
              DLM(:,3)=0

              GOTO 2000
           ENDIF
        ENDIF

        IF (IU0>=0) WRITE(IU0,'(A,I2)') ' magnetization potential read from file',ISP-1

      ENDDO

 2000 CONTINUE
      ! one centre corrections are passed through mixer and added now
      CALL RETRIVE_RHO_PAW(WDES, P, T_INFO, WDES%LOVERL, LMDIM, CDIJ, DLM)

      RETURN

 1000 CONTINUE
      IF (IU0>=0) &
       WRITE(IU0,*) 'WARNING: potential file is incomplete'
      ICHARG=0
      RETURN
      END SUBROUTINE


!*************************SUBROUTINE OPENGAMMA *************************
!
! the followin subroutines open and close the file GAMMA for reading
! and writing
! furthermore they write and read the header of the file
!
!***********************************************************************

      !> open the file GAMMA for reading and writing
      SUBROUTINE OPENGAMMA
        USE main_mpi
        OPEN(UNIT=IU_GAMMA  ,FILE=DIR_APP(1:DIR_LEN)//'GAMMA',STATUS='UNKNOWN')

      END SUBROUTINE OPENGAMMA


      !> Close the file GAMMA
      SUBROUTINE CLOSEGAMMA

        CLOSE(UNIT=IU_GAMMA)

      END SUBROUTINE CLOSEGAMMA

      !> Write the header to the file GAMMA
      SUBROUTINE WRITEGAMMA_HEAD( NKPTS , NBANDS, IO)
        USE base
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NKPTS   !< number of k-points
        INTEGER :: NBANDS  !< number of bands

        IF (IO%IU0>=0) THEN
           WRITE(IO%IU0,*) 'writing the density matrix to GAMMA'
        ENDIF
        IF (IO%IU6>=0) THEN
           WRITE(IU_GAMMA,'(2I10," ! number of k-points and bands")') NKPTS, NBANDS
        ENDIF

      END SUBROUTINE WRITEGAMMA_HEAD


      !> Write the header to the file GAMMA
      SUBROUTINE WRITEGAMMA_HEAD_HAM( NKPTS , NBANDS, IO)
        USE base
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NKPTS   !< number of k-points
        INTEGER :: NBANDS  !< number of bands

        IF (IO%IU0>=0) THEN
           WRITE(IO%IU0,*) 'writing the density matrix to GAMMA'
        ENDIF
        IF (IO%IU6>=0) THEN
           WRITE(IU_GAMMA,'("Hamiltonian",2I10," ! number of k-points and bands")') NKPTS, NBANDS
        ENDIF

      END SUBROUTINE WRITEGAMMA_HEAD_HAM

      !> Read the header from the file GAMMA
      SUBROUTINE READGAMMA_HEAD( IMODE, NKPTS , NBANDS, IO)
        USE base
        USE tutor, ONLY: vtutor, isError, Gamma
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NKPTS    !< number of k-points
        INTEGER :: NBANDS   !< number of bands
        CHARACTER (LEN=80) ::STRING ! header of Gamma file
        CHARACTER (LEN=1)  ::MODE   !  mode selector
        INTEGER :: IMODE    !< mode: IMODE=1 density matrix ; IMODE=2 Hamiltonian
     ! local
        INTEGER :: NKPTS_   ! number of k-points
        INTEGER :: NBANDS_  ! number of bands
        INTEGER :: IERR

     ! read first line
        READ(IU_GAMMA,'(A)') STRING
     ! try old format first
        IMODE=1
        READ(STRING,*, IOSTAT = IERR) NKPTS_, NBANDS_

        IF (IERR/=0) THEN
     ! try new format
           READ(STRING,*, IOSTAT = IERR) MODE, NKPTS_, NBANDS_
           IF (MODE=='H' .OR. MODE=='h') IMODE=2
        ENDIF

        ! if you do not care, one may simply supply -1 in the file
        IF (NKPTS_==-1)  NKPTS_=NKPTS
        IF (NBANDS_==-1) NBANDS_=NBANDS

        IF (IERR/= 0 .OR. NBANDS /= NBANDS_ .OR. NKPTS_ /= NKPTS) THEN
           CALL vtutor%write(isError, Gamma)
        ENDIF

        IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,*) 'reading the density matrix from GAMMA'
        ENDIF

      END SUBROUTINE READGAMMA_HEAD


!*************************SUBROUTINE WRITEGAMMA *************************
!
!> write the change of the density matrix for one k-point
!>
!> the matrix MUST not be stored in a distributed fashion, and
!> only one core must call the routine
!
!***********************************************************************

      SUBROUTINE WRITEGAMMA( NK, NB_TOT, NBDIM, CHAM , LROT)
        USE prec
        IMPLICIT NONE
        INTEGER :: NK    !< k-point index
        INTEGER :: NB_TOT!< total number of bands
        INTEGER :: NBDIM !< leading dimension of CHAM
        GDEF :: CHAM(NBDIM, NB_TOT)
        !> one can also supply an approximate unitary rotation matrix
        !> for instance obtained by Loewdin perturbation theory
        LOGICAL :: LROT
        INTEGER :: N1, N2
        INTEGER :: NB_START, NB_END

        ! change sign in the upper triangle (convert unitary matrix to density matrix)
        IF (LROT) THEN
           DO N1=1, NB_TOT
              DO N2=N1+1, NB_TOT
                 CHAM(N1, N2)= -CHAM(N1, N2)
              ENDDO
              CHAM(N1, N1)=0
           ENDDO
        ENDIF

        ! write starting band, and final band
        NB_START = 1
        NB_END = NB_TOT
        WRITE(IU_GAMMA, '(3I10," ! k-point index, start and end band")') NK, NB_START, NB_END

        DO N1=NB_START, NB_END
           WRITE(IU_GAMMA,'(16E15.7)') CHAM(N1, NB_START:NB_END)
        ENDDO

        ! change back to rotation matrix)
        IF (LROT) THEN
           DO N1=1, NB_TOT
              DO N2=N1+1, NB_TOT
                 CHAM(N1, N2)= -CHAM(N1, N2)
              ENDDO
              CHAM(N1, N1)=1
           ENDDO
        ENDIF

      END SUBROUTINE WRITEGAMMA


      !> read the change of the density matrix for one k-point
      !>
      !> READGAMMA may be either called by all cores, or better
      !> called only by one core with a subsequent broadcast of the
      !> data to all cores
      SUBROUTINE READGAMMA( NK, NB_TOT, NBDIM, CHAM, IO )
        USE prec
        USE base
        USE tutor, ONLY: vtutor, isError, Gamma2
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NK    !< k-point index
        INTEGER :: NB_TOT!< total number of bands
        INTEGER :: NBDIM !< leading dimension of CHAM
        GDEF :: CHAM(NBDIM, NB_TOT)
      ! local
        INTEGER :: NK_
        INTEGER :: N1, N2
        INTEGER :: NB_START, NB_END
        REAL(q) :: BUFFER(NBDIM*2)

        ! read starting band, and final band
        NB_START = 1
        NB_END = NB_TOT
        READ(IU_GAMMA, *) NK_, NB_START, NB_END

        NB_START=MAX(1, NB_START)
        NB_END  =MIN(NB_TOT, NB_END)
        IF (NK_==-1) NK_=NK

        IF (NK /= NK_) THEN
            CALL vtutor%write(isError, Gamma2)
        ENDIF

        CHAM=0
        DO N1=NB_START,NB_END
#ifdef gammareal
           READ(IU_GAMMA,*) CHAM(N1, NB_START:NB_END)
#else
           ! not to read into real buffer (since complex entries must be in parantheses)
           READ(IU_GAMMA,*) BUFFER(1:(NB_END-NB_START+1)*2)
           ! store as complex
           DO N2=NB_START,NB_END
              CHAM(N1, N2)=CMPLX(BUFFER((N2-NB_START)*2+1),BUFFER((N2-NB_START)*2+2),q)
           ENDDO
#endif

        ENDDO
      END SUBROUTINE READGAMMA

      !> Write header of Magres-formatted output (https://www.ccpnc.ac.uk/docs/magres)
      !> for NMR (LCHIMAG) and EFG (LEFG) calculations
      SUBROUTINE WRITE_MAGRES_HEADER(P,LATT_CUR,INFO,T_INFO,XC,KPOINTS,IU_MAGRES)
         use base
         use prec
         use version
         use pseudo_struct_def
         use poscar
         use mkpoints_struct_def
         use setexm_struct_def, ONLY : xc_info

         IMPLICIT NONE
   
         TYPE (potcar), TARGET :: P(:)
         TYPE (latt) LATT_CUR
         TYPE (info_struct) INFO
         TYPE (type_info) T_INFO
         TYPE (kpoints_struct) KPOINTS
         TYPE (xc_info) XC
         REAL(q) ENMAX
         INTEGER IU_MAGRES
         ! local
         INTEGER NI, I, J
         REAL(q) POS_X, POS_Y, POS_Z

         WRITE(IU_MAGRES, '(A)') "#$magres-abinitio-v1.0"
   
         ! calculation
         WRITE(IU_MAGRES, '(A13)') "[calculation]"
         WRITE(IU_MAGRES, '(A14)') "calc_code VASP"
         WRITE(IU_MAGRES, '(A18,A80)') "calc_code_version ", vasp()

         WRITE(IU_MAGRES, '(2A)') "calc_name ", INFO%SZNAM1 !T_INFO%SZNAM2
         WRITE(IU_MAGRES, '(A12)') "calc_comment"

         ! Functional, Cutoff, Potentials
         WRITE(IU_MAGRES, '(5A)') "calc_xcfunctional ", XC%NAME(:) ! unformatted write
         WRITE(IU_MAGRES, '(A18,(F10.5), A3)') "calc_cutoffenergy ", INFO%ENMAX, " eV"      ! unformatted write

         DO I=1,T_INFO%NTYP
            WRITE(IU_MAGRES, '(10A)') "calc_pspot", P(I)%SZNAMP
         ENDDO

         ! KPOINTS info
         ! KPOINTS read -1 when not specified, only write when bigger than 0
         IF (KPOINTS%NKPX.NE.-1) THEN
            WRITE(IU_MAGRES, '(A20,L6)') "calc_kpoint_mp_gamma", KPOINTS%LGAMMA
            WRITE(IU_MAGRES, '(A19,3(I4))') "calc_kpoint_mp_grid", KPOINTS%NKPX, KPOINTS%NKPY, KPOINTS%NKPZ
            WRITE(IU_MAGRES, '(A21,3(F10.5))') "calc_kpoint_mp_offset", KPOINTS%SHIFT(:)
         ELSE
            WRITE(IU_MAGRES, '(A19,(F10.5))') "calc_kpoint_spacing", KPOINTS%SPACING
         ENDIF

         WRITE(IU_MAGRES, '(A14)') "[/calculation]"

         ! atoms
         WRITE(IU_MAGRES, '(A7)') "[atoms]"
         WRITE(IU_MAGRES, '(A22)') "units lattice Angstrom"
         WRITE(IU_MAGRES, '(A8,9(F10.5))') "lattice ", LATT_CUR%A
         WRITE(IU_MAGRES, '(A19)') "units atom Angstrom"
         ! no symmetry (should be done on the python side)
         !WRITE(IO, *) "symmetry"

         ! Loop structure needs to accomodate for numbering per type
         NI = 0
         DO I=1, T_INFO%NTYP
            DO J=1, T_INFO%NITYP(I)
               NI = NI + 1

               ! express POSION in cartesian coordinates
               POS_X=LATT_CUR%A(1,1)*T_INFO%POSION(1,NI)+ &
               &       LATT_CUR%A(1,2)*T_INFO%POSION(2,NI)+ &
               &        LATT_CUR%A(1,3)*T_INFO%POSION(3,NI)
               POS_Y=LATT_CUR%A(2,1)*T_INFO%POSION(1,NI)+ &
               &       LATT_CUR%A(2,2)*T_INFO%POSION(2,NI)+ &
               &        LATT_CUR%A(2,3)*T_INFO%POSION(3,NI)
               POS_Z=LATT_CUR%A(3,1)*T_INFO%POSION(1,NI)+ &
               &       LATT_CUR%A(3,2)*T_INFO%POSION(2,NI)+ &
               &        LATT_CUR%A(3,3)*T_INFO%POSION(3,NI)

               WRITE(IU_MAGRES, '(A4,A3,A3,I4,3(F10.5))') "atom", P(T_INFO%ITYP(NI))%ELEMENT, P(T_INFO%ITYP(NI))%ELEMENT, J, POS_X, POS_Y, POS_Z !vasp does not have "species label"
            ENDDO
         ENDDO

         WRITE(IU_MAGRES, '(A8)') "[/atoms]"

         RETURN
         END SUBROUTINE WRITE_MAGRES_HEADER

#ifdef VASP_DEBUG
         SUBROUTINE WRITE_STOPCAR(COMM, IU0)

            USE mpimy, ONLY: communic
            USE base, ONLY: in_struct
            TYPE(communic), INTENT(IN) :: COMM
            INTEGER, INTENT(IN)        :: IU0

            INTEGER :: NODE_ME
            INTEGER :: IONODE

#ifdef MPI
            NODE_ME = COMM%NODE_ME
            IONODE = COMM%IONODE
#else
            NODE_ME = 0
            IONODE = 0
#endif

            IF (IU0 >= 0) WRITE(IU0, *) 'DEBUG: Writing STOPCAR and waiting for 3 seconds...'
            io_begin
            OPEN(99, FILE="STOPCAR", STATUS="REPLACE")
            WRITE(99, *) "LSTOP = .TRUE."
            CLOSE(99)
            CALL SLEEP(3)
            io_end

         END SUBROUTINE WRITE_STOPCAR
#endif

!*************************SUBROUTINE OUTPOT ****************************
!>   write potential     to a specified unit
!>   HEADER is currently created in pot.F
!>   z is slowest, x is fastest index
!***********************************************************************

      SUBROUTINE OUTPOT(GRIDC, IU, LLONG, CVTOT)
      USE prec
      USE mgrid
      TYPE(grid_3d), INTENT(IN) :: GRIDC              !< grid for potential
      INTEGER, INTENT(IN)       :: IU                 !< unit
      LOGICAL, INTENT(IN)       :: LLONG              !< TRUE: E17.11, FALSE: G11.5
      RGRID, INTENT(IN)         :: CVTOT(GRIDC%RL%NP) !< real-valued potential
! local variables
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      REAL(q),ALLOCATABLE ::  WORK(:)
      INTEGER ISTAT, NODE_ME, IONODE , N
      CHARACTER (40) FORM

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
#endif

      NALLOC=GRIDC%NGX*GRIDC%NGY

      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the potential immediate exit

      IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
        CALL vtutor%bug("internal ERROR: OUTPOT GRIDC%NPLWV is not compatibel with GRIDC%NGX,&
           &GRIDC%NGY,NGZC \n " // str(GRIDC%NPLWV) // " " // str(GRIDC%NGX) // " " // str(GRIDC%NGY) &
           // " " // str(GRIDC%NGZ), __FILE__, __LINE__)
      ENDIF

      IF (LLONG) THEN
        FORM='(1(1X,E17.11))'
        NWRITE=5
      ELSE
        FORM='(1G11.5)'
        NWRITE=10
      ENDIF

      do_io WRITE(IU,'(3I5)') GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ

      NWRITTEN=0
      DO NZ=1,GRIDC%NGZ
         CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CVTOT, NZ)
         io_begin
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
         io_end
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(WORK)

      RETURN
      END SUBROUTINE OUTPOT


!*********************** SUBROUTINE WRITE_POT **************************
!> @brief: Writes potential to file and/or hdf5 file.
!***********************************************************************
      SUBROUTINE WRITE_OEP(WDES, GRIDC, INFO, T_INFO, LATT_CUR, P, LMDIM, N_MIX_PAW, CVTOT, CDIJ, IO, FILEN )
         USE base
         USE wave_struct_def
         USE mgrid
         USE poscar_struct_def
         USE pseudo_struct_def 
         USE paw
         USE pawm, ONLY: WRT_RHO_PAW
         USE poscar, ONLY: OUTPOS
         USE main_mpi
         USE iso_c_binding, ONLY: c_loc, c_f_pointer 
#ifdef VASP_HDF5
         USE vhdf5, ONLY: IH5WAVEFILEID, VASPWAVE_FOUND, VH5_WRITE_LATTICE_ION, VH5_WRITE_GRID_QUANTITY, &
              VH5_WRITE_PAW_OCCUPANCIES, VH5_DELETE, VH5_EXISTS, VH5_ERROR
         USE string
#endif
         IMPLICIT NONE
         TYPE (wavedes), INTENT(IN) :: WDES  !< wavefunction descriptor 
         TYPE (grid_3d),INTENT(IN) :: GRIDC  !< grid for potentials/charge
         TYPE (info_struct), INTENT(IN) :: INFO
         TYPE (type_info), INTENT(IN) :: T_INFO !< species information
         TYPE (latt), INTENT(IN) :: LATT_CUR !< lattice information
         TYPE (potcar), INTENT(IN) :: P(T_INFO%NTYP) !< potcar handle 
         INTEGER, INTENT(IN) :: LMDIM !< maximum number of all magnetic quantum numbers
         INTEGER, INTENT(IN) :: N_MIX_PAW !< 
         COMPLEX(q), INTENT(INOUT), TARGET :: CVTOT(:,:)!< local potential
         OVERLAP, INTENT(IN) :: CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
         TYPE (in_struct), INTENT(IN) :: IO
         CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FILEN 
         ! local variables
         REAL(q)            :: DLM(N_MIX_PAW,WDES%NCDIJ)
         INTEGER :: I, ISP
         CHARACTER(LEN=80) :: MY_FILE_NAME
         CHARACTER(LEN=40) :: DATASETNAME
         RGRID, POINTER :: RVTOT(:)
 
         MY_FILE_NAME = 'POT'
         IF( PRESENT( FILEN ) ) MY_FILE_NAME = FILEN

         IF (  .NOT. IO%LH5 ) THEN
            IF (IO%IU0>=0) THEN
               IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE=TRIM(MY_FILE_NAME),STATUS='UNKNOWN')
               REWIND IO%IUVTOT
               CALL OUTPOS(IO%IUVTOT,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,T_INFO%POSION)
            ENDIF
            
            ! at the moment the spin up and down potential is written to the file
            CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                 CDIJ, DLM)

            CALL C_F_POINTER( C_LOC( CVTOT(:,1) ), RVTOT, [ GRIDC%RL%NP ] )
            CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,RVTOT)
            CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM(:,1), GRIDC%COMM, IO%IUVTOT)
               
            DO ISP=2,WDES%NCDIJ
               IF (IO%IU0>=0) THEN
                  WRITE( IO%IUVTOT,'(5E20.12)') (T_INFO%ATOMOM(I),I=1,T_INFO%NIONS)
               ENDIF
               CALL C_F_POINTER( C_LOC( CVTOT(:,ISP) ), RVTOT, [ GRIDC%RL%NP ] )
               CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,RVTOT)
               CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM(:,ISP), GRIDC%COMM, IO%IUVTOT )
            ENDDO
            ! close unit
            IF(IO%LOPEN .AND. IO%IU0>=0 ) CLOSE(IO%IUVTOT)
#ifdef VASP_HDF5
         ELSEIF( IO%LH5 ) THEN
            IF( .NOT. VASPWAVE_FOUND ) THEN
               CALL vtutor%bug( "Failed writing "//LOWERCASE(TRIM( MY_FILE_NAME))//", vaspwave.h5 not found!",__FILE__, __LINE__)
            ENDIF

            ! delete existing dataset
            IF (VH5_EXISTS(IH5WAVEFILEID, LOWERCASE( TRIM(MY_FILE_NAME) ), GRIDC%COMM ) ) THEN
               VH5_CHECK( VH5_DELETE(IH5WAVEFILEID, LOWERCASE( TRIM(MY_FILE_NAME) ) ) )
            ENDIF

            ! dump lattice information
            CALL VH5_WRITE_LATTICE_ION(IH5WAVEFILEID," ",T_INFO,LATT_CUR%SCALE,LATT_CUR%A,&
                      .FALSE.,T_INFO%POSION,GROUP=lowercase( TRIM(MY_FILE_NAME) ) ,SUBGROUP="position")

            ! the potential is complex for noncollinear calculations unless represented
            ! in terms of the Pauli matrices
             !CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .FALSE.)

             CALL VH5_WRITE_GRID_QUANTITY(IH5WAVEFILEID, LOWERCASE( TRIM(MY_FILE_NAME) ), 'total', GRIDC, CVTOT)

             ! at the moment the spin up and down potential is written to the file
             CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                  CDIJ, DLM)

             DO ISP = 1, WDES%NCDIJ
                CALL VH5_WRITE_PAW_OCCUPANCIES(IH5WAVEFILEID, LOWERCASE( TRIM(MY_FILE_NAME) ), P, T_INFO, &
                   INFO%LOVERL, DLM(:,ISP), GRIDC%COMM, ISP)
             END DO
             !
             !CALL RC_FLIP_POTENTIAL(CVTOT, GRIDC, WDES%NCDIJ, .TRUE.)

#endif 
         ENDIF
      END SUBROUTINE WRITE_OEP
      END MODULE



!***********************************************************************
!>  write out initial header for PCDAT
!***********************************************************************

      SUBROUTINE PCDAT_HEAD(IU, T_INFO, LATT_CUR, DYN, PACO, SZNAM1)
      USE prec
      USE lattice
      USE poscar
      USE base
      IMPLICIT NONE

      INTEGER IU
      TYPE (type_info) :: T_INFO
      TYPE (latt)::       LATT_CUR
      TYPE (dynamics)  :: DYN
      TYPE (paco_struct)  PACO
      CHARACTER (40) SZNAM1
! local variables
      INTEGER I
      REAL(q) AOMEGA
      AOMEGA=LATT_CUR%OMEGA/T_INFO%NIONS

      WRITE(IU,'(4I4,2E15.7)')1,T_INFO%NIONS,1,0,AOMEGA,DYN%TEMP
      WRITE(IU,*) ' CAR '
      WRITE(IU,*) SZNAM1
      WRITE(IU,'(3I4)') 0,0,0
      WRITE(IU,'(2I4)') 1,DYN%NBLOCK
      WRITE(IU,'(3I4)') PACO%NPACO,PACO%NPACO,PACO%NPACO
      WRITE(IU,'(1I4)') PACO%NPACO
      WRITE(IU,'(1E15.7)') 1E-10
      WRITE(IU,'(1E15.7)') PACO%APACO*1E-10/PACO%NPACO
      WRITE(IU,'(1I4)') DYN%NSW/DYN%NBLOCK/DYN%KBLOCK
      WRITE(IU,'(4E15.7)') DYN%POTIM*1E-15,((LATT_CUR%ANORM(I)*1E-10),I=1,3)

      RETURN
      END SUBROUTINE

!***********************************************************************
!>  write out temperature  to PCDAT file
!***********************************************************************

      SUBROUTINE PCDAT_WRITE(IU, PACO, TMEAN0, TMEAN)
      USE prec
      USE base
      IMPLICIT NONE

      INTEGER IU
      TYPE (paco_struct)  PACO
      REAL(q) :: TMEAN0, TMEAN
! local variables
      INTEGER :: I

      WRITE(IU,'(3E15.7)') TMEAN0,TMEAN

      ! more then one species write all data
      IF (UBOUND(PACO%SIPACO,2)>1) THEN
         DO  I=0,PACO%NPACO-1
            WRITE(IU,'(100F7.3)') PACO%SIPACO(I,:)
         ENDDO
      ELSE
      ! write only total (entry 0 equal entry 1)
         DO  I=0,PACO%NPACO-1
            WRITE(IU,'(100F7.3)') PACO%SIPACO(I,0)
         ENDDO
      ENDIF
      CALL WFORCE(IU)

      END SUBROUTINE

!***********************************************************************
!
! handling of vasp.lock  file
! this file created (open and close)
! there are three functions
! CREATE_VASP_LOCK creates the vasp.lock file on the root node
! REMOVE_VASP_LOCK removes the vasp.lock file on the root node
! WAIT_VASP_LOCK   waits for the file to be created by an external process
!
!***********************************************************************

      !> creates the vasp.lock file on the root node
      SUBROUTINE CREATE_VASP_LOCK(MYCOMM)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM

#ifdef MPI
        IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN
#endif
           OPEN(UNIT=23,FILE=DIR_APP(1:DIR_LEN)//'vasp.lock',STATUS='UNKNOWN')
           WRITE(23,*)'vasp is running happily'
           CLOSE(UNIT=23)
#ifdef MPI
        ENDIF
#endif
      END SUBROUTINE


      !>  waits for the file to be created by an external process
      SUBROUTINE WAIT_VASP_LOCK(MYCOMM)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM

      ! local
        INTEGER :: ISTAT

        DO
           ISTAT=0
           ! try opening file on root node
#ifdef MPI
           IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN
#endif
              OPEN(UNIT=23,FILE=DIR_APP(1:DIR_LEN)//'vasp.lock',STATUS='OLD',IOSTAT=ISTAT)
#ifdef MPI
           ENDIF
#endif
           ! "propagate" error status to all nodes
           CALLMPI( M_sum_i(MYCOMM, ISTAT,1))
           IF (ISTAT==0) EXIT
           ! all nodes sleep for 2 seconds
           CALL SLEEP(1)
           ! all nodes wait for root node
           CALLMPI( M_barrier(MYCOMM))
        ENDDO

        CLOSE(UNIT=23)

      END SUBROUTINE WAIT_VASP_LOCK

!******************** SUBROUTINE FILE_EXISTS ***************************
!> @brief Check if a file exists in the current directory
!***********************************************************************
      SUBROUTINE FILE_EXISTS(MYCOMM, FILENAME, RES)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        LOGICAL, INTENT(OUT) :: RES
        ! local
        LOGICAL :: FILE_PRESENT
        INTEGER :: FILE_INT

#ifdef MPI
           IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN
#endif
           INQUIRE(FILE=TRIM(FILENAME), EXIST=FILE_PRESENT)
           ! cast to int to make broadcast easier
           FILE_INT = MERGE(1, 0, FILE_PRESENT)
#ifdef MPI
           ENDIF
#endif
         CALLMPI(M_bcast_i(COMM%MPI_COMM, FILE_INT , 1 ))
         ! cast back to logical for easier evaluation
         RES = MERGE(.TRUE., .FALSE., FILE_INT == 1)

      END SUBROUTINE FILE_EXISTS

      !> removes the vasp.lock file on the root node
      SUBROUTINE REMOVE_VASP_LOCK(MYCOMM)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM

      ! local
        INTEGER :: ISTAT
#ifdef MPI
        IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN
#endif
           OPEN(UNIT=23,FILE=DIR_APP(1:DIR_LEN)//'vasp.lock',STATUS='UNKNOWN')
           WRITE(23,*) 'file will be removed'
           CLOSE(UNIT=23,STATUS='DELETE')
#ifdef MPI
        ENDIF
#endif
        ! all nodes sleep for 2 seconds to be sure NFS has propagated this change
        CALL SLEEP(1)
        ! all nodes wait for root node
        CALLMPI( M_barrier(MYCOMM))

      END SUBROUTINE REMOVE_VASP_LOCK
