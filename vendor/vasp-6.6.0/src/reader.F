#include "symbol.inc"
      SUBROUTINE READER &
     &       (IU5,IU0,INTERACTIVE,SZNAM1,ISTART,LHALF_INIT,IALGO,IALGO_COMPAT,LROBUST,ISEARCH,IMIX,MAXMIX,MREMOVE, &
     &        AMIX,BMIX,AMIX_MAG,BMIX_MAG,AMIN, &
     &        WC,INIMIX,MIXPRE,MIXFIRST,LFOUND,LDIAG,LSUBROT,LREAL,LREALD, &
     &        LPDENS,IBRION,ICHARG,INIWAV,NELM,NELMALL,NELMIN,NELMDL,EDIFF, &
     &        EDIFFG,NSW,ISIF,IWAVPR,MAXDIS,ISYM,NBLOCK,KBLOCK,ENMAX,POTIM, &
     &        TEBEG,TEEND,NFREE,EFOR, &
     &        NPACO,APACO,NTYPIN,NTYPD,SMASS,SCALEE,POMASS, &
     &        DARWIN_V,DARWIN_R,VCA,LVCADER, &
     &        RWIGS,NELECT,NUP_DOWN,TIME,EMIN,EMAX,EFERMI,FERMI_METHOD,ISMEAR, &
     &        BANDGAP,SPACING,LGAMMA,LKBLOWUP, &
     &        PSTRESS,NDAV,NKRYLOV,DAVITER,NSTRICT,FBREAK, &
     &        SIGMA,EFERMI_NEDOS,LTET,WEIMIN,EBREAK,DEPER,NWRITE,LCORR, &
     &        IDIOT,NIONS,NTYPP,RADIAL_FDORDER,lmusic,LOPTICS,STM, &
     &        ISPIN,ATOMOM,NIOND,LWAVE,LDOWNSAMPLE,SZPREC, &
     &        ENAUG,LORBIT,LELF,ROPT,ENINI, &
     &        NGX,NGY,NGZ,NGXF,NGYF,NGZF,NBANDS,NBANDSLOW,NBANDSHIGH,NEDOS,NBLK,LATT_CUR, &
     &        LPLANE_WISE,LCOMPAT,LMAX_CALC,LMAX_MIX,NSIM,LPAW,LADDGRID, &
     &        LNONCOLLINEAR,LSORBIT,SAXIS, &
     &        LSPIRAL,LZEROZ,QSPIRAL,LORBITALREAL, &
     &        LASPH,TURBO,IRESTART,NREBOOT,NMIN,EREF, &
     &        NLSPLINE,FFTW_PLAN_EFFORT, &
     &        LH5,WRT_POTENTIAL,WRT_DENSITY,LWAP,VELOCITY, &
     &        KERNEL_TRUNCATE, &
     &        EWALD_CUTOFF &
#ifdef libbeef
     &       ,LBEEFENS,LBEEFBAS &
#endif
     &       )


      USE prec
      USE base
      USE sym_prec, ONLY: TINY
      USE ini
      USE lattice
      USE scala
      USE wave_mpi
      USE constant
      USE pseudo   ! for subroutine EXTYP
      USE mkpoints
      USE vaspxml
      USE reader_tags
      USE fermi_energy, ONLY: TOO_SMALL_SIGMA
      USE string, ONLY: str, lowercase
      USE tutor, ONLY: vtutor, never, isBug, isError, isAlert, isWarning, isAdvice
      USE pot_struct_def, only: kernel_truncation
      !USE dynconstr


      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (latt)        LATT_CUR
      TYPE (kernel_truncation) KERNEL_TRUNCATE
#ifdef LONGCHAR
      CHARACTER (32767)INPLIN
#else
      CHARACTER (255)  INPLIN
#endif
      CHARACTER (1)    CHARAC, ALGO
      CHARACTER (40)   SZNAM1
      CHARACTER (40)   SZNAM
      CHARACTER (12)   SZPREC
      CHARACTER (12)   FERMI_METHOD, BANDGAP
      LOGICAL   LDUM,MIXFIRST,LFOUND,LDIAG,LSUBROT,LREAL,LREALD,LPDENS,LTET,LOPTICS, &
     &          LCORR,LOPEN,lmusic,LWAVE,LDOWNSAMPLE, &
     &          LORBIT_,LELF,LCOMPAT,LPAW,LADDGRID, &
     &          LNONCOLLINEAR,LSORBIT, &
     &          LBEEFENS,LBEEFBAS, &
     &          LPLANE_WISE, &
     &          LASPH,INTERACTIVE,LORBITALREAL,LVCADER, &
     &          LCONT, VELOCITY, LH5, LWAP, LHALF_INIT, &
     &          TRUNCATE_COULOMB_KERNEL, LROBUST
      TYPE(write_pot) WRT_POTENTIAL
      TYPE(write_dens) WRT_DENSITY
      DIMENSION POMASS(NTYPD),RWIGS(NTYPP), &
     &          ROPT(NTYPD),DARWIN_V(NTYPD),DARWIN_R(NTYPD),VCA(NTYPD)
      DIMENSION ATOMOM(3*NIOND), ATOMFOR(3*NIOND), EFOR(3,NIOND), TOTFOR(3)
      REAL(q)   SAXIS(3)
      REAL(q)   NELECT,NUP_DOWN
      REAL(q)   STM(7)
      INTEGER   TURBO,IRESTART,NREBOOT,NMIN
      INTEGER   DAVITER
      INTEGER   RADIAL_FDORDER
      REAL(q)   EREF
      REAL(q)   SPACING
      INTEGER   EFERMI_NEDOS
      LOGICAL   LGAMMA
      LOGICAL   LKBLOWUP
      LOGICAL   NLSPLINE
!-MM- Spin spiral stuff
      LOGICAL   LSPIRAL,LZEROZ
      REAL(q)   QSPIRAL(3)
!-MM- end of addition
      INTEGER   FFTW_PLAN_EFFORT
      REAL(q)   MAXDIS

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

! select which levels of error the code should stop on
      INPLIN="--"
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'STOP_ON', INPLIN, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(INPLIN,N,'L')
      CALL LOWER(INPLIN)
      SELECT CASE(INPLIN)
      CASE ('e', 'error')
          vtutor%stopOn = isError
      CASE ('a', 'alert')
          vtutor%stopOn = isAlert
      CASE ('--')
          ! do nothing if stop_on is not set
#ifdef debug
      CASE ('never')
          vtutor%stopOn = never
      CASE ('bug')
          vtutor%stopOn = isBug
      CASE ('warning')
          vtutor%stopOn = isWarning
      CASE ('advice')
          vtutor%stopOn = isAdvice
#endif
      CASE DEFAULT
          CALL vtutor%alert("Couldn't process STOP_ON = " // TRIM(INPLIN) // &
              ", please use STOP_ON = Error or STOP_ON = Alert")
      END SELECT

! 'title'-string (defaults to 'unknown system'), keyword 'SYSTEM'
      SZNAM='unknown system'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SYSTEM', SZNAM, 40, IERR, WRITEXMLINCAR)

      CALL STRIP(SZNAM,N,'L')
      SZNAM1=SZNAM
! start flag ISTART: a default value ISTART=1 should do the best job!
      ISTART=1
! ... of course if 'WAVECAR' doesnt exist --> take ISTART=0 ...
      IF (.NOT.LFOUND) ISTART=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISTART', ISTART, IERR, WRITEXMLINCAR)
      LHALF_INIT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LHALF_INIT', LHALF_INIT, IERR, WRITEXMLINCAR)
      IF (LHALF_INIT) ISTART=0
! the 'idiot flag' (for VTUTOR ...), defaults to 3 ('complete idiot')
      IDIOT=3
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IDIOT', IDIOT, IERR, WRITEXMLINCAR)
      IF (IDIOT<0) IDIOT=0
      IF (IDIOT>3) IDIOT=3
! ... read in the required precision (low - medium - high)
      SZNAM='NORMAL'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PREC', SZNAM, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(SZNAM,N,'L')
      CALL LOWER(SZNAM)
      SZPREC=SZNAM(1:12)
! algorithm: default is 8 (prec. CG)
      IALGO=38
      IALGO_COMPAT=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IALGO', IALGO, IERR, WRITEXMLINCAR)
! algorithm: tag ALGO overwrites IALGO
      INPLIN="--"
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ALGO', INPLIN, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(INPLIN,N,'L')
      CALL LOWER(INPLIN)
      ALGO=INPLIN(1:1)

      IF ( INPLIN(1:4)=='fast') THEN
         IALGO=68
      ELSE IF ( INPLIN(1:5)=='old f') THEN
         IALGO=68
         IALGO_COMPAT=1
      ELSE IF ( INPLIN(1:2)=='of') THEN
         IALGO=68
         IALGO_COMPAT=1
      ELSE IF ( INPLIN(1:2)=='fo') THEN
         IALGO=68
         IALGO_COMPAT=1
      ELSE IF ( INPLIN(1:1)=='f') THEN
         IALGO=68
      ELSE IF ( INPLIN(1:4)=='very') THEN
         IALGO=48
      ELSE IF ( INPLIN(1:5)=='old v') THEN
         IALGO=48
         IALGO_COMPAT=1
      ELSE IF ( INPLIN(1:2)=='ov') THEN
         IALGO=48
         IALGO_COMPAT=1
      ELSE IF ( INPLIN(1:2)=='vo') THEN
         IALGO=48
         IALGO_COMPAT=1
      ELSE IF ( INPLIN(1:1)=='v') THEN
         IALGO=48
      ELSE IF ( INPLIN(1:4)=='none') THEN
         IALGO=2
      ELSE IF ( INPLIN(1:7)=='nothing') THEN
         IALGO=2
      ELSE IF ( INPLIN(1:6)=='normal') THEN
         IALGO=38
      ELSE IF ( INPLIN(1:1)=='n') THEN
         IALGO=38
      ELSE IF ( INPLIN(1:4)=='diag') THEN
         IALGO=90
      ELSE IF ( INPLIN(1:5)=='exact') THEN
         IALGO=90
      ELSE IF ( INPLIN(1:6)=='davidi') THEN
         IALGO=88
      ELSE IF ( INPLIN(1:3)=='dav') THEN
         IALGO=119
      ELSE IF ( INPLIN(1:2)=='cf' .or. INPLIN(1:2)=='af' .or. INPLIN(1:4)=='allf'  ) THEN
         IALGO=108
      ! prevent using IALGO=38 for ALGO = CRPA and NBANDS < 0
      ELSE IF ( INPLIN(1:1)=='c' .AND. .NOT. INPLIN(1:4)=='crpa' ) THEN
         IALGO=58
      ELSE IF ( INPLIN(1:1)=='a' .AND. INPLIN(1:2)/='ac' ) THEN
         IALGO=58
      ELSE IF ( INPLIN(1:6)=='damped') THEN
         IALGO=53
      ELSE IF ( INPLIN(1:1)=='d') THEN
         IALGO=53
      ELSE IF ( INPLIN(1:8)=='eigenval') THEN
         IALGO=3
      ELSE IF ( INPLIN(1:6)=='subrot') THEN
         IALGO=4
      ELSE IF ( INPLIN(1:3)=='sub') THEN
         IALGO=4
      ELSE IF ( INPLIN(1:3)=='jdh' .OR. INPLIN(1:1)=='i') THEN
         IALGO=78
      ENDIF
! read in flag LROBUST
      LROBUST=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LROBUST', LROBUST, IERR, WRITEXMLINCAR)
      IF (LROBUST) THEN
        ISEARCH=1
      ELSE
        ISEARCH=0
      ENDIF
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISEARCH', ISEARCH, IERR, WRITEXMLINCAR)
! factor that controls the maximum size of the reduced basis of the non-blocked Davidson
      NKRYLOV=4
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKRYLOV', NKRYLOV, IERR, WRITEXMLINCAR)
! max. number of iterations in non-blocked Davidson
      DAVITER=25
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'DAVITER', DAVITER, IERR, WRITEXMLINCAR)
! max. number of iterations in non-blocked Davidson
      NSTRICT=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NSTRICT', NSTRICT, IERR, WRITEXMLINCAR)
! break-off criterion w.r.t. residual size in non-blocked Davidson
      FBREAK=1.E-3_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'FBREAK', FBREAK, IERR, WRITEXMLINCAR)
! max. number of iterations NRMM in RMM-DIIS (NDAV), default usually 4
      NDAV=4
      IF (IALGO>=70 .AND. IALGO<90) NDAV=40
      IF (IALGO>=110) NDAV = DAVITER
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NRMM', NDAV, IERR, WRITEXMLINCAR)
! band blocking in RMM-DIIS and Davidson (and some other subroutines)
      NSIM=4
      ! deep iterations, do only few bands in one go
      ! otherwise the diagonalization steps become very expensive
      IF (IALGO>=70 .AND. IALGO<90) NSIM=2
      ! for the Davidson it is advisable to increase blocking
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NSIM', NSIM, IERR, WRITEXMLINCAR)
! LDIAG -- use subspace diagonalization or not (default is TRUE):
      LDIAG=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDIAG', LDIAG, IERR, WRITEXMLINCAR)
! LSUBROT -- use subspace diagonalization or not (default is FALSE):
      LSUBROT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSUBROT', LSUBROT, IERR, WRITEXMLINCAR)
! LADDGRID -- use an additional grid for the calculation of the US-PP
      LADDGRID=.FALSE.
!      IF (SZPREC(1:1)=='a') LADDGRID=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ADDGRID', LADDGRID, IERR, WRITEXMLINCAR)
! read in flag LSORBIT
      LSORBIT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSORBIT', LSORBIT, IERR, WRITEXMLINCAR)
! read in flag LNONCOLLINEAR
      LNONCOLLINEAR=LSORBIT
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LNONCOLLINEAR', LNONCOLLINEAR, IERR, WRITEXMLINCAR)
      IF (LSORBIT) LNONCOLLINEAR=LSORBIT

! ... read spin quantisation axis
      SAXIS(1)=0
      SAXIS(2)=0
      SAXIS(3)=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SAXIS', SAXIS, 3, IERR, WRITEXMLINCAR)
      IF (SQRT(SAXIS(1)**2+SAXIS(2)**2+SAXIS(3)**2)<1E-10) THEN
        CALL vtutor%alert("The magnitude of SAXIS = " // str(SAXIS(1)) // &
             " "// str(SAXIS(2)) // " "// str(SAXIS(3)) // &
             " is too small. Using the default (SAXIS = 0 0 1) instead.")
        SAXIS(1)=0
        SAXIS(2)=0
        SAXIS(3)=1
      ENDIF
#if defined(NGXhalf) ||  defined(NGZhalf)
      IF (LNONCOLLINEAR) THEN
         CALL vtutor%error("ERROR: Noncollinear calculations require "//&
             "that VASP is compiled without the flag -DNGXhalf and"//&
             "-DNGZhalf.\n"//&
             "You probably tried to run a noncollinear calculation "//&
             "(LNONCOLLINEAR=T) with the executable vasp_std. "//&
             "Use vasp_ncl instead!")
      ENDIF
#endif
! spin polarized calculation? (1 is default)
      ISPIN=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISPIN', ISPIN, IERR, WRITEXMLINCAR)
      IF (ISPIN>=2) ISPIN=2
      IF (ISPIN<=1) ISPIN=1
      IF (LNONCOLLINEAR.AND.ISPIN==2) THEN
         ISPIN = 1
         CALL vtutor%warning("You requested both a collinear spin-polarized calculation (ISPIN=2) "//&
             "as well as a noncollinear one (LNONCOLLINEAR=T and/or LSORBIT=T). The noncollinear one "//&
             "will take precedence (defaulting back to ISPIN=1 internally).")
      ENDIF

! Mixing parameters: by default use IMIX=4 (Broyden) with AMIX=0.8,
! BMIX=1.0 (should work almost always ...), WC=100, INIMIX=1, MIXPRE=1
      IMIX=4
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IMIX', IMIX, IERR, WRITEXMLINCAR)

! MIXFIRST mix before diagonalization
      MIXFIRST=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MIXFIRST', MIXFIRST, IERR, WRITEXMLINCAR)

      MAXMIX=-45
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MAXMIX', MAXMIX, IERR, WRITEXMLINCAR)

      MREMOVE=5
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MREMOVE', MREMOVE, IERR, WRITEXMLINCAR)

      AMIX=0.8_q; IF (ISPIN == 2) AMIX = 0.4_q
      IF (LPAW) AMIX=0.4_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AMIX', AMIX, IERR, WRITEXMLINCAR)

      BMIX=1.0_q;
      IF (LPAW) BMIX=1.0_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BMIX', BMIX, IERR, WRITEXMLINCAR)

      AMIX_MAG=AMIX*4
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AMIX_MAG', AMIX_MAG, IERR, WRITEXMLINCAR)

      AMIN=MIN(0.1_q, AMIX, AMIX_MAG)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AMIN', AMIN, IERR, WRITEXMLINCAR)

      BMIX_MAG=BMIX
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BMIX_MAG', BMIX_MAG, IERR, WRITEXMLINCAR)

      WC=100._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WC', WC, IERR, WRITEXMLINCAR)

      INIMIX=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'INIMIX', INIMIX, IERR, WRITEXMLINCAR)

      MIXPRE=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MIXPRE', MIXPRE, IERR, WRITEXMLINCAR)

! initial charge density ICHARG (default 0, if startjob: default 2)
      ICHARG=0
      IF (ISTART==0) ICHARG=2
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ICHARG', ICHARG, IERR, WRITEXMLINCAR)

      LPDENS=.FALSE.
      IF (ICHARG<0) THEN
         ICHARG=0
         LPDENS=.TRUE.
      ENDIF
! initial wavefunctions (defaults is 1, warning: keyword is 'INIWAV')
      INIWAV=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'INIWAV', INIWAV, IERR, WRITEXMLINCAR)

! max/min. number of electronic minimization steps, delay ... (default
! shall be NELM=60, NELMIN=2, NELMDL=-5 if ISTART=0 and
! NELM=60, NELMIN=2, NELMDL=0 if ISTART/=0 ...):
      NELM=60
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELM', NELM, IERR, WRITEXMLINCAR)

      NELMALL=NELM
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELMALL', NELMALL, IERR, WRITEXMLINCAR)

      NELMDL=0
      IF (ISTART==0 .AND. INIWAV==1) THEN
         NELMDL=-5
         IF (IALGO>=40 .AND. IALGO<=50) NELMDL=-12
      ENDIF
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELMDL', NELMDL, IERR, WRITEXMLINCAR)

      NELMIN=2
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELMIN', NELMIN, IERR, WRITEXMLINCAR)

! conjugate gradient or quasi-Newton method? (default IBRION=0)
      IBRION=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IBRION', IBRION, IERR, WRITEXMLINCAR)

! number of degrees of freedom
      IF (IBRION==2) THEN
         NFREE=1
      ELSE
         NFREE=0
      ENDIF
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NFREE', NFREE, IERR, WRITEXMLINCAR)

! energy tolerances (defaults: EDIFF=1E-4, EDIFFG=1E-3)
      EDIFF=1E-4_q
      IF (IBRION==5 .OR. IBRION==6 .OR. IBRION==7 .OR. IBRION==8) EDIFF=1E-6_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EDIFF', EDIFF, IERR, WRITEXMLINCAR)
      ! for reasons of safety (crazy user are present all over the world):
      EDIFF=MAX(ABS(EDIFF),1.E-12_q)

      EDIFFG=EDIFF*10
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EDIFFG', EDIFFG, IERR, WRITEXMLINCAR)
! number of ionic steps, calculate stresses? (default NSW=0, ISIF=2):
      NSW=0
      IF (IBRION==5 .OR. IBRION==6 .OR. IBRION==7 .OR. IBRION==8) NSW=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NSW', NSW, IERR, WRITEXMLINCAR)
! IBRION is 'useless' if NSW=0, set this flag to -1 in this case ...
      IF (NSW==0) IBRION=-1
      ISIF=2
      ! if MD is selected dont calculate stress
      IF (IBRION==0) ISIF=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISIF', ISIF, IERR, WRITEXMLINCAR)
! prediction of wavefunction:
      IWAVPR=0
! MDs
      IF (IBRION==0) IWAVPR=2
#ifdef MPI
      IF (IWAVPR > 0) THEN
         IWAVPR=IWAVPR+10
      ENDIF
#endif
! relaxation: IWAVPR=1
      IF (IBRION>0) IWAVPR=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IWAVPR', IWAVPR, IERR, WRITEXMLINCAR)

      IF (IWAVPR==10) THEN
        ! MD: IWAVPR=12
        IF (IBRION==0) IWAVPR=12
        ! relaxation: IWAVPR=11
        IF (IBRION>0) IWAVPR=11
      ENDIF
      IF (IWAVPR==1) IWAVPR=11 ! makes the same but requires less memory :->

#ifdef MPI
      IWAVPR=MOD(IWAVPR,10)+10
#endif
! read MAXDIS, defaults to 0; no reset ever
      MAXDIS=0.0_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MAXDIS', MAXDIS, IERR, WRITEXMLINCAR)

! switch on symmetry (default ISYM=1):
      ISYM=1 ; IF (LPAW) ISYM=2
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISYM', ISYM, IERR, WRITEXMLINCAR)

! for reasons of safety (crazy user are present all over the world):
      TINY=1E-5
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SYMPREC', TINY, IERR, WRITEXMLINCAR)

! how often to write some data; defaults to KBLOCK=NSW, NBLOCK=1:
      NBLOCK=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBLOCK', NBLOCK, IERR, WRITEXMLINCAR)

      KBLOCK=NSW
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KBLOCK', KBLOCK, IERR, WRITEXMLINCAR)

! plane wave cutoff energy for wavefunctions ..., no default!!!!
      ENMAX=-1._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENMAX', ENMAX, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUT', ENMAX, IERR, WRITEXMLINCAR)

! plane wave cutoff energy for wavefunctions ..., no default!!!!
      ENINI=-1._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENINI', ENINI, IERR, WRITEXMLINCAR)

! cutoff for augmentation charge
      ENAUG=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENAUG', ENAUG, IERR, WRITEXMLINCAR)

! read in NGX, NGY, NGZ, NBANDS
      NGX=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NGX', NGX, IERR, WRITEXMLINCAR)

      NGY=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NGY', NGY, IERR, WRITEXMLINCAR)

      NGZ=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NGZ', NGZ, IERR, WRITEXMLINCAR)

      NGXF=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NGXF', NGXF, IERR, WRITEXMLINCAR)

      NGYF=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NGYF', NGYF, IERR, WRITEXMLINCAR)

      NGZF=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NGZF', NGZF, IERR, WRITEXMLINCAR)

      NBANDS=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDS', NBANDS, IERR, WRITEXMLINCAR)

      NBANDSLOW=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSLOW', NBANDSLOW, IERR, WRITEXMLINCAR)

      NBANDSHIGH=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSHIGH', NBANDSHIGH, IERR, WRITEXMLINCAR)

! ionic time step, default is POTIM=0.5. for IBRION/=0, else no default!
      POTIM=0.5_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'POTIM', POTIM, IERR, WRITEXMLINCAR)

! if IBRION=0 (MD) then POTIM must be given, otherwise error ... !
      IF (((IERR==3).OR.(N<1)).AND.(IBRION==0)) THEN
         CALL vtutor%error('Fatal error! IBRION=0, but no entry for POTIM &
            &on file INCAR. MUST be specified!')
      ENDIF
! start temperature and end temperature (default is 1E-4 for both),
      TEBEG=1.E-4_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'TEBEG', TEBEG, IERR, WRITEXMLINCAR)
      TEEND=TEBEG

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'TEEND', TEEND, IERR, WRITEXMLINCAR)

! pair-correlation functions ..., defaults are NPACO=256, APACO=10
      NPACO=256
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NPACO', NPACO, IERR, WRITEXMLINCAR)

      APACO=10._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'APACO', APACO, IERR, WRITEXMLINCAR)

! NEDOS subdivisions for DOS
      NEDOS=301
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NEDOS', NEDOS, IERR, WRITEXMLINCAR)

! NBLK blocking for some DGEMM commands
      NBLK=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBLK', NBLK, IERR, WRITEXMLINCAR)

! default for SMASS is -3 (micro canonical MD)
      SMASS=-3
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SMASS', SMASS, IERR, WRITEXMLINCAR)

! plane wave cutoff energy for wavefunctions ..., no default!!!!
! default for SMASS is -3 (micro canonical MD)
      SCALEE=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SCALEE', SCALEE, IERR, WRITEXMLINCAR)

! Well, we supply the atomic masses on file POTCAR, but in some cases
! one might wish to change them artificially (--> for example trying
! some kind of 'pre-conditioning' by hand for relaxation runs ...):
! by default we set all masses to negative numbers (this shall be the
! 'signal' to take the values from file POTCAR ...).
! same applies to VCA parameter which are defaulted from POTCAR files
! VCA parameters allow to weigh the potentials by a number supplied as
! VCA (Virtual Crystal Approximation) parameter (usually between 0 and 1)
      POMASS=-1._q
      RWIGS=-1._q
      ROPT=0
      VCA=-1.0_q

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'POMASS', POMASS, NTYPIN, IERR, WRITEXMLINCAR)

! "Cutoff radii" (Wigner-Seitz-radii) for l-projections (default is -1.)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'RWIGS', RWIGS, NTYPP, IERR, WRITEXMLINCAR)
!
! atom weight in virtual crystal approximation (VCA)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'VCA', VCA, NTYPP, IERR, WRITEXMLINCAR)

!
! which order to use for computing derivatives inside the radial grid
      RADIAL_FDORDER = 2
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'RADIAL_FDORDER', RADIAL_FDORDER, IERR, WRITEXMLINCAR)

! LVCADER -- calculate derivative with respect to VCA parameter
! for all ions with VCA not equal 1
      LVCADER=.FALSE.

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LVCADER', LVCADER, IERR, WRITEXMLINCAR)
! read in DARWIN_V and DARWIN_R
      DARWIN_R=0
      DARWIN_V=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'DARWINR', DARWIN_R, NTYPP, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'DARWINV', DARWIN_V, NTYPP, IERR, WRITEXMLINCAR)

! number of up down electrons
      NUP_DOWN=-1._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NUPDOWN', NUP_DOWN, IERR, WRITEXMLINCAR)

! Initial magnetic moments for each atom (default is 1. for all ions)
      AINI=1
      IF (NUP_DOWN >=0) THEN
        AINI=NUP_DOWN/ NIONS
      ENDIF

      NMAGMOM=NIONS
      IF (LNONCOLLINEAR) NMAGMOM=3*NIONS
      DO NI=1,NMAGMOM
         ATOMOM(NI)=AINI
      ENDDO
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MAGMOM', ATOMOM, NMAGMOM, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE., FOUNDNUMBER=N)
      IF (LNONCOLLINEAR .AND. N<3*NIONS .AND. N/=0) THEN
         CALL vtutor%error("ERROR: You have set "//str(N)//" value(s) for MAGMOM; however "//&
                           "for noncollinear calculations (LNONCOLLINEAR=T), "//&
                           "MAGMOM has to be a list of 3*NIONS values. "//&
                           "NIONS="//str(NIONS)//" in your calculation." )
      ENDIF
      IF (ISPIN==2 .AND. N<NIONS .AND. N/=0) THEN
         CALL vtutor%error("ERROR: You have set "//str(N)//" value(s) for MAGMOM; however "//&
                           "for spin-polarized calculations (ISPIN=2), "//&
                           "MAGMOM has to be a list of NIONS values. "//&
                           "NIONS="//str(NIONS)//" in your calculation." )
      ENDIF
      IF ((ISPIN==2.OR.LNONCOLLINEAR) .AND. N==0) THEN
         CALL vtutor%alert("You requested a magnetic or noncollinear calculation, but did not specify the "//&
                           "initial magnetic moment with the MAGMOM tag. Note that a default of "//str(AINI)// &
                           " will be used for all atoms. This ferromagnetic setup may break the symmetry of the "//&
                           "crystal, in particular it may rule out finding an antiferromagnetic solution. Thence, "//&
                           "we recommend setting the initial magnetic moment manually or verifying carefully that "//&
                           "this magnetic setup is desired.")
      END IF

! external forces on each atom..., default is 0. for all ions
      ATOMFOR=0.0_q
      EFOR=0.0_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EFOR', ATOMFOR, NIOND*3, IERR, WRITEXMLINCAR)
      EFOR=RESHAPE(ATOMFOR, SHAPE(EFOR))
      TOTFOR=0.0_q
      DO I=1,NIONS
         TOTFOR=TOTFOR+EFOR(:,I)
      ENDDO
      IF (ANY(ABS(TOTFOR) > 1.E-4_q)) THEN
            CALL vtutor%error('The external forces set &
            &do not add up to zero. Correct this or your system will drift. &
            &If you have fixed atoms to avoid drift using selective dynamics, &
            &please make sure to balance the total forces using the fixed atoms.')
      ENDIF

! number of electrons ..., default is NELECT=0 (= neutral cell)
      NELECT=0._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELECT', NELECT, IERR, WRITEXMLINCAR)

! Real-space projection: default should be POTCAR-dependent ... (if
! one finds 'optimization flag' then set LREAL=.TRUE., else .FALSE.)
      LREAL=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LREAL', LREAL, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE., FOUNDNUMBER=N)

! no input --> remind it and choose later the appropriate value ...
      LREALD=(IERR==3)
      IF (IERR==5) THEN
        INPLIN="--"
        CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LREAL', INPLIN, 40, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
        LREAL=.TRUE.

        CALL STRIP(INPLIN,N,'L')
        IF (INPLIN(1:1)=='O' .OR. INPLIN(1:1)=='o' .OR. &
            INPLIN(1:1)=='A' .OR. INPLIN(1:1)=='a' ) THEN
          IF ( INPLIN(1:1)=='A' .OR. INPLIN(1:1)=='a' ) THEN
            ROPTV=-2E-3
            IF  (SZPREC(1:1)=='l') ROPTV=-1E-2
            IF  (SZPREC(1:1)=='n') ROPTV=-5E-4
            IF  (SZPREC(1:1)=='s') ROPTV=-5E-4
            IF  (SZPREC(1:1)=='h') ROPTV=-4E-4
            IF  (SZPREC(1:1)=='a') ROPTV=-2.5E-4
          ELSE
             ROPTV=1.0_q
             IF  (SZPREC(1:1)=='l') ROPTV=1/1.5
             IF  (SZPREC(1:1)=='h') ROPTV=1.5
          ENDIF
          DO NTYP=1,NTYPIN
            ROPT(NTYP)=ROPTV
          ENDDO
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ROPT', ROPT, NTYPIN, IERR, WRITEXMLINCAR)

          IF ( INPLIN(1:1)=='A' .OR. INPLIN(1:1)=='a' ) THEN
            ROPT=-ABS(ROPT)
          ELSE
            ROPT=ABS(ROPT)
          ENDIF
        ELSE
          IERR=5
        ENDIF
      ENDIF
      LCONT=.FALSE.; CALL CHECK_ERROR_AND_NUMBER('LREAL', IU0, IERR, N, 1, LCONT)
! plane by plane distribution of data
      LPLANE_WISE=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPLANE', LPLANE_WISE, IERR, WRITEXMLINCAR)

! LCOMPAT .TRUE. means full compatibility
      LCOMPAT = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LCOMPAT', LCOMPAT, IERR, WRITEXMLINCAR)
! electronic timestep
      TIME=0.4_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'TIME', TIME, IERR, WRITEXMLINCAR)

! energy range for DOS (default is EMIN=10.,EMAX=-10. = automatic):
      EMIN=10._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EMIN', EMIN, IERR, WRITEXMLINCAR)

      EMAX=-10._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EMAX', EMAX, IERR, WRITEXMLINCAR)

! reference energy read from INCAR
      EREF=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EREF', EREF, IERR, WRITEXMLINCAR)
! determine the method to calculate the Fermi energy
! this can be either LEGACY, MIDGAP or FORCE (which is set by specifying a numeric value)
      INPLIN = "@"
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EFERMI', INPLIN, 12, IERR, WRITEXMLINCAR)
      SELECT CASE (INPLIN(1:1))
      CASE ("m", "M")
         FERMI_METHOD = "MIDGAP"
      CASE ("@", "l", "L")
         FERMI_METHOD = "LEGACY"
      CASE ("0":"9", "-", ".", "f", "F")
         FERMI_METHOD = "FORCE"
      CASE DEFAULT
         CALL VTUTOR%ERROR("Unknown method '" // TRIM(INPLIN) // "' to determine the &
            &Fermi energy. Only the following options are valid: MIDGAP, LEGACY, or &
            &a numeric value.")
      END SELECT
! forced Fermi level read from INCAR (defaults to EREF)
      IF (FERMI_METHOD == "FORCE") THEN
      EFERMI=EREF
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EFERMI', EFERMI, IERR, WRITEXMLINCAR)
      END IF
! z range for STM data
      STM=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'STM', STM, 7, IERR, WRITEXMLINCAR)

! BZ-integration type, default is ISMEAR=1 and SIGMA=0.2 ...
      ISMEAR=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISMEAR', ISMEAR, IERR, WRITEXMLINCAR)

      SIGMA=0.2_q
! If we provide fermi-weights on file INCAR the main intention is mostly
! to do calculations at given fixed occupancies -> this requires SIGMA=0
      IF (ISMEAR==-2) SIGMA=0._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SIGMA', SIGMA, IERR, WRITEXMLINCAR)
! ISMEAR<=-4 and ISMEAR>=30 means tetrahedron method for DOS ...,
! ISMEAR==-4,-5 and <=-7: also tetrahedron method for occ. numbers
! ISMEAR==-14,-15 new tetrahedron method with smearing
      LTET=((ISMEAR<=-4).OR.(ISMEAR>=30))
      IF (ISMEAR==-6) ISMEAR=-1
      IF (ISMEAR>=0) ISMEAR=MOD(ISMEAR,30)
      IF ((SIGMA < TOO_SMALL_SIGMA).AND.(ISMEAR /= -2).AND.(ISMEAR /= -14).AND.(ISMEAR /= -15)) THEN
         CALL vtutor%error("A very small value of SIGMA = " // str(SIGMA) // " was used. &
            &This will lead to problems in the determination of the Fermi level. Please &
            &increase SIGMA to a value larger than " // str(TOO_SMALL_SIGMA) // ".")
      END IF

! Number of points for the Gauss-Legendre integration grid for the determination of the fermi level
      EFERMI_NEDOS = 21
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EFERMI_NEDOS', EFERMI_NEDOS, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BANDGAP', BANDGAP, LEN(BANDGAP), IERR, WRITEXMLINCAR)
      SELECT CASE(BANDGAP(1:1))
      CASE ('c', 'C')
          BANDGAP = 'COMPACT'
      CASE ('w', 'W')
          BANDGAP = 'WEIGHT'
      CASE ('k', 'K')
          BANDGAP = 'KPOINT'
      END SELECT

! k-point spacing
      SPACING=0.5_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KSPACING', SPACING, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
      IF (N>0) LKSPACING_WAS_READ=.TRUE.
! include gamma point in k-points
      LGAMMA=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KGAMMA', LGAMMA, IERR, WRITEXMLINCAR)

! symmetrize k-lattice (the blow-up phase)
      LKBLOWUP=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KBLOWUP', LKBLOWUP, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
      IF (N>0) LKBLOWUP_WAS_READ=.TRUE.

! min. occupation number for 'high quality update' (default: WEIMIN=0)
      WEIMIN=0._q
! MD and relaxation: take WEIMIN=0.001
      IF (IBRION>=0) WEIMIN=0.001_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WEIMIN', WEIMIN, IERR, WRITEXMLINCAR)

! break condition for intra-band min. (default: 0.25*EDIFF/NBANDS)
! because we allow to be EDIFFG smaller than EDIFF also consider
! EDIFFG
      EBREAK=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EBREAK', EBREAK, IERR, WRITEXMLINCAR)

! relative break condition for intra-band minimization (default is 0.3)
      DEPER=0.3_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'DEPER', DEPER, IERR, WRITEXMLINCAR)

! 'verbosity' (default: 2):
      NWRITE=2
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NWRITE', NWRITE, IERR, WRITEXMLINCAR)
   ! allowed range is 0...4, if <0 assume 0, if >4 assume 4 ...
      IF (NWRITE<0) NWRITE=0
      IF (NWRITE>4) NWRITE=4

! Harris corrections for Hellman-Feynman forces ... (default yes):
      LCORR=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LCORR', LCORR, IERR, WRITEXMLINCAR)
! Pullay pressure (one could also say external pressure), default 0.
      PSTRESS=0._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PSTRESS', PSTRESS, IERR, WRITEXMLINCAR)
! max. L for onsite charge expansion in PAW method
      LMAX_CALC=-100
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMAXPAW', LMAX_CALC, IERR, WRITEXMLINCAR)
! max. L for the mixing and CHGCAR for the onsite charge expansion in PAW method
      LMAX_MIX=2
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMAXMIX', LMAX_MIX, IERR, WRITEXMLINCAR)

! some "music" ? (--> default is no ...)
      LMUSIC=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMUSIC', LMUSIC, IERR, WRITEXMLINCAR)
! general h5 io-switch
      LH5 = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LH5', LH5, IERR, WRITEXMLINCAR)
#ifndef VASP_HDF5
      IF (LH5) THEN
        CALL vtutor%alert("The INCAR tag LH5 was set to .TRUE., but this version of VASP was &
            &not compiled with HDF5 support. Please recompile VASP with HDF5 support or avoid &
            &setting LH5 to .TRUE., in order to get rid of this warning. Continuing with LH5=.FALSE.!")
        LH5 = .FALSE.
      ENDIF
#endif

! Sometimes we not interested in any WAVECAR file at all ...
      LWAVE=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWAVE', LWAVE, IERR, WRITEXMLINCAR)

! Sometimes we need to read a WAVECAR of a denser k-grid ...
      LDOWNSAMPLE=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDOWNSAMPLE', LDOWNSAMPLE, IERR, WRITEXMLINCAR)

! ... and maybe not even in any CHGCAR / CHG file ...
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LCHARG', WRT_DENSITY%LCHARG, IERR, WRITEXMLINCAR)

! ... write CHGCAR to h5
      WRT_DENSITY%LCHARGH5 = LH5 .AND. WRT_DENSITY%LCHARG
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LCHARGH5', WRT_DENSITY%LCHARGH5, IERR, WRITEXMLINCAR)
#ifndef VASP_HDF5
      IF (WRT_DENSITY%LCHARGH5) THEN
        CALL vtutor%alert("The INCAR tag LCHARGH5 was set to .TRUE., but this version of VASP was &
            &not compiled with HDF5 support. Please recompile VASP with HDF5 support or avoid &
            &setting LCHARGH5 to .TRUE., in order to get rid of this warning. Continuing with LCHARGH5=.FALSE.!")
        WRT_DENSITY%LCHARGH5 = .FALSE.
      ENDIF
#endif

! ... densities to the HDF5
      INPLIN = '--'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WRT_DENSITY', INPLIN, LEN(INPLIN), IERR, WRITEXMLINCAR)
      INPLIN = LOWERCASE(INPLIN)
      WRT_DENSITY%LDENSITY = INDEX(INPLIN, 'density') > 0
      WRT_DENSITY%LGRADIENT = INDEX(INPLIN, 'gradient') > 0
      WRT_DENSITY%LLAPLACIAN = INDEX(INPLIN, 'laplacian') > 0

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWRT_AUGMENTED_DENSITY', WRT_DENSITY%LWRT_AUGMENTED_DENSITY, IERR, WRITEXMLINCAR)

! ... write ion velocities to hdf5
      CALL PROCESS_INCAR(LOPEN, IO0, IU5, 'VELOCITY', VELOCITY, IERR, WRITEXMLINCAR)

! ... interested in partial charge density ?
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPARD', WRT_DENSITY%LPARD, IERR, WRITEXMLINCAR)
! ... write PARCHG to vaspout.h5 ?
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPARDH5', WRT_DENSITY%LPARCHGH5, IERR, WRITEXMLINCAR)

! ... a WAVECAR must exist
      IF (.NOT.LFOUND.AND.WRT_DENSITY%LPARD) THEN
        CALL vtutor%alert("The value NPARD=T was set, but no WAVECAR file was found. Hence, no partial &
               &charge densities will be written. Continue with LPARD=F.")
        WRT_DENSITY%LPARD = .FALSE.
      ENDIF
! ... but maybe in the total potential?
! ... write the potentials for electron-phonon calculation
      LWAP=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWAP', LWAP, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LVTOT', WRT_POTENTIAL%LVTOT, IERR, WRITEXMLINCAR)

! ... but maybe in the Hartree potential?
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LVHAR', WRT_POTENTIAL%LVHAR, IERR, WRITEXMLINCAR)
      ! if one request the Hartree potential LVHAR supercedes the LVTOT
      IF (WRT_POTENTIAL%LVHAR) WRT_POTENTIAL%LVTOT=.FALSE.

! ... potential to the HDF5
      INPLIN = '--'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WRT_POTENTIAL', INPLIN, LEN(INPLIN), IERR, WRITEXMLINCAR)
      INPLIN = LOWERCASE(INPLIN)
      WRT_POTENTIAL%TOTAL = INDEX(INPLIN, 'total') > 0
      WRT_POTENTIAL%HARTREE = INDEX(INPLIN, 'hartree') > 0
      WRT_POTENTIAL%IONIC = INDEX(INPLIN, 'ionic') > 0
      WRT_POTENTIAL%XC = INDEX(INPLIN, 'xc') > 0
      WRT_POTENTIAL%XCMU = INDEX(INPLIN, 'xcmu') > 0
      WRT_POTENTIAL%PAW = INDEX(INPLIN, 'paw') > 0

! ... check if truncating the Coulomb kernel is needed (open boundary conditions)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/LTRUNCATE', &
                         KERNEL_TRUNCATE%ACTIVE, IERR, WRITEXMLINCAR)
      IF (KERNEL_TRUNCATE%ACTIVE) THEN
          KERNEL_TRUNCATE%PADDING_FACTOR = HUGE(1)
          KERNEL_TRUNCATE%TRUNCATION_FACTOR = HUGE(1._q)
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/IDIMENSIONALITY', &
                         KERNEL_TRUNCATE%DIMENSIONALITY, IERR, WRITEXMLINCAR)
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/LCOARSEN', &
                         KERNEL_TRUNCATE%COARSEN_BEFORE_PAD, IERR, WRITEXMLINCAR)
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/ISURFACE', &
                         KERNEL_TRUNCATE%SURFACE_NORMAL_DIRECTION, IERR, WRITEXMLINCAR)
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/LSIMPLE2D', &
                         KERNEL_TRUNCATE%USE_SIMPLIFIED_2D_KERNEL, IERR, WRITEXMLINCAR)
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/FACTOR', &
                         KERNEL_TRUNCATE%TRUNCATION_FACTOR, IERR, WRITEXMLINCAR)
          CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KERNEL_TRUNCATION/IPAD', &
                         KERNEL_TRUNCATE%PADDING_FACTOR, IERR, WRITEXMLINCAR)
      ENDIF

! read in flag LORBIT
      LORBIT_=.FALSE.
      LORBIT=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LORBIT', LORBIT_, IERR, LWRITEXML=.FALSE., LCONTINUE=.TRUE., FOUNDNUMBER=N)
      IF (IERR==5) THEN
        CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LORBIT', LORBIT, IERR, LWRITEXML=.FALSE., FOUNDNUMBER=N)
      ELSE
         IF (LORBIT_) LORBIT=5
      ENDIF

      LCONT=.FALSE.; CALL CHECK_ERROR_AND_NUMBER('LORBIT', IU0, IERR, N, 1, LCONT)

! read in flag LELF
      LELF=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LELF', LELF, IERR, WRITEXMLINCAR)

! read in flag LOPTICS
      LOPTICS=.FALSE.
#ifndef noPAW
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LOPTICS', LOPTICS, IERR, WRITEXMLINCAR)
#endif
! if scaLAPACK is used it can be switched of in the INCAR file
      IF (LscaLAPACK) THEN
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSCALAPACK', LSCALAPACK, IERR, WRITEXMLINCAR)
      IF (IU0>0) THEN
         IF (LscaLAPACK) THEN
#ifndef ACC_OFFLOAD
            WRITE(IU0,*) 'scaLAPACK will be used'
#else
            WRITE(IU0,*) 'scaLAPACK will be used selectively (only on CPU)'
#endif
         ELSE
            WRITE(IU0,*) 'scaLAPACK is switched off'
         ENDIF
      ENDIF
      ENDIF

      LSCALU= .FALSE.
! the parallel LU decomposition might be slower than the serial
! one, hence we can switch it off
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSCALU', LSCALU, IERR, WRITEXMLINCAR)

      LSCAAWARE=LscaLAPACK
      LscaAWARE_read=.FALSE.
! the parallel LU decomposition might be slower than the serial
! one, hence we can switch it off
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSCAAWARE', LSCAAWARE, IERR, WRITEXMLINCAR)
      IF (IERR==0) LscaAWARE_read=.TRUE.

! try to overlap communication with calculations ?
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LASYNC', LASYNC, IERR, WRITEXMLINCAR)

! read in flag LASPH
      LASPH=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LASPH', LASPH, IERR, WRITEXMLINCAR)
! read in flag LORBITALREAL
      LORBITALREAL=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LORBITALREAL', LORBITALREAL, IERR, WRITEXMLINCAR)

!-MM- spin spiral stuff
! if LSPIRAL
      LSPIRAL=.FALSE.
      LZEROZ =.FALSE.
      QSPIRAL=0._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSPIRAL', LSPIRAL, IERR, WRITEXMLINCAR)

      IF (LSPIRAL.AND..NOT.LNONCOLLINEAR) THEN
        CALL vtutor%error("Spin spirals require LNONCOLLINEAR=.TRUE. \n Exiting VASP.")
      ENDIF
      IF (LSPIRAL.AND.LSORBIT) THEN
        CALL vtutor%error("Spin spirals (LSPIRAL=True) cannot be combined with "//&
           "spin-orbit coupling (LSORBIT=True). \n Exiting VASP.")
      ENDIF

! if LSPIRAL=.TRUE. we also need QSPIRAL, and possibly LZEROZ
      IF (LSPIRAL) THEN
! ... read propagation vector of spin spiral
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'QSPIRAL', QSPIRAL, 3, IERR, WRITEXMLINCAR)
! ... look for LZEROZ
         LZEROZ=.TRUE.
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LZEROZ', LZEROZ, IERR, WRITEXMLINCAR)
      ENDIF
! read in flag INTERACTIVE
      INTERACTIVE=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'INTERACTIVE', INTERACTIVE, IERR, WRITEXMLINCAR)
      IF (INTERACTIVE) IBRION=11
! read in flag TURBO
      TURBO=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'TURBO', TURBO, IERR, WRITEXMLINCAR)
! read in flag IRESTART
      IRESTART=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IRESTART', IRESTART, IERR, WRITEXMLINCAR)
! read in flag NREBOOT
      NREBOOT=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NREBOOT', NREBOOT, IERR, WRITEXMLINCAR)
! read in flag NMIN
      NMIN=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NMIN', NMIN, IERR, WRITEXMLINCAR)
! read in flag NLSPLINE
      NLSPLINE=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NLSPLINE', NLSPLINE, IERR, WRITEXMLINCAR)

! read the FFTWMAKEPLAN flag
      INPLIN="--"
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'FFTWMAKEPLAN', INPLIN, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(INPLIN,N,'L'); CALL LOWER(INPLIN)
      IF (INPLIN(1:1)=='e') FFTW_PLAN_EFFORT=0
      IF (INPLIN(1:1)=='m') FFTW_PLAN_EFFORT=1

! read EWALD_CUTOFF, the outer bound for the spline function used in FEWALD (ewald_params%max_argument)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EWALD_CUTOFF', EWALD_CUTOFF, IERR, WRITEXMLINCAR)

#ifdef libbeef
! read in flag LBEEFENS
      LBEEFENS=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LBEEFENS', LBEEFENS, IERR, WRITEXMLINCAR)

! read in flag LBEEFBAS
      LBEEFBAS=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LBEEFBAS', LBEEFBAS, IERR, WRITEXMLINCAR)
#endif
! Thats all from INCAR (for the first ...):
      CALL CLOSE_INCAR_IF_FOUND(IU5)
      RETURN

  150 CONTINUE
      CALL vtutor%error("Error code was IERR= " // str(IERR) // " ... . Found N= " // str(N) // " data.")

    END SUBROUTINE READER

    !
    !> Duplicated reading of mpi flags, since writing to hdf5 and xml not
    !> possible at the position in code were they are already needed.
    !
    SUBROUTINE MPI_READER(IU5,IU0)

      USE reader_tags
      USE prec
      USE base
      USE main_mpi, ONLY: COMM_KIN, COMM_WORLD, NCORE_ACT => NCORE
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      use incar_reader, ONLY: defined_in_image
      IMPLICIT NONE

      LOGICAL LOPEN
      INTEGER IU0, IU5
      INTEGER IERR, IDUM
      INTEGER NPAR, NPAR_ACT, FOUND_NPAR
      INTEGER NCORE, FOUND_NCORE, FOUND_NCORE_BAND
      LOGICAL LREAD_VAR
      INTEGER IMAGES
      INTEGER KIMAGES
      INTEGER FOURORBIT
      INTEGER KPAR
      LOGICAL LVCAIMAGES
      REAL(q) VCAIMAGES
#if defined(use_shmem) && !defined(_OPENMP)
      INTEGER NCSHMEM
#endif
#ifdef USENCCL
      LOGICAL LUSENCCL
#endif

#if defined(MPI) || defined(MPI_CHAIN)
      NPAR_ACT = COMM_KIN%NCPU / NCORE_ACT
      NCORE = NCORE_ACT
      NPAR = NPAR_ACT
#endif

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

      VCAIMAGES=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'VCAIMAGES', VCAIMAGES, IERR, WRITEXMLINCAR)
      IF (VCAIMAGES==-1) THEN
         LVCAIMAGES=.FALSE.
      ELSE
         LVCAIMAGES=.TRUE.
      ENDIF

      IF (LVCAIMAGES) THEN
         IMAGES=2

         ! Here we only check for strings that are forbidden in IMAGE_* in the INCAR file
         CALL FORBID_TAG_IN_IMAGE('NCORE_IN_IMAGE1')
         CALL FORBID_TAG_IN_IMAGE('ML_OUTBLOCK')
         ! ML_FF_OUTBLOCK is a deprecated alternative to ML_OUTBLOCK.
         CALL FORBID_TAG_IN_IMAGE('ML_FF_OUTBLOCK')
         CALL FORBID_TAG_IN_IMAGE('NBLOCK')   
         CALL FORBID_TAG_IN_IMAGE('IBRION')
         CALL FORBID_TAG_IN_IMAGE('ISIF')
         CALL FORBID_TAG_IN_IMAGE('MD_ALGO')
         CALL FORBID_TAG_IN_IMAGE('NSW')
         CALL FORBID_TAG_IN_IMAGE('KBLOCK')
         CALL FORBID_TAG_IN_IMAGE('TEBEG')
         CALL FORBID_TAG_IN_IMAGE('TEEND')
         CALL FORBID_TAG_IN_IMAGE('POTIM')
         CALL FORBID_TAG_IN_IMAGE('ML_ESTBLOCK')
         ! ML_IERR is a deprecated alternative to ML_ESTBLOCK.
         CALL FORBID_TAG_IN_IMAGE('ML_IERR')

      ELSE
         IMAGES=0
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IMAGES', IMAGES, IERR, WRITEXMLINCAR)
      ENDIF


      KIMAGES=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KIMAGES', KIMAGES, IERR, WRITEXMLINCAR)
      IF (KIMAGES>0) THEN
         IDUM=0
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'FOURORBIT', IDUM, IERR, WRITEXMLINCAR)
      ENDIF

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPAR', KPAR, IERR, WRITEXMLINCAR)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NCORE', NCORE, IERR, WRITEXMLINCAR, FOUNDNUMBER=FOUND_NCORE)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NCORES_PER_BAND', NCORE, IERR, WRITEXMLINCAR, FOUNDNUMBER=FOUND_NCORE_BAND)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NPAR', NPAR, IERR, WRITEXMLINCAR, FOUNDNUMBER=FOUND_NPAR)

#if defined(MPI) || defined(MPI_CHAIN)
#if !defined(_OPENMP) && !defined(ACC_OFFLOAD)
      IF ((NCORE /= NCORE_ACT).AND.(FOUND_NCORE + FOUND_NCORE_BAND > 0)) THEN
         CALL vtutor%alert("The value NCORE = " // str(NCORE) // " specified in &
            &the INCAR file was overwritten, because it was not compatible with &
            &the " // str(COMM_KIN%NCPU) // " processes available:\n&
            &NCORE = " // str(NCORE_ACT) // &
            "\nwas used instead, please check that this makes sense for &
            &your machine.")
      END IF
      IF ((NPAR /= NPAR_ACT).AND.(FOUND_NPAR > 0)) THEN
         CALL vtutor%alert("The value NPAR = " // str(NPAR) // " specified in &
            &the INCAR file was overwritten, because it was not compatible with &
            &the " // str(COMM_KIN%NCPU) // " processes available:\n&
            &NPAR = " // str(NPAR_ACT) // &
            "\nwas used instead, please check that this makes sense for &
            &your machine.")
      END IF
#endif
#endif
#if defined(use_shmem) && !defined(_OPENMP)
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NCSHMEM', NCSHMEM, IERR, WRITEXMLINCAR)
#endif
#ifdef USENCCL
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LUSENCCL', LUSENCCL, IERR, WRITEXMLINCAR)
#endif
      CALL CLOSE_INCAR_IF_FOUND(IU5)
      RETURN

      CONTAINS
      
         SUBROUTINE FORBID_TAG_IN_IMAGE(TAGNAME)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: TAGNAME
            ! Local variables
            LOGICAL :: PRESENT_IN_IMAGE
            PRESENT_IN_IMAGE = defined_in_image(INCAR_F, TAGNAME)
            CALLMPI( M_or(COMM_WORLD, PRESENT_IN_IMAGE, 1) )
            IF (PRESENT_IN_IMAGE) CALL vtutor%error( TAGNAME // &
               " must not be set in IMAGE block.")
         END SUBROUTINE FORBID_TAG_IN_IMAGE
    END SUBROUTINE MPI_READER

    !
    !> Duplicated reading of OpenMP flags, since writing to hdf5 and xml not
    !> possible at the position in code were they are already needed.
    !
    SUBROUTINE OMP_READER(IU5,IU0)
      USE reader_tags
      USE prec
      USE base

      IMPLICIT NONE

      INTEGER IU0, IU5
    ! local
      LOGICAL LOPEN
      INTEGER IERR
      INTEGER OMP_NTHREADS
      INTEGER OMP_NTHREADS_ALLTOALL
      INTEGER OMP_NTHREADS_NONLR_RSPACE

#ifdef _OPENMP
      INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS
      OMP_NTHREADS=OMP_GET_MAX_THREADS()
#else
      OMP_NTHREADS=1
#endif

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

      OMP_NTHREADS_ALLTOALL=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NTHREADS_A2A', OMP_NTHREADS_ALLTOALL, IERR, WRITEXMLINCAR)

      OMP_NTHREADS_NONLR_RSPACE=OMP_NTHREADS
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NTHREADS_RPROJ', OMP_NTHREADS_NONLR_RSPACE, IERR, WRITEXMLINCAR)

      CALL CLOSE_INCAR_IF_FOUND(IU5)

      RETURN
    END SUBROUTINE OMP_READER

    SUBROUTINE POTCAR_H5_READER(IU5, IU0, POTENTIALTYP)
      USE reader_tags
      IMPLICIT NONE
      LOGICAL :: LOPEN
      INTEGER, INTENT(IN) :: IU5, IU0
      CHARACTER(LEN=4), INTENT(OUT) :: POTENTIALTYP
      INTEGER :: IERR, N

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
      POTENTIALTYP = 'PBE'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'POTENTIALTYP', POTENTIALTYP, 40, IERR, LWRITEXML=.FALSE.)
      CALL CLOSE_INCAR_IF_FOUND(IU5)
    END SUBROUTINE POTCAR_H5_READER
