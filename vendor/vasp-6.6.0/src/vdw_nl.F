#include "symbol.inc"

MODULE vdw_df_helper
  INTEGER :: IU0=6,IU6=6,ISIF=0
CONTAINS
  SUBROUTINE SET_vdW_IO_UNITS(IU0_in, IU6_in, ISIF_in)
    INTEGER IU0_in, IU6_in, ISIF_in
    IU0=IU0_in
    IU6=IU6_in
    ISIF=ISIF_in
  END SUBROUTINE SET_vdW_IO_UNITS
END MODULE vdw_df_helper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> This is used to evaluate the non-local
!> correlation functional of Dion et al., PRL 2004.
!>
!> The fast procedure of Roman-Perez and Soler PRL 2010 is used
!> which is available in Siesta and GPAW which this implementation follows.
!> Written by jik 2009-2011
!> Please cite \n
!>   G. Rom\'{a}n-P\'{e}rez and J. Soler, Phys. Rev. Lett. 103, 096102 (2009) \n
!>   J. Klime\v{s}, D. R. Bowler, and A. Michaelides, PRB 83, 195131 (2011) \n
!> and \n
!>   J. Klime\v{s}, D. R. Bowler, and A. Michaelides, J. Phys.: Cond. Matt. 22, 022201 (2010) \n
!> if you use optPBE-vdW or optB88-vdW functionals.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE VDW_LL
  USE prec
  USE lattice
  USE mpimy
  USE mgrid
  USE constant
  USE vdw_df_helper, ONLY : IU0,IU6,ISIF

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: VDW_NONLOC, VDW_NONLOC_spin

  INTEGER, parameter :: Nalp=30, Nr=4096

  REAL(q) :: spl(Nalp,Nalp,4)
  REAL(q) :: phik(0:Nr,Nalp,Nalp), d2phidk2(0:Nr,Nalp,Nalp)
  REAL(q) :: phir(0:Nr,Nalp,Nalp)
  REAL(q) :: qs_vdw(Nalp)
  REAL(q), parameter :: q0cut=10.0_q, lambda=1.2_q
  REAL(q), parameter :: qAlp_ratio=31.949333280885_q
  REAL(q), parameter :: rcut = 100._q
  REAL(q), parameter :: asymptot_phi=12._q*(4._q*pi/9._q)**3
  REAL(q) :: dr
  LOGICAL :: INIT_VDW_LL=.TRUE. !do we need to initialise the stuff?
  LOGICAL :: LREAD_KERN = .FALSE.
  LOGICAL :: LHAVE_PHI_TABLE = .FALSE.
  logical :: lhave_alpha_table = .FALSE.

  !dmesh parameters, normal settings
  !integer, parameter :: nd = 20
  !real(q), parameter :: dcut = 30._q
  !real(q), parameter :: dratio = 20._q
  !dmesh parameters, tight
  integer, parameter :: nd = 35
  real(q), parameter :: dcut = 50._q
  real(q), parameter :: dratio = 20._q
  logical :: luse_asymptotic = .true.
  !careful !!! the kernel needs to be recalculated when
  !luse_asymptotic is changed (the highest coefficients are
  !zero for .false.)
  !normal settings, double nintg for tight
  INTEGER, PUBLIC :: nintg
  REAL(q), parameter :: maxin =200._q
  !REAL(q), parameter :: maxin =100._q


  real(q) :: dmesh(nd)
  logical :: lhave_nd_points = .FALSE.
  real(q) :: kcut
  real(q), parameter :: dsoft = 0.5_q, phi_0 =0.8_q

  real(q) :: phi_table(0:3,0:3,nd,nd)

CONTAINS

  SUBROUTINE VDW_NONLOC(DCHARG,DWORK,DWORKG,DWORK4,EXC,GRIDC,LATT_CUR,stress,ISPIN,IVDW_NL,Zab_vdW,GAMMA_VDW,ALPHA_VDW)

__GPU USE moffload_struct_def

    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    REAL(q) Zab_vdW,GAMMA_VDW,ALPHA_VDW
    REAL(q) DCHARG(GRIDC%RL%NP) !charge density on input and output
    REAL(q) DWORKG(GRIDC%RL%NP) !df/drho input, add the same non-local here
    REAL(q) DWORK4(GRIDC%RL%NP) !gradient on input !df/d(drho)/(drho)potential on output
    REAL(q) EXC
    REAL(q) stress(3,3)
    REAL(q) E_corr
    COMPLEX(q), allocatable :: CINP(:) !, THETA(:,:)
    REAL(q) ,allocatable :: k_k(:)
    RGRID   DWORK(GRIDC%RL%NP)
    INTEGER ISPIN, IVDW_NL

    IF (INIT_VDW_LL) THEN
      call MIN_K(GRIDC,LATT_CUR)
      call GET_SPLINES
      call GET_KERNELS(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
      INIT_VDW_LL=.FALSE.
    ENDIF

    ALLOCATE(K_K(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(K_K) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:K_K) IF(OFFLOAD_ON)

    call GET_K_IND(k_k,GRIDC,LATT_CUR)

    ALLOCATE(CINP(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CINP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CINP) IF(OFFLOAD_ON)

    CALL VDW_CALC(DCHARG,DWORK,DWORKG,DWORK4,EXC,GRIDC,LATT_CUR,CINP,CINP,spl, phik, d2phidk2,qs_vdw,k_k,stress,ISIF,E_corr,ISPIN,Zab_vdW)

GPACC !$ACC EXIT DATA DELETE(CINP,K_K) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CINP, K_K) IF(OFFLOAD_ON)
    DEALLOCATE(CINP)
    DEALLOCATE(k_k)

    IF (IU6>=0) THEN
       WRITE(IU6,'(A,F14.7)') 'Total vdW correction in eV: ',E_corr
     ENDIF
  END SUBROUTINE VDW_NONLOC

  SUBROUTINE VDW_NONLOC_spin(DCHARG,DCHARGS,DWORK,DWORKS,DWORKG,DWORKGS, &
                 & DWORK4,DWORK4S,EXC,GRIDC,LATT_CUR,stress,IVDW_NL,Zab_vdW,GAMMA_VDW,ALPHA_VDW)

__GPU USE moffload_struct_def

    !spin polarised version of the previous routine
    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    REAL(q) Zab_vdW,GAMMA_VDW,ALPHA_VDW
    REAL(q) DCHARG(GRIDC%RL%NP), DCHARGS(GRIDC%RL%NP) !charge density on input and output
    REAL(q) DWORKG(GRIDC%RL%NP), DWORKGS(GRIDC%RL%NP) !df/drho input, add the same non-local here
                    !don't touch this one yet
    REAL(q) DWORK4(GRIDC%RL%NP), DWORK4S(GRIDC%RL%NP) !gradient on input
                    !df/d(drho)/(drho)potential on output
    REAL(q) EXC
    REAL(q) stress(3,3) !this guy is not yet calculated and passed to xcspin.F
    REAL(q) E_corr
    COMPLEX(q), allocatable :: CINP(:) !, THETA(:,:)
    REAL(q) ,allocatable :: k_k(:)
    RGRID   DWORK(GRIDC%RL%NP), DWORKS(GRIDC%RL%NP)
    INTEGER IVDW_NL

    IF (INIT_VDW_LL) THEN
      call MIN_K(GRIDC, LATT_CUR)
      call GET_SPLINES
      call GET_KERNELS(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
      INIT_VDW_LL=.FALSE.
    ENDIF

    ALLOCATE(K_K(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(K_K) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:K_K) IF(OFFLOAD_ON)

    call GET_K_IND(k_k,GRIDC,LATT_CUR)

    ALLOCATE(CINP(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CINP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CINP) IF(OFFLOAD_ON)

    CALL VDW_CALC_spin(DCHARG,DCHARGS,DWORK,DWORKS,DWORKG,DWORKGS,DWORK4,DWORK4S, &
              EXC,GRIDC,LATT_CUR,CINP,CINP,spl, phik, d2phidk2,qs_vdw,k_k,stress,ISIF,E_corr,Zab_vdW)

GPACC !$ACC EXIT DATA DELETE(CINP,K_K) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CINP, K_K) IF(OFFLOAD_ON)

    DEALLOCATE(CINP)
    DEALLOCATE(k_k)

    IF (IU6>=0) THEN
       WRITE(IU6,'(A,F14.7)') 'Total vdW correction in eV: ',E_corr
     ENDIF
  END SUBROUTINE VDW_NONLOC_spin

  SUBROUTINE GET_SPLINES
    !construct the spline coefficients for interpolation
    REAL(q) :: a(Nalp,Nalp), c(Nalp,Nalp), l(Nalp,Nalp), mi(Nalp,Nalp), z(Nalp,Nalp)
    REAL(q) :: b(Nalp-1,Nalp), d(Nalp-1,Nalp), h(Nalp-1), al(Nalp-1,Nalp)
    INTEGER :: i,j

    spl=0._q
    a=0._q
    b=0._q
    d=0._q
    h=0._q
    al=0._q
    c=0._q
    l=0.0_q
    mi=0._q
    z=0._q

    call GEN_ALPHA_POINTS()

    DO i=1,Nalp
      a(i,i)=1._q
    ENDDO
    DO i=1,Nalp-1
      h(i)=qs_vdw(i+1)-qs_vdw(i)
    ENDDO
    DO i=2,Nalp-1
      DO j=1,Nalp
        al(i,j)=3._q/h(i)*(a(i+1,j)-a(i,j))-3._q/h(i-1)*(a(i,j)-a(i-1,j))
      ENDDO
    ENDDO
    l(1,:)=1._q
    DO i=2,Nalp-1
      l(i,:)=2._q*(qs_vdw(i+1)-qs_vdw(i-1))-h(i-1)*mi(i-1,:)
      mi(i,:)=h(i)/l(i,:)
      z(i,:)=(al(i,:)-h(i-1)*z(i-1,:))/l(i,:)
    ENDDO
    l(Nalp,:)=1._q
    z(Nalp,:)=0._q
    c(Nalp,:)=0._q
    DO i=Nalp-1,1,-1
      c(i,:)=z(i,:)-mi(i,:)*c(i+1,:)
      b(i,:)=(a(i+1,:)-a(i,:))/h(i)-h(i)/3._q*(c(i+1,:)+2._q*c(i,:))
      d(i,:)=(c(i+1,:)-c(i,:))/(3._q*h(i))
    ENDDO

    spl(:,:,1)=a(:,:)
    spl(:,:,3)=c(:,:)
    DO i=1,Nalp-1
      spl(i,:,2)=b(i,:)
      spl(i,:,4)=d(i,:)
    ENDDO

  END SUBROUTINE GET_SPLINES

    !> This procedure calculates the fourier transformed vdw kernels.
  SUBROUTINE GET_KERNELS(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
    IMPLICIT NONE
    INTEGER :: alf,bet,i,j,irad, nk, ik, IVDW_NL
    REAL (q) :: q1,q2, dk, krad
    REAL (q) :: r(0:nr)
    REAL (q) :: phir_work(0:nr), phik_work(0:nr)
    REAL (q) :: d2phik_work(0:nr)
    REAL (q) :: x_work(0:nr), u_work(0:nr)
    REAL (q) :: GAMMA_VDW,ALPHA_VDW
    TYPE (grid_3d) GRIDC

    !setting up the mesh of interpolating points
    !call GEN_ALPHA_POINTS()
    dr= rcut/Nr
    dk= pi/rcut
    nk= int(kcut/dk)+1

    !find the interpolation coefficients for the alpha points
    DO alf = 1, Nalp
      DO bet = 1, alf
        phik_work=0._q
        phir_work=0._q
        !desaturate q interpolation points
        q1=DESATURATE(qs_vdw(alf))
        q2=DESATURATE(qs_vdw(bet))

        !get the kernel on a real grid by using the precalculated interpolation
        !data
        DO irad = 0, Nr
          r(irad) = irad * dr
          phir(irad,alf,bet) = PHI_INTERPOLATE(GRIDC,q1*r(irad),q2*r(irad),IVDW_NL,GAMMA_VDW,ALPHA_VDW)
        ENDDO

        !smoothen for large r
        DO irad =0, Nr
          phir(irad,alf,bet)=phir(irad,alf,bet)*S_CUT(r(irad)/rcut)
        ENDDO

        !fourier transform it and smoothen
        phir_work(:)=phir(:,alf,bet)
        call PHIK_FFT(phir_work(:), phik_work(:))
        phik(:,alf,bet)=phik_work*(2._q*pi)**1.5_q
        phik(nk:nr,alf,bet)=0._q
        DO ik =1,nk
          krad=ik*dk
          phik(ik,alf,bet)=phik(ik,alf,bet)*S_CUT(krad/kcut)
        ENDDO

        !get coefficients for cubic splines
        !phik_work to splines
        call SPLINE(dk,phik(:,alf,bet),Nr+1,0._q,0._q,d2phik_work,x_work,u_work)

        !symmetrise
        phik(:,bet,alf)=phik(:,alf,bet)
        d2phidk2(:,bet,alf)=d2phik_work
        d2phidk2(:,alf,bet)=d2phik_work
      ENDDO
    ENDDO

  END SUBROUTINE GET_KERNELS

    !> Find the interpolation points for the kernel interpolation q_0 grid
    !> i.e. fill the qs_vdw array.
  SUBROUTINE GEN_ALPHA_POINTS()
    REAL(q) :: a,b,e
    INTEGER :: i

    qs_vdw(:)=0._q
    a=log(qAlp_ratio)/(Nalp-1._q)
    b=q0cut/(qAlp_ratio-1._q)

    DO i=2,Nalp-1
      qs_vdw(i)=b*exp(a*(i-1._q))-b
    ENDDO
    qs_vdw(Nalp)=q0cut
    lhave_alpha_table = .true.

  END SUBROUTINE GEN_ALPHA_POINTS

    !> Find q0 corresponding to saturated q0 by interval halving.
  FUNCTION DESATURATE(x)
    REAL(q) :: desaturate, x
    REAL(q) :: xdo, xup, xme, xmesat
    REAL(q), parameter :: xtol=1.e-14_q

    !check first
    IF ((x<0._q).OR.(x>q0cut)) then
      CALL vtutor%error("DESATURATE: out of range " // str(x))
    endif
    !this is done for the original q0 values before saturation
    xdo =  0._q
    xup = q0cut*5._q
    DO
      xme = (xdo+xup)*0.5_q
      xmesat = SATURATE(xme)
      !if smaller
      IF (abs(xmesat-x).LE.xtol ) THEN
        desaturate=xup
        RETURN
      ELSE IF (x.LE.xmesat) THEN
        xup=xme
      ELSE
        xdo=xme
      ENDIF
    ENDDO

  END FUNCTION DESATURATE

    !> The same that we use down in the main routine
    !> R-P&S' equation No 4 or so.
  FUNCTION SATURATE(x)
    REAL(q):: saturate, x, coe
    INTEGER :: iu

    coe=0._q
    DO iu=1,12
      coe=coe+(x/q0cut)**iu/real(iu,kind=q)
    ENDDO
    saturate=q0cut*(1._q-exp(-coe))
  END FUNCTION SATURATE

    !> Finds the kernel using phi_generate
    !> returns the value of kernel for d1, d2.
  FUNCTION PHI_INTERPOLATE(GRIDC, d1, d2, IVDW_NL, GAMMA_VDW, ALPHA_VDW)
    REAL (q) :: phi_interpolate
    REAL (q) :: d1, d2
    REAL (q) :: dd1(0:3), dd2(0:3)
    REAL (q) :: GAMMA_VDW, ALPHA_VDW
    INTEGER :: id1, id2
    INTEGER :: n1, n2
    INTEGER :: IVDW_NL
    logical,save :: lfirst_G = .TRUE.
    TYPE (grid_3d) GRIDC

    !the interpolation table needs to be generated before we access it
    IF (.NOT. lhave_phi_table) call PHI_GENERATE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
    !don't go over the max
      IF ((d1.GE.dcut).OR.(d2.GE.dcut)) THEN
      !or use C/d^6
        IF (.NOT.luse_asymptotic) THEN
          phi_interpolate=0._q
        ELSE
          !if d1 or d2 are 0 then the thing doesn't work,
          !but as my boss of long time ago (Fiala, UFA) said:
          !"cut-off at 1, it's close to infinity"
          !the terms here are on the order of 10^-6 compared to the big ones
          !so that the errors should be small (when d1 or d2 are not that big)
          !and presumably not cutting the thing off will help? we'll see...
          IF ((d1.GE.dmesh(2)) .AND. (d2.GE.dmesh(2))) THEN
            phi_interpolate=-asymptot_phi/(d1*d1*d2*d2*(d1*d1+d2*d2))
          ELSE
            phi_interpolate=0._q
          ENDIF
        ENDIF
      RETURN
    ENDIF

    !find the indices and interpolate the value
    id1 = IOFD(d1)
    id2 = IOFD(d2)
    dd1(0) = 1._q
    dd2(0) = 1._q
    dd1(1) = (d1 - dmesh(id1)) / ( dmesh(id1+1) - dmesh(id1))
    dd2(1) = (d2 - dmesh(id2)) / ( dmesh(id2+1) - dmesh(id2))
    dd1(2) = dd1(1)**2
    dd2(2) = dd2(1)**2
    dd1(3) = dd1(1)**3
    dd2(3) = dd2(1)**3

    phi_interpolate = 0._q
    DO n2 = 0,3
      DO n1 = 0,3
        phi_interpolate = phi_interpolate + phi_table(n1,n2,id1,id2)*dd1(n1)*dd2(n2)
      ENDDO
    ENDDO

  END FUNCTION PHI_INTERPOLATE

    !> Find the interpolation points for the dmesh grid
    !> i.e. fill the dmesh array.
  SUBROUTINE GEN_ND_POINTS()
    REAL(q) :: a,b,e
    INTEGER :: i

    dmesh(:)=0._q
    a=log(dratio)/(nd-1._q)
    b=dcut/(dratio-1._q)

    DO i=1,nd-1
      dmesh(i)=b*exp(a*(i-1._q))-b
    ENDDO
    dmesh(nd)=dcut
    lhave_nd_points = .true.

  END SUBROUTINE GEN_ND_POINTS

    !> Finds closest smallest index point for the kernel interpolation
    !> ok so the point is that we have numbers b*(exp(a*(i-1))-1).
  FUNCTION IOFD(d)
    INTEGER :: iofd
    REAL (q) :: d
    REAL (q), save :: a,b
    LOGICAL, save :: lfirst = .true.

    IF (lfirst) THEN
      a = log(dratio)/(nd-1._q)
      a = max( a,1.e-13_q)
      b = dcut/(dratio-1._q)
      lfirst = .false.
    ENDIF

    iofd = 1+log(1+(d-dmesh(1))/b)/a
    iofd = max ( 1, iofd)
    iofd = min (nd-1, iofd)

  END FUNCTION IOFD

    !> Generate the interpolation table for the kernel using phi.
  SUBROUTINE PHI_GENERATE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
    logical :: lkern,lkern_new
    integer :: IVDW_NL,IVDW_NL_TMP,NINTG_TMP,nmesh,IOstatus
    real(q) :: GAMMA_VDW,GAMMA_VDW_TMP,ALPHA_VDW,ALPHA_VDW_TMP
    TYPE (grid_3d) GRIDC

    !ok if the kernel table is here read it
    inquire( file='vdw_kernel.bindat', exist=lkern)
    IF (lkern) then
       open (unit=10, file='vdw_kernel.bindat', form='unformatted')
       read(10,IOSTAT=IOstatus) IVDW_NL_TMP,GAMMA_VDW_TMP,ALPHA_VDW_TMP,NINTG_TMP,nmesh
       IF (IOstatus.EQ.0) THEN
          lkern_new=.true.
       ELSE
          rewind(10)
          read(10,IOSTAT=IOstatus) nmesh
          IF (IOstatus.EQ.0) lkern_new=.false.
          IF (IVDW_NL.NE.1) call vtutor%error('PHI_GENERATE: The distributed vdw_kernel.bindat file can only be used with IVDW_NL=1.')
       ENDIF
       IF (IOstatus.EQ.0) THEN
          IF (lkern_new) THEN
             IF ((IVDW_NL_TMP.NE.IVDW_NL).OR.(GAMMA_VDW_TMP.NE.GAMMA_VDW).OR.(ALPHA_VDW_TMP.NE.ALPHA_VDW).OR.(NINTG_TMP.NE.NINTG).OR.(nmesh.NE.nd)) THEN
                IF (IU0>=0) WRITE(IU0,*) 'Parameters in vdw_kernel.bindat: IVDW_NL,GAMMA_VDW,ALPHA_VDW,NINTG,NMESH=',IVDW_NL_TMP,GAMMA_VDW_TMP,ALPHA_VDW_TMP,NINTG_TMP,nmesh
                IF (IU0>=0) WRITE(IU0,*) 'Parameters of required vdW kernel: IVDW_NL,GAMMA_VDW,ALPHA_VDW,NINTG,NMESH=',IVDW_NL,GAMMA_VDW,ALPHA_VDW,NINTG,nd
                IF (IU0>=0) WRITE(IU0,*) 'Parameters differ. Recalculating the vdW kernel (this might take time).'
                close(10)
                CALL GET_KERNEL_TABLE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
                lhave_phi_table=.true.
                RETURN
             ENDIF
          ELSE
             IF (nmesh.NE.nd) THEN
                IF (IU0>=0) WRITE(IU0,*) 'Parameter in vdw_kernel.bindat: NMESH=',nmesh
                IF (IU0>=0) WRITE(IU0,*) 'Parameter of required vdW kernel: NMESH=',nd
                IF (IU0>=0) WRITE(IU0,*) 'Parameters differ. Recalculating the vdW kernel (this might take time).'
                close(10)
                CALL GET_KERNEL_TABLE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
                lhave_phi_table=.true.
                RETURN
             ENDIF
          ENDIF
          read(10,IOSTAT=IOstatus) dmesh
          read(10,IOSTAT=IOstatus) phi_table
          IF (IOstatus.NE.0) THEN
             IF (IU0>=0) WRITE(IU0,*) 'vdw_kernel.bindat is not legible. Recalculating the vdW kernel (this might take time).'
             close(10)
             CALL GET_KERNEL_TABLE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
          ELSE
             IF (IU6>=0) WRITE(IU6,'(A/,(10F8.4))') 'vdW kernel read from vdw_kernel.bindat', dmesh(:) !, phi_table(0,2,:,5)
             IF (IU6>=0) WRITE(IU6,*)
          ENDIF
       ELSE
          IF (IU0>=0) WRITE(IU0,*) 'vdw_kernel.bindat is not legible. Recalculating the vdW kernel (this might take time).'
          close(10)
          CALL GET_KERNEL_TABLE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
       ENDIF
       close(10)
    ELSE
       IF (IU0>=0) WRITE(IU0,*) 'No vdw_kernel.bindat file found. Calculating vdW kernel (this might take time).'
       CALL GET_KERNEL_TABLE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
    ENDIF
    lhave_phi_table=.true.

  END SUBROUTINE PHI_GENERATE

    !> Here we want to find bicubic interpolation coefficients for the vdW kernel
    !> and store them in variable
    !> phi_table(0:3,0:3,nd,nd).
    !> That is we need to get the values, derivations and cross derivs for the phi(d1,d2)
    !> function as defined in Dion et al and RPetS papers.
  SUBROUTINE GET_KERNEL_TABLE(GRIDC,IVDW_NL,GAMMA_VDW,ALPHA_VDW)
    INTEGER :: id1,id2,id,id_tmp,nd2,id1vec((nd*nd-nd)/2+nd),id2vec((nd*nd-nd)/2+nd),IVDW_NL
    INTEGER :: NODE_ME, IONODE
    REAL (q) :: d1, d2, dtot, d1m, d1p, d2m, d2p, dd
    REAL (q) :: val1
    REAL (q) :: phi(nd,nd), phi_der_d1(nd,nd), phi_der_d2(nd,nd), phi_der_d1d2(nd,nd)
    REAL (q) :: ker_d1p_d2, ker_d1m_d2, ker_d1_d2p, ker_d1_d2m
    REAL (q) :: ker_d1p_d2p, ker_d1m_d2p, ker_d1p_d2m, ker_d1m_d2m
    REAL (q) :: by(4), by_d1(4), by_d2(4), by_d1d2(4), coefs(0:3,0:3)
    REAL (q) :: GAMMA_VDW, ALPHA_VDW
    TYPE (grid_3d) GRIDC

#ifdef MPI
    NODE_ME=GRIDC%COMM%NODE_ME
    IONODE =GRIDC%COMM%IONODE
#else
    NODE_ME=0
    IONODE =0
#endif

    !dmesh is read from the kernel file if that exists and is legible
    !in that case this routine is not called and dmesh needs to be done
    IF (.NOT. lhave_nd_points) THEN
      CALL GEN_ND_POINTS()
    ENDIF

    do id1=3,nd-1
      d1=dmesh(id1)
      phi(id1,id1)=GET_KER_VAL(d1,d1,IVDW_NL,gamma_vdw,alpha_vdw)
!      IF (IU0>=0) WRITE(IU0,*)  d1, phi(id1,id1)
    enddo

    id=0
    DO id1=1,nd
      DO id2=1,id1
         id=id+1
         id1vec(id)=id1
         id2vec(id)=id2
      ENDDO
    ENDDO

    phi=0._q
    phi_der_d1=0._q
    phi_der_d2=0._q
    phi_der_d1d2=0._q
    dd=0.01_q
    nd2=(nd*nd-nd)/2+nd
    !we'll call GET_KER_VAL(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw) to get the phi(d1,d2) value
    DO id_tmp=1, nd2, GRIDC%COMM%NCPU
      id=id_tmp+GRIDC%COMM%NODE_ME-1
      IF (id.LE.nd2) THEN
        id1=id1vec(id)
        id2=id2vec(id)
        !IF (IU0>=0) WRITE(IU0,*)  id1, id2
        !the thing is symmetric so we can do half the work
        d1 = dmesh(id1)
        d2 = dmesh(id2)
        !value of phi
        phi(id1,id2)=GET_KER_VAL(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw)
        !derivative
        !now dd should be chosen in a sensible way, the minimal values are 0.00 0.27  0.585
        !and the largest 21.5 25.4 30.0
        d1p=d1+dd
        d1m=d1-dd
        ker_d1p_d2=GET_KER_VAL(d1p,d2,IVDW_NL,gamma_vdw,alpha_vdw)
        ker_d1m_d2=GET_KER_VAL(d1m,d2,IVDW_NL,gamma_vdw,alpha_vdw)
        phi_der_d1(id1,id2)=(ker_d1p_d2 -ker_d1m_d2 )/(2.d0*dd)
        d2p=d2+dd
        d2m=d2-dd
        ker_d1_d2p=GET_KER_VAL(d1,d2p,IVDW_NL,gamma_vdw,alpha_vdw)
        ker_d1_d2m=GET_KER_VAL(d1,d2m,IVDW_NL,gamma_vdw,alpha_vdw)
        phi_der_d2(id1,id2)=(ker_d1_d2p -ker_d1_d2m )/(2.d0*dd)
        !symmetrise
        ker_d1p_d2p=GET_KER_VAL(d1p,d2p,IVDW_NL,gamma_vdw,alpha_vdw)
        ker_d1m_d2p=GET_KER_VAL(d1m,d2p,IVDW_NL,gamma_vdw,alpha_vdw)
        ker_d1p_d2m=GET_KER_VAL(d1p,d2m,IVDW_NL,gamma_vdw,alpha_vdw)
        ker_d1m_d2m=GET_KER_VAL(d1m,d2m,IVDW_NL,gamma_vdw,alpha_vdw)
        phi_der_d1d2(id1,id2)=(ker_d1p_d2p+ker_d1m_d2m-ker_d1m_d2p-ker_d1p_d2m)/(4.d0*dd*dd)
      ENDIF
    ENDDO

#ifdef MPI
    CALL M_sum_d(GRIDC%COMM, phi, nd*nd)
    CALL M_sum_d(GRIDC%COMM, phi_der_d1, nd*nd)
    CALL M_sum_d(GRIDC%COMM, phi_der_d2, nd*nd)
    CALL M_sum_d(GRIDC%COMM, phi_der_d1d2, nd*nd)
#endif

    !The phi-arrays are symmetric
    DO id1=1,nd
      DO id2=1,id1-1
        phi(id2,id1)= phi(id1,id2)
        phi_der_d2(id2,id1)= phi_der_d1(id1,id2)
        phi_der_d1(id2,id1)= phi_der_d2(id1,id2)
        phi_der_d1d2(id2,id1)= phi_der_d1d2(id1,id2)
      ENDDO
    ENDDO

    IF (IU0>=0) WRITE(IU0,*) 'The calculation of the vdW kernel is finished'
    !zero the border components
    !this actually should be avoided when the asymptotic form of the kernel
    !is included in PHI_INTERPOLATE, for the "nd" parts of the arrays
    IF (.NOT.luse_asymptotic) THEN
      phi(:,nd)=0._q
      phi(nd,:)=0._q
      phi_der_d1(:,nd)=0._q
      phi_der_d1(nd,:)=0._q
      phi_der_d2(:,nd)=0._q
      phi_der_d2(nd,:)=0._q
      phi_der_d1d2(:,nd)=0._q
      phi_der_d1d2(nd,:)=0._q
    ENDIF
    phi(1,nd)=0._q
    phi(nd,1)=0._q
    phi_der_d1(1,nd)=0._q
    phi_der_d1(nd,1)=0._q
    phi_der_d2(1,nd)=0._q
    phi_der_d2(nd,1)=0._q
    phi_der_d1d2(1,nd)=0._q
    phi_der_d1d2(nd,1)=0._q

    phi_der_d1(1,:)=0._q
    phi_der_d2(:,1)=0._q
    phi_der_d1d2(:,1)=0._q
    phi_der_d1d2(1,:)=0._q

    DO id1=1,nd-1
      DO id2 = 1,id1
        by(:) =(/phi(id1,id2),phi(id1+1,id2),phi(id1+1,id2+1),phi(id1,id2+1) /)
        by_d1(:) =(/phi_der_d1(id1,id2),phi_der_d1(id1+1,id2),phi_der_d1(id1+1,id2+1),phi_der_d1(id1,id2+1) /)
        by_d2(:) =(/phi_der_d2(id1,id2),phi_der_d2(id1+1,id2),phi_der_d2(id1+1,id2+1),phi_der_d2(id1,id2+1) /)
        by_d1d2(:) =(/phi_der_d1d2(id1,id2),phi_der_d1d2(id1+1,id2),phi_der_d1d2(id1+1,id2+1),phi_der_d1d2(id1,id2+1) /)
        CALL BCUCOF(by,by_d1,by_d2,by_d1d2,dd,dd,coefs)
        phi_table(:,:,id1,id2)=coefs
        IF (id1.NE.id2) phi_table(:,:,id2,id1)=coefs
      ENDDO
    ENDDO
    !save the data
    IF (NODE_ME == IONODE) THEN
      open(10, file='vdw_kernel.bindat',form='unformatted',status='replace')
      write(10) IVDW_NL,GAMMA_VDW,ALPHA_VDW,NINTG,nd
      write(10) dmesh
      write(10) phi_table
      close(10)
    ENDIF

  END SUBROUTINE GET_KERNEL_TABLE

    !> This function integrates the Dion et al kernel using NR p 130. remark
    !> the point is to get the phi in some sensible points, fit with spline
    !> and then use the integrated 3.3.3[NR].
  FUNCTION GET_KER_INT(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw)
__GPU USE moffload

    REAL (q) :: d1, d2, gamma_vdw, alpha_vdw, get_ker_int
    !we need to make grid for a and b, get the integrand on the grid
    !then fit it with splines and integrate using the pita but nice formula
    REAL (q) :: intgrid(nintg), dina(nintg), dinb(nintg), ders(nintg)
    REAL (q) :: x_work(nintg), u_work(nintg)
    REAL (q) :: a,b,sa,sb,ca,cb,b2,a2,W,T,din
    REAL (q) :: dint1
    INTEGER :: ia,ib,ii,IVDW_NL

GPACC !$ACC ROUTINE(NY) SEQ
GPACC !$ACC ROUTINE(SPLINE) SEQ
GPACC !$ACC ROUTINE(INTSPL) SEQ

__GPU CALL OFFLOAD_ON_PUSH(.TRUE.)

    dint1=maxin/nintg
    DO ii=1,nintg
      intgrid(ii)=(ii-1)*dint1+1.d-13
    ENDDO

GPACC !$ACC ENTER DATA COPYIN(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,dina,intgrid,nintg,dint1) __IF_ASYNC__
GPACC !$ACC ENTER DATA CREATE(dinb,ders,x_work,u_work) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(TO:d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,dina,intgrid,nintg,dint1) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:dinb,ders,x_work,u_work) IF(OFFLOAD_ON)

__CPU !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
__CPU !$OMP SHARED(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,dina,intgrid,nintg,dint1) &
__CPU !$OMP PRIVATE(ia,ib,a,b,sa,sb,ca,cb,a2,b2,W,T,dinb,ders,x_work,u_work)

GPACC !$ACC PARALLEL LOOP PRESENT(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,dina,intgrid,nintg,dint1) &
GPACC !$ACC PRIVATE(ib,a,b,sa,sb,ca,cb,a2,b2,W,T,dinb,ders,x_work,u_work) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,dina,intgrid,nintg,dint1) &
GPOMP !$OMP PRIVATE(ib,a,b,sa,sb,ca,cb,a2,b2,W,T,dinb,ders,x_work,u_work) IF(OFFLOAD_ON)

    DO ia=1,nintg
GPACC !$ACC LOOP SEQ
      DO ib=1,nintg
        a=intgrid(ia)
        b=intgrid(ib)
        sa=sin(a)
        sb=sin(b)
        ca=cos(a)
        cb=cos(b)
        a2=a*a
        b2=b*b
        IF ((ia.GT.0).AND.(ib.GT.0)) THEN
          W = 2._q*((3._q-a2)*b*cb*sa+(3._q-b2)*a*ca*sb+(a2+b2-3._q)*sa*sb-3._q*a*b*ca*cb )/(a*b)
        ELSE IF ((ia.EQ.0).AND.(ib.EQ.0)) THEN
          W = 2._q*((3._q-a2)*b*cb+(3._q-b2)*ca*sb+(a2+b2-3._q)*sb-3._q*b*ca*cb )/(b)
        ELSE IF ((ia.GT.0).AND.(ib.EQ.0)) THEN
          W = 2._q*((3._q-a2)*cb*sa+(3._q-b2)*a*ca+(a2+b2-3._q)*sa-3._q*a*ca*cb )/(a)
        ELSE
          W = 2._q*((3._q-a2)*cb+(3._q-b2)*ca+(a2+b2-3._q)-3._q*ca*cb )
        ENDIF
        T = 0.5_q*(1._q/(NY(a,d1,IVDW_NL,gamma_vdw,alpha_vdw)+NY(b,d1,IVDW_NL,gamma_vdw,alpha_vdw))+1._q/(NY(a,d2,IVDW_NL,gamma_vdw,alpha_vdw) +NY(b,d2,IVDW_NL,gamma_vdw,alpha_vdw))) &
         *(1._q/((NY(a,d1,IVDW_NL,gamma_vdw,alpha_vdw)+NY(a,d2,IVDW_NL,gamma_vdw,alpha_vdw))*(NY(b,d1,IVDW_NL,gamma_vdw,alpha_vdw)+NY(b,d2,IVDW_NL,gamma_vdw,alpha_vdw))) + &
         1._q/((NY(a,d1,IVDW_NL,gamma_vdw,alpha_vdw)+NY(b,d2,IVDW_NL,gamma_vdw,alpha_vdw))*(NY(a,d2,IVDW_NL,gamma_vdw,alpha_vdw)+NY(b,d1,IVDW_NL,gamma_vdw,alpha_vdw))))
        dinb(ib) =  W*T
      ENDDO
      call SPLINE(dint1, dinb, nintg, 0._q, 0._q, ders, x_work, u_work)
      dina(ia)=INTSPL(intgrid, dinb, nintg, ders)
    ENDDO

__CPU !$OMP END PARALLEL DO

GPACC !$ACC EXIT DATA COPYOUT(dina) DELETE(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,intgrid,nintg,dint1,dinb,ders,x_work,u_work) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:dina) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:d1,d2,IVDW_NL,gamma_vdw,alpha_vdw,intgrid,nintg,dint1,dinb,ders,x_work,u_work) IF(OFFLOAD_ON)

__GPU CALL OFFLOAD_ON_POP

    call SPLINE(dint1, dina, nintg, 0._q, 0._q, ders, x_work, u_work)
    get_ker_int=INTSPL(intgrid, dina, nintg, ders)*2._q/pi/pi

  END function GET_KER_INT

  FUNCTION INTSPL(x, y, n, y2)
GPACC !$ACC ROUTINE SEQ
    REAL (q) :: x(n), y(n), y2(n), intspl
    INTEGER :: n, ii
GPOMP !$OMP DECLARE TARGET(INTSPL)
    intspl=0._q
!    DO ii=1,n-1
!      intspl=intspl+0.5_q*(y(ii)+y(ii+1))*(x(ii+1)-x(ii)) -1._q/24._q*(y2(ii)+y2(ii+1))*(x(ii+1)-x(ii))**3.
!    ENDDO
    DO ii=n-1,1,-1
      intspl=intspl+0.5_q*(y(ii)+y(ii+1))*(x(ii+1)-x(ii)) -1._q/24._q*(y2(ii)+y2(ii+1))*(x(ii+1)-x(ii))**3.
    ENDDO


  END FUNCTION INTSPL

  !> If dtot is smaller then dsoft the kernel is gently cut.
  !> This is done by using a function phi_0+dtot**2*phi_2+dtot**4*phi_4.
  !> The phi_0 is a parameter defined above (0.8), dsoft=0.5
  !> so that only phi_2 and phi_4 need to be found.
  !> This can be done by calculating a value and derivative of the kernel for phi(d1',d2')
  !> where d1'/d2'=d1/d2 and sqrt(d1'**2+d2'**2)=dtot, this is done numerically.
  !> The phi_2 and phi_4 formulae are then straightforward
  FUNCTION GET_KER_VAL(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw)
  REAL (q) :: d1,d2,gamma_vdw,alpha_vdw,dtot, phi_2, phi_4, phi_val, phi_der
  REAL (q) :: phi_pl, phi_mi, d1pl, d2pl, d1mi, d2mi, ddelta
  REAL (q) :: get_ker_val
  INTEGER :: IVDW_NL

  ddelta = dsoft*0.02_q
  !some sensible value
  dtot = sqrt(d1*d1+d2*d2)
  IF ((dtot .LE. 0._q).OR.(d1.LT.0._q).OR.(d2.LT.0._q)) THEN
    get_ker_val=phi_0
    RETURN
  ELSE IF (dtot .GE. dsoft) THEN
    get_ker_val=GET_KER_INT(d1,d2,IVDW_NL,gamma_vdw,alpha_vdw)
    RETURN
  ELSE
    d1pl = d1*(dsoft+ddelta)/dtot
    d2pl = d2*(dsoft+ddelta)/dtot
    d1mi = d1*(dsoft-ddelta)/dtot
    d2mi = d2*(dsoft-ddelta)/dtot
    phi_pl = GET_KER_INT(d1pl,d2pl,IVDW_NL,gamma_vdw,alpha_vdw)
    phi_mi = GET_KER_INT(d1mi,d2mi,IVDW_NL,gamma_vdw,alpha_vdw)
    phi_val = 0.5_q*(phi_pl+phi_mi)
    phi_der = (phi_pl-phi_mi)/(2._q*ddelta)
    phi_2 = (4._q*(phi_val-phi_0)-dsoft*phi_der)/(2._q*dsoft*dsoft)
    phi_4 = (phi_0 - phi_val+0.5_q*dsoft*phi_der)/(dsoft**4.)
    get_ker_val=phi_0+phi_2*dtot*dtot +phi_4*dtot**4.
  ENDIF
  END function GET_KER_VAL

  FUNCTION NY(ar,d,IVDW_NL,gamma_vdw,alpha_vdw)
    REAL(q) :: ar,d,gamma_vdw,alpha_vdw,ny,y,y2,y4,gamm2
    INTEGER :: IVDW_NL
GPOMP !$OMP DECLARE TARGET(NY)
    IF ((ar<1.e-10_q).AND.(d>0._q)) THEN
      ny = 0.5_q*d*d/gamma_vdw
    ELSEIF ((ar>0._q).AND.(d>0._q)) THEN
      IF (IVDW_NL==1) THEN
        y=ar/d
        ny = ar*ar*0.5_q/(1._q-exp(-gamma_vdw*y*y))
      ELSEIF ((IVDW_NL==3).OR.(IVDW_NL==4)) THEN
        gamm2=gamma_vdw*gamma_vdw
        y=ar/d
        y2=y*y
        y4=y2*y2
        ny = ar*ar*0.5_q/( 1._q - 1._q / ( 1._q + gamma_vdw*y2 + gamm2*y4 + alpha_vdw*y4*y4 ) )
      ENDIF
    ELSE
      ny = 0.5_q*ar*ar
    ENDIF
    RETURN
  END function NY

  !> Bicubic spline from NR.
  SUBROUTINE BCUCOF(y,y1,y2,y12,d1,d2,c)
    IMPLICIT NONE
    REAL(q), INTENT(IN) :: d1,d2
    REAL(q), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
    REAL(q), DIMENSION(4,4), INTENT(OUT) :: c
    REAL(q), DIMENSION(16) :: x
    REAL(q), DIMENSION(16,16) :: wt
    DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
    8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
    2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
    2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
    -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
    -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
    -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
    x(1:4)=y
    x(5:8)=y1*d1
    x(9:12)=y2*d2
    x(13:16)=y12*d1*d2
    x=matmul(wt,x)
    c=reshape(x,(/4,4/),order=(/2,1/))
  END SUBROUTINE BCUCOF

    !> This performs the FFT of the kernel which is spherically symmetric
    !> and thus one just needs to perform sin(kr) transformation in principle
    !> times the $k^{-1} r$ factors.
    !> The sin transforms are used in setlocalpp.F so I guess I copy FOUR1 from there
    !> the radial array is [2^N] anyway.
    !> We need to know Nr, rcut.
  SUBROUTINE PHIK_FFT(phir_work, phik_work)
    REAL (q) :: phir_work(0:nr), phik_work(0:nr)
    REAL (q) :: work(2,0:2*nr)
    REAL (q) :: dr, dk, fac, kval, rval
    INTEGER :: ik, ir, ii
    dr=rcut/Nr
    fac=dr/sqrt(2._q*pi)
    dk=pi/rcut
    work=0._q

    do ii=1, 2*Nr-1
      !this can be loop over half but I'm tired
      if (ii.lt.Nr) then
        ir = ii
        rval = dr*ir
        work(2,ii) = -fac*phir_work(ii)*rval
      else if (ii.eq.Nr) then
        work(2,ii) = 0._q
      else
        ir = 2*Nr-ii
        rval = -dr*ir
        work(2,ii) = -fac*phir_work(ir)*rval
      endif
    enddo

    call FOUR1(work, 2*Nr, 1)
    do ik=1, Nr
      kval=ik*dk
      phik_work(ik)=work(1,ik)/kval
    enddo

    !calculate phik_work(0)
    phik_work(0)=0._q
    do ii=1,Nr
      rval=dr*ii
      phik_work(0)=phik_work(0)+rval*rval*phir_work(ii)
    enddo
    phik_work(0)=phik_work(0)*2._q*fac

  END SUBROUTINE PHIK_FFT

  SUBROUTINE FOUR1(DATA,NN,ISIGN)
      USE prec
      IMPLICIT NONE
      REAL(q) WR,WI,WPR,WPI,WTEMP,THETA,TEMPI,TEMPR
      REAL(q) DATA(*)
      INTEGER M,ISTEP,MMAX,ISIGN,NN,N,I,J
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959_q/(ISIGN*MMAX)
        WPR=-2._q*DSIN(0.5_q*THETA)**2
        WPI=DSIN(THETA)
        WR=1._q
        WI=0._q
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
  END SUBROUTINE FOUR1

    !> cut off function for the kernels,
    !> personally this looks quite harsh.
  FUNCTION S_CUT(val)
    REAL(q) :: val, s_cut
    if (val .le. 0._q) then
      s_cut=1._q
    else if (val .ge. 1) then
      s_cut=0._q
    else
      s_cut=(1._q-val**8)**4
    endif

  END FUNCTION S_CUT

      !> From Num recipes.
  SUBROUTINE SPLINE(dx,y,n,yp1,ypn,y2,x,u)
GPACC !$ACC ROUTINE SEQ
      INTEGER n
      REAL(q) :: yp1,ypn,x(n),y(n),y2(n),dx
      INTEGER i,k
      REAL(q) :: p,qm,sig,un,u(n)
GPOMP !$OMP DECLARE TARGET(SPLINE)

      do i=1,n
       x(i)=dx*i
      enddo
      if (yp1.gt..99e30_q) then
        y2(1)=0._q
        u(1)=0._q
      else
        y2(1)=-0.5_q
        u(1)=(3._q/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do  i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2._q
        y2(i)=(sig-1._q)/p
        u(i)=(6._q*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30_q) then
        qm=0._q
        un=0._q
      else
        qm=0.5_q
        un=(3._q/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qm*u(n-1))/(qm*y2(n-1)+1._q)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
  END  SUBROUTINE SPLINE


  SUBROUTINE MIN_K(GRIDC,LATT_CUR)
    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    INTEGER :: I, N1, N2, N3, NC, ix, iy, iz
    REAL(q) :: cell(3,3)
    REAL(q) :: k_x,k_y,k_z,kte
    !the latt_cur%B has rec vector x stored as 11 21 31
    !and similarly for the rest
    !to get total coordinate x in k-space sum is done over
    !11 12 13
    N1=GRIDC%NGX
    N2=GRIDC%NGY
    N3=GRIDC%NGZ
    kcut = 10.e10_q
    cell=latt_cur%B*AUTOA
    !now we want to know the smallest boundary k-vector
    !IF (IU0>=0) WRITE(IU0,*)  'reciprocal cell'
    !IF (IU0>=0) WRITE(IU0,*)  cell
    !IF (IU0>=0) WRITE(IU0,*)  N1*cell(1,1), N2*cell(1,2), N3*cell(1,3)
    !IF (IU0>=0) WRITE(IU0,*) N1*cell(2,1),N2*cell(2,2), N3*cell(2,3)
    !IF (IU0>=0) WRITE(IU0,*) N1*cell(3,1), N2*cell(3,2), N3*cell(3,3)

    do ix = -1,1
      do iy = -1,1
        do iz = -1,1
          if ((abs(ix)+abs(iy)+abs(iz)).ne.0 ) then
            k_x=N1*cell(1,1)*ix+N2*cell(1,2)*iy+N3*cell(1,3)*iz
            k_y=N1*cell(2,1)*ix+N2*cell(2,2)*iy+N3*cell(2,3)*iz
            k_z=N1*cell(3,1)*ix+N2*cell(3,2)*iy+N3*cell(3,3)*iz
            kte=2._q*pi*sqrt(k_x*k_x+k_y*k_y+k_z*k_z)*0.5_q
            kcut=min(kcut, kte)
          endif
        enddo
      enddo
    enddo

  END SUBROUTINE MIN_K


  SUBROUTINE GET_K_IND(k_k,GRIDC,LATT_CUR)

__GPU USE moffload_struct_def

    !get the k-space interpolation indices
    !can be done once
    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    REAL(q) k_k(GRIDC%MPLWV)
    INTEGER :: I, N1, N2, N3, NC
    REAL(q) :: cell(3,3)
    REAL(q) :: k_x, k_y, k_z

    cell=latt_cur%B*AUTOA
GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,K_K) COPYIN(CELL) PRIVATE(N1,NC,N2,N3,K_X,K_Y,K_Z) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,K_K) MAP(TO:CELL) &
GPOMP !$OMP PRIVATE(N1,NC,N2,N3,K_X,K_Y,K_Z) IF(OFFLOAD_ON)

    do i=1,gridc%rc%np
      n1=Mod((i-1),gridc%rc%nrow)+1
      nc=(i-1)/gridc%rc%nrow+1
      n2=gridc%rc%i2(nc)
      n3=gridc%rc%i3(nc)
      k_x=gridc%lpctx(n1)*cell(1,1)+gridc%lpcty(n2)*cell(1,2)+gridc%lpctz(n3)*cell(1,3)
      k_y=gridc%lpctx(n1)*cell(2,1)+gridc%lpcty(n2)*cell(2,2)+gridc%lpctz(n3)*cell(2,3)
      k_z=gridc%lpctx(n1)*cell(3,1)+gridc%lpcty(n2)*cell(3,2)+gridc%lpctz(n3)*cell(3,3)
      k_k(i)=2._q*pi*sqrt(k_x*k_x+k_y*k_y+k_z*k_z)
    enddo

  END SUBROUTINE GET_K_IND


END MODULE VDW_LL



SUBROUTINE VDW_CALC(DCHARG,DWORK,DWORKG,DWORK4,EXC,GRIDC,LATT_CUR,CINP,DINP,spl,phik,d2phidk2,qs_vdw,k_k,stress,ISIF,E_corr,ISPIN,Zab_vdW)

__GPU USE moffload

  USE prec
  USE lattice
  USE mpimy
  USE mgrid
  USE constant

GPACC !$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

  IMPLICIT NONE
  INTEGER NP
!I guess LATT_CURR is the one that we need not OMEGA only
  TYPE (latt) LATT_CUR
  TYPE (grid_3d) GRIDC
  REAL(q) OMEGA
  REAL(q) Zab_vdW
  REAL(q) EXC
  REAL(q) DCHARG(GRIDC%RL%NP) !charge density on input and output
  REAL(q) DWORKG(GRIDC%RL%NP) !df/drho input, add the same non-local here
                    !don't touch this one yet
  REAL(q) DWORK4(GRIDC%RL%NP) !gradient on input
                    !df/d(drho)/(drho)potential on output
  REAL(q) stress(3,3), E_corr
  INTEGER :: I, IU, alf, bet, ik_k, iofq, ii, jj, ISIF, ISPIN
  INTEGER :: NODE_ME, IONODE, IDUMP
  INTEGER,parameter :: Nalp=30, Nr=4096
  COMPLEX(q) :: CINP(GRIDC%MPLWV)
  COMPLEX(q), allocatable::   THETA(:,:),f_k(:)

  RGRID :: DINP(DIMREAL(GRIDC%MPLWV))
  REAL(q),parameter :: lda_a=0.0310907_q, lda_a1=0.21370_q, lda_b1=7.5957_q
  REAL(q),parameter :: lda_b2=3.5876_q, lda_b3=1.6382_q, lda_b4=0.49294_q
  !REAL(q),parameter :: Zab=-0.8491_q !stored in setexm.inc to use vdw_df2 easily
  REAL(q), parameter :: rcut=100._q
  REAL(q),parameter :: dk_k=pi/rcut
  REAL(q),parameter :: soft_cor=-0.0000000000000_q
  REAL(q) :: s, nrs, rnrs, ex_LDA, ec_LDA, vc_LDA
  REAL(q) :: diff, coef(4)
  REAL(q),parameter :: q0cut=10.0_q, lambda=1.2_q
  REAL(q) :: spl(Nalp,Nalp,4)
  REAL(q) :: phik(0:Nr,Nalp,Nalp), d2phidk2(0:Nr,Nalp,Nalp)
  REAL(q) , allocatable :: q0(:),dq0(:),pa(:),f_g(:,:)
  REAL(q) :: qs_vdw(Nalp)
  REAL(q) :: k_k(GRIDC%MPLWV)
  REAL(q),allocatable :: dj_k(:)
  REAL(q) :: energ,soft_pot,soft_der
  REAL(q) :: a_k, b_k, a2_k, b2_k, a3_k, b3_k, pref
  REAL(q) :: ro,gr,dq0dn, dq0dgr, dpadq0, dtheadn, dtheadgr
  REAL(q) :: dhdx, pa_tm, vld1, vld2, vld3, vld4, phi, dphidk
  REAL(q) :: k_x(3), k_y, k_z, k_a, cell(3,3)
  REAL(q) :: max_n, max_gr, a_iofq, b_iofq
  REAL(q) :: RINPL, THRD
  REAL(q), allocatable :: k_x_pre(:,:)
  REAL(q) :: str11,str12,str13,str21,str22,str23,str31,str32,str33
  RGRID   DWORK(GRIDC%RL%NP)
  INTEGER,allocatable :: j_k(:),alp(:)
  INTEGER :: n1,n2,n3,nc
  INTEGER ::  i_max_n, i_max_gr

  NP=GRIDC%RL%NP
  !  IF (IU0>=0) WRITE(IU0,*)  "MPLWV: ", GRIDC%MPLWV, DIMREAL(GRIDC%MPLWV), GRIDC%RC%NP
  ALLOCATE(q0(NP),dq0(NP),alp(NP),pa(NP),f_g(NP,Nalp))
  allocate(THETA(GRIDC%MPLWV,Nalp),f_k(GRIDC%MPLWV))
  allocate(dj_k(GRIDC%MPLWV),j_k(GRIDC%MPLWV))
  ALLOCATE(k_x_pre(4,GRIDC%RC%NP))
GPACC !$ACC ENTER DATA CREATE(q0,dq0,alp,pa,f_g,theta,f_k,dj_k,j_k) __IF_ASYNC__
GPACC !$ACC ENTER DATA CREATE(k_x_pre) __IF_ASYNC__

GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:Q0,DQ0,ALP,PA,F_G,THETA,F_K,DJ_K,J_K) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:K_X_PRE) IF(OFFLOAD_ON)

  cell=latt_cur%B*AUTOA
GPACC !$ACC ENTER DATA COPYIN(qs_vdw,spl,phik,d2phidk2,cell) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:QS_VDW,SPL,PHIK,D2PHIDK2,CELL) IF(OFFLOAD_ON)
  str11=0._q
  str21=0._q
  str31=0._q
  str12=0._q
  str22=0._q
  str32=0._q
  str13=0._q
  str23=0._q
  str33=0._q

  RINPL = 1._q/(real(GRIDC%NGX,kind=q)*GRIDC%NGY*GRIDC%NGZ)
  THRD = 1._q/3._q

  energ=0.0_q

#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
      IDUMP=0
#ifdef debugNODE_ME
      IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
      NODE_ME=0
      IONODE =0
      IDUMP=0
#ifdef debug
      IDUMP=1
#endif
#endif
  !calculate q0 in each point and it's saturated version (Soler's paper)
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(DCHARG,DWORK4,q0) &
GPACC !$ACC& PRIVATE(ro,gr,s,nrs,rnrs,ex_LDA,ec_LDA) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:DCHARG,DWORK4,Q0) &
GPOMP !$OMP PRIVATE(RO,GR,S,NRS,RNRS,EX_LDA,EC_LDA) IF(OFFLOAD_ON)

  DO I=1,NP
    IF (DCHARG(i).LT. 1e-11_q) THEN
      ro=1e-15_q
      gr=0._q
    ELSE
      ro=DCHARG(i)
      gr=DWORK4(i)
    ENDIF
    s=gr*AUTOA4*0.5_q/ro/(3._q*pi*pi*ro)**THRD/(AUTOA4)
    nrs=(3._q/(4._q*pi*ro*AUTOA3))**THRD
    rnrs=sqrt(nrs)
    ex_LDA=-(3._q*0.25_q/pi)*(9._q*pi*0.25_q)**THRD/nrs
    ec_LDA=-2.0_q*lda_a*(1._q+lda_a1*nrs)*log(1._q+1._q/(2._q*lda_a*(lda_b1*rnrs+lda_b2*nrs+lda_b3*rnrs*nrs+lda_b4*nrs*nrs)))
    q0(i)=(3._q*pi*pi*ro*AUTOA3)**THRD*(1._q+ec_LDA/ex_LDA-Zab_VDW/9._q*s*s)
    !now saturated
    nrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO iu=1,12
      nrs=nrs+(q0(i)/q0cut)**iu/real(iu,kind=q)
    ENDDO
    q0(i)=q0cut*(1._q-exp(-nrs))
  ENDDO



  !calculate alp and dq0 in each point
  !that is which interpolating point we will be using and how far from it are we
  a_iofq=log((qs_vdw(Nalp)-qs_vdw(Nalp-1))/(qs_vdw(2)-qs_vdw(1)))/(Nalp-2.0_q)
  a_iofq=max(a_iofq,1e-12_q)
  b_iofq=(qs_vdw(2)-qs_vdw(1))/(exp(a_iofq)-1.0_q)
GPACC !$ACC PARALLEL PRESENT(q0,alp,dq0,qs_vdw) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(iofq)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:Q0,ALP,QS_VDW) PRIVATE(IOFQ) IF(OFFLOAD_ON)

  DO I=1,NP
    iofq=1+log(1+(q0(i)-qs_vdw(1))/b_iofq)/a_iofq
    iofq=max(1,iofq)
    iofq=min(Nalp,iofq)
    IF(q0(i).GT.(q0cut-0.00015_q)) iofq=Nalp
    alp(i)=iofq
  ENDDO

  !this is the difference from interpolating point
GPACC !$ACC LOOP GANG VECTOR
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:Q0,DQ0,ALP,QS_VDW) IF(OFFLOAD_ON)
  DO I=1,NP
    dq0(i)=q0(i)-qs_vdw(alp(i))
  ENDDO
GPACC !$ACC END PARALLEL

  !for each interpolating point create p_\alpha (Soler's paper)
  !and Fourier transform it, this is stored in \Theta array
  DO alf=1,Nalp
GPACC !$ACC PARALLEL PRESENT(spl,alp,dq0,pa,DCHARG,dinp) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(coef,diff)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:SPL,ALP,DQ0,PA) PRIVATE(coef,diff) IF(OFFLOAD_ON)
    DO I=1,NP
      coef(:)=spl(alp(i),alf,:)
      diff=dq0(i)
      pa(i)=coef(1)+diff*(coef(2)+diff*(coef(3)+diff*coef(4)))
    ENDDO
GPACC !$ACC LOOP GANG VECTOR
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DCHARG,DINP,PA) IF(OFFLOAD_ON)
    DO I=1,NP
      IF (DCHARG(i).ge.0._q) THEN
        dinp(i)=pa(i)*DCHARG(i)*AUTOA3
      ELSE
        dinp(i)=0._q
      ENDIF
    ENDDO
GPACC !$ACC END PARALLEL

    CALL OPSYNC(DINP,CINP,GRIDC%NPLWV)
    CALL FFT3D(CINP,GRIDC,-1)

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(THETA,CINP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:THETA,CINP) IF(OFFLOAD_ON)
    DO I=1,GRIDC%RC%NP
      THETA(i,alf)=CINP(i)
    ENDDO
  ENDDO

  !find values of the kernel on the reciprocal grid
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(dj_k,k_k,j_k) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DJ_K,K_K,J_K) IF(OFFLOAD_ON)

  DO i=1,GRIDC%RC%NP
    dj_k(i)=0.5_q*k_k(i)/pi*rcut
    j_k(i)=floor(dj_k(i))
    dj_k(i)=dj_k(i)-j_k(i)
    dj_k(i)=dj_k(i)*tpi/rcut
  ENDDO

  IF (ISIF .ge. 1) THEN
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(k_x_pre,gridc,cell) PRIVATE(n1,nc,n2,n3) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:K_X_PRE,GRIDC,CELL) &
GPOMP !$OMP PRIVATE(N1,NC,N2,N3) IF(OFFLOAD_ON)

     DO i=1,GRIDC%RC%NP
        n1=Mod((i-1),gridc%rc%nrow)+1
        nc=(i-1)/gridc%rc%nrow+1
        n2=gridc%rc%i2(nc)
        n3=gridc%rc%i3(nc)
        k_x_pre(2,i)=tpi*(gridc%lpctx(n1)*cell(1,1)+gridc%lpcty(n2)*cell(1,2)+gridc%lpctz(n3)*cell(1,3))
        k_x_pre(3,i)=tpi*(gridc%lpctx(n1)*cell(2,1)+gridc%lpcty(n2)*cell(2,2)+gridc%lpctz(n3)*cell(2,3))
        k_x_pre(4,i)=tpi*(gridc%lpctx(n1)*cell(3,1)+gridc%lpcty(n2)*cell(3,2)+gridc%lpctz(n3)*cell(3,3))
        k_x_pre(1,i)=sqrt(k_x_pre(2,i)*k_x_pre(2,i)+k_x_pre(3,i)*k_x_pre(3,i)+k_x_pre(4,i)*k_x_pre(4,i))  !*2._q*pi factor is here
     ENDDO
  ENDIF

  !now we are going to calculate the energy
  !this requires two nested loops over \alpha and \beta
  !first we do the sum over \Theta_\beta * kernel
  DO alf=1,Nalp
GPACC !$ACC KERNELS PRESENT(f_K) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:F_K) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    f_k=0._q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRESENT(k_k,phik,d2phidk2,f_k,theta,gridc,k_x_pre) &
GPACC !$ACC& PRIVATE(ik_k,phi,dphidk,a_k,b_k,a2_k,b2_k,a3_k,b3_k,n1,nc,n3,pref) &
GPACC !$ACC& REDUCTION(+:str11,str12,str13,str21,str22,str23,str31,str32,str33) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) MAP(PRESENT,TO:K_K,PHIK,D2PHIDK2,F_K,THETA,GRIDC,K_X_PRE) &
GPOMP !$OMP PRIVATE(IK_K,PHI,DPHIDK,A_K,B_K,A2_K,B2_K,A3_K,B3_K,N1,NC,N3,PREF) &
GPOMP !$OMP REDUCTION(+:STR11,STR12,STR13,STR21,STR22,STR23,STR31,STR32,STR33) IF(OFFLOAD_ON)

    DO bet=1,Nalp
      DO i=1,GRIDC%RC%NP
        ik_k=k_k(i)/dk_k
        IF (ik_k.GE.Nr) THEN
          !I really don't like the if here
          phi=0._q
          dphidk=0._q
        ELSE
          a_k=((ik_k+1._q)*dk_k-k_k(i))/dk_k
          b_k=1._q-a_k
          a2_k=(3._q*a_k**2.-1._q)*dk_k/6._q
          b2_k=(3._q*b_k**2.-1._q)*dk_k/6._q
          a3_k=(a_k**3.-a_k)*dk_k**2._q/6._q
          b3_k=(b_k**3.-b_k)*dk_k**2._q/6._q
          phi=a_k*phik(ik_k, bet, alf)+b_k*phik(ik_k+1,bet,alf) &
           +a3_k*d2phidk2(ik_k,bet,alf)+b3_k*d2phidk2(ik_k+1,bet,alf)
          dphidk=(-phik(ik_k, bet, alf)+phik(ik_k+1,bet,alf))/dk_k &
                 -a2_k*d2phidk2(ik_k,bet,alf)+b2_k*d2phidk2(ik_k+1,bet,alf)

        ENDIF
#ifdef OFFLOADING
        CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(f_k(i),theta(i,bet)*phi)
#else // OFFLOADING
        f_k(i)=f_k(i)+theta(i,bet)*phi
#endif // OFFLOADING
        IF (ISIF .ge. 1) THEN
        if (k_x_pre(1,i) .gt. 1.0e-13_q) then
          pref=dphidk/k_x_pre(1,i)*theta(i,bet)*conjg(theta(i,alf))
#ifdef NGZhalf
          nc=(i-1)/gridc%rc%nrow+1
          n3=gridc%rc%i3(nc)
          if ((n3.ne.1).and.(n3.ne.(GRIDC%NGZ/2+1))) pref=2._q*pref
#elif defined(NGXhalf)
          n1=Mod((i-1),gridc%rc%nrow)+1
          if ((n1.ne.1).and.(n1.ne.(GRIDC%NGX/2+1))) pref=2._q*pref
#endif
          str11=str11- pref*k_x_pre( 1+1,i) * k_x_pre( 1+1,i)
          str12=str12- pref*k_x_pre( 1+1,i) * k_x_pre( 2+1,i)
          str13=str13- pref*k_x_pre( 1+1,i) * k_x_pre( 3+1,i)
          str21=str21- pref*k_x_pre( 2+1,i) * k_x_pre( 1+1,i)
          str22=str22- pref*k_x_pre( 2+1,i) * k_x_pre( 2+1,i)
          str23=str23- pref*k_x_pre( 2+1,i) * k_x_pre( 3+1,i)
          str31=str31- pref*k_x_pre( 3+1,i) * k_x_pre( 1+1,i)
          str32=str32- pref*k_x_pre( 3+1,i) * k_x_pre( 2+1,i)
          str33=str33- pref*k_x_pre( 3+1,i) * k_x_pre( 3+1,i)
        endif
        ENDIF !only calculate stress when required

      ENDDO !q points
    ENDDO  !beta

GPACC !$ACC PARALLEL PRESENT(gridc,f_k,theta,CINP) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(nc,n1,n3,pref) REDUCTION(+:energ)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,F_K,THETA) &
GPOMP !$OMP PRIVATE(NC,N1,N3,PREF) REDUCTION(+:ENERG) IF(OFFLOAD_ON)

    DO i=1,GRIDC%RC%NP
      pref=real(f_k(i)*conjg(theta(i,alf)),kind=q)
#ifdef NGZhalf
!for MPI version reduced in Z direction
!and then we sum it with \Theta_\alpha
      nc=(i-1)/gridc%rc%nrow+1
      n3=gridc%rc%i3(nc)
      if ((n3.ne.1).and.(n3.ne.(GRIDC%NGZ/2+1))) pref=2._q*pref
#else
#ifdef NGXhalf
!for serial version reduced in X direction
      n1=Mod((i-1),gridc%rc%nrow)+1
      if ((n1.ne.1).and.(n1.ne.(GRIDC%NGX/2+1))) pref=2._q*pref
#endif
#endif
      energ=energ+pref
    ENDDO

    !now we fourier transform the \Theta * kernel to real space to reuse it for potential
GPACC !$ACC LOOP GANG VECTOR
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:F_K,CINP) IF(OFFLOAD_ON)
    DO i=1,GRIDC%RC%NP
      CINP(I)=f_k(i)
    ENDDO
GPACC !$ACC END PARALLEL

    CALL OPSYNC(DINP,CINP,GRIDC%NPLWV)
    CALL FFT3D(CINP,GRIDC,1)

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(dinp,f_g) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DINP,F_G) IF(OFFLOAD_ON)
    DO i=1,GRIDC%RL%NP
      f_g(i,alf)=dinp(i)*RINPL
    ENDDO
  ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
  !add the energy to the xc energy, include volume factor
  !  number of grid points since we've done FFT, and factors to get eV
  EXC=EXC + 0.5_q*energ*LATT_CUR%OMEGA*RINPL/AUTOA3*2._q*RYTOEV
  !compared to vander (my local version (jik)),
  !we divide by (M)PLWV only once since it's done in xcgrad for the second time
  !so we can say that we get rid of the FFT 1/N factor and keep the \Omega/N for xcgrad
  E_corr=0.5_q*energ*LATT_CUR%OMEGA*RINPL**2./AUTOA3*2._q*RYTOEV

  max_n=0._q
  max_gr=0._q

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(DCHARG,DWORK4,q0,spl,alp,dq0,f_g,DWORKG,DWORK) &
GPACC !$ACC& PRIVATE(ro,gr,nrs,rnrs,vld1,vld2,vld3,vld4,ec_LDA,vc_LDA,dq0dn,dq0dgr,dhdx, &
GPACC !$ACC& dtheadn,dtheadgr,coef,diff,pa_tm,dpadq0) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:DCHARG,DWORK4,Q0,SPL,ALP,DQ0,F_G,DWORKG,DWORK) &
GPOMP !$OMP PRIVATE(RO,GR,NRS,RNRS,VLD1,VLD2,VLD3,VLD4,EC_LDA,VC_LDA,DQ0DN,DQ0DGR,DHDX, &
GPOMP !$OMP DTHEADN,DTHEADGR,COEF,DIFF,PA_TM,DPADQ0) IF(OFFLOAD_ON)

  DO i=1,gridc%RL%NP
    !calc the potential
    !here we calculate the lda_c potential, and derivatives of q0 wrt density and density gradient
    IF (DCHARG(i).LT. 1e-11_q) THEN
      ro=1e-15_q
      gr=0._q
    ELSE
      ro=DCHARG(i)*AUTOA3
      gr=DWORK4(i)*AUTOA4
    ENDIF
    nrs=(3._q/(4._q*pi*ro))**THRD
    rnrs=sqrt(nrs)
    vld1=lda_b1*rnrs+lda_b2*nrs+lda_b3*rnrs*nrs+lda_b4*nrs*nrs
    vld2=2._q*lda_a*vld1
    vld3=0.5_q*lda_b1/rnrs+lda_b2+1.5_q*lda_b3*rnrs+2._q*lda_b4*nrs
    vld4=log(1._q+1._q/vld2)
    ec_LDA=-2.0_q*lda_a*(1._q+lda_a1*nrs)*vld4
    vc_LDA=-2.0_q*lda_a*lda_a1*vld4+2.0_q*lda_a*(1._q+lda_a1*nrs)*vld3/(vld1*(vld2+1._q))
    vc_LDA=(ec_LDA+ro*vc_LDA*(-THRD)*(3._q/(4._q*pi))**THRD*ro**(-4._q*THRD))
    ec_LDA=ec_LDA*ro
    dq0dn=(pi/(3._q*ro))**(2._q*THRD)+4._q*pi*THRD*(ec_LDA/ro-vc_LDA)/ro+7._q/108._q*Zab_VDW*gr*gr*ro**(-10._q/3._q)/(pi*pi*3._q)**THRD
    dq0dgr=-Zab_VDW/(18._q*(3._q*pi*pi*ro)**THRD)*gr/(ro*ro)
    nrs=0._q
    rnrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:RNRS) IF(OFFLOAD_ON)
    DO iu=1,12
      rnrs=rnrs+(q0(i)/q0cut)**iu/real(iu,kind=q)
    ENDDO
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO iu=0,11
      nrs=nrs+(q0(i)/q0cut)**iu
    ENDDO
    dhdx=nrs*exp(-rnrs)
    dtheadn=0._q
    dtheadgr=0._q
    !one loop is needed to calculate the potential
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:DTHEADN,DTHEADGR) IF(OFFLOAD_ON)
    DO alf=1,Nalp
      coef=spl(alp(i),alf,:)
      diff=dq0(i)
      pa_tm=coef(1)+diff*(coef(2)+diff*(coef(3)+diff*coef(4)))
      dpadq0=coef(2)+diff*(2._q*coef(3)+3._q*diff*coef(4))
      dtheadn=dtheadn+(pa_tm+ro*dpadq0*dq0dn*dhdx)*f_g(i,alf)
      dtheadgr=dtheadgr+ro*dpadq0*dhdx*dq0dgr*f_g(i,alf)
    ENDDO

    !add the potentials to the GGA ones
    !this is d/dn
    DWORKG(i)=DWORKG(i)+dtheadn*2._q*RYTOEV
    !this is d/ddn
    IF (ISPIN==1) THEN
       DWORK(i)=DWORK(i)+dtheadgr*2._q*RYTOEV/MAX(DWORK4(I),1.E-10_q)*AUTOA
    ELSE
       DWORK(i)=DWORK(i)+dtheadgr*2._q*RYTOEV*AUTOA
    ENDIF
  ENDDO !over points on real space grid

GPACC !$ACC EXIT DATA DELETE(qs_vdw,spl,phik,d2phidk2,cell) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(q0,dq0,alp,pa,f_g,theta,f_k,dj_k,j_k) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(k_x_pre) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:QS_VDW,SPL,PHIK,D2PHIDK2,CELL) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:Q0,DQ0,ALP,PA,F_G,THETA,F_K,DJ_K,J_K) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:K_X_PRE) IF(OFFLOAD_ON)

  DEALLOCATE(q0,dq0,alp,pa)
  DEALLOCATE(dj_k,theta,f_k,f_g,j_k)
  DEALLOcATE(k_x_pre)
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
  stress(1,1)=str11
  stress(2,1)=str21
  stress(3,1)=str31
  stress(1,2)=str12
  stress(2,2)=str22
  stress(3,2)=str32
  stress(1,3)=str13
  stress(2,3)=str23
  stress(3,3)=str33
  stress=0.5_q*stress*LATT_CUR%OMEGA*RINPL**2.*2._q*RYTOEV/AUTOA3


#ifdef MPI
 CALLMPI( M_sum_d(GRIDC%COMM,e_corr,1))
#endif

! IF (NODE_ME == IONODE) THEN
!   WRITE(*,'(A,F14.7)') 'Total vdW correction in eV: ', e_corr
! ENDIF

END SUBROUTINE VDW_CALC

SUBROUTINE VDW_CALC_spin(DCHARG,DCHARGS,DWORK,DWORKS,DWORKG,DWORKGS,DWORK4,DWORK4S, &
                       & EXC,GRIDC,LATT_CUR,CINP,DINP,spl,phik,d2phidk2,qs_vdw, & 
                       & k_k,stress,ISIF,E_corr,Zab_vdW)

__GPU USE moffload

  USE prec
  USE lattice
  USE mpimy
  USE mgrid
  USE constant
  USE ldalib, only: EX, VX, CORPBE_LDA

GPACC !$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

  IMPLICIT NONE
  INTEGER NP
  TYPE (latt) LATT_CUR
  TYPE (grid_3d) GRIDC
  REAL(q) OMEGA
  REAL(q) Zab_vdW
  REAL(q) EXC
  !spin vdw mod these arrays have two instances, one for spin up, other for spin down
  REAL(q) DCHARG(GRIDC%RL%NP), DCHARGS(GRIDC%RL%NP) !charge density on input and output
  REAL(q) DWORKG(GRIDC%RL%NP), DWORKGS(GRIDC%RL%NP) !df/drho input, add the same non-local here
  REAL(q) DWORK4(GRIDC%RL%NP), DWORK4S(GRIDC%RL%NP) !gradient on input
                    !df/d(drho)/(drho)potential on output
  RGRID   DWORK(GRIDC%RL%NP), DWORKS(GRIDC%RL%NP)
  REAL(q) stress(3,3),E_corr
  INTEGER :: I, IU, alf, bet, ik_k, iofq, ii, jj, ISIF
  INTEGER :: NODE_ME, IONODE, IDUMP
  INTEGER,parameter :: Nalp=30, Nr=4096
  COMPLEX(q) :: CINP(GRIDC%MPLWV)
  COMPLEX(q), allocatable::   THETA(:,:),f_k(:)

  RGRID :: DINP(DIMREAL(GRIDC%MPLWV))
  REAL(q),parameter :: lda_a=0.0310907_q, lda_a1=0.21370_q, lda_b1=7.5957_q
  REAL(q),parameter :: lda_b2=3.5876_q, lda_b3=1.6382_q, lda_b4=0.49294_q
! REAL(q),parameter :: Zab=-0.8491_q !read from INCAR for easy-of-use vdw-df2
  REAL(q),parameter :: dk_k=pi/100._q
  REAL(q),parameter :: soft_cor=-0.0000000000000_q
  REAL(q), parameter :: THRD=1._q/3._q
  !spin vdw mod
  REAL(q) :: s, ss, nrs, zeta, rnrs, ex_LDA, ec_LDA, vc_LDA, qx
  REAL(q) :: KBEX, KBEC, KBVU, KBVD
  REAL(q) :: diff, coef(4)
  REAL(q),parameter :: q0cut=10.0_q, lambda=1.2_q, rcut=100._q
  REAL(q) :: spl(Nalp,Nalp,4)
!  REAL(q) :: phi_ab(Nalp,Nalp,0:Nr/2-1,4)
  REAL(q) :: phik(0:Nr,Nalp,Nalp), d2phidk2(0:Nr,Nalp,Nalp)
  REAL(q) , allocatable :: q0(:),dq0(:),pa(:),f_g(:,:)
  REAL(q) :: qs_vdw(Nalp)
  REAL(q) :: k_k(GRIDC%MPLWV)
  REAL(q),allocatable :: dj_k(:)
  REAL(q) :: energ,soft_pot,soft_der
  REAL(q) :: a_k, b_k, a2_k, b2_k, a3_k, b3_k, pref
  REAL(q) :: ro,gr,dq0dn, dq0dgr, dpadq0, dtheadn, dtheadgr
  !spin vdw mod
  REAL(q) :: ros,grs,dq0dns, dq0dgrs, dtheadns, dtheadgrs
  REAL(q) :: dqup_dq, dqdn_dq
  REAL(q) :: dhdx, pa_tm, vld1, vld2, vld3, vld4, phi, dphidk
  REAL(q) :: k_x(3), k_y, k_z, k_a, cell(3,3)
  REAL(q) :: max_n, max_gr, a_iofq, b_iofq
  INTEGER,allocatable :: j_k(:),alp(:)
  INTEGER :: n1,n2,n3,nc
  INTEGER ::  i_max_n, i_max_gr
  !previous mod, optimization?
  REAL(q), allocatable :: k_x_pre(:,:)
  REAL(q) :: str11,str12,str13,str21,str22,str23,str31,str32,str33
  !spin vdw mod (jk: not sure, needs to be checked)
  REAL(q) :: q0_temp
  ! fixing bad convergence for one spin density too low
  ! while the other OK
  REAL(q), parameter :: rhomin=2.e-10_q
  REAL(q), parameter :: spin_cut=4._q
  LOGICAL :: LDO_UP, LDO_DN
  REAL(q) :: qx_up, qx_dn, k_F
! jik not helped with "support procedure error"
  REAL(q) :: DQZXDN, DQZXDDN

GPACC !$ACC ROUTINE(CORPBE_LDA) SEQ

  NP=GRIDC%RL%NP
  !jk spin: none of these has a spin version
  ALLOCATE(q0(NP),dq0(NP),alp(NP),pa(NP),f_g(NP,Nalp))
  allocate(THETA(GRIDC%MPLWV,Nalp),f_k(GRIDC%MPLWV))
  allocate(dj_k(GRIDC%MPLWV),j_k(GRIDC%MPLWV))
  ALLOCATE(k_x_pre(4,GRIDC%RC%NP))
GPACC !$ACC ENTER DATA CREATE(q0,dq0,alp,pa,f_g,theta,f_k,dj_k,j_k) __IF_ASYNC__
GPACC !$ACC ENTER DATA CREATE(k_x_pre) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:Q0,DQ0,ALP,PA,F_G,THETA,F_K,DJ_K,J_K) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:K_X_PRE) IF(OFFLOAD_ON)

  cell=latt_cur%B*AUTOA
GPACC !$ACC ENTER DATA COPYIN(qs_vdw,spl,phik,d2phidk2,cell) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:QS_VDW,SPL,PHIK,D2PHIDK2,CELL) IF(OFFLOAD_ON)
  str11=0._q
  str21=0._q
  str31=0._q
  str12=0._q
  str22=0._q
  str32=0._q
  str13=0._q
  str23=0._q
  str33=0._q

  energ=0.0_q

#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
      IDUMP=0
#ifdef debugNODE_ME
      IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
      NODE_ME=0
      IONODE =0
      IDUMP=0
#ifdef debug
      IDUMP=1
#endif
#endif
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(DCHARG,DCHARGS,DWORK4,DWORK4S,q0) &
GPACC !$ACC& PRIVATE(ro,ros,gr,grs,s,ss,nrs,rnrs,zeta,ex_LDA,ec_LDA,vc_LDA,LDO_UP,LDO_DN,qx_up,qx_dn,KBEC,KBVU,KBVD,k_F) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:DCHARG,DCHARGS,DWORK4,DWORK4S,Q0) &
GPOMP !$OMP PRIVATE(RO,ROS,GR,GRS,S,SS,NRS,RNRS,ZETA,EX_LDA,EC_LDA,VC_LDA,LDO_UP,LDO_DN,QX_UP,QX_DN,KBEC,KBVU,KBVD,K_F) IF(OFFLOAD_ON)

!!!
  DO I=1,NP
    !spin vdw mod, need to take care of noisy densities in one or both spins
    q0(i)=q0cut
    ro=DCHARG(i)*AUTOA3  !convert to a.u.
    gr=DWORK4(i)*AUTOA4
    ros=DCHARGS(i)*AUTOA3  !convert to a.u.
    grs=DWORK4S(i)*AUTOA4

    LDO_UP=.TRUE.
    LDO_DN=.TRUE.

    ! skip this grid point if the density is too low
    ! default the q0 to q0cut
    !-! new way
    !actually crashes if commented out... for problematic cases
    ! as the total density can be negative and things like that
    IF ( (DCHARG(i)+DCHARGS(i) ) .LT. rhomin) CYCLE

    ! check each spin individually for low density
    ! if the density is too low, the q0x contribution
    ! will not be calculated for that spin

    IF (DCHARG(i) .LT. rhomin*0.5_q) LDO_UP=.FALSE.
    IF (DCHARGS(i) .LT. rhomin*0.5_q) LDO_DN=.FALSE.

    !calculate the qx contributions for both spins    
    qx_up=0._q
    IF (LDO_UP) THEN
       !difference to Q
       !put twice density and grad due to use of scaling relation
       !jik could not get this to work with acc
       !qx_up=QZX(2._q*ro,2._q*gr)
       !so inline these routines multiplying the density and grad twice  
       k_F=(3._q*pi*pi*2*ro)**(1./3.)
       s=abs(2*gr)/(2._q*k_F*2*ro)
       qx_up=(1.0_q-Zab_vdW/9._q*s*s)*k_F
       
       nrs=0._q
! jik? is the line below appropriate? assumed from similar directive below
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
       DO iu=1,12
         nrs=nrs+(qx_up/(spin_cut*q0cut))**iu/real(iu,q)
       ENDDO
       qx_up=spin_cut*q0cut*(1._q-exp(-nrs))
    ENDIF

    qx_dn=0._q
    IF (LDO_DN) THEN
       !put twice density and grad due to use of scaling relation
       !qx_dn=QZX(2._q*ros,2._q*grs)
       !line above works without acc
       k_F=(3._q*pi*pi*2*ros)**(1./3.)
       s=abs(2*grs)/(2._q*k_F*2*ros)
       qx_dn=(1.0_q-Zab_vdW/9._q*s*s)*k_F

       nrs=0._q
! jik? is the line below appropriate?
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
       DO iu=1,12
         nrs=nrs+(qx_dn/(spin_cut*q0cut))**iu/real(iu,q)
       ENDDO
       qx_dn=spin_cut*q0cut*(1._q-exp(-nrs))
    ENDIF

    ! this should make sense unless problematic rho is too negative
    ! then the ro+ros would be negative, in which case we cycle
    nrs=(3._q/(4._q*pi*(ro+ros)))**(1._q/3._q)
    rnrs=sqrt(nrs)

    ! LDA correlation is calculated for every point 
    ! with density larger than minimal
    zeta=(ro-ros)/(ro+ros)
    ! avoid going unphysical
    IF (ABS(zeta)>1.0_q) zeta = SIGN(1.0_q, zeta)

    CALL CORPBE_LDA(nrs, zeta, KBEC, KBVU, KBVD)
    ! KBEC, KBVU, KBVD returned in Rydbergs,  convert to Ha
    KBEC=KBEC*0.5_q; KBVU=KBVU*0.5_q; KBVD=KBVD*0.5_q

    ! eqs 5 to 8 in Thonhauser et al. PRL 115, 136402 (2015)
    q0(i)=(-4._q*pi/3._q)*KBEC + (qx_up*ro + qx_dn *ros)/(ro+ros)

    !now saturated
    nrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO iu=1,12
      nrs=nrs+(q0(i)/q0cut)**iu/real(iu,q)
    ENDDO
      q0(i)=q0cut*(1._q-exp(-nrs))
  ENDDO

  !calculate alp and dq0 in each point
  a_iofq=log((qs_vdw(Nalp)-qs_vdw(Nalp-1))/(qs_vdw(2)-qs_vdw(1)))/(Nalp-2.0_q)
  a_iofq=max(a_iofq,1e-12_q)
  b_iofq=(qs_vdw(2)-qs_vdw(1))/(exp(a_iofq)-1.0_q)
GPACC !$ACC PARALLEL PRESENT(q0,alp,dq0,qs_vdw) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(iofq)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:Q0,ALP,QS_VDW) PRIVATE(IOFQ) IF(OFFLOAD_ON)
  DO I=1,NP
    iofq=1+log(1+(q0(i)-qs_vdw(1))/b_iofq)/a_iofq
    iofq=max(1,iofq)
    iofq=min(Nalp,iofq)
    !!!
    IF(q0(i).GT.9.99985_q) iofq=Nalp
    alp(i)=iofq
  ENDDO
GPACC !$ACC LOOP GANG VECTOR
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DQ0,Q0,ALP,QS_VDW) IF(OFFLOAD_ON)
  DO I=1,NP
    dq0(i)=q0(i)-qs_vdw(alp(i))
  ENDDO
GPACC !$ACC END PARALLEL


  DO alf=1,Nalp
GPACC !$ACC PARALLEL PRESENT(spl,alp,dq0,pa,DCHARG,DCHARGS,dinp) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(coef,diff)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:SPL,ALP,DQ0,PA) &
GPOMP !$OMP PRIVATE(COEF,DIFF) IF(OFFLOAD_ON)
    DO I=1,NP
      coef(:)=spl(alp(i),alf,:)
      diff=dq0(i)
      pa(i)=coef(1)+diff*(coef(2)+diff*(coef(3)+diff*coef(4)))
    ENDDO
!spin vdw mod
!jik: not sure if the VECTOR will be still doing something here efficiently
GPACC !$ACC LOOP GANG VECTOR PRIVATE(nrs)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DCHARG,DCHARGS,DINP) REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO I=1,NP
      nrs=0._q
      IF (DCHARG(i).GE.rhomin/2._q) nrs=DCHARG(i)
      IF (DCHARGS(i).GE.rhomin/2._q) nrs=nrs+DCHARGS(i)

      dinp(i)=pa(i)*nrs*AUTOA3
    ENDDO
GPACC !$ACC END PARALLEL

    CALL OPSYNC(DINP,CINP,GRIDC%NPLWV)
    CALL FFT3D(CINP,GRIDC,-1)

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(THETA,CINP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:THETA,CINP) IF(OFFLOAD_ON)
    DO I=1,GRIDC%RC%NP
      THETA(i,alf)=CINP(i)
    ENDDO
  ENDDO

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(dj_k,k_k,j_k) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DJ_K,K_K,J_K) IF(OFFLOAD_ON)
  DO i=1,GRIDC%RC%NP
    dj_k(i)=0.5_q*k_k(i)/pi*rcut
    j_k(i)=floor(dj_k(i))
    dj_k(i)=dj_k(i)-j_k(i)
    dj_k(i)=dj_k(i)*2._q*pi/rcut
  ENDDO

!prev mod: prepare k_x array if needed, no spin obviously involved
  IF (ISIF .ge. 1) THEN
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(k_x_pre,gridc,cell) PRIVATE(n1,nc,n2,n3) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:K_X_PRE,GRIDC,CELL) &
GPOMP !$OMP PRIVATE(N1,NC,N2,N3) IF(OFFLOAD_ON)
     DO i=1,GRIDC%RC%NP
        n1=Mod((i-1),gridc%rc%nrow)+1
        nc=(i-1)/gridc%rc%nrow+1
        n2=gridc%rc%i2(nc)
        n3=gridc%rc%i3(nc)
        k_x_pre(2,i)=tpi*(gridc%lpctx(n1)*cell(1,1)+gridc%lpcty(n2)*cell(1,2)+gridc%lpctz(n3)*cell(1,3))
        k_x_pre(3,i)=tpi*(gridc%lpctx(n1)*cell(2,1)+gridc%lpcty(n2)*cell(2,2)+gridc%lpctz(n3)*cell(2,3))
        k_x_pre(4,i)=tpi*(gridc%lpctx(n1)*cell(3,1)+gridc%lpcty(n2)*cell(3,2)+gridc%lpctz(n3)*cell(3,3))
        k_x_pre(1,i)=sqrt(k_x_pre(2,i)*k_x_pre(2,i)+k_x_pre(3,i)*k_x_pre(3,i)+k_x_pre(4,i)*k_x_pre(4,i))  !*2._q*pi factor is here
     ENDDO
  ENDIF ! (ISIF .ge. 1)

  DO alf=1,Nalp
GPACC !$ACC KERNELS PRESENT(f_K) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:F_K) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    f_k=0._q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS


! jk: do the 2D thing with theta, no spin probably involved
GPACC !$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRESENT(k_k,phik,d2phidk2,f_k,theta,gridc,k_x_pre) &
GPACC !$ACC& PRIVATE(ik_k,phi,dphidk,a_k,b_k,a2_k,b2_k,a3_k,b3_k,n1,nc,n3,pref) &
GPACC !$ACC& REDUCTION(+:str11,str12,str13,str21,str22,str23,str31,str32,str33) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:K_K,PHIK,D2PHIDK2,F_K,THETA,GRIDC,K_X_PRE) &
GPOMP !$OMP PRIVATE(IK_K,PHI,DPHIDK,A_K,B_K,A2_K,B2_K,A3_K,B3_K,N1,NC,N3,PREF) &
GPOMP !$OMP REDUCTION(+:STR11,STR12,STR13,STR21,STR22,STR23,STR31,STR32,STR33) IF(OFFLOAD_ON)

    DO bet=1,Nalp
      DO i=1,GRIDC%RC%NP
        ik_k=k_k(i)/dk_k
        IF (ik_k.GE.Nr) THEN
          !I really don't like the if here
          phi=0._q
          dphidk=0._q
        ELSE
          a_k=((ik_k+1._q)*dk_k-k_k(i))/dk_k
          b_k=1._q-a_k
          a2_k=(3._q*a_k**2.-1._q)*dk_k/6._q
          b2_k=(3._q*b_k**2.-1._q)*dk_k/6._q
          a3_k=(a_k**3.-a_k)*dk_k**2._q/6._q
          b3_k=(b_k**3.-b_k)*dk_k**2._q/6._q
          phi=a_k*phik(ik_k, bet, alf)+b_k*phik(ik_k+1,bet,alf)&
              +a3_k*d2phidk2(ik_k,bet,alf)+b3_k*d2phidk2(ik_k+1,bet,alf)
          dphidk=(-phik(ik_k, bet, alf)+phik(ik_k+1,bet,alf))/dk_k &
                 -a2_k*d2phidk2(ik_k,bet,alf)+b2_k*d2phidk2(ik_k+1,bet,alf)

        ENDIF
#ifdef OFFLOADING
        CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(f_k(i),theta(i,bet)*phi)
#else // OFFLOADING
        f_k(i)=f_k(i)+theta(i,bet)*phi
#endif // OFFLOADING
        IF (ISIF .ge. 1) THEN
        if (k_x_pre(1,i) .gt. 1.0e-13_q) then
          pref=dphidk/k_x_pre(1,i)*theta(i,bet)*conjg(theta(i,alf))
#ifdef NGZhalf
          nc=(i-1)/gridc%rc%nrow+1
          n3=gridc%rc%i3(nc)
          if ((n3.ne.1).and.(n3.ne.(GRIDC%NGZ/2+1))) pref=2._q*pref
#elif defined(NGXhalf)
          n1=Mod((i-1),gridc%rc%nrow)+1
          if ((n1.ne.1).and.(n1.ne.(GRIDC%NGX/2+1))) pref=2._q*pref
#endif
          str11=str11- pref*k_x_pre( 1+1,i) * k_x_pre( 1+1,i)
          str12=str12- pref*k_x_pre( 1+1,i) * k_x_pre( 2+1,i)
          str13=str13- pref*k_x_pre( 1+1,i) * k_x_pre( 3+1,i)
          str21=str21- pref*k_x_pre( 2+1,i) * k_x_pre( 1+1,i)
          str22=str22- pref*k_x_pre( 2+1,i) * k_x_pre( 2+1,i)
          str23=str23- pref*k_x_pre( 2+1,i) * k_x_pre( 3+1,i)
          str31=str31- pref*k_x_pre( 3+1,i) * k_x_pre( 1+1,i)
          str32=str32- pref*k_x_pre( 3+1,i) * k_x_pre( 2+1,i)
          str33=str33- pref*k_x_pre( 3+1,i) * k_x_pre( 3+1,i)
        endif
        ENDIF !only calculate stress when required

      ENDDO
    ENDDO !beta

GPACC !$ACC PARALLEL PRESENT(gridc,f_k,theta,CINP) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(nc,n1,n3,pref) REDUCTION(+:energ)

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,F_K,THETA) &
GPOMP !$OMP PRIVATE(NC,N1,N3,PREF) REDUCTION(+:ENERG) IF(OFFLOAD_ON)

    DO i=1,GRIDC%RC%NP
      pref=real(f_k(i)*conjg(theta(i,alf)),kind=q)
#ifdef NGZhalf
!for MPI version reduced in Z direction
!and then we sum it with \Theta_\alpha
      nc=(i-1)/gridc%rc%nrow+1
      n3=gridc%rc%i3(nc)
      if ((n3.ne.1).and.(n3.ne.(GRIDC%NGZ/2+1))) pref=2._q*pref
#else
#ifdef NGXhalf
!for serial version reduced in X direction
      n1=Mod((i-1),gridc%rc%nrow)+1
      if ((n1.ne.1).and.(n1.ne.(GRIDC%NGX/2+1))) pref=2._q*pref
#endif
#endif
      energ=energ+pref
    ENDDO

GPACC !$ACC LOOP GANG VECTOR
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CINP,F_K) IF(OFFLOAD_ON)
    DO i=1,GRIDC%RC%NP
      CINP(I)=f_k(i)
    ENDDO
GPACC !$ACC END PARALLEL

    CALL OPSYNC(DINP,CINP,GRIDC%NPLWV)
    CALL FFT3D(CINP,GRIDC,1)

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(dinp,f_g,gridc) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DINP,F_G,GRIDC) IF(OFFLOAD_ON)
    DO i=1,GRIDC%RL%NP
      f_g(i,alf)=dinp(i)/(real(GRIDC%NGX,q)*GRIDC%NGY*GRIDC%NGZ)
    ENDDO
  ENDDO !alf

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
  EXC=EXC + 0.5_q*energ*LATT_CUR%OMEGA/(real(GRIDC%NGX,q)*GRIDC%NGY*GRIDC%NGZ)/AUTOA3*2._q*RYTOEV
  !compared to vander, we divide by (M)PLWV only once since it's done in xcgrad for the second time
  !so we can say that we get rid of the FFT 1/N factor and keep the \Omega/N for xcgrad
  E_corr=0.5_q*energ*LATT_CUR%OMEGA/(real(GRIDC%NGX)*GRIDC%NGY*GRIDC%NGZ)**2./AUTOA3*2._q*RYTOEV

  max_n=0.d0
  max_gr=0.d0

!vdw spin mod whole loop
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(DCHARG,DCHARGS,DWORK4,DWORK4S,q0,spl,alp,dq0,f_g,DWORKG,DWORKGS,DWORK,DWORKS) &
GPACC !$ACC& PRIVATE(ro,ros,gr,grs,nrs,rnrs,vld1,vld2,vld3,vld4,ec_LDA,vc_LDA,dq0dn,dq0dns,dq0dgr,dhdx, &
GPACC !$ACC& dtheadn,dtheadns,dtheadgr,dtheadgrs,coef,diff,pa_tm,dpadq0,LDO_UP,LDO_DN,dqup_dq,dqdn_dq, &
GPACC !$ACC& qx_up,qx_dn,dq0dgrs,KBEC,KBVU,KBVD,k_F,s,DQZXDN,zeta) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:DCHARG,DCHARGS,DWORK4,DWORK4S,Q0,SPL,ALP,DQ0, &
GPOMP !$OMP F_G,DWORKG,DWORKGS,DWORK,DWORKS) PRIVATE(RO,ROS,GR,GRS,NRS,RNRS,VLD1,VLD2,VLD3,VLD4,EC_LDA, &
GPOMP !$OMP VC_LDA,DQ0DN,DQ0DNS,DQ0DGR,DHDX,DTHEADN,DTHEADNS,DTHEADGR,DTHEADGRS,COEF,DIFF,PA_TM,DPADQ0, &
GPOMP !$OMP LDO_UP,LDO_DN,DQUP_DQ,DQDN_DQ,QX_UP,QX_DN,DQ0DGRS,KBEC,KBVU,KBVD,K_F,S,DQZXDN,ZETA) IF(OFFLOAD_ON)

  DO i=1,gridc%RL%NP
    !calc the potential
    LDO_UP=.TRUE.
    LDO_DN=.TRUE.
    ro=DCHARG(i)*AUTOA3  !convert to a.u.
    gr=DWORK4(i)*AUTOA4
    ros=DCHARGS(i)*AUTOA3
    grs=DWORK4S(i)*AUTOA4

    !crashes for potential if commented out
    !as the total density can be negative and things like that
    IF ( (DCHARG(i)+DCHARGS(i) ) .LT. rhomin) CYCLE
    !don't do the exchange contributions when density is too small
    !for either of the spins
    IF (DCHARG(i) .LT. rhomin*0.5_q) LDO_UP=.FALSE.
    IF (DCHARGS(i) .LT. rhomin*0.5_q) LDO_DN=.FALSE.
 
    dqup_dq=0._q
    dqdn_dq=0._q

    IF (LDO_UP .eqv. .TRUE.) THEN
       !qx_up=QZX(2._q*ro,2._q*gr)
       k_F=(3._q*pi*pi*2*ro)**(1./3.)
       s=abs(2*gr)/(2._q*k_F*2*ro)
       qx_up=(1.0_q-Zab_vdW/9._q*s*s)*k_F

       nrs=0._q
       rnrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
       DO iu=1,12
         nrs=nrs+(qx_up/(spin_cut*q0cut))**iu/real(iu,q)
       ENDDO
       qx_up=spin_cut*q0cut*(1._q-exp(-nrs))
       
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:RNRS) IF(OFFLOAD_ON)
       DO iu=0,11
         rnrs=rnrs+(qx_up/(spin_cut*q0cut))**iu
       ENDDO
       dqup_dq=rnrs*exp(-nrs)
    ENDIF

    IF (LDO_DN .eqv. .TRUE.) THEN
       !qx_dn=QZX(2._q*ros,2._q*grs)
       k_F=(3._q*pi*pi*2*ros)**(1./3.)
       s=abs(2*grs)/(2._q*k_F*2*ros)
       qx_dn=(1.0_q-Zab_vdW/9._q*s*s)*k_F

       nrs=0._q
       rnrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
       DO iu=1,12
         nrs=nrs+(qx_dn/(spin_cut*q0cut))**iu/real(iu,q)
       ENDDO
       qx_dn=spin_cut*q0cut*(1._q-exp(-nrs))
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:RNRS) IF(OFFLOAD_ON)
       DO iu=0,11
         rnrs=rnrs+(qx_dn/(spin_cut*q0cut))**iu
       ENDDO
       dqdn_dq=rnrs*exp(-nrs)
    ENDIF

    dq0dn=0._q
    dq0dns=0._q

    IF (LDO_UP) THEN
       ! add qx_up terms to the up derivative
       DQZXDN=(pi/3._q/(2.*ro))**(2./3.)+(7._q*Zab_vdW/108._q)*(2.*gr)*(2.*gr)*(2.*ro)**(-10./3._q)/(3.*pi*pi)**(1./3.)
       dq0dn =  ros/(ro+ros)**2.*qx_up + ro/(ro+ros)*2._q*DQZXDN*dqup_dq
       ! add the qx_up term to the down derivative
       dq0dns= -ro/(ro+ros)**2.*qx_up
    ENDIF
    IF (LDO_DN) THEN
       DQZXDN=(pi/3._q/(2.*ros))**(2./3.)+(7._q*Zab_vdW/108._q)*(2.*grs)*(2.*grs)*(2.*ros)**(-10./3._q)/(3.*pi*pi)**(1./3.)
       dq0dns=dq0dns+ro/(ro+ros)**2.*qx_dn + ros/(ro+ros)*2._q*DQZXDN*dqdn_dq
       ! add the qx_down term to the up derivative
       dq0dn=dq0dn-qx_dn*ros/(ro+ros)**2.
    ENDIF
    
    ! correlation is always calculated
    nrs=(3._q/(4._q*pi*(ro+ros)))**(1._q/3._q)
    rnrs=sqrt(nrs)
    zeta=(ro-ros)/(ro+ros)
    IF (ABS(zeta)>1.0_q) zeta = SIGN(1.0_q, zeta)
    !use LDA routines
    CALL CORPBE_LDA(nrs,zeta, KBEC, KBVU, KBVD)   ! KBEC is LDAc energy density
    KBEC=KBEC*0.5_q; KBVU=KBVU*0.5_q; KBVD=KBVD*0.5_q
 
    !add the correlation contributions to the potentials
    dq0dn=dq0dn+(4.*pi/3.)*(KBEC-KBVU)/(ro+ros)
    dq0dns=dq0dns+(4.*pi/3.)*(KBEC-KBVD)/(ro+ros)

    ! calculate the derivatives only when the density is above threshold
    dq0dgr=0._q
    dq0dgrs=0._q
    IF (LDO_UP) THEN
       DQZXDDN=-Zab_vdW/18._q*(2*gr)/(3._q*pi*pi*(2*ro))**(1./3.)/(4*ro*ro)
       dq0dgr=ro/(ro+ros)*DQZXDDN*2._q*dqup_dq
    ENDIF
    IF (LDO_DN) THEN
       DQZXDDN=-Zab_vdW/18._q*(2*grs)/(3._q*pi*pi*(2*ros))**(1./3.)/(4*ros*ros)
       !dq0dgrs=ros/(ro+ros)*DQZXDDN(2*ros,2*grs)*2._q*dqdn_dq
       dq0dgrs=ros/(ro+ros)*DQZXDDN*2._q*dqdn_dq
    ENDIF

    nrs=0._q
    rnrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:RNRS) IF(OFFLOAD_ON)
    DO iu=1,12
      rnrs=rnrs+(q0(i)/q0cut)**iu/real(iu,q)
    ENDDO
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO iu=0,11
      nrs=nrs+(q0(i)/q0cut)**iu
    ENDDO
    dhdx=nrs*exp(-rnrs)

    !vdw spin mod
    dtheadn=0._q
    dtheadgr=0._q
    dtheadns=0._q
    dtheadgrs=0._q

GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:DTHEADN,DTHEADGR,DTHEADNS,DTHEADGRS) IF(OFFLOAD_ON)
    DO alf=1,Nalp
      coef=spl(alp(i),alf,:)
      diff=dq0(i)
      pa_tm=coef(1)+diff*(coef(2)+diff*(coef(3)+diff*coef(4)))
      dpadq0=coef(2)+diff*(2._q*coef(3)+3._q*diff*coef(4))
      dtheadn=dtheadn+(pa_tm+(ro+ros)*dpadq0*dq0dn*dhdx)*f_g(i,alf)
      dtheadgr=dtheadgr+(ro+ros)*dpadq0*dhdx*dq0dgr*f_g(i,alf)
      dtheadns=dtheadns+(pa_tm+(ro+ros)*dpadq0*dq0dns*dhdx)*f_g(i,alf)
      dtheadgrs=dtheadgrs+(ro+ros)*dpadq0*dhdx*dq0dgrs*f_g(i,alf)
    ENDDO

    !this is d/dn
    DWORKG(i)=DWORKG(i)+dtheadn*2._q*RYTOEV
    DWORKGS(i)=DWORKGS(i)+dtheadns*2._q*RYTOEV
    !this is d/ddn
    DWORK(i)=DWORK(i)+dtheadgr*2._q*RYTOEV* AUTOA
    DWORKS(i)=DWORKS(i)+dtheadgrs*2._q*RYTOEV* AUTOA
    ! in the spin polarized case the gradient is used later MAX(DWORK4(I),1.E-10_q)

  ENDDO

GPACC !$ACC EXIT DATA DELETE(qs_vdw,spl,phik,d2phidk2,cell) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(q0,dq0,alp,pa,f_g,theta,f_k,dj_k,j_k) __IF_ASYNC__
GPACC !$ACC EXIT DATA DELETE(k_x_pre) __IF_ASYNC__

GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:QS_VDW,SPL,PHIK,D2PHIDK2,CELL) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:Q0,DQ0,ALP,PA,F_G,THETA,F_K,DJ_K,J_K) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:K_X_PRE) IF(OFFLOAD_ON)

  DEALLOCATE(q0,dq0,alp,pa)
  DEALLOCATE(dj_k,theta,f_k,f_g,j_k)
  DEALLOcATE(k_x_pre)
GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
  stress(1,1)=str11
  stress(2,1)=str21
  stress(3,1)=str31
  stress(1,2)=str12
  stress(2,2)=str22
  stress(3,2)=str32
  stress(1,3)=str13
  stress(2,3)=str23
  stress(3,3)=str33
  stress=0.5_q*stress*LATT_CUR%OMEGA/(real(GRIDC%NGX,kind=q)*GRIDC%NGY*GRIDC%NGZ)**2.*2._q*RYTOEV/AUTOA3

#ifdef MPI
  CALLMPI( M_sum_d(GRIDC%COMM,e_corr,1))
#endif

! IF (NODE_ME == IONODE) then
!   WRITE(*,*) 'Total vdW correction in eV: ', e_corr
! ENDIF

! CONTAINS

! jik not sure how to deal with this, I keep getting 
! NVFORTRAN-S-1062-Support procedure called within a compute region -
! __kmpc_threadprivate_cached
! !   FUNCTION QZX(rho, grad)
! !!$ACC ROUTINE SEQ
! !    REAL(q) :: rho, grad
! !    REAL(q) :: QZX
! !   !internal 
! !    REAL(q) :: k_F, s
! !
! !    k_F=(3._q*pi*pi*rho)**(1./3.)
! !    s=abs(grad)/(2._q*k_F*rho)   ! rho is set to +delta when lower than that
! !    QZX=(1.0_q-Zab_vdW/9._q*s*s)*k_F
! !   END FUNCTION QZX
! !
! !   FUNCTION DQZXDN(rho, grad)
! !!$ACC ROUTINE SEQ
! !    REAL(q) :: rho, grad
! !    REAL(q) :: DQZXDN
! !   !internal 
! !
! !    DQZXDN=(pi/3._q/rho)**(2./3.)+(7._q*Zab_vdW/108._q)*grad*grad*rho**(-10./3._q)/(3.*pi*pi)**(1./3.)
! !   END FUNCTION DQZXDN
! !
! !   FUNCTION DQZXDDN(rho, grad)
! !!$ACC ROUTINE SEQ
! !    REAL(q) :: rho, grad
! !    REAL(q) :: DQZXDDN
! !
! !    DQZXDDN=-Zab_vdW/18._q*grad/(3._q*pi*pi*rho)**(1./3.)/(rho*rho)
! !   END FUNCTION DQZXDDN

END SUBROUTINE VDW_CALC_spin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> rVV10
!>
!> 2015/05/14 \n
!> Modified by Zenghui Yang, and Haowei Peng \n
!> to evaluate the rVV10 non-local correlation functional of Phys. Rev. B 87, 041108(R) (2013),
!> which is a simplified/revised version of the VV10 functional of J. Chem. Phys. 133, 244103 (2010).
!>
!> Please cite \n
!>     Haowei Peng*, Zeng-Hui Yang, Jianwei Sun, and John P. Perdew,
!>     "Versatile van der Waals Density Functional Based on a Meta-Generalized
!>      Gradient Approximation", PRX 6, 041005 (2016).
!>
!>     Haowei Peng* and John P. Perdew,
!>     "Rehabilitation of the Perdew-Burke-Ernzerhof generalized gradient approximation
!>      for layered materials", PRB 95, 081105(R) (2017).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE VDW_rVV10
  USE prec
  USE lattice
  USE mpimy
  USE mgrid
  USE constant
  USE vdw_df_helper, ONLY : IU0,IU6,ISIF

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: VDW_NONLOC_RVV10,Nalp,Nr,q0cut,rcut

! REAL(q),parameter :: bparam = 6.3_q !rVV10, not the 5.9 of VV10
! REAL(q),parameter :: betaparam = 0.03125_q * (3.0_q / (bparam**2.0_q) )**(0.75_q) * 2._q*RYTOEV  ! Already in eV
! REAL(q),parameter :: Cparam = 0.0093_q

  INTEGER, parameter :: Nalp=24, Nr=2048
  REAL(q), parameter :: qAlp_ratio=2048._q
  REAL(q), parameter :: q0cut=0.5_q
  REAL(q), parameter :: rcut = 100._q
  REAL(q) :: spl(Nalp,Nalp,4)
  REAL(q) :: phik(0:Nr,Nalp,Nalp), d2phidk2(0:Nr,Nalp,Nalp)
  REAL(q) :: phir(0:Nr,Nalp,Nalp)
  REAL(q) :: qs_vdw(Nalp)
  REAL(q) :: dr

  LOGICAL :: INIT_VDW_LL=.TRUE. !do we need to initialise the stuff?
  LOGICAL :: LREAD_KERN = .FALSE.
  logical :: lhave_alpha_table = .FALSE.

  REAL(q) :: kcut

CONTAINS

  SUBROUTINE VDW_NONLOC_RVV10(DCHARG,DWORK,DWORKG,DWORK4,EXC,GRIDC,LATT_CUR,STRESS,ISPIN,bparam,cparam)

__GPU USE moffload_struct_def

    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    REAL(q) bparam,cparam
    REAL(q) DCHARG(GRIDC%RL%NP) !charge density on input and output
    REAL(q) DWORKG(GRIDC%RL%NP) !df/drho input, add the same non-local here
    REAL(q) DWORK4(GRIDC%RL%NP) !gradient on input !df/d(drho)/(drho)potential on output
    REAL(q) EXC
    REAL(q) STRESS(3,3)
    COMPLEX(q), ALLOCATABLE :: CINP(:) !, THETA(:,:)
    REAL(q), ALLOCATABLE :: K_K(:)
    RGRID   DWORK(GRIDC%RL%NP)
    INTEGER ISPIN

    IF (INIT_VDW_LL) THEN
      CALL MIN_K(GRIDC,LATT_CUR)
      CALL GET_SPLINES
      CALL GET_KERNELS
      INIT_VDW_LL=.FALSE.
    ENDIF

    ALLOCATE(K_K(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(K_K) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:K_K) IF(OFFLOAD_ON)

    CALL GET_K_IND(K_K,GRIDC,LATT_CUR)

    ALLOCATE(CINP(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(CINP) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:CINP) IF(OFFLOAD_ON)

    CALL VDW_CALC_RVV10(DCHARG,DWORK,DWORKG,DWORK4,EXC,GRIDC,LATT_CUR,CINP,CINP,SPL,PHIK,D2PHIDK2,QS_VDW,K_K,STRESS,ISIF,ISPIN,bparam,cparam)

GPACC !$ACC EXIT DATA DELETE(CINP,K_K) __IF_ASYNC__
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:CINP,K_K) IF(OFFLOAD_ON)
    DEALLOCATE(CINP)
    DEALLOCATE(k_k)
  END SUBROUTINE VDW_NONLOC_RVV10

    !> Construct the spline coefficients for interpolation.
  SUBROUTINE GET_SPLINES
    REAL(q) :: a(Nalp,Nalp), c(Nalp,Nalp), l(Nalp,Nalp), mi(Nalp,Nalp), z(Nalp,Nalp)
    REAL(q) :: b(Nalp-1,Nalp), d(Nalp-1,Nalp), h(Nalp-1), al(Nalp-1,Nalp)
    INTEGER :: i,j

    spl=0._q
    a=0._q
    b=0._q
    d=0._q
    h=0._q
    al=0._q
    c=0._q
    l=0.0_q
    mi=0._q
    z=0._q

    call GEN_ALPHA_POINTS()

    DO i=1,Nalp
      a(i,i)=1._q
    ENDDO
    DO i=1,Nalp-1
      h(i)=qs_vdw(i+1)-qs_vdw(i)
    ENDDO
    DO i=2,Nalp-1
      DO j=1,Nalp
        al(i,j)=3._q/h(i)*(a(i+1,j)-a(i,j))-3._q/h(i-1)*(a(i,j)-a(i-1,j))
      ENDDO
    ENDDO
    l(1,:)=1._q
    DO i=2,Nalp-1
      l(i,:)=2._q*(qs_vdw(i+1)-qs_vdw(i-1))-h(i-1)*mi(i-1,:)
      mi(i,:)=h(i)/l(i,:)
      z(i,:)=(al(i,:)-h(i-1)*z(i-1,:))/l(i,:)
    ENDDO
    l(Nalp,:)=1._q
    z(Nalp,:)=0._q
    c(Nalp,:)=0._q
    DO i=Nalp-1,1,-1
      c(i,:)=z(i,:)-mi(i,:)*c(i+1,:)
      b(i,:)=(a(i+1,:)-a(i,:))/h(i)-h(i)/3._q*(c(i+1,:)+2._q*c(i,:))
      d(i,:)=(c(i+1,:)-c(i,:))/(3._q*h(i))
    ENDDO

    spl(:,:,1)=a(:,:)
    spl(:,:,3)=c(:,:)
    DO i=1,Nalp-1
      spl(i,:,2)=b(i,:)
      spl(i,:,4)=d(i,:)
    ENDDO

  END SUBROUTINE GET_SPLINES

    !> This procedure calculates the fourier transformed vdw kernels.
  SUBROUTINE GET_KERNELS
    IMPLICIT NONE
    INTEGER :: alf,bet,i,j,irad, nk, ik
    REAL (q) :: q1,q2, dk, krad
    REAL (q) :: r(0:nr)
    REAL (q) :: phir_work(0:nr), phik_work(0:nr)
    REAL (q) :: d2phik_work(0:nr)

    !setting up the mesh of interpolating points
    !call GEN_ALPHA_POINTS()
    dr= rcut/Nr
    dk= pi/rcut
    nk= int(kcut/dk)+1

    !find the interpolation coefficients for the alpha points
    DO alf = 1, Nalp
      DO bet = 1, alf
        phik_work=0._q
        phir_work=0._q
        !desaturate q interpolation points
        q1=DESATURATE(qs_vdw(alf))
        q2=DESATURATE(qs_vdw(bet))

        !get the kernel on a real grid by using the precalculated interpolation
        !data
        DO irad = 0, Nr
          r(irad) = irad * dr
          phir(irad,alf,bet) = PHI_VALUE(q1*(r(irad)**2), q2*(r(irad)**2))
        ENDDO

        !smoothen for large r
        DO irad =0, Nr
          phir(irad,alf,bet)=phir(irad,alf,bet)*S_CUT(r(irad)/rcut)
        ENDDO

        !fourier transform it and smoothen
        phir_work(:)=phir(:,alf,bet)
        call PHIK_FFT(phir_work(:), phik_work(:))
        phik(:,alf,bet)=phik_work*(2._q*pi)**1.5_q
        phik(nk:nr,alf,bet)=0._q
        DO ik =1,nk
          krad=ik*dk
          phik(ik,alf,bet)=phik(ik,alf,bet)*S_CUT(krad/kcut)
        ENDDO

        !get coefficients for cubic splines
        !phik_work to splines
        call SPLINE(dk,phik(:,alf,bet),Nr+1,0._q,0._q,d2phik_work)

        !symmetrise
        phik(:,bet,alf)=phik(:,alf,bet)
        d2phidk2(:,bet,alf)=d2phik_work
        d2phidk2(:,alf,bet)=d2phik_work
      ENDDO
    ENDDO

  END SUBROUTINE GET_KERNELS

    !> Find the interpolation points for the kernel interpolation q_0 grid
    !> i.e. fill the qs_vdw array.
  SUBROUTINE GEN_ALPHA_POINTS()
    REAL(q) :: a,b,e
    INTEGER :: i

    qs_vdw(:)=0._q
    a=log(qAlp_ratio)/(Nalp-1._q)
    b=q0cut/(qAlp_ratio-1._q)

    DO i=2,Nalp-1
      qs_vdw(i)=b*exp(a*(i-1._q))-b
    ENDDO
    qs_vdw(Nalp)=q0cut
    lhave_alpha_table = .true.

  END SUBROUTINE GEN_ALPHA_POINTS

    !> Find q0 corresponding to saturated q0 by interval halving.
  FUNCTION DESATURATE(x)
    REAL(q) :: desaturate, x
    REAL(q) :: xdo, xup, xme, xmesat
    REAL(q), parameter :: xtol=1.e-14_q

    !check first
    IF ((x<0._q).OR.(x>q0cut)) then
      CALL vtutor%error("DESATURATE: out of range " // str(x))
    endif
    !this is done for the original q0 values before saturation
    xdo =  0._q
    xup = q0cut*5._q
    DO
      xme = (xdo+xup)*0.5_q
      xmesat = SATURATE(xme)
      !if smaller
      IF (abs(xmesat-x).LE.xtol ) THEN
        desaturate=xup
        RETURN
      ELSE IF (x.LE.xmesat) THEN
        xup=xme
      ELSE
        xdo=xme
      ENDIF
    ENDDO

  END FUNCTION DESATURATE

    !> The same that we use down in the main routine
    !> R-P&S' equation No 4 or so.
  FUNCTION SATURATE(x)
    REAL(q):: saturate, x, coe
    INTEGER :: iu

    coe=0._q
    DO iu=1,12
      coe=coe+(x/q0cut)**iu/real(iu,kind=q)
    ENDDO
    saturate=q0cut*(1._q-exp(-coe))
  END FUNCTION SATURATE

    !> Finds the kernel using phi_generate.
    !> Returns the value of kernel for d1, d2.
  FUNCTION PHI_VALUE(d1, d2)
    REAL (q) :: phi_value
    REAL (q) :: d1, d2

    phi_value = -1.5_q/((d1+1._q)*(d2+1._q)*(d1+d2+2._q))

  END FUNCTION PHI_VALUE

    !> This performs the FFT of the kernel which is spherically symmetric
    !> and thus one just needs to perform sin(kr) transformation in principle
    !> times the $k^{-1} r$ factors.
    !> The sin transforms are used in setlocalpp.F so I guess I copy FOUR1 from there.
    !> The radial array is [2^N] anyway.
    !> We need to know Nr, rcut.
  SUBROUTINE PHIK_FFT(phir_work, phik_work)
    REAL (q) :: phir_work(0:nr), phik_work(0:nr)
    REAL (q) :: work(2,0:2*nr)
    REAL (q) :: dr, dk, fac, kval, rval
    INTEGER :: ik, ir, ii
    dr=rcut/Nr
    fac=dr/sqrt(2._q*pi)
    dk=pi/rcut
    work=0._q

    do ii=1, 2*Nr-1
      !this can be loop over half but I'm tired
      if (ii.lt.Nr) then
        ir = ii
        rval = dr*ir
        work(2,ii) = -fac*phir_work(ii)*rval
      else if (ii.eq.Nr) then
        work(2,ii) = 0._q
      else
        ir = 2*Nr-ii
        rval = -dr*ir
        work(2,ii) = -fac*phir_work(ir)*rval
      endif
    enddo

    call FOUR1(work, 2*Nr, 1)
    do ik=1, Nr
      kval=ik*dk
      phik_work(ik)=work(1,ik)/kval
    enddo

    !calculate phik_work(0)
    phik_work(0)=0._q
    do ii=1,Nr
      rval=dr*ii
      phik_work(0)=phik_work(0)+rval*rval*phir_work(ii)
    enddo
    phik_work(0)=phik_work(0)*2._q*fac

  END SUBROUTINE PHIK_FFT

  SUBROUTINE FOUR1(DATA,NN,ISIGN)
      USE prec
      IMPLICIT NONE
      REAL(q) WR,WI,WPR,WPI,WTEMP,THETA,TEMPI,TEMPR
      REAL(q) DATA(*)
      INTEGER M,ISTEP,MMAX,ISIGN,NN,N,I,J
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959_q/(ISIGN*MMAX)
        WPR=-2._q*DSIN(0.5_q*THETA)**2
        WPI=DSIN(THETA)
        WR=1._q
        WI=0._q
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
  END SUBROUTINE FOUR1

    !> Cut off function for the kernels,
    !> personally this looks quite harsh.
  FUNCTION S_CUT(val)
    REAL(q) :: val, s_cut
    if (val .le. 0._q) then
      s_cut=1._q
    else if (val .ge. 1) then
      s_cut=0._q
    else
      s_cut=(1._q-val**8)**4
    endif

  END FUNCTION S_CUT

      !> From Num recipes.
  SUBROUTINE SPLINE(dx,y,n,yp1,ypn,y2)
      INTEGER n
      REAL(q) :: yp1,ypn,x(n),y(n),y2(n),dx
      INTEGER i,k
      REAL(q) :: p,qm,sig,un,u(n)

      do i=1,n
       x(i)=dx*i
      enddo
      if (yp1.gt..99e30_q) then
        y2(1)=0._q
        u(1)=0._q
      else
        y2(1)=-0.5_q
        u(1)=(3._q/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do  i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2._q
        y2(i)=(sig-1._q)/p
        u(i)=(6._q*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30_q) then
        qm=0._q
        un=0._q
      else
        qm=0.5_q
        un=(3._q/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qm*u(n-1))/(qm*y2(n-1)+1._q)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
  END  SUBROUTINE SPLINE


  SUBROUTINE MIN_K(GRIDC,LATT_CUR)
    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    INTEGER :: I, N1, N2, N3, NC, ix, iy, iz
    REAL(q) :: cell(3,3)
    REAL(q) :: k_x,k_y,k_z,kte
    !the latt_cur%B has rec vector x stored as 11 21 31
    !and similarly for the rest
    !to get total coordinate x in k-space sum is done over
    !11 12 13
    N1=GRIDC%NGX
    N2=GRIDC%NGY
    N3=GRIDC%NGZ
    kcut = 10.e10_q
    cell=latt_cur%B*AUTOA
    !now we want to know the smallest boundary k-vector
    !IF (IU0>=0) WRITE(IU0,*) 'reciprocal cell'
    !IF (IU0>=0) WRITE(IU0,*) cell
    !IF (IU0>=0) WRITE(IU0,*) N1*cell(1,1), N2*cell(1,2), N3*cell(1,3)
    !IF (IU0>=0) WRITE(IU0,*) N1*cell(2,1),N2*cell(2,2), N3*cell(2,3)
    !IF (IU0>=0) WRITE(IU0,*) N1*cell(3,1), N2*cell(3,2), N3*cell(3,3)

    do ix = -1,1
      do iy = -1,1
        do iz = -1,1
          if ((abs(ix)+abs(iy)+abs(iz)).ne.0 ) then
            k_x=N1*cell(1,1)*ix+N2*cell(1,2)*iy+N3*cell(1,3)*iz
            k_y=N1*cell(2,1)*ix+N2*cell(2,2)*iy+N3*cell(2,3)*iz
            k_z=N1*cell(3,1)*ix+N2*cell(3,2)*iy+N3*cell(3,3)*iz
            kte=2._q*pi*sqrt(k_x*k_x+k_y*k_y+k_z*k_z)*0.5_q
            kcut=min(kcut, kte)
          endif
        enddo
      enddo
    enddo

  END SUBROUTINE MIN_K


  SUBROUTINE GET_K_IND(k_k,GRIDC,LATT_CUR)

__GPU USE moffload_struct_def

    !get the k-space interpolation indices
    !can be done once
    TYPE (latt) LATT_CUR
    TYPE (grid_3d) GRIDC
    REAL(q) k_k(GRIDC%MPLWV)
    INTEGER :: I, N1, N2, N3, NC
    REAL(q) :: cell(3,3)
    REAL(q) :: k_x, k_y, k_z

    cell=latt_cur%B*AUTOA

GPACC !$ACC PARALLEL LOOP PRESENT(GRIDC,K_K) COPYIN(CELL) PRIVATE(N1,N2,N3,NC,K_X,K_Y,K_Z) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,K_K) PRIVATE(N1,N2,N3,NC,K_X,K_Y,K_Z) MAP(TO:CELL) IF(OFFLOAD_ON)

    do i=1,gridc%rc%np
      n1=Mod((i-1),gridc%rc%nrow)+1
      nc=(i-1)/gridc%rc%nrow+1
      n2=gridc%rc%i2(nc)
      n3=gridc%rc%i3(nc)
      k_x=gridc%lpctx(n1)*cell(1,1)+gridc%lpcty(n2)*cell(1,2)+gridc%lpctz(n3)*cell(1,3)
      k_y=gridc%lpctx(n1)*cell(2,1)+gridc%lpcty(n2)*cell(2,2)+gridc%lpctz(n3)*cell(2,3)
      k_z=gridc%lpctx(n1)*cell(3,1)+gridc%lpcty(n2)*cell(3,2)+gridc%lpctz(n3)*cell(3,3)
      k_k(i)=2._q*pi*sqrt(k_x*k_x+k_y*k_y+k_z*k_z)
    enddo

  END SUBROUTINE GET_K_IND


END MODULE VDW_rVV10



SUBROUTINE VDW_CALC_RVV10(DCHARG,DWORK,DWORKG,DWORK4,EXC,GRIDC,LATT_CUR,CINP,DINP,spl,phik,d2phidk2,qs_vdw,k_k,stress,ISIF,ISPIN,bparam,cparam)

__GPU USE moffload_struct_def

  USE prec
  USE lattice
  USE mpimy
  USE mgrid
  USE constant
  USE VDW_rVV10, ONLY : Nalp, Nr, q0cut, rcut

GPACC !$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

  IMPLICIT NONE
  INTEGER NP
  REAL(q) RINPL, NELECT
!I guess LATT_CURR is the one that we need not OMEGA only
  TYPE (latt) LATT_CUR
  TYPE (grid_3d) GRIDC
  REAL(q) OMEGA
  REAL(q) bparam,cparam
  REAL(q) EXC
  REAL(q) DCHARG(GRIDC%RL%NP) !charge density on input and output
  REAL(q) DWORKG(GRIDC%RL%NP) !df/drho input, add the same non-local here
                    !don't touch this one yet
  REAL(q) DWORK4(GRIDC%RL%NP) !gradient on input
                    !df/d(drho)/(drho)potential on output
  REAL(q) stress(3,3)
  INTEGER :: I, IU, alf, bet, ik_k, iofq, ii, jj, ISIF, ISPIN
  INTEGER :: NODE_ME, IONODE, IDUMP
  COMPLEX(q) :: CINP(GRIDC%MPLWV)
  COMPLEX(q), allocatable::   THETA(:,:),f_k(:)

  RGRID :: DINP(DIMREAL(GRIDC%MPLWV))
  REAL(q),parameter :: dk_k=pi/rcut
  REAL(q) :: nrs, rnrs
  REAL(q) :: diff, coef(4), E_corr
  REAL(q) :: spl(Nalp,Nalp,4)
  REAL(q) :: phik(0:Nr,Nalp,Nalp), d2phidk2(0:Nr,Nalp,Nalp)
  REAL(q) , allocatable :: q0(:),dq0(:),pa(:),f_g(:,:)
  REAL(q) :: qs_vdw(Nalp)
  REAL(q) :: k_k(GRIDC%MPLWV)
  REAL(q) :: energ
  REAL(q) :: a_k, b_k, a2_k, b2_k, a3_k, b3_k, pref
  REAL(q) :: ro,gr,dq0dn, dq0dgr, dpadq0, dtheadn, dtheadgr
  REAL(q) :: dhdx, pa_tm, phi, dphidk
  REAL(q) :: k_x(3), k_y, k_z, k_a, cell(3,3)
  REAL(q) :: max_n, max_gr, a_iofq, b_iofq
  REAL(q) :: DCHARGSUM
  RGRID   DWORK(GRIDC%RL%NP)
  INTEGER,allocatable :: alp(:)
  INTEGER :: n1,n2,n3,nc
  INTEGER ::  i_max_n, i_max_gr

  REAL(q) :: omega_g2, omega_p2, kval, omega_0
  REAL(q) :: betaparam

  NP=GRIDC%RL%NP
  RINPL = 1._q/(real(GRIDC%NGX,kind=q)*GRIDC%NGY*GRIDC%NGZ)

  !  IF (IU0>=0) WRITE(IU0,*)  "MPLWV: ", GRIDC%MPLWV, DIMREAL(GRIDC%MPLWV), GRIDC%RC%NP
  ALLOCATE(q0(NP),dq0(NP),alp(NP),pa(NP),f_g(NP,Nalp))
  ALLOCATE(THETA(GRIDC%MPLWV,Nalp),f_k(GRIDC%MPLWV))
GPACC !$ACC ENTER DATA CREATE(Q0,DQ0,ALP,PA,F_G,THETA,F_K) COPYIN(QS_VDW,SPL,PHIK,D2PHIDK2) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(ALLOC:Q0,DQ0,ALP,PA,F_G,THETA,F_K) IF(OFFLOAD_ON)
GPOMP !$OMP TARGET ENTER DATA MAP(TO:QS_VDW,SPL,PHIK,D2PHIDK2) IF(OFFLOAD_ON)

  cell=latt_cur%B*AUTOA
  stress=0._q
GPACC !$ACC ENTER DATA COPYIN(CELL,STRESS) __IF_ASYNC__
GPOMP !$OMP TARGET ENTER DATA MAP(TO:CELL,STRESS) IF(OFFLOAD_ON)

  energ=0.0_q

  betaparam=0.03125_q * (3.0_q / (bparam**2.0_q) )**(0.75_q) * 2._q*RYTOEV

#ifdef MPI
      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE
      IDUMP=0
#ifdef debugNODE_ME
      IF (NODE_ME==IONODE) IDUMP=1
#endif
#else
      NODE_ME=0
      IONODE =0
      IDUMP=0
#ifdef debug
      IDUMP=1
#endif
#endif
  !calculate q0 in each point and it's saturated version (Soler's paper)
GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(DCHARG,DWORK4,Q0) &
GPACC !$ACC& PRIVATE(RO,GR,OMEGA_P2,OMEGA_G2,KVAL,OMEGA_0,NRS) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:DCHARG,DWORK4,Q0) &
GPOMP !$OMP PRIVATE(RO,GR,OMEGA_P2,OMEGA_G2,KVAL,OMEGA_0,NRS) IF(OFFLOAD_ON)

  DO I=1,NP
    IF (DCHARG(i).LT. 1e-11_q) THEN
      ro=1e-15_q
      gr=0._q
    ELSE
      ro=DCHARG(i)*AUTOA3
      gr=DWORK4(i)*AUTOA4
    ENDIF

    omega_p2 = 4.0_q*pi*ro
    omega_g2 = Cparam*((gr/ro)**4)
    kval=1.5_q*pi*bparam*((ro/(9.0_q*pi))**(1.0_q/6.0_q))
    omega_0 = sqrt(omega_g2 + omega_p2/3.0_q)
    q0(i) = omega_0 / kval

    !now saturated
    nrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO iu=1,12
      nrs=nrs+(q0(i)/q0cut)**iu/real(iu,kind=q)
    ENDDO
    q0(i)=q0cut*(1._q-exp(-nrs))
  ENDDO

  !calculate alp and dq0 in each point
  !that is which interpolating point we will be using and how far from it are we
  a_iofq=log((qs_vdw(Nalp)-qs_vdw(Nalp-1))/(qs_vdw(2)-qs_vdw(1)))/(Nalp-2.0_q)
  a_iofq=max(a_iofq,1e-12_q)
  b_iofq=(qs_vdw(2)-qs_vdw(1))/(exp(a_iofq)-1.0_q)
GPACC !$ACC PARALLEL PRESENT(Q0,QS_VDW,ALP,DQ0) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(IOFQ)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:Q0,QS_VDW,ALP) PRIVATE(IOFQ) IF(OFFLOAD_ON)
  DO I=1,NP
    iofq=1+log(1+(q0(i)-qs_vdw(1))/b_iofq)/a_iofq
    iofq=max(1,iofq)
    iofq=min(Nalp,iofq)
    IF(q0(i).GT.(q0cut-0.00001_q)) iofq=Nalp
    alp(i)=iofq
  ENDDO
  !this is the difference from interpolating point
GPACC !$ACC LOOP GANG VECTOR
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:Q0,QS_VDW,ALP,DQ0) IF(OFFLOAD_ON)
  DO I=1,NP
    dq0(i)=q0(i)-qs_vdw(alp(i))
  ENDDO
GPACC !$ACC END PARALLEL

  !for each interpolating point create p_\alpha (Soler's paper)
  !and Fourier transform it, this is stored in \Theta array
  DO alf=1,Nalp
GPACC !$ACC PARALLEL PRESENT(spl,alp,dq0,pa,DCHARG,dinp) __IF_ASYNC__
GPACC !$ACC LOOP GANG VECTOR PRIVATE(diff,coef)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:SPL,ALP,DQ0,PA) PRIVATE(DIFF,COEF) IF(OFFLOAD_ON)
    DO I=1,NP
      coef(:)=spl(alp(i),alf,:)
      diff=dq0(i)
      pa(i)=coef(1)+diff*(coef(2)+diff*(coef(3)+diff*coef(4)))
    ENDDO
GPACC !$ACC LOOP GANG VECTOR PRIVATE(ro,kval)
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DCHARG,DINP,PA) PRIVATE(RO,KVAL) IF(OFFLOAD_ON)
    DO I=1,NP
      IF (DCHARG(i).ge.0._q) THEN
        ro=DCHARG(i)*AUTOA3

        kval=1.5_q*pi*bparam*((ro/(9.0_q*pi))**(1.0_q/6.0_q))
        dinp(i)=pa(i)*ro/(kval**1.5_q)
      ELSE
        dinp(i)=0._q
      ENDIF
    ENDDO
GPACC !$ACC END PARALLEL

    CALL OPSYNC(DINP,CINP,GRIDC%NPLWV)
    CALL FFT3D(CINP,GRIDC,-1)

GPACC !$ACC PARALLEL LOOP PRESENT(THETA,CINP) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:THETA,CINP) IF(OFFLOAD_ON)
    DO I=1,GRIDC%RC%NP
      THETA(i,alf)=CINP(i)
    ENDDO
  ENDDO

  !now we are going to calculate the energy
  !this requires two nested loops over \alpha and \beta
  !first we do the sum over \Theta_\beta * kernel
  DO alf=1,Nalp
GPACC !$ACC KERNELS PRESENT(f_k) __IF_ASYNC__
GPOMP !$OMP TARGET MAP(PRESENT,TO:F_K) IF(OFFLOAD_ON)
IFXW8 GPOMP !$OMP TEAMS WORKDISTRIBUTE IF(OFFLOAD_ON)
    f_k=0._q
IFXW8 GPOMP !$OMP END TEAMS WORKDISTRIBUTE
GPOMP !$OMP END TARGET
GPACC !$ACC END KERNELS

GPACC !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(k_k,phik,d2phidk2,f_k,theta,gridc,cell,stress) &
GPACC !$ACC& PRIVATE(ik_k,phi,dphidk,a_k,b_k,a2_k,b2_k,a3_k,b3_k,n1,n2,n3,nc,k_x,k_a,pref) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
GPOMP !$OMP MAP(PRESENT,TO:K_K,PHIK,D2PHIDK2,F_K,THETA,GRIDC,CELL,STRESS) &
GPOMP !$OMP PRIVATE(IK_K,PHI,DPHIDK,A_K,B_K,A2_K,B2_K,A3_K,B3_K,N1,N2,N3,NC,K_X,K_A,PREF) IF(OFFLOAD_ON)

    DO bet=1,Nalp
      DO i=1,GRIDC%RC%NP
        ik_k=k_k(i)/dk_k
        IF (ik_k.GE.Nr) THEN
          !I really don't like the if here
          phi=0._q
          dphidk=0._q
        ELSE
          a_k=((ik_k+1._q)*dk_k-k_k(i))/dk_k
          b_k=1._q-a_k
          a2_k=(3._q*a_k**2.-1._q)*dk_k/6._q
          b2_k=(3._q*b_k**2.-1._q)*dk_k/6._q
          a3_k=(a_k**3.-a_k)*dk_k**2._q/6._q
          b3_k=(b_k**3.-b_k)*dk_k**2._q/6._q
          phi=a_k*phik(ik_k, bet, alf)+b_k*phik(ik_k+1,bet,alf) &
           +a3_k*d2phidk2(ik_k,bet,alf)+b3_k*d2phidk2(ik_k+1,bet,alf)
          dphidk=(-phik(ik_k, bet, alf)+phik(ik_k+1,bet,alf))/dk_k &
                 -a2_k*d2phidk2(ik_k,bet,alf)+b2_k*d2phidk2(ik_k+1,bet,alf)

        ENDIF
#ifdef OFFLOADING
        CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(f_k(i),theta(i,bet)*phi)
#else // OFFLOADING
        f_k(i)=f_k(i)+theta(i,bet)*phi
#endif // OFFLOADING
        IF (ISIF .ge. 1) THEN
           n1=Mod((i-1),gridc%rc%nrow)+1
           nc=(i-1)/gridc%rc%nrow+1
           n2=gridc%rc%i2(nc)
           n3=gridc%rc%i3(nc)
           k_x(1)=tpi*(gridc%lpctx(n1)*cell(1,1)+gridc%lpcty(n2)*cell(1,2)+gridc%lpctz(n3)*cell(1,3))
           k_x(2)=tpi*(gridc%lpctx(n1)*cell(2,1)+gridc%lpcty(n2)*cell(2,2)+gridc%lpctz(n3)*cell(2,3))
           k_x(3)=tpi*(gridc%lpctx(n1)*cell(3,1)+gridc%lpcty(n2)*cell(3,2)+gridc%lpctz(n3)*cell(3,3))
           k_a=sqrt(k_x(1)*k_x(1)+k_x(2)*k_x(2)+k_x(3)*k_x(3))  !*2._q*pi factor is here
           if (k_a .gt. 1.0e-13_q) then
              pref=dphidk/k_a*theta(i,bet)*conjg(theta(i,alf))
#ifdef NGXhalf
              if ((n1.ne.1).and.(n1.ne.(GRIDC%NGX/2+1))) pref=2._q*pref
#endif
#ifdef NGZhalf
              if ((n3.ne.1).and.(n3.ne.(GRIDC%NGZ/2+1))) pref=2._q*pref
#endif
              do ii=1,3
                 do jj=1,3
GPACC !$ACC ATOMIC UPDATE
GPOMP !$OMP ATOMIC UPDATE
                    stress(ii,jj) = stress(ii,jj) - pref*k_x(ii)*k_x(jj)
                 enddo
              enddo
           endif
        ENDIF !only calculate stress when required

      ENDDO !q points
    ENDDO  !beta

GPACC !$ACC PARALLEL LOOP PRESENT(gridc,f_k,theta) PRIVATE(nc,n1,n3,pref) REDUCTION(+:energ) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:GRIDC,F_K,THETA) &
GPOMP !$OMP PRIVATE(NC,N1,N3,PREF) REDUCTION(+:ENERG) IF(OFFLOAD_ON)

    DO i=1,GRIDC%RC%NP
      pref=real(f_k(i)*conjg(theta(i,alf)),kind=q)
#ifdef NGZhalf
!for MPI version reduced in Z direction
!and then we sum it with \Theta_\alpha
      nc=(i-1)/gridc%rc%nrow+1
      n3=gridc%rc%i3(nc)
      if ((n3.ne.1).and.(n3.ne.(GRIDC%NGZ/2+1))) pref=2._q*pref
#else
#ifdef NGXhalf
!for serial version reduced in X direction
      n1=Mod((i-1),gridc%rc%nrow)+1
      if ((n1.ne.1).and.(n1.ne.(GRIDC%NGX/2+1))) pref=2._q*pref
#endif
#endif
      energ=energ+pref
    ENDDO

    !now we fourier transform the \Theta * kernel to real space to reuse it for potential
GPACC !$ACC PARALLEL LOOP PRESENT(CINP,f_k) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:CINP,F_K) IF(OFFLOAD_ON)
    DO i=1,GRIDC%RC%NP
      CINP(I)=f_k(i)
    ENDDO

    CALL OPSYNC(DINP,CINP,GRIDC%NPLWV)
    CALL FFT3D(CINP,GRIDC,1)

GPACC !$ACC PARALLEL LOOP PRESENT(f_g,dinp) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:F_G,DINP) IF(OFFLOAD_ON)
    DO i=1,GRIDC%RL%NP
      f_g(i,alf)=dinp(i)*RINPL
    ENDDO
  ENDDO

  DCHARGSUM=0
GPACC !$ACC PARALLEL LOOP PRESENT(DCHARG) REDUCTION(+:DCHARGSUM) __IF_ASYNC__
GPOMP !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO MAP(PRESENT,TO:DCHARG) REDUCTION(+:DCHARGSUM) IF(OFFLOAD_ON)
  DO i=1,NP
    DCHARGSUM=DCHARGSUM+DCHARG(i)
  ENDDO

GPACC !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
  !add the energy to the xc energy, include volume factor
  !  number of grid points since we've done FFT, and factors to get eV
  E_corr=0.5_q*energ*LATT_CUR%OMEGA*RINPL/AUTOA3*2._q*RYTOEV + betaparam * DCHARGSUM*LATT_CUR%OMEGA
  EXC=EXC + E_corr
  !compared to vander (my local version (jik)),
  !we divide by (M)PLWV only once since it's done in xcgrad for the second time
  !so we can say that we get rid of the FFT 1/N factor and keep the \Omega/N for xcgrad
  E_corr=E_corr*RINPL

  max_n=0._q
  max_gr=0._q

GPACC !$ACC PARALLEL LOOP GANG VECTOR PRESENT(GRIDC,DCHARG,DWORK4,DWORKG,DWORK,q0,spl,alp,dq0,f_g) &
GPACC !$ACC& PRIVATE(ro,gr,dq0dn,dq0dgr,nrs,rnrs,dhdx,dtheadn,dtheadgr,coef,diff,pa_tm,dpadq0) __IF_ASYNC__

GPOMP !$OMP TARGET TEAMS DISTRIBUTE MAP(PRESENT,TO:DCHARG,DWORK4,DWORKG,DWORK,Q0,SPL,ALP,DQ0,F_G) &
GPOMP !$OMP PRIVATE(RO,GR,DQ0DN,DQ0DGR,NRS,RNRS,DHDX,DTHEADN,DTHEADGR,COEF,DIFF,PA_TM,DPADQ0) IF(OFFLOAD_ON)

  DO i=1,gridc%RL%NP
    !calc the potential
    IF (DCHARG(i).LT. 1e-11_q) THEN
      ro=1e-15_q
      gr=0._q
    ELSE
      ro=DCHARG(i)*AUTOA3
      gr=DWORK4(i)*AUTOA4
    ENDIF

    dq0dn=(-39._q*Cparam*(gr**4)+8._q*pi*(ro**5)) / (9._q*(3._q**(1._q/6._q))*bparam*(pi**(5._q/6._q))*(ro**(19._q/6._q))*sqrt(3._q*Cparam*(gr**4)+4._q*pi*(ro**5)))
    dq0dgr=4._q*Cparam*(gr**3)/((3._q**(1._q/6._q))*bparam*(pi**(5._q/6._q))*(ro**(13._q/6._q))*sqrt(3._q*Cparam*(gr**4)+4._q*pi*(ro**5)))

    nrs=0._q
    rnrs=0._q
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:RNRS) IF(OFFLOAD_ON)
    DO iu=1,12
      rnrs=rnrs+(q0(i)/q0cut)**iu/real(iu,kind=q)
    ENDDO
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:NRS) IF(OFFLOAD_ON)
    DO iu=0,11
      nrs=nrs+(q0(i)/q0cut)**iu
    ENDDO
    dhdx=nrs*exp(-rnrs)
    dtheadn=0._q
    dtheadgr=0._q
    !one loop is needed to calculate the potential
GPACC !$ACC LOOP SEQ
GPOMP !$OMP PARALLEL DO REDUCTION(+:DTHEADN,DTHEADGR) IF(OFFLOAD_ON)
    DO alf=1,Nalp
      coef=spl(alp(i),alf,:)
      diff=dq0(i)
      pa_tm=coef(1)+diff*(coef(2)+diff*(coef(3)+diff*coef(4)))
      dpadq0=coef(2)+diff*(2._q*coef(3)+3._q*diff*coef(4))
      dtheadn=dtheadn+(pa_tm/(4._q*(bparam**1.5_q)*(pi**1.25_q)*(ro**0.25_q))+(ro**0.75_q)*dpadq0*dq0dn*dhdx/(3._q*(pi**1.25_q)*(bparam**1.5_q))) *f_g(i,alf)*2._q**1.5_q
      dtheadgr=dtheadgr+ 2._q**1.5_q*f_g(i,alf) * (ro**0.75_q)*dpadq0*dhdx*dq0dgr/(3._q*(pi**1.25_q)*(bparam**1.5_q))
    ENDDO

    !add the potentials to the GGA ones
    !this is d/dn
    DWORKG(i)=DWORKG(i)+dtheadn*2._q*RYTOEV  + betaparam
    !this is d/ddn
    IF (ISPIN==1) THEN
       DWORK(i)=DWORK(i)+dtheadgr*2._q*RYTOEV/MAX(DWORK4(I),1.E-10_q)*AUTOA
    ELSE
       DWORK(i)=DWORK(i)+dtheadgr*2._q*RYTOEV*AUTOA
    ENDIF
  ENDDO !over points on real space grid


GPACC !$ACC EXIT DATA DELETE(Q0,DQ0,ALP,PA,F_G,THETA,F_K,CELL,QS_VDW,SPL,PHIK,D2PHIDK2) COPYOUT(STRESS) IF(OFFLOAD_ON) WAIT(ACC_ASYNC_Q)
GPOMP !$OMP TARGET EXIT DATA MAP(RELEASE:Q0,DQ0,ALP,PA,F_G,THETA,F_K,CELL,QS_VDW,SPL,PHIK,D2PHIDK2) IF(OFFLOAD_ON) 
GPOMP !$OMP TARGET EXIT DATA MAP(FROM:STRESS) IF(OFFLOAD_ON)
  DEALLOCATE(q0,dq0,alp,pa)
  DEALLOCATE(theta,f_k,f_g)

  stress=0.5_q*stress*LATT_CUR%OMEGA*RINPL**2._q *2._q*RYTOEV/AUTOA3

  NELECT = DCHARGSUM*LATT_CUR%OMEGA*RINPL

  CALLMPI( M_sum_d(GRIDC%COMM,E_corr,1))
  CALLMPI( M_sum_d(GRIDC%COMM,NELECT,1))

  IF (NODE_ME == IONODE) THEN
    WRITE(*,'(A,2F14.7)') 'Total rVV10 correction in eV: ', E_corr, NELECT
  ENDIF

END SUBROUTINE VDW_CALC_RVV10
