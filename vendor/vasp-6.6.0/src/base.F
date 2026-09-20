#include "symbol.inc"
!************************************************************************
! RCS:  $Id: base.F,v 1.2 2001/02/20 14:44:56 kresse Exp $
!
!> this module contains some control data structures for VASP
!
!***********************************************************************
      MODULE PREC
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
#ifndef noQuadPrecision
      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(30)
#else
      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(10)
#endif
      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
      INTEGER, PARAMETER :: qi4=SELECTED_INT_KIND(8)
      INTEGER, PARAMETER :: qi8=SELECTED_INT_KIND(15)
!
!> this parameter controls the step width in for numerical differentiation
!> in some VASP routines
!> 1E-5 is very reliable yielding at least 7 digits in all contributions
!> to the forces
!> 1E-4, however, is better suited for second derivatives
!> for reasons of consistentcy with previous versions 1E-5
      REAL(q), PARAMETER :: fd_displacement=1E-5
      END MODULE

      MODULE BASE
      USE prec
      IMPLICIT NONE
!
!> header type information
!
      TYPE PLUGINS_SETTINGS
        !> Configuration of the interface to plugins
        CHARACTER(LEN=:), ALLOCATABLE :: MODE !< are the plugins executed in serial or parallel
        LOGICAL :: FORCE_AND_STRESS = .FALSE. !< call the force and stress plugins
        LOGICAL :: LOCAL_POTENTIAL = .FALSE. !< call the local potential plugins
        LOGICAL :: FORCE_FIELD = .FALSE.  !< augment the VASP run by input from an external force field
        LOGICAL :: OCCUPANCIES = .FALSE. !< allow modification of occupation related parameters
        LOGICAL :: STRUCTURE = .FALSE. !< call the structure update plugins
        REAL(q) :: NEIGHBOR_CUTOFF = -1.0_q !< cutoff for construction of neighbor list, set to < 0 to skip
        INTEGER :: OUTBLOCK = 1 !< Step size of the
        INTEGER :: OUTPUT_MODE = 0 !< Controls whether PCDAT etc. files are calculated and written or not
        CHARACTER(LEN=:), ALLOCATABLE :: ML_MODE !< Execution mode for ML plugin
      END TYPE PLUGINS_SETTINGS

      TYPE info_struct
! only INFO
! mode information
        LOGICAL LREAL                !< real space projection/ reciprocal space proj.
        LOGICAL LOVERL               !< vanderbilt type PP read in ?
        LOGICAL LCORE                !< any partial core read in ?
        LOGICAL LCHCON               !< charge density constant during run
        LOGICAL LCHCOS               !< allways same as above in cur. impl.
        LOGICAL LONESW               !< use all bands simultaneous
        LOGICAL LONESW_AUTO          !< switch automatically between LONESW and DIIS
        LOGICAL LPRECONDH            !< precondition the subspace rotation matrix
        LOGICAL LDAVID               !< use block davidson
        LOGICAL LDAVID_FULL          !< use non-blocked davidson
        LOGICAL LEXACT_DIAG          !< use exact diagonalization
        LOGICAL LRMM                 !< use RMM-DIIS algorithm
        LOGICAL LORTHO               !< orthogonalize
        LOGICAL LHARD                !< hard stop: stop after this electronic minimization
        !> Signals that a soft stop (early graceful exit of the ionic loop) was requested.
        !>
        !> @warning Do not expect that setting this variable to .TRUE. actually stops the ionic
        !> loop. If you intend to trigger a soft stop from an internal routine use the routine
        !> TRIGGER_SOFT_STOP() and LSOFT_TRIGGER.
        LOGICAL :: LSOFT = .FALSE.
        !> Variable triggering a soft stop from an internal routine (not a STOPCAR file).
        !>
        !> Setting this to .TRUE. anywhere in the ionic loop triggers a graceful exit after the
        !> current ionic step or block of ionic steps (NBLOCK, ML_OUTBLOCK). This variable should
        !> not be set directly but rather updated via the TRIGGER_SOFT_STOP subroutine.
        LOGICAL :: LSOFT_TRIGGER = .FALSE.
        LOGICAL LABORT, LSTOP        !< stop with this iteration
        LOGICAL LPOTOK               !< local potential ok ?
        LOGICAL LMIX                 !< mixing done
        LOGICAL LASPH                !< aspherical radial PAW
! exchange correlation spin
        LOGICAL LXLDA                !< calculate LDA exchange only in POTLOK
        INTEGER ISPIN                !< spin 1=no 2 =yes
        REAL(q) RSPIN                !< 2 for spinpolar. 1 for non spinpol.
! electronic relaxation
        INTEGER ISTART               !< how to start up
        LOGICAL LHALF_INIT           !< initialize orbitals through DeePAW-HALF
        LOGICAL LHALF_API            !< obtain the HALF density from DeePAW-eSCN
        CHARACTER(255) HALF_ESCN_URL !< DeePAW-eSCN endpoint selected in INCAR
        INTEGER ICHARG               !< initial charge density
        INTEGER INIWAV               !< how to initialize wavefunctions
        INTEGER INICHG               !< how to initialize charge
        INTEGER NELM                 !< maximal number of el-steps in SCC
        INTEGER NELMALL              !< maximal number of el-steps in direct optimization
        INTEGER NELMIN               !< minimal number of el-steps
        INTEGER NELMDL               !< number of delay el-steps
        INTEGER IALGO                !< algorithm for el-relax
        INTEGER IALGO_COMPAT         !< compatibility flag for IALGO
                                     !< MOD(IALGO_COMPAT,2) = 1 use old DIIS algorithm
        LOGICAL LROBUST              !< specify how line search is performed in CG
        INTEGER ISEARCH              !< specify how line search is performed in CG
        INTEGER NDAV                 !< number of steps in RMM-DIIS
        INTEGER NKRYLOV              !< maximum size of reduced basis in the non-blocked Davidson
        INTEGER DAVITER              !< maximum number of iterations in the non-blocked Davidson
        INTEGER NSTRICT              !< the number of orbitals that will be converged tightly
        REAL(q) FBREAK               !< break-off criterion w.r.t. residual size
        REAL(q) TIME                 !< timestep for el (IALGO>50)
        REAL(q) WEIMIN,EBREAK,DEPER  !< control tags for elm
        REAL(q) EDIFF                !< accuracy for electronic relaxation
        LOGICAL LDIAG                !< level reordering allowed or not
        LOGICAL LSUBROT              !< sub space rotation to optimize rotation matrix
        LOGICAL LPDIAG               !< sub space rotation before iterat. diag.
        LOGICAL LCDIAG               !< recalculate eigenvalues after iterat. diag.
! cutoff information
        REAL(q) ENMAX                !< cutoff for calculations
        REAL(q) ENINI                !< cutoff during delay
        REAL(q) ENAUG                !< cutoff for augmentation charges
! some important things
        REAL(q) EALLAT               !< total energy of all atoms
        REAL(q) NELECT               !< number of electrons
        REAL(q) NUP_DOWN             !< spin multiplicity
        INTEGER NBANDTOT             !< total number of bands
        INTEGER MCPU                 !< max number of proc (dimensioned)
        INTEGER NCPU                 !< actual number of proc
        LOGICAL LCORR                !< correction to forces
        INTEGER IQUASI               !< eigenvalue/occupation number corrections
        INTEGER TURBO                !< turbo mode
        INTEGER IFOLD                !< folded eigenproblem: 1=LFOLD,2=LFOLDHpsi
        INTEGER IRESTART             !< whether to restart: 2=restart with 2 optimized vectors
        INTEGER IHARMONIC            !< harmonic Ritz values
        INTEGER NREBOOT              !< number of reboots
        INTEGER NMIN                 !< reboot dimension
        REAL(q) EREF                 !< reference energy to select bands
        LOGICAL NLSPLINE             !< use spline interpolation to construct projection operator
        INTEGER RADIAL_FDORDER       !< order of the derivative stencil on the radial grids
! GW algorithms
        INTEGER ICHIREAL             !< which GW algorithm , > 0 selects space time algos
        LOGICAL LGW                  !< perform GW calculations
        LOGICAL LCHI                 !< calculate response functions
        LOGICAL LscQPGW              !< iterate wavefunctions (selfconsistent GW)
        LOGICAL LGW0                 !< leave initial wavefunctions unmodified for calculations of W
        LOGICAL LG0W0                !< leave initial wavefunctions unmodified for calculations of W
        LOGICAL LCRPA                !< swithes on CRPA in GW routines
        LOGICAL LBSE                 !< solve BSE or Cassida equations (TD-DFT)
        LOGICAL LACFDT               !< ACFDT calculation (adiabatic connection, fluctuation dissipation theorem)
        LOGICAL LOEP                 !< OEP calculation: calculate the RPA-OEP potential
        LOGICAL LEXX                 !< OEP calculation: calculate the EXX-OEP potential
        LOGICAL LHFCORRECT           !< HFCORRECT calculation (adiabatic connection, fluctuation dissipation theorem)
        INTEGER EXXOEP               !< Energy cutoff settings for OEP/LHF methods
        LOGICAL LGWNO                !< GW natural orbitals (optimized for description of Coloumb hole)
        LOGICAL LCORBSE              !< solve BSE or Cassida equations and calculate correlation energy
        LOGICAL LQPBSE               !< solve BSE or Cassida equations and calculate self-energy
        LOGICAL L2E4W                !< Calculate 2-electron 4-wannier-orbital integrals
        LOGICAL L2E4W_ALL            !< all terms not just one-center terms
        TYPE(PLUGINS_SETTINGS) PLUGIN  !< which quantities may be altered by plugins
! characters allways last
        CHARACTER*40 SZNAM1          !< header of INCAR
        CHARACTER*12 SZPREC          !< precision information
      END TYPE

      !> This type controls which potential are written
      TYPE write_pot
        LOGICAL :: LVTOT = .FALSE.    !< legacy flag to write total local potential
        LOGICAL :: LVHAR = .FALSE.    !< legacy flag to write Hartree potential
        LOGICAL :: TOTAL = .FALSE.    !< write total local potential
        LOGICAL :: HARTREE = .FALSE.  !< write hartree potential
        LOGICAL :: IONIC = .FALSE.    !< write ionic potential
        LOGICAL :: XC = .FALSE.       !< write local xc potential (for LDA and GGA)
        LOGICAL :: XCMU = .FALSE.     !< write xc potential associated to kinetic energy density (for meta-GGA)
        LOGICAL :: PAW = .FALSE.      !< write paw strengths
      END TYPE write_pot

      !> This type controls which densities are written
      TYPE write_dens
        LOGICAL :: LCHARG = .TRUE.         !< write charge density (to CHGCAR or if LCHARGH5 to vaspwave.h5)
        LOGICAL :: LCHARGH5 = .FALSE.      !< write charge density to vaspwave.h5
        LOGICAL :: LTAU = .FALSE.          !< write kinetic energy density (to TAUCAR or if LCHARGH5 to vaspwave.h5)
        LOGICAL :: LPARD = .FALSE.         !< write partial charge density (to PARCHG or if LPARCHGH5 to vaspout.h5)
        LOGICAL :: LPARCHGH5 = .FALSE.     !< write partial charge densities to vaspout.h5
        LOGICAL :: LDENSITY = .FALSE.      !< write total densities (core + valence as used by the XC functional) to vaspout.h5
        LOGICAL :: LGRADIENT = .FALSE.     !< write gradient of the density (LDENSITY) to vaspout.h5
        LOGICAL :: LLAPLACIAN = .FALSE.    !< write laplacian of the density (LDENSITY) to vaspout.h5
        LOGICAL :: LWRT_AUGMENTED_DENSITY = .TRUE. !< write augmented density (if applicable) for LDENSITY/LGRADIENT/LLAPLACIAN
      END TYPE write_dens

      TYPE in_struct
! only IO
        LOGICAL LOPEN                !< files open at startup
        INTEGER IU0                  !< unit for error
        INTEGER IU6                  !< unit for stdout
        INTEGER IU5                  !< unit for stdin
        INTEGER NWRITE               !< how much information is written out
        INTEGER IDIOT                !< how much information is written out
        INTEGER ICMPLX               !< size of a complex item upon IO
        INTEGER MRECL                !< maximal size of record length
        LOGICAL LREALD               !< no LREAL read in
        LOGICAL LMUSIC               !< jF (just a joke)
        LOGICAL LFOUND               !< WAVECAR exists ?
        LOGICAL LWAVE                !< write WAVECAR
        LOGICAL LPDENS               !< write partial density (charge density for one band)
        INTEGER INITAU               !< 0: tau in file incompatible, -1: tau_m not read
        INTEGER LORBIT               !< write orbit/dos
        LOGICAL LELF                 !< write elf
        LOGICAL LOPTICS              !< calculate/write optical matrix elements
        LOGICAL LPETIM               !< timing information
        INTEGER IUVTOT               !< unit for local potential
        LOGICAL INTERACTIVE          !< vasp runs interactive
        INTEGER IRECLW               !< record lenght for WAVECAR
        LOGICAL LDOWNSAMPLE          !< read WAVECAR of denser k-grid
        LOGICAL LWAVEDERF            !< write file WAVEDERF
        LOGICAL :: DRY_RUN = .FALSE. !< execute dry run of VASP
        LOGICAL :: VELOCITY = .FALSE.!< write velocities to HDF5 file
        LOGICAL LWAP                 !< activate writting of the potential for electron-phonon calculation
#ifdef GDBDEBUG
        LOGICAL :: HALT = .FALSE.    !< wait for debugger to attach
#endif
        LOGICAL :: LH5 = .FALSE.     !< redirect restart information to vaspwave
        LOGICAL :: LSYNCH5 = .FALSE. !< synchronize the vaspout.h5 during run to allow SWMR access
        TYPE(write_pot) WRT_POTENTIAL !< which potentials are written
        TYPE(write_dens) WRT_DENSITY  !< which densities are written
      END TYPE

      TYPE mixing
! only MIX
        INTEGER IUBROY               !< unit for broyden mixer
        REAL(q) AMIX                 !< mixing parameter A
        REAL(q) BMIX                 !< mixing parameter B
        REAL(q) AMIX_MAG             !< mixing parameter A for magnetization
        REAL(q) BMIX_MAG             !< mixing parameter B for magnetization
        REAL(q) AMIN                 !< minimal mixing parameter A
        REAL(q) WC                   !< weight factor for Johnsons method
        INTEGER IMIX                 !< type of mixing
        INTEGER INIMIX               !< initial mixing matrix
        INTEGER MIXPRE               !< form of metric for mixing
        LOGICAL LRESET               !< reset mixer on next call (set when ions move)
        LOGICAL HARD_RESET           !< force hard reset of mixer (force full reset regardless of MAXMIX)
        INTEGER MAXMIX               !< maximum number of mixing steps (if positive LRESET does not apply)
        INTEGER NEIG                 !< number of eigenvalues
        INTEGER MREMOVE              !< how many vectors are removed once the iteration depth is reached
        REAL(q) EIGENVAL(512)        !< eigenvalues of dielectric matrix
        REAL(q) AMEAN                !< mean eigenvalue
        LOGICAL MIXFIRST             !< mix before diagonalization (or after)
      END TYPE

      TYPE symmetry
!only SYMM
        INTEGER, POINTER :: ROTMAP(:,:,:) !
        REAL(q), POINTER :: TAU(:,:)      ! jF
        REAL(q), POINTER :: TAUROT(:,:)   ! jF
        REAL(q), POINTER :: WRKROT(:)     ! jF
        REAL(q), POINTER :: PTRANS(:,:)   ! jF
        REAL(q), POINTER :: MAGROT(:,:)   ! jF
        INTEGER, POINTER :: INDROT(:)     ! jF
        INTEGER ISYMOP(3,3,48)            !< Space group symmetry operation
        INTEGER IGRPOP(3,3,48)            !< Rotation part of the sapce group symmetry operation
        REAL(q) GTRANS(3,48)              !< Translation part of the space group symmetry operation
        INTEGER INVMAP(48)                !< Map to the inverse elements of each group element
        REAL(q) AP(3,3)                   !< Lattice parameters of the primitive cell
        INTEGER NPCELL                    !< Number of primitive cells
        INTEGER ISYM                      !< symmetry on/of
        INTEGER NROT                      !< number of rotations
        INTEGER NPTRANS                   !< number of primitive translations
      END TYPE

      TYPE prediction
!only PRED
        INTEGER IWAVPR               !< prediction of wavefunctions
        INTEGER INIPRE               !< initialized yes/no
        INTEGER IPRE                 !< what was done in wavefunction predic.
        INTEGER IUDIR                !< unit for prediction of wavefunction
        INTEGER ICMPLX               !< size of complex word
        REAL(q) ALPHA,BETA           !< parameters of extrapolation
        REAL(q) MAXDIS               !< maximum distance; if an atom moves beyond, reset density
      END TYPE

      TYPE dipol
!only DIP
        INTEGER IDIPCO               !< direction (0 no dipol corrections)
        LOGICAL LCOR_DIP             !< correct potential
        REAL(q) POSCEN(3)            !< position of center
        REAL(q) DIPOLC(3)            !< calculated dipol
        REAL(q) QUAD                 !< trace of quadrupol
        INTEGER INDMIN(3)            !< position of minimum
        REAL(q) EDIPOL,EMONO,E_ION_EXTERN
        REAL(q), POINTER :: FORCE(:,:)
        REAL(q) VACUUM(2)            !< vacuum level
      END TYPE

      TYPE smear_struct
!only SMEAR_LOOP
        INTEGER              :: ISMCNT          !
        REAL(q), allocatable :: SMEARS(:)       !< (ismear, sigma, ...)
      END TYPE


      TYPE paco_struct
!only PACO
        INTEGER NPACO                !< number of grid points for pair corr.
        REAL(q) APACO                !< cutoff
        REAL(q),ALLOCATABLE :: SIPACO(:,:) ! accumulated partial pair correlation function, second index is flattened
        INTEGER SMEANP               !< number of configurations stored
        REAL(q) OMEGA                !< average volume times number of accumulated timesteps
      END TYPE

      TYPE energy
        REAL(q)    :: TOTENASPH  = 0.0_q  !< total energy for aspherical GGA
        REAL(q)    :: EBANDSTR   = 0.0_q  !< bandstructure energy
        REAL(q)    :: DENC       = 0.0_q  !< -1/2 hartree (d.c.)
        REAL(q)    :: XCENC      = 0.0_q  !< -V(xc)+E(xc) (d.c.)
        REAL(q)    :: EXCG       = 0.0_q  !< E(xc) (LDA+GGA)
        REAL(q)    :: EXLDA      = 0.0_q  !< LDA excchange energy
        REAL(q)    :: ECLDA      = 0.0_q  !< LDA correlation energy
        REAL(q)    :: EXGGA      = 0.0_q  !< GGA exchange energy
        REAL(q)    :: ECGGA      = 0.0_q  !< GGA correlation energy
        REAL(q)    :: EXHF       = 0.0_q  !< Hartree-Fock exchange energy
        REAL(q)    :: EXHF_ACFDT = 0.0_q  !< difference between HF energy, and exchange energy in ACFDT
        REAL(q)    :: EDOTP      = 0.0_q  !< Electric field \dot Polarization
        REAL(q)    :: TEWEN      = 0.0_q  !< Ewald energy
        REAL(q)    :: PSCENC     = 0.0_q  !< alpha Z (V(q->0) Z)
        REAL(q)    :: EENTROPY   = 0.0_q  !< Entropy term
        REAL(q)    :: PAWPS      = 0.0_q  !< paw double counting corrections
        REAL(q)    :: PAWAE      = 0.0_q  !< paw double counting corrections
        REAL(q)    :: PAWPSG     = 0.0_q  !< paw xc energies (LDA+GGA)
        REAL(q)    :: PAWAEG     = 0.0_q  !< paw xc energies (LDA+GGA)
        REAL(q)    :: PAWCORE    = 0.0_q  !< exchange correlation energy of core (LDA+GGA)
        REAL(q)    :: PAWPSAS    = 0.0_q  !< paw xc energies (aspherical)
        REAL(q)    :: PAWAEAS    = 0.0_q  !< paw xc energies (aspherical)
        COMPLEX(q) :: CVZERO     = (0.0_q, 0.0_q)  !< average local potential
        REAL(q)    :: ECGWGM     = 0.0_q  !< GW Galitskii-Migdal correlation energy
        REAL(q)    :: EKLGG0     = 0.0_q  !< Klein-contribution: Tr( G_0^-1 . G - 1 )
        REAL(q)    :: EKLLOG     = 0.0_q  !< Klein-contribution: Tr[Ln( G . G_0^-1 )]
        REAL(q)    :: ETRGHF     = 0.0_q  !< trace of HF hamiltonian with denstiy matrix
        REAL(q)    :: ELOG1G     = 0.0_q  !< Tr Log ( 1 - G S )
        REAL(q)    :: ELOGG0     = 0.0_q  !< Tr Log ( G_0 )
        REAL(q)    :: ERPA_CUT   = 0.0_q  !< RPA corr. energy with highest cutoff caculated
        REAL(q)    :: ERPA_INF   = 0.0_q  !< RPA corr. energy extrapolated to infinte basis set cutoff
        REAL(q)    :: ESCPC      = 0.0_q  !< self consistent potential correction double counting
        REAL(q)    :: EPLUGINS   = 0.0_q  !< contributions from plugins (python)
        REAL(q)    :: POTHARZERO = 0.0_q  !< g=0 component for the hartree potential
        REAL(q)    :: POTIONZERO = 0.0_q  !< g=0 component for the ion potential
      END TYPE

      CHARACTER (LEN=10) :: INCAR='INCAR'

      !> Create a pointer to complex array reinterpreting it as a real array
      INTERFACE CAST_COMPLEX_TO_REAL
         MODULE PROCEDURE C_D_CONVERT1, C_D_CONVERT2
      END INTERFACE CAST_COMPLEX_TO_REAL

      !> Create a pointer to complex array reinterpreting it as an RGRID array
      INTERFACE CAST_COMPLEX_TO_RGRID
         MODULE PROCEDURE C_RG_CONVERT1, C_RG_CONVERT2
      END INTERFACE

      !> Create a pointer to RGRID array reinterpreting it as a complex
      INTERFACE CAST_RGRID_TO_COMPLEX
         MODULE PROCEDURE RG_C_CONVERT1, RG_C_CONVERT2
      END INTERFACE

      INTERFACE RGRID2C
         MODULE PROCEDURE FRG2C, FRG2C2
      END INTERFACE

      PRIVATE RG_C_CONVERT1, C_RG_CONVERT1, C_RG_CONVERT2, C_D_CONVERT1, C_D_CONVERT2

      CONTAINS
!
!> small subroutine which tries to give good dimensions for 1 dimension
!
      SUBROUTINE MAKE_STRIDE (N)
      INTEGER N,NEW
      INTEGER, PARAMETER :: NGOOD=16

      NEW=(N+NGOOD)/NGOOD
      NEW=NEW*NGOOD+1
      N=NEW

      END SUBROUTINE

!****************************************************************************
! DESCRIPTION:
!>returns a double-precision real from a possible emulated quad precision
!>real
!> @param[in]  A  (argument)
!****************************************************************************
   FUNCTION TOREAL( A )
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      REAL(q) :: TOREAL
      QDPREAL :: A
      TOREAL = A
   END FUNCTION TOREAL

      SUBROUTINE C_D_CONVERT1(C, R)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         COMPLEX(Q), INTENT(IN), TARGET :: C(:)
         REAL(Q), INTENT(INOUT), POINTER :: R(:)

         CALL C_F_POINTER(C_LOC(C), R, SHAPE=[2*SIZE(C)])

      END SUBROUTINE C_D_CONVERT1

      SUBROUTINE C_D_CONVERT2(C, R)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         COMPLEX(Q), INTENT(IN), TARGET :: C(:,:)
         REAL(Q), INTENT(INOUT), POINTER :: R(:,:)

         CALL C_F_POINTER(C_LOC(C), R, SHAPE=[2*SIZE(C,1), SIZE(C,2)])

      END SUBROUTINE C_D_CONVERT2

      SUBROUTINE C_RG_CONVERT1(C, RG)
         COMPLEX(Q), INTENT(IN), TARGET :: C(:)
         RGRID, INTENT(INOUT), POINTER :: RG(:)
#ifdef realmode
         CALL CAST_COMPLEX_TO_REAL(C, RG)
#else
         RG => C
#endif
      END SUBROUTINE C_RG_CONVERT1

      SUBROUTINE C_RG_CONVERT2(C, RG)
         COMPLEX(Q), INTENT(IN), TARGET :: C(:,:)
         RGRID, INTENT(INOUT), POINTER :: RG(:,:)
#ifdef realmode
         CALL CAST_COMPLEX_TO_REAL(C, RG)
#else
         RG => C
#endif
      END SUBROUTINE C_RG_CONVERT2

      SUBROUTINE RG_C_CONVERT1(RG, C)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         RGRID, INTENT(IN), TARGET :: RG(:)
         COMPLEX(q), INTENT(INOUT), POINTER :: C(:)
#ifdef realmode
         CALL C_F_POINTER(C_LOC(RG), C, SHAPE=[SIZE(RG)/2])
#else
         C => RG
#endif
      END SUBROUTINE RG_C_CONVERT1

      SUBROUTINE RG_C_CONVERT2(RG, C)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         RGRID, INTENT(IN), TARGET :: RG(:,:)
         COMPLEX(q), INTENT(INOUT), POINTER :: C(:,:)
#ifdef realmode
         CALL C_F_POINTER(C_LOC(RG), C, SHAPE=[SIZE(RG,1)/2, SIZE(RG,2)])
#else
         C => RG
#endif
      END SUBROUTINE RG_C_CONVERT2

      FUNCTION FRG2C(RG) RESULT(C)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         RGRID, INTENT(IN), TARGET :: RG(:)
         COMPLEX(q), POINTER :: C(:)
#ifdef realmode
         CALL C_F_POINTER(C_LOC(RG), C, SHAPE=[SIZE(RG,1)/2])
#else
         C => RG
#endif
      END FUNCTION FRG2C

      FUNCTION FRG2C2(RG) RESULT(C)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         RGRID, INTENT(IN), TARGET :: RG(:,:)
         COMPLEX(q), POINTER :: C(:,:)
#ifdef realmode
         CALL C_F_POINTER(C_LOC(RG), C, SHAPE=[SIZE(RG,1)/2, SIZE(RG,2)])
#else
         C => RG
#endif
      END FUNCTION FRG2C2

      !*********************************************************************************************
      !> Update given soft stop trigger value with given request.
      !>
      !> Makes sure the trigger value stays .TRUE. once triggered, even if later requests are
      !> .FALSE..
      !*********************************************************************************************
      SUBROUTINE TRIGGER_SOFT_STOP(CURRENT_TRIGGER_VALUE, REQUEST)

         LOGICAL, INTENT(INOUT) :: CURRENT_TRIGGER_VALUE
         LOGICAL, INTENT(IN)    :: REQUEST

         CURRENT_TRIGGER_VALUE = CURRENT_TRIGGER_VALUE .OR. REQUEST

      END SUBROUTINE TRIGGER_SOFT_STOP

      END MODULE
