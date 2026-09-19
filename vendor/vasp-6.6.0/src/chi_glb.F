#include "symbol.inc"

!> contains the GW reader, most important data structures 
!> and most important global variables for all GW and RPA/ACFDT/BSE algorithms

MODULE chi_glb
  USE base 
  USE esf_struct_def
  USE fock_glb, ONLY: LMAX_FOCKAE
  USE poscar_struct_def, ONLY: type_info, latt
  USE wave_struct_def, ONLY: wavedes, wavedes1
  USE nonlr_struct_def, ONLY: nonlr_struct
  USE nonl_struct_def, ONLY: nonl_struct
  USE mgrid_struct_def, ONLY: grid_3d
  USE mkpoints_struct_def, ONLY: skpoints_trans
  USE mpimy
#ifdef use_shmem
  USE mpi_shmem, ONLY : M_shmem_lock,M_shmem
#endif
#ifdef _OPENMP
  USE omp_lib_kinds, ONLY : OMP_LOCK_KIND
#endif
  IMPLICIT NONE

  INTEGER,SAVE :: NBSEEIG          !< Number of BSE eigenfunctions written out in BSEEIGEN file
                                   !< starting from the lowest exciton, ie. LAMBDA=1,NBSEEIG.

!
!> use sprectral function
!
  LOGICAL, SAVE :: LSPECTRALGW=.FALSE.

!> perform GW calculations
!
  LOGICAL,SAVE :: LGW
!
!> calculate response functions
!
  LOGICAL,SAVE :: LCHI
!
!> calculate response functions in real space using G G
!> 1) using old version and folding between k-points
!> 2) using new version and entirely in real space
!
  INTEGER,SAVE :: ICHIREAL=0

!> use exact one-centere terms in HF routines for space-time GW? 
  LOGICAL, SAVE :: LFOCKSTD=.FALSE.

!
!> use old distribution scheme
!
   LOGICAL :: LBSEOLD
!
!> solve BSE or Cassida equations (TD-DFT)
!
  LOGICAL,SAVE :: LBSE
!
!> solve BSE or Cassida equations and calculate correlation energy 
!> selected using ALGO=CORBSE or CORTDHF
!
  LOGICAL,SAVE :: LCORBSE
!> solve BSE or Cassida equations and calculate self-energy 
!> selected using ALGO=QPBSE or QPRPA 
!
  LOGICAL,SAVE :: LQPBSE
!
!> perform test charge-test charge GW selfenergy evaluation
!> default is .FALSE.
!
  LOGICAL :: LTCTC=.FALSE.
!
!> use LTIME_EVOLUTION
!
  LOGICAL,SAVE :: LTIME_EVOLUTION=.FALSE.
!
!> how to solve BSE equations
!> IBSE=0   conventional matrix diagonalization (old driver)
!> IBSE=1   time evolution (new driver)
!> IBSE=2   conventional matrix diagonalization (new driver)
!> IBSE=3   Lanczos iterative matrix diagonalization (new driver)
!
  INTEGER,SAVE :: IBSE=2
!
! precision variable for IBSE=1 and IBSE=3
  CHARACTER(40) :: BSEPREC
!> IBSE=1 controls the length of the time-step
!> IBSE=3 controls the RMS for the dielectric function
!
!> use Gaussian or Lorentzian functions when ploting the dielectric function
! LOGICAL,SAVE :: LBSEGAUSS=.FALSE.
  INTEGER,SAVE :: ISMEARBSE=1
!> ISMEARBSE = 1   use Lorentzian functions
!> ISMEARBSE = 0   use Gaussian functions
!
!> hidden tag, controls whether or not the Lanczos algorithm also decomposes the Lanczos matrix into
!> eigenvalues and eigenvectors. Not documented for users
!> Default is TRUE to force the code to print compute the transition amplitudes and energies
  LOGICAL,SAVE :: LANCZOSEIGEN=.TRUE.
! 
!> number of TAU groups can be set by NTAUPAR
  INTEGER,SAVE :: NTAUPAR=-1
!> number of OMEGA groups can be set by NOMEGAPAR
  INTEGER,SAVE :: NOMEGAPAR=1
!
!> total number of frequencies
!
  INTEGER,SAVE :: NOMEGA                 
 
!
!> NOMEGA_DUMP selects which point is printed to WFULL
!
  INTEGER,SAVE   :: NOMEGA_DUMP
!
!> number of frequencies for QPBSE calculation
!
  INTEGER,SAVE :: NOMEGABSE                 
!> number of frequencies along real axis
!> the remaining frequencies are chosen to run along the
!> imaginary frequency axis
!
  INTEGER,SAVE :: NOMEGAR                
!
!> number of bands used for GW calculations
!> second entry in NBANDS lines
! 
  INTEGER,SAVE :: NBANDSEXACT
!
!> select all-in-one mode
!> 
  LOGICAL,SAVE :: LALL_IN_ONE
!
!> select specific all-in-one mode 
! 
  INTEGER,SAVE :: IALL_IN_ONE
! 
!> force writing a specific number of bands written to WAVECAR
! 
  INTEGER,SAVE :: NBANDS_WAVE
!
!> number of bands for which quasiparticle shift is evaluated
!> this is usually simply twice the number of valence bands
! 
  INTEGER,SAVE :: NBANDSGW
!
!> lowest band included in the calculation of the response function
!> defaults to 1
! 
  INTEGER,SAVE :: NBANDSGWLOW
!
!> number of bands included for the determination of the ladder diagrams
! 
  INTEGER,SAVE :: NBANDSO
  INTEGER,SAVE :: NBANDSV
!
!> calculate RPA contribution to forces
!
   LOGICAL :: LRPAFORCE=.FALSE.
!
!> use LADDER diagrams in the dielectric matrix
!
   LOGICAL :: LADDER=.FALSE.
!
!>  LFXC used DFT xc kernel
!
   LOGICAL :: LFXC=.FALSE.
!
!> use Hartree diagrams in the dielectrix matrix (presently only available in BSE)
!
   LOGICAL :: LHARTREE=.TRUE.
!
!> use test charge-test electron dielectric function
!> default is .FALSE.
!
    LOGICAL :: LTCTE=.FALSE.
!
!> use only second order terms
!> default is .FALSE.
!
    LOGICAL :: L2ORDER=.FALSE.
!
!> laplace transformed MP2
!
    LOGICAL :: LMP2LT=.FALSE.
!
!> stochastic orbitals laplace transformed MP2
!
    LOGICAL :: LSMP2LT=.FALSE.
!
!> use test electron-test electron dielectric function
!> default is true
!
    LOGICAL :: LTETE=.FALSE.
!
!> calculate triplet solution in BSE (default is singlet)
!
    LOGICAL :: LTRIPLET=.FALSE.
!
!> calculate a model fxc = epsilon-1(w=0) Xi(w=0)^-1
!
    LOGICAL :: LFXCEPS=.FALSE.
!
!> include local field effects corresponding to homogenous electron gas
!
    LOGICAL :: LFXHEG=.FALSE.
!
!> determines how natural orbitals are calculated (default rotate entire space)
!
    INTEGER :: NATURALO=-2
!
!> use holes for self-consistent GW 
!
   LOGICAL :: LHOLEGF = .FALSE.

!
!> include antiresonant term
!> BSE: 
!> ANTIRES = 0 corresponds to the Tamm Dancoff approximation, only direct terms
!> ANTIRES = 1 applies a fancy approximation that includes
!>             resonant-antiresonant coupling yielding eact  w=0 results
!> ANTIRES = 2 full version
!
!> NANOQUANTA-kernel
!> ANTIRES = 0 use only direct terms
!> ANTIRES = 1 include resonant-antiresonant exchange coupling
!>             but add it to direct coupling (exact at w=0)
!> ANTIRES = 2 include resonant-antiresonant exchange coupling
!>             exactly
    INTEGER :: ANTIRES=0

  !> dirty: to handle resonant-resonant coupling correctly the resonant part
  !> of the response function is required and stored in CHIR
  !> if this is not know/calculated a backup is implemented that
  !> gets either the exact results at w=0 or at large frequencies
    REAL(q) :: SCALE_TBSE=0.5_q  ! for exact w=0 results
  !  REAL(q) :: SCALE_TBSE=1.0_q  ! for exact high frequency results
!
!> maximum L for one center terms (per default not included)
!
    INTEGER :: LMAXMP2=-1
!
!> maximum frequency OMEGAMAX
!> up to OMEGAMAX the dielectric function is determined using
!> a rather dense grid 
!
  REAL(q),SAVE :: OMEGAMAX=0

!
!> minimum frequency OMEGAMIN
!> lowest frequency required: typically the band gap
!
  REAL(q),SAVE :: OMEGAMIN=0
!
!> maximum frequency OMEGATL
!> up to OMEGATL a rather course grid is used
!
  REAL(q),SAVE :: OMEGATL=0
!
!>  OMEGAGRID specifies the grid applied for the frequency integration
!
  INTEGER,SAVE :: OMEGAGRID=0
!
!> complex shift used in the evaluation of response functions
!> the shift is used to shift the poles slightly away from the real axis
!
  REAL(q),SAVE :: SHIFT=0
!
!> frequency dependent self energy required
!
  LOGICAL,SAVE :: LSELFENERGY=.FALSE.
!> variables determining spacing and max and min energy for the frequency-dependent self-energy 
!> calculated in CALC_SELFENERGY 
  INTEGER,SAVE :: NOMEGA_SIGMA=1001
  REAL(q),SAVE :: OMEGA_SIGMA_MIN=-150._q
  REAL(q),SAVE :: OMEGA_SIGMA_MAX= 150._q
!
!> energy cutoff for responsefunctions
!> the response function is developed in plane waves
!> and usually the cutoff ENMAX (ENCUT) is imposed 
!> the cutoff can be reduced to around 100-150 eV for
!> almost all applications without loss of precision
! 
  REAL(q),SAVE :: ENCUTGW
!
!> soft energy cutoff for Coulomb kernel GW related routines
! 
  REAL(q),SAVE :: ENCUTGWSOFT_OLD ! not used except to check for compatibility 
  REAL(q),SAVE :: ENCUTGWSOFT
!
!> TELESCOPE k-point sampling
! 
  INTEGER,SAVE :: TELESCOPE
!
!> energy cutoff for local field effects
!> I think this flag is hardly ever required
!
  REAL(q),SAVE :: ENCUTLF
!
!> scissor correction
!
  REAL(q),SAVE :: SCISSOR
!
!> k-point at which response function is determined (presently BSE only)
!> first entry specified the k-point in the KPOINT grid
!> KPOINT_BSE(2:4) specify an additional shift
!> G= KPOINT_BSE(2) b_1 + KPOINT_BSE(3) b_2  + KPOINT_BSE(4) b_3

  INTEGER,SAVE :: KPOINT_BSE(4)

!> block size for band parallelization in BSE
  INTEGER,SAVE :: NBSEBLOCKV, NBSEBLOCKO
!> position of fixed hole in fractional coordnates
  REAL(q),SAVE :: BSEHOLE(3)

!> position of fixed electron in fractional coordnates
  REAL(q),SAVE :: BSEELECTRON(3)

!> add augmenation charges in exciton WF
  LOGICAL,SAVE :: LBSEWFAUG


!> Hartree and kinetic energy contribution to eigenvalues
   COMPLEX(q), POINTER :: CELTOT_HARTREE_KINETIC(:,:,:)=>NULL()
!> Hartree and kinetic energy contribution to eigenvalues distributed
   COMPLEX(Q), POINTER :: CELEN_HARTREE_KINETIC(:,:,:)=>NULL()
!> exchange energy contribution to eigenvalues
   COMPLEX(q), POINTER :: CELTOT_X(:,:,:)=>NULL()
!> exchange energy contribution to eigenvalues distributed
   COMPLEX(q), POINTER :: CELEN_X(:,:,:)=>NULL()

!
!> subtract exact exchange from W and calculate sigma = G (W-v)
!
   LOGICAL :: LFOCK_SUBTRACT=.TRUE.
!
!> flag that allows to avoid the gamma point when the response 
!> function is determined,
!> can be used in combination with high symmetry cells
!> to sample only "odd k-points"
!> it can be used to mimic Monkhorst Pack k-point grids which are
!> otherwise not possible with GW within VASP
!> i.e. 8x8x8 k-point grid, NKRED = 2, ODDONLYGW=.TRUE.
!> results in the same results as a 4x4x4 Monkhorst Pack grid
!
   LOGICAL  :: ODDONLYGW=.FALSE.
   LOGICAL  :: EVENONLYGW=.FALSE.

!> these variable allow to reduce the grid used for the evaluation
!> of the local field effects
  INTEGER, SAVE :: NKREDLFX=1, NKREDLFY=1, NKREDLFZ=1

!
!> flag determines whether spectral functions are used when calculating 
!> response functions and screened two electron integrals
!
   LOGICAL  :: LSPECTRAL=.TRUE.

!
!> flag determines how many electronic iterations are performed
!
   INTEGER :: NELMGW=1

!
!>flag determines how many electronic iterations are performed inside the HF part
!
   INTEGER :: NELMHF=1

!> DIM in INCAR file
!> flag sets the special mode for low dimensional objects
!> DIMGW=0     0 dimensional 
!> DIMGW=1     1 dimensional line (periodic in z direction)
!> DIMGW=2     2 dimensional slab (periodic in x,y direction)
!> DIMGW=3     3 dimensional bulk
!> NOTE: not yet fully supported
!
   INTEGER :: DIMGW=3
   INTEGER :: IDIR_MAX  ! loops over direction can be restricted

!> IEPSILON determines which components of the dielectric function are calculated
!> presently only supported by time evolution code
!> IEPSILON=1-3   response to field in cartesian direction x, y or z
!> IEPSILON=4     response to field in x, y and z direction
   INTEGER :: IEPSILON=4

!> damping parameter for Newton step in the RPA routines

   REAL(q) :: DAMP_NEWTON
!
!> available memory many machines now have typically 3-4 Gbyte
!> this should be save on 3 Gbyte machines
!
   INTEGER  :: MAXMEM=2800
!
!> autoset MAXMEM, is not an INCAR tag 
   LOGICAL,SAVE :: LAUTOSET_MAXMEM = .FALSE. 

!
!> LSPECTRALCHI implies that the spectral method (Kramers-Kronig from
!>  imaginary part) is used for the evaluation of the response functions 
!
   LOGICAL  :: LSPECTRALCHI=.FALSE.
!
!> use W potentials for exchange like diagrams
!
   LOGICAL :: LGWLF=.FALSE.
!
!> iterate wavefunctions (selfconsistent GW)
!
   LOGICAL :: LscQPGW=.FALSE.
!
!> GW natural orbitals (optimized for description of Coloumb hole)
!
   LOGICAL :: LGWNO=.FALSE.
!
!> leave initial wavefunctions unmodified for calculations of W
!
   LOGICAL :: LGW0=.FALSE.
!
!> leave initial wavefunctions unmodified for calculations of W
!
   LOGICAL :: LG0W0=.FALSE.
!
!> ACFDT calculation (adiabatic connection, fluctuation dissipation theorem)
!
   LOGICAL :: LACFDT=.FALSE.
!
!> use squeezed Coulomb kernel
!
  LOGICAL :: LSCK=.FALSE.
!
!> Laplace transformed direct MP2
!
   LOGICAL :: LLTDMP2=.FALSE.
!
!> Determine linear RPA term (MP1), i.e. integrate CHI.V not (CHI.V)^2 in acfdt routines 
!
   LOGICAL :: LDMP1=.FALSE.
!
!> use Matsubara formalism in RPA and GW, i.e.
!> finite temperature formalism of many-body perturbation theory.
!> works only if ISMEAR=-1 (Fermi-smearing) is set
!
   LOGICAL :: LFINITE_TEMPERATURE=.FALSE.
!>
!> calculate Luttinger Ward grand canonical potential
!>
   LOGICAL :: LUTTINGER_WARD = .FALSE. 

!> use HF only in GWR (for testing purposes)
   LOGICAL :: LGREENHF=.FALSE.
!
!> CRPA calculations 
!
   LOGICAL :: LCRPA=.FALSE.
!
!> calculate singles contributions
!
   LOGICAL :: LSINGLES=.FALSE.
!
!> update Fermi-energy in G 
!
   LOGICAL :: LFERMIGW=.FALSE.
!
!> when the Nanoquanta kernel is used this flag determines wether 
!> f_xc = X^-1 G G W G G X^-1 is calculated (TRUE) or TBSE simply stores G G W G G
!
   LOGICAL :: LINVXI=.TRUE.
!
!> range-separated RPA (total - short range) in ACFDT 
!
   LOGICAL :: LRSRPA=.FALSE.
!
!> HFCORRECT calculation (adiabatic connection, fluctuation dissipation theorem)
!
   LOGICAL :: LHFCORRECT=.FALSE.
!
!> HFSCREEN_ORIG (simply read from INCAR)
!
   REAL(q) :: HFSCREEN_ORIG
!
!> OEP calculation: calculate the RPA-OEP potential
!
   LOGICAL :: LOEP=.FALSE.
!
!> OEP calculation: calculate the EXX-OEP potential
!
   LOGICAL :: LEXX=.FALSE.
!
!> MP2
!
   LOGICAL :: LMP2=.FALSE.
!> MP2 KPAR
   LOGICAL :: LMP2KPAR=.FALSE.
!> MP2 natual orbitals
   LOGICAL :: LMP2NO=.FALSE.
!> dumps interaction parameter for Full CI 
   LOGICAL :: LFCIDUMP=.FALSE.
!> RPA + second order 
   LOGICAL :: LRPAX=.FALSE.
!> Coupled cluster singles and doubles
   LOGICAL :: LCCSD=.FALSE.
!> 
   LOGICAL :: LBRACKETST=.FALSE.
!
!> Calculate 2-electron 4-wannier-orbital integrals
!
   LOGICAL, SAVE :: L2E4W=.FALSE.
!
!> Calculate 2-electron 4-wannier-orbital integrals
!> all terms not just one-center terms
!
   LOGICAL, SAVE :: L2E4W_ALL=.FALSE.
!
!> coupling constant strenght lambda
!> this allows to scale down the Coulomb interaction by a factor LAMBDA
!
   REAL(q), SAVE :: LAMBDA=1.0_q
!
!> use tricubic spline fit of electronic strcuture factor (ESF) in ACFDT[R][RK]
!> to reduce energy cutoff error, substitude for k-p pertubation theory
!> which uses WAVEDER and Taylor expansion around the q->0 limit of chi
!
   LOGICAL :: LESF_SPLINES=.TRUE.
!
!> same for exact exchange
   LOGICAL :: LESFX_SPLINES=.TRUE.
!
!> maximum number of blocks done at the same time
!
   INTEGER, SAVE :: NSTRIP_TOTAL=64
!
!> logical variable that decided whether parallelization is over inner 
!> k loop or other q loop
!
   LOGICAL,PARAMETER :: INNER_LOOP_PARALLEL=.TRUE.
   LOGICAL,PARAMETER :: OUTER_LOOP_PARALLEL=.NOT. INNER_LOOP_PARALLEL

!> mixed stochastic-deterministic compression of virtual states, see PRL 132, 086401 (2024)
!> https://doi.org/10.1103/PhysRevLett.132.086401
!> constant energy ratio in per cent 
   REAL(q), SAVE :: MSDGW_F = -1._q
!> constant energy ratio in per centfor occuppied states
   REAL(q), SAVE :: MSDGW_F_OCC = -1._q
!> number of stochastic subspaces at each k-point (N_s in paper)
   INTEGER, SAVE :: MSDGW_NS = -1 
!> number of protected bands above the fermi energy  (N_p in paper)
   INTEGER, SAVE :: MSDGW_NP = -1 
!> number of protected bands below the fermi energy
   INTEGER, SAVE :: MSDGW_NP_OCC = -1 
!> number of pseudo bands per stochastic subspace (N_xi in paper)
   INTEGER, SAVE :: MSDGW_NXI = 2
!> number of pseudo bands per stochastic subspace in occuppied manifold 
   INTEGER, SAVE :: MSDGW_NXI_OCC = 2
!> mersenne Twister seed for stochastic compression 
   INTEGER, SAVE :: MSDGW_SEED = 4357
!> dump WAVECAR and exit 
   LOGICAL, SAVE :: MSDGW_LWAVE = .FALSE.

! from chi_base.F
  !
  !> XI_EMPTY_THRESHHOLD sets a threshhold for emtpy orbitals
  !> 0.5 implies that the bands are divided into empty and occupied bands
  !> smaller values imply that more bands are included in the sum over
  !> occupied states and unoccupied states
  !> for metals a value of 0.001 is required to get values that are
  !> converged to 1 meV
  !> maybe we need to be even stricter for very high accuracy !?
  ! 
  REAL(q), PARAMETER :: XI_EMPTY_THRESHHOLD=0.00001
  !REAL(q), PARAMETER :: XI_EMPTY_THRESHHOLD=1E-6_q

  !> neglect orbitals that yield values smaller than -50 in the exponent
  !> this is presently only implemented for unoccupied states, since the condition
  !> is rarely met for occupied states
  REAL(q), PARAMETER :: GPOS_THRESHHOLD=50._q

  !> threshhold wether a wavevector G^2 is zero or not
  !> comparisons are done using ABS(GX**2+GY**2+GZ**2)<G2ZERO
  REAL(q), PARAMETER :: G2ZERO=1E-12_q
  !
  !> structure to store the irreducable polarizability tensor
  !
  !> @details @ref openmp :
  !> Under OpenMP the members
  !>    responsefunction::storage, responsefunction::weight,
  !>    responsefunction::ncache,  responsefunction::stpos, and
  !>    responsefunction::w_inter
  !> acquire an additional dimension.
  !> Furthermore the OpenMP version uses two additional arrays:
  !>    responsefunction::omp_lck1
  !>    responsefunction::omp_lck2
  !
  TYPE responsefunction
     INTEGER NOMEGA                      !< number of frequencies of the responsefunctions
     REAL(q)    :: SHIFT                 !< complex shift used in evaluation
     REAL(q)    :: BETA                  !< inverse temperature in units of 1/eV 
     REAL(q)    :: OMEGA0                   !< zero frequency limit for finite temperature 
     COMPLEX(q), POINTER:: COMEGA(:)     !< frequencies
     REAL(q), POINTER   :: OMEGA(:)      !< absolute frequencies
     INTEGER            :: NOMEGA_LOW    !< low frequency in global OMEGA array
     INTEGER            :: NOMEGA_HIGH   !< highest frequency in global OMEGA array
     INTEGER            :: NQ            !< q-vectors in WDES at which response function is required
     REAL(q)            :: VKPT(3)       !< q-vectors at which response function is required
     COMPLEX(q),POINTER, CONTIGUOUS :: RESPONSEFUN(:,:,:)  !< responsefunctions xi_q(G',G,w)
     REAL(q),POINTER, CONTIGUOUS    :: RESPONSER  (:,:,:)  !< real valued responsefunction (LREALSTORE)
     COMPLEX(q),POINTER, CONTIGUOUS :: RESPONSEONE(:,:)    !< one center responsefunction  xi_q(kappa,nu,w)
     REAL(q),POINTER, CONTIGUOUS    :: RESPONSEONER(:,:)   !< real valued one center responsefunction (LREALSTORE)
     INTEGER            :: NP            !< number of G vectors
     INTEGER            :: NP1           !< number of complex words along first dimension
     INTEGER            :: NP2           !< number of complex words along second dimension  
     COMPLEX(q),POINTER :: WING(:,:,:)   !< q->0 wings of responsefunctions of (:,1:3,NOMEGA) divided by 1/q
                                         !< for each cartesian direction one column 
                                         !< WING(:,IDIR,NOMEGA)=RESPONSEFUN(:,1,NOMEGA)
     COMPLEX(q),POINTER :: CWING(:,:,:)  !< q->0 wings of responsefunctions of (:,1:3,NOMEGA) divided by 1/q
                                         !< for each cartesian direction one column
                                         !< CWING(:,IDIR,NOMEGA)=RESPONSEFUN(1,:,NOMEGA)
     REAL(q),POINTER    :: WINGR(:,:,:)  !< real valued wing (LREALSTORE)
     REAL(q),POINTER    :: CWINGR(:,:,:) !< real valued wing (LREALSTORE)
     COMPLEX(q),POINTER :: HEAD(:,:,:)   !< q->0 head of responsefunctions of (1:3,1:3,NOMEGA) divided by 1/q
     LOGICAL            :: LGAMMA        !< special treatment of head and wing required
     LOGICAL            :: LREAL         !< special mode using sin consine transforms
                                         !< applicable if complex to real FFTs are used
                                         !< usually equivalent to WGW%LGAMMA
     LOGICAL            :: LREALSTORE    !< special mode using real valued response function
                                         !< reduces storage requirement by a factor 2
                                         !< but works only along imaginary axis or for spectral rep
                                         !< or for w=0
     ! shmem related
     LOGICAL            :: LSHMEM        !< shmem is used for this response function (always set)
#ifdef use_shmem
     TYPE(M_shmem_lock) :: SHMEMLOCK     !< id of semaphore for shm write operations
     TYPE(M_shmem)      :: RESPONSEFUN_shmem
#endif
     INTEGER            :: SHMID         !< id of shmem
     LOGICAL            :: LLEAD         !< lead node, all operations need to be done on this node (always set)
     TYPE(communic), POINTER :: COMM_SHMEM       => NULL() !< shmem communicator
     TYPE(communic), POINTER :: COMM_INTER_SHMEM => NULL() !< communicator between different shmem arrays

     ! cache structure 
     INTEGER            :: MAXCACHE      !< maximum number of cache lines
     COMPLEX(q),POINTER, CONTIGUOUS :: STRONE (:,:,:)              !< for one center terms
     COMPLEX(q),POINTER, CONTIGUOUS :: STORAGE(:,:,: __omp_arg(:)) !< to be added to responsefunction (:,NCACHE_DIM,NOMEGA)
     COMPLEX(q),POINTER :: WEIGHT (:,:   __omp_arg(:)) !< weight when added to responsefunction
     INTEGER,   POINTER :: NCACHE (:     __omp_arg(:)) !< how many elements are stored in cache bin (NOMEGA)
     INTEGER,   POINTER :: STPOS  (:,:,: __omp_arg(:)) !< storage postion of integral in screened two electron table
     REAL(q),   POINTER :: W_INTER(:,:,: __omp_arg(:)) !< weights for interpolation of selfenergy and derivative

!$   INTEGER(KIND=OMP_LOCK_KIND), POINTER :: OMP_LCK1(:),OMP_LCK2(:)

  END TYPE responsefunction

!> Bloch integrals of the type
!>\f[
!> \sum_{k_1} \sum_{k_2} u^*_{k_1}(r) u^*_{k_2}(r) F(r,r') u_{k_1}(r') u_{k_2}(r')
!>\f]
!>
!> can be performed
!> ) either by restricting k1 to the IRZ and looping over all k2
!> ) or by restricting q= k1-k2 to the IRZ and looping over all k1
!>   (with k2=k1-q).
!> If eigenvalues are required the first version is easier to use,
!> but the SCREENED_TWO_ELECTRON_INTEGRAL and ADD_XI use the second
!> variant.
!> The screened_2e_handle stores all required quantities to map loops
!> of the first kind to the second kind.
  TYPE screened_2e_handle
!>  number of q-points in IRZ
     INTEGER :: NUMBER_OF_NQ
!>  total number of q-points in full BZ
     INTEGER :: NUMBER_OF_NQ_FULL
!>  stores the index NK1 mapped into the IRZ
     INTEGER, POINTER :: NQ(:)
!>   K1_IN_IRZ(NQ,NK1)  stores the index NK1 mapped into the IRZ
     INTEGER, POINTER :: K1_IN_IRZ(:,:)
!>   K2_IN_IRZ(NQ,NK1)  stores the corresponding index NK2,
!>                      when NK1 is mapped into the IRZ
     INTEGER, POINTER :: K2_IN_IRZ(:,:)
!>   REQUIRED(NQ,NK1)   determines whether this integral needs
!>                      to be calculated, due to symmetry some
!>                      integrals are redundant
!>                      .T. integral needs to be calculated
!>                      .F. integral does not need to be calculated
     LOGICAL, POINTER :: REQUIRED(:,:)
!>   WTKPT(K1_IN_IRZ,K2_IN_IRZ)
!>                      stores the weight of the two electron integral
!>                      this takes care of the fact that some
!>                      integrals are symmetry equivalent
!>                      and properly weights the contributions
     REAL(q), POINTER :: WTKPT(:,:)
!>   NUMBER_OF_REDUNDANT(K1_IN_IRZ,K2_IN_IRZ)
!>                      counter for the number of in principle
!>                      symmetry inequivalent two electron integrals
!>                      usually 1 or 0
     INTEGER, POINTER :: NUMBER_OF_REDUNDANT(:,:)
!>   K2_STORE_INDEX(K1_IN_IRZ,K2_IN_IRZ)
!>                      since not all K2_IN_IRZ indices might be required
!>                      (only those with WTKPT(K1_IN_IRZ,K2_IN_IRZ)/=0
!>                       are required)
!>                      the compact index allows to use a more compact
!>                      storage mode for screened to electron integrals
!>                      it stores the actual storage position NK2_STORE
!>                      for each K1, K2 pair
     INTEGER, POINTER :: K2_STORE_INDEX(:,:)
  END TYPE screened_2e_handle


  !> handle used for spline fit of the frequency dependent screened two electron integrals along the frequency axis
  TYPE screened_2e_spline
  !> number of frequency points
     INTEGER :: NOMEGA
  !> stores the spline coefficients of the real part 
     REAL(q), POINTER :: REAL_PART(:,:)
  !> stores the spline coefficients of the imaginary part 
     REAL(q), POINTER :: IMAG_PART(:,:)
  END TYPE screened_2e_spline

  !> handle that contains information about required memory for GW jobs
  TYPE mem_gw_handle
     !> storage requirement for one Green's function at one time point:
     !> deduced from ALLOCATE_G_RECIPROCAL for G%GG, G%G_PROJ G%PROJ_PROJ
     REAL(q) :: M_GTAU = 0._q
     !> storage requirement for one Green's function at one time point:
     !> deduced from ALLOCATE_RESPONSEFUN_DISTRI 
     REAL(q) :: M_CHITAU = 0._q
     !> storage requirement for one Green's function at all frequ point:
     !> deduced from ALLOCATE_RESPONSEFUN_DISTRI 
     REAL(q) :: M_CHIOMEGA = 0._q
     !> storage requirement in RESPONSE_SUPER and co
     REAL(q) :: M_GTAU_G = 0._q
     !> storage requirement in SIGMA_SUPER and co
     REAL(q) :: M_GTAU_WTAU = 0._q
     !> storage for orbtial basis
     REAL(q) :: M_ORB = 0._q
     !> already allocated part
     REAL(q) :: M_REST = 0._q
     !> part comming from screened 2 electron integrals
     REAL(q) :: M_SCREENED2E = 0._q
     !> additionally requried wavefunctionsl
     REAL(q) :: M_WAVEFUN = 0._q
     !> total memory predicted 
     REAL(q) :: M_TOTAL = 0._q
  END TYPE mem_gw_handle

!> This descriptor contains all contributions of the full Green's function.
!> In the PAW method, we need three distinct entities to reconstruct
!> the full Green's function. That is, the plane-wave contribution
!>\f[
!>  G({\bf r}',{\bf r}) =     <{\bf r}'|a> <a| {\bf r}> , 
!>\f]
!> the plane-wave - augmenation sphere contribution 
!>\f[
!>  G({\bf r},\alpha) = <{\bf r} |a> <a| p_{\alpha} >, 
!>\f]
!> another plane-wave - augmenation sphere contribution 
!>\f[
!>  G(\alpha,{\bf r}) = <p_{\alpha} |a> <a| {\bf r} >, 
!>\f]
!> and the augmentation sphere - augmenation sphere contribution 
!>\f[
!>  G(\alpha,\beta) = <p_{\alpha} |a> <a| p_{\beta} >,
!>\f]
!> note \f$ |a>\f$ is the cell-periodic part of a orbital, \f$\alpha,\beta\f$ ionic indices and \f$ |p_{\alpha}>\f$ 
!> a projector inside the augmentation sphere of ion \f$ \alpha \f$.
  TYPE greensf
!
!> determines the data distribtution
!> 1   data distribution over second index
!
     INTEGER :: DATA_DISTRIBUTION 
     GDEF, POINTER, CONTIGUOUS :: RR(:,:)=>NULL()           !<  \f$G(r,r') \f$
     GDEF, POINTER, CONTIGUOUS :: R_PROJ(:,:)=>NULL()       !<  \f$G(r,\alpha) \f$
     GDEF, POINTER, CONTIGUOUS :: PROJ_R(:,:)=>NULL()       !<  \f$G(\alpha,r) \f$
     GDEF, POINTER, CONTIGUOUS :: PROJ_PROJ(:,:)=>NULL()    !<  \f$G(\alpha,\beta) \f$

     GDEF, POINTER, CONTIGUOUS :: GG(:,:)=>NULL()           !<  \f$ G(G,G')     \f$
     GDEF, POINTER, CONTIGUOUS :: G_PROJ(:,:)=>NULL()       !<  \f$ G(G,\alpha) \f$
     GDEF, POINTER, CONTIGUOUS :: PROJ_G(:,:)=>NULL()       !<  \f$ G(\alpha,G')\f$

     GDEF, POINTER, CONTIGUOUS :: R_G(:,:)=>NULL()          !<  \f$ G(R,G) \f$
     GDEF, POINTER, CONTIGUOUS :: G_R(:,:)=>NULL()          !<  \f$ G(G,R) \f$
  END TYPE greensf
  
!> the following data structure determine the data distribution
!> of the orbital character (CPROJ) and plane wave part
!> as well as real space part.
!> ions are distributed in a round robin distribution identical to 
!> wave_mpi.F. 
!> there are some nasty details, though. For instance
!> an FFT of a real valued function yields NRPLWV_ROW complex coefficients
!> for the Gamma only version these are interpreted as sine and cosine transforms
!> with real coefficients.
!> for the Gamma only version NRPLWV_ROW_DATA_POINTS = 2* NRPLWV_ROW
!> for complex version        NRPLWV_ROW_DATA_POINTS = NRPLWV_ROW
  TYPE greensfdes
     REAL(q)    RSPIN              !< spin multiplicity
     INTEGER :: NRSPINORS          !< number of spinors (1 for collinear, 2 for non collinear)
     LOGICAL :: LNONCOLLINEAR      !< noncollinear calculations

     INTEGER :: NRPLWV_ROW         !< total number of rows for plane wave coefficients in reciprocal space for orbitals (always dividable by cores)
     INTEGER :: NRPLWV_COL         !< number of columns for plane wave coefficients in reciprocal space for orbitals (distributed in parallel version)
     INTEGER :: NRPLWV_POS         !< global (or row) index of first column on a core in reciprocal space for orbitals 

     INTEGER :: NRPLWV_ROW_DATA_POINTS  !< number of data points in reciprocal space for orbitals along row on each node
     INTEGER :: NRPLWV_COL_DATA_POINTS  !< number of data points in reciprocal space for orbitals along column on each node
     INTEGER :: NRPLWV_COL_MAX_DATA_POINTS  !< corresponding maximum value in reciprocal space for orbitals evaluated over all nodes

     INTEGER :: RES_NRPLWV_ROW         !< total number of rows for plane wave coefficients in reciprocal space for response function (always dividable by cores)
     INTEGER :: RES_NRPLWV_COL         !< number of columns for plane wave coefficients in reciprocal space for response function(distributed in parallel version)
     INTEGER :: RES_NRPLWV_POS         !< global (or row) index of first column on a core in reciprocal space for response function
     INTEGER :: RES_NRPLWV_ROW_DATA_POINTS  !< number of data points along row on each node in reciprocal space for response function
     INTEGER :: RES_NRPLWV_COL_DATA_POINTS  !< number of data points along column on each node in reciprocal space for response function
     INTEGER :: RES_NRPLWV_COL_MAX_DATA_POINTS  !< corresponding maximum value evaluated over all nodes in reciprocal space for response function

     INTEGER :: MPLWV_ROW          !< total number of rows = total number of grid points in real space
     INTEGER :: MPLWV_COL          !< number of columns (grid points) on each core in real space
     INTEGER :: MPLWV_COL_MAX      !< maximum columns (grid points) over all cores in real space
     INTEGER :: MPLWV_POS          !< global (or row) index of first column on a core in real space

     INTEGER :: NTYP               !< data layout in column direction (round robin fashion distributed): number of types in the column direction 
     INTEGER :: NIONS              !< data layout in column direction (round robin fashion distributed): number of ions in the column direction on each core
     INTEGER,POINTER :: NIONS_CORE(:) =>NULL() !< data layout in column direction (round robin fashion distributed): number of ions on each core
     INTEGER,POINTER :: NITYP(:)      =>NULL() !< data layout in column direction (round robin fashion distributed): number of ions per type along column
     INTEGER,POINTER :: ITYP(:)       =>NULL() !< data layout in column direction (round robin fashion distributed): local type index of the ions treated locally
     INTEGER,POINTER :: NT_GLOBAL(:)  =>NULL() !< data layout in column direction (round robin fashion distributed): global index for this type
     INTEGER,POINTER :: NI_GLOBAL(:)  =>NULL() !< data layout in column direction (round robin fashion distributed): global index of each ion treated locally
     REAL(q),POINTER :: POSION(:,:)   =>NULL() !< data layout in column direction (round robin fashion distributed): ionic positions considered locally
     INTEGER,POINTER :: NT_LOCAL(:)   =>NULL() !< data layout in column direction (round robin fashion distributed): local type index of the ions treated locally

     INTEGER :: NPRO_ROW           !< total number of rows for orbital character = total number of orbital characters
     INTEGER :: NPRO_COL           !< number of columns for orbital character on each node
     INTEGER :: NPRO_COL_MAX       !< number of columns for orbital character on each node
     INTEGER,POINTER :: NPRO_LMMAX(:)     =>NULL() !< possible LM quantum numbers for each type
     INTEGER,POINTER :: NPRO_LMBASE(:)    =>NULL() !< NPRO_LMBASE(I)+1 is the local index of the first projection onto locally owned atom I
     INTEGER,POINTER :: NPRO_POS(:,:)     =>NULL() !< global (or row) index for each ion along column; attention zero based
     INTEGER,POINTER :: NPRO_ENTRIES(:,:) =>NULL() !< number of entries for each ion same as LMMAX

     !  description of density related quantities 
     ! these are usually related by Clebsch Gordon transformation to the orbital
     ! related quantities
     INTEGER :: NLM_ROW            !< total number of row indices for CG transformed augmentation charges
     INTEGER :: NLM_COL            !< total number of column indices for CG transf. charges
     INTEGER :: NLM_COL_MAX        !< maximum columns over all cores
     INTEGER,POINTER :: NLM_LMMAX(:)     =>NULL() !< total CG transformed LM quantum numbers for each type
     INTEGER,POINTER :: NLM_POS(:,:)     =>NULL() !< global (or row) index for each ion along column; attention zero based
     INTEGER,POINTER :: NLM_ENTRIES(:,:) =>NULL() !< number of entries for each ion same as NLM_LMMAX

     !> MPI communicator 
     TYPE(communic) :: COMM

     !> at some k-points inversion symmetry can be used to reduce the number of PW coefficients
     !> this is possible since G(g+k,g'+k) = G*(-g-k,-g'-k)
     !> the same applies to response functions, but for response function this property is not yet used
     !> there are routines to go from the full set of g vectors to the inversion symmetry reduced set, 
     !> where the reduced set is used for the second index [G_k(g,g') - >  G_k(g,g_red')]
     !> GDES elements LUSEINV and NGVECTOR_INV are set up by INIT_SUPERCELL_FFT and linked to entries in supercell S
     LOGICAL,POINTER :: LUSEINV(:)=>NULL()                       !< specifies that inversion symmetry can be used at this k-point
     INTEGER,POINTER :: NGVECTOR_INV(:)=>NULL()                  !< number of actual PW components if inversion symmetry is used
     ! the following arrays replicate the previous data structures, but are now k-point dependent
     INTEGER,POINTER :: NRPLWV_ROW_DATA_POINTS_NK(:)=>NULL()     !< number of data points along row on each node for each k-point
     INTEGER,POINTER :: NRPLWV_COL_DATA_POINTS_NK(:)=>NULL()     !< number of data points along column on each node for each k-point
     INTEGER,POINTER :: NRPLWV_COL_MAX_DATA_POINTS_NK(:)=>NULL() !< corresponding maximum value evaluated over all nodes for each k-point
     INTEGER,POINTER :: MAP_TO_FULL(:,:,:)=>NULL()               !< map from half-grid mode to full grid mode
!< for instance,  MAP_TO_FULL(:,1,NK) index the PW G+k has in the standard mode without inversion symmetry
!< for instance,  MAP_TO_FULL(:,2,NK) index of the PW -G-k has has in the standard mode without inversion symmetry
   END TYPE greensfdes

  !> defines a data structure for the long wave limit of the polarizability  
  TYPE head_handle
     !> head, G=G'=q=0 component of polarizability, i.e. \f$ \chi_{\bf 0 0}({\bf 0},\omega) \f$
     COMPLEX(q),POINTER :: HEAD(:,:,:)  => NULL()
     !> wing, G'=0 component of polarizability, i.e. \f$ \chi_{\bf G 0}({\bf 0},\omega) \f$
     COMPLEX(q),POINTER :: WING(:,:,:)  => NULL()
     !> wing, G=0 component of polarizability, i.e. \f$ \chi_{\bf 0 G'}({\bf 0},\omega) \f$
     COMPLEX(q),POINTER :: CWING(:,:,:) => NULL() 
     !> real valued wing (used for gamma-only)
     REAL(q),POINTER    :: WINGR(:,:,:)  => NULL()
     !> real valued cwing (used for gamma-only)
     REAL(q),POINTER    :: CWINGR(:,:,:) => NULL()
     !> frequency points of polarizability
     REAL(q), POINTER ::  OMEGA(:) => NULL() 
     !>contains information of distribution of frequency points into MPI groups 
     INTEGER, POINTER :: DISTRIBUTION(:,:) => NULL()
     !> number of complex row data points (in the response function)
     INTEGER :: NROWS
     !> number of frequency points
     INTEGER :: NOMEGA
     !> global communicator
     TYPE(communic) :: COMM
     !> communicator between frequency groups
     TYPE(communic) :: COMMBET
     !> communicator in frequency groups
     TYPE(communic) :: COMMIN      
     !> used for constrained RPA, where polarizability is calculated once for entire Fock space and 
     !> once for the correlated sub-space and the contributions are accumulated in the same instance. 
     !> However, in the latter case FACTOR is set to -1 and to subtract the correlated part
     REAL(q) :: FACTOR = 1._q
  ENDTYPE head_handle

  !> handle to store OEP related charges and potentials
  TYPE oep_handle 
  !> stores soft part of potential
     RGRID, POINTER :: SV(:,:) => NULL()
  !> stores the DIJ charges that appear in PAW
     OVERLAP, POINTER :: CDIJ(:,:,:,:) => NULL()
  !> stores the total OEP potential 
     COMPLEX(q), POINTER :: CVTOT(:,:) => NULL()

  !> this shifts the OEP at \f$ r\to \infty\f$ 
  !> to a specific value 
     REAL(q), POINTER :: ESHIFT(:) => NULL()
  END TYPE oep_handle

  TYPE correlation 
     INTEGER :: NE                                  !< number of energy points
     REAL(q), POINTER ::    ENCUTGW(:)              !< selected cutoff for which correlation energy is calculated
     REAL(q), POINTER ::    ENCUTGWSOFT(:)          !< selected cutoff for which correlation energy is calculated
     COMPLEX(q), POINTER :: CORRELATION(:)          !< total correlation energy
     COMPLEX(q), POINTER :: CORRSOSEX(:)            !< total correlation energy
     COMPLEX(q), POINTER :: CORRELATION_LAMBDA(:,:) !< lambda dependency of correlation energy
                                                    !< total correlation energy (lambda dependency)
     COMPLEX(q), POINTER :: CORRELATION_K(:)        !< correlation energy at specific k-point
     COMPLEX(q), POINTER :: CORRSOSEX_K(:)          !< correlation energy at specific k-point
     COMPLEX(q), POINTER :: CORRMP2DIR(:)           !< second order direct correlation energy
     COMPLEX(q), POINTER :: CORRMP2EX(:)            !< second order exchange correlation energy
     COMPLEX(q), POINTER :: CORRMP2DIR_K(:)         !< total correlation energy at specific k-point
     COMPLEX(q), POINTER :: CORRMP2EX_K(:)          !< total correlation energy at specific k-point

     TYPE( elstufac_handle ), POINTER   :: ESF(:)=>NULL()    !< electronic structure factors

  END TYPE correlation



CONTAINS 

!*********************************************************************
!
!> helper routine that sets the most important tags from reading the
!> ALGO tag
!
!*********************************************************************

  SUBROUTINE SET_GW_FROM_ALGO( TEXT, IU0,IU6, INFO, IO )
    USE fock_glb
    USE base 
    USE tutor, ONLY: VTUTOR, isAlert, isAdvice, ISERROR, ALGO_GW, ALGO_GW0, ALGO_G0W0, ALGO_SCGW, ALGO_SCGW0, &
                     FockStdGWR, &
#if defined( ACC_OFFLOAD ) && defined( gammareal )
                     ACFDTR_OpenACC, &
#endif
                     RPAFORCES_COMPAT, NCOREAndGWR, NonCollinearAndGWR, newArgument
    CHARACTER (LEN=*) :: TEXT
    INTEGER :: IU0, IU6
    TYPE( info_struct ) :: INFO
    TYPE( in_struct ) :: IO
    ! local 
    INTEGER :: N
    INTEGER, EXTERNAL :: LENGTH
    ! only supported codes 
    CALL STRIP(TEXT,N,'L')
    CALL LOWER(TEXT)
    N=LENGTH(TEXT)

    LGW=.FALSE.      ! GW algorithm is chosen
    LCHI=.FALSE.     ! chi is calculated
    LscQPGW=.FALSE.  ! self-consistent QP GW 
    LGW0 =.FALSE.    ! update of Green's function, Dyson equation for G is solved 
    LG0W0=.FALSE.    ! eigenvalue GW, only diagonal of self-energy is used to shift poles of G_0 
    LCRPA=.FALSE.
    LBSE=.FALSE.
    LACFDT=.FALSE.
    LOEP=.FALSE.
    LEXX=.FALSE.
    LHFCORRECT=.FALSE.
    EXXOEP=0
    LGWNO=.FALSE.
    ICHIREAL=0
    LCORBSE=.FALSE.
    LQPBSE=.FALSE.
    LLTDMP2=.FALSE.

    IF (TEXT(1:N)=='bse') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='corbse') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.TRUE.
       LCORBSE=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='qpbse') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.TRUE.
       LCORBSE=.FALSE.
       LQPBSE=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='qprpa') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.FALSE.
       LCORBSE=.FALSE.
       LQPBSE=.TRUE.
       LADDER=.FALSE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ! tdhf and cassida are synonyms
    ELSE IF (TEXT(1:MIN(N,6))=='timeev') THEN
       LCHI=.FALSE.
       LTIME_EVOLUTION=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ! tdhf and cassida are synonyms
    ELSE IF (TEXT(1:N)=='tdhf') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.FALSE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='cortdhf') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.FALSE.
       LCORBSE=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='cassida') THEN
       LCHI=.TRUE.
       LBSE=.TRUE.
       LGWLF=.FALSE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='chi') THEN
       LCHI=.TRUE.
       LGW=.FALSE.
    ELSE IF (TEXT(1:N)=='gw' .OR. TEXT(1:N)=='evgw' ) THEN
       LCHI=.TRUE.
       LGW=.TRUE.
       IF ( TEXT(1:N)=='gw' ) THEN
          CALL vtutor%write(isAlert, ALGO_GW)
       ENDIF
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='gw0' .OR. TEXT(1:N)=='evgw0' ) THEN
       LCHI=.TRUE.
       LGW =.TRUE.
       LGW0=.TRUE.
       IF ( TEXT(1:N)=='gw0' ) THEN
          CALL vtutor%write(isAlert, ALGO_GW0)
       ENDIF
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='g0w0') THEN
       LCHI=.TRUE.
       LGW =.TRUE.
       LGW0=.TRUE.
       LG0W0=.TRUE.
       IF ( TEXT(1:N)=='g0w0' ) THEN
          CALL vtutor%write(isAlert, ALGO_G0W0)
       ENDIF
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='g0w0rk'.OR. TEXT(1:N)=='evgw0rk') THEN
       LCHI=.TRUE.
       ICHIREAL=1
       LGW=.TRUE.
       LG0W0=.TRUE. ! do not solve Dyson equation for G, only pole shift in G 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='g0w0r' .OR. TEXT(1:N)=='evgw0r' )  THEN
       LCHI=.TRUE.
       ICHIREAL=2
       LGW=.TRUE.
       LG0W0=.TRUE. ! do not solve Dyson equation for G, only pole shift in G 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='gw0rk' .OR. TEXT(1:N)=='scgw0rk' ) THEN
       LCHI=.TRUE.
       ICHIREAL=1
       LGW=.TRUE.
       LGW0=.TRUE. ! Dyson equation for G is solved 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='gw0r' .OR. TEXT(1:N)=='scgw0r' ) THEN
       LCHI=.TRUE.
       ICHIREAL=2
       LGW=.TRUE.
       LGW0=.TRUE. ! Dyson equation for G is solved 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='gwrk'.OR. TEXT(1:N)=='scgwrk' ) THEN
       LCHI=.TRUE.
       ICHIREAL=1
       LGW=.TRUE. ! Dyson equation for G and W is solved 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='gwr' .OR. TEXT(1:N)=='scgwr' ) THEN
       LCHI=.TRUE.
       ICHIREAL=2
       LGW=.TRUE. ! Dyson equation for G and W is solved 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='qpgw0rk') THEN
       LCHI=.TRUE.
       ICHIREAL=1
       LGW=.TRUE.
       LG0W0=.TRUE.  ! storage wise behaves like LG0W0
       LscQPGW=.TRUE.! self-energy is hermitized 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='qpgw0r' ) THEN
       LCHI=.TRUE.
       ICHIREAL=2
       LGW=.TRUE.
       LG0W0=.TRUE.  ! storage wise behaves like LG0W0
       LscQPGW=.TRUE.! self-energy is hermitized 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='qpgw' .OR. TEXT(1:N)=='scgw' ) THEN
       LCHI=.TRUE.
       LGW =.TRUE.
       LscQPGW=.TRUE.
       IF ( TEXT(1:N)=='scgw' ) THEN
          CALL vtutor%write(isAlert, ALGO_scGW)
       ENDIF
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='qpgw0' .OR. TEXT(1:N)=='scgw0') THEN
       LCHI=.TRUE.
       LGW =.TRUE.
       LGW0=.TRUE.
       LscQPGW=.TRUE.
       IF ( TEXT(1:N)=='scgw0' ) THEN
          CALL vtutor%write(isAlert, ALGO_scGW0)
       ENDIF
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='gwno') THEN
       LCHI=.TRUE.
       LGW =.TRUE.
       LGW0=.TRUE.
       LGWNO=.TRUE.
       LscQPGW=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='crpa') THEN
       LCHI=.TRUE.
       LCRPA=.TRUE.
    ELSE IF (TEXT(1:N)=='crpar') THEN
       ICHIREAL=2
       LCHI=.TRUE.
       LCRPA=.TRUE.
       LGW=.FALSE.
       LACFDT=.FALSE.
    ELSE IF (TEXT(1:N)=='crpark') THEN
       ICHIREAL=1
       LCHI=.TRUE.
       LCRPA=.TRUE.
       LGW=.FALSE.
       LACFDT=.FALSE.
    ELSE IF (TEXT(1:N)=='ltdmp2') THEN
       LCHI=.TRUE.
       LGW=.FALSE.
       LACFDT=.TRUE.
       LLTDMP2=.TRUE.
    ELSE IF (TEXT(1:N)=='ltdmp2r') THEN
       LCHI=.TRUE.
       LGW=.FALSE.
       LACFDT=.TRUE.
       ICHIREAL=2
       LLTDMP2=.TRUE.
    ELSE IF (TEXT(1:N)=='ltdmp2rk') THEN
       LCHI=.TRUE.
       LGW=.FALSE.
       LACFDT=.TRUE.
       ICHIREAL=1
       LLTDMP2=.TRUE.
    ELSE IF (TEXT(1:N)=='rpa'.OR. TEXT(1:N)=='acfdt' .OR. TEXT(1:N)=='acdft') THEN
       LCHI=.TRUE.
       LGW=.FALSE.
       LACFDT=.TRUE.
    ELSE IF (TEXT(1:N)=='rpar'.OR. TEXT(1:N)=='acfdtr' .OR. TEXT(1:N)=='acdftr') THEN
       LCHI=.TRUE.
       ICHIREAL=2
       LGW=.FALSE.
       LACFDT=.TRUE.
#if defined(ACC_OFFLOAD) && defined(gammareal)
       ! warn user that code path is not ported for gamma only 
       CALL vtutor%write(isAlert, ACFDTR_OpenACC )
#endif 
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='rpark'.OR. TEXT(1:N)=='acfdtrk' .OR. TEXT(1:N)=='acdftrk') THEN
       LCHI=.TRUE.
       ICHIREAL=1
       LGW=.FALSE.
       LACFDT=.TRUE.
       IF( .NOT. LFOCKSTD ) CALL vtutor%write(isAdvice, FockStdGWR )
    ELSE IF (TEXT(1:N)=='hfc') THEN
       LCHI=.TRUE.
       LGW=.FALSE.
       LACFDT=.TRUE.
       LHFCORRECT=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='oep') THEN
       LCHI=.TRUE.
       LGW=.TRUE.
       LscQPGW=.TRUE.
       LOEP=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='oeprk') THEN
       LCHI=.TRUE.
       LG0W0=.TRUE.
       LGW=.TRUE.
       LOEP=.TRUE.
       ICHIREAL=1
    ELSE IF (TEXT(1:N)=='oepr') THEN
       LCHI=.TRUE.
       LG0W0=.TRUE.
       LGW=.TRUE.
       LOEP=.TRUE.
       ICHIREAL=2
    ELSE IF (TEXT(1:N)=='exxr') THEN
       LCHI=.TRUE.
       LG0W0=.TRUE.
       LGW=.TRUE.
       LOEP=.TRUE.
       LEXX=.TRUE.
       ICHIREAL=2
    ELSE IF (TEXT(1:N)=='exxrk') THEN
       LCHI=.TRUE.
       LG0W0=.TRUE.
       LGW=.TRUE.
       LOEP=.TRUE.
       LEXX=.TRUE.
       ICHIREAL=1
    ELSE IF (TEXT(1:N)=='exx') THEN
       LCHI=.TRUE.
       LGW=.TRUE.
       LscQPGW=.TRUE.
       LOEP=.TRUE.
       LEXX=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='mp2') THEN
       LCHI=.TRUE.
       LMP2=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='mp2kpar') THEN
       LCHI=.TRUE.
       LMP2KPAR=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='mp2no') THEN
       LCHI=.TRUE.
       LMP2NO=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='rpax') THEN
       LCHI=.TRUE.
       LRPAX=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='ccsd') THEN
       LCHI=.TRUE.
       LCCSD=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='(t)') THEN
       LCHI=.TRUE.
       LBRACKETST=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='fcidump') THEN
       LCHI=.TRUE.
       LFCIDUMP=.TRUE.
       IF (LFINITE_TEMPERATURE) GOTO 99    ! not supported yet 
    ELSE IF (TEXT(1:N)=='2e4w') THEN
       LCHI=.TRUE.
       L2E4W=.TRUE.
    ELSE IF (TEXT(1:5)=='2e4wa') THEN
       LCHI=.TRUE.
       L2E4W=.TRUE.
       L2E4W_ALL=.TRUE.
    ELSE
      ! mark TEXT that no GW algorithm was selected
       TEXT='NOT_A_GW_ALGORITHM'
    ENDIF

    ! concerning RPA forces: 
    IF (LRPAFORCE) THEN
       IF (LACFDT) THEN
          ! Green's functions are required, so select G0W0 as well 
          LGW=.TRUE.
          LG0W0=.TRUE.
          LACFDT=.FALSE.
       ENDIF 
    ENDIF 
   
! store info about GW algorithms into INFO
! this is not really used, probably best to get rid of it 
    INFO%LGW       = LGW       
    INFO%LCHI      = LCHI      
    INFO%LscQPGW   = LscQPGW   
    INFO%LGW0      = LGW0      
    INFO%LG0W0     = LG0W0     
    INFO%LCRPA     = LCRPA     
    INFO%LBSE      = LBSE      
    INFO%LACFDT    = LACFDT    
    INFO%LOEP      = LOEP      
    INFO%LEXX      = LEXX      
    INFO%LHFCORRECT= LHFCORRECT
    INFO%EXXOEP    = EXXOEP    
    INFO%LGWNO     = LGWNO     
    INFO%ICHIREAL  = ICHIREAL  
    INFO%LCORBSE   = LCORBSE   
    INFO%LQPBSE    = LQPBSE    
    INFO%L2E4W     = L2E4W
    INFO%L2E4W_ALL = L2E4W_ALL

    RETURN 

    ! might not be supported 
99  CONTINUE 
    CALL VTUTOR%ERROR("LFINITE_TEMPERATURE supported only for ALGO = &
     & ACFDT[R][RK] | RPA[R][RK] | [QP]GW0[R|RK] | [QP]GW[R|RK]")

  END SUBROUTINE SET_GW_FROM_ALGO

!> helper routine used to dump read ALGO tag to OUTCAR
  FUNCTION ALGO_FROM_GW()
    CHARACTER (LEN=7) :: ALGO_FROM_GW

    IF (L2E4W_ALL) THEN
       ALGO_FROM_GW='2E4WA '
    ELSE IF (L2E4W) THEN
       ALGO_FROM_GW='2E4W  '
    ELSE IF (LMP2) THEN
       ALGO_FROM_GW='MP2   '
    ELSE IF (LMP2KPAR) THEN
       ALGO_FROM_GW='MP2KPAR'
    ELSE IF (LMP2NO) THEN
       ALGO_FROM_GW='MP2NO '
    ELSE IF (LRPAX) THEN
       ALGO_FROM_GW='RPAX  '
    ELSE IF (LCCSD) THEN
       ALGO_FROM_GW='CCSD  '
    ELSE IF (LBRACKETST) THEN
       ALGO_FROM_GW='(T)  '
    ELSE IF (LFCIDUMP) THEN
       ALGO_FROM_GW='FCIDUMP'
    ELSE IF (LHFCORRECT) THEN
       ALGO_FROM_GW='HFC'
    ELSE IF (LACFDT) THEN
       ALGO_FROM_GW='ACFDT '
    ELSE IF (LG0W0.AND. LscQPGW .AND. ICHIREAL>=2) THEN
       ALGO_FROM_GW='QPGW0R'
    ELSE IF (LG0W0.AND. LscQPGW .AND. ICHIREAL>=1) THEN
       ALGO_FROM_GW='QPGW0RK'
    ELSE IF (LG0W0.AND. ICHIREAL>=2) THEN
       ALGO_FROM_GW='EVG0W0R'
    ELSE IF (LG0W0.AND. ICHIREAL>=1) THEN
       ALGO_FROM_GW='EVG0W0RK'
    ELSE IF (LGW0.AND. ICHIREAL>=2) THEN
       ALGO_FROM_GW='EVGW0R'
    ELSE IF (LGW0.AND. ICHIREAL>=1) THEN
       ALGO_FROM_GW='EVGW0RK'
    ELSE IF (LACFDT.AND. ICHIREAL>=2) THEN
       ALGO_FROM_GW='ACFDTRK'
    ELSE IF (LACFDT.AND. ICHIREAL>=1) THEN
       ALGO_FROM_GW='ACFDTR'
    ELSE IF (ICHIREAL>=2) THEN
       ALGO_FROM_GW='scGWR'
    ELSE IF (ICHIREAL>=1) THEN
       ALGO_FROM_GW='scGWRK'
    ELSE IF (LTIME_EVOLUTION) THEN
       ALGO_FROM_GW='TIMEEV'
    ELSE IF (LBSE) THEN
       IF (LCORBSE) THEN
          IF (LGWLF) THEN
             ALGO_FROM_GW='CORBSE'
          ELSE
             ALGO_FROM_GW='CORTDHF'
          ENDIF
       ELSE
          IF (LGWLF) THEN
             ALGO_FROM_GW='BSE   '
          ELSE
             ALGO_FROM_GW='TDHF  '
          ENDIF
       ENDIF
       IF (LQPBSE) THEN
          IF (LGWLF) THEN
             ALGO_FROM_GW='QPBSE '
          ELSE
             ALGO_FROM_GW='QPRPA '
          ENDIF
       ENDIF
    ELSE IF (LEXX) THEN
       ALGO_FROM_GW='EXX   '
    ELSE IF (LOEP) THEN
       ALGO_FROM_GW='OEP   '
    ELSE IF (LEXX.AND.ICHIREAL>=1) THEN
       ALGO_FROM_GW='EXXRK '
    ELSE IF (LEXX.AND.ICHIREAL>=2) THEN
       ALGO_FROM_GW='EXXR  '
    ELSE IF (LOEP.AND.ICHIREAL>=1) THEN
       ALGO_FROM_GW='OEPRK '
    ELSE IF (LOEP.AND.ICHIREAL>=2) THEN
       ALGO_FROM_GW='OEPR  '
    ELSE IF (LGWNO .AND. LscQPGW .AND. LGW0) THEN
       ALGO_FROM_GW='GWNO  '
    ELSE IF (LscQPGW .AND. LGW0) THEN
       ALGO_FROM_GW='QPGW0 '
    ELSE IF (LCRPA) THEN
       ALGO_FROM_GW='CRPA '
    ELSE IF (LCRPA) THEN
       ALGO_FROM_GW='CRPAR '
    ELSE IF (LscQPGW) THEN
       ALGO_FROM_GW='QPGW  '
    ELSE IF (LG0W0) THEN
       ALGO_FROM_GW='G0W0  '
    ELSE IF (LGW0) THEN
       ALGO_FROM_GW='EVGW0 '
    ELSE IF (LGW) THEN
       ALGO_FROM_GW='EVGW  '
    ELSE IF (LCHI) THEN
       ALGO_FROM_GW='CHI   '
    ELSE
       ALGO_FROM_GW='none  '
    ENDIF
  END FUNCTION ALGO_FROM_GW

  SUBROUTINE MAP_SELF_ENERGY(SIGMAU_TAU, SIGMAO_TAU, IMAG_GRIDS, NBANDSGW, SIGMA)
     USE minimax_struct , ONLY :imag_grid_handle
     IMPLICIT NONE
     GDEF, INTENT(IN)                   :: SIGMAU_TAU(:,:)  ! negative time self-energy
     GDEF, INTENT(IN)                   :: SIGMAO_TAU(:,:)  ! positive time self-energy
     TYPE(imag_grid_handle), INTENT(IN) :: IMAG_GRIDS        !time and frequency grid handle
     INTEGER, INTENT(IN)                :: NBANDSGW        !number of bands to be written to file
     COMPLEX(q), INTENT(INOUT)          :: SIGMA(:,:)   !self-energy on imaginary frequency axis 
     !local variables 
     REAL(q), ALLOCATABLE :: RE(:,:)
     REAL(q), ALLOCATABLE :: IM(:,:)

     PROFILING_START( 'map_self_energy' )

     ALLOCATE(RE(IMAG_GRIDS%NOMEGA, NBANDSGW))
     ALLOCATE(IM(IMAG_GRIDS%NOMEGA, NBANDSGW))
  
     ! build real part of self-energy
     SIGMA = 0

     IM( :, 1:NBANDSGW ) =   REAL( SIGMAO_TAU( :, 1:NBANDSGW) - SIGMAU_TAU( :, 1:NBANDSGW),q)
     CALL DGEMM('N','N', IMAG_GRIDS%NOMEGA, NBANDSGW, IMAG_GRIDS%NOMEGA, one, &
          IMAG_GRIDS%TO_FER_RE(1,1), IMAG_GRIDS%NOMEGA, IM(1,1), IMAG_GRIDS%NOMEGA, & 
          zero, RE(1,1), IMAG_GRIDS%NOMEGA )
 
     SIGMA = RE

     ! build imaginary part of self-energy
     RE( :, 1:NBANDSGW) = -REAL( SIGMAU_TAU( :, 1:NBANDSGW) + SIGMAO_TAU( :, 1:NBANDSGW),q)
     CALL DGEMM('N','N', IMAG_GRIDS%NOMEGA, NBANDSGW, IMAG_GRIDS%NOMEGA, one, &
          IMAG_GRIDS%TO_FER_RE_CONJG(1,1), IMAG_GRIDS%NOMEGA, RE(1,1), IMAG_GRIDS%NOMEGA, & 
          zero, IM(1,1), IMAG_GRIDS%NOMEGA )
  
     SIGMA = SIGMA + CMPLX(0._q,1._q)*IM 

     DEALLOCATE(RE, IM)
     PROFILING_STOP( 'map_self_energy' )
  END SUBROUTINE MAP_SELF_ENERGY


END MODULE chi_glb



!**********************************************************************
!
!> reads all variables related to the responsefunctions
!
!**********************************************************************

  SUBROUTINE RESPONSE_READER( INFO, IO )
   
    USE base
    USE reader_tags
    USE chi_glb
    USE string, ONLY: str
    USE tutor, ONLY: vtutor, isError, isAlert, NOMEGA_DUMPtooLarge, LSPECTRALGWset,&
       ImagGridParForLTMP2, SinglesNaturalo, LuttingerWardNaturalo, newArgument
    USE mkpoints_struct_def, ONLY: LSHIFT_KPOINTS

    IMPLICIT NONE
    TYPE (info_struct) INFO
    TYPE( in_struct ) :: IO
! local
    INTEGER IU5, IU6, IU0
    INTEGER IDUM, N, IERR
    REAL(q) RDUM
    COMPLEX(q) CDUM
    LOGICAL LOPEN, LDUM
    CHARACTER (1) :: CHARAC
    CHARACTER (40) :: TEXT

    IU5=IO%IU5
    IU6=IO%IU6
    IU0=IO%IU0

    CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
! calculate RPA forces

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LRPAFORCE', LRPAFORCE, IERR, WRITEXMLINCAR)
    IF ( LRPAFORCE .AND. INFO%LREAL ) THEN
       CALL vtutor%error( "LRPAFORCE=T and LREAL=Auto (True) is currently not implemented")
    ENDIF
    ! switch on finite temperature in RPA and GW routines
    ! this forces all (fermionic) bosonic correlation functions to have a 
    ! (anti) period of 1/(k_b T) = beta 
    LFINITE_TEMPERATURE=.FALSE.
    ! anyone sane will go beyond RPA for BSE
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFINITE_TEMPERATURE', LFINITE_TEMPERATURE, IERR, WRITEXMLINCAR )

    LGREENHF=.FALSE.
    ! use only Hartree-Fock approximation in GWR algos (used for testing purposes)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LGREENHF', LGREENHF, IERR, WRITEXMLINCAR )

    !> use one-centre terms in HF routines for euclidean GW/RPA 
    LFOCKSTD=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFOCKSTD', LFOCKSTD, IERR, WRITEXMLINCAR)

! algorithm: ALGO tag
    TEXT="--"
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ALGO', TEXT, 40, IERR, WRITEXMLINCAR)
    CALL SET_GW_FROM_ALGO( TEXT, IU0, IU6, INFO, IO )

    ! compute two-centre Coulomb integrals 
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LTWO_CENTRE', L2E4W_ALL, IERR, WRITEXMLINCAR )

! read number of frequencies NOMEGA
    NOMEGA=0
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NOMEGA', NOMEGA, IERR, WRITEXMLINCAR)
    ! space time routines require decent number of imaginary frequency and time points 
    IF ( NOMEGA<=0 .AND. ( LACFDT .OR. LCRPA .OR. ICHIREAL>=1 ) ) THEN
       NOMEGA=12
    ENDIF
    ! exact exchange OEP useless to use more than one frequency
    ! and even that is already too much :)
    IF ((LEXX .AND. ICHIREAL <1) .OR. ( LCRPA .AND. ICHIREAL < 1 ) ) THEN
       NOMEGA=1
    ENDIF
    ! exact exchange in euclidean QFT formalism
    IF (NOMEGA<=0 .AND. LCHI .AND. .NOT. LBSE .AND. ICHIREAL < 1 ) THEN
       NOMEGA=100
    ENDIF
    !
    ! QPBSE, set the number of frequency points to evaluate self-energy
    ! 
    IF (LBSE .AND. LQPBSE) THEN
       NOMEGABSE=NOMEGA
       NOMEGA=0
    ENDIF
    !
    ! BSE, no need for frequency grid (just eats up storage possibly)
    ! 
    IF (LBSE) THEN
       NOMEGA=0
    ENDIF

    LDMP1=.FALSE.
    ! compute linear term in acfdt routine instead of direct MP2 term (which is
    ! of second order )
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDMP1', LDMP1, IERR, WRITEXMLINCAR)

    ! the distribution of tau points among the CPUs is governed by the NTAUPAR flag
    ! NTAUPAR determines how many TAU groups should be created
    NTAUPAR=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NTAUPAR', NTAUPAR, IERR, WRITEXMLINCAR)
    ! the distribution of omega points among the CPUs is governed by the NOMEGAPAR flag
    ! NOMEGAPAR determines how many OMEGA groups should be created
    ! default: same number of tau and omega groups
    NOMEGAPAR=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NOMEGAPAR', NOMEGAPAR, IERR, WRITEXMLINCAR)

    ! NOMEGA_DUMP selects global frequency point of screened potential 
    ! which is written to WFULL???? files for each q-point in the 1.BZ 
    ! negative value means WFULL files are not written 
    ! NOMEGA_DUMP = 0 means W(omega=0) is written
    NOMEGA_DUMP = -1 
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NOMEGA_DUMP', NOMEGA_DUMP, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    IF ( NOMEGA_DUMP > NOMEGA ) THEN
       CALL vtutor%write(isAlert, NOMEGA_DUMPtooLarge)
       NOMEGA_DUMP = NOMEGA 
    ENDIF
    !> set NOMEGA_DUMP = 0 for cRPA 
    IF( N < 1 .AND. LCRPA ) THEN
       NOMEGA_DUMP = 0 
    ENDIF
    !> finite temperature grid contains omega=0 point,
    !> if user sets NOMEGA_DUMP = 0 no need to integrate over tau
    IF( LFINITE_TEMPERATURE .AND. NOMEGA_DUMP == 0  ) THEN
       NOMEGA_DUMP = MAX( NOMEGA_DUMP, 1 ) 
    ENDIF 
! use IP approximation in BSE (default is .FALSE. here)
! this is for testing only
    LHARTREE=.TRUE.
! anyone sane will go beyond RPA for BSE
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LHARTREE', LHARTREE, IERR, WRITEXMLINCAR)

! use LADDER diagrams (default is .FALSE.)
    LADDER=.FALSE.
! anyone sane will go beyond RPA for BSE
    IF (LBSE) LADDER=.TRUE.
    LADDER=.NOT. LADDER
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LRPA', LADDER, IERR, WRITEXMLINCAR)
    LADDER=.NOT. LADDER

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LADDER', LADDER, IERR, WRITEXMLINCAR)

! use approximate ladder diagrams (default is .FALSE.)
    LFXC=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFXC', LFXC, IERR, WRITEXMLINCAR)

! use RSRPA (default is .FALSE. here)

    LRSRPA=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LRSRPA', LRSRPA, IERR, WRITEXMLINCAR)

! calculate HF singles
    LSINGLES=.FALSE.

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSINGLES', LSINGLES, IERR, WRITEXMLINCAR)

! iterate Fermi-level
    LFERMIGW=.FALSE.

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFERMIGW', LFERMIGW, IERR, WRITEXMLINCAR)

!
! IBSE flag (how to do BSE)
!
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IBSE', IBSE, IERR, WRITEXMLINCAR)
! BSEPREC flag for precision in IBSE=2 and IBSE=3
    BSEPREC = 'NORMAL'
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BSEPREC', BSEPREC, 40, IERR, WRITEXMLINCAR)
! LBSEGAUSS option for ploting the dielectric function
!   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LBSEGAUSS', LBSEGAUSS, IERR, WRITEXMLINCAR)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISMEARBSE', ISMEARBSE, IERR, WRITEXMLINCAR)
! LANCZOSEIGEN option for ploting the dielectric function
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LANCZOSEIGEN', LANCZOSEIGEN, IERR, WRITEXMLINCAR)    
!

! use LTCTC (default is .FALSE.)

    LTCTC=.FALSE.

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LTCTC', LTCTC, IERR, WRITEXMLINCAR)

! read in the k-point at which BSE should be done
    KPOINT_BSE = 0 ; KPOINT_BSE(1)=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPOINT_BSE', KPOINT_BSE, 4, IERR, WRITEXMLINCAR,  lcontinue=.true.)
    IF (IERR/=0) CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPOINT_BSE', KPOINT_BSE, 1, IERR, WRITEXMLINCAR)

! NBSEBLOCKO and NBSEBLOCKV blocking valence and conduction bands

    NBSEBLOCKV=-1; NBSEBLOCKO=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBSEBLOCKO', NBSEBLOCKO, IERR, WRITEXMLINCAR)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBSEBLOCKV', NBSEBLOCKV, IERR, WRITEXMLINCAR)


! switch to old BSE driver
    LBSEOLD=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LBSEOLD', LBSEOLD, IERR, WRITEXMLINCAR)
! read in the position at which hole should be fixed
    BSEHOLE = -123456.0_q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BSEHOLE', BSEHOLE, 3, IERR, WRITEXMLINCAR)

! read in the position at which electron should be fixed
    BSEELECTRON = -123456.0_q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BSEELECTRON', BSEELECTRON, 3, IERR, WRITEXMLINCAR)

! use augmenation charges in exciton WF
    LBSEWFAUG=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LBSEWFAUG', LBSEWFAUG, IERR, WRITEXMLINCAR)

! NBSEEIG number of BSE wave functions written

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBSEEIG', NBSEEIG, IERR, WRITEXMLINCAR)

! use LTCTE (default is .FALSE.)

    LTCTE=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LTCTE', LTCTE, IERR, WRITEXMLINCAR)

! use LTETE (default is .TRUE.)

    LTETE=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LTETE', LTETE, IERR, WRITEXMLINCAR)

! use LTRIPLET (default is .TRUE.)

    LTRIPLET=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LTRIPLET', LTRIPLET, IERR, WRITEXMLINCAR)

! use LFXCEPS (default is .TRUE.)

    LFXCEPS=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFXCEPS', LFXCEPS, IERR, WRITEXMLINCAR)

! use LFXHEG (default is .TRUE.)

    LFXHEG=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFXHEG', LFXHEG, IERR, WRITEXMLINCAR)

! compute Luttinger Ward grand potential
    LUTTINGER_WARD=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LUTTINGER_WARD', LUTTINGER_WARD, IERR, WRITEXMLINCAR )

! use NATURALO (default is 2, i.e. rotation in unoccupied space)

    NATURALO=-2
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NATURALO', NATURALO, IERR, WRITEXMLINCAR, FOUNDNUMBER=N) 
    IF ( N < 1 ) THEN
       ! use natural orbitals at T>0 or OEP by default
       IF( (LFINITE_TEMPERATURE .AND. LUTTINGER_WARD) .OR. LOEP ) THEN
          NATURALO = 0 
       ! use DFT orbitals for EVGW0[R|RK] 
       ELSE IF( LG0W0 .AND. .NOT. LGW0 ) THEN
          NATURALO = 4 
       ! otherwise use NO only for occ-unoccupied block for T=0 
       ELSE
          NATURALO = 2 
       ENDIF
    ENDIF

    !> check consistency with LSINGLES 
    IF ( NATURALO /= 2 .AND. LSINGLES ) THEN
       CALL vtutor%write(isAlert, SinglesNaturalo, newArgument( IVAL = [ NATURALO ] ) )
    ENDIF

    !> check consistency with LUTTINGER_WARD
    IF ( NATURALO /= 0 .AND. LUTTINGER_WARD ) THEN
       CALL vtutor%write(isError, LuttingerWardNaturalo, newArgument( IVAL = [ NATURALO ] ) )
    ENDIF


! use LHOLEGF (default is .FALSE.)

    LHOLEGF=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LHOLEGF', LHOLEGF, IERR, WRITEXMLINCAR)

! use L2ORDER (default is .FALSE.)

    L2ORDER=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'L2ORDER', L2ORDER, IERR, WRITEXMLINCAR)
    IF (L2ORDER) THEN
       LTETE=.FALSE.
       LTCTE=.FALSE.
    ENDIF

! use LMP2LT (default is .FALSE.)

    LMP2LT=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMP2LT', LMP2LT, IERR, WRITEXMLINCAR)

! use LSMP2LT (default is .FALSE.)

    LSMP2LT=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSMP2LT', LSMP2LT, IERR, WRITEXMLINCAR)

    ! safety measure , there is no omega-parallelization available for this method
    IF ( LMP2LT .OR. LSMP2LT ) THEN
       IF ( NOMEGAPAR > 1 .OR. NTAUPAR > 1 ) THEN
           CALL vtutor%write(isAlert, ImagGridParForLTMP2)
       ENDIF
       NOMEGAPAR = 1
       NTAUPAR = 1
    ENDIF

! solve BSE beyond TDA; TDA: ANTIRES=0  full BSE: ANTIRES=2
    ANTIRES=0
    IF (LQPBSE) ANTIRES=2
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ANTIRES', ANTIRES, IERR, WRITEXMLINCAR)
    IF (ANTIRES==2 .AND. IBSE==3) THEN
      CALL VTUTOR%ERROR("Lanczos algorithm is not implemented to go beyond the Tamm-Dancoff approximation (TDA), i.e. ANTIRES=2. If you &
           &want to go beyond TDA then set IBSE=0, 1, or 2. If you wish use the Lanczos diagonalization, please set ANTIRES=0.")
    ENDIF
! number of frequencies along real axis
    NOMEGAR=0
    IF (LCHI) NOMEGAR=NOMEGA

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NOMEGAR', NOMEGAR, IERR, WRITEXMLINCAR)

    ! force NOMEGAR to be zero for ACFDT
    IF (LACFDT) NOMEGAR=0
! read entry NBANDSEXACT
    NBANDSEXACT=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSEXACT', NBANDSEXACT, IERR, WRITEXMLINCAR)

! switch on all-in-one mode 
! can be used with NBANDS set in INCAR 
    LALL_IN_ONE=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LALL_IN_ONE', LALL_IN_ONE, IERR, WRITEXMLINCAR)

! switch on all-in-one mode 
    IALL_IN_ONE=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IALL_IN_ONE', IALL_IN_ONE, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    ! select default mode 
    IF ( LALL_IN_ONE .AND. N<1 ) THEN
       IALL_IN_ONE = 1
    ELSE IF ( IALL_IN_ONE > 0 ) THEN
    ! overwrite LALL_IN_ONE 
       LALL_IN_ONE = .TRUE. 
    ENDIF

! switch on all-in-one mode 
    NBANDS_WAVE=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDS_WAVE', NBANDS_WAVE, IERR, WRITEXMLINCAR)

! read number of bands for which the shifts are calculated
    NBANDSGW=-1

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSGW', NBANDSGW, IERR, WRITEXMLINCAR)

! read number of bands for which the shifts are calculated
    NBANDSGWLOW=1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSLOW', NBANDSGWLOW, IERR, WRITEXMLINCAR)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSGWLOW', NBANDSGWLOW, IERR, WRITEXMLINCAR)
    NBANDSGWLOW=MAX(NBANDSGWLOW,1)
!
! read number of bands which are included in the calculation
! of local field effects (vertex corrections)
!
    NBANDSO=-1

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSLF', NBANDSO, IERR, WRITEXMLINCAR)

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSO', NBANDSO, IERR, WRITEXMLINCAR)

    NBANDSV=NBANDSO
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSV', NBANDSV, IERR, WRITEXMLINCAR)

! use LGWLF (default is .FALSE.)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LGWLF', LGWLF, IERR, WRITEXMLINCAR)
!
! OMEGAMAX from INCAR
!
    OMEGAMAX=-30
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'OMEGAMAX', OMEGAMAX, IERR, WRITEXMLINCAR)
!
! OMEGAMIN from INCAR
!
    OMEGAMIN=-30
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'OMEGAMIN', OMEGAMIN, IERR, WRITEXMLINCAR)
!
! OMEGATL from INCAR
!
    OMEGATL=-200
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'OMEGATL', OMEGATL, IERR, WRITEXMLINCAR)
!
! type of grid for frequency integration
!
    IF (LACFDT) THEN
       ! 140 is the new default: much more reliable and fully automatic
       OMEGAGRID=140
       ! time and frequency grid for space-time algorithm
       IF (ICHIREAL >=1) THEN
          OMEGAGRID=145
       ENDIF
       ! time grid for lapalace transformed direct MP2 
       IF ( LLTDMP2 ) OMEGAGRID=141
    ELSE IF (ICHIREAL >=1) THEN
       ! real space/ imaginary time routines always work in the complex time/ frequency domain
       ! same grid as for ACFDT
       OMEGAGRID=145
       ! pL grid 145 also creates a sine-grid with an optimized setting for
       ! the sine frequencies, this is important for GWr code
       OMEGAGRID=145
    ELSE
       OMEGAGRID=0 !(default for old GW code, NOMRGA=NOMEGAR=120)
    ENDIF

    !pL if NOMEGA>0, NOMEGAR=0, then OMEGAGRID=145
    ! no real frequency points use optimized GW(cos+sin) grids
    ! this is important for GW0 with pade fit
    IF ( NOMEGAR == 0 .AND. LGW) THEN
       OMEGAGRID=145  
    ENDIF 

    IF ( LFINITE_TEMPERATURE .AND. OMEGAGRID < 160 ) OMEGAGRID = OMEGAGRID + 10 

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'OMEGAGRID', OMEGAGRID, IERR, WRITEXMLINCAR)

    ! in case ALGO is a GW algorithm 
    ! TEXT has been set by GW_FROM_ALGO if no GW was selected
    IF ( TEXT(1:18) /= 'NOT_A_GW_ALGORITHM' ) THEN
       !quick check
       IF ((ICHIREAL>=1 .AND. LGW) .OR. (NOMEGAR==0 .AND. NOMEGA>0 )) THEN
          IF (OMEGAGRID < 140 .OR. OMEGAGRID > 169 ) THEN
             CALL vtutor%error("OMEGAGRID should be set between 141 and 169 " // str(OMEGAGRID))
          ENDIF
       ENDIF
    ENDIF

! CSHIFT from INCAR
    SHIFT=-0.1_q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CSHIFT', SHIFT, IERR, WRITEXMLINCAR)

! use squeezed Coulomb kernel
    LSCK = .FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSCK', LSCK, IERR, WRITEXMLINCAR)

! use linear response functions
    LSELFENERGY=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSELFENERGY', LSELFENERGY, IERR, WRITEXMLINCAR)
! spacing of frequency depended self-energy (in screened_2e.F) 
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NOMEGA_SIGMA', NOMEGA_SIGMA, IERR, WRITEXMLINCAR)
! minimum frequency of frequency dependent self-energy calculated in screened_2e.F
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'OMEGA_SIGMA_MIN',OMEGA_SIGMA_MIN, IERR, WRITEXMLINCAR)
! maximum frequency of frequency dependent self-energy calculated in screened_2e.F
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'OMEGA_SIGMA_MAX',OMEGA_SIGMA_MAX, IERR, WRITEXMLINCAR)

! avoid the gamma point when calculating response function
    ODDONLYGW=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ODDONLYGW', ODDONLYGW, IERR, WRITEXMLINCAR)

! avoid the odd points when calculating response function
    EVENONLYGW=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EVENONLYGW', EVENONLYGW, IERR, WRITEXMLINCAR)

    IF (LGWNO) THEN
       NOMEGA=1
    ENDIF

! switch on spectral method for more than 24 grid points
    IF (NOMEGA > 24) THEN
       LSPECTRAL=.TRUE.
    ELSE
       LSPECTRAL=.FALSE.
    ENDIF
! no real frequency points, then LSPECTRAL = .FALSE.
    IF (NOMEGAR==0) THEN
       LSPECTRAL=.FALSE.
    ENDIF

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSPECTRAL', LSPECTRAL, IERR, WRITEXMLINCAR)

! LSPECTRALGW is a refined methods to calculate sigma(w) = \int W(w') G(w'-w) dw'
! it should be slightly more accurate in particular for selfenergies
    LSPECTRALGW=LSPECTRAL
    LSPECTRALGW=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSPECTRALGW', LSPECTRALGW, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)

    IF (N<1 .AND. LSPECTRALGW  .AND. LGW) THEN
       CALL vtutor%write(isAlert, LSPECTRALGWset)
    ENDIF

!
! read in ENCUTGW cutoff for Fock exchange
!
    ENCUTGW=-2
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUTGW', ENCUTGW, IERR, WRITEXMLINCAR)

!
! read in TELESCOPE tag (allows to use less k-points at high energies)
!
    TELESCOPE=0
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'TELESCOPE', TELESCOPE, IERR, WRITEXMLINCAR)
!
! read in ENCUTGWSOFT cutoff for Fock exchange
!
    ENCUTGWSOFT=ENCUTGW
    ENCUTGWSOFT_OLD=ENCUTGW
! default for ENCUTGWSOFT
! ENCUTGWSOFT now consistenly defaults to ENCUTGW*0.8
! this increases compatabilty between different calculations
! also RPA force calculations converge slower than they need to when ENCUTGWSOFT = ENCUTGW
    IF ((LACFDT .OR. LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR. &
    & LCCSD .OR. LBRACKETST) .AND. ENCUTGW /= -1 .AND. ENCUTGW /=-2 ) THEN
       ENCUTGWSOFT_OLD=ENCUTGW*0.8
    ENDIF
    IF ( ENCUTGW > 0  ) THEN
       ENCUTGWSOFT=ENCUTGW*0.8
    ENDIF
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUTGWSOFT', ENCUTGWSOFT, IERR, WRITEXMLINCAR)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUTGWSOFT', ENCUTGWSOFT_OLD, IERR, WRITEXMLINCAR)
!
! read in ENCUTLF cutoff for Fock exchange
!
    ENCUTLF=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUTLF', ENCUTLF, IERR, WRITEXMLINCAR)
!
! read in SCISSOR correction for GW or BSE calculations
! shifts the unoccupied states by a correction SCISSOR
!
    SCISSOR=0
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SCISSOR', SCISSOR, IERR, WRITEXMLINCAR)
!
! read available memory amount per MPI rank on one node 
!
    MAXMEM=2800
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MAXMEM', MAXMEM, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    ! in case not set in INCAR and low-scaling GW selected
    ! try to set this automatically in chi_super, 
    IF ( N<1 ) THEN
       LAUTOSET_MAXMEM = .TRUE. 
    ENDIF
!
! read number of electronic iterations (NELM)
!
    NELMGW=1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELMGW', NELMGW, IERR, WRITEXMLINCAR)

! read number of electronic iterations for the Hartree-Fock part (NELMHF)
!
    NELMHF=1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NELMHF', NELMHF, IERR, WRITEXMLINCAR)
!
! read dimension
!
    DIMGW=3
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'DIM', DIMGW, IERR, WRITEXMLINCAR)
    DIMGW=MIN(3,MAX(DIMGW,0))
    IF (DIMGW==0) THEN
       IDIR_MAX=1
    ELSE
       IDIR_MAX=3
    ENDIF
!
! NKREDLF reduce k-points for local field effects
!
    NKREDLFX=1
    NKREDLFY=1
    NKREDLFZ=1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDLF', NKREDLFX, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)

    IF (N>=1) THEN
       NKREDLFY=NKREDLFX
       NKREDLFZ=NKREDLFX
    ENDIF

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDLFX', NKREDLFX, IERR, WRITEXMLINCAR)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDLFY', NKREDLFY, IERR, WRITEXMLINCAR)
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDLFZ', NKREDLFZ, IERR, WRITEXMLINCAR)

    LMAXMP2=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMAXMP2', LMAXMP2, IERR, WRITEXMLINCAR)

! set default for LMAX_FOCKAE
    LMAX_FOCKAE=-1
    IF (LCHI .AND. LMAXMP2==-1) LMAX_FOCKAE=4

    HFSCREEN_ORIG=0._q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFSCREEN', HFSCREEN_ORIG, IERR, WRITEXMLINCAR)

    IF (LCHI.OR.LMP2 .OR. LMP2KPAR .OR.L2E4W .OR. LCRPA) THEN
       ! full k-point grid contains also all k-points 
       ! that are difference vectors between two other k-points
       LSHIFT_KPOINTS=.TRUE.
    ENDIF

!
! LAMBDA  from INCAR
!
    ! \fixme[vasp]: in constrmag.F Lambda defaults to 0.
    LAMBDA=1.0_q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LAMBDA', LAMBDA, IERR, WRITEXMLINCAR)
! use tricubic splines to fit electronic structure factor S(q+G) in ACFDT
! routines
    LESF_SPLINES = .FALSE. 
! anyone sane will go beyond RPA for BSE
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ESF_SPLINES', LESF_SPLINES, IERR, WRITEXMLINCAR)

    IF ( LACFDT .AND. LADDER .AND. LESF_SPLINES ) THEN
       CALL vtutor%error('tricubic spline fit for ladder diagrams in ACFDT'//&
          ' routines not implemented yet' )
    ENDIF

    LESFX_SPLINES = .FALSE. 
! exact exchange interpolation with tricubic splines
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ESFX_SPLINES', LESFX_SPLINES, IERR, WRITEXMLINCAR)

    IF ( (LESF_SPLINES .OR. LESFX_SPLINES ) .AND. LGW .AND. ICHIREAL==0  ) THEN
       CALL vtutor%error('Tricubic spline fitting (ESF_SPLINES/ESFX_SPLINES) is not supported for real-frequency GW.')

    ENDIF
!
! IEPSILON from INCAR
!
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IEPSILON', IEPSILON, IERR, WRITEXMLINCAR)
!
! damping parameter for Newton step in the RPA routines
!
    DAMP_NEWTON=0.8
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'DAMP_NEWTON', DAMP_NEWTON, IERR, WRITEXMLINCAR)

!
! mixed stochastic-deterministic compression of virtual states
!
! constant energy ratio in per cent 
    MSDGW_F = -1._q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_F'  , MSDGW_F, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    !> requested energy ratio must be a percentage between 0 and 1
    IF( N>=1 ) THEN
       IF( MSDGW_F > 1._q .OR. MSDGW_F< 0._q ) THEN
          call vtutor%error( 'requested constant energy ratio within pseudo subspaces must be a'//&
             ' real number between 0 and 1! I found a value for MSDGW_F='//str(MSDGW_F))
       ENDIF
    ENDIF
! constant energy ratio in per cent for occuppied states
    MSDGW_F_OCC = -1._q
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_F_OCC'  , MSDGW_F_OCC, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    !> requested energy ratio must be a percentage between 0 and 1
    IF( N>=1 ) THEN
       IF( MSDGW_F_OCC > 1._q .OR. MSDGW_F_OCC< 0._q ) THEN
          call vtutor%error( 'requested constant energy ratio within pseudo subspaces must be a'//&
             ' real number between 0 and 1! I found a value for MSDGW_F_OCC='//str(MSDGW_F_OCC))
       ENDIF
    ENDIF
! number of stochastic subspaces at each k-point (N_s in paper)
    MSDGW_NS=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_NS' , MSDGW_NS, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    IF( N>=1 ) THEN
       IF( MSDGW_NS < 1 ) THEN
          call vtutor%error( 'requested number of pseudo subspaces must be a'//&
             ' positive number! I found a value for MSDGW_NS='//str(MSDGW_NS))
       ENDIF
    ENDIF
! number of protected bands above the fermi energy  (N_p in paper)
    MSDGW_NP=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_NP' , MSDGW_NP, IERR, WRITEXMLINCAR)
! number of protected bands above below fermi energy 
    MSDGW_NP_OCC=-1
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_NP_OCC' , MSDGW_NP_OCC, IERR, WRITEXMLINCAR)

    IF( MSDGW_NP_OCC < 0 .AND. MSDGW_F_OCC > 0 ) THEN
       IF( IU0>=0 ) THEN
          WRITE(*,*) 'setting MSDGW_NP_OCC=0.5 *NELECT = ', INT( INFO%NELECT*0.5)
          MSDGW_NP_OCC = INT( INFO%NELECT*0.5)
       ENDIF
    ENDIF

! number of pseudo bands per stochastic subspace (N_xi in paper)
    MSDGW_NXI=2
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_NXI', MSDGW_NXI, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    IF( N>=1 ) THEN
       IF( MSDGW_NXI < 1 ) THEN
          call vtutor%error( 'requested number of pseudo bands per degenerate subspace must e a'//&
             ' positive number! I found a value for MSDGW_NXI='//str(MSDGW_NXI))
       ENDIF
    ENDIF
! number of pseudo bands per stochastic subspace (N_xi in paper)
    MSDGW_NXI_OCC=2
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_NXI_OCC', MSDGW_NXI_OCC, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    IF( N>=1 ) THEN
       IF( MSDGW_NXI_OCC < 1 ) THEN
          call vtutor%error( 'requested number of pseudo bands per degenerate subspace must e a'//&
             ' positive number! I found a value for MSDGW_NXI_OCC='//str(MSDGW_NXI_OCC))
       ENDIF
    ENDIF

    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_SEED', MSDGW_SEED, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
    IF( N > 0 ) THEN
       IF( MSDGW_SEED < 0 ) THEN
          CALL VTUTOR%ERROR( "I found a negative seed MSDGW_SEED="//str(MSDGW_SEED)//". Only positive values are allowed")
       ENDIF
    ENDIF
    
    !> if true, dumps wavecar and stops job
    MSDGW_LWAVE =.FALSE.
    CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MSDGW_LWAVE', MSDGW_LWAVE, IERR, WRITEXMLINCAR)

    CALL CLOSE_INCAR_IF_FOUND(IU5)
  END SUBROUTINE RESPONSE_READER
