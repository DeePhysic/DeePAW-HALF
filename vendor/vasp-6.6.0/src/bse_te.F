!jF:  This here is a fully reworked version! The original version was partly
!     broken and incomplete in many senses: There were inconsistencies with
!     bse.F concerning single/double precision, independent-particle spectra
!     and "gammareal" were in principle implemented (Tobias Sander) but not
!     activated, "beyond Tamm-Dancoff" which already existed in a very early
!     (although not really cleaned-up) version was fully broken/destroyed (but
!     not only here: the mess started already in bse.F setting up BSEMATRIX),
!     only diagonal elements of the dieletric tensor were calculated (although
!     calculation of all elements was already implemented -- but effectively
!     switched off by "XYZINFO=0" in the old code and switching it on had then
!     required to set for example "ZYZINFO=-1" in bse_te.F; but who else except
!     me should ever have known that in this formerly widely uncommented code?),
!     and finally also some few other "bugs"/improper changes were identified.
!     For that reason I decided to make a "big clean-up" and partial re-write,
!     and in addition I also added many comments and explanations of the code
!     (nobody could understand the undocumented old code except the authors).

#include "symbol.inc"
#include "bse.inc"

!******************************************************************************
!
!> This module provides routines to solve the BSE equation (and to obtain BSE
!> spectra) using the time-evolution method (here also referenced with label "TD"
!> like "Time Development") developed in Friedhelm Bechstedt's group in Jena (by
!> Gero Schmidt, Stephan Glutsch, Patrick Hahn, Frank Fuchs, and Claudia Roedl).
!>
!> Extenions to work with the Gamma point version contributed by Tobias Sander.
!> MAJOR CLEANUP done by jF 2020 (throw away all outdated code zombies ...)!
!> Number of comment lines and explanations was vastly increased ... (also jF).
!>
!> The time-evolution method coded here is fully described in the paper of
!>
!>   W. G. Schmidt, S. Glutsch, P.H. Hahn, and F. Bechstedt
!>   "Efficient O(N^2) method to solve the Bethe-Salpeter-equation"
!>   Physical Review B 67, 085307 (2003)
!
!******************************************************************************
MODULE bse_te


  USE prec
  USE constant
  USE wave_high
  USE bse_struct
  USE bse_base
  USE c2f_interface, ONLY : VTIME

  IMPLICIT NONE

  INTEGER,SAVE :: FIRST_ROW_INDEX=0
  INTEGER,SAVE :: NCV_LOCAL

! the flag "single_prec_bse" allows to select single precision storage of the
! BSE matrix and thus saves a factor 2 in storage but usually even a factor 2
! in CPU time since on modern CPUs usually single precision performance is
! twice the double precision performance (!) -- therefore its use is very
! highly recommeded (and still delivers high accuracy!); this flag must be set
! identical to what is set in/for bse.F -- otherwise the code will crash!!
! the default in bse.F is (or should at least be) to set this flag if one
! did not define "double_prec_bse" (or did not define anything ...) and is
! currently set to "double_prec_bse" -- if this default is changed in bse.F
! don't forget to change it here as well or otherwise the code will be broken!!
!
! Best practice is anyway to define either precompiler flag "single_prec_bse"
! or precompiler flag "double_prec_bse" explicitly in "makefile.include" ...

! for simplicity we define an precision-independent real kind "mq" and also
! define virtual BLAS and MPI routines (BDOTC, BGEMM, M_sumbse_v) -- this will
! save a lot of "#ifdef ..." blocks and code doubling later in order to always
! distinguish the two different cases "single_prec_bse" and "double_prec_bse";
! here is what we need to define for the case "single_prec_bse"

! output shall be written only on the head node ...
#define infwrite        IF (node_me==node_lead) write(*,*)
#define infwritef       IF (node_me==node_lead) write

!> restore the behaviour of the old original implementation if set to true
!> (however, this is not really what you would like to do ... !!); this is
!> at best for "educational purposes" to see whether the new implementation
!> has improved speed and accuracy with respect to the old version ...
  LOGICAL, PRIVATE :: LOLD_BSE_TE=.FALSE.
!> Use a simplified (2nd order) version of the iteration scheme (or not) which is
!> closer to the old scheme than the new corrected (effectively 3rd order) scheme
  LOGICAL, PRIVATE :: BSE_TE_ORDER2=.FALSE.
!> Use fixed time step by default
  LOGICAL, PRIVATE :: LFIXED=.TRUE.

!------------------------------------------------------------------------

! TD-BSE control flags -- will actually be set by the calling routines or
! will be defined depending on given parameters in the calling routine ...
!
!> (1) BSE, (0) IP -> according to settings of LADDER and LHARTREE (and LTRIPLET)
      INTEGER,PRIVATE  :: EXCITON
!> (1) going beyond ("ANTIRES==2"), (0) sticking with Tamm-Dancoff approximation
      INTEGER,PRIVATE  :: BEYONDTD
!> Number of spin components (2 for collinear spin-polarization and 1 else); just
!> needed for a norm factor (containing a spin-degeneracy factor ...)
      INTEGER,PRIVATE  :: ISPIN
!> Number of k-points -> what is NKPTS in VASP ...
      INTEGER,PRIVATE  :: MAXK
!> Energy cutoff used for the calculation of the spectra (the OMEGAMAX from VASP)
      REAL(mq),PRIVATE :: ENCUT_SP
!> Unit cell volume in bohr^3 (OMEGA handed over by VASP divided by a_bohr**3)
      REAL(mq),PRIVATE :: CELLVOL
!> Broadening for the spectra in eV (will be set to CSHIFT from INCAR)
      REAL(mq),PRIVATE :: GAM

!  Some further (local) variables

!> Maximum number of time steps for the TD (a kind of soft limit in this version)
      INTEGER,PRIVATE  :: NTDSTEPS
!> Number of the last actual time step in the TD
      INTEGER,SAVE,PRIVATE  :: LAST_TDSTEP
!> Matrix dimension used for the calculation (-> rank of BSE matrix from bse::CALCULATE_BSE)
      INTEGER,PRIVATE  :: MATDIM
!> Number of rows at each node (BSE matrix distributed for MPI jobs ... !)
      INTEGER,PRIVATE, ALLOCATABLE  :: MATDIM_NODE(:)
!> Corresponding "global index" of first local matrix element at each node
      INTEGER,PRIVATE, ALLOCATABLE  :: MATSTART(:)
!> Matrix dimension used for the calculation in "beyond Tamm-Dancoff" mode
!> (will be 2 times MATDIM for "beyond Tamm-Dancoff" and MATDIM otherwise)
      INTEGER,PRIVATE  :: MATDIM_BTD


! Some (local) arrays

!> Respective k-point weight for each transition (should come from bse::CALCULATE_BSE)
!> -- currently not yet really used but "we would be prepared" ...
      REAL(mq),PRIVATE, ALLOCATABLE    :: KPTWEIGHT(:)
!> Time steps and cumulated time for TD
      REAL(mq),PRIVATE, ALLOCATABLE    :: TIMESTEP_V(:,:)
!> Scalar products ("projections") for each time step (key quantity in TD !)
      COMPLEX(mq),PRIVATE, ALLOCATABLE  :: SCALAR(:,:,:)

! some local constants
      REAL(mq),PRIVATE, PARAMETER     :: HARTREE  = 2._mq*REAL(RYTOEV,mq)
      COMPLEX(mq),PRIVATE, PARAMETER  :: IMUN     = (0.0_mq,1.0_mq)
      COMPLEX(mq),PRIVATE, PARAMETER  :: INV_IMUN = (0.0_mq,-1.0_mq)
      COMPLEX(mq),PRIVATE, PARAMETER  :: CONE     = (1.0_mq,0.0_mq)
      COMPLEX(mq),PRIVATE, PARAMETER  :: CMONE    = (-1.0_mq,0.0_mq)
      COMPLEX(mq),PRIVATE, PARAMETER  :: CZERO    = (0.0_mq,0.0_mq)
#ifdef gammareal
      REAL(mq),PRIVATE,PARAMETER      :: RONE     =  1.0_mq
      REAL(mq),PRIVATE,PARAMETER      :: RMONE    =  -1.0_mq
      REAL(mq),PRIVATE,PARAMETER      :: RZERO    =  0.0_mq
#endif


      INTEGER,PRIVATE :: NODE_LEAD=0    !< Number of the leading process.
      INTEGER,PRIVATE :: NODE_ME        !< Number of my process (set by VASP).
      INTEGER,PRIVATE :: NODES          !< Total number of processes (set by VASP).
      TYPE(communic),PRIVATE :: MY_COMM !< Communicator to be used (set by VASP).
      INTEGER, PRIVATE  :: NPROW, NPCOL,MYROW,MYCOL,NQ,NP, COLS
      INTEGER, EXTERNAL    ::     NUMROC

    CONTAINS

#ifdef scaLAPACK


!**********************************************************************
!
! the remaining routines have been initially added by Claudia Rödl (as
! taken over from Patrick Hahn, André Schleife, and Frank Fuchs) but
! were then re-written by jF (using some contributions of Tobias Sander)
!
!**********************************************************************

!************* SUBROUTINE CALCULATE_BSE_TIME_EVOLUTION ****************
!
!> this is the main driver routine called from outside (bse::CALCULATE_BSE) and must
!> properly provide all necessary parameters to be used in the TD scheme
!
!**********************************************************************

      SUBROUTINE CALCULATE_BSE_TIME_EVOLUTION(WHF,SHIFT,OMEGA,ISPIN_IN,MAXK_IN,NCV,KPTWEIGHT_IN, ENCUT_SP_IN, LEXCITON, LBEYONDTD, NEDOS_IN, IO, QPOINT, BSEMATRIX, OPTMAT)

      USE base
      USE constant
      IMPLICIT NONE
      TYPE (wavespin) ::   WHF
      INTEGER :: ISPIN_IN,MAXK_IN,NCV
! mind: everything coming from "outside" has always kind type "q" ... !!
! this is exactly the reason why we copy them to other variables with a
! potentially different kind type ("mq") as defined in this module above
      REAL(q) :: SHIFT
      REAL(q) :: OMEGA
      REAL(q) :: KPTWEIGHT_IN(MAXK_IN)
      REAL(q) :: ENCUT_SP_IN
      LOGICAL :: LEXCITON
      LOGICAL :: LBEYONDTD
      INTEGER :: NEDOS_IN
      TYPE (in_struct) :: IO
      REAL(q) :: QPOINT(3)
#ifdef single_prec_bse
      GDEFS, INTENT(IN) :: BSEMATRIX(:)
      GDEFS, INTENT(IN) :: OPTMAT(:,:)
#endif

#ifdef double_prec_bse
      GDEF, INTENT(IN) :: BSEMATRIX(:)
      GDEF, INTENT(IN) :: OPTMAT(:,:)
#endif

      PROFILING_START('calculate_bse_time_evolution')
! the strip size is ceil(ncv/nodes), some nodes have less or no data
      NODE_ME = WHF%WDES%COMM%NODE_ME-1

! copy the MPI communicator for later use
      MY_COMM=WHF%WDES%COMM

! calculate excitonic spectra or independent-particle spectra only?
      IF (LEXCITON) THEN
         EXCITON=1
      ELSE
! Special note here: "EXCITON=0" employs an array "LDAGW" which has already
! been allocated in subroutine TW4O_STORE_STRIP_PREPARE in this module and is
! then set in routine INIT_BSE_OSZI_STRENGTH which can be found in file bse.F
! (which includes this module here and hence knows all definitions ... ;-)).
! aT disable special treatment of the IP case
         EXCITON=1
      ENDIF

! go beyond Tamm-Dancoff approximation (ANTIRES==2) or not (ANTIRES/=2)?
      IF (LBEYONDTD) then
         BEYONDTD=1
      ELSE
         BEYONDTD=0
      ENDIF

! if only IP spectra are calculated then the setting of ANTIRES (LBEYONDTD) is
! ignored/senseless (and therefore we silently reset BEYONDTD here to "0" ...)
! aT disable special treatment of the IP case
      !IF (.NOT.LEXCITON) BEYONDTD=0

! default for ("soft") maximum of NTDSTEPS (warning if number will be exceeded);
! usual CSHIFT ("GAM") values (>= 0.1) and not too large OMEGAMAX ("ENCUT_SP")
! values (smaller about 70-80 eV) should always result in less than 20000 steps
! and if you exceed this number then it seems that either CSHIFT is very small
! or/and OMEGAMAX was very large (and that you have either set up an excessively
! "accurate" calculation -- but since small CSHIFT only make sense for really
! dense k-meshes or giant unit cells (what is both hard to impossible to do ...)
! one should most likely check CSHIFT (one can hardly go below 0.1-0.2 eV in
! most cases). Just a very large OMEGAMAX (only possible for small systems and
! a rather small number of k-points since OMEGAMAX lets quickly explode the
! rank of the BSE matrix) could result in a bit more than 20000 steps ...;
      NTDSTEPS=20000

! maximum transition energy (OMEGAMAX in INCAR) and broadening (CSHIFT in INCAR)
      ENCUT_SP=ENCUT_SP_IN
      GAM=SHIFT

! cell volume in a.u. (Hartree units) -- internally we usually work in a.u.
! (avoids all these hbars, elementary charges, electron masses and so on ...)
      CELLVOL=OMEGA/AUTOA**3
      ISPIN=ISPIN_IN
      MAXK=MAXK_IN

! this is the global matrix dimension (rank of the BSE matrix) as set in bse.F
      MATDIM=NCV
      MATDIM_BTD=NCV
! in the case of "beyond Tamm-Dancoff" we made a simple trick: we doubled one
! dimension of "BSEMATRIX" and of the vector to be propagated where 1:NCV will
! then contain the part corresponding "AMAT" for IBSE=0 and NCV+1:2*NCV will
! contain the ("off-diagonal") part corresponding to "BMAT" for IBSE=0 ...
      IF (BEYONDTD==1) MATDIM_BTD=2*NCV

      ALLOCATE(KPTWEIGHT(MAXK))
! at the moment we do not yet support symmetry-reduced k-sets with different
! weights for each k (each transition), just non-reduced k-meshes for which
! all weights are identical (to 1/NKPTS); hence we currently set all elements
! to the same value (which is defined by KPTWEIGHT_IN(1) of the Gamma point) ...
!!!  THIS MAY NEED AN ADAPTION LATER IF VASP WOULD SUPPORT REDUCED MESHES  !!!
! (and I actually don't even know yet how and at which place it needs to enter)
      KPTWEIGHT(:)=KPTWEIGHT_IN(:)

! At this point all relevant input quantities should be known / set up;
! now we can call the "kernel routine" which actually lets set up the time
! grid, lets then perform the TD scheme and constructs (and writes out) then
! the dielectric function based on the results returned by the TD scheme ... :
      CALL CALC_AND_WRITE_SPECTRA(NEDOS_IN,LEXCITON,IO,QPOINT,BSEMATRIX,OPTMAT)

      PROFILING_STOP('calculate_bse_time_evolution')

      END SUBROUTINE CALCULATE_BSE_TIME_EVOLUTION


!******************************************************************
!> Set up the time grid for the iteration scheme. Time steps have to
!> be smaller than the inverse max. eigenvalue, otherwise the scheme
!> will diverge. These values are a good choice, with good speed and
!> good convergence, so best is to leave them that way ... .
!
! This was by the way a never-ending story: In the main driver routine
! a certain ("large") maximum number of times steps NTDSTEPS is set.
! It was always a "soft limit" insofar that exceeding this limit did
! not stop execution -- but effectively reduced TMAX what is a certain
! risk insofar that accuracy could suffer more than one might like --
! where the final issue is that by setting a finite TMAX effectively
! an extra broadening of order 1/TMAX would be added to the broadening
! set by CHSIFT in file INCAR and we know that in particular very small
! CSHIFT give rise to an excessive number of time steps (resulting in
! such a reduction of TMAX); but then the effective TMAX reached for
! step 20000 could be so small that at the end the "extra broadening"
! 1/TMAX is much larger than CSHIFT and the spectra rather look like
! spectra calculated with a very large CSHIFT (>(>) 0.2 eV); hmmm ... .
! For that reason the compromise is now following: We definitely allow
! now to exceed this limit and set the  *true*  value needed to reach
! TMAX, no matter how large it would be! But if we exceed the value of
! NTDSTEPS set before (seriously) a (serious) warning shall be issued!
!******************************************************************

      SUBROUTINE CALC_TIMESTEPS(LSET,BROADENING,IO)
      USE base
      USE vaspxml
      USE reader_tags

      IMPLICIT NONE
! We have to allocate some arrays with some dimension of number of time steps;
! initially we do not yet know this number and therefore we have first to call
! this routine   *without*   accessing/filling array TIMESTEP_V (LSET=.FALSE.)
! but just counting what we really need -- then we can make a second call that
! actually fills the arrays (which are properly allocated after the first call)
      LOGICAL   :: LSET
      REAL(mq)  :: BROADENING
      TYPE (in_struct)     IO
! local
      INTEGER   :: STEP,IERR,N
      REAL(mq)  :: TMAX,TIMESTEP,SUMTIME,TFIXED,TSCAL,DAMP,TMAXFAC
      LOGICAL   :: LOPEN,BSE_SCALE_DAMP
      CHARACTER(39) :: STRING
      CHARACTER(1) :: BSEPREC


! read some (test) parameter(s) from INCAR (for "fine-tuning" -> this is  NEW! )
      CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)
! BSEPREC selects different settings for the time step setup; only the first
! character counts (upper/lower case does not matter) and default is "N(ormal)",! optional selections are "A(ccurate)" [equivalent to "H(igh)"], "F(ast)", and
! "L(ow)". For "backward-compatibilty" also the mode "O(ld/riginal)" and (only
! partly equivalent to this) the mode "G(eorg)"/"gK" is defined. Both set up
! the original parameters for the time grids -- however "o" will then use the
! new corrected iteration scheme and also not put any restrictions while "g"
! will use the old sligthly wrong ("fixed-step") update equations (see below!)
! but also strictly limits the number of time steps to NTDSTEPS (= 20000) ... .
      STRING='Normal'
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'BSEPREC', STRING, 39, IERR, LWRITEXML=.FALSE., LCONTINUE=.TRUE.)
! something went wrong; don't stop, just restore the default and continue ...
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*) &
            'Error reading item ''BSEPREC'' from file INCAR. Using default (''Normal'') ...'
         STRING='Medium'
      ENDIF
      CALL STRIP(STRING,N,'L')
      CALL LOWER(STRING)
      BSEPREC=STRING(1:1)
! Unknown value -> mapped to 'm' like 'M(edium)' ...
      IF ((BSEPREC/='a').AND.(BSEPREC/='h').AND.(BSEPREC/='n').AND.(BSEPREC/='m').AND.(BSEPREC/='f').AND. &
          (BSEPREC/='l').AND.(BSEPREC/='o').AND.(BSEPREC/='g')) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*) &
            'Undefined value for item ''BSEPREC'' in INCAR. Using default (''Medium'') ...'
         BSEPREC='m'
      ENDIF
! see below for an explanation of these flags ...
      BSE_SCALE_DAMP=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'BSE_SCALE_DAMP', BSE_SCALE_DAMP, IERR, LWRITEXML=.FALSE., LCONTINUE=.TRUE.)
! something went wrong; don't stop, just restore the default and continue ...
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*) &
            'Error reading item ''BSE_SCALE_DAMP'' from file INCAR. Using default (''.F.'') ...'
         BSE_SCALE_DAMP=.FALSE.
      ENDIF
      BSE_TE_ORDER2=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'BSE_TE_ORDER2', BSE_TE_ORDER2, IERR, LWRITEXML=.FALSE., LCONTINUE=.TRUE.)
! something went wrong; don't stop, just restore the default and continue ...
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*) &
            'Error reading item ''BSE_TE_ORDER2'' from file INCAR. Using default (''.F.'') ...'
         BSE_TE_ORDER2=.FALSE.
      ENDIF
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LFIXED', LFIXED, IERR, LWRITEXML=.FALSE., LCONTINUE=.TRUE.)
      CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

! Now do the setting of some key parameters defining the grid; the following
! defines the 'o'riginal settings of 'g'K from the old original code of previous
! VASP versions (along with an explanation of the meaning of each quantity):
      LOLD_BSE_TE=.FALSE.
      IF ((BSEPREC=='o') .OR. (BSEPREC=='g')) THEN
         LFIXED=.FALSE.
! restore the behaviour of the old code (including errors in the update scheme;
! see below for a discussion of these -- fortunately not too large -- errors);
! for comparison this is only done for BSEPREC=g while BSEPREC=o will set the
! original mesh but will utilize then the new corrected update equations ...
         IF (BSEPREC=='g') LOLD_BSE_TE=.TRUE.
! Here some well-working settings as choosen by gK; note there was an additional
! intermediate quantitity "EVMAX" which was set to "1.5*ENCUT_SP" -- I have now
! killed "EVMAX" and did directly rescale "ENCUT_SP" (which defines "TFIXED")
! and I also have rewritten the expression for "TMAX" introducing "TMAXFAC" ...
         TMAXFAC=5._mq
! this factor appears as a kind of "damping parameter" in the scaling function
         DAMP=6._mq
! "scaling parameter" for gK's dynamic time step scaling; for usage see below
! this is the default for "LFIXED=.FALSE." and reproduces gK's original scaling
! if set to "59" -- finding a reasonable value is not so trivial: The larger
! this value the faster the increase (the less steps you need) but the danger
! to run into instabilities will also grow with this value since the final step
! may exceed the "critical stability boundary" resulting in full divergence
! what is also discussed again in more detail a bit below (so, "59" turned out
! to be a "well-balanced" choice in terms of "efficiency" and "stability" ...)
         TSCAL=59._mq
! this is the "fixed time step" (as used in the original code of Schmidt et al.)
         TFIXED=(1._mq/(ENCUT_SP*90._mq))
! The values given above are the original ones as given in the original code
! (re-written in a different way with the scaling parameters introduce above);
! but I fealt that there is still room for playing around with it ... :

! ===== HERE NOW SOME MODIFIED SETUP PARAMETERS =====
!
! First of all some general remarks and explanations concerning all parameters:
! -----------------------------------------------------------------------------
!
! The crucial point about gK's dynamic scaling is to start with a small step
! for small total times and then to increase it step by step for the long-time
! "tail" (where an "EXP(-GAM*SUMTIME)" factor appears in the time integral
! for the spectra --> reduced accuracy for long times should not harm overall
! accuracy ...). Anyway, the construction of TIMESTEP(SUMTIME) must happen in
! such a way that NEVER EVER(!) the critical time step where the scheme starts
! to diverge will be exceeded at the end (the critical time step is given by
! 1/"maximum eigenvalue of BSEMATRIX" --> something of the order "1/ENCUT_SP").
! For safety I recommend a final value of maximum about (0.85...0.95)/ENCUT_SP.
!
! Despite the exponential factor for long times empirically increasing TMAX
! increases overall accuracy; for TMAX=n/BROADENING we get a factor EXP(-n) at
! t=TMAX. Originally, gK has chosen n=5, resulting in a factor EXP(-5)=0.006738
! but larger values bring down the integrand much further and hence the residual
! error (given by the unknown integral from TMAX to infinity, where the error
! scales somehow proportional to EXP(-n)/n ...). So, let's be more generous!
! This "n" from above is what I called now "TMAXFAC" in my re-written code
! and as explained above it also determines final accuracy by adding a small
! artificial (numerical fake) extra broadening of the order 1/TMAX ... -- so
! clearly increasing TMAX should help to increase also accuracy ... !
!
! In general one can start with a much larger initial time step without really
! loosing accuracy; empirically the approximate limit for sufficiently accurate
! results is about 10 times larger than what has been chosen by gK ... . And so,
! we will increase the initial time step "TFIXED" in the following settings. Of
! course, still a certain limitation holds since despite full stability of the
! scheme for any time step smaller than the "critical time step", approaching
! this "critical time step" closer and closer will finally increase numerical
! errors very quickly and so for reasonable accuracy one should still stay
! "well below" this critical stability boundary -- only for large times (by
! the help of the exp(-GAM*t) factor) we may steadily approach the critical
! limit more and more (but not too fast! -- what is also meaningful for "DAMP"
! which actually should then not be too small ...)
!
! Finally, we have to set reasonable parameters "DAMP" and "TSCAL": Starting
! with a larger initial step but not being allowed to exceed a certain critical
! final time step at TMAX we have now to adapt the increase rate of the time
! steps accordingly, decreasing or increasing "DAMP" and "TSCAL"; the following
! is a proposal (and one might try to play around with it a bit -- as a general
! rule it must hold  TFIXED*(1+TSCAL*(1-EXP(-TMAXFAC/DAMP))) <= "0.9"/ENCUT_SP
! where the optimum is to approach the equal sign ... [and "0.9" could also be
! marginally larger, e.g. 0.92...0.95, but around "1" we definitely DIVERGE!]).
! For that reason, I also introduced then an "automatic setting" of "TSCAL"
! below, maximizing its value for a given "DAMP" (and "TMAXFAC" and "TFIXED").
!
! Finally, as a rule of thumb one should note that increasing "TMAXFAC" will
! of course increase the number of time steps. Keeping the number of time steps
! on low level requires then an increased "TFIXED" (can be increased by some
! order of 5-10 compared to gK's choice) and finally a  *decrease*  of "DAMP".
! Thereby, be aware that a large DAMP will slow down the increase of the time
! intervals (driving us more and more towards almost fixed time steps) while a
! small "DAMP" will quickly increase the time intervals from "TFIXED" towards
! (asymptotically) "TFIXED*(1+TSCAL)" -- where for the new settings the values
! of "TSCAL" (if set automatically) are much smaller than gK's value of 59
! (a typical order of magnitude will be about "10" and it decreases a bit for
! decreasing values of "DAMP" -- i.e., small "DAMP" will more quickly increase
! the time intervals but the asymptotic increase will be smaller); notice that
! for small "DAMP" we obtain again almost fixed times steps at longer times
! corresponding to time steps close to the "critical time step" -- what may
! impose serious accuracy problems if we reach too large time steps too early
! (so decreasing DAMP gives a performance boost but decreases accuracy ... !)
!
! Let us start with the most generous (and nominally most 'A'ccurate) setting;
! It is approximately as slow as the 'O'riginal or 'g'K scheme of the old code.
! However, accuracy (average and even more maximum errors) is definitly higher.
      ELSE IF (BSEPREC=='a') THEN
         TFIXED=(1._mq/(ENCUT_SP*18._mq))
         TMAXFAC=15.0_mq
         DAMP=15.0_mq
         TSCAL=-1._mq
! This is a still 'H'igh accuracy setting which is about a factor 1,4 faster
! than 'A'ccurate (or the old code) but still a bit more accurate than the old
! code (with significantly smaller maximum errors and similar average errors)
! Compared to 'A'ccurate errors are typically about a factor 2 larger ...
      ELSE IF (BSEPREC=='h') THEN
         TFIXED=1._mq/(ENCUT_SP*15._mq)
         TMAXFAC=12.0_mq
         DAMP=6.0_mq
         TSCAL=-1._mq
! This is a 'M'edium (or 'N'ormal) accuracy setting (it shall be the default).
! It is approximately a factor 2 faster than the old code or 'A'ccurate but
! accuracy is still on a good level. The maximum errors are still a bit lower
! than for the old code and the average errors are still small (although not
! anymore smaller than for the old code). So, it's not worse than the old code.
! Compared to 'A'ccurate errors are typically about a factor 4 larger ...
      ELSE IF ((BSEPREC=='m').OR.(BSEPREC=='n')) THEN
         TFIXED=(1._mq/(ENCUT_SP*12.0_mq))
         TMAXFAC=10.0_mq
         DAMP=2.5_mq
         TSCAL=-1._mq
! Finally, the 'F'ast lane with a 'L'ow(er) precision setting (which is however
! still acceptable). It is approximately a factor 3 faster than the old code
! or 'A'ccurate (or approximately a factor 1,5 faster than the default setting).
! Errors are somewhat larger than for the old code but not yet serious -- this
! setting still produces "publishable results" and marginal errors ("order of
! line thickness") become only visible (if compared to IBSE=0) for very steep
! increases or decreases of very sharp peaks as they may occur for very small
! broadenings CSHIFT -- which also produce extreme peak heights which might
! then also become slightly wrong. But there sooner or later all settings show
! increasing errors crudely scaling like 1/CSHIFT, just with smaller pre-factor.
! Compared to 'A'ccurate errors are typically about a factor 8-10 larger ...
      ELSE IF ((BSEPREC=='l').OR.(BSEPREC=='f')) THEN
         TFIXED=(1._mq/(ENCUT_SP*9._mq))
         TMAXFAC=7.5_mq
         DAMP=1.5_mq
         TSCAL=-1._mq
! and this here should never ever be reached at all ...
      ELSE
         infwritef(IO%IU0,*) 'BUG! Undefined value for BSEPREC (not correctly mapped to default value)!'
         infwritef(IO%IU0,*) '     Continuing now with the default value (medium precision) ...'
! but for safety just continue with the default setting ...
         TFIXED=(1._mq/(ENCUT_SP*12.0_mq))
         TMAXFAC=10.0_mq
         DAMP=2.5_mq
         TSCAL=-1._mq
      ENDIF

      IF (LFIXED .AND. .NOT. BSEPREC=='o' .AND. .NOT. BSEPREC=='g' ) THEN
         SELECT CASE (BSEPREC)
            CASE('a')
               TFIXED=(1._mq/(ENCUT_SP*4._mq))
               TMAXFAC=10.0_mq
            CASE('h')
               TFIXED=(1._mq/(ENCUT_SP*3._mq))
               TMAXFAC=7.5_mq
            CASE('m','n')
               TFIXED=(1._mq/(ENCUT_SP*2.5_mq))
               TMAXFAC=6.25_mq
            CASE('l','f')
               TFIXED=(1._mq/(ENCUT_SP*2._mq))
               TMAXFAC=5.0_mq
            CASE DEFAULT
               TFIXED=(1._mq/(ENCUT_SP*2.5_mq))
               TMAXFAC=6.25_mq
         END SELECT
      ENDIF

! This here is another very heuristic gimmick and could be needed or not
! (needed or not in such a sense that major deviations seem to occur in very
! "steep" parts of the many narrow and high peaks occurring for small CSHIFT
! if just looking at the bare numbers compared to results obtained with IBSE=0
! -- anyway "visual inspection" of plots still reveals a quasi perfect ("within
! line thickness") reproduction of IBSE=0 (diagonalization) reference results!
! As a consequence, I decided to  NOT  use this by default but to allow at
! least to switch this on by setting a (test) flag "BSE_SCALE_DAMP " in INCAR.
! But that's mainly for (very patient) "accuracy hunters" and for testing ...
!
! For very small broadenings (large TMAX) in particular for a small DAMP (as
! used for BSEPREC='m' or 'l') not only TMAX and hence the number of steps will
! increase (proportional to 1/CSHIFT) but in addition also an increase of
! average absolute deviations from IBSE=0 (diagonalization) reference results
! can be observed (also almost proportional 1/CHSIFT). Now remember the impact
! of a small DAMP: It quickly increases the step size from TFIXED (initial
! step size) towards the critical step size reached at TMAX -- by that the
! number of steps becomes minimal but the scheme comes closer and closer to
! a "fixed step" scheme using the "critical step size" -- but such step sizes
! (though the scheme remains stable) somehow maximize numerical errors (and
! these errors propgate and increase from step to step with a cumulated final
! error which scales then roughly proportional to the number of steps ...).
! On a lower error level (with finally similar increases of numerical errors)
! this holds by the way for any accuracy setting, i.e., also for 'h' or 'a'.
! The only possible help is then to increase DAMP to keep step sizes a bit
! longer smaller (increasing accuracy but unfortunately also the number of
! steps ... !). The choice of the "enhancement function" below is rather
! arbitrary -- the choice done shall only increase DAMP for small broadenings
! and should not make DAMP yet too much smaller for large broadenings and it
! is constructed such that DAMP takes the value from above for CHSIFT=0.2
! (and it shall not be applied to the 'o'riginal (or 'g'K) scheme ...).

! Fortunately, it only has a significant influence for smearings below 0.1 eV
! (increasing the number of steps and hence CPU times then more proportional
! to an order of (1/CSHIFT)**(3/2) instead of 1/CSHIFT -- but it improves
! accuracy and for smearings CSHIFT>0.2 it provides even a small extra speedup);
! anyway: for standard purposes I do not consider this as necessary since even
! for BSE_SCALE_DAMP=.FALSE. the formally reduced accuracy for small CSHIFT
! values still produces "publishable spectra" (plots) and is reliable enough
      IF (BSE_SCALE_DAMP.AND.((BSEPREC/='o').OR.(BSEPREC/='g'))) &
         DAMP=DAMP*(1._mq+exp(0.2_mq/BROADENING))/(1._mq+exp(1._mq))

! this is a recommended "automatic setting" for having a "just stable" scheme;
! you may set TSCAL as an independent free parameter (like in the 'o'ld scheme)
! but before you undergo a long search for each combination of TFIXED, TMAXFAC,
! and DAMP for an "optimum value" better use this (giving a minimum number of
! time steps but still safe stability over all time steps) -- the only drawback
! is that it also gives us the "steepest increase" of time steps ... . Thereby,
! one could still play a bit with the "0.93" below, e.g., also "0.95" still
! works but about here inaccuracies start to increase and at/beyond "1." there
! is then full divergence (so, I think "0.95" is the absolute maximum ...).
      IF (TSCAL < 0._mq) &
        TSCAL=(0.93_mq/(TFIXED*ENCUT_SP)-1._mq)/(1._mq-EXP(-1._mq*TMAXFAC/DAMP))
! setting parameter "TSCAL" to "0" leaves you with a fixed time step ...
      IF (LFIXED) TSCAL=0._mq

! the maximum time used should finally scale proportional to 1/BROADENING ... !
! this has also to do with the fact that cutting all down to effectively zero at
! t>TMAX gives rise to an extra broadening (being a numerical artefact) of the
! order 1/TMAX -- which should be always smaller than the "external" broadening;
! with the following setting the ratio between BROADENING (wanted broadening)
! and 1/TMAX (unintended extra broadening) would be by the way TMAXFAC and
! ideally this ratio should be close to infinity (but of course, a too large
! TMAXFAC will finally let explode the number of time steps -- so one has to
! find some compromise between accuracy [numerical artefacts] and effort ...)
      TMAX=TMAXFAC/BROADENING

      SUMTIME=0._mq
      STEP=0

! this branch sets array TIMESTEP_V (requires a call with LSET=.FALSE. before)
      IF (LSET) THEN
! for the trapezoidal rule integration later we need a minimum of two steps ;-)
         DO WHILE ((SUMTIME<=TMAX).OR.(STEP<2))
            STEP=STEP+1
! here comes gK's "dynamic rescaling" of time step "TFIXED" (re-written in some
! tricky jF-way using "TSCAL" in order to avoid an "IF (LFIXED) ...") statement:
! for TSCAL=59 (LFIXED=.FALSE.) the factor is  1+59*(1-exp(...))=60-59*exp(...)
! which is the original factor of gK and for TSCAL=0 you simply get a factor 1
!original:    TIMESTEP=TFIXED*(60._mq-59._mq*EXP(-BROADENING*SUMTIME/6._mq))
            TIMESTEP=TFIXED*(1._mq+TSCAL*(1._mq-EXP(-BROADENING*SUMTIME/DAMP)))
! store the timestep (element "1") and the integral time (element "2") ...
            TIMESTEP_V(1,STEP)=TIMESTEP
            TIMESTEP_V(2,STEP)=SUMTIME
! ... and since above SUMTIME shall be the sum of all    previous    timesteps
! we just add up the current timestep  HERE  (after storing the current SUMTIME)
! --> be aware that by incrementing SUMTIME here it is guaranteed that element
! TIMESTEP_V(2,1)=0., i.e., we correctly start at t=0 as needed later ... !
            SUMTIME=SUMTIME+TIMESTEP
! here we fully restore the original 'g(K)' behaviour: never allow more than
! NTDSTEPS time steps! -- exit then and issue a warning / "error" message:
! (the "SUMTIME<=TMAX" just checks whether we would continue at all or have
! accidently reached the final number of steps ... -- then stay silent ...)
            IF ((BSEPREC=='g').AND.(STEP==NTDSTEPS).AND.(SUMTIME<=TMAX)) THEN
               IF (IO%IU0>=0) THEN
! Too many time steps ... -- never ever publish such a run! If this has happened
! with the old code then you should definitely re-check your results and might
! need to publish errata (I hope everybody took the printed warnings serious!).
! The basic problem is that implicitly TMAX is cut down to SUMTIME(step 20000)
! and as already noted several times a too small TMAX would introduce a (very)
! large artificial extra broadening (which becomes indeed visible in plots ...)
! and if the effective broadening becomes even much larger than CSHIFT then the
! spectra will not be what you initially intended to calculate (maybe not yet
! complete nonsense but fully "smeared out" -- so, better don't publish this!)
                  WRITE(IO%IU0,'(1x,A)') 'Urgent    W A R N I N G    from routine CALC_TIMESTEPS:'
                  WRITE(IO%IU0,'(1x,A,I8,A)') 'You have selected BSEPREC=''g'' and will exceed the maximum of ',NTDSTEPS,' steps!'
                  WRITE(IO%IU0,'(1x,A)') 'Be warned that for compatibility the behavior of the old code will be mimicked!'
                  WRITE(IO%IU0,'(1x,A)') 'This means that implicitly we will then reset TMAX to the actual integral time'
                  WRITE(IO%IU0,'(1x,A,F10.4,A)') 'reached so far which is currently ',SUMTIME,' while the original value of TMAX'
                  WRITE(IO%IU0,'(1x,A,F10.4,A)') 'was set to ',TMAX,' (if the difference is too large you will be in trouble!).'
                  WRITE(IO%IU0,*)
                  WRITE(IO%IU0,'(1x,A)') '!! RESULTING ERRORS MIGHT BE DEVASTATING!! Do not publish this result ... !!'
                  WRITE(IO%IU0,*)
               ENDIF
               IF (IO%IU6>=0) THEN
! ... and this should also go to OUTCAR!
                  WRITE(IO%IU6,'(1x,A)') 'Urgent    W A R N I N G    from routine CALC_TIMESTEPS:'
                  WRITE(IO%IU6,'(1x,A,I8,A)') 'You have selected BSEPREC=''g'' and will exceed the maximum of ',NTDSTEPS,' steps!'
                  WRITE(IO%IU6,'(1x,A)') 'Be warned that for compatibility the behavior of the old code will be mimicked!'
                  WRITE(IO%IU6,'(1x,A)') 'This means that implicitly we will then reset TMAX to the actual integral time'
                  WRITE(IO%IU6,'(1x,A,F10.4,A)') 'reached so far which is currently ',SUMTIME,' while the original value of TMAX'
                  WRITE(IO%IU6,'(1x,A,F10.4,A)') 'was set to ',TMAX,' (if the difference is too large you will be in trouble!).'
                  WRITE(IO%IU6,*)
                  WRITE(IO%IU6,'(1x,A)') '!! RESULTING ERRORS MIGHT BE DEVASTATING!! Do not publish this result ... !!'
                  WRITE(IO%IU6,*)
               ENDIF
               EXIT
            ENDIF
         ENDDO !while

! Internal consistency check -- should never fail if nobody ruins this code ...
         IF (LAST_TDSTEP/=STEP) THEN
            IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Internal ERROR in routine CALC_TIMESTEPS. Number of time steps not reproducible.'
            STOP
         ENDIF

! Information output and warnings only issued upon second call with LSET=.TRUE.
         IF (LAST_TDSTEP.GT.NTDSTEPS) THEN
! if NTDSTEPS was set large enough one should not arrive here; if it happens
! then issue a warning and advice to check input parameters but do not stop
            IF (IO%IU0>=0) THEN
               WRITE(IO%IU0,'(1x,A)') 'Urgent    W A R N I N G    from routine CALC_TIMESTEPS:'
! Too many time steps ...
               IF (LAST_TDSTEP.LE.2*NTDSTEPS) WRITE(IO%IU0,'(1x,A)') 'You will need a suspiciously large number of time steps for the time evolution!'
! ... or a dramatically excessive number of time steps
               IF (LAST_TDSTEP.GT.2*NTDSTEPS) WRITE(IO%IU0,'(1x,A)') 'You will need an excessively large number of time steps for the time evolution!'
               WRITE(IO%IU0,'(1x,A,I8,A,I8,A)') 'Usually this should not take more than ',NTDSTEPS,' steps but you need ',LAST_TDSTEP,' !'
               WRITE(IO%IU0,'(1x,A)') 'I bet something is wrong -- in most cases a too small CSHIFT may cause this.'
               WRITE(IO%IU0,'(1x,A)') 'Re-check your input parameters please. The job will continue without failures'
               WRITE(IO%IU0,'(1x,A)') 'but be advised that this job might take an excessively long execution time ...'
               WRITE(IO%IU0,*)
            ENDIF
         ENDIF

         IF (IO%IU0>=0) THEN
            WRITE(IO%IU0,'(1x,(A),I8)') 'Number of time steps: ',LAST_TDSTEP
            WRITE(IO%IU0,'(1x,(A),E13.6,(A))') 'First time step: ',TIMESTEP_V(1,1),' 1/eV'
            WRITE(IO%IU0,'(1x,(A),E13.6,(A))') 'Last  time step: ',TIMESTEP_V(1,LAST_TDSTEP),' 1/eV'
            WRITE(IO%IU0,*)
         ENDIF
         IF (IO%IU6>=0) THEN
            WRITE(IO%IU6,'(1x,(A),I8)') 'Number of time steps: ',LAST_TDSTEP
            WRITE(IO%IU6,*)
            WRITE(IO%IU6,'(1x,(A),E13.6,(A))') 'First time step: ',TIMESTEP_V(1,1),' 1/eV'
            WRITE(IO%IU6,'(1x,(A),E13.6,(A))') 'Last  time step: ',TIMESTEP_V(1,LAST_TDSTEP),' 1/eV'
            WRITE(IO%IU6,*)
         ENDIF
      ELSE ! lset
! information output (shall only appear upon first call)
         IF (IO%IU0>=0) THEN
            WRITE(IO%IU0,'(1x,A,F10.4)') ' Omega_max:',ENCUT_SP
            WRITE(IO%IU0,'(1x,A,F10.4)') 'Broadening:',BROADENING
            WRITE(IO%IU0,'(1x,A,F10.4)') '     T_max',TMAX
            IF (BSEPREC=='o') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Old setup (but corrected update scheme and >20000 steps allowed)'
            IF (BSEPREC=='g') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Old setup (full backward-compatibility mode, NOT RECOMMENDED!)'
            IF (BSEPREC=='a') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Accurate'
            IF (BSEPREC=='h') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  High'
            IF (BSEPREC=='m') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Medium (default)'
            IF (BSEPREC=='n') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Medium (default)'
            IF (BSEPREC=='l') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Low (Fast)'
            IF (BSEPREC=='f') WRITE(IO%IU0,'(1x,A)') '   BSEPREC:  Low (Fast)'
            IF (.NOT.LOLD_BSE_TE) THEN
               IF (BSE_TE_ORDER2) WRITE(IO%IU0,'(1x,A)') '    Scheme:  2nd order update equations'
               IF (.NOT.BSE_TE_ORDER2) WRITE(IO%IU0,'(1x,A)') '    Scheme:  3rd order update equations'
               IF (BSE_SCALE_DAMP.AND.((BSEPREC/='o').OR.(BSEPREC/='g'))) &
                  WRITE(IO%IU0,'(1x,A)') '             (with experimental broadening-dependent precision enhancement)'
            ELSE
               WRITE(IO%IU0,'(1x,A)') '    Scheme:  Old slightly incorrect update scheme (quality like BSEPREC=Low)'
            ENDIF
            WRITE(IO%IU0,*)
         ENDIF
         IF (IO%IU6>=0) THEN
            WRITE(IO%IU6,'(1x,A,F10.4)') ' Omega_max:',ENCUT_SP
            WRITE(IO%IU6,'(1x,A,F10.4)') 'Broadening:',BROADENING
            WRITE(IO%IU6,'(1x,A,F10.4)') '     T_max:',TMAX
            IF (BSEPREC=='o') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Old setup (but corrected update scheme and >20000 steps allowed)'
            IF (BSEPREC=='g') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Old setup (full backward-compatibility mode, NOT RECOMMENDED!)'
            IF (BSEPREC=='a') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Accurate'
            IF (BSEPREC=='h') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  High'
            IF (BSEPREC=='m') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Medium (default)'
            IF (BSEPREC=='n') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Medium (default)'
            IF (BSEPREC=='l') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Low (Fast)'
            IF (BSEPREC=='f') WRITE(IO%IU6,'(1x,A)') '   BSEPREC:  Low (Fast)'
            IF (.NOT.LOLD_BSE_TE) THEN
               IF (BSE_TE_ORDER2) WRITE(IO%IU6,'(1x,A)') '    Scheme:  2nd order update equations'
               IF (.NOT.BSE_TE_ORDER2) WRITE(IO%IU6,'(1x,A)') '    Scheme:  3rd order update equations'
               IF (BSE_SCALE_DAMP.AND.((BSEPREC/='o').OR.(BSEPREC/='g'))) &
                   WRITE(IO%IU6,'(1x,A)') '             (with experimental broadening-dependent precision enhancement)'
            ELSE
               WRITE(IO%IU6,'(1x,A)') '    Scheme:  Old slightly incorrect update scheme (quality like BSEPREC=Low)'
            ENDIF
            WRITE(IO%IU6,*)
         ENDIF

! this just counts how many steps we will need (same as above without storing):
         DO WHILE ((SUMTIME<=TMAX).OR.(STEP<2))
            STEP=STEP+1
            TIMESTEP=TFIXED*(1._mq+TSCAL*(1._mq-EXP(-BROADENING*SUMTIME/DAMP)))
            SUMTIME=SUMTIME+TIMESTEP
! test output I only want to see (or not) for the first call with LSET=.FALSE.
            IF (MOD(STEP,1000)==0) THEN
               IF (IO%IU0>=0) THEN
! in case you really want to see some more information uncomment this ...
!                  WRITE(IO%IU0,'(I6,2(1x,E13.6))') STEP,TIMESTEP,SUMTIME
               ENDIF
            ENDIF
            IF ((BSEPREC=='g').AND.(STEP==NTDSTEPS).AND.(SUMTIME<=TMAX)) EXIT
         ENDDO !while
! save the final value of the counter (number of time steps) and return
         LAST_TDSTEP=STEP
      ENDIF !lset

      END SUBROUTINE CALC_TIMESTEPS


!========================================================================
!>     "Kernel" routine performing the time development algorithm for BSE:
!>  
!>     We solve the complex differential equation     i dv(t)/dt = H v(t).
!>     Important note: all routines here assume "Hartree units" and hence
!>     there will be no "hbar" and all other factors need to be converted!
!>     Due to the leap-frog scheme used here we need two initial vectors.
!>     The second vector is calculated using a double-(half-)time step.
!========================================================================

      SUBROUTINE TIME_DEVELOP_EXCI(FULLDIM,IO,BSEMATRIX,OPTMAT)

GPACC USE moffload

      USE base

      implicit none
      TYPE (in_struct)        :: IO
      INTEGER                 :: FULLDIM
#ifdef single_prec_bse
    GDEFS, INTENT(IN) :: BSEMATRIX(:)
    GDEFS, INTENT(IN) :: OPTMAT(:,:)
#endif

#ifdef double_prec_bse
    GDEF, INTENT(IN) :: BSEMATRIX(:)
    GDEF, INTENT(IN) :: OPTMAT(:,:)
#endif
! local
      INTEGER                 :: IDIR,KDIR,LDIR,I,J,K,NBTD,IBTD,IV(3)
      REAL(q)                 :: VPUT0,CPUT0,VPUT1,CPUT1,VPUT,CPUT
      REAL(mq)                :: POM,POM2
! vectors to be iterated: every step needs to calculate a third vector out of
! two vectors of two previous time steps and all that has to be done for three
! different directions (polarizations); all these vectors are stored in a single
! array "VECTOR(1:FULLDIM,1:3,1:3)" where the last index counts the different
! iteration steps, the mid index the three directions, and the first dimension
! essentially holds the vector data itself (FULLDIM complex numbers) ...
      COMPLEX(mq),ALLOCATABLE :: VECTOR(:,:,:)
! some further helper arrays for intermediate results or inputs
      COMPLEX(mq),ALLOCATABLE :: VECTOR_0(:,:)         ! initial vectors (usually optical matrix elements)
      COMPLEX(mq),ALLOCATABLE :: VECTORELEMENTS(:,:,:) ! holds the local blocks of v' = bsematrix v
      COMPLEX(mq),ALLOCATABLE :: VECTOR_D(:,:)         ! merged result vector v' = bsematrix v
      COMPLEX(mq),ALLOCATABLE :: VECTOR_S(:,:,:)       ! "super vector" for BEYONDTD==1
#ifdef gammareal
! needed in addition as helper vectors for the "gammareal" case
      REAL(mq),ALLOCATABLE    :: VECTOR0_REAL(:,:)    ! real counterpart of "VECTOR_0"
      REAL(mq),ALLOCATABLE    :: VECTORELEMENTS_REAL(:,:,:)  ! real counterpart of "VECTORELEMENTS"
      REAL(mq),ALLOCATABLE    :: VECTOR_REAL(:,:)     ! real counterpart of "VECTOR_D" 
      REAL(mq),ALLOCATABLE    :: VECTORS_REAL(:,:,:)  ! real counterpart of "VECTOR_S"
#endif
      INTEGER                 :: IERR
! BLAS functions (as re-defined before, like "mq", depending on precision)
      COMPLEX(mq)             :: BDOTC
      INTEGER                 :: IFAIL
      INTEGER                 :: MATDIM, MATSTRT
      INTEGER                 :: NODE, NPROW, NPCOL, MYROW, MYCOL
      INTEGER, EXTERNAL       :: NUMROC

      IFAIL = 0

      PROFILING_START('time_develop_exci')
      CALL BLACS_GRIDINFO(BSE_DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)

      ! global index of the first local element
      MATSTRT=1
      DO I=1,NODE_ME
         MATSTRT=MATSTRT + NUMROC(BSE_DESC(M_),BSE_DESC(MB_),I-1,0,NPROW)
      ENDDO

      ! local number of rows
      MATDIM = NUMROC(BSE_DESC(M_),BSE_DESC(MB_),NODE_ME,0,NPROW)

! some (local) dimension parameter for handling "beyond Tamm-Dancoff" (or not)
      NBTD=1
      IF (BEYONDTD==1) NBTD=2

! allocations of local arrays; "FULLDIM" should be "MATDIM_BTD" from outside
! and all the "3"s are the three cartesian coordinates (x, y, z) -- where in
! the Gamma-only case it is also sometimes "6" (3x real + 3x imaginary part!)
! and as explained above for VECTOR there is a third coordinate counting the
! three different iteration steps involved in the update equation; some of
! the helper arrays also carry a third dimension "NBTD" depending on "BEYONDTD"
      ALLOCATE(VECTOR(FULLDIM,3,3))
      ALLOCATE(VECTOR_0(FULLDIM,3))
      ALLOCATE(VECTOR_D(FULLDIM,3))
      ALLOCATE(VECTORELEMENTS(MATDIM,3,NBTD))
      !$ACC ENTER DATA CREATE(VECTOR,VECTOR_0,VECTOR_D,VECTORELEMENTS) __IF_ASYNC__
#ifdef gammareal
      ALLOCATE(VECTOR0_REAL(FULLDIM,3))
      ALLOCATE(VECTOR_REAL(FULLDIM,6))
      ALLOCATE(VECTORELEMENTS_REAL(MATDIM,6,NBTD))
      !$ACC ENTER DATA CREATE(VECTOR_REAL,VECTOR0_REAL,VECTORELEMENTS_REAL) __IF_ASYNC__
#endif
      IF (BEYONDTD==1) THEN
         ALLOCATE(VECTOR_S(FULLDIM,3,NBTD))
         !$ACC ENTER DATA CREATE(VECTOR_S) __IF_ASYNC__
#ifdef gammareal
         ALLOCATE(VECTORS_REAL(FULLDIM,6,NBTD))
         !$ACC ENTER DATA CREATE(VECTORS_REAL) __IF_ASYNC__
#endif
      ELSE
         ALLOCATE(VECTOR_S(1,1,1))
         !$ACC ENTER DATA CREATE(VECTOR_S) __IF_ASYNC__
#ifdef gammareal
         ALLOCATE(VECTORS_REAL(1,1,1))
         !$ACC ENTER DATA CREATE(VECTORS_REAL) __IF_ASYNC__
#endif
      ENDIF

! "SCALAR" (global module array, comes from outside ...) finally will hold all
! projections of current vectors v(t) onto the initial vectors v(0); initialize
! since at the end we have to sum up partial projections of "CPU-local stripes"!
      SCALAR=CZERO
! "IV(1:3)" is now a kind of tricky "address array" for the three vectors that
! need to be iterated later (for transparency denoted as "VECTOR_1", "VECTOR_2",
! and "VECTOR_3" belonging to iterations i-1, i, and i+1 in all later comments
! -- in the code you have then to identify "VECTOR_i <=> VECTOR(:,:,IV(i)) ...);
! at the end from iteration to iteration we have to "shift the index" in such
! a way that having determined VECTOR_3 = f(VECTOR_1,VECTOR_2) before entering
! the next step we can "discard" the original (oldest) VECTOR_1 and VECTOR_2
! has then to take over the role of VECTOR_1 and VECTOR_3 the role of VECTOR_2
! (i.e., in pseudo-code we needed to perform copy operations VECTOR_1=VECTOR_2
! and VECTOR_2=VECTOR_3 -- and could then overwrite again VECTOR_3 in the next
! iteration); but instead of a physical copy of one vector to another vector
! we will basically only (cyclically) shift the "addresses" of the vectors so
! that at a later stage the following initial setting will change every step!
      IV(1)=1 ; IV(2)=2 ; IV(3)=3

!=============================================================================
! starting point for leap frog scheme (t=0) according to Eq. (20) of the paper
!=============================================================================
#ifdef ACC_OFFLOAD
      ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
      !$ACC ENTER DATA COPYIN(OPTMAT) __IF_ASYNC__

      IF (BEYONDTD==0) THEN
#ifdef gammareal
! OPTMAT is REAL for Gamma-only (!)
         !$ACC KERNELS PRESENT(VECTOR_0,OPTMAT,VECTOR0_REAL) __IF_ASYNC__
         VECTOR0_REAL(1:FULLDIM,1:3)=OPTMAT(1:FULLDIM,1:3)
         VECTOR_0(1:FULLDIM,1:3)=CMPLX(VECTOR0_REAL(1:FULLDIM,1:3),0._mq,mq)
         !$ACC END KERNELS
#else
! but in general OPTMAT is a complex quantity
         !$ACC KERNELS PRESENT(VECTOR_0,OPTMAT) __IF_ASYNC__
         VECTOR_0(1:FULLDIM,1:3)=OPTMAT(1:FULLDIM,1:3)
         !$ACC END KERNELS
#endif
      ELSEIF (BEYONDTD==1) THEN
! furth: I was digging around a lot in many old code(s) and found all possible
!        combinations (CONJG or not, plus or minus, ...) but that's definitely
!        what we need (and delivers exactly the same result as IBSE=0 ... ;-))
#ifdef gammareal
! OPTMAT is REAL for Gamma-only (!)
         !$ACC KERNELS PRESENT(VECTOR_0,OPTMAT,VECTOR0_REAL) __IF_ASYNC__
         VECTOR0_REAL(1:FULLDIM/2,1:3)=OPTMAT(1:FULLDIM/2,1:3)
! for the "coupling" part we should just start with minus the "resonant" part
         VECTOR0_REAL(FULLDIM/2+1:FULLDIM,1:3)=RMONE*OPTMAT(1:FULLDIM/2,1:3)
         VECTOR_0(1:FULLDIM,1:3)=CMPLX(VECTOR0_REAL(1:FULLDIM,1:3),0._mq,mq)
         !$ACC END KERNELS
#else
! but in general OPTMAT is a complex quantity
         !$ACC KERNELS PRESENT(VECTOR_0,OPTMAT) __IF_ASYNC__
         VECTOR_0(1:FULLDIM/2,1:3)=OPTMAT(1:FULLDIM/2,1:3)
! for the "coupling" part we should just start with minus the "resonant" part
         VECTOR_0(FULLDIM/2+1:FULLDIM,1:3)=CMONE*(OPTMAT(1:FULLDIM/2,1:3))
         !$ACC END KERNELS
#endif
      ENDIF

! "VECTOR_0" is now left untouched forever and is the initial vector we need to
! project onto ("initial value"); for time propagation the general propagation
! step is explained below in detail and as you will see three vectors called
! VECTOR_1, VECTOR_2, and VECTOR_3 are introduced representing the vectors at
! timesteps "i-1", "i", and "i+1" -- at this point here (i=0) we have now simply
! to set VECTOR_1 equal to VECTOR_0 and then initialization step one is complete
#ifdef gammareal
      !$ACC KERNELS PRESENT(VECTOR_REAL,VECTOR0_REAL,OPTMAT) __IF_ASYNC__
      VECTOR_REAL(1:FULLDIM,1:3)=VECTOR0_REAL(1:FULLDIM,1:3)
      VECTOR_REAL(1:FULLDIM,4:6)=RZERO
      !$ACC END KERNELS
#else
      !$ACC KERNELS PRESENT(VECTOR,VECTOR_0) __IF_ASYNC__
      VECTOR(:,:,IV(1))=VECTOR_0(:,:)
      !$ACC END KERNELS
#endif

! start some timing
      IF (NODE_ME==NODE_LEAD) CALL VTIME(VPUT0,CPUT0)

!==============================================================================
! Now we determine the second initial vector required for the leap-frog scheme:
! In Eq. 21 we needed two vectors ("i-1" and "i") to determine the third one
! ("i+1"). At the current point (this is now actually the determination of the
! vector at t=dt ("i=1")!) we needed now an "i=0" (this is "VECTOR_1=VECTOR_0)
! but also some "i=-1" -- which is not known/defined at all. So like all central
! difference methods determining a new vector out of two others at t-dt and
! t one has now to replace Eq. 21 by some "forward-approximation" for the very
! first step where only one vector is present (known). What is used here is
! not really described in the paper of Schmidt et al. but it corresponds to
!
!    v(1) = v(0) + dt/(i*hbar) H v(0) + 0.5*[dt/(i*hbar)]^2 H H v(0)
!
! One can use different philosophies to arrive at this formula: Basically one
!   i) can use a second-order expansion of the exact propagator exp(-i*dt*H)
!      to approximate the formal exact solution v(1) = exp(-i*dt*H).v(0)
!  ii) can start from  v(1) = v(0) - i*H.v(dt/2)  as a kind of "midpoint"
!      scheme for the first "forward step" and then we can
!      a) either insert v(dt/2)=v(0)-i*(dt/2)*H.v(0) [and apply H to it] or
!      b) start from v(dt)=v(0)-i*dt*H.v(0), calculate then H.v(dt), and
!         use the approximation/interpolation H.v(dt/2)=0.5*[H.v(0)+H.v(1)]
! No matter what you do, you always end up with the same formula (given above).
!==============================================================================

! In the following we always have to calculate products of the BSE matrix times
! some vectors -- and this for three directions (x, y,and z). Mind that the BSE
! matrix is locally stored in "stripes" (DIMENSION(stripe,full)) and that the
! result is then also just a stripe of the resulting vector. This result will
! be stored in "VECTORELEMENTS"; therefore, the final result vector (with full
! dimension) has then to be constructed (on all nodes) out of these stripes
! by gathering all the stripes from the different nodes and to store them at
! the appropriate places in the global vector_<out>(full) ... (see below).
!
! In principle, "matrix times vector" is a classical *GEMV (level 2 BLAS) task.
! Hence, in the "historical code" of Schmidt et al. all was done with "*GEMV".
! However, mind that we apply all to three vectors, so we can also build up a
! "matrix" vector(dim,3) (or vectorelements(stripe,3)) and formally turn all
! into a matrix-matrix multiplication, i.e., a typical *GEMM task (level 3
! BLAS) -- which is formally much more efficient than *GEMV. However, note
! that the very small matrix dimension "3" for the "matrix of vectors" makes
! *GEMM still rather inefficient compared to very large matrices (the full
! power of *GEMM is unleashed only for large dimensions ...). So, one may now
! ask whether such "*GEMM" calls really can outperform three "*GEMV" calls. The
! answer is: The gain is limited but there is a small gain (unless there is a
! desastrous problem with the implementation of the "*GEMM" routines -- but then
! you should rather exchange your BLAS library and not turn back to "*GEMV").
! There are by the way situations (Gamma-only or "beyond Tamm-Dancoff") where
! the extra matrix dimension is even larger (6 or even 12 instead of 3) and in
! these cases there is definitely a substantial visible performance boost ... .
! FOR THAT REASON ANY OLD "*GEMV" BRANCH WAS ELIMINATED FROM THIS VERSION ... !

! The basic structure of the following appears three times (this "first step",
! then for the "second step", and inside the loop for all the remaining steps;
! just mind that inside the loop we have just to process VECTOR_2 then ...):
      IF (EXCITON==1) THEN
         IF (BEYONDTD==0) THEN
#ifdef gammareal
! case gammareal (BGEMM=DGEMM/SGEMM); here the additional dimension is even 6(!)
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,6,FULLDIM,RONE, &
                 BSEMATRIX,MATDIM,VECTOR_REAL,FULLDIM, &
                 RZERO,VECTORELEMENTS_REAL,MATDIM)
#else
! general complex case (BGEMM=CGEMM/ZGEMM)
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,3,FULLDIM,CONE, &
                 BSEMATRIX,MATDIM,VECTOR(1,1,IV(1)),FULLDIM, &
                 CZERO,VECTORELEMENTS,MATDIM)
#endif

         ELSEIF (BEYONDTD==1) THEN
! "beyond Tamm-Dancoff" is now the "tricky" part: formally we have to perform
! four matrix multiplications (two parts of the BSEMATRIX multiplied with two
! parts of the vectors) -- all with dimensions MATDIM=FULLDIM/2. But some trick
! (building up some "supervector" containing the corresponding parts of our
! vectors arranged in an appropriate way) allows to rewrite this then into one
! *single*   matrix multiplication again: Basically the "supervector" (VECTOR_S)
! contains each part of the vectors twice (!) and the storage order must be
!            "part 1      part 2      -part2      -part1"
! This "supervector" with this peculiar order has to be multiplied now with
! matrix "BSEMATRIX" (where basically for the second part one had formally to
! exchange the storage-order of the two parts of BSEMATRIX -- but it is
! equivalent to exchangeing the two parts of the vectors ... !). Again I
! was digging in various old code (and some papers, also Claudia Roedl's
! thesis) and came to the conclusion that this must be what we need:
!
! Note: Despite the fact that formally this will now take four times more
! operations the "small dimension" of 6 of the "matrix of (super)vectors"
! gives already such a significant performance boost over dimension 3 that
! effectively (for a highly tuned BLAS like MKL or OpenBLAS) the CPU time
! approximately increases by a factor about 2 ... max. 2.4 only (!) -- like
! for IBSE=0 where two diagonalizations instead of one need to be performed
! (i.e., no real disadvantage / "slow-down" will be observed for IBSE=1 in
! competition with IBSE=0 -- leaving the "cross-over" matrix rank beyond which
! IBSE=1 will be faster than IBSE=0 approximately the same as for the case
! of the Tamm-Dancoff approximation ...)
#ifdef gammareal
! build a "supervector" of different parts of VECTOR_1: first VECTOR_1 ...
            !$ACC KERNELS PRESENT(VECTOR_REAL,VECTORS_REAL) __IF_ASYNC__
            VECTORS_REAL(1:FULLDIM,1:6,1)=VECTOR_REAL(1:FULLDIM,1:6)
! ... and then as "second element" the two parts of VECTOR_1 swapped times -1
            VECTORS_REAL(1:FULLDIM/2,1:6,2)=RMONE*VECTOR_REAL(FULLDIM/2+1:FULLDIM,1:6)
            VECTORS_REAL(FULLDIM/2+1:FULLDIM,1:6,2)=RMONE*VECTOR_REAL(1:FULLDIM/2,1:6)
            !$ACC END KERNELS
! now the big "all-in-one" multiplication (here even with dimension "12" now!)
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,12,FULLDIM,RONE, &
                 BSEMATRIX,MATDIM,VECTORS_REAL,FULLDIM, &
                 RZERO,VECTORELEMENTS_REAL,MATDIM)
#else
! build a "supervector" of different parts of VECTOR_1: first VECTOR_1 ...
            !$ACC KERNELS PRESENT(VECTOR_S,VECTOR) __IF_ASYNC__
            VECTOR_S(1:FULLDIM,1:3,1)=VECTOR(1:FULLDIM,1:3,IV(1))
! ... and then as "second element" the two parts of VECTOR_1 swapped times -1
            VECTOR_S(1:FULLDIM/2,1:3,2)=CMONE*VECTOR(FULLDIM/2+1:FULLDIM,1:3,IV(1))
            VECTOR_S(FULLDIM/2+1:FULLDIM,1:3,2)=CMONE*VECTOR(1:FULLDIM/2,1:3,IV(1))
            !$ACC END KERNELS
! now the big "all-in-one" multiplication (here we still get a dimension "6"!)
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,6,FULLDIM,CONE, &
                 BSEMATRIX,MATDIM,VECTOR_S,FULLDIM, &
                 CZERO,VECTORELEMENTS,MATDIM)
#endif
         ENDIF
      ENDIF     !exciton

!===========================================================================
! rebuild VECTORELEMENTS & VECTOR_1 in complex-type for the case "gammareal"
!===========================================================================
#ifdef gammareal
      !$ACC KERNELS PRESENT(VECTORELEMENTS_REAL,VECTORELEMENTS) __IF_ASYNC__
      DO IBTD=1,NBTD
       DO IDIR=1,3
         DO I=1,MATDIM
            VECTORELEMENTS(I,IDIR,IBTD)=CMPLX(VECTORELEMENTS_REAL(I,IDIR,IBTD),VECTORELEMENTS_REAL(I,IDIR+3,IBTD),mq)
         ENDDO
       ENDDO
      ENDDO
      !$ACC END KERNELS

      !$ACC KERNELS PRESENT(VECTOR,VECTOR_REAL) __IF_ASYNC__
      DO IDIR=1,3
         DO I=1,FULLDIM
            VECTOR(I,IDIR,IV(1))=CMPLX(VECTOR_REAL(I,IDIR),VECTOR_REAL(I,IDIR+3),mq)
         ENDDO
      ENDDO
      !$ACC END KERNELS
#endif

!!========================================================================
!! merge VECTORELEMENTS into some temporary work array VECTOR_D = "H v(0)"
!!========================================================================


       !$ACC KERNELS PRESENT(VECTOR_D) __IF_ASYNC__
       VECTOR_D=MAGIC
       !$ACC END KERNELS

#ifndef CUBLASMP
       !$ACC UPDATE SELF(VECTORELEMENTS)  __IF_ASYNC__
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC  CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif
       CALL PGGEMR2Dtebse(FULLDIM/NBTD,3,VECTORELEMENTS(:,:,1), 1, 1, DVEC_DESC, VECTOR_D, 1, 1, VEC_DESC, DVEC_DESC(2))
       IF (BEYONDTD==1 .AND. EXCITON==1) THEN
         CALL PGGEMR2Dtebse(FULLDIM/NBTD,3,VECTORELEMENTS(:,:,2), 1, 1, DVEC_DESC, VECTOR_D, 1+FULLDIM/NBTD, 1, VEC_DESC, DVEC_DESC(2))
       ENDIF
#ifndef CUBLASMP
GPACC  CALL OFFLOAD_ON_POP
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       !$ACC UPDATE DEVICE(VECTOR_D)  __IF_ASYNC__
#endif

       CALL M_bcastbsete_g_from(MY_COMM, VECTOR_D, SIZE(VECTOR_D), 1)


       !$ACC KERNELS PRESENT(VECTOR_D) __IF_ASYNC__
       IFAIL = COUNT(VECTOR_D == MAGIC)
       !$ACC END KERNELS

       IF (IFAIL>0) THEN
          CALL vtutor%bug("values in the polarization vector are not set correctly", __FILE__, __LINE__)
       ENDIF
! Finally build an "intermediate" result VECTOR_2 (here v(0)-i*dt*H*v(0)/hbar;
! mind: we work in atomic units with "hbar=1" and so you won't find any hbar!)
      !$ACC KERNELS PRESENT(VECTOR_D,VECTOR) __IF_ASYNC__
      VECTOR(:,:,IV(2))=VECTOR(:,:,IV(1))+(INV_IMUN*TIMESTEP_V(1,1))*VECTOR_D(:,:)
      !$ACC END KERNELS
      IF (NODE_ME==NODE_LEAD) THEN
         CALL VTIME(VPUT,CPUT)
         infwritef(*,'(1x,(a),f13.6,(a))') 'Time for vector 1: ',CPUT-CPUT0,' sec'
         CALL VTIME(VPUT1,CPUT1)
      ENDIF

!=======================================================================
! now copy VECTOR_D to VECTOR_1 (to be used as input for the next step)
!=======================================================================
      !$ACC KERNELS PRESENT(VECTOR,VECTOR_D) __IF_ASYNC__
      VECTOR(:,:,IV(1))=VECTOR_D(:,:)
      !$ACC END KERNELS


!===============================================================================
! apply the BSE matrix ("H") once again to VECTOR_1 (containing here "H*v(0)",
! i.e., the result will then be H*(H*v(0)) -- exactly what we need next ...)
! and it is just the same code (no comments) as above ("copied from above" ...)
!===============================================================================
      IF (EXCITON==1) THEN
#ifdef gammareal
         !$ACC PARALLEL LOOP PRESENT(VECTOR,VECTOR_REAL) __IF_ASYNC__
         DO IDIR=1,3
               VECTOR_REAL(:,IDIR)=REAL(VECTOR(:,IDIR,IV(1)),mq)
               VECTOR_REAL(:,IDIR+3)=AIMAG(VECTOR(:,IDIR,IV(1)))
         ENDDO

         IF (BEYONDTD==0) THEN
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,6,FULLDIM,RONE, &
                 BSEMATRIX,MATDIM,VECTOR_REAL,FULLDIM, &
                 RZERO,VECTORELEMENTS_REAL,MATDIM)
         ELSEIF (BEYONDTD==1) THEN
            !$ACC KERNELS PRESENT(VECTORS_REAL,VECTOR_REAL) __IF_ASYNC__
            VECTORS_REAL(1:FULLDIM,1:6,1)=VECTOR_REAL(1:FULLDIM,1:6)
            VECTORS_REAL(1:FULLDIM/2,1:6,2)=RMONE*VECTOR_REAL(FULLDIM/2+1:FULLDIM,1:6)
            VECTORS_REAL(FULLDIM/2+1:FULLDIM,1:6,2)=RMONE*VECTOR_REAL(1:FULLDIM/2,1:6)
            !$ACC END KERNELS
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,12,FULLDIM,RONE, &
                 BSEMATRIX,MATDIM,VECTORS_REAL,FULLDIM, &
                 RZERO,VECTORELEMENTS_REAL,MATDIM)
         ENDIF
#else
         IF (BEYONDTD==0) THEN
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,3,FULLDIM,CONE, &
                 BSEMATRIX,MATDIM,VECTOR(1,1,IV(1)),FULLDIM, &
                 CZERO,VECTORELEMENTS,MATDIM)
         ELSEIF (BEYONDTD==1) THEN
            !$ACC KERNELS PRESENT(VECTOR_S,VECTOR) __IF_ASYNC__
            VECTOR_S(1:FULLDIM,1:3,1)=VECTOR(1:FULLDIM,1:3,IV(1))
            VECTOR_S(1:FULLDIM/2,1:3,2)=CMONE*VECTOR(FULLDIM/2+1:FULLDIM,1:3,IV(1))
            VECTOR_S(FULLDIM/2+1:FULLDIM,1:3,2)=CMONE*VECTOR(1:FULLDIM/2,1:3,IV(1))
            !$ACC END KERNELS
            IF(MATDIM>0) &
            CALL BGEMM('N','N', MATDIM,6,FULLDIM,CONE, &
                 BSEMATRIX,MATDIM,VECTOR_S,FULLDIM, &
                 CZERO,VECTORELEMENTS,MATDIM)
         ENDIF

#endif
      ENDIF

#ifdef gammareal
! after this we do not need to use anymore VECTORELEMENTS_REAL for EXCITON=0
      IF (EXCITON==1) THEN
       !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(VECTOR_S,VECTOR) __IF_ASYNC__
       DO IBTD=1,NBTD
        DO IDIR=1,3
         DO I=1,MATDIM
            VECTORELEMENTS(I,IDIR,IBTD)=CMPLX(VECTORELEMENTS_REAL(I,IDIR,IBTD),VECTORELEMENTS_REAL(I,IDIR+3,IBTD),mq)
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif

       !$ACC KERNELS PRESENT(VECTOR_D) __IF_ASYNC__
       VECTOR_D=MAGIC
       !$ACC END KERNELS

#ifndef CUBLASMP
       !$ACC UPDATE SELF(VECTORELEMENTS)  __IF_ASYNC__
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC  CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif
       CALL PGGEMR2Dtebse(FULLDIM/NBTD,3,VECTORELEMENTS(:,:,1), 1, 1, DVEC_DESC, VECTOR_D, 1, 1, VEC_DESC, DVEC_DESC(2))
       IF (BEYONDTD==1 .AND. EXCITON==1) THEN
         CALL PGGEMR2Dtebse(FULLDIM/NBTD,3,VECTORELEMENTS(:,:,2), 1, 1, DVEC_DESC, VECTOR_D, 1+FULLDIM/NBTD, 1, VEC_DESC, DVEC_DESC(2))
       ENDIF
#ifndef CUBLASMP
GPACC  CALL OFFLOAD_ON_POP
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       !$ACC UPDATE DEVICE(VECTOR_D)  __IF_ASYNC__
#endif

       CALL M_bcastbsete_g_from(MY_COMM, VECTOR_D, SIZE(VECTOR_D), 1)

       !$ACC KERNELS PRESENT(VECTOR_D) __IF_ASYNC__
       IFAIL = COUNT(VECTOR_D == MAGIC)
       !$ACC END KERNELS

       IF (IFAIL>0) THEN
          CALL vtutor%bug("values in the polarization vector are not set correctly", __FILE__, __LINE__)
       ENDIF

!===============================================================================
! Now we can build the final result VECTOR_2 (reflecting "v(1)" at this point!)
! Then all "initialization" is finished and we can finally enter the "big loop"
!===============================================================================
      !$ACC KERNELS PRESENT(VECTOR) __IF_ASYNC__
      VECTOR(:,:,IV(2))=VECTOR(:,:,IV(2))+(0.5_mq*(INV_IMUN*TIMESTEP_V(1,1))**2)*VECTOR_D(:,:)
      !$ACC END KERNELS

! before entering the main iteration loop (mind at the moment we are at "i=1")
! we have now to restore again VECTOR_1 to be VECTOR_0 which is v(0) and then
! with VECTOR_2 being v(1) at this point the first iteration of the following
! time-step loop will properly determine VECTOR_3 being v(2) out of it ... .
#ifdef gammareal
      !$ACC KERNELS PRESENT(VECTOR,VECTOR0_REAL) __IF_ASYNC__
      VECTOR(1:FULLDIM,1:3,IV(1))=CMPLX(VECTOR0_REAL(1:FULLDIM,1:3),0._mq,mq)
      !$ACC END KERNELS
#else
      !$ACC KERNELS PRESENT(VECTOR_0,VECTOR) __IF_ASYNC__
      VECTOR(:,:,IV(1))=VECTOR_0(:,:)
      !$ACC END KERNELS
#endif

      IF (NODE_ME==NODE_LEAD) THEN
         CALL VTIME(VPUT,CPUT)
         infwritef(*,'(1x,(a),f13.6,(a))') 'Time for vector 2: ',CPUT-CPUT1,' sec'
! somehow I got the feeling the second step alone gives a more precise forecast
         infwritef(*,'(1x,(a),f13.2,(a))') 'The full time propagation will take approx. ',(CPUT-CPUT1)*REAL(LAST_TDSTEP,mq),' sec'
         infwrite
      ENDIF

!==============================================================================
! The calculation of all timesteps follows. According to Schmidt et al.
! the Schroedinger-like Eq. 19 is solved using the method of central
! differences (Eq.21):     i*hbar (v(i+1)-v(i-1))/(2*dt) = H v(i)
! This yields (Eq.21 resolved):   v(i+1) = v(i-1) + 2*dt/(i*hbar) H v(i)
! And at the end we need also   scalar(i) = <v(0)|v(i)>   entering Eq. 18
!
! "H" corresponds of course to the BSE matrix (the locally stored stripe)
! Just for confusion: as it is implemented here "VECTOR_2" takes the role
! of "v(i)" which will be multiplied with matrix "H" (times 2*dt/(i*hbar))
! and "VECTOR_1" represents "v(i-1)", and "VECTOR_3" is then "v(i+1)", while
! "VECTOR_0" remains untouched and is "v(0)" onto which we have to project.
! However, also mind that the current loop counter below will effectively
! correspond to "i+1" (which is the time index of VECTOR_3) so that our
! loop counting and the time index count "i" in the comments differ by one!
!
! There is now just one nasty detail left: Schmidt et al. applied UNIFORM
! (means equally-spaced) grids for everything but as you may remember some
! comments in the routine setting up the time grid, gK had now introduced
! some NON-UNIFORM grid ("dynamical increase of time step for large times t")
! and immediately one may raise doubts that the equations given in the paper
! of Schmidt et al. can be naively used 1:1 also for non-uniform grids.
!
! So, some important notes about "central differences" for first derivatives:
! For  non-uniform  grids one has to be a bit cautious; normally the central
! differences formula is derived by taking a WEIGHTED average of a forward
! difference representing Psi'(t+dt_p/2) and a backward difference which
! represents Psi'(t-dt_m/2) with weights that reflect the position of "t"
! within the interval [t-dt_m/2,t+dt_p/2] which is not(!) a central position
! if dt_p /= dt_m but some "asymmetric" position t /= [dt_m+dt_p]/2 ... :
!
!     Psi'(t) =    dt_m/(dt_p+dt_m) * (Psi(t+dt_p) - Psi(t))     /dt_p)
!                + dt_p/(dt_p+dt_m) * (Psi(t)      - Psi(t-dt_m))/dt_m)
!
! Obviously (and logically) the weighting factors must be dt_m/(dt_p+dt_m)
! and dt_p/(dt_p+dt_m) -- and correctly sum up to "one" as it should be.
!
! With dt_p = dt_m = dt this reduces to Psi'(t) = (Psi(t+dt)-Psi(t-dt)/(2*dt)
! but if dt_p and dt_p are different we need to rewrite it a bit (with a common
! nominator (dt_m+dt_p) as generalization of the "2*dt") which gives us finally
!
! Psi'(t)=   {  (dt_m/dt_p)* [Psi(t+dt_p) - Psi(t)    ]
!             + (dt_p/dt_m)* [Psi(t) -     Psi(t-dt_m)] } / (dt_m+dt_p)
!
! Now applied to the numerical approximation of "i*Psi'(t)=H.Psi(t)" with a
! representation of Psi'(t) according to what is written above the update rule
! also needs to be generalized into a much more complicated update equation:
!
!   Psi(t+dt_p) = [1-(dt_p/dt_m)^2] * Psi(t) + (dt_p/dt_m)^2 * Psi(t-dt_m)
!                 -i*(dt_p+dt_m)*(dt_p/dt_m) * H.Psi(t)
!
! For dt_p=dt_m=dt we recover again the simple (well-known) update equation
!
!    Psi(t+dt) = Psi(t-dt) - i*H.Psi(t)*(2*dt)
!
!
! Note: From a mathematical point of view (and there are also strict derivations
!       in mathematics literature, see for example M.K. Bowen and Ronald Smith,
!       Proceedings of The Royal Society A Mathematical Physical and Engineering
!       Sciences 461, 1975-1997 (2005) Appendix A, Eq. (A 3b) ...) this is an
!       effectively third order approximation for the first derivative based on
!       three function values; as an alternative a second order approximation
!       (in a mathematician's view) would by the way yield the simple result
!
!            Psi(t+dt) = Psi(t-dt) - i*H.Psi(t)*(dt_p+dt_m)
!
!       which is almost the old formula for uniform grids but with "2*dt"
!       (which actually can be identified as "2*dt_p" in the original code)
!       replaced by its generalization "dt_p+dt_m" -- so the uniform mesh
!       code is not really comletely wrong, at the end it contains just a
!       wrong "dt" (with dt_p>dt_m in our case it is always a bit too large);
!       there are by the way somehow disputes (depends also on the structure
!       of the mathematical problem) whether the third order approximation
!       is always better than the simple second order approximation -- for
!       me empirically the second order approximation does not make a visible
!       difference (frequently no difference at all ...) and so the dispute
!       may now go on but currently the default will be the third order scheme
!       (both schemes require the same quantities and hence identical effort).
!       Anyway, the user may have the choice: I introduced an INCAR (test) flag
!       "BSE_TE_ORDER2" which allows to switch on this 2nd order approximation.
!
!==============================================================================

! loop over all (remaining) time steps (next step to be "2" now ...)
      DO I=2,LAST_TDSTEP
! From now on we have to apply the BSE matrix ("H") to VECTOR_2 (being "v(I-1)",
! the actual vector of the previous time step) -- code is again same as above!
         PROFILING_START('te_main_loop_out')
         IF (EXCITON==1) THEN
#ifdef gammareal
            !$ACC PARALLEL LOOP PRESENT(VECTOR,VECTOR_REAL) __IF_ASYNC__
            DO IDIR=1,3
               VECTOR_REAL(:,IDIR)=REAL(VECTOR(:,IDIR,IV(2)),mq)
               VECTOR_REAL(:,IDIR+3)=AIMAG(VECTOR(:,IDIR,IV(2)))
            ENDDO

            IF (BEYONDTD==0) THEN
               IF(MATDIM>0) &
               CALL BGEMM('N','N', MATDIM,6,FULLDIM,RONE, &
                    BSEMATRIX,MATDIM,VECTOR_REAL,FULLDIM, &
                    RZERO,VECTORELEMENTS_REAL,MATDIM)
            ELSEIF (BEYONDTD==1) THEN
               !$ACC KERNELS PRESENT(VECTORS_REAL,VECTOR_REAL)  __IF_ASYNC__
               VECTORS_REAL(1:FULLDIM,1:6,1)=VECTOR_REAL(1:FULLDIM,1:6)
               VECTORS_REAL(1:FULLDIM/2,1:6,2)=RMONE*VECTOR_REAL(FULLDIM/2+1:FULLDIM,1:6)
               VECTORS_REAL(FULLDIM/2+1:FULLDIM,1:6,2)=RMONE*VECTOR_REAL(1:FULLDIM/2,1:6)
               !$ACC END KERNELS
               IF(MATDIM>0) &
               CALL BGEMM('N','N', MATDIM,12,FULLDIM,RONE, &
                    BSEMATRIX,MATDIM,VECTORS_REAL,FULLDIM, &
                    RZERO,VECTORELEMENTS_REAL,MATDIM)
            ENDIF
#else
            IF (BEYONDTD==0) THEN
               PROFILING_START('te_main_loop')
               IF(MATDIM>0) &
               CALL BGEMM('N','N', MATDIM,3,FULLDIM,CONE, &
                    BSEMATRIX,MATDIM,VECTOR(1,1,IV(2)),FULLDIM, &
                    CZERO,VECTORELEMENTS,MATDIM)
               PROFILING_STOP('te_main_loop')
            ELSEIF (BEYONDTD==1) THEN
               !$ACC KERNELS PRESENT(VECTOR_S,VECTOR)  __IF_ASYNC__
               VECTOR_S(1:FULLDIM,1:3,1)=VECTOR(1:FULLDIM,1:3,IV(2))
               VECTOR_S(1:FULLDIM/2,1:3,2)=CMONE*VECTOR(FULLDIM/2+1:FULLDIM,1:3,IV(2))
               VECTOR_S(FULLDIM/2+1:FULLDIM,1:3,2)=CMONE*VECTOR(1:FULLDIM/2,1:3,IV(2))
               !$ACC END KERNELS
               IF(MATDIM>0) &
               CALL BGEMM('N','N', MATDIM,6,FULLDIM,CONE, &
                    BSEMATRIX,MATDIM,VECTOR_S,FULLDIM, &
                    CZERO,VECTORELEMENTS,MATDIM)
            ENDIF
#endif
         ENDIF

#ifdef gammareal
         IF (EXCITON==1) THEN
          !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(VECTORELEMENTS,VECTORELEMENTS_REAL) __IF_ASYNC__
          DO IBTD=1,NBTD
           DO IDIR=1,3
            DO K=1,MATDIM
               VECTORELEMENTS(K,IDIR,IBTD)=CMPLX(VECTORELEMENTS_REAL(K,IDIR,IBTD),VECTORELEMENTS_REAL(K,IDIR+3,IBTD),mq)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
#endif
!========================================================================
! merge VECTORELEMENTS into some temporary work array VECTOR_D = "H v(0)"
!========================================================================
       !$ACC KERNELS PRESENT(VECTOR_D) __IF_ASYNC__
       VECTOR_D=MAGIC
       !$ACC END KERNELS

#ifndef CUBLASMP
       !$ACC UPDATE SELF(VECTORELEMENTS)  __IF_ASYNC__
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
GPACC  CALL OFFLOAD_ON_PUSH(.FALSE.)
#endif
       CALL PGGEMR2Dtebse(FULLDIM/NBTD,3,VECTORELEMENTS(:,:,1), 1, 1, DVEC_DESC, VECTOR_D, 1, 1, VEC_DESC, DVEC_DESC(2))
       IF (BEYONDTD==1 .AND. EXCITON==1) THEN
         CALL PGGEMR2Dtebse(FULLDIM/NBTD,3,VECTORELEMENTS(:,:,2), 1, 1, DVEC_DESC, VECTOR_D, 1+FULLDIM/NBTD, 1, VEC_DESC, DVEC_DESC(2))
       ENDIF
#ifndef CUBLASMP
GPACC  CALL OFFLOAD_ON_POP
       !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       !$ACC UPDATE DEVICE(VECTOR_D)  __IF_ASYNC__
#endif

       CALL M_bcastbsete_g_from(MY_COMM, VECTOR_D, SIZE(VECTOR_D), 1)

       !$ACC KERNELS PRESENT(VECTOR_D) __IF_ASYNC__
       IFAIL = COUNT(VECTOR_D == MAGIC)
       !$ACC END KERNELS

       IF (IFAIL>0) THEN
          CALL vtutor%bug("values in the polarization vector are not set correctly", __FILE__, __LINE__)
       ENDIF


!==============================================================================
! Now do the "projections" for the current time step (to be stored in "SCALAR");
! we have to project VECTOR_1 ("v(I-2)", the vector in the  second-last step)
! what might make you wondering a bit but mind that  no  projection has been
! done so far and that for "I=2" (first projection!) VECTOR_1 still contains
! "v(0)"; so we are somehow two indices "behind". But that is what we actually
! want since then we really start at "(t=)0"! And since the array "SCALAR" is
! starting with index "1" and we want to fill element "1" for loop counter "I=2"
! we actually have to store everything into array element "SCALAR(I-1,...)" ...
! At the end "SCALAR" should be filled with all these projections which appear
! in Eq. 18 as "<mu^j|xi^j(t)(t)>" -- just with the generalization that now not
! only diagonal elements are calculated but the full matrix <mu^i|xi^j(t)(t)>
! where "I-1" is the time (counter) ["t"], "i" is "KDIR", "j" is "LDIR and the
! "mu" represents "VECTOR_0" and "xi(t)" our current "VECTOR_1" for each "t"
! (not yet confused? -- congratulations, then you got the "iteration count")
!==============================================================================

! here we project "stripes" only -- which are then all summed   *after*   the
! main loop (and the final projection) for further use in calculating epsilon
         IF (BEYONDTD==0) THEN
            DO LDIR=1,3
             DO KDIR=1,3
! this shortcut here is explained at the very end ...
               IF (KDIR==MODULO(LDIR,3)+1) CYCLE
               IF(MATDIM>0) &
                  SCALAR(I-1,KDIR,LDIR)=BDOTC(MATDIM,VECTOR_0(MATSTRT,KDIR),1,VECTOR(MATSTRT,LDIR,IV(1)),1)
             ENDDO
            ENDDO
         ELSE IF (BEYONDTD==1) THEN
!jF: again a search for the correct combinations of +/- etc. led me to this ...
            DO LDIR=1,3
             DO KDIR=1,3
               IF (KDIR==MODULO(LDIR,3)+1) CYCLE
               IF(MATDIM>0) &
                  SCALAR(I-1,KDIR,LDIR)=BDOTC(MATDIM,VECTOR_0(MATSTRT,KDIR),1,VECTOR(MATSTRT,LDIR,IV(1)),1) &
                          -BDOTC(MATDIM,VECTOR_0(MATSTRT+FULLDIM/2,KDIR),1,VECTOR(MATSTRT+FULLDIM/2,LDIR,IV(1)),1)
             ENDDO
            ENDDO
         ENDIF
        !$ACC UPDATE DEVICE(SCALAR)  __IF_ASYNC__



!====================================================================
! now propagate the vector in time and store the result in VECTOR_3 !
!====================================================================

! uniform grid/paper:  VECTOR_3=VECTOR_1+2._mq*INV_IMUN*TIMESTEP_V(1,I)*VECTOR_D
! corrected TD scheme by jF:  for  NON-UNIFORM  time grids this here is correct
! (identifying TIMESTEP_V(1,I) with "dt_p" and TIMESTEP_V(1,I-1) with "dt_m" in
! the lengthy update expression derived above for non-uniform time grids; as you
! may recognize basically there always appears a factor "dt_p/dt_m" ("POM") or
! its square ("POM2") and finally a sum of both timesteps (instead of "2*dt")!
         IF (LOLD_BSE_TE) THEN
! the "wrong scheme" (only used for BSEPREC='g') ...
            !$ACC KERNELS PRESENT(VECTOR,VECTOR_D)  __IF_ASYNC__
            VECTOR(:,:,IV(3))=VECTOR(:,:,IV(1))+2._mq*INV_IMUN*TIMESTEP_V(1,I)*VECTOR_D(:,:)
            !$ACC END KERNELS
         ELSE
! the new improved (3rd order) scheme (used in all other cases, default)
            POM=TIMESTEP_V(1,I)/TIMESTEP_V(1,I-1)
! setting this one can by the way go down to the second order approximation of
! the first derivative (set test flag BSE_TE_ORDER2=.TRUE. in INCAR to use it)
            IF (BSE_TE_ORDER2) POM=1._mq
            POM2=POM*POM
            !$ACC KERNELS PRESENT(VECTOR,VECTOR_D)  __IF_ASYNC__
            VECTOR(:,:,IV(3))=(1._mq-POM2)*VECTOR(:,:,IV(2))+POM2*VECTOR(:,:,IV(1))+INV_IMUN*POM*(TIMESTEP_V(1,I-1)+TIMESTEP_V(1,I))*VECTOR_D(:,:)
            !$ACC END KERNELS
         ENDIF

! and now prepare everything for the following iteration and advance the index
! of the vector count by VECTOR_1 now taking over the role of VECTOR_2 and by
! VECTOR_2 now taking over the role of VECTOR_3 ...; this means symbolically
!         VECTOR_1=VECTOR_2
!         VECTOR_2=VECTOR_3
! but actually (see comment on "IV" at the beginning) we do not copy physically
! but will just cyclically shift the roles (counting up all indices IV(i) by one
! and taking all modulo 3), so that 1 points now to 2, 2 to 3, and 3 to 1 ...;
! confused? -- never mind, it's exactly the way it works ...
         IV(1)=MODULO(IV(1),3)+1 ; IV(2)=MODULO(IV(2),3)+1 ; IV(3)=MODULO(IV(3),3)+1

         IF (NODE_ME==NODE_LEAD) THEN
            IF (MOD(I,500)==0) THEN
! some "progress monitoring" (but each 500 steps only => usually 20-40 lines)
               !CALL CPU_TIME(Z2)
               CALL VTIME(VPUT,CPUT)
               infwritef(*,'(1x,(a),i6,(a),f13.2,(a))') 'Step ',i,' done: ',CPUT-CPUT1,' sec'
            ENDIF
         ENDIF
         PROFILING_STOP('te_main_loop_out')
      ENDDO       ! time_step(i)



!==============================================================================
! after leaving the time step loop we still have one more vector that can be
! projected: the counter is now at "LAST_TDSTEP+1" and still the preparation
! step for the next iteration (which is not executed anymore) was done and so
! VECTOR_1 contains now former VECTOR2 (which now contains former VECTOR_3).
! So finally we project once again VECTOR_1 onto VECTOR_0 obtaining the element
! "SCALAR(LAST_TDSTEP,...)" which is still missing ... :
!==============================================================================

      I=LAST_TDSTEP+1    ! "just for consistency with the loop counter above"
! it is the same code as in the loop above ...
      IF (BEYONDTD==0) THEN
         DO LDIR=1,3
          DO KDIR=1,3
            IF (KDIR==MODULO(LDIR,3)+1) CYCLE
            IF(MATDIM>0) &
               SCALAR(I-1,KDIR,LDIR)=BDOTC(MATDIM,VECTOR_0(MATSTRT,KDIR),1,VECTOR(MATSTRT,LDIR,IV(1)),1)
          ENDDO
         ENDDO
      ELSE IF (BEYONDTD==1) THEN
         DO LDIR=1,3
          DO KDIR=1,3
            IF (KDIR==MODULO(LDIR,3)+1) CYCLE
            IF(MATDIM>0) &
               SCALAR(I-1,KDIR,LDIR)=BDOTC(MATDIM,VECTOR_0(MATSTRT,KDIR),1,VECTOR(MATSTRT,LDIR,IV(1)),1) &
                       -BDOTC(MATDIM,VECTOR_0(MATSTRT+FULLDIM/2,KDIR),1,VECTOR(MATSTRT+FULLDIM/2,LDIR,IV(1)),1)
          ENDDO
         ENDDO
      ENDIF
      !$ACC UPDATE DEVICE(SCALAR)  __IF_ASYNC__
      !$ACC EXIT DATA COPYOUT(OPTMAT) __IF_ASYNC__


! and now collect (sum up) all data for SCALAR (sum all stripes ...)
     CALLMPI(M_sumbsete_v(MY_COMM,SCALAR,9*LAST_TDSTEP))




     IF (NODE_ME==NODE_LEAD) THEN
! and the final timing ...
         CALL VTIME(VPUT,CPUT)
         infwritef(*,'(1x,(a),f13.2,(a))') 'All time steps done. Time was ',(CPUT - CPUT1)/1._mq,' sec'
         infwrite
         IF (IO%IU6>=0) WRITE(IO%IU6,'(1x,A,F13.2,A)') 'CPU time (time steps)',(CPUT - CPUT1)/1._mq,'sec'
         IF (IO%IU6>=0) WRITE(IO%IU6,*)
      ENDIF
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
      !$ACC EXIT DATA DELETE(VECTORELEMENTS,VECTOR,VECTOR_0,VECTOR_D) __IF_ASYNC__
      !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

! but now we are really finished -- clean up and ciao ...
      DEALLOCATE(VECTORELEMENTS)
      DEALLOCATE(VECTOR_D)
      DEALLOCATE(VECTOR)
      DEALLOCATE(VECTOR_0,VECTOR_S)
#ifdef gammareal
      DEALLOCATE(VECTOR_REAL)
      DEALLOCATE(VECTORELEMENTS_REAL)
      DEALLOCATE(VECTOR0_REAL,VECTORS_REAL)
#endif

      PROFILING_STOP('time_develop_exci')

      END SUBROUTINE TIME_DEVELOP_EXCI


!     *******************************************************************
!>     "Driver routine" to calculate and write out the spectra ...
!     *******************************************************************

      SUBROUTINE CALC_AND_WRITE_SPECTRA(NEDOS_IN,LEXCITON,IO, QPOINT,BSEMATRIX,OPTMAT)
      USE base
      USE chi_glb, ONLY: LGWLF
#ifdef VASP_HDF5
      USE vhdf5
#endif

      IMPLICIT NONE
      INTEGER              :: NEDOS_IN
      TYPE (in_struct)     :: IO
#ifdef single_prec_bse
      GDEFS, INTENT(IN) :: BSEMATRIX(:)
      GDEFS, INTENT(IN) :: OPTMAT(:,:)
#endif

#ifdef double_prec_bse
      GDEF, INTENT(IN) :: BSEMATRIX(:)
      GDEF, INTENT(IN) :: OPTMAT(:,:)
#endif
! local
      INTEGER              :: DIR,IDIR,I,J
      INTEGER              :: MAXOMEGA,KDIR,LDIR,NEDOS
      REAL(mq)             :: NORM,W,DELTAW,DT
      COMPLEX(mq)          :: EXPFAC,ALPHA
      COMPLEX(mq), POINTER  :: EPSILON_TOT(:,:,:)
      COMPLEX(q), ALLOCATABLE :: EPSILON(:,:,:)
      REAL(mq)             :: WEIGHT, TIMESTEP
! needed for writing vasprun.xml (must be always of precision "q" !!!)
      REAL(q),POINTER      :: EPSILON_RE(:,:,:), EPSILON_IM(:,:,:)
      LOGICAL              :: LEXCITON
      REAL(q)              :: QPOINT(3)
#ifdef VASP_HDF5
      CHARACTER(LEN=:), ALLOCATABLE :: PREFIX
#endif


      PROFILING_START('calc_and_write_spectra')

      IF (IO%IU0>=0) THEN
         WRITE(IO%IU0,'(1x,(A))') 'Calculation of spectra started ...'
         WRITE(IO%IU0,*)
         IF (.NOT. LEXCITON) WRITE(IO%IU0,'(1x,A)')'Calculating IP spectra only'
         IF (LEXCITON) THEN
            WRITE(IO%IU0,'(1x,A)',ADVANCE='NO')'Calculating BSE spectra using'
            IF (BEYONDTD==0) WRITE(IO%IU0,'(1x,A)')'the Tamm-Dancoff approximation'
            IF (BEYONDTD==1) WRITE(IO%IU0,'(1x,A)')'full res.-antires. coupling'
         ENDIF
         WRITE(IO%IU0,'(1x,A,I10)')'Matrix rank: ',MATDIM
         WRITE(IO%IU0,*)
      ENDIF
      IF (IO%IU6>=0) THEN
         WRITE(IO%IU6,*)
         WRITE(IO%IU6,'(1x,A)')'Time evolution parameters used:'
         WRITE(IO%IU6,'(1x,A)')'-------------------------------'
         WRITE(IO%IU6,*)
         IF (.NOT. LEXCITON) WRITE(IO%IU6,'(1x,A)')'Calculating IP spectra only'
         IF (LEXCITON) THEN
            WRITE(IO%IU6,'(1x,A)',ADVANCE='NO')'Calculating BSE spectra using'
            IF (BEYONDTD==0) WRITE(IO%IU6,'(1x,A)')'the Tamm-Dancoff approximation'
            IF (BEYONDTD==1) WRITE(IO%IU6,'(1x,A)')'full res.-antires. coupling'
         ENDIF
         WRITE(IO%IU6,'(1x,A,I10)')'Matrix rank: ',MATDIM
         WRITE(IO%IU6,*)
      ENDIF

! fake allocation only in order to avoid an uninitialized pointer
      ALLOCATE(TIMESTEP_V(2,1))
! call setup routine for time step parameters: first to count what we need ...
      CALL CALC_TIMESTEPS(.FALSE.,GAM,IO)
! now we can do the (re-)allocations with correct dimension "LAST_TDSTEP" ...
      DEALLOCATE(TIMESTEP_V)
      ALLOCATE(TIMESTEP_V(2,LAST_TDSTEP))
      ALLOCATE(SCALAR(LAST_TDSTEP,3,3))
      !$ACC ENTER DATA CREATE(SCALAR)  __IF_ASYNC__
! ... and now we can call it a second time actually setting up TIMESTEP_V ...
      CALL CALC_TIMESTEPS(.TRUE.,GAM,IO)

! default frequency step; it should be ENCUT_SP/(NEDOS_IN-1) but minimum 0.01
      DELTAW=MIN(0.01_mq,(ENCUT_SP/MAX((NEDOS_IN-1),1)))
! set NEDOS
      NEDOS=MAX(NEDOS_IN,1+NINT(ENCUT_SP/DELTAW))
! number of frequencies shall be determined by DELTAW set above; we add few eV
! to ENCUT_SP (in the original code it was 4 eV = 400 additional frequencies)
      MAXOMEGA=MAX(NINT(MIN(ENCUT_SP+4,1.1_mq*ENCUT_SP)/DELTAW),NEDOS-1)
      IF (IO%IU0>=0) THEN
         WRITE(IO%IU0,'(1x,A,1x,F8.5)') 'dOmega  :',DELTAW
         WRITE(IO%IU0,'(1x,A,1x,I6)') 'N_Omega :',NEDOS
         WRITE(IO%IU0,*)
      ENDIF
      IF (IO%IU6>=0) THEN
         WRITE(IO%IU6,'(1x,A,1x,F8.5)') 'dOmega  :',DELTAW
         WRITE(IO%IU6,'(1x,A,1x,I6)') 'N_Omega :',NEDOS
         WRITE(IO%IU6,*)
      ENDIF

! now we can allocate the array(s) storing the dielectric matrix
      ALLOCATE(EPSILON_TOT(0:MAXOMEGA,3,3))
      ALLOCATE(EPSILON_RE(1:NEDOS,3,3))
      ALLOCATE(EPSILON_IM(1:NEDOS,3,3))
      ALLOCATE(EPSILON(1:NEDOS,3,3))

! norm factor appearing in front of the integral in Eq. 18 (in Hartree units!!)
      ! assuming that all k-point have the same weight, i.e.
      NORM=4._mq*REAL(PI,mq)*(HARTREE**3)/CELLVOL*REAL(KPTWEIGHT(1),kind=mq)
! additional spin-degeneracy factor 2 for non-spin polarized calculations ...
      IF (ISPIN==1) NORM=2._mq*NORM
! if you like to see the value of the norm factor uncomment these two lines:
!      IF(IO%IU6>=0) WRITE(IO%IU6,*)'NORM:',NORM
!      IF(IO%IU0>=0) WRITE(IO%IU0,*)'NORM:',NORM

! initialize all to zero
      EPSILON_TOT=CZERO
      EPSILON_RE=0.0_q
      EPSILON_IM=0.0_q

!==============================================================================
! Call the "kernel routine" doing the iterations and projections of the vectors;
! here one has to use MATDIM_BTD (which is either MATDIM [for BEYONDTD=0] or
! 2*MATDIM [for BEYONDTD=1] and therefore always has the correct dimension)
!==============================================================================
      CALL TIME_DEVELOP_EXCI(MATDIM_BTD,IO,BSEMATRIX,OPTMAT)

!=====================
! now calculate the DF
!=====================
      TIMESTEP=TIMESTEP_V(1,1)


! According to Eq. 18 we have now to perform a Fourier transformation in order
! to get spectra from the "real-time signal" generated by "TIME_DEVELOP_EXCI"
! (actually it is then numerically a sum instead of an integral ... ;-)).
! Note: It's a "slow" Fourier transform -- fortunately (non-uniform time grid!);
!       currently this is not (yet) parallelized (but it's order(N^0) ... [!])
      IF (NODE_ME==NODE_LEAD) THEN
         !$ACC ENTER DATA COPYIN(TIMESTEP_V,EPSILON_TOT,EPSILON_RE,EPSILON_IM) __IF_ASYNC__
         !$ACC KERNELS PRESENT(TIMESTEP_V,SCALAR,EPSILON_TOT) __IF_ASYNC__
         DO LDIR=1,3
          DO KDIR=1,3
            DO I=0,MAXOMEGA
               IF (KDIR==MODULO(LDIR,3)+1) CYCLE
! actual frequency "omega"
               IF (LFIXED) THEN
                  W=ASIN(REAL(I,mq)*DELTAW*TIMESTEP)/TIMESTEP
                  WEIGHT=1._mq/SQRT(1._mq - TIMESTEP**2*(REAL(I,mq)*DELTAW)**2)
               ELSE
                  W=DELTAW*REAL(I,mq)
                  WEIGHT=1._mq
               ENDIF
! polarizability alpha (needs to be initialized to zero here for the summation)
               ALPHA=CZERO
               DO J=1,LAST_TDSTEP
! Set EXPFAC to the exponential-factor appearing in Eq. 18
! (with TIMESTEP(2,..) containing the integral time t)
                  EXPFAC=EXP((IMUN*W-GAM)*TIMESTEP_V(2,J))
! the "dt", set like for trapezoidal integration rule (for non-uniform grids);
! in the original code it was implicitly "dt=timestep_v(1,j)=dt_p" (= uniform
! grid) but now we should properly account for the non-uniformity of the grid
!
! Additional note: Various tests with alternative integration schemes (e.g.
! some generalized Simpson rule, adapted to the case of non-uniform grids)
! uncovered that at the end we can do whatever we like: the very primitive
! trapezoidal rule always gives "identical results" (crudely spoken plus/minus
! machine accuracy ...) -- so, at the end the decision was to stick to the
! trapezoidal rule (adapted for non-uniform meshes) as "fully sufficient" ...

! first set up the "dt" for the "left boundary" (t=0=t(J=1))
                  IF (J==1) THEN
                     DT=0.5_mq*TIMESTEP_V(1,1)
! "right boundary" (t=t(J=last_tdstep)) -- mind: timestep "J" holds the time
! difference t(J+1)-t(J) and therefore we have to use argument "J-1" here ... !
                  ELSE IF (J==LAST_TDSTEP) THEN
                     DT=0.5_mq*TIMESTEP_V(1,J-1)
! any point "inside" -- again we need to replace "2*dt" by the sum "dt_m+dt_p"
                  ELSE
                     DT=0.5_mq*(TIMESTEP_V(1,J-1)+TIMESTEP_V(1,J))
                  ENDIF

! Now something tricky: Read Eq. 18 carefully and think about what it says! The
! quantity   [ <v_j|Xi_j> - <v_j|Xi_j>^* ]   cancels out the real part (due to
! the minus sign) but since "- -" is "+" it gives two times the imaginary part!
! But this multiplied by "i" and the resulting "i^2" finally gives a "-" sign!
! Hence we can simply use   -2._mq*AIMAG(SCALAR(...))"   now; but now caution!
!
! Looking closely to Eq. 18 you will recognize that it just states something
! about the diagonal elements of polarizability "alpha" (and of "EPSILON_TOT");
! but of course, we do not only want the diagonal elements -- we want  "ALL"  !
! It should be obvious to everybody that the simple generalization of Eq. 18
! for  NON-DIAGONAL  elements "alpha_ij" is just to replace the expression
! [ <v_j|Xi_j> - <v_j|Xi_j>^* ]   by    [ <v_i|Xi_j> - <v_i|Xi_j>^* ]  (!).
! But then still the same "construct" applies as for the diagonal elements!
                  IF (BEYONDTD==0) THEN
                     ALPHA=ALPHA-2.0_mq*AIMAG(SCALAR(J,KDIR,LDIR))*EXPFAC*DT*WEIGHT
                  ELSE IF (BEYONDTD==1) THEN
! And again a little puzzle -- here things go a bit different: coming back to
! Eq. 18 from the paper, the "<v_j|Xi_j> - <v_i|Xi_j>^*" above contains now
! coupling terms (it is the expression "BDOTC(first part) - BDOTC(second part)"
! in the projections above). But only in Tamm-Dancoff approximation for the
! second term (antiresonant part) "anti-resonant" = "minus conjg(resonant)"
! holds (or in the language of "beyond Tamm-Dancoff" we have only calculated
! "BDOTC(first term)" and have complememented the expression just here with the
! yet missing "- BDOTC(second term) = - CONJG(BDOTC(first term))" for the
! Tamm-Dancoff approximation. Including all couplings ("beyond Tamm-Dancoff")
! one has now to understand that the second term "- BDOTC(second part)" in the
! projections above is actually what corresponds then to "- <v_j|Xi_j>^*" in
! Tamm-Dancoff approximation (but different from "- CONJG(BDOTC(first term))".
! Hence, we just have to take the plain projection (SCALAR) already containing
! all we need and have just to do the multiplicatio with imaginary unit "i"
! (no chance for further simplification / reductions like above ... !), i.e.,
                     ALPHA=ALPHA+IMUN*SCALAR(J,KDIR,LDIR)*EXPFAC*DT*WEIGHT
                  ENDIF
               ENDDO ! time steps
! Finally, apply the norm factor in front of the integral in Eq. 18 to "alpha"
! (and store the result directly into array "EPSILON_TOT" ...):
               EPSILON_TOT(I,KDIR,LDIR)=NORM*ALPHA
! Going from "alpha" to the dielectric function "EPSILON_TOT" we have to add a
! unit matrix (means to add a value of "one" to all   diagonal   elements).
            ENDDO    !frequencies
          ENDDO      !kdir
         ENDDO       !ldir
         !$ACC END KERNELS

! at this point we will now only write to vasprun.xml like for IBSE=0, the code
! writing out to the old "spectra file" (some relict from pre-vasprun.xml times)
! is now definitely gone with this version (not needed and not used anymore!)

         IF (ANY(ABS(QPOINT)>1E-6_q)) THEN
            !$ACC KERNELS PRESENT(EPSILON_TOT) __IF_ASYNC__
            EPSILON_TOT=EPSILON_TOT/(SUM(QPOINT**2)*TPI**2)
            !$ACC END KERNELS
            !$ACC PARALLEL LOOP PRESENT(EPSILON_TOT) __IF_ASYNC__
            DO I = 1,3
               EPSILON_TOT(:,I,I) = 1/(1-EPSILON_TOT(:,I,I))
            ENDDO
         ELSE
            !$ACC PARALLEL LOOP PRESENT(EPSILON_TOT) __IF_ASYNC__
            DO I=1,3
               EPSILON_TOT(:,I,I)=EPSILON_TOT(:,I,I)+CONE
            ENDDO
         ENDIF

!================================
! prepare epsilon for vasprun.xml
!================================
         !$ACC KERNELS PRESENT(EPSILON_TOT,EPSILON_RE,EPSILON_IM) __IF_ASYNC__
         DO LDIR=1,3
          DO KDIR=1,3
! just to explain it at least here (this construct was already used previously
! several times ...): the dielectric tensor is symmetric and hence for output
! we only need the element 11, 22, 33, 12, 23, and 31 but NOT 21, 32, or 13
! and this here is constructed in such a way that KDIR,LDIR = 21, 32, or 13
! will be cycled (filtered out as "unnecessary elements") -- and at the end we
! do not even need to "complete" the arrays handed over to the output routine
! (since the yet "missing" elements still containing zeros are anyway ignored!)
           IF (KDIR==MODULO(LDIR,3)+1) CYCLE
           DO I=0,NEDOS-1
! for single_prec_bse we would have all quantities in single precision but the
! output routine strictly requires double precision (KIND "_q") arguments ...;
! just for that reason we create now copies into arrays EPSILON_RE/EPSILON_IM
! which have been declared "REAL(q)" and can be handed over to "XML_EPSILON_W"
            EPSILON_RE(I+1,KDIR,LDIR)=REAL(EPSILON_TOT(I,KDIR,LDIR),mq)
            EPSILON_IM(I+1,KDIR,LDIR)=AIMAG(EPSILON_TOT(I,KDIR,LDIR))
           ENDDO
          ENDDO
         ENDDO
         !$ACC END KERNELS

! write the dielectric function into "vasprun.xml" -- but only NEDOS elements!
! note jF: the factor "1._q" enforces conversion to double precision if DELTAW
! is single precision only (XML_EPSILON_W requires first argument of type "_q"!)
! SO DO NOT TRY TO "OPTIMIZE" THIS (i.e., do NOT remove the "1._q") OR YOU FAIL!
         !$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
         !$ACC EXIT DATA COPYOUT(EPSILON_RE,EPSILON_IM) DELETE(TIMESTEP_V,EPSILON_TOT) __IF_ASYNC__
         IF (ANY(ABS(QPOINT)>1E-6_q)) THEN
            CALL XML_EPSILON_SCALAR_W(1._q*DELTAW, EPSILON_RE(:,1,1), EPSILON_IM(:,1,1), NEDOS)
         ELSE
            CALL XML_EPSILON_W(1._q*DELTAW, EPSILON_RE(:,:,:), EPSILON_IM(:,:,:), NEDOS)
         ENDIF
#ifdef VASP_HDF5
         IF (LGWLF) THEN
           PREFIX = "bse"
         ELSE
           PREFIX = "tdhf"
         END IF
         EPSILON = CMPLX(EPSILON_RE,EPSILON_IM)
         CALL VH5_WRITE_DIELECTRIC_ENERGIES(IH5OUTFILEID, 1._q*DELTAW, NEDOS, PREFIX)
         IF (ANY(ABS(QPOINT)>1E-6_q)) THEN
            CALL VH5_WRITE_DIELECTRIC_SCALAR_DYNAMIC(IH5OUTFILEID, EPSILON(1:NEDOS,1,1), PREFIX)
         ELSE
            CALL VH5_WRITE_DIELECTRIC_DYNAMIC(IH5OUTFILEID, EPSILON, PREFIX)
         ENDIF
#endif
      ENDIF    ! node_me

      infwritef(*,'(1x,(a))') 'Calculation of spectra finished successfully. Output written to ''vasprun.xml'''

! clean up everything and ciao ...
      !$ACC EXIT DATA DELETE(SCALAR)  __IF_ASYNC__
      DEALLOCATE(EPSILON_TOT)
      DEALLOCATE(EPSILON)
      DEALLOCATE(EPSILON_RE)
      DEALLOCATE(EPSILON_IM)
      DEALLOCATE(TIMESTEP_V)
      DEALLOCATE(SCALAR)

      PROFILING_STOP('calc_and_write_spectra')

      END SUBROUTINE CALC_AND_WRITE_SPECTRA
#endif
END MODULE bse_te
