#include "symbol.inc"
!> This module provides routines to determine the Fermi energy.
!>
!> The basic algorithm to determine the Fermi energy works as follows:
!>
!> * We consider a single k-point and evaluate the number of occupied states
!>   per spin channel there. Here, we need to account for the fact that the
!>   states may be split over multiple CPUs. We use a bisection algorithm with
!>   step function smearing to do so.
!>
!> * Then we compute the band gap of the material. To this end, we assume the
!>   number of occupied states is the same at all k points and evaluate the
!>   difference in the eigenvalues between the valence band maximum and the
!>   conduction band minimum. Note that for a metal the assumption of constant
!>   occupation at all k points is false and we will obtain a negative band gap
!>   that we set to zero.
!>
!> * If the band gap is not zero, we put the Fermi energy in the middle of the
!>   gap and check the number of electrons. If the number of electrons agrees
!>   with the expected value, we found the correct Fermi energy and return.
!>
!> * For metals or when the smearing lead to an incorrect number of electrons,
!>   we determine the Fermi energy with a bisection method between the band
!>   extrema. We compute the number of electrons for a Fermi energy in the
!>   middle of the interval. When the result is smaller or larger than the
!>   expected number of electrons, we replace the lower or upper interval
!>   boundary, respectively.
!>
!>
!> The alternative method to evaluate the Fermi energy traditionally used by
!> VASP worked differently:
!>
!> * We compute the integrated DOS using the smearing method of choice.
!>
!> * We set the Fermi energy to where the lowest energy where the integrated
!>   DOS approximately agrees with the number of electrons or is larger than it.
!>
!> * Then we start a bisection algorithm in the vicinity of this approximate
!>   Fermi energy. We calculate the Fermi weights for all bands and k-points and
!>   the number of electrons as sum of all of those. Then we update the interval
!>   analogously to the algorithm above.
module fermi_energy

    use prec, only: q
    use wave, only: wavedes

    implicit none

    !> Derive new smearings from this base class so that the determination of
    !> the Fermi energy works.
    type, abstract, public :: smearing_base
        !> used to reduce the data over the communicators
        type(wavedes), pointer :: wdes => null()
        !> number of electrons fitting into a single band
        real(q) spin_degeneracy
        !> the electronic eigenvalues
        real(q), allocatable :: band_structure(:,:,:)
        !> tightest attainable convergence in the determination
        !> of the number of electrons
        real(q) :: attainable_convergence_estimate = 0
    contains
        procedure :: initialize
        procedure :: comm => get_communicator
        procedure(number_electrons_interf), deferred :: number_electrons
    end type smearing_base

    abstract interface
        !> Every smearing subclass must define a function like this.
        !>
        !> The bisection algorithm terminates when the return value of this
        !> function agrees with the desired number of electrons
        real(q) function number_electrons_interf(this, fermi_energy)
            import q, smearing_base
            class(smearing_base), intent(in) :: this
            real(q), intent(in) :: fermi_energy  !< trial value for the Fermi energy
        end function number_electrons_interf
    end interface

    !> Implements the Gaussian and Methfessel-Paxtion smearing.
    type, extends(smearing_base), public :: default_smearing
        integer ismear
        real(q) sigma
    contains
        procedure :: number_electrons => default_number_electrons
    end type default_smearing

    private
    public calculate_fermi_energy, calculate_fermi_energy_legacy, calculate_fermi_weights
    public delstp

    integer, parameter :: min_ = 1, max_ = 2
    integer, parameter :: vbm = max_, cbm = min_
    real(q), parameter :: meV_precision = 1e-3_q
    real(q), parameter :: converged_number_electrons = 1e-10_q

    !> If the smearing is too small, the code may be numerically instable.
    real(q), parameter, public :: too_small_sigma = 1e-5_q

contains

    !> Initialize the common fields that all smearing subclasses have.
    subroutine initialize(this, band_structure, spin_degeneracy, wdes)
        class(smearing_base), intent(inout) :: this
        !> The electronic eigenvalues.
        real(q), intent(in) :: band_structure(:,:,:)
        !> Number of electrons per band.
        real(q), intent(in) :: spin_degeneracy
        !> If provided the algorithms will run in parallel over the communicators
        !> defined in the descriptor.
        type(wavedes), intent(in), target, optional :: wdes
        if (present(wdes)) this%wdes => wdes
        allocate(this%band_structure, source=band_structure)
        this%spin_degeneracy = spin_degeneracy
    end subroutine initialize


    type(communic) function get_communicator(smearing) result (comm)
        use mpimy
        class(smearing_base), intent(in) :: smearing
        if (associated(smearing%wdes)) then
            comm = smearing%wdes%comm
        else
#if defined(MPI) || defined(MPI_CHAIN)
            comm%mpi_comm = MPI_COMM_SELF
            call M_initc(comm)
#endif
        end if
    end function get_communicator


    !> Determine a Fermi energy, for which the desired number of electrons fit
    !> into the band structure.
    !>
    !> Initially, we try to find a semiconducting solution by setting the occupancy
    !> to the same values as for the first **k** point. If this leads to a Fermi
    !> energy in the gap producing the correct number of electrons, we are done.
    !> Otherwise, we bisect the available spectrum to find the energy at which
    !> the correct number of electrons are in the bands.
    real(q) function calculate_fermi_energy(smearing, number_electrons) result(fermi_energy)
        !> smearing method used to calculate the number of electrons
        class(smearing_base), intent(in) :: smearing
        !> desired number of electrons
        real(q), intent(in) :: number_electrons
        integer number_occupied(size(smearing%band_structure, 3))
        real(q) band_edges(2), error_number_electrons
        !
        number_occupied = first_kpoint_occupations(smearing, number_electrons)
        band_edges = find_band_edges(smearing, number_occupied)
        fermi_energy = 0.5 * sum(band_edges)
        error_number_electrons = number_electrons - smearing%number_electrons(fermi_energy)
        if (error_is_acceptable(smearing, error_number_electrons)) &
            return
        fermi_energy = bisect_band_structure(smearing, number_electrons)
    end function calculate_fermi_energy


    function first_kpoint_occupations(smearing, number_electrons) result (number_occupied)
        class(smearing_base), intent(in) :: smearing
        real(q), intent(in) :: number_electrons
        integer number_occupied(size(smearing%band_structure, 3))
        integer num_elec
        real(q) bounds(2), fermi_energy
        integer spin
        !
        if (no_bisection_required(smearing)) then
            ! trivial case - only one spin channel on a single cpu
            number_occupied = nint(number_electrons / smearing%spin_degeneracy)
            return
        end if
        !
        bounds(min_) = minval(smearing%band_structure)
        bounds(max_) = maxval(smearing%band_structure)
        call reduce_over_comm(smearing, bounds)
        !
        do
            fermi_energy = 0.5_q * sum(bounds)
            do spin = 1, size(number_occupied)
               number_occupied(spin) = count(smearing%band_structure(:, 1, spin) < fermi_energy)
            end do
            num_elec = smearing%spin_degeneracy * sum(number_occupied)
            call sum_over_comm(smearing, num_elec)
            if (num_elec == nint(number_electrons)) then
                exit
            else if (bounds(max_) - bounds(min_) < meV_precision) then
                exit
            else if (num_elec < nint(number_electrons)) then
                bounds(min_) = fermi_energy
            else
                bounds(max_) = fermi_energy
            end if
        end do
    end function first_kpoint_occupations


    function find_band_edges(smearing, number_occupied) result (band_edges)
        class(smearing_base), intent(in) :: smearing
        integer, intent(in) :: number_occupied(:)
        real(q) band_edges(2), band_edges_spin(2, size(number_occupied))
        real(q), parameter :: not_set = -huge(1.0_q) ! set too large negative so that maxval works
        integer spin
        !
        band_edges_spin = not_set
        do spin = 1, size(number_occupied)
            if (number_occupied(spin) > 0) &
                band_edges_spin(vbm, spin) = maxval(smearing%band_structure(number_occupied(spin), :, spin))
            if (number_occupied(spin) < size(smearing%band_structure, 1)) &
                band_edges_spin(cbm, spin) = minval(smearing%band_structure(number_occupied(spin) + 1, :, spin))
        end do
        band_edges(vbm) = maxval(band_edges_spin(vbm, :))
        band_edges(cbm) = minval(band_edges_spin(cbm, :))
        call reduce_over_comm(smearing, band_edges)
    end function find_band_edges


    real(q) function bisect_band_structure(smearing, number_electrons) result (fermi_energy)
        use string, only: str
        use tutor, only: vtutor
        class(smearing_base), intent(in) :: smearing
        real(q), intent(in) :: number_electrons
        real(q) error_number_electrons
        real(q) bounds(2)
        !
        bounds(min_) = minval(smearing%band_structure)
        bounds(max_) = maxval(smearing%band_structure)
        call reduce_over_comm(smearing, bounds)
        do
            fermi_energy = 0.5_q * sum(bounds)
            error_number_electrons = number_electrons - smearing%number_electrons(fermi_energy)
            if (error_is_acceptable(smearing, error_number_electrons)) &
                return
            if ((bounds(2) - bounds(1)) / (abs(fermi_energy) + converged_number_electrons) < 1E-14_q) &
                exit
            if (error_number_electrons > 0.0_q) then
                bounds(min_) = fermi_energy
            else
                bounds(max_) = fermi_energy
            end if
        end do
        !
        call vtutor%warning("Midgap Fermi energy: can't reach specified precision\n&
            &Number of Electrons is NELECT = " // str(number_electrons - error_number_electrons))
        !
    end function bisect_band_structure


    logical function error_is_acceptable(smearing, error_number_electrons)
        class(smearing_base), intent(in) :: smearing
        real(q), intent(in) :: error_number_electrons
        real(q) acceptable_error
        acceptable_error = max(smearing%attainable_convergence_estimate * 2, converged_number_electrons)
        error_is_acceptable = abs(error_number_electrons) <= acceptable_error
    end function error_is_acceptable


    logical function no_bisection_required(smearing)
        class(smearing_base), intent(in) :: smearing
        logical distributed, spin_polarized
        distributed = .false.
        if (associated(smearing%wdes)) distributed = smearing%wdes%comm%ncpu > 1
        spin_polarized = size(smearing%band_structure, 3) > 1
        no_bisection_required = .not.distributed .and. .not.spin_polarized
    end function no_bisection_required


    subroutine sum_over_comm(smearing, num_elec)
        class(smearing_base), intent(in) :: smearing
        integer, intent(inout) :: num_elec
        if (.not.associated(smearing%wdes)) return
        CALLMPI( M_sum_i(smearing%wdes%comm, num_elec, 1) )
        num_elec = num_elec / (smearing%wdes%comm_inb%ncpu * smearing%wdes%comm_kinter%ncpu)
    end subroutine sum_over_comm


    subroutine reduce_over_comm(smearing, bounds)
        class(smearing_base), intent(in) :: smearing
        real(q), intent(inout) :: bounds(2)
        if (.not.associated(smearing%wdes)) return
        ! multiply lower bound with -1 so that maximum across CPUs works
        bounds(min_) = -bounds(min_)
        CALLMPI( M_max_d(smearing%wdes%comm, bounds, size(bounds)) )
        bounds(min_) = -bounds(min_)
    end subroutine reduce_over_comm


    real(q) function default_number_electrons(this, fermi_energy) result (number_electrons)
       use tutor, only: vtutor
       class(default_smearing), intent(in) :: this
       real(q), intent(in) :: fermi_energy
       real(q), allocatable :: weights(:,:,:)
       integer kpoint, kfirst, kstep
       !
       if (.not.associated(this%wdes)) call vtutor%bug("Assertion failed that WDES is not null pointer", __FILE__, __LINE__)
       allocate(weights, mold=this%band_structure)
       call calculate_fermi_weights(this, fermi_energy, weights)
       !
       kfirst = this%wdes%comm_kinter%node_me
       kstep = this%wdes%comm_kinter%ncpu
       number_electrons = 0
       do kpoint = kfirst, this%wdes%nkpts, kstep
          number_electrons = number_electrons + sum(weights(:,kpoint,:)) * this%wdes%wtkpt(kpoint)
       end do
       CALLMPI( M_sum_d( this%wdes%comm, number_electrons, 1) )
       number_electrons = number_electrons * this%spin_degeneracy / this%wdes%comm_inb%ncpu
    end function default_number_electrons


    subroutine calculate_fermi_weights(this, fermi_energy, weights)
       use tutor, only: vtutor
       class(default_smearing), intent(in) :: this
       real(q), intent(in) :: fermi_energy
       real(q), intent(out) :: weights(:,:,:)
       integer spin, kpoint, kfirst, kstep, band
       real(q) xx, ignored, weight
       !
       if (.not.associated(this%wdes)) call vtutor%bug("Assertion failed that WDES is not null pointer", __FILE__, __LINE__)
       weights = 0
       kfirst = this%wdes%comm_kinter%node_me
       kstep = this%wdes%comm_kinter%ncpu
       do spin = 1, size(this%band_structure, 3)
          do kpoint = kfirst, this%wdes%nkpts, kstep
             do band = 1, this%wdes%nbands
                xx = (fermi_energy - this%band_structure(band, kpoint, spin)) / this%sigma
                call delstp(this%ismear, xx, ignored, weights(band, kpoint, spin))
             end do
          end do
       end do
    end subroutine calculate_fermi_weights


    !> Determine a Fermi energy, for which the desired number of electrons fit
    !> into the band structure.
    !>
    !> Evaluate a first guess for the Fermi energy based on the integrated density
    !> of states. Then increase or decrease the Fermi energy until the error until
    !> the error in the number of electrons changes sign. Finally, bisect in this
    !> interval to find the correct Fermi energy.
    real(q) function calculate_fermi_energy_legacy(smearing, number_electrons, &
            emin, delta_e, int_dos, sigma_, efermi_force, routine) result (fermi_energy)
        use string, only: str
        use tutor, only: vtutor
        use wave_struct_def, only: wavedes
        class(smearing_base), intent(in) :: smearing
        real(q), intent(in) :: emin, delta_e, int_dos(:,:), number_electrons, sigma_, efermi_force
        character(len=*), intent(in) :: routine
        !
        real(q) bounds(2), num_elec
        logical enlarge_bounds(2)
        integer ii
        !
        if (efermi_force /= HUGE(1.0_q)) then
            fermi_energy = efermi_force
            return
        end if
        !
        ! calculate approximated fermi energy
        do ii = 1, size(int_dos, 1)
           num_elec = sum(int_dos(ii,:))
           if (abs(num_elec - number_electrons) < 0.01_q &
              .or. num_elec > number_electrons) exit
        end do
        fermi_energy = emin + (ii - 1) * delta_e
        if (sigma_ < too_small_sigma) return
        !
        ! search now for exact Fermi-level using bisectioning
        bounds(min_) = emin + (ii - 2) * delta_e
        bounds(max_) = emin + (ii - 1) * delta_e
        enlarge_bounds = .true.
        do
            fermi_energy = 0.5 * sum(bounds)
            num_elec = smearing%number_electrons(fermi_energy)
            ! compare with number of electrons
            if (error_is_acceptable(smearing, num_elec - number_electrons)) &
                return
            if ((bounds(2) - bounds(1)) / (abs(fermi_energy) + converged_number_electrons) < 1E-14_q) &
                exit
            !
            if (num_elec > number_electrons) then
                if (enlarge_bounds(min_)) bounds(min_) = bounds(min_) - delta_e
                bounds(max_) = fermi_energy
                enlarge_bounds(max_) = .false.
            else
                if (enlarge_bounds(max_)) bounds(max_) = bounds(max_) + delta_e
                bounds(min_) = fermi_energy
                enlarge_bounds(min_) = .false.
            end if
        end do
        !
        call vtutor%warning(routine // ": can't reach specified precision\n&
            &Number of Electrons is NELECT = " // str(num_elec))
        !
    end function calculate_fermi_energy_legacy


!******************** DELSTP    ****************************************
!
!> Returns generalised delta and step functions (Methfessel & Paxton)
!>
!>  Input:
!>      n >=0 : order of approximant; x : argument
!>      for n=-1 Fermi-function is used
!>      for n<-1 Gaussian broadening is used
!>  Output:
!>      D_n (x) ,  S_n (x)
!>  Remarks:
!>      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
!>      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
!>      where H is a Hermite polynomial and
!>      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
!
!***********************************************************************

      SUBROUTINE DELSTP(N,X,D,S)
GPACC !$ACC ROUTINE SEQ
GPACC USE moffload_struct_def, ONLY : ACC_ON_DEVICE,ACC_DEVICE_NOT_HOST
GPOMP USE omp_lib, ONLY : OMP_IS_INITIAL_DEVICE

      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF

      IMPLICIT REAL(q) (A-H,O-Z)
GPOMP !$OMP DECLARE TARGET (DELSTP)

      INTEGER :: N, I, K

      IF (N == -1) THEN
         OUT_OF_BOUNDS = 50.0_q
      ELSE
         OUT_OF_BOUNDS = 20.0_q
      END IF

      IF (X < -OUT_OF_BOUNDS) THEN
         D = 0._q
         S = 0._q
         RETURN
      END IF
      IF (X > OUT_OF_BOUNDS) THEN
         D = 0._q
         S = 1._q
         RETURN
      END IF
!=======================================================================
!  If n == -1: use Fermi-Dirac statistics
!=======================================================================
      IF (N==-1) THEN
!         S=  1/(1 + EXP(-X))
!         ! hyperbolic function is nummerically safer
!         !D= EXP(-X)/(1 + EXP(-X))**2
!         ! is
!         D = 1/(4*COSH(X/2)**2)
         !
         ! use a robust determination of S = f(-x)=1/(1+EXP(-x))
         !
         S = 1 / ( 1 + EXP(-X ) )
         D = S*(1-S)
         RETURN
!=======================================================================
!  If n < 0 : assume Gaussian type smearing
!  (must return  same as N=0 or ... )
!=======================================================================
      ELSE IF (N<0) THEN
         D=EXP(-(X*X))/SQRT(PI)
#ifdef OFFLOADING
#ifdef IFXWAR5
         !Compiler bug for IFX (<= 2025.3.2) ERF does not compile. custom function works.
         S=0.5_q+0.5_q*ERRF(X)
#else // IFXWAR5
GPACC    IF(ACC_ON_DEVICE(ACC_DEVICE_NOT_HOST)) THEN
GPOMP    IF(.NOT.OMP_IS_INITIAL_DEVICE()) THEN
            S=0.5_q+0.5_q*ERF(X)
         ELSE
            S=0.5_q+0.5_q*ERRF(X)
         ENDIF
#endif // IFXWAR5
#else // OFFLOADING
         S=0.5_q+0.5_q*ERRF(X)
#endif // OFFLOADING
         RETURN
      END IF
!=======================================================================
! Methfessel & Paxton
!=======================================================================
      EX2=EXP(-(X*X))
#ifdef OFFLOADING
#ifdef IFXWAR5
      !Compiler bug for IFX (<= 2025.3.2) ERF does not compile. custom function works.
      S0=0.5_q*ERRF(X)
#else // IFXWAR5
GPACC IF(ACC_ON_DEVICE(ACC_DEVICE_NOT_HOST)) THEN
GPOMP IF(.NOT.OMP_IS_INITIAL_DEVICE()) THEN
          S0=0.5_q*ERF(X)
      ELSE
          S0=0.5_q*ERRF(X)
      ENDIF
#endif // IFXWAR5
#else // OFFLOADING
      S0=0.5_q*ERRF(X)
#endif // OFFLOADING
      A=1._q/SQRT(PI)
      K=0
      H1=1._q
      H2=2._q*X
      S=0._q
      D=A
      DO I=1,N
         A=A/((-4._q)*I)
         K=K+1
         H3=H1
         H1=H2
         H2=2._q*X*H2-2*K*H3
         S=S+A*H1
         K=K+1
         H3=H1
         H1=H2
         H2=2._q*X*H2-2*K*H3
         D=D+A*H1
      ENDDO
      D=D*EX2
      S=0.5_q+S0-S*EX2
      RETURN
      END SUBROUTINE DELSTP

end module fermi_energy
