!NEC$ options "-no-ftrace"

#include "symbol.inc"
!===========================================================================
! RCS:  $Id: random.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!> This random number generator originally appeared in Toward a Universal
!> Random Number Generator by George Marsaglia and Arif Zaman.
!> Florida State University Report: FSU-SCRI-87-50 (1987)
!>
!> It was later modified by F. James and published in A Review of Pseudo-
!> random Number Generators
!>
!> Some final small modifications have been done by J. Furthmueller
!> Technical University of Vienna, November 1993
!>
!> THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
!>       (However, a newly discovered technique can yield
!>         a period of 10^600. But that is still in the development stage.)
!>
!> It passes ALL of the tests for random number generators and has a period
!>   of 2^144, is completely portable (gives bit identical results on all
!>   machines with at least 24-bit mantissas in the floating point
!>   representation).
!>
!> The algorithm is a combination of a Fibonacci sequence (with lags of 97
!>   and 33, and operation "subtraction plus one, modulo one") and an
!>   "arithmetic sequence" (using subtraction).
!>
! On a Vax 11/780, this random number generator can produce a number in
!    13 microseconds.
! (Note by J. Furthmueller: in 2.5 microseconds on a IBM RS6000/Model 580)
!========================================================================
module random_seeded
    use prec, only: q, qi4, qi8
    use constant

    implicit none

    real(q), parameter :: twopi = 6.283185307179586_q
    real(q), parameter :: cd = 7654321.0_q / 16777216.0_q
    real(q), parameter :: cm = 16777213.0_q / 16777216.0_q
    integer, parameter :: max_ij = 31328, max_kl = 30081
    integer, parameter :: countlarge = 1000000
    integer, parameter :: seed1_max = 900000000, seed2_max = 1000000

    type, public :: random_number_generator
        logical, private :: initialized = .false.
        real(q), private :: u(97)
        real(q), private :: c = 362436.0_q / 16777216.0_q
        integer, private :: i97 = 97, j97 = 33

    contains
        generic, public :: random => random_real_scalar, random_real_array
        generic, public :: uniform => uniform_real_scalar, uniform_real_array, &
            uniform_complex_scalar, uniform_complex_array
        procedure, private :: initialize
        procedure, private :: random_real_scalar, random_real_array
        procedure, private :: uniform_real_scalar, uniform_real_array, &
            uniform_complex_scalar, uniform_complex_array
    end type random_number_generator

    type(random_number_generator), public :: gen_wave, gen_ion


    integer(qi8), parameter :: pcg_multiplier = 6364136223846793005_qi8
    real(q), parameter :: resolution = 2.3283064365386963e-10_q  ! 1 / 2**32

    type pcg_generator
        integer(qi8) :: state = 0, increment
    end type pcg_generator


    integer, parameter, public :: default_rng = 1
    integer, parameter, public :: pcg_32bit_rng = 2

    private
    public :: rane_ion, rang_ion, rane, rang, rane_array, pcg_generator, new_pcg_generator, random_number, &
        random_reader, write_random

    interface random_number
        module procedure pcg_integer_scalar, pcg_integer_array, pcg_real_scalar, pcg_real_array, &
            pcg_complex_scalar, pcg_complex_array
    end interface random_number

contains

    subroutine random_reader(incar, comm, seed, generator, pcg_seed)
        use incar_reader, only: incar_file, process_incar, not_found_error
        use mpimy, only: communic
        use string, only: lowercase
        use tutor, only: vtutor, isError, RandomSeed, newArgument
        !
        type(incar_file), intent(inout) :: incar
        type(communic), intent(in) :: comm
        integer, intent(out) :: seed(3), generator
        integer(qi8), intent(out) :: pcg_seed
        character(len=:), allocatable :: random_generator
        integer(qi4) pcg_seed_copy(2)
        integer ierr
        !
        seed = 0
#ifdef tbdyn
        call process_incar(incar, 'RANDOM_SEED', seed, ierr)
        if (ierr == not_found_error) then
            ! with this choice we are sure not to exceed SEED1_MAX (maximal
            ! value we can get here is 893581000)
            call seed_from_time(seed(1))
        else if ((seed(1) < 0) .or. (seed(1) > seed1_max) .or. &
                 (seed(2) < 0) .or. (seed(2) > seed2_max) .or. &
                 (seed(3) < 0)) then
            call vtutor%write(isError, RandomSeed, newArgument(ival = seed))
        end if
        CALLMPI( M_bcast_i(comm, seed, size(seed)) )
#endif
        random_generator = 'default'
        call process_incar(incar, 'RANDOM_GENERATOR', random_generator)
        random_generator = lowercase(random_generator)
        select case (random_generator)
        case ('default')
            generator = default_rng
        case ('pcg_32')
            generator = pcg_32bit_rng
        case default
            call vtutor%alert("RANDOM_GENERATOR = " // trim(random_generator) // " could not be mapped onto &
                &a known random number. The default will be used instead. Please make sure there is no &
                &spelling error in the INCAR file.")
            generator = default_rng
        end select
        !
        if (generator == pcg_32bit_rng) then
            pcg_seed_copy = 0
            call process_incar(incar, 'PCG_SEED', pcg_seed_copy, ierr)
            if (ierr == not_found_error) &
                call seed_from_dev_random(pcg_seed_copy)
            CALLMPI( M_bcast_i(comm, pcg_seed_copy, size(pcg_seed_copy)) )
            pcg_seed = transfer(pcg_seed_copy, 1_qi8)
        end if
        !
    end subroutine random_reader

    subroutine seed_from_dev_random(seed)
        integer(qi4), intent(out) :: seed(2)
        integer iunit
        open(newunit=iunit, file='/dev/urandom', access='stream', form='unformatted', err=999)
        read(iunit, err=998) seed
        close(iunit)
        return
        ! fall back to reading from time
998     close(iunit)
999     continue
        call seed_from_time(seed(1))
    end subroutine seed_from_dev_random

    subroutine seed_from_time(seed)
        integer, intent(out) :: seed
        integer values(8)
        call date_and_time(values=values)
        seed = values(3)*24*60*60*300 + values(5)*60*60*1000 + values(6)*60*1000 + values(7)*1000 + values(8)
        seed = mod(seed, seed1_max)
    end subroutine seed_from_time

    subroutine write_random(iunit, generator, pcg_seed)
        use tutor, only: vtutor
        integer, intent(in) :: iunit, generator
        integer(qi8), intent(in) :: pcg_seed
        integer(qi4) pcg_seed_copy(2)
        !
        if (iunit < 0) return
        !
        write(iunit,'(a)') ' Random number generation:'
        select case (generator)
        case (default_rng)
            write(iunit,'(a)') '   RANDOM_GENERATOR = DEFAULT'
            write(iunit,'(a)') '   PCG_SEED         = not used'
        case (pcg_32bit_rng)
            write(iunit,'(a)') '   RANDOM_GENERATOR = PCG_32'
            pcg_seed_copy = transfer(pcg_seed, 1_qi4, size(pcg_seed_copy))
            write(iunit,'(a,2i12)') '   PCG_SEED         = ', pcg_seed_copy
        case default
            call vtutor%bug("Unknown random number generator", __FILE__, __LINE__)
        end select
        !
    end subroutine write_random

    real(q) function rane(reset)
        logical, intent(in), optional :: reset
        if (present(reset)) gen_wave%initialized = .not.reset
        call gen_wave%random(rane)
    end function rane

    subroutine rane_array(array, draw_n, reset)
        real(q), intent(out) :: array(:)
        integer, intent(in) :: draw_n
        logical, intent(in), optional :: reset
        if (present(reset)) gen_wave%initialized = .not.reset
        call gen_wave%random(array, draw_n)
    end subroutine rane_array

    !> This subroutine produces a normal distribution (Gaussian distribution)
    real(q) function rang(center, width)
        real(q), intent(in) :: center, width
        call gen_wave%uniform(rang, center, width)
    end function rang

    subroutine rane_ion(rng_ion, put, get)
        real(q), intent(out) :: rng_ion
        integer, intent(in), optional :: put(3)
        integer, intent(out), optional :: get(3)
        integer, save :: rng_state(3) = [54217137, 0, 0]
        integer, parameter :: seed1_max = 900000000
            ! The first state stores the seeds ij * (max_kl + 1) + kl, the other two
            ! the number of random numbers generated so far:
            ! rng_state(3) * countlarge + rng_state(2)
        if (present(get)) then
            ! just return the RNG seed (and initialize if needed)
            get = rng_state
        else if (present(put)) then
            ! ensure that seed(1) doesn't exceed max. sensible value
            rng_state(1) = mod(put(1), seed1_max)
            rng_state(2) = mod(put(2), countlarge)
            rng_state(3) = put(3) + put(2) / countlarge
            ! reseed and reinitialise if SEED is provided
            call reset(gen_ion, rng_state)
        else
            ! return random number if no special action is taken
            rng_state(2) = mod(rng_state(2) + 1, countlarge)
            if (rng_state(2) == 0) rng_state(3) = rng_state(3) + 1
            call gen_ion%random(rng_ion)
        end if
    end subroutine rane_ion

    !> This subroutine produces a normal distribution (Gaussian distribution)
    real(q) function rang_ion(center, width)
        real(q), intent(in) :: center, width
        call gen_ion%uniform(rang_ion, center, width)
    end function rang_ion

    subroutine reset(generator, rng_state)
        class(random_number_generator), intent(out) :: generator
        integer, intent(in) :: rng_state(3)
        integer seed_ij, seed_kl, ii, jj
        real(q) skip
        !c regenerate pair of seeds used originally
        seed_ij = rng_state(1) / (max_kl + 1)
        seed_kl = rng_state(1) - (max_kl + 1) * seed_ij
        call generator%initialize(seed_ij, seed_kl)
        ! skip rng_state3 * countlong + rng_state2 numbers (presumably already used
        ! in a previous run) on initialization
        do ii = 1, rng_state(3)
            do jj = 1, countlarge
                call generator%random(skip)
            end do
        end do
        do ii = 1, rng_state(2)
            call generator%random(skip)
        end do
    end subroutine reset

    !> This is the initialization routine for the random number generator RANMAR()
    !>
    !> NOTE: The seed variables can have values between:    0 <= IJ <= 31328
    !>                                                      0 <= KL <= 30081
    !> The random number sequences created by these two seeds are of sufficient
    !> length to complete an entire calculation with. For example, if several
    !> different groups are working on different parts of the same calculation,
    !> each group could be assigned its own IJ seed. This would leave each group
    !> with 30000 choices for the second seed. That is to say, this random
    !> number generator can create 900 million different subsequences -- with
    !> each subsequence having a length of approximately 10^30.
    !>
    !> Use IJ = 1802 & KL = 9373 to test the random number generator. The
    !> subroutine RANMAR should be used to generate 20000 random numbers.
    !> Then display the next six random numbers generated multiplied by 4096*4096
    !> If the random number generator is working properly, the random numbers
    !> should be:
    !>           6533892.0  14220222.0  7275067.0
    !>           6172232.0  8354498.0   10633180.0
    subroutine initialize(generator, seed_ij, seed_kl)
        use tutor, only: vtutor
        class(random_number_generator), intent(out) :: generator
        integer, intent(in) :: seed_ij, seed_kl
        if ((seed_ij < 0).or.(seed_ij > max_ij).or.(seed_kl < 0).or.(seed_kl > max_kl)) then
            call vtutor%error("The first random number seed must have a value between 0 and 31328 \n The &
               &second seed must have a value between 0 and 30081")
        end if
        call initialize_u(seed_ij, seed_kl, generator%u)
        generator%initialized = .true.
    end subroutine initialize

    pure subroutine initialize_u(seed_ij, seed_kl, u)
        integer, intent(in) :: seed_ij, seed_kl
        real(q), intent(out) :: u(97)
        integer i, j, k, l, m, ii, jj
        real(q) s, t
        !
        i = mod(seed_ij / 177, 177) + 2
        j = mod(seed_ij, 177) + 2
        k = mod(seed_kl / 169, 178) + 1
        l = mod(seed_kl, 169)
        !
        do ii = 1, size(u)
            s = 0.0_q
            t = 0.5_q
            do jj = 1, 24
                m = mod(mod(i * j, 179) * k, 179)
                i = j
                j = k
                k = m
                l = mod(53 * l + 1, 169)
                if (mod(l * m, 64) >= 32) then
                    s = s + t
                end if
                t = 0.5_q * t
            end do
            u(ii) = s
        end do
    end subroutine initialize_u

    !> This is the random number generator proposed by George Marsaglia in
    !> Florida State University Report: FSU-SCRI-87-50
    !> It was slightly modified by F. James to produce an array of pseudorandom
    !> numbers.
    subroutine random_real_scalar(generator, scalar)
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: scalar
        integer, parameter :: seed_ij = 1802, seed_kl = 9373
        real(q) uni
        if (.not.generator%initialized) call generator%initialize(seed_ij, seed_kl)
        uni = generator%u(generator%i97) - generator%u(generator%j97)
        if (uni < 0.0_q) uni = uni + 1.0_q
        generator%u(generator%i97) = uni
        generator%i97 = generator%i97 - 1
        if (generator%i97 == 0) generator%i97 = 97
        generator%j97 = generator%j97 - 1
        if (generator%j97 == 0) generator%j97 = 97
        generator%c = generator%c - cd
        if (generator%c < 0.0_q) generator%c = generator%c + cm
        uni = uni - generator%c
        if (uni < 0.0_q) uni = uni + 1.0_q
        scalar = uni
    end subroutine random_real_scalar

    subroutine random_real_array(generator, array, draw_n)
    class(random_number_generator), intent(inout) :: generator
    real(q), intent(inout) :: array(:)
    integer, intent(in) :: draw_n
    integer, parameter :: seed_ij = 1802, seed_kl = 9373
    integer :: ii, i97, j97
    real(q) :: uni, c

    if (.not. generator%initialized) call generator%initialize(seed_ij, seed_kl)

    ! Cache state locally
    i97 = generator%i97
    j97 = generator%j97
    c   = generator%c

    ! inline random_real_scalar for speed
    ! explicitly give size draw_n to function for flexibility
    ! but also confuses some compilers and optimize too aggressively if size(array) is used
    do ii = 1, draw_n
        uni = generator%u(i97) - generator%u(j97)
        if (uni < 0.0_q) uni = uni + 1.0_q

        generator%u(i97) = uni

        i97 = i97 - 1
        if (i97 == 0) i97 = 97

        j97 = j97 - 1
        if (j97 == 0) j97 = 97

        c = c - cd
        if (c < 0.0_q) c = c + cm

        uni = uni - c
        if (uni < 0.0_q) uni = uni + 1.0_q

        array(ii) = uni
    end do

    ! Write back final state
    generator%i97 = i97
    generator%j97 = j97
    generator%c   = c

    end subroutine random_real_array

    !> This subroutine produces a normal distribution (Gaussian distribution)
    !> using the built in fortran random number generator
    !> which is initialized with a random seed for tbdyn
    subroutine uniform_real_scalar(generator, scalar, center, width)
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: scalar
        real(q), intent(in), optional :: center, width
        real(q), parameter :: threshold = 1.0e-8_q
        real(q) first, second
        !
        call generator%random(first)
        do
            ! the random number generator includes zero, which leads to a problem
            ! in the logarithm, hence we replace it with a nonzero value
            call generator%random(second)
            if (second > threshold) exit
        end do
        scalar = cos(twopi * first) * sqrt(-2.0_q * log(second))
        if (present(width)) scalar = scalar * width
        if (present(center)) scalar = scalar + center
    end subroutine uniform_real_scalar

    subroutine uniform_real_array(generator, array, center, width)
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: array(:)
        real(q), intent(in), optional :: center, width
        integer ii
        do ii = 1, size(array)
            call generator%uniform(array(ii), center, width)
        end do
    end subroutine uniform_real_array

    subroutine uniform_complex_scalar(generator, scalar, center, width)
        class(random_number_generator), intent(inout) :: generator
        complex(q), intent(out) :: scalar
        real(q), intent(in), optional :: center, width
        real(q) re_part, im_part
        call generator%uniform(re_part, center, width)
        call generator%uniform(im_part, width=width) ! only real part shifted
        scalar = cmplx(re_part, im_part, kind=q)
    end subroutine uniform_complex_scalar

    subroutine uniform_complex_array(generator, array, center, width)
        class(random_number_generator), intent(inout) :: generator
        complex(q), intent(out) :: array(:)
        real(q), intent(in), optional :: center, width
        integer ii
        do ii = 1, size(array)
            call generator%uniform(array(ii), center, width)
        end do
    end subroutine uniform_complex_array


    type(pcg_generator) function new_pcg_generator(seed, stream) result (generator)
        integer(qi8), intent(in) :: seed, stream
        integer(qi4) ignored
        generator%increment = ior(ishft(stream, 1), 1_qi8)
        call random_number(generator, ignored)
        generator%state = generator%state + seed
        call random_number(generator, ignored)
    end function new_pcg_generator

    subroutine pcg_integer_scalar(generator, scalar)
        type(pcg_generator), intent(inout) :: generator
        integer(qi4), intent(out) :: scalar
        integer(qi8) oldstate
        integer(qi4) xorshifted, rotation
        oldstate = generator%state
        generator%state = oldstate * pcg_multiplier + generator%increment
        xorshifted = ishft(ieor(ishft(oldstate, -18), oldstate), -27)
        rotation = -ishft(oldstate, -59)
        ! The line below is equivalent to
        !    scalar = ishftc(xorshifted, rotation)
        ! but the intel compiler does not seem to recognize the cyclic bit shift
        ! unless we use the more complicated
        scalar = ior(ishft(xorshifted, rotation), ishft(xorshifted, iand(rotation, 31_qi4)))
    end subroutine pcg_integer_scalar

    subroutine pcg_integer_array(generator, array)
        type(pcg_generator), intent(inout) :: generator
        integer(qi4), intent(out) :: array(:)
        integer ii
        do ii = 1, size(array)
            call random_number(generator, array(ii))
        end do
    end subroutine pcg_integer_array

    subroutine pcg_real_scalar(generator, scalar)
        type(pcg_generator), intent(inout) :: generator
        real(q), intent(out) :: scalar
        integer(qi4) sample
        call random_number(generator, sample)
        scalar = sample * resolution + 0.5_q
    end subroutine pcg_real_scalar

    subroutine pcg_real_array(generator, array)
        type(pcg_generator), intent(inout) :: generator
        real(q), intent(out) :: array(:)
        integer ii
        do ii = 1, size(array)
            call random_number(generator, array(ii))
        end do
    end subroutine pcg_real_array

    subroutine pcg_complex_scalar(generator, scalar)
        type(pcg_generator), intent(inout) :: generator
        complex(q), intent(out) :: scalar
        real(q) parts(2)
        call random_number(generator, parts)
        scalar = cmplx(parts(1), parts(2), kind=q)
    end subroutine pcg_complex_scalar

    subroutine pcg_complex_array(generator, array)
        type(pcg_generator), intent(inout) :: generator
        complex(q), intent(out) :: array(:)
        integer ii
        do ii = 1, size(array)
            call random_number(generator, array(ii))
        end do
    end subroutine pcg_complex_array

end module random_seeded
