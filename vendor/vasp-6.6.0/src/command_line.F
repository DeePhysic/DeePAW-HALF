#include "symbol.inc"
!> Read arguments from the command line
!>
!> The point of command line arguments is that they facilitate getting some quicker
!> answers then reading the INCAR file would, because it is processed before the
!> INCAR file.
!>
!> Currently the following command line arguments are implemented
!>
!> # --version | -v
!>
!> Return some information about the version of VASP that is used.
!>
!> # --dry-run | -n
!>
!> halt program and wait to attach debugger 
!>
!> # --debug | -g
!>
!> Execute only a part of the code to test whether the setup would fit into
!> memory and the input can be processed.
module command_line

    use base, only: in_struct
    use tutor, only: vtutor, isAlert, isAdvice, isError
    use build_info, only : cpp_options, link_line

    implicit none

    integer, parameter :: default_status = -1
    character(len=*), parameter :: dry_run = "VASP will execute in dry-run mode. &
        &This means that only a small subset of the code will be run through to &
        &get some insight whether the setup is reasonable. The input files will &
        &be read and some sanity checks are conducted on them. After addressing &
        &any potential issues rerun VASP removing the dry-run command line option."
#ifdef GDBDEBUG
    character(len=*), parameter :: debug = "VASP will write all running PIDs to &
        &PAUSE_VASP and pause execution, i.e. waiting for debugger to attach to it. &
        &Once attached, continue execution from debugger. For optimal usage, &
        &compile VASP with debugging symbols and backtrace information and use &
        &the tmpi-debug tool."
#endif
    type parsed_argument
        character(len=:), allocatable :: description
        logical :: should_stop = .false.
        integer :: stop_status = 0
        integer :: status_ = default_status
    end type parsed_argument

contains

    !> Parse any potential command line arguments and react to them.
    subroutine parse_command_line(io, should_write)
        !> Some command line arguments modify the settings.
        type(in_struct), intent(inout) :: io
        !> Switch writing to the standard output on or off. Useful to limit output to a single rank.
        logical, intent(in) :: should_write
        character(len=256) text
        type(parsed_argument) argument
        integer index_
        do index_ = 1, command_argument_count()
            call get_command_argument(index_, text)
            call parse_argument(trim(text), io, argument)
            if (should_write) call log_description(argument%description, argument%status_)
            if (argument%should_stop) call vtutor%stopCode(argument%stop_status)
        end do
    end subroutine parse_command_line

    pure subroutine parse_argument(text, io, argument)
        use version
        character(len=*), intent(in) :: text
        type(in_struct), intent(inout) :: io
        type(parsed_argument), intent(out) :: argument
        select case(text)
        case ('--version', '-v')
            argument%description = vasp()
            argument%should_stop = .true.
        case ('--cpp-options', '-c')
            argument%description = cpp_options
            argument%should_stop = .true.
        case ('--link-line', '-l')
            argument%description = link_line
            argument%should_stop = .true.
        case ('--dry-run', '-n')
            argument%description = dry_run
            argument%status_ = isAdvice
            io%dry_run = .true.
#ifdef GDBDEBUG
        case ('--debug', '-g')
            argument%description = debug
            argument%status_ = isAdvice
            io%halt = .true.
#endif
        case default
            argument%description = "Command line argument '" // trim(text) // "' was not understood."
            argument%should_stop = .true.
            argument%stop_status = isError
        end select
    end subroutine parse_argument

    subroutine log_description(description, status_)
        use iso_fortran_env, only: output_unit
        character(len=*), intent(in) :: description
        integer, intent(in) :: status_
        select case (status_)
        case (default_status)
            write(output_unit, '(a)') description
        case (isAlert)
            call vtutor%alert(description)
        case (isAdvice)
            call vtutor%advice(description)
        case (isError)
            call vtutor%error(description)
        end select
    end subroutine log_description

end module command_line
