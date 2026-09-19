#include "symbol.inc"
!> Error handling module for VASP
!>
!> # Overview
!>
!> All noticable events in Vasp should be dealt with by this module, which defines
!> the [different categories](#Error severities) of events that the code recognizes. If
!> an event is severe enough the code will be stopped, otherwise just the event message
!> will be written. The user can adjust which severity of errors they want to stop the code
!> on. Funneling all error handling through this single interface allows a consistent cleanup
!> behavior and error code returned on exit.
!>
!> In addition this module defines a list of recognized error flags, which are
!> associated with a predefined error message. The severity of these messages is still
!> given by the calling code.
!>
!> # Usage
!>
!> To use the error handling you need to import the global logger from the module
!> `use tutor, only: vtutor`. When an event occured that you want to notify the user about,
!> you should call the [appropriate event subroutine](#Error severities). Note that for most
!> cases only the root node will produce any output so make sure that you call the error
!> processing routines from all nodes. Generally, you want to put the error message inside
!> the code with the following exception: If a solution to the issue exists that is difficult
!> to implement but possible to fix for the user, you may introduce an elaborate message
!> explaining how to overcome the issue in this module by defining a new error flag.
!> Then you use `call vtutor%write` with the appropriate severity and the newly defined error
!> flag.
!>
!> # Error severities
!>
!> ## Bug
!>
!> You should use the *bug* severity to check whether something unexpected happened. It
!> should not be possible to reach this statement by any setting in the Incar file. So
!> it does not indicate a usage error of the code, but that something was implemented
!> in a wrong way. An example would be that your implementation assumes all values of
!> an array are positive, then you can `call vtutor%bug` if the input contains any
!> negative values. Note that this routine behaves differently in two significant ways.
!> It has a different interface taking also the filename and line number as an argument
!> and it will always write even if you are not on the ionode.
!>
!> ## Error
!>
!> The *error* severity indicates an unrecoverable usage error, typically caused by an
!> incorrectly formed input file. Ideally your error message should state, why the error
!> occured but also give some hints into how to fix it.
!>
!> ## Alert
!>
!> You can raise an *alert* message if something strange happened. The idea behind this
!> message is that it is possibly a usage error, but there is a relatively straight
!> forward way to fix it. In this case, you may choose to perform this fix, but you should
!> raise an error message of this severity to notify the user that you overwrote their
!> input. If the user specifies the conservative mode `STOP_ON = Alert`, the code will
!> stop if it encounters this message.
!>
!> ## Warning
!>
!> You can notify the user that something unexpected happened by a *warning*. Please use
!> this kind of message sparingly, because most users will just ignore it. So carefully
!> assess whether it is something you can fix by a more robust coding or warrants a more
!> severe error message.
!>
!> ## Advice
!>
!> You can *advice* the user to change their setup to improve their performance. This
!> kind of messages should not change the result of the calculation for typical requirements.
!> Another use case are to provide literature references for code contributed by external
!> developers.
module ml_ff_tutor
    use ML_FF_PREC, only: q
    use iso_fortran_env, only: output_unit, error_unit

    implicit none

    enum, bind(c)
        enumerator :: never, isBug, isError, isAlert, isWarning, isAdvice
    end enum
!$acc declare copyin(isBug, isError)

    enum, bind(c)
        ! NOTE: add new error codes to the end to avoid screwing up the tests
        enumerator :: UndefinedCode
        enumerator :: STORAGE_LRC
        enumerator :: numTopic
    end enum

    type logger
        integer :: unitOut = output_unit, unitErr = error_unit, unitBug = error_unit
        integer :: stopOn = isError
        logical :: discardOutputToNegativeUnit = .false.
    contains
        procedure :: bug => logBug
        procedure :: error => logError
        procedure :: alert => logAlert
        procedure :: warning => logWarning
        procedure :: advice => logAdvice
        procedure shouldStop
        procedure stopCode
        procedure write
        procedure printMessageToUnits
    end type logger

    ! Note that the global logger vtutor needs to be initialized before its first use
    ! otherwise it will use the default output and error unit
    type(logger) :: ml_tutor

    character(*), parameter :: headerBug = '&
        &                _     ____    _    _    _____     _ \n&
        &               | |   |  _ \  | |  | |  / ____|   | |\n&
        &               | |   | |_) | | |  | | | |  __    | |\n&
        &               |_|   |  _ <  | |  | | | | |_ |   |_|\n&
        &                _    | |_) | | |__| | | |__| |    _ \n&
        &               (_)   |____/   \____/   \_____|   (_)\n\n'
    character(*), parameter :: footerBug = '\n\n&
        &If you are not a developer, you should not encounter this problem.\n&
        &Please submit a bug report.\n\n'

    character(*), parameter :: headerError = '\n&
        &EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###\n&
        &E        R     R  R     R  O     O  R     R     ###     ###     ###\n&
        &E        R     R  R     R  O     O  R     R     ###     ###     ###\n&
        &EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       # \n&
        &E        R   R    R   R    O     O  R   R                          \n&
        &E        R    R   R    R   O     O  R    R      ###     ###     ###\n&
        &EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###\n\n'
    character(*), parameter :: footerError = '\n\n&
        &  ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! <----\n\n'

    character(*), parameter :: headerAlert = '\n&
        &      W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!\n&
        &      W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!\n&
        &      W    W  A    A  R    R  N N  N  II  N N  N  G       !!!\n&
        &      W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   ! \n&
        &      WW  WW  A    A  R   R   N   NN  II  N   NN  G    G     \n&
        &      W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!\n\n'
    character(*), parameter :: footerAlert = '\n\n'

    character(*), parameter :: headerWarning = ' WARNING: '
    character(*), parameter :: footerWarning = ''

    character(*), parameter :: headerAdvice = '\n&
        &          ----> ADVICE to this user running VASP <----\n\n'
    character(*), parameter :: footerAdvice = '\n\n'


    character(*), parameter :: horizontalLineSegment = '-', verticalLineSegment = '|'

    type messageFormatter
        integer :: tabWidth = 0, padding = 5, totalWidth = 77
        logical :: useBox = .true.
    contains
        procedure formatMessage
        procedure horizontalLineOfBox
        procedure wrapMessage
        procedure determineWidth
        procedure formatLine
    end type messageFormatter

    type argument
        integer, allocatable :: ival(:)
        logical, allocatable :: lval(:)
        real(q), allocatable :: rval(:)
    end type argument

    type messageSplitter
        integer :: first = 1, last = 0
        character(len=:), allocatable :: message
    contains
        procedure finished
        procedure nextPart
    end type messageSplitter

    private
    public logger, argument, ml_tutor, never, isBug, isError, isAlert, isWarning, isAdvice, &
        numTopic, STORAGE_LRC

#ifdef MPI
    include 'mpif.h'
#endif

contains

    !> Print a bug error message.
    !>
    !> Calling this routine will stop the code with an error message. The intent of
    !> this routine is for coding errors. So use when a certain assertion about the
    !> passed quantities is not fulfilled. It should not be possible to raise these
    !> error by the user and indicate a coding error, which is why this routine is
    !> special in that it takes the filename and linenumber to help with the
    !> debugging. You can use the precompiler flags `__FILE__` and `__LINE__` for that.
    subroutine logBug(this, message, filename, linenumber)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message  !< error message printed to the output
        character(len=*), intent(in) :: filename  !< file in which the error occured
        integer, intent(in) :: linenumber  !< line in which the error occured
        type(messageFormatter) formatter
        character(len=:), allocatable :: whereInCode, formattedMessage
        whereInCode = createWhereInCodeInfo(filename, linenumber)
        formattedMessage = formatter%formatMessage(headerBug // whereInCode // message // footerBug)
        call this%printMessageToUnits(formattedMessage, [this%unitBug])
        if (this%stopOn >= isBug) call this%stopCode(isBug)
    end subroutine logBug

    function createWhereInCodeInfo(filename, linenumber) result(res)
        use string, only: str
        character(len=*), intent(in) :: filename
        integer, intent(in) :: linenumber
        character(len=:), allocatable :: res
        if (linenumber > 0) then
           res = "internal error in: " // trim(filename) // "  at line: " // str(linenumber) // "\n\n"
        else
           res = "internal error in: " // trim(filename) // "  at line: " // "UnknownLine" // "\n\n"
        endif
    end function createWhereInCodeInfo

    !> Print a user error message.
    !>
    !> Calling this routine will stop the code with the given error message. Use this
    !> e.g. when some settings are incompatible with a particular branch of the code.
    !> If possible give advice to the user what to change so that s*he can run the
    !> calculation.
    subroutine logError(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message  !< error message printed to the output
        logical, intent(in), optional :: dont_stop  !< use this to enforce the code will not stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formattedMessage = formatter%formatMessage(headerError // message // footerError)
        call this%printMessageToUnits(formattedMessage, [this%unitOut, this%unitErr])
        if (this%shouldStop(isError, dont_stop)) call this%stopCode(isError)
    end subroutine logError

    !> Print a severe warning.
    !>
    !> Calling this routine will print an error message, but not stop the code unless
    !> the user selects `STOP_ON = Alert` in the INCAR file. Use this when something goes
    !> wrong, but the code may still recover.
    subroutine logAlert(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message  !< warning message printed to the output
        logical, intent(in), optional :: dont_stop  !< use this to enforce the code will not stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formattedMessage = formatter%formatMessage(headerAlert // message // footerAlert)
        call this%printMessageToUnits(formattedMessage, [this%unitOut, this%unitErr])
        if (this%shouldStop(isAlert, dont_stop)) call this%stopCode(isAlert)
    end subroutine logAlert

    !> Print a warning message to the OUTCAR file.
    !>
    !> Because the warning message is only present in the OUTCAR file it will not rise
    !> to the users attention typically. Avoid this if possible.
    subroutine logWarning(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message  !< warning message printed to the OUTCAR file
        logical, intent(in), optional :: dont_stop  !< use this to enforce the code will not stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formatter = messageFormatter(tabWidth = 5, padding = 0, totalWidth = 78, &
            useBox = .false.)
        formattedMessage = formatter%formatMessage(headerWarning // message // footerWarning)
        call this%printMessageToUnits(formattedMessage, [this%unitOut])
        if (this%shouldStop(isWarning, dont_stop)) call this%stopCode(isWarning)
    end subroutine logWarning

    !> Print advice how to improve the calculation.
    !>
    !> This message is very prominently featured in the standard output, but does not
    !> terminate the code, even if `STOP_ON = Alert` is used. Use this to warn about
    !> issues that can impact the performance but should not alter the final result.
    subroutine logAdvice(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message  !< advice message printed to the output
        logical, intent(in), optional :: dont_stop  !< use this to enforce the code will not stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formattedMessage = formatter%formatMessage(headerAdvice // message // footerAdvice)
        call this%printMessageToUnits(formattedMessage, [this%unitOut, this%unitErr])
        if (this%shouldStop(isAdvice, dont_stop)) call this%stopCode(isAdvice)
    end subroutine logAdvice

    function formatMessage(formatter, message) result(res)
        use string, only: parseNewLine
        class(messageFormatter), intent(in) :: formatter
        character(len=*), intent(in) :: message
        character(len=:), allocatable :: res
        res = formatter%horizontalLineOfBox()
        res = res // formatter%wrapMessage(parseNewline(message))
        res = res // formatter%horizontalLineOfBox()
    end function formatMessage

    function horizontalLineOfBox(formatter) result (res)
        class(messageFormatter), intent(in) :: formatter
        character(len=:), allocatable :: res
        if (formatter%useBox) then
            res = ' ' // repeat(horizontalLineSegment, formatter%totalWidth) // new_line('n')
        else
            res = ''
        end if
    end function horizontalLineOfBox

    function wrapMessage(formatter, message) result(res)
        class(messageFormatter), intent(in) :: formatter
        character(len=*), intent(in) :: message
        character(len=:), allocatable :: res, part
        type(messageSplitter) splitter
        logical firstLine
        integer width
        res = ''
        firstLine = .true.
        allocate(character(len=len(message)) :: splitter%message)
        splitter%message = message
        do while (.not.splitter%finished())
            width = formatter%determineWidth(firstLine)
            part = splitter%nextPart(width)
            res = res // formatter%formatLine(firstLine, part)
            firstLine = .false.
        end do
    end function wrapMessage

    logical function finished(splitter)
        class(messageSplitter), intent(in) :: splitter
        finished = splitter%first > len(splitter%message)
    end function finished

    function nextPart(splitter, width) result (res)
        class(messageSplitter), intent(inout) :: splitter
        integer, intent(in) :: width
        character(len=:), allocatable :: res
        integer, parameter :: skipSpace = 2
        integer ii
        character char
        do ii = splitter%first, splitter%first + width
            if (ii > len(splitter%message)) then
                splitter%last = len(splitter%message)
                exit
            end if
            char = splitter%message(ii:ii)
            if (char == ' '.or.char == new_line('n')) splitter%last = ii - 1
            if (char == new_line('n')) exit
        end do
        res = splitter%message(splitter%first:splitter%last)
        splitter%first = splitter%last + skipSpace
    end function nextPart

    integer function determineWidth(formatter, firstLine) result(width)
        class(messageFormatter), intent(in) :: formatter
        logical, intent(in) :: firstLine
        width = formatter%totalWidth - 2 * formatter%padding
        if (.not.firstLine) width = width - formatter%tabWidth
    end function determineWidth

    function formatLine(formatter, firstLine, part) result(res)
        class(messageFormatter), intent(in) :: formatter
        logical, intent(in) :: firstLine
        character(len=*), intent(in) :: part
        character(len=:), allocatable :: res
        integer numSpace
        if (formatter%useBox) then
            res = repeat(' ', formatter%padding) // part
            numSpace = formatter%totalWidth - len(res)
            res = verticalLineSegment // res // repeat(' ', numSpace) // verticalLineSegment
        else
            res = part
            if (.not.firstLine) res = repeat(' ', formatter%tabWidth) // res
        end if
        res = res // new_line('n')
    end function formatLine

    subroutine printMessageToUnits(this, message, units)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message
        integer, intent(in) :: units(:)
        integer iunit
        do iunit = 1, size(units)
            if (this%discardOutputToNegativeUnit .and. units(iunit) < 0) cycle
            if (unitNotOpen(units(iunit))) cycle
            !$omp critical(omp_wrt_stdout)
            write(units(iunit), '(a)') message
            flush(units(iunit))
            !$omp end critical(omp_wrt_stdout)
        end do
    end subroutine printMessageToUnits

    logical function unitNotOpen(iunit)
        integer, intent(in) :: iunit
        logical unitOpen
        inquire(iunit, opened = unitOpen)
        unitNotOpen = .not.unitOpen
    end function unitNotOpen
#ifdef MPI
    subroutine stopMpi(errCode)
        integer, intent(in), optional :: errCode
        integer errCode_, request, ierr
        logical complete
        errCode_ = 0
        if (present(errCode)) errCode_ = errCode
        if (errCode_ /= 0) then
            call MPI_ibarrier(MPI_comm_world, request, ierr)
            complete = tryWaitingForOtherCPU(request, ierr)
            if (ierr /= MPI_success.or..not.complete) &
                call MPI_abort(MPI_comm_world, 1, ierr)
        end if
        call MPI_finalize(ierr)
        if (ierr /= MPI_success) &
            call MPI_abort(MPI_comm_world, 1, ierr)
    end subroutine stopMpi

    logical function tryWaitingForOtherCPU(request, ierr) result(complete)
        integer, intent(in) :: request
        integer, intent(inout) :: ierr
        integer try
        if (ierr /= MPI_success) return
        do try = 1, 30
            call MPI_test(request, complete, MPI_STATUS_IGNORE, ierr)
            if (ierr /= MPI_success.or.complete) exit
            call sleep(1)
        end do
    end function tryWaitingForOtherCPU
#endif

    logical function shouldStop(this, severity, dont_stop)
        class(logger), intent(in) :: this
        integer, intent(in) :: severity
        logical, intent(in), optional :: dont_stop
        shouldStop = this%stopOn >= severity
        if (present(dont_stop)) then
            shouldStop = .not.dont_stop.and.shouldStop
        end if
    end function shouldStop

    !> Stop the code in a clean manner finalizing the necessary accesses.
    subroutine stopCode(this, errCode)
        class(logger), intent(in) :: this
        integer, intent(in), optional :: errCode  !< when an error code is provided the code will stop in error
        logical errStop
        errStop = .false.
        if (present(errCode)) errStop = errCode /= never
        CALLMPI( stopMpi(errCode) )
        if (errStop) stop 1
        stop
    end subroutine stopCode

    !> Write error message predefined in this module
    !>
    !> Use this for messages that may be called from multiple routines or are too
    !> length so that it would impact the readability of the other code. Introduce
    !> a new error message to the generateMessage routine below.
    subroutine write(this, severity, topic, arg, filename, linenumber)
        class(logger), intent(in) :: this
        integer, intent(in) :: severity  !< the level of severity of the error
        integer, intent(in) :: topic  !< specify what kind of error occured
        type(argument), intent(in), optional :: arg  !< can be used to provide arguments to the error message
        character(len=*), intent(in), optional :: filename  !< for bug messages only: file in which the error occured
        integer, intent(in), optional :: linenumber  !< for bug messages only: line in which the error occured
        character(len=:), allocatable :: message, fnam
        integer line
        message = generateMessage(topic, arg)
        select case (severity)
        case (isBug)
            fnam = "UnkownFile" ; if (present(filename)) fnam = filename
            line = -1 ; if (present(linenumber)) line = linenumber
            call this%bug(message, fnam, line)
        case (isError)
            call this%error(message)
        case (isAlert)
            call this%alert(message)
        case (isWarning)
            call this%warning(message)
        case (isAdvice)
            call this%advice(message)
        end select
    end subroutine write

    function generateMessage(topic, arg) result(res)
        use string, only: str
        integer, intent(in) :: topic
        type(argument), intent(in), optional :: arg
        character(len=:), allocatable :: res, work
        select case (topic)
        case(STORAGE_LRC)
            res = "Not enough storage reserved for local &
               &reference configurations, please increase ML_MB. If you &
               &intend to keep the current storage size you may use &
               &ML_LBASIS_DISCARD=.TRUE. to enable automatic discarding. &
               &Alternatively, copy ML_ABN to ML_AB and continue with a 10 times &
               &increased ML_EPS_LOW (however, keep EPS_LOW<1E-7). This may &
               &yield a more memory-efficient but potentially less accurate &
               &force field."
        case default
            res = 'undefined tutor message ' // str(topic)
        end select
    end function generateMessage

end module ml_ff_tutor
