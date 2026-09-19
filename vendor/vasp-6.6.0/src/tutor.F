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
module tutor
    use prec, only: q
    use iso_fortran_env, only: output_unit, error_unit

    implicit none

    enum, bind(c)
        enumerator :: never, isBug, isError, isAlert, isWarning, isAdvice
    end enum
!$acc declare copyin(isBug, isError)

    enum, bind(c)
        ! NOTE: add new error codes to the end to avoid screwing up the tests
        enumerator :: UndefinedCode
        enumerator :: RealSpaceWithoutOptimization, SplineInterpolateProjectors, &
            VASP_4_4, NoRealSpaceAndYouCould, NoRealSpaceAndYouShould, IALGO8, &
            RealSpaceNoMoreRecommended, WrongOptimizationRealSpace, LongLattice, &
            DifferentXCGradTypes, NumberOfElectrons, NBANDSchanged, &
            DifferentLDAXCTypes, DifferentXCFiles, AM05spin, NoLDAU, &
            DoubleCounting, NoWannier, LREALA, DifferentSlaterXCTypes, &
            PartialDOS, NoOptics, Position, CenterOfMassDrift, EnforcedLDA, &
            DipolEpsilon, DipolCubic, ALGO_A_ISMEAR, ModelGW, HFUSPP, &
            LHFLinearResponse, LinearReponsePOSCAR, LinearResponseLREAL, &
            FockForce, HF_NPAR, NPAREfficiency, NPAR, LHFOneCenter, HFParam, &
            OEPdouble, NMAXFOCKAE, ISYM2, LSPECTRALinefficient, LSPECTRALGWset, &
            LPOTOK, XiSingular, AEXXzero, GWoptics, GW_HFoptics, ENCUTGW_LHF, &
            DeprENCUTFOCK, PRECFOCK, KPOINTS_HF, KPOINTS_PEAD, KPOINTS_INTER, &
            GWkweight, GWnoreal, ACFDTRnotSupported, PEADnoVirtuals, PEADlargeField, &
            PEAD_NCOREnot1, PEADPolNoMetal, PEADEpsNoMetal, &
            LINTET, GRDSHFT_NOTET, GammaK, KPOINTS_IO, NoInversion, &
            FFTGridIsNotSufficient, MetaGGAandForces, MetaGGAResponse, &
            SmallCSHIFT, TetNoMetal, POTIMlarge, GRDSHFTnotInIBZ, GRDSHFTK, &
            NoOcc, NoOddGrid, GammaC, ISYM_MD, BSEantires, BSEconjg, BSEncv, &
            Wmissing, AuxWAVECAR, SpinSpiral, HighestBandsOcccupied, LVTOT, &
            vdWKlimes, NoKPAR, FERWEOptics, OmegaPar, TauPar, &
            TauParCPU, ImagGridParForLTMP2, chiGG_NOmega, OmegaGridOMEGAMIN, WannNotFound, &
            TargetNotSet, NTARGET_STATES, NCRPALOW, CRPADisentangle, LCAO_INCAR, &
            SHMEM_CRPA, OmegaGridOmegaTL, NOmega32, OmegaGrid140, NUPDOWN, &
            NPAR_IALGO8, SHMEMError, KProjection, Gamma, Gamma2, DifferentRelXCTypes, &
            Alternant, RandomSeed, PEADnowIsMetal, &
            NBANDlessNWAN, ELPHtesting, &
            NOMEGA_DUMPtooLarge, SigmaTMBPT, ALLOC_WAVEFUN, &
            FermiFiniteTemp, FermiIntegration, BlochBasisUsed, TargetOut, NCRPA_BANDSout, &
            LOCPROJnotSet, ALGO_GW, ALGO_GW0, ALGO_G0W0, ALGO_scGW0, &
            ALGO_scGW, POTCARmissing, SymprecDependence, &
            ALLOC_RESPONSE, ALLOC_RESPONSE_DISTR, MAXMEM_LOW, AUTOSET_MEM, &
            RPAFORCES_COMPAT, RPAENCUTGWSOFT, ORTHONORMALIZATION_FAILED, &
            ACFDTR_OpenACC, DISTRIBUTE_RESPONSE_SUPER_MU_NGVECTOR, LOPTICS_SPLINES, NoVGVapplied,&
            SinglesNaturalo,BSEoptics, BSEnocder, ElectronPhononRequiresHDF5, NCOREAndGWR, NonCollinearAndGWR, &
            FockStdGWR, LuttingerWardNaturalo, LOPTICSAndGW, LIBXCWrongName, state_order, &
            CheckpointContinue, CheckpointSingle
        enumerator :: numTopic
    end enum
!$acc declare copyin(DISTRIBUTE_RESPONSE_SUPER_MU_NGVECTOR)

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
    type(logger) :: vtutor

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
        character(20), allocatable :: cval(:)
    end type argument

    type messageSplitter
        integer :: first = 1, last = 0
        character(len=:), allocatable :: message
    contains
        procedure finished
        procedure nextPart
    end type messageSplitter

    private
    public logger, argument, newArgument, vtutor, never, isBug, isError, isAlert, isWarning, isAdvice, &
        numTopic, RealSpaceWithoutOptimization, SplineInterpolateProjectors, VASP_4_4, &
        NoRealSpaceAndYouCould, NoRealSpaceAndYouShould, IALGO8, DifferentXCFiles, &
        RealSpaceNoMoreRecommended, WrongOptimizationRealSpace, LongLattice, &
        DifferentXCGradTypes, NumberOfElectrons, NBANDSchanged, DifferentLDAXCTypes, &
        AM05spin, NoLDAU, DoubleCounting, NoWannier, LREALA, DifferentSlaterXCTypes, &
        PartialDOS, NoOptics, Position, CenterOfMassDrift, EnforcedLDA, &
        DipolEpsilon, DipolCubic, ALGO_A_ISMEAR, ModelGW, HFUSPP, ISYM_MD, &
        LHFLinearResponse, LinearReponsePOSCAR, LinearResponseLREAL, FockForce, &
        HF_NPAR, NPAREfficiency, NPAR, LHFOneCenter, HFParam, OEPdouble, NMAXFOCKAE, &
        ISYM2, LSPECTRALinefficient, LSPECTRALGWset, LPOTOK, XiSingular, AEXXzero, &
        GWoptics, GW_HFoptics, ENCUTGW_LHF, DeprENCUTFOCK, PRECFOCK, KPOINTS_HF, &
        KPOINTS_PEAD, KPOINTS_INTER, GWkweight, GWnoreal, ACFDTRnotSupported, &
        PEADnoVirtuals, PEADlargeField, PEAD_NCOREnot1, PEADPolNoMetal, PEADEpsNoMetal, &
        PEADnowIsMetal, LINTET, GRDSHFT_NOTET, GammaK, KPOINTS_IO, NoInversion, &
        FFTGridIsNotSufficient, MetaGGAandForces, MetaGGAResponse, SmallCSHIFT, &
        TetNoMetal, POTIMlarge, GRDSHFTnotInIBZ, GRDSHFTK, NoOcc, NoOddGrid, &
        GammaC, BSEantires, BSEconjg, BSEncv, Wmissing, AuxWAVECAR, SpinSpiral, &
        HighestBandsOcccupied, LVTOT, vdWKlimes, NoKPAR, FERWEOptics, &
        OmegaPar, TauPar, TauParCPU, ImagGridParForLTMP2, chiGG_NOmega, &
        OmegaGridOMEGAMIN, WannNotFound, TargetNotSet, NTARGET_STATES, NCRPALOW, &
        CRPADisentangle, LCAO_INCAR, SHMEM_CRPA, OmegaGridOmegaTL, NOmega32, &
        OmegaGrid140, NUPDOWN, NPAR_IALGO8, SHMEMError, KProjection, Gamma, Gamma2, &
        DifferentRelXCTypes, Alternant, RandomSeed, &
        NBANDlessNWAN, ELPHtesting, &
        NOMEGA_DUMPtooLarge, SigmaTMBPT, FermiFiniteTemp, FermiIntegration, &
        BlochBasisUsed, TargetOut, NCRPA_BANDSout, LOCPROJnotSet, ALGO_GW, &
        ALGO_GW0, ALGO_G0W0, ALGO_scGW0, ALGO_scGW,POTCARmissing, SymprecDependence, &
        vtutor_error, ALLOC_RESPONSE, ALLOC_RESPONSE_DISTR, MAXMEM_LOW, AUTOSET_MEM, &
        ALLOC_WAVEFUN, RPAFORCES_COMPAT, RPAENCUTGWSOFT, deactivateVtutorOutputForTest, vtutor_bug, &
        ORTHONORMALIZATION_FAILED, ACFDTR_OpenACC, &
        DISTRIBUTE_RESPONSE_SUPER_MU_NGVECTOR, LOPTICS_SPLINES, NoVGVapplied, &
        SinglesNaturalo,BSEoptics, BSEnocder, ElectronPhononRequiresHDF5, NCOREAndGWR, NonCollinearAndGWR, &
        FockStdGWR, LuttingerWardNaturalo, LOPTICSAndGW, LIBXCWrongName, state_order, &
        CheckpointContinue, CheckpointSingle

#ifdef MPI
    include 'mpif.h'
#endif

contains

    !> Use this function where `call vtutor%bug` is not allowed.
    subroutine vtutor_bug(message, filename, linenumber)
        character(len=*), intent(in) :: message, filename
        integer, intent(in) :: linenumber
        call vtutor%bug(message, filename, linenumber)
    end subroutine vtutor_bug

    !> Use this function where `call vtutor%error` is not allowed
    subroutine vtutor_error(message)
        character(len=*), intent(in) :: message
        call vtutor%error(message)
    end subroutine vtutor_error

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
        logical set_last
        integer ii
        character char
        set_last = .false.
        do ii = splitter%first, splitter%first + width
            if (ii > len(splitter%message)) then
                set_last = .true.
                splitter%last = len(splitter%message)
                exit
            end if
            char = splitter%message(ii:ii)
            if (char == ' '.or.char == new_line('n')) then
                set_last = .true.
                splitter%last = ii - 1
            end if
            if (char == new_line('n')) exit
        end do
        if (set_last) then
            res = splitter%message(splitter%first:splitter%last)
            splitter%first = splitter%last + skipSpace
        else
            splitter%last = splitter%first + width
            res = splitter%message(splitter%first:splitter%last)
            splitter%first = splitter%last + 1
        end if
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
        use vhdf5_struct_def
        class(logger), intent(in) :: this
        integer, intent(in), optional :: errCode  !< when an error code is provided the code will stop in error
        logical errStop
        integer :: h5ierr
        errStop = .false.
        if (present(errCode)) errStop = errCode /= never
#ifdef VASP_HDF5
        !close h5
        if( ih5wavefileid>=0 ) then
            call h5fclose_f(ih5wavefileid, h5ierr)
        endif
#endif
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
        case(RealSpaceWithoutOptimization)
            res = 'The real-space projection scheme for the treatment of the non-local &
                &pseudopotentials has been switched on but on file POTCAR I have &
                &not found any entries, signalling me that you have never performed a &
                &real-space optimization of all non-local projectors! BE WARNED that &
                &a calculation using the real-space projection scheme together with &
                &non-local projectors which were not real-space optimized (according &
                &to the plane-wave cutoff used in this calculation) might give more &
                &or less inaccurate results. I hope you know what you are doing.'
        case(SplineInterpolateProjectors)
            res = 'The projection operators are constructed using a spline interpolation. &
                &This is needed for the calculation of chemical shifts by means of &
                &linear response but yields slightly different total energies compared &
                &to the default setup of the projectors (NLSPLINE=.FALSE.). '
        case(VASP_4_4)
            res = 'You are running vasp.4.5 in the vasp.4.4 compatibility mode. &
                &vasp.4.5 has some numerical improvements, which are not applied in this &
                &mode (charges at unbalanced lattice vectors are no longer zeroed, &
                &PAW augmentation charges are integrated more accurately).'
        case(NoRealSpaceAndYouCould)
            res = 'You have a (more or less) ''large supercell'' and for larger cells &
                &it might be more efficient to use real-space projection operators. &
                &Therefore, try LREAL= Auto in the INCAR file.\n&
                &Mind: If you want to do very accurate calculations, keep the &
                &reciprocal projection scheme (i.e. LREAL=.FALSE.).'
        case(NoRealSpaceAndYouShould)
            res = 'You have a (more or less) ''large supercell'' and for larger cells &
                &it might be more efficient to use real-space projection operators. &
                &Therefore, try LREAL= Auto in the INCAR file.\n&
                &Mind: For very accurate calculation, you might also keep the &
                &reciprocal projection scheme (i.e. LREAL=.FALSE.).'
        case(IALGO8)
            res = 'Recently, Corning got a patent for the Teter Allan Payne algorithm. &
                &Therefore, VASP.4.5 does not support IALGO=8 any longer. &
                &A much faster algorithm, IALGO=38, is now implemented in VASP. &
                &This algorithm is a blocked Davidson-like method and as reliable as &
                &IALGO=8 used to be.\n&
                &For ultimate performance, IALGO=48 is still the method of choice &
                &-- SO MUCH ABOUT PATENTS :)'
        case(RealSpaceNoMoreRecommended)
            res = 'You have a (more or less) ''small supercell'' and for smaller cells &
                &it is recommended to use the reciprocal-space projection scheme! &
                &The real-space optimization is not efficient for small cells and it &
                &is also less accurate ...\n&
                &Therefore, set LREAL=.FALSE. in the INCAR file.'
        case(WrongOptimizationRealSpace)
            res = 'One real-space projector is optimized for E = ' // str(arg%rval(1)) // &
                ' eV but you are using a cutoff of ENMAX = ' // str(arg%rval(2)) // &
                ' eV (QCUT = ' // str(arg%rval(3)) // ' a.u.). This makes no sense. &
                &Reoptimize the projector with the a.u. value given above.'
        case(LongLattice)
            res = 'One of the lattice vectors is very long (>50 A), but AMIN is rather large. &
                &This can spoil convergence since charge sloshing might occur &
                &along the long lattice vector. If problems with convergence are &
                &observed, try to decrease AMIN to a smaller value (e.g. 0.01).\n&
                &Note: This warning only applies if the self-consistency cycle is used.'
        case(DifferentXCGradTypes, DifferentLDAXCTypes)
            res = 'You have built up your multi-ion-type POTCAR file out of POTCAR &
                &files with incompatible specifications for the XC-types used to &
                &generate the pseudopotential. This makes no sense at all!! What &
                &I found is LEXCH = ' // TRIM(ADJUSTL(arg%cval(1))) // ' for atom types < ' // &
                &TRIM(ADJUSTL(arg%cval(3))) // ' but LEXCH = ' // TRIM(ADJUSTL(arg%cval(2))) // ' was &
                &found for atom type = ' // TRIM(ADJUSTL(arg%cval(3))) // '. Use identical &
                &XC-functionals for the pseudopotential generation for all atom types, please ... !'
        case(NumberOfElectrons)
            res = 'The number of bands is not sufficient to hold all electrons. &
                &I am sorry, but this calculation does not make sense. &
                &If you want to calculate only selected states, you need to set &
                &either EFERMI or EREF. EREF specifies around which energy states &
                &are supposed to be calculated.\nI found '// TRIM(ADJUSTL(arg%cval(1))) // ' = ' // &
                str(arg%ival(2)) // ', NELECT = ' // str(arg%ival(1)) // '.'
        case(NBANDSchanged)
            res = 'The number of bands has been changed from the values supplied &
                &in the INCAR file. This is a result of running the parallel version. &
                &The orbitals not found in the WAVECAR file will be initialized with &
                &random numbers, which is usually adequate. For correlated calculations, &
                &however, you should redo the groundstate calculation.\nI found NBANDS = ' // &
                &str(arg%ival(1)) // '. Now, NBANDS = ' // str(arg%ival(2)) // '.'
        case(DifferentXCFiles)
            res = 'The XC type on the POTCAR file is not compatible with the XC type on the EXHCAR file.'
        case(AM05spin)
            res = 'The AM05 functional for spin-polarized calculations is not supported. &
                &Maybe in a future release...'
        case(NoLDAU)
            res = 'VASP supports LDA+U only for PAW potentials and not for US or NC &
                &pseudopotentials. Please restart with appropriate potentials.'
        case(DoubleCounting)
            res = 'LDA+U in combination with GW calculations causes a double &
                &counting of unidentifiable many-body perturbation terms.'
        case(NoWannier)
            res = 'VASP supports the calculation of maximally-localized Wannier &
                &functions only in the Gamma-only version of the code.'
        case(LREALA)
            res = 'LREAL=A is not well tested yet, please use LREAL=O &
                &or test your results very carefully.'
        case(DifferentSlaterXCTypes)
            res = 'You have built up your multi-ion-type POTCAR file out of POTCAR &
                &files with incompatible specifications for the XC-types used to &
                &generate the pseudopotential. This makes no sense at all!! What &
                &I found is slater parameter = ' // str(arg%rval(1)) // ' for atom &
                &types <= ' // str(arg%ival(1)-1) // ' but slater parameter = ' // &
                str(arg%rval(2)) // ' was found for atom type = ' // str(arg%ival(1)) // &
                '. Use identical slater parameters in the exchange functionals for &
                &the pseudopotential generation for all atom types, please ... !'
        case(DifferentRelXCTypes)
            res = 'You have built up your multi-ion-type POTCAR file out of POTCAR &
                &files with incompatible specifications for the XC-types used to &
                &generate the pseudopotential. This makes no sense at all!! What &
                &I found is that the flag which switches on/off the relativistic &
                &corrections has been set to ' // str(arg%lval(1), compact = .false.) // &
                ' for atom types <= ' // str(arg%ival(1) - 1) // ' but it was &
                &set to ' // str(arg%lval(2), compact = .false.) // ' for atom type no. ' // &
                str(arg%ival(1)) // '. Use identical XC-functionals for ' // &
                'the pseudopotential generation for all atom types, please ...!'
        case(PartialDOS)
            res = 'The partial DOS and the PROCAR file are not evaluated for NPAR/=1. &
                &Please rerun with NPAR=1.'
        case(NoOptics)
            res = 'The file OPTICS can be written only for NPAR=1. &
                &Optical properties can be found in vasprun.xml and OUTCAR, however. &
                &If you need the OPTICS file, rerun with NPAR=1.'
        case(Position)
            res = 'The distance between some ions is very small. &
                &Please check the nearest-neighbor list in the OUTCAR file.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING! '
        case(CenterOfMassDrift)
            res = 'The initial velocities result in a center-of-mass drift but &
                &there must be no drift. The drift will be removed! '
        case(EnforcedLDA)
            res = 'You enforced a specific xc type in the INCAR file but &
                &a different type was found in the POTCAR file.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(DipolEpsilon)
            res = 'You specified a dielectric constant different from 1 (EPSILON). &
                &This indicates that you perform calculations for the bulk. &
                &At the same time, LDIPOL is set to .TRUE., a mode which does not &
                &work properly for bulk systems!\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(DipolCubic)
            res = 'LDIPOL = .TRUE. must be selected only for cubic supercells, since &
                &the quadrupole corrections are currently only implemented for this &
                &specific geometry.\n&
                &Please change your box, or set LDIPOL=.FALSE.'
        case(ALGO_A_ISMEAR)
            res = 'ALGO=A and IALGO=5X tend to fail with the tetrahedron method &
                &(e.g. Bloechl''s method, ISMEAR=-5, is not variational). &
                &Please switch to IMSEAR=0-n, except for DOS calculations. &
                &For DOS calculations use IALGO=53 after pre-converging with ISMEAR>=0.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(ModelGW)
            res = 'Model GW is not a variational method. &
                &The conjugate-gradient algorithm therefore does not work. &
                &Please switch to IALGO=53-55.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(HFUSPP)
            res = 'HF is only implemented for the PAW method and for norm-conserving &
                &pseudopotentials.\n&
                &One of the POTCAR data sets corresponds to an US pseudopotential.'
        case(LHFLinearResponse)
            res = 'HF is not implemented in the linear-response routines. &
                &Please use the finite-difference versions.'
        case(LinearReponsePOSCAR)
            res = 'The POSCAR file is problematic for the linear-response routines &
                &since the displacement is not along Cartesian directions. &
                &Please send this POSCAR to vasp support, we will try to fix the issue. &
                &You can use finite differences for the time being.'
        case(LinearResponseLREAL)
            res = 'Linear response in principle supports LREAL = Auto (or LREAL=.TRUE.). &
                &However, in some cases we have found that the real-space projection &
                &is not sufficiently accurate.\n&
                &I strongly recommend that you switch it off by setting LREAL=.FALSE.'
        case(FockForce)
            res = 'ALLOCATE for the force calculations using Hartree-Fock routines failed. &
                &Your forces are not correct!!!\n&
                &However, the WAVECAR file might be useful for continuation.'
        case(HF_NPAR)
            res = 'HF is only implemented for NPAR = number of nodes. ' // &
                'Otherwise, communication between nodes would be required during the &
                &calculation of the exchange potential, which would be utterly slow.'
        case(NPAREfficiency)
            res = 'For optimal performance we recommend to set \n&
                &  NCORE = 2 up to number-of-cores-per-socket\n' // &
                'NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE). &
                &This setting can greatly improve the performance of VASP for DFT. &
                &The default, NCORE=1 might be grossly inefficient &
                &on modern multi-core architectures or massively parallel machines. &
                &Do your own testing! More info at https://www.vasp.at/wiki/index.php/NCORE\n&
                &Unfortunately you need to use the default for GW and RPA calculations &
                &(for HF NCORE is supported but not extensively tested yet).'
        case(NPAR)
            res = 'VASP internal routines have requested a change of the k-point set. ' // &
                'Unfortunately, this is only possible if NPAR=number of nodes. &
                &Please remove the tag NPAR from the INCAR file and restart the &
                &calculation.'
        case(LHFOneCenter)
            res = 'HF one-center treatment and a full HF treatment are not compatible. &
                &Please either set LHFCALC to .FALSE. or remove the LHFONE tag &
                &from the INCAR file.'
        case(HFParam)
            res = 'A combination of HF parameters was found in the INCAR file that is &
                &not supported, e.g. HFRCUT and HFKIDENT.'
        case(LIBXCWrongName)
            res = 'One of the selected functionals for Libxc does not exist. &
                &Either you mistyped the name or your installed version of &
                &Libxc is not recent enough to include this functional.'
        case(state_order)
          res = 'You selected a constrained occupation calculation &
                &(ISMEAR=-2), but band ordering could not be preserved.&
                & Switch to ALGO=Damped or ALGO=All and set LDIAG=.FALSE.'
        case(CheckpointContinue)
            res = 'You are resuming a calculation from a checkpoint file. &
                &Make sure that you use the same calculational parameters (ENCUT, POTIM, &
                &KPOINTS, etc.) as in the original checkpointed calculation. &
                &Otherwise, the results might be incorrect.'
        case(CheckpointSingle)
            res = 'You are running a single independent calculation in a split &
                &calculation scheme. Make sure that you use the same calculational &
                &parameters (ENCUT, POTIM, KPOINTS, etc.) as in the preparation run. &
                &Otherwise, the results might be incorrect.'
        case(OEPdouble)
            res = 'The EXX-OEP and LHF methods do not support a different grid for HF &
                &and the local Hamiltonian.\nPRECFOCK was set to Normal.'
        case(NMAXFOCKAE)
            res = 'You use NMAXFOCKAE=2 in the INCAR file. Be careful:\n&
                &Adding two functions to mimic the all-electron charge density on the &
                &plane-wave grid requires very fine FFT grids. Please set NGX manually &
                &or reconsider your choice.'
        case(ISYM2)
            res = 'You have selected ISYM=2 for a HF type calculation. This will &
                &symmetrize the charge density but not the exchange potential. &
                &I suggest to use ISYM=3 instead. This uses symmetry to obtain ' // &
                'the orbitals at all k-points in the Brillouin zone, but does not &
                &apply symmetry to the Hartree potential directly. &
                &The resultant charge might have lower symmetry than the crystal &
                &but, at least, Hartree and exchange are fully compatible.'
        case(ISYM_MD)
            res = 'ISYM>0 for a molecular dynamics simmulation is not recommended.'
        case(LSPECTRALinefficient)
            res = 'LSPECTRAL is very inefficient if not sufficient memory is available. &
                &The available memory per core should be set by MAXMEM = XXXX. &
                &The default for MAXMEM is presently XXXX=2800 (MByte). ' // &
                'If this has already been done and if the warning prevails, &
                &you should either decrease NOMEGA or decrease ENCUTGW, or distribute ' // &
                'the calculations onto more nodes.\n&
                &Note: The present calculations might lead to excessive memory paging!!'
        case(LSPECTRALGWset)
            res = 'The default for LSPECTRALGW is now LSPECTRALGW=LSPECTRAL. &
                &This will change QP energies but increases the numerical precision &
                &of the QP slightly (for instance, Z factors are usually more accurate). &
                &You can use the old version by specifying LSPECTRALGW=.FALSE.'
        case(LPOTOK)
            res = 'Presently the RPA-OEP routine requires reading in a POT file &
                &containing an initial guess for the local potential. This could be &
                &generated by using the KLI or EXX-OEP routines (EXXOEP=1 or 2) and &
                &setting LVTOT=.TRUE. in order to write the POT file. &
                &Supply this file by setting ICHARG=4.'
        case(XiSingular)
            res = 'The response function can not be inverted. &
                &This indicates that the number of linearly independent transitions &
                &is smaller than the rank of the response function. &
                &Try to increase the number of k-points, or try to decrease ENCUTGW.'
        case(AEXXzero)
            res = 'You calculate the electron-hole interaction using AEXX=0. &
                &This is hardly what you want to do.\n&
                &Maybe you have forgotten to set AEXX in the INCAR file &
                &(For LHFCALC=.TRUE., the default is AEXX=0.25 but if LHFCALC is not &
                &set, the default is AEXX=0.0).'
        case(BSEoptics)
            res = 'The number of bands in WAVEDER ('// str(arg%ival(1)) //') &
                &is different from NBANDS ('// str(arg%ival(2)) //').&
                & Reading '// str(arg%ival(3)) // ' bands from WAVEDER.'
        case(BSEnocder)
            res = 'WAVEDER is missing or contains insufficient number of bands &
                &for this calculation, which will result in incorrect transitions and &
                &dielectric function. You should redo the groundstate &
                &calculation with LOPTICS = .TRUE. If you wish to proceed anyway, &
                &you can set IDIOT = 0.'
        case(GWoptics)
            res = 'The derivative of the wavefunctions with respect to k (WAVEDER) &
                &can not be found. You should redo the groundstate calculation &
                &using LOPTICS=.TRUE. in order to write the WAVEDER dataset. &
                &However for metals, the present setting is ok.'
        case(GW_HFoptics)
            res = 'Presently LREAL is not supported for the calculation of optical &
                &properties and HF type Hamiltonians.\nPlease use LREAL=.FALSE.'
        case(ENCUTGW_LHF)
            res = 'ENCUTGW is set for a local HF calculation. &
                &This might cause very slow convergence with respect to &
                &ENCUTGW, since even the Hartree potential is cut off. &
                &Use with extreme care (or set ENCUTGW=0 in the INCAR file).'
        case(DeprENCUTFOCK)
            res = 'ENCUTFOCK is no longer supported!\n&
                &Use PRECFOCK to select the mode for HF type calculations:\n&
                &PRECFOCK= L low (coarse grid for HF, normal augmentation charge)\n&
                &PRECFOCK= M medium (normal grid for HF, normal augmentation charge)\n&
                &PRECFOCK= F fast (coarse grid for HF, soft augmentation charge)\n&
                &PRECFOCK= N normal (PREC=N grid for HF, soft augmentation charge)\n&
                &PRECFOCK= A accurate (PREC=A grid for HF, soft augmentation charge)\n\n&
                &L is equivalent to ENCUTFOCK=0 in vasp.5.2.2.\n&
                &M is equivalent to vasp.5.2.2 if ENCUTFOCK is not set.\n&
                &M&L cause significant noise in the forces and are no longer recommended.\n&
                &N is recommended for routine calculations, as there is little noise to be expected.\n&
                &A is now recommended for very accurate calculations.\n&
                &F is now recommended for quick calculations (even phonons often ok, &
                & speedup between 2 and 4).'
        case(PRECFOCK)
            res = 'Use PRECFOCK to select the mode for HF type calculations:\n&
                &PRECFOCK= L low (coarse grid for HF, normal augmentation charge)\n&
                &PRECFOCK= M medium (normal grid for HF, normal augmentation charge)\n&
                &PRECFOCK= F fast (coarse grid for HF, soft augmentation charge)\n&
                &PRECFOCK= N normal (PREC=N grid for HF, soft augmentation charge)\n&
                &PRECFOCK= A accurate (PREC=A grid for HF, soft augmentation charge)\n\n&
                &L is equivalent to ENCUTFOCK=0 in vasp.5.2.2.\n&
                &M is equivalent to vasp.5.2.2 if ENCUTFOCK is not set.\n&
                &M&L cause significant noise in the forces and are no longer recommended.\n&
                &N is recommended for routine calculations, as there is little noise to be expected.\n&
                &A is now recommended for very accurate calculations.\n&
                &F is now recommended for quick calculations (even phonons often ok, &
                & speedup between 2 and 4).'
        case(KPOINTS_HF)
            res = 'Your generating k-point grid is not commensurate with the symmetry &
                &of the lattice. This can cause slow convergence with respect &
                &to k-points for HF-type calculations.\n&
                &Suggested SOLUTIONS:\n&
                & ) If not already the case, use automatic k-point generation,\n&
                & ) and/or shift your grid to Gamma (G) (e.g. required for\n&
                &   hexagonal, fcc, and fcc-orthorhombic lattices).\n&
                & ) and/or set KBLOWUP = .FALSE.'
        case(KPOINTS_PEAD)
            res = 'Your generating k-point grid is not commensurate with the symmetry &
                &of the lattice. This does not sit well in combination with the &
                &PEAD routines, sorry ...\n&
                &Suggested SOLUTIONS:\n&
                & ) If not already the case, use automatic k-point generation,\n&
                & ) and/or shift your grid to Gamma (G) (e.g. required for\n&
                &   hexagonal, fcc, and fcc-orthorhombic lattices).\n&
                & ) and/or set KBLOWUP = .FALSE.'
        case(KPOINTS_INTER)
            res = 'Your generating k-point grid is not commensurate with the symmetry &
                &of the lattice. Presently, it is not possible to write the velocity ' // &
                'operator in the extended zone scheme for such k-point grids.\n&
                &Suggested SOLUTIONS:\n&
                & ) If not already the case, use automatic k-point generation,\n&
                & ) and/or shift your grid to Gamma (G) (e.g. required for\n&
                &   hexagonal, fcc, and fcc-orthorhombic lattices).\n&
                & ) and/or set KBLOWUP = .FALSE.'
        case(GWkweight)
            res = 'The weights do not sum up to one for the summation over q in the &
                &GW routine. To save your life, they have been renormalized to 1. &
                &This circumstance can occur under certain conditions. &
                &Please check your results carefully.'
        case(GWnoreal)
            res = 'Presently, the real-space GW code only supports Gamma-point-only &
                &calculations. Please use the standard version for arbitrary k-point &
                &grids.'
        case(ACFDTRnotSupported)
            res = 'Presently, the real-space GW code requires scaLAPACK. &
                &This code is designed for massively parallel computations. &
                &Please use the standard code for routine calculations.'
        case(PEADnoVirtuals)
            res = 'Since you included no unoccupied states, the PEAD routines can not &
                &estimate whether EFIELD_PEAD exceeds the critical field strength. &
                &Please check your results carefully.'
        case(PEADlargeField)
            res = 'One or more components of EFIELD_PEAD are too large for comfort. In all &
                &probability, you are too near to the onset of Zener tunneling.\n' // &
                '\n     e |E dot A_1| = ' // str(arg%rval(1), '(F8.5)') // ' > 1/10, E_g/N_1 = ' // str(arg%rval(4), '(F8.5)') // &
                '\n     e |E dot A_2| = ' // str(arg%rval(2), '(F8.5)') // ' > 1/10, E_g/N_2 = ' // str(arg%rval(5), '(F8.5)') // &
                '\n     e |E dot A_3| = ' // str(arg%rval(3), '(F8.5)') // ' > 1/10, E_g/N_3 = ' // str(arg%rval(6), '(F8.5)') // &
                '\n\nPossible SOLUTIONS:\n&
                & ) Choose a smaller electric field.\n&
                & ) Use a less dense grid of k-points.'
        case(PEAD_NCOREnot1)
            res = 'The PEAD routines do not work for NCORE /= 1 (or NPAR /= NCPU). &
                &Several features internally rely on the PEAD routines, so even if you &
                &have not specified LPEAD=.TRUE., try restarting your job with NCORE=1.'
        case(PEADPolNoMetal)
            res = 'The calculation of the macroscopic polarization by means of the &
                &Berry-phase expressions (LCALCPOL=.TRUE.) requires your system to be &
                &insulating. This does not seem to be the case.\n&
                &VASP will skip this part of your job, sorry...'
        case(PEADEpsNoMetal)
            res = 'The calculation of changes in the macroscopic polarization due to the &
                &application of a finite electric field, by means of the PEAD method &
                &(LCALCEPS=.TRUE. or EFIELD_PEAD/=0), requires your system to be &
                &insulating. This does not seem to be the case.\n&
                &VASP will skip this part of your job, sorry...'
        case(PEADnowIsMetal)
            res = 'The application of a finite electric field, by means of the PEAD &
                &method (LCALCEPS=.TRUE. or EFIELD_PEAD/=0), requires your system to be &
                &insulating. This does not seem to be the case (any longer). You might try &
                &to reduce the size of the electric field and/or the smearing width.'
        case(LINTET)
            res = 'Changing the smearing to ISMEAR=1 (the default) for those parts &
                &of this job that use k-points specified in line-mode.'
        case(GRDSHFT_NOTET)
            res = 'The linear tetrahedron method can not be used together with the &
                &grid-shift interpolation routine.'
        case(GammaK)
            res = 'You are using the Gamma-point-only version with more than one k-point &
                &or some other non-Gamma k-point.'
        case(KPOINTS_IO)
            res = 'Error reading KPOINTS file.\n&
                &The error occurred at line: ' // str(arg%ival(1)) // '.'
        case(NoInversion)
            res = 'Full k-point grid generated. Inversion symmetry is not applied.'
        case(FFTGridIsNotSufficient)
            res = 'Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accurate &
                &calculation. Thus, the results might be wrong. Good settings for NGX, &
                &NGY and NGZ are ' // str(arg%ival(1)) // ', ' // str(arg%ival(2)) // &
                ' and ' // str(arg%ival(3)) //  ', respectively.\n&
                &Mind: This setting results in a small but reasonable wrap-around error. &
                &It is also necessary to adjust these values to the FFT routines you use.'
        case(MetaGGAandForces)
            res = 'You have switched METAGGA and/or ASPHERICAL radial calculations on &
                &but ions are moved based on forces. &
                &METAGGA and LASPH change only the exchange and correlation energy. &
                &You might want to reconsider your choices.\n&
                &Current values: LASPH=' // str(arg%lval(1), compact = .false.) // &
                ' and LMETAGGA=' // str(arg%lval(2), compact = .false.) // '.'
        case(MetaGGAResponse)
            res = 'You have switched on METAGGA and linear-response routines. &
                &This combination is currently not supported. &
                &If you have selected LEPSILON=.TRUE., you can try using LCALCEPS=.TRUE. instead. &
                &In all other cases, VASP will stop. SOOO sorry...'
        case(SmallCSHIFT)
            res = 'You are using a small complex shift. This can yield unreliable results &
               &when the Kramers-Kronig relation is used for systems with spin-orbit &
               &interactions.\nPlease increase CSHIFT in the INCAR file.'
        case(TetNoMetal)
            res = 'Integration of the imaginary part of the dielectric function by means &
               &of the tetrahedron method is not supported for metallic systems. &
               &VASP will skip the calculation of the optical properties.\n&
               &Please rerun with ISMEAR > -4.'
        case(POTIMlarge)
            res = 'Your timestep is larger than 0.1 Angst.\n&
                &For finite differences, this really does not make sense. I will &
                &reset POTIM to 0.015. I recommend to use 0.01 to 0.02 for finite &
                &differences.'
        case(GRDSHFTnotInIBZ)
            res = 'The grid-shift interpolation procedure is not able to find this point in the IBZ.'
        case(GRDSHFTK)
            res = 'The grid-shift interpolation procedure only works for odd KINTER.'
        case(NoOcc)
            res = 'The energy of the highest-lying unoccupied state is less than the &
                &lowest lying occupied states. Are you sure your system contains &
                &unoccupied states?'
        case(NoOddGrid)
            res = 'The grid-shift interpolation procedure does not work for even grid &
                &divisions.'
        case(GammaC)
            res = 'The k-points need to be Gamma centered.'
        case(BSEantires)
            res = 'BSE calculations beyond the Tamm-Dankoff approximation require you &
                &to set the flag LORBITALREAL=.TRUE. in the groundstate as well as &
                &in the BSE calculation.\n&
                &Please redo the groundstate calculation with this flag set.'
        case(BSEconjg)
            res = 'BSE calculations beyond the Tamm-Dankoff approximation require you &
                &to set the flag LORBITALREAL=.TRUE. in the groundstate as well as &
                &in the BSE calculation.\n&
                &Please redo the groundstate calculation with this flag set.'
        case(BSEncv)
            res = 'The number of conduction band and valence band pairs is zero. &
                &This indicates that NBANDSV or NBANDSO are not properly set up. &
                &The VASP default is to set them to the number of occupied orbitals. &
                &Try to set them manually.\n&
                &It is also possible that OMEGAMAX is too small; please increase it.'
        case(Wmissing)
            res = 'The file WXXXX.tmp for XXXX=' // str(arg%ival(1)) // ' can not be read. &
                &Rerun the GW calculations with LRPA= .TRUE. and do not change &
                &ENCUTGW or NBANDS between runs.'
        case(AuxWAVECAR)
            res = 'An auxiliary WAVECAR file was read.\n&
                &W will be calculated using these auxiliary wavefunctions. &
                &The WAVEDER.XXX file must match the WAVECAR file!'
        case(SpinSpiral)
            res = 'To represent the spin spiral you requested, with a kinetic &
                &energy cutoff of ENINI= ' // str(arg%rval(1)) // ' eV, &
                &choose ENMAX > ' // str(arg%rval(3)) // ' eV. &
                &Currently, ENMAX= ' // str(arg%rval(2)) // ' eV.'
        case(HighestBandsOcccupied)
            res = 'Your highest band is occupied at some k-points! Unless you are &
                &performing a calculation for an insulator or semiconductor, without &
                &unoccupied bands, you have included TOO FEW BANDS!! Please increase &
                &the parameter NBANDS in file INCAR to ensure that the highest band &
                &is unoccupied at all k-points. It is always recommended to &
                &include a few unoccupied bands to accelerate the convergence of &
                &molecular-dynamics runs (even for insulators or semiconductors), &
                &since the presence of unoccupied bands improves wavefunction &
                &prediction and helps to suppress ''band-crossings''.' // &
                fermiWeights(arg%ival, arg%rval)
        case(LVTOT)
            res = 'For LVTOT=.TRUE., VASP.5.x writes the TOTAL local potential to &
                &the file LOCPOT. If you want only the Hartree contributions, use &
                &LVHAR=.TRUE. instead. '
        case(vdWKlimes)
            res = 'You have switched on vdW-DFT.\n&
                &This routine was written and supplied by Jiri Klimes. &
                &We recommend that you carefully read and cite the following &
                &publications:\n&
                &  J. Klimes, D.R. Bowler, A. Michelides\n&
                &     J. Phys.: Cond Matt. 22 022201 (2010)\n&
                &  J. Klimes, D.R. Bowler, A. Michelides\n&
                &     Phys. Rev. B. 83, 195131 (2011)\n&
                &and references therein.'
        case(NoKPAR)
            res = 'You have enabled k-point parallelism (KPAR>1) but this does not &
                &work for the interpolation routines.'
        case(FERWEOptics)
            res = 'Integration of the imaginary part of the dielectric function by &
                &means of fixed occupancies (ISMEAR=-2) is not supported. VASP will &
                &skip the calculation of the optical properties.'
        case(OmegaPar)
            res = 'NOMEGAPAR must be smaller than NOMEGA and the used number of MPI &
                & ranks. I am setting NOMEGAPAR to '//str(arg%ival(1))//', the largest possible value.'
        case(TauPar)
            res = 'NTAUPAR must be smaller than NOMEGA and the used number of MPI &
                & ranks. I am setting NTAUPAR to '//str(arg%ival(1))//', the largest possible value.'
        case(TauParCPU)
            res = 'The number of MPI ranks ('//str(arg%ival(1))//') is not&
                  & divisible by NTAUPAR ('//str(arg%ival(2))//') and/or NOMEGAPAR &
                  &('//str(arg%ival(3))//').'
        case(ImagGridParForLTMP2)
            res = 'The present laplace transformed MP2 implentation does not&
                & support parallelization over imaginary time and frequency &
                & points. NTAUPAR and NOMEGAPAR are set to 1 in the following.'
        case(chiGG_NOmega)
            res = 'Presently, the real-time/real-space code requires a minimum of &
                &3 imaginary times.\nPlease set NOMEGA to a sensible number. &
                &Typically 8 imaginary times suffice.'
        case(OmegaGridOMEGAMIN)
            res = 'You have a very small band gap. This results in slow convergence &
                &of the energy with respect to NOMEGA. &
                &Usually, VASP defaults OMEGAMIN to the bandgap which guarantees an &
                &optimal grid for OMEGAGRID=140. &
                &The gap was found to be OMEGAMIN=' // str(arg%rval(1)) // '. &
                &Now, OMEGAMIN is set to 0.02 instead, which effectively removes some &
                &excitations from the spectrum.\n&
                &You can set OMEGAMIN manually in the INCAR, however.'
        case(WannNotFound)
            res = 'The file wannier90.win can not be found.'
        case(TargetNotSet)
            res = 'You have not defined target states for the projection onto the &
                &subspace. VASP will continue and use all states of the localized &
                &basis as targets.\n&
                &You can set target states in the INCAR file with NTARGET_STATES'
        case(NTARGET_STATES)
            res = 'NTARGET_STATES contains more non-zero values than target states &
                &available. VASP will continue and neglect screening effects of &
                &ALL Wannier states! If this is not desired NTARGET_STATES has to &
                &be set properly in the INCAR.'
        case(NCRPALOW)
            res = 'NCRPALOW is not set in INCAR.\n&
                &NCRPALOW should correspond to the lowest band represented by &
                &the localized basis.'
        case(CRPADisentangle)
            res = 'CRPA-disentanglement method needs a minimal basis with the same size &
                &as the target space. Your chosen targets span a ' // str(arg%ival(2)) // &
                '-dimensional sub space in a space with ' // str(arg%ival(1)) // &
                ' dimensions. Change either NTARGET_STATES or the localized basis set. &
                &Alternatively, use the standard CRPA or the weighted method.'
        case(LCAO_INCAR)
            res = 'You have not specified an ion number WANPROJ_I and/or an angular &
                &momentum WANPROJ_L in the INCAR file. This is mandatory if you &
                &want to use LCAOs!!'
        case(SHMEM_CRPA)
            res = 'NCSHMEM>1 is not supported for CRPA calculations.\n&
                &Please set NCSHMEM=1 for all CRPA calculations.'
        case(OmegaGridOmegaTL)
            res = 'The maximum frequency in the frequency grid is smaller than the &
                &largest excitation energy. This can cause errors in the &
                &correlation energy. Please increase NOMEGA to make sure that the &
                &error is sufficiently small.'
        case(NOmega32)
            res = 'OMEGAGRIRD = 140 does not allow for more than 32 grid points. &
                &For this method, the error decays exponentially with the number of &
                &grid points. Typically, NOMEGA=8-16 is sufficient even for small-gap &
                &systems. Please do your own testing.'
        case(OmegaGrid140)
            res = 'Presently, the real-time/real-space code requires either &
                &OMEGAGRID=140 (minimax grids) or OMEGAGRID=150 (least square grids).'
        case(NUPDOWN)
            res = 'NUPDOWN is presently not supported because it can result in a &
                &different EFERMI for up and down electrons. Switching it off now &
                &and continuing... Please check that occupancies are still ok by, e.g., ' // &
                'performing a one-step KS calculation.'
        case(NPAR_IALGO8)
            res = 'IALGO = 8 supports only parallelism over k-points (KPAR) and &
                &distribution of orbitals over cores.\n&
                &Specifically, NPAR must be set to 1 and NCORE should not be set in &
                &the INCAR file.'
        case(SHMEMError)
            res = 'Can not allocate shared memory. &
                &The simple solution is to run using NCSHMEM=1. &
                &However, this can increase the memory requirements significantly. &
                &If possible, increase the available shared memory. &
                &The available memory can be set in the kernel. &
                &/sbin/sysctl kernel{max,all,mni} can be used to determine the &
                &available memory. ''max'' often needs to be increased (half the physical &
                &memory is often a good choice).'
        case(KProjection)
            res = 'The projection LKPROJ is presently only implemented for NCORE=1. &
                &Please set NCORE=1 (or NPAR to the total number of cores) &
                &and restart the calculations.'
        case(Gamma)
            res = 'The file GAMMA is not compatible with the present calculation. &
                &The first line must specify the number of ions, bands and k-points.'
        case(Gamma2)
            res = 'The file GAMMA can not be read. Each set must start with the k-point &
                &number the start and final band for which the correction of the &
                &density matrix is supplied. NK=' // str(arg%ival(1)) // '.'
        case(Alternant)
            res = 'Could not find all extrema of error function. Try to increase the &
                &number of sampling points. Number of points used in this run: &
                &NMAXALT=' // str(arg%ival(1)) // '.'
        case(NoVGVapplied)
            res = 'You have set LVGVCALC=.FALSE. and LVGVAPPL=.TRUE.. The latter will &
                &take precedence, i.e., LVGVCALC will now be set to .TRUE.'
        case(RandomSeed)
            res = 'The Random Number Generator may be seeded as follows:\n&
                &  0 <= SEED1 < 900000000\n&
                &  0 <= SEED2 < 1000000\n&
                &  0 <= SEED3\n&
                &Your current RANDOM_SEED:' // str(arg%ival(1)) // ' ' // &
                str(arg%ival(2)) // ' ' // str(arg%ival(3)) // '.'
        case(NBANDlessNWAN)
            res = 'The total number of bands (' // str(arg%ival(1)) // ') &
                &is smaller than the available number of local projection centers (' // &
                str(arg%ival(2)) // '). For generating Wannier functions, this &
                &does not make sense.\n Please increase the total number of &
                &bands by setting a higher value for NBANDS in the INCAR file. &
                &Alternatively, you might want to lower the number of local &
                &functions. Consult the documentation of the LOCPROJ flag for &
                &further information.'
        case(ELPHtesting)
            res = 'Electron-Phonon testing mode is active (LELPH_TESTING = .TRUE.). &
                &Unwanted behavior might occur. I highly recommend disabling it.'
        case(NOMEGA_DUMPtooLarge)
            res = 'I found a value for NOMEGA_DUMP > NOMEGA, will use NOMEGA_DUMP &
                &= NOMEGA instead.'
        case(SigmaTMBPT)
            res = 'The temperature for the Matsubara axis differs from the one found &
                &in the WAVECAR. This might cause severe integration errors for &
                &GW and RPA calculations. TMBPT should not be used in production &
                &runs.\n&
                &The inverse temperature (beta) is set via the k-point smearing &
                &beta=1/SIGMA.'
        case(FermiFiniteTemp)
            res = 'Finite temperature many-body pertubation algorithms work only &
                &with Fermi smearing weights (ISMEAR=-1).'
        case(FermiIntegration)
            res = 'Number of electrons determined from Greens function with an error of&
                & '//str( arg%rval(1) ) //'. I recommend to increase NOMEGA.'
        case(BlochBasisUsed)
            res = 'You have not specified any basis functions for your CRPA run. &
                &VASP will use Bloch functions. You can define a basis set with the ' // &
                'LOCPROJ tag or using the VASP2WANNIER90 interface.'
        case(TargetOut)
            res = 'You selected more target states than present for this run.\n&
                &Check your NTARGET_STATES setting.'
        case(NCRPA_BANDSout)
            res = 'You selected more CRPA bands than present for this run.\n&
                &Check your NCRPA_BANDS setting.'
        case(LOCPROJnotSet)
            res = 'You have not specified any basis functions with the LOCRPOJ tag &
                &for your CRPA calculation!'
        case(ALGO_GW)
            res = 'ALGO = GW will be deprecated in future releases. Please use &
                &ALGO = EVGW from now on.'
        case(ALGO_GW0)
            res = 'ALGO = GW0 will be deprecated in future releases. Please use &
                &ALGO = EVGW0 from now on.'
        case(ALGO_G0W0)
            res = 'ALGO = G0W0 will be deprecated in future releases. Please use &
                &ALGO = EVGW0 ; NELM=1 from now on.'
        case(ALGO_scGW0)
            res = 'ALGO = scGW0 will be deprecated in future releases. Please use &
                &ALGO = QPGW0 from now on.'
        case(ALGO_scGW)
            res = 'ALGO = scGW will be deprecated in future releases. Please use &
                &ALGO = QPGW from now on.'
        case(POTCARmissing)
            res = 'Cannot find POTCAR file!'
        case(SymprecDependence)
            res = 'The determination of the symmetry of your systems shows a strong &
                &dependence on the tolerance parameter SYMPREC. This is often the result &
                &of "noise" in the structural parameters (lattice vectors and/or &
                &ionic positions). Set NWRITE=3 for more details.'
        case(ALLOC_RESPONSE)
            res = 'Could not allocate body of response function on mpi rank 0 of size: '&
                &//str(arg%ival(1))//' MB. Reducing NOMEGAPAR or using more computing ' // &
                'nodes might solve this problem.'
        case(ALLOC_RESPONSE_DISTR)
            res = 'Could not allocate body of response function on mpi rank 0 of size: '&
                &//str(arg%ival(1))//' MB. Reducing NTAUPAR or using more computing ' // &
                'nodes might solve this problem.'
        case(ALLOC_WAVEFUN)
            res = 'Could not allocate wavefunction on mpi rank 0 of size: '&
                &//str(arg%ival(1))//' MB. Reducing ENCUT and/or using more ' // &
                'computing nodes might solve this problem.'
        case(MAXMEM_LOW)
            res = 'This job will probably crash, due to insufficient memory &
                &available. Available memory per mpi rank: '//str(arg%ival(1))//' MB, &
                &required memory: '//str(arg%ival(2))//' MB. Reducing ' // &
                'NTAUPAR or using more computing nodes might solve this problem.'
        case(AUTOSET_MEM)
            res = 'Failed to automatically set available memory for this job. &
                &Please set MAXMEM by hand, because the default value of '//&
                &str(arg%ival(1))//' MB might be grossly inefficient.'
        case(RPAFORCES_COMPAT)
            res = 'Currently RPA forces (LRPAFORCES = .TRUE.) can be calculated &
                &only with the formalism of euclidean quantum field theory &
                &(imaginary time Greens functions). That is, only ALGO = &
                &RPAR[K] | G0W0R[K] or similar options are supported.'
        case(RPAENCUTGWSOFT)
            res = 'The default for ENCUTGWSOFT is different in this version of VASP. &
                & If you wish to obtain identical results as using pre vasp.6.3 &
                & please manually set ENCUTGWSOFT =' // str(arg%rval(1)) // ' in the INCAR file.'
        case(ORTHONORMALIZATION_FAILED)
            select case (arg%ival(1))
            case (1)
                work = 'LAPACK: Routine ZPOTRF failed!'
            case (2)
                work = 'LAPACK: Routine ZTRTRI failed!'
            case (3)
                work = 'scaLAPACK: Routine ZPOTRF ZTRTRI failed!'
            case default
                work = 'unknown argument'
            end select
            res = 'Orbital orthonormalization failed in the inversion of matrix\n' &
                // work // ' kpoint: ' // str(arg%ival(2)) // ' spin: ' // str(arg%ival(3)) // '\n\n&
                &Possible solutions for this issue are:\n&
                &- Please check whether the atoms are too close to each other.\n&
                &- Decreasing POTIM might help during relaxations or MDs.\n&
                &- Deleting CHGCAR or WAVECAR if incompatible with the structure.\n&
                &- Try calculating without symmetrizing (ISYM = 0 or ISYM = -1).\n&
                &- As last resort: using a different PAW for some of the atoms.'
        case(DISTRIBUTE_RESPONSE_SUPER_MU_NGVECTOR)
            res = 'DISTRIBUTE_RESPONSE_SUPER_MU: NGVECTOR not correct'
        case(LOPTICS_SPLINES)
            res = 'You have selected k-p perturbation theory (LOPTICS = T or WAVEDER present)&
                  & as well as the spline interpolation method (ESF_SPLINES = T) to calculate&
                  & the long-wave limit q->0 of the electronic structure factor.&
                  & Selecting both methods is not recommened, because the former deteriorates&
                  & k-point convergence of the spline method and increases the&
                  & computational cost of the job.'
        case(ACFDTR_OpenACC)
            res = 'The gamma-only version for ALGO = ACFDTR|RPAR is not fully ported to GPU yet. &
                &Performance degradation can be expected. You might consider running the complex &
                &version of vasp (vasp_std) instead.'
        case(SinglesNaturalo)
            res = 'Renormalized singles correction LSINGLES=T is consistently calculated only &
                &for NATURALO = 2. I have found NATURALO='//str(arg%ival(1))//' in INCAR.'
        case(NCOREAndGWR)
            res = 'Low-scaling GW/ACFDT calculations support only NCORE=1. I found a value of &
               &'//str(arg%ival(1))//' in INCAR.'
        case(NonCollinearAndGWR)
            res = 'Sorry, non-collinear low-scaling GW/ACFDT calculations are not supported yet.'
        case(ElectronPhononRequiresHDF5)
            res = 'Electron-phonon calculations require VASP to be compiled with HDF5 file support'
        case(FockStdGWR)
            res = 'For accurate exact exchange energies (FHF) set LFOCKSTD=T.'
        case(LuttingerWardNaturalo)
            res = 'Luttinger-Ward functional is consistently calculated only for NATURALO = 0. &
                &I have found NATURALO='//str(arg%ival(1))//' in INCAR.'
        case(LOPTICSAndGW)
            res = 'Calculation of the dielectric matrix after GW requires LPEAD=T.'
        case default
            res = 'undefined tutor message ' // str(topic)
        end select
    end function generateMessage

    function fermiWeights(kptIndex, weight) result(res)
        use string, only: str
        integer, allocatable, intent(in) :: kptIndex(:)
        real(q), allocatable, intent(in) :: weight(:)
        character(len=:), allocatable :: res
        integer iband
        logical writeWeights
        writeWeights = .false.
        if (allocated(kptIndex) .and. allocated(weight)) then
            writeWeights = size(kptIndex) == size(weight)
        end if
        if (writeWeights) then
            res = '\n\nIn the following, all k-points are listed with the &
                &Fermi weights of the highest band given in paranthesis:\n'
            do iband = 1, size(kptIndex) - 1
                res = res // str(kptIndex(iband), '(3X,I5)') // ' (' // &
                    str(weight(iband), '(F8.5)') // ')\n'
            end do
            iband = size(kptIndex)
            res = res // 'The total occupancy of band no. ' // str(kptIndex(iband)) // &
                ' is ' // str(weight(iband)) // ' electrons...'
        else
            res = ''
        end if
    end function fermiWeights

    subroutine deactivateVtutorOutputForTest()
        vtutor%discardOutputToNegativeUnit = .true.
        vtutor%unitOut = -1
        vtutor%unitErr = -1
        vtutor%unitBug = -1
    end subroutine deactivateVtutorOutputForTest

    type(argument) function newArgument(ival, lval, rval, cval)
        integer, intent(in), optional :: ival(:)
        logical, intent(in), optional :: lval(:)
        real(q), intent(in), optional :: rval(:)
        character(20), intent(in), optional :: cval(:)
        if (present(ival)) &
            allocate(newArgument%ival, source=ival)
        if (present(lval)) &
            allocate(newArgument%lval, source=lval)
        if (present(rval)) &
            allocate(newArgument%rval, source=rval)
        if (present(cval)) &
            allocate(newArgument%cval, source=cval)
    end function newArgument

end module tutor
