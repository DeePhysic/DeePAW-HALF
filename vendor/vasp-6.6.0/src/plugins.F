#include "symbol.inc"
!> Interface VASP calculation with external plugins.
!>
!> This module contains the tools to convert data from VASP types to C structs
!> and calls a external C code that can then interface with other libraries or
!> codes.
!>
!> There are a few design choices here that came up during the development
!> process:
!>
!> * All preprocessor statements are limited to this file so that the rest of
!>   the code can remain without modification. There is a constant `use_plugins`
!>   defined in this module that can be used in place of the preprocessor
!>   statements in other module if it turns out to be necessary. To facilitate
!>   this design all public procedures will turn to empty placeholders if the
!>   code is compiled without the preprocessor flag.
!>
!> * The interface routines map closely to current functionality of VASP. For
!>   example, the local-potential update is called from within `POTLOK` and the
!>   force update from within `FORCE_AND_STRESS`. With this choice the plugins
!>   are less likely to break existing VASP functionality.
!>
!> * Return values of plugins are additions to the current values. This choice
!>   allows multiple plugins running on the same interface at the same time
!>   without conflict (unless the plugins themselves break this functionality).
!>
!> Adding a new interface requires the following steps:
!>
!> * Create a new function in this module to convert the VASP data to C structs.
!> * Define the C structs in both the C header and the Fortran file (folder plugins)
!> * Define the dataclasses in the Python package
!> * Please keep the interfaces in alphabetical order.
!>
!> All interface routines have a common logic to them:
!>
!> * They perform a sanity check.
!> * They allocate the necessary memory in a buffer.
!> * They encode the immutable VASP data in a object called constants.
!> * They encode the mutable VASP data in a object called additions.
!> * They call the interface to the plugin.
!> * They handle the error if it is nonzero.
!> * They decode the returned additions object and apply the changes to the VASP data.
!>
module plugins

    use base, only: q, plugins_settings, info_struct
    use iso_c_binding, only: c_int, c_ptr
    use mgrid_struct_def, only: grid_3d
    use mpimy, only: communic
    use poscar_struct_def, only: latt, type_info
    use pseudo_struct_def, only: potcar
    use verlet_cell_list_mod, only: verlet_cell_list_type

#ifdef PLUGINS
    use plugins_interface
#endif

    implicit none

    !> Use this in code outside of this module instead of the preprocessor flags.
#ifdef PLUGINS
    logical, parameter, public :: use_plugins = .true.
#else
    logical, parameter, public :: use_plugins = .false.
#endif

    ! internal variables to handle errors when the module is initialized/finalized twice
    logical :: initialized = .false., finalized = .false.

    !
    ! private helper types to contain the memory during the communication
    !

    type, abstract :: buffer_base
        type(communic) comm
        logical broadcast, active_rank
        type(verlet_cell_list_type), pointer :: verlet_list => null()
#ifdef PLUGINS
        type(neighbor_list_data), allocatable :: neighbor_list(:)
#endif
    end type buffer_base

    type, extends(buffer_base) :: buffer_force_and_stress
        integer(c_int), allocatable :: atomic_numbers(:)
        real(q), allocatable :: ZVAL(:)
        real(q), allocatable :: charge_density(:,:,:)
        real(q), allocatable :: forces(:,:)
        real(q), allocatable :: stress(:,:)
    end type buffer_force_and_stress

    type, extends(buffer_base) :: buffer_local_potential
        integer(c_int), allocatable :: atomic_numbers(:)
        real(q), allocatable :: ZVAL(:)
        real(q), allocatable :: charge_density(:,:,:)
        real(q), allocatable :: hartree_potential(:,:,:)
        real(q), allocatable :: ion_potential(:,:,:)
        real(q), allocatable :: potential(:,:,:)
        real(q), allocatable :: dipole_moment(:)
    end type buffer_local_potential

    type, extends(buffer_base) :: buffer_structure
        real(q), allocatable :: lattice_vectors(:,:)
        real(q), allocatable :: positions(:,:)
        integer(c_int), allocatable :: atomic_numbers(:)
        real(q), allocatable :: charge_density(:,:,:)
    end type buffer_structure

    type, extends(buffer_base) :: buffer_occupancies
    end type buffer_occupancies

    private
    public start_plugins, stop_plugins, plugins_reader, write_plugins, &
        plugins_force_and_stress, plugins_local_potential, plugins_machine_learning, &
        plugins_occupancies, plugins_structure

contains

    !> This routine must be called before any of the interface routines. You cannot call
    !> this routine a second time to reinitialize after a prior `stop_plugins` call.
    subroutine start_plugins()
        use tutor, only:vtutor
        integer(c_int) :: error
        if (initialized) then
            call vtutor%bug("Plugins are already initialized, do not initialize them another time.", &
                __FILE__, __LINE__)
        end if
#ifdef PLUGINS
        error = plugins_initialize()
        if (error /= NO_ERROR) &
            call vtutor%error("Error initializing plugins")
#endif
        initialized = .true.
    end subroutine start_plugins


    !> This routine should be called at the end of the code. You must not call further
    !> interface routines after calling this function.
    subroutine stop_plugins()
        use tutor, only: vtutor
        integer(c_int) :: error
        if (.not.initialized) then
            call vtutor%bug("Finalize called before initialize.", __FILE__, __LINE__)
        end if
        if (finalized) then
            call vtutor%bug("Plugins are already finalized, do not finalize them another time.", &
                __FILE__, __LINE__)
        end if
#ifdef PLUGINS
        error = plugins_finalize()
        if (error /= NO_ERROR) &
            call vtutor%error("Error finalizing plugins")
#endif
    end subroutine stop_plugins


    !> Read the INCAR file and detect which interfaces should be activated.
    !> You may use existing INCAR tags to activate plugins where this makes
    !> sense like it is done for IBRION.
    subroutine plugins_reader(incar, ibrion, settings)
        use incar_reader, only: incar_file, process_incar
        use string, only: lowercase
        use tutor, only: vtutor
        type(incar_file), intent(inout) :: incar
        integer, intent(inout) :: ibrion
        type(plugins_settings), intent(out) :: settings
        logical update_structure, machine_learning
        !
        settings%mode = "--"
        call process_incar(incar, "PLUGINS/MODE", settings%mode)
        settings%mode = lowercase(settings%mode)
        call read_flag(incar, "PLUGINS/FORCE_AND_STRESS", settings%force_and_stress)
        call read_flag(incar, "PLUGINS/LOCAL_POTENTIAL", settings%local_potential)
#ifndef NO_MACHINE_LEARNING
        machine_learning = .false.
        call read_flag(incar, "PLUGINS/MACHINE_LEARNING", machine_learning)
        ! Set default of PLUGINS/ML_OUTPUT_MODE and PLUGINS/ML_MODE
        if (machine_learning) then
            settings%output_mode = 0
            settings%ml_mode = "run-legacy"
        else
            settings%output_mode = 1
            settings%ml_mode = "none"
        endif
        call process_incar(incar, "PLUGINS/ML_OUTPUT_MODE", settings%output_mode)
        call process_incar(incar, "PLUGINS/ML_OUTBLOCK", settings%outblock)
        call process_incar(incar, "PLUGINS/ML_MODE", settings%ml_mode)
        settings%ml_mode = lowercase(settings%ml_mode)
#endif
        call read_flag(incar, "PLUGINS/OCCUPANCIES", settings%occupancies)
        update_structure = .false.
        call read_flag(incar, "PLUGINS/STRUCTURE", update_structure)
        settings%structure = ibrion == 12 .or. update_structure
        call process_incar(incar, "PLUGINS/NEIGHBOR_CUTOFF", settings%neighbor_cutoff)
        !
        ! sanity checks
        select case(settings%mode)
        case ("--")
            settings%mode = "serial"
        case ("serial", "parallel", "vasp")
            if (.not.use_plugins) &
                call error_no_plugin_support("PLUGINS/MODE")
        case default
            call vtutor%error("PLUGINS/MODE = " // settings%mode // " is not " // &
                "implemented. Please check the spelling and correct the selection.")
        end select
        !
        select case (settings%ml_mode)
        case ("run", "run-legacy", "delta")
            settings%force_field = .true.
        case ("none")
            settings%force_field = .false.
        case default
            call vtutor%error("plugins_reader: unknown mode selected for PLUGINS/ML_MODE.")
        end select
        if (machine_learning) then
            call vtutor%alert("Please note that PLUGINS/MACHINE_LEARNING is deprecated because of a &
                &unit mismatch for the stress between the underlying ASE implementation and VASP. You &
                &can now set PLUGINS/FORCE_AND_STRESS = T and PLUGINS/ML_MODE = RUN for an analogous &
                &setup with correct treatment of the stress. Please check out\n&
                &https://vasp.at/wiki/PLUGINS/ML_MODE and\n&
                &https://vasp.at/wiki/PLUGINS/FORCE_AND_STRESS\n&
                &for more details.")
        end if
        if (settings%outblock < 1) then
            call vtutor%error("plugins_reader: PLUGINS/OUTBLOCK must be larger than 0.")
        endif
        !
        if (.not.use_plugins.and.settings%structure) then
            call vtutor%error("IBRION = 12 requires PLUGINS. Please recompile " // &
                "the code with -DPLUGINS and rerun the calculation!")
        end if
        if (update_structure .and. (ibrion /= 0 .and. ibrion /= 12)) &
            call vtutor%error("IBRION must not be set if 'PLUGINS/STRUCTURE' is set")
        !
        if (settings%structure) ibrion = 12
        if (.not.initialized) &
            call start_plugins()
        !
    end subroutine plugins_reader

    subroutine read_flag(incar, key, flag)
        use incar_reader, only: incar_file, process_incar
        type(incar_file), intent(inout) :: incar
        character(len=*), intent(in) :: key
        logical, intent(inout) :: flag
        !
        call process_incar(incar, key, flag)
        if (.not.use_plugins.and.flag) &
            call error_no_plugin_support(key)
    end subroutine read_flag

    subroutine error_no_plugin_support(key)
        use tutor, only: vtutor
        character(len=*), intent(in) :: key
        call vtutor%error(key // " was set in the INCAR file but the code " // &
            "was not compiled with PLUGINS support. Please recompile with " // &
            "-DPLUGINS and rerun the calculation!")
    end subroutine error_no_plugin_support

    !> Report the settings to the OUTCAR file
    subroutine write_plugins(unit_, settings)
        integer, intent(in) :: unit_
        type(plugins_settings), intent(in) :: settings
        character(len=*), parameter :: format_ = "(3x,a22,l6,4x,a)"
        character(len=:), allocatable :: explanation
        !
        write(unit_, *) "Python interface"
        if (use_plugins) then
            write(unit_, format_) "PLUGINS active      = ", .true., "the code was linked to PLUGINS"
            select case (settings%mode)
            case ("serial")
                explanation = "run plugins only on root process"
            case ("parallel")
                explanation = "run plugins in parallel but do not distribute data"
!            case ("vasp")
!                explanation = "run plugins in parallel with VASP's data distribution"
            end select
            write(unit_, '(3x,a22,a,t36,a)') "  -mode             = ", settings%mode, explanation
            write(unit_, format_) "  -force_and_stress = ", settings%force_and_stress, "update the forces and the stress externally"
#ifndef NO_MACHINE_LEARNING
            select case (settings%ml_mode)
            case ("none")
                explanation = "add to ab initio forces and stress"
            case ("run", "delta")
                explanation = "replace ab initio forces with plugin"
            case ("run-legacy")
                explanation = "legacy consistency with PLUGINS/MACHINE_LEARNING"
            end select
            write(unit_, '(3x,a22,a,t36,a)') "  -ml_mode          = ", settings%ml_mode, explanation
#endif
            write(unit_, format_) "  -local_potential  = ", settings%local_potential, "update the potential externally"
            write(unit_, format_) "  -occupancies      = ", settings%occupancies, "update quantities related to the occupancies"
            write(unit_, format_) "  -structure        = ", settings%structure, "overwrite the structure externally"
            if (settings%neighbor_cutoff <= 0) then
                explanation = "neighbor list not used"
            else
                explanation = "cutoff radius for neighbor list"
            end if
            write(unit_, '(3x,a22,f8.2,t36,a)') "  -neighbor_cutoff  = ", settings%neighbor_cutoff, explanation
        else
            write(unit_, format_) "PLUGINS active      = ", .false., "recompile with -DPLUGINS to activate"
        end if
    end subroutine write_plugins

    !
    ! The interfaces to the plugins
    !

    !> This routine allows to modify the force and the stress before the ionic update
    !>
    !> Note: if the charge density is not present, ML interface is called instead
    subroutine plugins_force_and_stress(info, t_info, latt_cur, p, gridc, &
            energy, forces, stress, charge_density)
        type(info_struct), intent(in) :: info
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        real(q), intent(inout) :: energy
        real(q), intent(inout) :: forces(:,:)
        real(q), intent(inout) :: stress(:,:)
        complex(q), intent(in), optional :: charge_density(:)
        integer(c_int) error
        !
#ifdef PLUGINS
#ifndef NO_FORCE_AND_STRESS
        type(constants_force_and_stress) constants
        type(additions_force_and_stress) additions
        type(buffer_force_and_stress), volatile :: buffer
        !   volatile because it may be modified by the external program
        !
        call sanity_check_on_data_types_and_status()
        call allocate_buffer_force_and_stress(info%plugin%mode, t_info%nions, p, gridc, buffer, charge_density)
        call create_neighbor_list(info%plugin%neighbor_cutoff, t_info, latt_cur, buffer)
        if (buffer%active_rank) then
            call encode_constants_force_and_stress(buffer, gridc, info, t_info, latt_cur, forces, stress, constants)
            call encode_additions_force_and_stress(buffer, additions)
            error = interface_force_and_stress(constants, additions)
        end if
        call error_handling(buffer, "plugins_force_and_stress", error)
        call decode_additions_force_and_stress(buffer, additions, energy, forces, stress)
#endif
#endif
    end subroutine plugins_force_and_stress


    !> This interface allows to add an external local potential in real space
    subroutine plugins_local_potential(info, t_info, latt_cur, p, gridc, charge_density, &
            potential, hartree_potential, ion_potential, dipole_moment, energy)
        type(info_struct), intent(in) :: info
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        complex(q), intent(in) :: charge_density(:)
        complex(q), intent(inout) :: potential(:)
        complex(q), intent(in) :: hartree_potential(:)
        complex(q), intent(in) :: ion_potential(:)
        real(q), intent(in)    :: dipole_moment(:)
        real(q), intent(inout) :: energy
#ifdef PLUGINS
#ifndef NO_LOCAL_POTENTIAL
        type(constants_local_potential) :: constants
        type(additions_local_potential) :: additions
        type(buffer_local_potential), volatile :: buffer
        integer(c_int) error
        !
        call sanity_check_on_data_types_and_status()
        call allocate_buffer_local_potential(info%plugin%mode, t_info%nions, p, gridc, charge_density, &
                                             hartree_potential, ion_potential, dipole_moment, buffer)
        call create_neighbor_list(info%plugin%neighbor_cutoff, t_info, latt_cur, buffer)
        if (buffer%active_rank) then
            call encode_constants_local_potential(buffer, info, t_info, latt_cur, gridc, constants)
            call encode_additions_local_potential(buffer, additions)
            error = interface_local_potential(constants, additions)
        end if
        call error_handling(buffer, "plugins_local_potential", error)
        call decode_additions_local_potential(buffer, additions, gridc, energy, potential)
#endif
#endif
    end subroutine plugins_local_potential


    !> This interface allows to provide an external machine learned force field
    subroutine plugins_machine_learning(info, t_info, latt_cur, p, gridc, dyn, &
        symm, isymop, nrotk, ldo_ab_initio, energy, forces, stress, press)
        use base, only: symmetry
        use constant, only: EVTOJ
        use poscar_struct_def, only: dynamics
        !
        type(info_struct), intent(in) :: info
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        type(dynamics), intent(in) :: dyn
        type(symmetry), intent(in) :: symm
        integer, intent(in) :: isymop(:,:,:), nrotk
        logical, intent(in) :: ldo_ab_initio
        real(q), intent(out) :: energy, forces(:,:), stress(:,:), press
         !
#ifdef PLUGINS
#ifndef NO_MACHINE_LEARNING
         !
         energy = 0.0_q
         forces = 0.0_q
         stress = 0.0_q
         !
         PROFILING_START("plugin_machine_learning")
         call plugins_force_and_stress(info, t_info, latt_cur, p, gridc, energy, forces, stress)
         PROFILING_STOP("plugin_machine_learning")
         !
         if (symm%isym > 0) then
            call forsym(forces, symm%rotmap(1,1,1), t_info%ntyp, t_info%nityp, t_info%niond, &
                symm%taurot, symm%wrkrot, latt_cur%a)
            if (dyn%isif /= 0) call tsym(stress, isymop, nrotk, latt_cur%a)
            press = (stress(1,1) + stress(2,2) + stress(3,3)) / 3._q &
                - dyn%pstress / (EVTOJ * 1e22_q) * latt_cur%omega
         end if
#endif
#endif
         !
    end subroutine plugins_machine_learning


    !> This interface allow to change quantities related to the occupancies at the end of each ionic step
    subroutine plugins_occupancies(info, comm, nelect, efermi, efermi_set, sigma, ismear, emin, emax, nupdown)
        type(info_struct), intent(in) :: info
        type(communic), intent(in) :: comm
        real(q), intent(inout) :: nelect
        real(q), intent(in) :: efermi
        real(q), intent(inout) :: efermi_set
        real(q), intent(inout) :: sigma
        integer, intent(inout) :: ismear
        real(q), intent(inout) :: emin
        real(q), intent(inout) :: emax
        real(q), intent(inout) :: nupdown
        integer(c_int)         :: error
#ifdef PLUGINS
#ifndef NO_OCCUPANCIES
        type(constants_occupancies) :: constants
        type(additions_occupancies) :: additions
        type(buffer_occupancies), volatile :: buffer
        call sanity_check_on_data_types_and_status()
        call allocate_buffer_occupancies(info%plugin%mode, comm, buffer)
        if (buffer%active_rank) then
            call encode_constants_occupancies(nelect, efermi, sigma, ismear, emin, emax, nupdown, constants)
            call encode_additions_occupancies(additions)
            error = interface_occupancies(constants, additions)
        end if
        call error_handling(buffer, "plugins_occupancies", error)
        call decode_additions_occupancies(buffer, additions, nelect, efermi, efermi_set, sigma, ismear, emin, emax, nupdown)
#endif
#endif
    end subroutine plugins_occupancies


    !> This interface allows to update the structure at the end of an ionic step
    !> replacing the _normal_ update procedure.
    subroutine plugins_structure(gridc, info, t_info, latt_cur, p, comm, energy, forces, stress, density)
        type(grid_3d), intent(in) :: gridc
        type(info_struct), intent(in) :: info
        type(type_info), intent(inout) :: t_info
        type(latt), intent(inout) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(communic), intent(in) :: comm
        real(q), intent(in) :: energy
        real(q), intent(inout) :: forces(:,:)
        real(q), intent(inout) :: stress(:,:)
        complex(q), allocatable :: density(:,:)
        complex(q), allocatable :: cwork(:)
        integer(c_int) error
#ifdef PLUGINS
#ifndef NO_STRUCTURE
        type(constants_structure) :: constants
        type(additions_structure) :: additions
        type(buffer_structure), volatile    :: buffer
        real(q), allocatable :: positions(:,:)
        real(q), allocatable :: lattice_vectors(:,:)
        !
        allocate(positions, source=t_info%posion)
        allocate(lattice_vectors, source=latt_cur%a)
        positions = 0.0_q ; lattice_vectors = 0.0_q
        call sanity_check_on_data_types_and_status()
        if (allocated(density)) then
            allocate(cwork(gridc%MPLWV))
            cwork = 0.0_q
            call rc_add(density(1,1), 1.0_q, &
                        density(1,1), 0.0_q, &
                        cwork(1), gridc)
            call FFT3D(cwork, gridc, 1)
        else
            allocate(cwork(0))
        end if
        call allocate_buffer_structure(gridc, info%plugin%mode, comm, p, positions, lattice_vectors, cwork, buffer)
        call create_neighbor_list(info%plugin%neighbor_cutoff, t_info, latt_cur, buffer)
        if (buffer%active_rank) then
            call encode_constants_structure(buffer, gridc, info, t_info, latt_cur, energy, forces, stress, constants)
            call encode_additions_structure(buffer, additions)
            error = interface_update_structure(constants, additions)
        end if
        call error_handling(buffer, "plugins_update_structure", error)
        call decode_additions_structure(buffer, latt_cur, t_info)
#endif
#endif
    end subroutine plugins_structure

    !
    ! common helper routines used by all interfaces
    !

    subroutine sanity_check_on_data_types_and_status
        use iso_c_binding, only: c_double
        use tutor, only: vtutor
        real(q) double_fortran
        real(c_double) double_cpp
        integer int_fortran
        integer(c_int) int_cpp
        !
        if (.not.initialized) then
            call vtutor%bug("Please call `start_plugins` before calling any interface to a plugin", &
                __FILE__, __LINE__)
        end if
        if (finalized) then
            call vtutor%bug("Cannot call any interface after finalizing the plugins with `stop_plugins`", &
                __FILE__, __LINE__)
        end if
        if (kind(double_fortran) /= kind(double_cpp)) then
            call vtutor%bug("The bytesize of floating-point numbers is different. " &
                // "Please check the setup of the VASP datatypes!", __FILE__, __LINE__)
        end if
        if (kind(int_fortran) /= kind(int_cpp)) then
            call vtutor%bug("The bytesize of integer numbers is different. " &
                // "Please check the setup of the VASP datatypes!", __FILE__, __LINE__)
        end if
    end subroutine sanity_check_on_data_types_and_status

    subroutine error_handling(buffer, routine, error_this_rank)
        use string, only: str
        use tutor, only: vtutor
        class(buffer_base), intent(in) :: buffer
        character(len=*), intent(in) :: routine
        integer(c_int), intent(in) :: error_this_rank
        integer error_any_rank
        !
#ifdef PLUGINS
        if (buffer%active_rank) then
            error_any_rank = error_this_rank
        else
            error_any_rank = NO_ERROR
        end if
        CALLMPI( M_max_i(buffer%comm, error_any_rank, 1) )
        !
        select case (error_any_rank)
        case (NO_ERROR)
            return
        case (BUG_IN_CPP)
            call vtutor%bug("When calling Python from C++ in " // routine // " an error " // &
                "occurred. Please check the output for possible hints. The most common " // &
                "issues are typos or missing arguments for the Python functions.", &
                __FILE__, __LINE__)
        case (IMPORT_ERROR)
            call vtutor%error("Could not run 'import vasp' in Python. Please make sure that " // &
                "you installed the package and it is in your active environment")
        case (PYTHON_EXCEPTION)
            call vtutor%error("One of the plugins raised an exception in " // routine // &
                ", please check the output for more details. For further clarification, " // &
                "contact the author of the plugin to resolve the issue.")
        case default
            call vtutor%error("Calling the interface " // routine // " returned a unknown " // &
                "error code " // str(error_any_rank) // ". Perhaps the output may give advice on " // &
                "how to proceed. Otherwise your VASP version may be incompatible with the " // &
                "plugin system.")
        end select
#endif
    end subroutine error_handling


    subroutine merge_grid_quantity(mode, gridc, distributed, gathered)
        use tutor, only: vtutor
        character(len=*), intent(in) :: mode
        type(grid_3d), intent(in) :: gridc
        complex(q), intent(in) :: distributed(:)
        real(q), intent(out), allocatable :: gathered(:,:,:)
        RGRID, allocatable :: work(:)
        real(q), allocatable :: rwork(:,:)
        integer :: nz
        !
        select case (mode)
        case ("serial", "parallel")
            allocate(work(DIMREAL(gridc%mplwv)))
            allocate(gathered(gridc%ngz, gridc%ngy, gridc%ngx))
            allocate(rwork(gridc%ngx,gridc%ngy))
            work = 0._q
            call rl_add(distributed, 1.0_q, distributed, 0.0_q, work, gridc)
            do nz = 1, gridc%ngz
               call mrg_grid_rl_plane(gridc, rwork, work, nz)
               gathered(nz,:,:) = transpose(rwork)
            end do
        case default
            call vtutor%bug("Merging grid quantities is not implemented for mode '" // mode // "'.", __FILE__, __LINE__)
        end select
    end subroutine merge_grid_quantity

    subroutine setup_buffer_base(buffer, mode, comm)
        class(buffer_base), intent(inout) :: buffer
        character(len=*), intent(in) :: mode
        type(communic), intent(in) :: comm
        buffer%comm = comm
        buffer%broadcast = mode == "serial"
        buffer%active_rank = mode /= "serial" .or. comm%node_me == 1
    end subroutine setup_buffer_base

#ifdef PLUGINS
    subroutine create_neighbor_list(cutoff, t_info, latt_cur, buffer)
        use constant, only: TPI
        use lattice, only: dirkar
        use neighbor_list_mod, only: neighbor_list_type
        use tutor, only: vtutor
        real(q), intent(in) :: cutoff
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        class(buffer_base), intent(inout) :: buffer
        type(verlet_cell_list_type), allocatable, save, target :: neighbor_list
        integer ion, number_ions
        real(q) inverse_lattice(3,3)
        real(q), allocatable :: positions(:,:)
        !
        if (cutoff <= 0) &
            return
        !
#if defined(MPI) || defined(MPI_CHAIN)
        number_ions = size(t_info%posion, 2)
        allocate(positions, source=t_info%posion)
        call dirkar(number_ions, positions, latt_cur%a)
        inverse_lattice = TPI * latt_cur%b
        if (.not.allocated(neighbor_list)) then
            allocate(neighbor_list)
            call neighbor_list%init_nn_cell_list(cutoff, latt_cur%a, positions, t_info%ntyp)
            call neighbor_list%forced_update(positions, t_info%ityp, latt_cur%a, inverse_lattice)
        end if
        call neighbor_list%compute_nn_list(positions, t_info%ityp, latt_cur%a, inverse_lattice)  ! TODO: MPI
        !
        buffer%verlet_list => neighbor_list
        allocate(buffer%neighbor_list(number_ions))
        do ion = 1, number_ions
            buffer%neighbor_list(ion)%size = size(neighbor_list%nn_list(ion)%neighbors)
            buffer%neighbor_list(ion)%neighbors = c_loc(neighbor_list%nn_list(ion)%neighbors)
            buffer%neighbor_list(ion)%distances = c_loc(neighbor_list%nn_list(ion)%neighbor_dist)
            buffer%neighbor_list(ion)%directions = c_loc(neighbor_list%nn_list(ion)%neighbor_connect)
        end do
#else
        call vtutor%error("Constructing a neighbor list requires parallel version of the code. &
            &Please recompile with MPI support.")
#endif`
        !
    end subroutine create_neighbor_list

    type(c_ptr) function encode_neighbor_list(buffer) result (neighbor_list)
        use iso_c_binding, only: c_null_ptr
        class(buffer_base), intent(in), target :: buffer
        if (.not.allocated(buffer%neighbor_list)) then
            neighbor_list = c_null_ptr
        else
            neighbor_list = c_loc(buffer%neighbor_list)
        end if
    end function encode_neighbor_list


#ifndef NO_FORCE_AND_STRESS
    !
    ! implementation details of the force and stress interface
    !

    subroutine allocate_buffer_force_and_stress(mode, number_ions, p, gridc, buffer, charge_density)
        character(len=*), intent(in) :: mode
        integer, intent(in) :: number_ions
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        type(buffer_force_and_stress), intent(out) :: buffer
        complex(q), intent(in), optional :: charge_density(:)
        !
        call setup_buffer_base(buffer, mode, gridc%comm)
        if (present(charge_density)) then
            call merge_grid_quantity(mode, gridc, charge_density, buffer%charge_density)
        else
            allocate(buffer%charge_density(0,0,0))
        end if
        allocate(buffer%atomic_numbers, source=atomic_number(p%element))
        allocate(buffer%ZVAL(size(p)))
        allocate(buffer%forces(3, number_ions), source=0.0_q)
        allocate(buffer%stress(3, 3), source=0.0_q)
        buffer%ZVAL(:) = p(:)%ZVALF
    end subroutine allocate_buffer_force_and_stress

    subroutine encode_constants_force_and_stress(buffer, gridc, info, t_info, latt_cur, forces, stress, constants)
        type(buffer_force_and_stress), intent(in), target :: buffer
        type(grid_3d), intent(in) :: gridc
        type(info_struct), intent(in) :: info
        type(type_info), intent(in), target :: t_info
        type(latt), intent(in), target :: latt_cur
        real(q), intent(in), target :: forces(:,:), stress(:,:)
        type(constants_force_and_stress), intent(out) :: constants
        !
        constants%ENCUT = info%enmax
        constants%NELECT = info%nelect
        if (size(buffer%charge_density) == 0) then
            ! The charge density if passed as an array with size zero to the C++
            constants%shape_grid = [0, 0, 0]
        else
            constants%shape_grid = [gridc%ngx, gridc%ngy, gridc%ngz]
        end if
        constants%number_ions = t_info%nions
        constants%number_ion_types = t_info%ntyp
        constants%ion_types = c_loc(t_info%ityp)
        constants%atomic_numbers = c_loc(buffer%atomic_numbers)
        constants%lattice_vectors = c_loc(latt_cur%a)
        constants%positions = c_loc(t_info%posion)
        constants%ZVAL = c_loc(buffer%ZVAL)
        constants%POMASS = c_loc(t_info%POMASS)
        constants%forces = c_loc(forces)
        constants%stress = c_loc(stress)
        constants%charge_density = c_loc(buffer%charge_density)
#ifndef NO_MACHINE_LEARNING
        if (info%plugin%ml_mode == "run-legacy") then
            constants%mode = MACHINE_LEARNING
        else
            constants%mode = FORCE_AND_STRESS
        end if
#endif
        constants%neighbor_list = encode_neighbor_list(buffer)
    end subroutine encode_constants_force_and_stress

    subroutine encode_additions_force_and_stress(buffer, additions)
        type(buffer_force_and_stress), intent(in), target :: buffer
        type(additions_force_and_stress), intent(out) :: additions
        additions%total_energy = 0.0_q
        additions%forces = c_loc(buffer%forces)
        additions%stress = c_loc(buffer%stress)
    end subroutine encode_additions_force_and_stress

    subroutine decode_additions_force_and_stress(buffer, additions, energy, forces, stress)
        type(buffer_force_and_stress), intent(in) :: buffer
        type(additions_force_and_stress), intent(in) :: additions
        real(q), intent(inout) :: energy, forces(:,:), stress(:,:)
        !
        if (buffer%broadcast) then
            CALLMPI( M_bcast_d(buffer%comm, additions%total_energy, 1) )
            CALLMPI( M_bcast_d(buffer%comm, buffer%forces, size(buffer%forces)) )
            CALLMPI( M_bcast_d(buffer%comm, buffer%stress, size(buffer%stress)) )
        end if
        energy = energy + additions%total_energy
        forces = forces + buffer%forces
        stress = stress + buffer%stress
    end subroutine decode_additions_force_and_stress
#endif

#ifndef NO_LOCAL_POTENTIAL
    !
    ! implementation details of the local potential interface
    !

    subroutine allocate_buffer_local_potential(mode, number_ions, p, gridc, charge_density, &
                                               hartree_potential, ion_potential, &
                                               dipole_moment, buffer)
        character(len=*), intent(in) :: mode
        integer, intent(in) :: number_ions
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        complex(q), intent(in) :: charge_density(:)
        complex(q), intent(in) :: hartree_potential(:)
        complex(q), intent(in) :: ion_potential(:)
        real(q), intent(in)    :: dipole_moment(:)
        type(buffer_local_potential), intent(out) :: buffer
        !
        call setup_buffer_base(buffer, mode, gridc%comm)
        call merge_grid_quantity(mode, gridc, charge_density, buffer%charge_density)
        call merge_grid_quantity(mode, gridc, hartree_potential, buffer%hartree_potential)
        call merge_grid_quantity(mode, gridc, ion_potential, buffer%ion_potential)
        allocate(buffer%dipole_moment, source=dipole_moment)
        allocate(buffer%atomic_numbers, source=atomic_number(p%element))
        allocate(buffer%ZVAL(size(p)))
        allocate(buffer%potential(gridc%ngz, gridc%ngy, gridc%ngx))
        buffer%ZVAL(:) = p(:)%ZVALF
        buffer%potential = 0
    end subroutine allocate_buffer_local_potential

    subroutine encode_constants_local_potential(buffer, info, t_info, latt_cur, gridc, &
            constants)
        type(buffer_local_potential), intent(in), target :: buffer
        type(info_struct), intent(in) :: info
        type(type_info), intent(in), target :: t_info
        type(latt), intent(in), target :: latt_cur
        type(grid_3d), intent(in) :: gridc
        type(constants_local_potential), intent(out) :: constants
        !
        constants%ENCUT = info%enmax
        constants%NELECT = info%nelect
        constants%shape_grid = [gridc%ngx, gridc%ngy, gridc%ngz]
        constants%number_ions = t_info%nions
        constants%number_ion_types = t_info%ntyp
        constants%ion_types = c_loc(t_info%ityp)
        constants%atomic_numbers = c_loc(buffer%atomic_numbers)
        constants%lattice_vectors = c_loc(latt_cur%a)
        constants%positions = c_loc(t_info%posion)
        constants%ZVAL = c_loc(buffer%ZVAL)
        constants%charge_density = c_loc(buffer%charge_density)
        constants%hartree_potential = c_loc(buffer%hartree_potential)
        constants%ion_potential = c_loc(buffer%ion_potential)
        ! Currently the dipole moment is sent only if LDIPOL is set
        constants%dipole_moment_size = size(buffer%dipole_moment)
        constants%dipole_moment = c_loc(buffer%dipole_moment)
        constants%neighbor_list = encode_neighbor_list(buffer)
    end subroutine encode_constants_local_potential

    subroutine encode_additions_local_potential(buffer, additions)
        type(buffer_local_potential), intent(in), target :: buffer
        type(additions_local_potential), intent(out) :: additions
        additions%total_energy = 0.0_q
        additions%potential = c_loc(buffer%potential)
    end subroutine encode_additions_local_potential

    subroutine decode_additions_local_potential(buffer, additions, gridc, &
            energy, potential)
        type(buffer_local_potential), intent(in) :: buffer
        type(additions_local_potential), intent(in) :: additions
        type(grid_3d), intent(in) :: gridc
        real(q), intent(inout) :: energy
        complex(q), intent(inout) :: potential(:)
        complex(q), allocatable :: added_potential(:)
        real(q), allocatable :: work(:,:)
        integer :: ni
        !
        allocate(work(gridc%ngx,gridc%ngy))
        allocate(added_potential(gridc%mplwv))
        added_potential = 0._q
        if (buffer%broadcast) &
            CALLMPI( M_bcast_d(buffer%comm, additions%total_energy, 1) )
        energy = additions%total_energy
        do ni = 1, gridc%ngz
           !
           ! note that the x/z axis is swapped for fortran/c ordering
           !
           work = transpose(buffer%potential(ni,:,:))
           call dis_grid_rl_plane(gridc, work, added_potential, buffer%broadcast, ni)
        end do
        call rl_add(potential, 1._q, added_potential, 1._q, potential, gridc)
    end subroutine decode_additions_local_potential
#endif

#ifndef NO_STRUCTURE
    !
    ! implementation details of the update structure interface
    !

    subroutine allocate_buffer_structure(gridc, mode, comm, p, positions, lattice_vectors, charge_density, buffer)
        type(grid_3d), intent(in) :: gridc
        character(len=*), intent(in) :: mode
        type(communic), intent(in) :: comm
        type(potcar), intent(in) :: p(:)
        real(q), intent(in) :: positions(:,:)
        real(q), intent(in) :: lattice_vectors(:,:)
        complex(q), intent(in) :: charge_density(:)
        type(buffer_structure), intent(out) :: buffer
        !
        call setup_buffer_base(buffer, mode, comm)
        allocate(buffer%lattice_vectors, source=lattice_vectors)
        allocate(buffer%positions, source=positions)
        allocate(buffer%atomic_numbers, source=atomic_number(p%element))
        if (size(charge_density) > 0) then
            call merge_grid_quantity(mode, gridc, charge_density, buffer%charge_density)
        else
            allocate(buffer%charge_density(0,0,0))
        end if
    end subroutine allocate_buffer_structure

    subroutine encode_constants_structure(buffer, gridc, info, t_info, latt_cur, energy, forces, &
            stress, constants)
        type(buffer_structure), intent(in), target :: buffer
        type(grid_3d), intent(in) :: gridc
        type(info_struct), intent(in) :: info
        type(type_info), intent(in), target :: t_info
        type(latt), intent(in), target :: latt_cur
        real(q), intent(in), target :: energy
        real(q), intent(in), target :: forces(:,:), stress(:,:)
        type(constants_structure), intent(out) :: constants
        !
        constants%number_ions = t_info%nions
        constants%number_ion_types = t_info%ntyp
        constants%ion_types = c_loc(t_info%ityp)
        constants%atomic_numbers = c_loc(buffer%atomic_numbers)
        constants%lattice_vectors = c_loc(latt_cur%a)
        constants%positions = c_loc(t_info%posion)
        constants%POMASS = c_loc(t_info%POMASS)
        constants%total_energy = energy
        constants%forces = c_loc(forces)
        constants%stress = c_loc(stress)
        if (size(buffer%charge_density) == 0) then
            ! The charge density if passed as an array with size zero to the C++
            constants%shape_grid = [0, 0, 0]
        else
            constants%shape_grid = [gridc%ngx, gridc%ngy, gridc%ngz]
        end if
        constants%charge_density = c_loc(buffer%charge_density)
        constants%neighbor_list = encode_neighbor_list(buffer)
    end subroutine encode_constants_structure

    subroutine encode_additions_structure(buffer, additions)
        type(buffer_structure), intent(in), target :: buffer
        type(additions_structure), intent(out) :: additions
        additions%lattice_vectors = c_loc(buffer%lattice_vectors)
        additions%positions = c_loc(buffer%positions)
    end subroutine encode_additions_structure

    subroutine decode_additions_structure(buffer, latt_cur, t_info)
        type(buffer_structure), intent(in) :: buffer
        type(latt), intent(inout) :: latt_cur
        type(type_info), intent(inout) :: t_info
        !
        if (buffer%broadcast) then
            CALLMPI( M_bcast_d(buffer%comm, buffer%lattice_vectors, size(buffer%lattice_vectors)) )
            CALLMPI( M_bcast_d(buffer%comm, buffer%positions, size(buffer%positions)) )
        end if
        latt_cur%A = latt_cur%A + buffer%lattice_vectors
        t_info%posion = t_info%posion + buffer%positions
    end subroutine decode_additions_structure
#endif

#ifndef NO_OCCUPANCIES
    subroutine allocate_buffer_occupancies(mode, comm, buffer)
        character(len=*), intent(in) :: mode
        type(communic), intent(in) :: comm
        type(buffer_occupancies) :: buffer
        call setup_buffer_base(buffer, mode, comm)
    end subroutine

    subroutine encode_constants_occupancies(nelect, efermi, sigma, ismear, emin, emax, nupdown, &
                                            constants)
        real(q), intent(in) :: nelect
        real(q), intent(in) :: efermi
        real(q), intent(in) :: sigma
        integer, intent(in) :: ismear
        real(q), intent(in) :: emin
        real(q), intent(in) :: emax
        real(q), intent(in) :: nupdown
        type(constants_occupancies), intent(out) :: constants

        constants%nelect = nelect
        constants%efermi = efermi
        constants%sigma = sigma
        constants%ismear = ismear
        constants%emin = emin
        constants%emax = emax
        constants%nupdown = nupdown
    end subroutine encode_constants_occupancies

    subroutine encode_additions_occupancies(additions)
        type(additions_occupancies), intent(out) :: additions

        additions%nelect = 0.0_q
        additions%efermi = 0.0_q
        additions%sigma = 0.0_q
        additions%ismear = 0
        additions%emin = 0.0_q
        additions%emax = 0.0_q
        additions%nupdown = 0.0_q
    end subroutine encode_additions_occupancies

    subroutine decode_additions_occupancies(buffer, additions, nelect, efermi, efermi_set, sigma, ismear, emin, emax, &
                                            nupdown)
        type(buffer_occupancies), intent(in) :: buffer
        type(additions_occupancies), intent(out) :: additions
        real(q), intent(inout) :: nelect
        real(q), intent(in) :: efermi
        real(q), intent(inout) :: efermi_set
        real(q), intent(inout) :: sigma
        integer, intent(inout) :: ismear
        real(q), intent(inout) :: emin
        real(q), intent(inout) :: emax
        real(q), intent(inout) :: nupdown

        if (buffer%broadcast) then
            CALLMPI( M_bcast_d(buffer%comm, additions%nelect, 1) )
            CALLMPI( M_bcast_d(buffer%comm, additions%efermi, 1) )
            CALLMPI( M_bcast_d(buffer%comm, additions%sigma, 1) )
            CALLMPI( M_bcast_d(buffer%comm, additions%ismear, 1) )
            CALLMPI( M_bcast_d(buffer%comm, additions%emin, 1) )
            CALLMPI( M_bcast_d(buffer%comm, additions%emax, 1) )
            CALLMPI( M_bcast_d(buffer%comm, additions%nupdown, 1) )
        end if
        nelect = nelect + additions%nelect
        if ((additions%efermi /= 0.0_q)) then
            efermi_set = efermi + additions%efermi
        end if
        sigma = sigma + additions%sigma
        ismear = ismear + additions%ismear
        emin = emin + additions%emin
        emax = emax + additions%emax
        nupdown = nupdown + additions%nupdown
    end subroutine decode_additions_occupancies
#endif
#endif

    ! this routine should perhaps be moved to poscar.F

    !> Map an ion type onto the atomic number
    elemental integer function atomic_number(element)
        use string, only: lowercase
        character(len=2), intent(in) :: element
        character(len=2), parameter :: elements(94) = ['h ', 'he', &
            'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne', &
            'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', &
            'k ', 'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', &
                'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', &
            'rb', 'sr', 'y ', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', &
                'cd', 'in', 'sn', 'sb', 'te', 'i ', 'xe', &
            'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', &
                'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w ', 're', 'os', 'ir', 'pt', &
                'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', &
            'fr', 'ra', 'ac', 'th', 'pa', 'u ', 'np', 'pu']
        integer, parameter :: unknown_atom = -1
        character(len=2) element_
        integer i
        element_ = lowercase(element)
        atomic_number = unknown_atom
        do i = 1, size(elements)
            if (element_ == elements(i)) then
                atomic_number = i
                exit
            end if
        end do
    end function atomic_number

end module
