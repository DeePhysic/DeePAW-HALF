module plugins_interface

    use iso_c_binding

    implicit none

    enum, bind(c)
        enumerator :: NO_ERROR, BUG_IN_CPP, IMPORT_ERROR, PYTHON_EXCEPTION
    end enum

    enum, bind(c)
        enumerator :: FORCE_AND_STRESS, MACHINE_LEARNING
    end enum

    type, bind(c) :: neighbor_list_data
        integer(c_int) size
        type(c_ptr) neighbors  ! size
        type(c_ptr) distances  ! size
        type(c_ptr) directions ! 3 x size
    end type neighbor_list_data

#ifndef NO_FORCE_AND_STRESS
    type, bind(c) :: constants_force_and_stress
        real(c_double) ENCUT
        real(c_double) NELECT
        integer(c_int) shape_grid(3)
        integer(c_int) number_ions
        integer(c_int) number_ion_types
        type(c_ptr) ion_types  ! number_ions integers
        type(c_ptr) atomic_numbers  ! number_ion_types integers
        type(c_ptr) lattice_vectors  ! 3 x 3
        type(c_ptr) positions  ! 3 x number_ions reals
        type(c_ptr) ZVAL  ! number_ion_types reals
        type(c_ptr) POMASS  ! number_ion_types reals
        type(c_ptr) forces  ! 3 x number_ions reals
        type(c_ptr) stress  ! 3 x 3
        type(c_ptr) charge_density  ! shape_grid reals/complexes
        integer(c_int) mode  ! force_and_stress or machine_learning
        type(c_ptr) :: neighbor_list  ! number_ions, pointer to neighbor_list_data, only set if requested
    end type constants_force_and_stress

    type, bind(c) :: additions_force_and_stress
        real(c_double) total_energy
        type(c_ptr) forces  ! 3 x number_ions reals
        type(c_ptr) stress  ! 3 x 3
    end type additions_force_and_stress
#endif

#ifndef NO_LOCAL_POTENTIAL
    type, bind(c) :: constants_local_potential
        real(c_double) ENCUT
        real(c_double) NELECT
        integer(c_int) shape_grid(3)
        integer(c_int) number_ions
        integer(c_int) number_ion_types
        type(c_ptr) ion_types  ! number_ions integers
        type(c_ptr) atomic_numbers  ! number_ion_types integers
        type(c_ptr) lattice_vectors  ! 3 x 3
        type(c_ptr) positions  ! 3 x number_ions reals
        type(c_ptr) ZVAL  ! number_ion_types reals
        type(c_ptr) charge_density  ! shape_grid reals/complexes
        type(c_ptr) hartree_potential ! shape_grid reals/complexes
        type(c_ptr) ion_potential ! shape_grid reals/complexes
        integer(c_int) dipole_moment_size
        type(c_ptr) dipole_moment ! 3 / 0 depending on if computed
        type(c_ptr) :: neighbor_list  ! number_ions, pointer to neighbor_list_data, only set if requested
    end type constants_local_potential

    type, bind(c) :: additions_local_potential
        real(c_double) total_energy
        type(c_ptr) :: potential ! shape_grid reals/complexes
    end type additions_local_potential
#endif

#ifndef NO_OCCUPANCIES
    type, bind(c) :: constants_occupancies
        real(c_double) :: NELECT
        real(c_double) :: EFERMI
        real(c_double) :: SIGMA
        integer(c_int) :: ISMEAR
        real(c_double) :: EMIN
        real(c_double) :: EMAX
        real(c_double) :: NUPDOWN
    end type constants_occupancies

    type, bind(c) :: additions_occupancies
        real(c_double) :: NELECT
        real(c_double) :: EFERMI
        real(c_double) :: SIGMA
        integer(c_int) :: ISMEAR
        real(c_double) :: EMIN
        real(c_double) :: EMAX
        real(c_double) :: NUPDOWN
    end type additions_occupancies
#endif

#ifndef NO_STRUCTURE
    type, bind(c) :: constants_structure
        integer(c_int) number_ions
        integer(c_int) number_ion_types
        type(c_ptr) ion_types  ! number_ions integers
        type(c_ptr) atomic_numbers  ! number_ion_types integers
        type(c_ptr) lattice_vectors  ! 3 x 3
        type(c_ptr) positions  ! 3 x number_ions reals
        type(c_ptr) POMASS  ! number_ion_types reals
        real(c_double) total_energy
        type(c_ptr) forces  ! 3 x number_ions reals
        type(c_ptr) stress  ! 3 x 3
        integer(c_int) shape_grid(3)
        type(c_ptr) charge_density ! shape_grid reals/complexes
        type(c_ptr) :: neighbor_list  ! number_ions, pointer to neighbor_list_data, only set if requested
    end type constants_structure

    type, bind(c) :: additions_structure
        type(c_ptr) :: lattice_vectors
        type(c_ptr) :: positions
    end type additions_structure
#endif

    interface
        integer(c_int) function plugins_initialize() bind(c, name="plugins_initialize")
            import c_int
        end function plugins_initialize

        integer (c_int) function plugins_finalize() bind(c, name="plugins_finalize")
            import c_int
        end function plugins_finalize

#ifndef NO_FORCE_AND_STRESS
        integer(c_int) function interface_force_and_stress(constants, additions) &
                bind(c, name="interface_force_and_stress")
            import c_int, constants_force_and_stress, additions_force_and_stress
            type(constants_force_and_stress), intent(in) :: constants
            type(additions_force_and_stress), intent(inout) :: additions
        end function interface_force_and_stress
#endif

#ifndef NO_LOCAL_POTENTIAL
        integer(c_int) function interface_local_potential(constants, &
                                                          additions) &
                bind(c, name="interface_local_potential")
            import c_int, constants_local_potential, additions_local_potential
            type(constants_local_potential), intent(in) :: constants
            type(additions_local_potential), intent(inout) :: additions
        end function interface_local_potential
#endif

#ifndef NO_STRUCTURE
        integer(c_int) function interface_update_structure(constants, &
                                                           additions) &
               bind(c, name="interface_update_structure")
           import c_int, constants_structure, additions_structure
           type(constants_structure), intent(in) :: constants
           type(additions_structure), intent(inout) :: additions
        end function interface_update_structure
#endif

#ifndef NO_OCCUPANCIES
        integer(c_int) function interface_occupancies(constants, &
                                                      additions) &
              bind(c, name="interface_occupancies")
           import c_int, constants_occupancies, additions_occupancies
           type(constants_occupancies), intent(in) :: constants
           type(additions_occupancies), intent(inout) :: additions
        end function interface_occupancies
#endif

    end interface

end module plugins_interface
