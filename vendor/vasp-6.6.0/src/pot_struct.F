#include "symbol.inc"

module pot_struct_def
    use prec, only: q
    use mgrid_struct_def, only: grid_3d, transit

    implicit none

    type kernel_truncation
        logical :: active                   = .false. !< Decides if the Coulomb truncation method is used
        integer :: dimensionality           = 3       !< Dimensionality of the system.
                                                      !< 0 -> atoms and molecules
                                                      !< 1 -> nanorods and wires
                                                      !< 2 -> surfaces, slabs and 2D materials
                                                      !< 3 -> bulk crystals
        integer :: surface_normal_direction = -1      !< Direction in which the surface normal is oriented
        integer :: padding_factor           = huge(1) !< Number of cells that will be padded in any direction
        real(q) :: truncation_factor        = huge(1._q)   !< Truncate the Coulomb kernel with this factor
        real(q) :: truncation_length        = -1._q   !< Truncation length based on the truncation factor
        logical :: coarsen_before_pad       = .false. !< Coarsen the grid before padding it with the zeros in real space.
        logical :: use_simplified_2d_kernel = .false. !< if r_c = L_z/2, use the simplified version of the kernel
        type(grid_3d), pointer :: fine_grid         => NULL() !< The fine grid used in the padding procedure
        type(grid_3d), pointer :: coarse_grid       => NULL() !< The coarse grid used in the padding procedure
        type(transit), pointer :: translation_table => NULL() !< Translation table between fine and coarse grids
    end type

    type coulomb_potential
        type(kernel_truncation) :: kernel_truncate    !< Parameters for truncating the kernel
    end type

end module pot_struct_def
