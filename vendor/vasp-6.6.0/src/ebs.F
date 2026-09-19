#include "symbol.inc"
  module ebs
    use prec, only: q
    use constant, only: pi, tpi, felect
    use tutor, only: vtutor
    use pot_struct_def, only: kernel_truncation, coulomb_potential
    use lattice, only: latt
    use mgrid, only: grid_3d, inilgrd, gen_rc_grid, gen_grid, gen_sub_grid_
    use mgrid_struct_def, only: transit
    use pseudo_struct_def, only: potcar
    use poscar_struct_def, only: type_info
    use mpimy, only: communic
    use coulomb_cutoff, only: index_to_gvector, &
                              kernel_for_gvector, &
                              phase_factor_for_index, &
                              generate_padded_lattice, &
                              generate_padded_t_info
    implicit none

    private

    type ewald_parameters
        integer :: spline_points  = 4000
        real(q) :: max_argument   = 4._q   !< Max value for spline
        real(q) :: min_argument   = 0.4_q  !< Min value for spline
    end type

    abstract interface
      function function_at_point(x) result(output)
          import q
          real(q), intent(in) :: x
          real(q) :: output
      end function
    end interface

    type(ewald_parameters) :: ewald_params

    public :: ewald_params, fewald, fewald_lib, fewald_data_writer

    ! these need to be public for the unit-test (ebs_test.pf)
    public :: ewald_parameters, ewald_truncation_parameter, &
              maximum_cells_real_space, maximum_cells_reciprocal_space, &
              erfc_dividedby_x, deriv_erfc_dividedby_x, &
              expsq_dividedby_xsq, deriv_expsq_divided_by_xsq, &
              spline_erfc_dividedby_x, spline_expsq_dividedby_xsq, &
              ewald_real_space, ewald_reciprocal_space, &
              ewald_self_correction, ewald_homogeneous_background, &
              reciprocal_space_two_dimensional_summand

    contains

    subroutine raise_issue_if_zero(x, calling_routine)
        real(q) :: x
        character(len=*) :: calling_routine
#ifdef DEBUG_EWALD
        if (x == 0._q) then
            call vtutor%error("Passing zero to the " // calling_routine &
                            // " routine")
        end if
#endif
    end subroutine

    function erfc_dividedby_x(x) result(output)
        real(q), intent(in) :: x
        real(q) :: output
        call raise_issue_if_zero(x, "erfc_dividedby_x")
        output = erfc(x) / x
    end function

    function deriv_erfc_dividedby_x(x) result(output)
        real(q), intent(in) :: x
        real(q) ::output
        call raise_issue_if_zero(x, "deriv_erfc_dividedby_x")
        output = -2._q * exp(-1._q * x**2._q) / ( x * sqrt(pi) ) - erfc(x) / x**2._q
    end function

    function expsq_dividedby_xsq(x) result(output)
        real(q), intent(in) :: x
        real(q) ::output
        call raise_issue_if_zero(x, "expsq_dividedby_xsq")
        output = exp(-x**2) / x**2._q
    end function

    function deriv_expsq_divided_by_xsq(x) result(output)
        real(q), intent(in) :: x
        real(q) ::output
        call raise_issue_if_zero(x, "deriv_expsq_dividedby_xsq")
        output = -2._q * exp(-x**2._q) * ( 1 + x**2._q) / x**3._q
    end function

    function spline_fit_quantity(quantity_generator, &
                                      derivative_quantity_generator, &
                                      spline_points, max_argument) &
                                 result(spline)
         procedure(function_at_point) :: quantity_generator
         procedure(function_at_point) :: derivative_quantity_generator
         integer :: spline_points
         real(q) :: max_argument
         real(q), allocatable :: spline(:,:)
         real(q) :: derivative_at_first_point
         integer :: i
         !
         allocate(spline(spline_points, 5))
         !
         spline = 0._q
         !
         do i = 1, spline_points
           spline(i,1) = max_argument * i / spline_points
           spline(i,2) = quantity_generator(spline(i,1))
           if (i == 1) derivative_at_first_point = &
               derivative_quantity_generator(spline(i,1))
         end do
         !
         call splcof(spline, spline_points, spline_points, &
                     derivative_at_first_point)
    end function

    function x_to_index(x, spline_points, max_argument) result(idx)
        real(q) :: x
        integer :: spline_points
        real(q) :: max_argument
        real(q) :: idx
        idx = int(x * spline_points / max_argument)
    end function

    function spline_value(x, spline, &
                          quantity_generator, &
                          spline_points,  &
                          max_argument, min_argument, is_gradient)
        real(q) :: x
        real(q), intent(in) :: spline(:,:)
        procedure(function_at_point) :: quantity_generator
        integer :: spline_points
        real(q) :: max_argument, min_argument
        integer :: idx
        real(q) :: rem
        real(q) :: spline_value
        logical :: is_gradient
        !
        if (x <= max_argument .and. x >= min_argument) then
           idx = x_to_index(x, spline_points, max_argument)
           rem = x - spline(idx,1)
           if (is_gradient) then
               spline_value = spline(idx, 3) &
                            + 2._q * rem * spline(idx,4) &
                            + 3._q * rem**2 * spline(idx,5)
           else
               spline_value = spline(idx,2) &
                            + rem * spline(idx,3) &
                            + rem**2 * spline(idx,4) &
                            + rem**3 * spline(idx,5)
           end if
        else if (x < min_argument) then
           spline_value = quantity_generator(x)
        else
           spline_value = 0._q
        end if
    end function

    subroutine check_if_initialized_properly(is_initialized, spline)
        logical :: is_initialized
        real(q), intent(in), allocatable :: spline(:,:)
#ifdef DEBUG_EWALD
        if (is_initialized .and. (.not. allocated(spline))) then
            call vtutor%error("The spline to fit either the real or &
                              &reciprocal space contributions is not &
                              properly allocated.")
         end if
#endif
    end subroutine

    subroutine spline_erfc_dividedby_x(x, ewald_params, &
                                       output, output_gradient, &
                                       initialize)
        real(q), intent(in)  :: x
        real(q), intent(out) :: output
        real(q), intent(out) :: output_gradient
        type(ewald_parameters), intent(in) :: ewald_params
        logical, optional :: initialize
        logical, save :: is_initialized = .false.
        real(q), allocatable, save :: spline(:,:)
        !
        if (present(initialize)) then
        if (initialize) then
            spline = spline_fit_quantity(erfc_dividedby_x, &
                                     deriv_erfc_dividedby_x, &
                                     ewald_params%spline_points, &
                                     ewald_params%max_argument)
            is_initialized = .true.
        end if
        end if
        !
        call check_if_initialized_properly(is_initialized, spline)
        !
        output = spline_value(x, spline, erfc_dividedby_x, &
                              ewald_params%spline_points, &
                              ewald_params%max_argument, &
                              ewald_params%min_argument, &
                              is_gradient=.false.)
        output_gradient = spline_value(x, spline, &
                                       deriv_erfc_dividedby_x, &
                                       ewald_params%spline_points, &
                                       ewald_params%max_argument, &
                                       ewald_params%min_argument, &
                                       is_gradient=.true.)
    end subroutine

    subroutine spline_expsq_dividedby_xsq(x, ewald_params, &
                                          output, output_gradient, &
                                          initialize)
        real(q), intent(in)  :: x
        real(q), intent(out) :: output
        real(q), intent(out) :: output_gradient
        type(ewald_parameters), intent(in) :: ewald_params
        logical, optional :: initialize
        logical, save :: is_initialized = .false.
        real(q), allocatable, save :: spline(:,:)
        !
        if (present(initialize)) then
        if (initialize) then
            spline = spline_fit_quantity(expsq_dividedby_xsq, &
                                     deriv_expsq_divided_by_xsq, &
                                     ewald_params%spline_points, &
                                     ewald_params%max_argument)
            is_initialized = .true.
        end if
        end if
        !
        call check_if_initialized_properly(is_initialized, spline)
        !
        output = spline_value(x, spline, expsq_dividedby_xsq, &
                              ewald_params%spline_points, &
                              ewald_params%max_argument, &
                              ewald_params%min_argument, &
                              is_gradient=.false.)
        output_gradient = spline_value(x, spline, &
                                       deriv_expsq_divided_by_xsq, &
                                       ewald_params%spline_points, &
                                       ewald_params%max_argument, &
                                       ewald_params%min_argument, &
                                       is_gradient=.true.)
    end subroutine

    function ewald_truncation_parameter(latt_cur) result(alpha)
        type(latt), intent(in)              :: latt_cur
        real(q) :: surface_area
        real(q) :: alpha

        alpha = pi / (latt_cur%omega)**(2._q / 3._q)
    end function

    subroutine maximum_cells(argument, norm_vector, &
                              cell1, cell2, cell3)
        real(q), intent(in)  :: argument
        real(q), intent(in)  :: norm_vector(:)
        integer, intent(out) :: cell1, cell2, cell3
        !
        cell1 = ceiling(argument * norm_vector(1))
        cell2 = ceiling(argument * norm_vector(2))
        cell3 = ceiling(argument * norm_vector(3))
#ifdef DEBUG_EWALD
        if (cell1 == 0 .or. cell2 == 0 .or. cell3 == 0) then
            call vtutor%error("Maximum cell is zero in some &
                    direction in the maximum_cells routine &
                    of the Ewald summation. Either your cell &
                    is particularly strange of the alpha &
                    parameter of the Ewald summation is not &
                    correct.")
        end if
#endif
    end subroutine

    subroutine maximum_cells_real_space(ewald_params, latt_cur, &
                                        coulomb_pot, alpha, &
                                        cell1, cell2, cell3)
        type(ewald_parameters), intent(in)  :: ewald_params
        type(latt), intent(in)              :: latt_cur
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        integer, intent(out)                :: cell1, cell2, cell3
        real(q)                             :: max_radial_distance

        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality, &
                  surface_normal_direction => coulomb_pot%kernel_truncate%surface_normal_direction)
        !
        ! The max radial distance in real space is given be the maximal
        ! value that r can taken given by the max_argument parameter in
        ! ewald_params Explicitly, it is max_params ~ sqrt(alpha) *
        ! R_{max}
        !
        max_radial_distance = ewald_params%max_argument / sqrt(alpha)
        !
        call maximum_cells(max_radial_distance, latt_cur%bnorm, &
                           cell1, cell2, cell3)
        !
        if (is_truncated .and. dimensionality == 2) then
            !
            ! The 2D Ewald summation does not include any summation
            ! in the surface normal direction
            !
            if (surface_normal_direction == 1) then
                cell1 = 0
            else if (surface_normal_direction == 2) then
                cell2 = 0
            else if (surface_normal_direction == 3) then
                cell3 = 0
            end if
        end if
        !
        end associate
    end subroutine

    function get_maximum_distance(gridc, latt_cur, t_info, coulomb_pot) &
                                  result(z_max)
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        type(type_info), intent(in)         :: t_info
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q)                             :: z_max
        real(q)                             :: z_ij
        integer                             :: nion_1, nion_2
        integer                             :: number_nodes, istart

#ifdef MPI
        number_nodes = min(t_info%nions, gridc%comm%ncpu)
        istart = gridc%comm%node_me
#else
        number_nodes = 1
        istart = 1
#endif

        associate(surface_normal_direction => coulomb_pot%kernel_truncate%surface_normal_direction)
        !
        z_max = 0.0_q
        !
        do nion_1 = istart, t_info%nions, number_nodes
        do nion_2 = nion_1, t_info%nions
           !
           z_ij = distance_two_atoms_axis(nion_1, nion_2, &
                                          t_info%posion, &
                                          latt_cur%A, &
                                          surface_normal_direction, &
                                          mode="c", &
                                          aperiodic_dimension=&
                                          surface_normal_direction)
           !
           ! only interested in the absolute value of z_ij to get
           ! the required g-vectors for 2D boundary conditions
           !
           if (abs(z_ij) > z_max) z_max = abs(z_ij)
        end do
        end do
        !
        CALLMPI( M_max_d(gridc%comm, z_max ,1))
        !
        end associate
    end function

    subroutine maximum_cells_reciprocal_space(ewald_params, gridc, &
                                              latt_cur, &
                                              t_info, &
                                              coulomb_pot, &
                                              alpha, &
                                              cell1, cell2, cell3)
        type(ewald_parameters), intent(in)  :: ewald_params
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        type(type_info), intent(in)         :: t_info
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        integer, intent(out)                :: cell1, cell2, cell3
        real(q)                             :: max_radial_distance
        real(q)                             :: z_max

        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality)
        !
        if (is_truncated .and. dimensionality == 0) then
           !
           ! There is no reciprocal space contribution for molecular
           ! systems and hence the required cells are also zero
           !
           cell1 = 0
           cell2 = 0
           cell3 = 0
        else if (is_truncated .and. dimensionality == 2) then
           !
           z_max = get_maximum_distance(gridc, latt_cur, t_info, &
                                        coulomb_pot)
           !
           call maximum_cells_2D(coulomb_pot, gridc, latt_cur, z_max, alpha, &
                                 cell1, cell2, cell3)
        else
           !
           ! The max radial distance in reciprocal space is given by the
           ! maximum argument that can be achieved by g / (2sqrt(alpha).
           ! Since there has been no 2pi normalization of g, the number
           ! of cells (n * b) will be given by maximum_argument *
           ! 2sqrt(alpha) / 2pi
           !
           max_radial_distance = ewald_params%max_argument &
                               * sqrt(alpha) / pi
           call maximum_cells(max_radial_distance, latt_cur%anorm, &
                              cell1, cell2, cell3)
        end if
        end associate
    end subroutine

    subroutine maximum_cells_2D(coulomb_pot, gridc, latt_cur, z_max, alpha, cell1, cell2, cell3)
        type(coulomb_potential), intent(in) :: coulomb_pot
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        real(q), intent(in)                 :: z_max
        real(q), intent(in)                 :: alpha
        integer, intent(out)                :: cell1, cell2, cell3
        integer                             :: nmax_a, nmax_b, ncell_a, ncell_b
        integer                             :: idx(3)
        real(q)                             :: g_vector(3), g_norm, largest_g_norm
        real(q)                             :: summand

        associate(surface_normal_direction => coulomb_pot%kernel_truncate%surface_normal_direction)
        !
        if (surface_normal_direction == 1) then
            !
            cell1 = 0
            nmax_a = gridc%ngy; nmax_b = gridc%ngz
        else if (surface_normal_direction == 2) then
            !
            cell2 = 0
            nmax_a = gridc%ngz; nmax_b = gridc%ngx
        else if (surface_normal_direction == 3) then
            !
            cell3 = 0
            nmax_a = gridc%ngx; nmax_b = gridc%ngy
        end if
        !
        ! Figure out the other two in-plane maximum g-vectors
        !
        ncell: do ncell_a = 1, nmax_a
        do ncell_b = 1, nmax_b
           !
           if (surface_normal_direction == 1) then
               !
               idx(:) = [cell1, ncell_a, ncell_b]
           else if (surface_normal_direction == 2) then
               !
               idx(:) = [ncell_b, cell2, ncell_a]
           else if (surface_normal_direction == 3) then
               !
               idx(:) = [ncell_a, ncell_b, cell3]
           end if
           !
           if (all(idx == 0)) cycle
           !
           call index_to_gvector(latt_cur, idx, g_vector)
           !
           g_norm = sqrt( g_vector(1)**2._q &
                        + g_vector(2)**2._q &
                        + g_vector(3)**2._q )
           !
           summand = reciprocal_space_two_dimensional_summand(&
                                              g_norm, alpha, z_max,&
                                              mode="p")
           !
           if (summand / g_norm < 1e-10_q) then
               !
               ! this value is chosen to be commensurate with the 3D
               ! ewald summation cutoffs for numerical precision
               !
               largest_g_norm = g_norm
               !
               exit ncell
           end if
        end do
        end do ncell
        !
        if (surface_normal_direction == 1) then
            !
            cell2 = largest_g_norm / norm2(latt_cur%b(:,2)) + 1
            cell3 = largest_g_norm / norm2(latt_cur%b(:,3)) + 1
        else if (surface_normal_direction == 2) then
            !
            cell3 = largest_g_norm / norm2(latt_cur%b(:,3)) + 1
            cell1 = largest_g_norm / norm2(latt_cur%b(:,1)) + 1
        else if (surface_normal_direction == 3) then
            !
            cell1 = largest_g_norm / norm2(latt_cur%b(:,1)) + 1
            cell2 = largest_g_norm / norm2(latt_cur%b(:,2)) + 1
        end if
        end associate
    end subroutine

    subroutine ewald_real_space(gridc, latt_cur, p, t_info, &
                                ewald_params, coulomb_pot, alpha, &
                                energy, force, stress)
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        type (type_info), intent(in)        :: t_info
        type (potcar), intent(in)           :: p(t_info%ntyp)
        type(ewald_parameters), intent(in)  :: ewald_params
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        real(q), intent(out)                :: energy
        real(q), allocatable, intent(out)   :: force(:,:), stress(:,:)
        real(q), allocatable                :: force_contribs(:,:,:)
        integer                             :: istart, number_nodes
        integer                             :: cell1, cell2, cell3
        integer                             :: n_cell1, n_cell2, n_cell3
        integer                             :: type_1, type_2
        real(q)                             :: double_counting_factor
        real(q), allocatable                :: dposion(:), rdiff(:)
        real(q)                             :: dummy, dummy_gradient
        real(q)                             :: summand, summand_gradient
        real(q) :: z1, z2, x1, x2, x3, xdiff, ydiff, zdiff, r_ij
        integer :: nion_1, nion_2, i, j
        !
        PROFILING_START('ewald_real_space')
        !
        allocate(dposion(3), rdiff(3))
        allocate(force_contribs(3, t_info%nions, t_info%nions))
        allocate(force(3, t_info%nions))
        allocate(stress(3,3))
        !
#ifdef MPI
        number_nodes = min(t_info%nions, gridc%comm%ncpu)
        istart = gridc%comm%node_me
#else
        number_nodes = 1
        istart = 1
#endif
        !
        energy = 0._q
        force_contribs = 0._q
        force = 0._q
        stress = 0._q
        !
        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality, &
                  surface_normal_direction => coulomb_pot%kernel_truncate%surface_normal_direction)
        !
        if (is_truncated .and. dimensionality == 0) then
          !
          ! Just the coulomb energy and force
          !
          do nion_1 = istart, t_info%nions, number_nodes
          do nion_2 = nion_1, t_info%nions
            !
            if (nion_1 == nion_2) cycle
            !
            type_1 = t_info%ityp(nion_1)
            type_2 = t_info%ityp(nion_2)
            !
            z1 = p(type_1)%ZVALF * t_info%vca(type_1)
            z2 = p(type_2)%ZVALF * t_info%vca(type_2)
            !
            dposion = t_info%posion(:,nion_1) - t_info%posion(:,nion_2)
            !
            x1 = modulo(dposion(1)+100.5_q, 1._q)-0.5_q
            x2 = modulo(dposion(2)+100.5_q, 1._q)-0.5_q
            x3 = modulo(dposion(3)+100.5_q, 1._q)-0.5_q
            !
            xdiff = x1 * latt_cur%A(1,1) &
                  + x2 * latt_cur%A(1,2) &
                  + x3 * latt_cur%A(1,3)
            ydiff = x1 * latt_cur%A(2,1) &
                  + x2 * latt_cur%A(2,2) &
                  + x3 * latt_cur%A(2,3)
            zdiff = x1 * latt_cur%A(3,1) &
                  + x2 * latt_cur%A(3,2) &
                  + x3 * latt_cur%A(3,3)
            !
            rdiff = [xdiff, ydiff, zdiff]
            !
            r_ij = sqrt(xdiff**2 + ydiff**2 + zdiff**2)
            !
            ! if r_ij is less than 1E-10 it likely means that VCA
            ! is activated and the energies must not be computed for two
            ! vca atoms
            !
            if (r_ij < 1e-10_q) cycle
            !
            energy = energy + (z1 * z2) / r_ij
            !
            do i = 1, size(rdiff)
                force_contribs(i,nion_2,nion_1) &
                              = force_contribs(i,nion_2,nion_1) &
                              + z1 * z2 * rdiff(i) / r_ij**3
            end do
          end do
          end do
          !
        else
          !
          call maximum_cells_real_space(ewald_params, latt_cur, &
                                        coulomb_pot, alpha, &
                                        cell1, cell2, cell3)
          !
          call spline_erfc_dividedby_x(ewald_params%max_argument, &
                                       ewald_params, dummy, &
                                       dummy_gradient, initialize=.true.)
          !
          do n_cell1 = -cell1, cell1
          do n_cell2 = -cell2, cell2
          do n_cell3 = -cell3, cell3
             !
             do nion_1 = istart, t_info%nions, number_nodes
             do nion_2 = nion_1, t_info%nions
                !
                if (nion_1 == nion_2 .and. n_cell1 == 0 .and. &
                    n_cell2 == 0 .and. n_cell3 == 0) cycle ! Same atom
                !
                type_1 = t_info%ityp(nion_1)
                type_2 = t_info%ityp(nion_2)
                !
                z1 = p(type_1)%ZVALF * t_info%vca(type_1)
                z2 = p(type_2)%ZVALF * t_info%vca(type_2)
                !
                if (nion_2 == nion_1) then
                    double_counting_factor = 0.5_q
                else
                    double_counting_factor = 1.0_q
                end if
                !
                dposion = t_info%posion(:,nion_1) - t_info%posion(:,nion_2)
                !
                x1 = modulo(dposion(1)+100.5_q, 1._q)-0.5_q + n_cell1
                x2 = modulo(dposion(2)+100.5_q, 1._q)-0.5_q + n_cell2
                x3 = modulo(dposion(3)+100.5_q, 1._q)-0.5_q + n_cell3
                !
                if (is_truncated) then
                    !
                    ! if truncated, make sure that the minimum image
                    ! convention is not used along the aperiodic
                    ! dimension
                    !
                    if (surface_normal_direction == 1) then
                        !
                        x1 = dposion(1)
                    else if (surface_normal_direction == 2) then
                        !
                        x2 = dposion(2)
                    else
                        !
                        x3 = dposion(3)
                    end if
                end if
                !
                xdiff = x1 * latt_cur%A(1,1) &
                      + x2 * latt_cur%A(1,2) &
                      + x3 * latt_cur%A(1,3)
                ydiff = x1 * latt_cur%A(2,1) &
                      + x2 * latt_cur%A(2,2) &
                      + x3 * latt_cur%A(2,3)
                zdiff = x1 * latt_cur%A(3,1) &
                      + x2 * latt_cur%A(3,2) &
                      + x3 * latt_cur%A(3,3)
                rdiff = [xdiff, ydiff, zdiff]
                !
                r_ij = sqrt(xdiff**2 + ydiff**2 + zdiff**2)
                !
                ! if r_ij is less than 1E-10 it likely means that VCA
                ! is activated and the energies must not be computed for two
                ! vca atoms
                !
                if (r_ij < 1e-10_q) cycle
                !
                call spline_erfc_dividedby_x(sqrt(alpha) * r_ij, &
                                             ewald_params, summand, &
                                             summand_gradient)
                !
                energy = energy &
                       + double_counting_factor &
                       * z1 * z2 * summand
                !
                do i = 1, size(rdiff)
                    force_contribs(i,nion_2,nion_1) &
                                  = force_contribs(i,nion_2,nion_1) &
                                  - rdiff(i) / r_ij &
                                  * z1 * z2 * summand_gradient
                    do j = 1, size(rdiff)
                        stress(i,j) = stress(i,j) &
                                    - rdiff(i) * rdiff(j) / r_ij &
                                    * double_counting_factor &
                                    * z1 * z2 * summand_gradient
                    end do
                end do
             end do
             end do
          end do
          end do
          end do
        end if
        !
        ! Multiply by sqrt(alpha) to correct for the argument of the
        ! spline fit and felect for unit conversion if not truncated
        !
        CALLMPI( M_sum_d(gridc%comm, energy, 1))
        !
        if (is_truncated .and. dimensionality == 0) then
            !
            energy = energy * felect
        else
            !
            energy = energy * felect * sqrt(alpha)
        end if
        !
        CALLMPI( M_sum_d(gridc%comm, force_contribs, 3*t_info%nions*t_info%nions))
        call midpoint_forces(t_info%nions, force_contribs, force)
        !
        ! Multiply by alpha to correct for the argument of the spline
        ! fine and felect for unit conversion (Angstrom already
        ! accounted for) if not truncated
        !
        if (is_truncated .and. dimensionality == 0) then
            force = force * felect
        else
            force = force * felect * alpha
        end if
        !
        stress(2,1) = stress(1,2)
        stress(3,1) = stress(1,3)
        stress(3,2) = stress(2,3)
        stress = stress * felect * alpha
        !
        CALLMPI( M_sum_d(gridc%comm, stress, 9))
        end associate
        !
        PROFILING_STOP('ewald_real_space')
    end subroutine

    function reciprocal_space_two_dimensional_summand(g_norm, &
                                                alpha, z_ij, mode) &
                                                result(summand)
        real(q), intent(in)          :: g_norm
        real(q), intent(in)          :: alpha
        real(q), intent(in)          :: z_ij
        character(len=1), intent(in) :: mode
        real(q)                      :: summand_1, summand_2
        real(q)                      :: summand
        real(q)                      :: operand_erfc_1, operand_erfc_2

        operand_erfc_1 = alpha**0.5_q * z_ij &
                       + g_norm / (2._q * alpha**0.5_q)
        !
        operand_erfc_2 = -1._q * alpha**0.5_q * z_ij &
                       + g_norm / (2._q * alpha**0.5_q)
        !
        ! summand_1 = exp(g_norm * z_ij) * erfc(operand_erfc_1)
        ! numerically stable summands
        ! using erfc(x) = exp(-x^2) * erfc_scaled(x)
        if (operand_erfc_1 >= 0._q) then
            summand_1 = exp(g_norm*z_ij - operand_erfc_1**2) * erfc_scaled(operand_erfc_1)
        else
            ! use erfc(-x)=2-erfc(x) to keep erfc_scaled argument nonnegative
            summand_1 = 2._q*exp(g_norm*z_ij) &
                      - exp(g_norm*z_ij - operand_erfc_1**2) * erfc_scaled(-operand_erfc_1)
        end if

        ! summand_2 = exp(-1._q * g_norm * z_ij) * erfc(operand_erfc_2)
        if (operand_erfc_2 >= 0._q) then
            summand_2 = exp(-1._q * g_norm*z_ij - operand_erfc_2**2) * erfc_scaled(operand_erfc_2)
        else
            summand_2 = 2._q*exp(-g_norm*z_ij) &
                      - exp(-1._q * g_norm*z_ij - operand_erfc_2**2) * erfc_scaled(-operand_erfc_2)
        end if

#ifdef FULL_GRID_REC_EWALD_2D
        if (operand_erfc_1 > 6._q) then
            !
            summand_1 = 0._q
        else
            !
            summand_1 = exp(g_norm * z_ij) * erfc(operand_erfc_1)
        end if
        !
        if (operand_erfc_2 > 6._q) then
            !
            summand_2 = 0._q
        else
            !
            summand_2 = exp(-1._q * g_norm * z_ij) &
                      * erfc(operand_erfc_2)
        end if
#endif
        !
        if (mode == "p") then
            !
            summand = (summand_1 + summand_2)
        else
            !
            summand = (summand_1 - summand_2)
        end if
    end function

    subroutine reciprocal_space_two_dimensional(gridc, latt_cur, p, &
                                                coulomb_pot, &
                                                ewald_params, &
                                                t_info, alpha, energy, &
                                                force_contribs)
        type(grid_3d), intent(in), target   :: gridc
        type(latt), intent(in)              :: latt_cur
        type (type_info), intent(in)        :: t_info
        type (potcar), intent(in)           :: p(t_info%ntyp)
        type(coulomb_potential), intent(in) :: coulomb_pot
        type(ewald_parameters), intent(in)  :: ewald_params
        real(q), intent(in)                 :: alpha
        real(q), intent(inout)              :: energy
        real(q), intent(inout)              :: force_contribs(:,:,:)
        integer                             :: nc, ni
        integer                             :: number_nodes, istart
        integer                             :: cell1, cell2, cell3
        integer                             :: n_cell1, n_cell2, n_cell3
        integer                             :: nion_1, nion_2
        integer                             :: idx(3)
        real(q), allocatable                :: r_ij(:)
        real(q)                             :: g_vector(3)
        real(q)                             :: g_norm, z_ij
        real(q)                             :: summand, summand_
        real(q)                             :: g_dot_r
        real(q)                             :: surface_area
        integer                             :: type_1, type_2
        real(q)                             :: z1, z2
        real(q)                             :: double_counting_factor
        integer                             :: i
        real(q)                             :: common_force_term
        type(grid_3d), pointer              :: relevant_grid
        type(grid_3d), target               :: trick_grid
        type(transit)                       :: translation_table

        allocate(r_ij(3))
        !
#ifdef MPI
        number_nodes = min(t_info%nions, gridc%comm%ncpu)
        istart = gridc%comm%node_me
#else
        number_nodes = 1
        istart = 1
#endif
        !
        call maximum_cells_reciprocal_space(ewald_params, gridc, &
                                            latt_cur,&
                                            t_info, &
                                            coulomb_pot, &
                                            alpha, cell1, cell2, &
                                            cell3)
        !
        call inilgrd(2*cell1+1, 2*cell2+1, 2*cell3+1, trick_grid)
        !
        trick_grid%comm => gridc%comm
        !
        trick_grid%NGX_rd = 2*cell1+1
        trick_grid%NGY_rd = 2*cell2+1
        trick_grid%NGZ_rd = 2*cell3+1
        trick_grid%lreal = .false.
        trick_grid%real2cplx = .false.
        call gen_sub_grid_(trick_grid, trick_grid, &
                           translation_table, .true., .false.)
        !
        relevant_grid => trick_grid
#ifdef FULL_GRID_REC_EWALD_2D
        relevant_grid => gridc
#endif
        !
        associate(dimensionality => coulomb_pot%kernel_truncate%dimensionality, &
                  surface_normal_direction => coulomb_pot%kernel_truncate%surface_normal_direction)
        !
        energy = 0._q
        !
        do nion_1 = 1, t_info%nions
        do nion_2 = nion_1, t_info%nions
           !
           if (nion_2 == nion_1) then
               !
               double_counting_factor = 0.5_q
           else
               !
               double_counting_factor = 1.0_q
           end if
           !
           type_1 = t_info%ityp(nion_1)
           type_2 = t_info%ityp(nion_2)
           !
           z1 = p(type_1)%ZVALF * t_info%vca(type_1)
           z2 = p(type_2)%ZVALF * t_info%vca(type_2)
           !
           do nc = 1, relevant_grid%rc%ncol
           n_cell2 = relevant_grid%rc%i2(nc)
           n_cell3 = relevant_grid%rc%i3(nc)
           do n_cell1 = 1, relevant_grid%rc%nrow
              !
              idx(:) = [relevant_grid%lpctx(n_cell1), &
                        relevant_grid%lpcty(n_cell2), &
                        relevant_grid%lpctz(n_cell3)]
              !
              ! the g=0 term is handled in homogeneous background charge
              !
              if (all(idx == 0)) cycle
              !
              call index_to_gvector(latt_cur, idx, g_vector)
              !
              if (g_vector(surface_normal_direction) > 1e-10_q) cycle
              !
              ! only interested in reciprocal lattice vectors which are
              ! zero along surface normal
              !
              g_norm = sqrt( g_vector(1)**2._q &
                           + g_vector(2)**2._q &
                           + g_vector(3)**2._q )
              !
              r_ij = distance_vector_two_atoms(nion_1, nion_2, &
                                               t_info%posion, &
                                               latt_cur%A, mode="d", &
                                               aperiodic_dimension=&
                                               surface_normal_direction)
              !
              z_ij = distance_two_atoms_axis(nion_1, nion_2, &
                                             t_info%posion, &
                                             latt_cur%A, &
                                             surface_normal_direction, &
                                             mode="c", &
                                             aperiodic_dimension=&
                                             surface_normal_direction)
              !
              summand = reciprocal_space_two_dimensional_summand(&
                                                g_norm, alpha, z_ij,&
                                                mode="p")
              !
              g_dot_r = tpi * ( idx(1) * r_ij(1) &
                              + idx(2) * r_ij(2) &
                              + idx(3) * r_ij(3) )
              !
              energy = energy + double_counting_factor &
                              * z1 * z2 * cos(g_dot_r) * summand &
                              / g_norm
              !
              do i = 1, 3
                 !
                 if (i == surface_normal_direction) then
                    !
                    summand_ = reciprocal_space_two_dimensional_summand(&
                                                g_norm, alpha, z_ij, &
                                                mode="m")
                    !
                    common_force_term = -1._q * double_counting_factor &
                                      * z1 * z2 * cos(g_dot_r) &
                                      * summand_
                 else
                    !
                    common_force_term = double_counting_factor &
                                      * z1 * z2 * sin(g_dot_r) &
                                      * summand * g_vector(i) / g_norm
                 end if
                 !
                 force_contribs(i,nion_2,nion_1) = &
                                        force_contribs(i,nion_2,nion_1)&
                                      + common_force_term
              end do
           end do
           end do
        end do
        end do
        !
        surface_area = latt_cur%omega &
                     / latt_cur%anorm(surface_normal_direction)
        !
        energy = energy * pi / surface_area * felect
        !
        force_contribs = force_contribs &
                       * pi / surface_area * felect
        !
        end associate
    end subroutine

    subroutine reciprocal_space_three_dimensional(gridc, latt_cur, p, &
                                                coulomb_pot, &
                                                ewald_params, &
                                                t_info, alpha, energy, &
                                                force_contribs, stress)
        type(grid_3d), intent(in), target   :: gridc
        type(latt), intent(in)              :: latt_cur
        type(type_info), intent(in)         :: t_info
        type(potcar), intent(in)            :: p(t_info%ntyp)
        type(coulomb_potential)             :: coulomb_pot
        type(ewald_parameters), intent(in)  :: ewald_params
        real(q), intent(in)                 :: alpha
        real(q), intent(out)                :: energy
        real(q), intent(inout)              :: force_contribs(:,:,:)
        real(q), intent(inout)              :: stress(:,:)
        integer                             :: number_nodes, istart
        integer                             :: cell1, cell2, cell3
        real(q)                             :: dummy, dummy_gradient
        integer                             :: idx(3)
        real(q)                             :: g_vector(3)
        real(q)                             :: g_norm
        integer                             :: type_1, type_2
        complex(q)                          :: phase_factor_1
        complex(q)                          :: phase_factor_2
        real(q)                             :: summand, summand_gradient
        real(q)                             :: kernel
        real(q)                             :: common_force_term
        real(q)                             :: common_stress_term
        real(q)                             :: z1, z2
        real(q), allocatable                :: stress_(:,:)
        real(q)                             :: structure_factor
        real(q)                             :: double_counting_factor
        integer                             :: i, j
        integer                             :: n_cell1,n_cell2,n_cell3
        integer                             :: nion_1, nion_2, nc
        type(grid_3d), pointer              :: relevant_grid
        type(grid_3d), target               :: trick_grid
        type(transit)                       :: translation_table
        real(q)                             :: factor_grid
        complex(q)                          :: phase_norm

        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality, &
                  kernel_truncate => coulomb_pot%kernel_truncate)
        !
        allocate(stress_, mold=stress)
        !
#ifdef MPI
        number_nodes = min(t_info%nions, gridc%comm%ncpu)
        istart = gridc%comm%node_me
#else
        number_nodes = 1
        istart = 1
#endif
        !
        if (is_truncated .and. dimensionality == 2) then
            !
            relevant_grid => gridc
            !
            factor_grid = 2._q
        else
            !
            call maximum_cells_reciprocal_space(ewald_params, gridc, &
                                                latt_cur,&
                                                t_info, &
                                                coulomb_pot, &
                                                alpha, cell1, cell2, &
                                                cell3)
            !
            call inilgrd(2*cell1+1, 2*cell2+1, 2*cell3+1, trick_grid)
            !
            trick_grid%comm => gridc%comm
            !
            trick_grid%NGX_rd = 2*cell1+1
            trick_grid%NGY_rd = 2*cell2+1
            trick_grid%NGZ_rd = 2*cell3+1
            trick_grid%lreal = .false.
            trick_grid%real2cplx = .false.
            call gen_sub_grid_(trick_grid, trick_grid, &
                               translation_table, .true., .false.)
            !
            relevant_grid => trick_grid
            !
            factor_grid = 1._q
        end if
        !
        call spline_expsq_dividedby_xsq(ewald_params%max_argument, &
                                        ewald_params, dummy, &
                                        dummy_gradient, &
                                        initialize=.true.)
        !
        do nc = 1, relevant_grid%rc%ncol
        n_cell2 = relevant_grid%rc%i2(nc)
        n_cell3 = relevant_grid%rc%i3(nc)
        do n_cell1 = 1, relevant_grid%rc%nrow
           !
           idx(:) = [relevant_grid%lpctx(n_cell1), &
                     relevant_grid%lpcty(n_cell2), &
                     relevant_grid%lpctz(n_cell3)]
           !
           if (all(idx == 0)) cycle
           !
           ! g=0 component handled in homogenous background charge
           !
           call index_to_gvector(latt_cur, idx, g_vector)
           !
           g_norm = sqrt( g_vector(1)**2._q &
                        + g_vector(2)**2._q &
                        + g_vector(3)**2._q )
           !
           call kernel_for_gvector(coulomb_pot, g_vector, kernel)
           !
           call spline_expsq_dividedby_xsq(&
                                        g_norm / 2._q / sqrt(alpha), &
                                        ewald_params, summand, &
                                        summand_gradient)
           !
           ! kernel takes care of the 1/g^2
           !
           summand = summand * ( g_norm**2._q ) / (4._q * alpha)
           !
           structure_factor = 0.0_q
           !
           do nion_1 = 1, t_info%nions
           do nion_2 = nion_1, t_info%nions
             !
             if (nion_2 == nion_1) then
                 !
                 double_counting_factor = 0.5_q
             else
                 !
                 double_counting_factor = 1.0_q
             end if
             !
             type_1 = t_info%ityp(nion_1)
             type_2 = t_info%ityp(nion_2)
             !
             z1 = p(type_1)%ZVALF * t_info%vca(type_1)
             z2 = p(type_2)%ZVALF * t_info%vca(type_2)
             !
             ! vca is accounted for in z1 and z2
             !
             call phase_factor_for_index(idx, t_info%posion(:,nion_1), &
                                   1._q, phase_factor_1)
             call phase_factor_for_index(idx, t_info%posion(:,nion_2), &
                                   1._q, phase_factor_2)
             !
             phase_norm = phase_factor_1 * conjg(phase_factor_2)
             !
             energy = energy &
                    + double_counting_factor &
                    * z1 * z2 * kernel * phase_norm &
                    * summand / latt_cur%omega
             !
             common_force_term = -1._q * double_counting_factor &
                             * z1 * z2 * kernel * aimag(phase_norm) &
                             * summand / latt_cur%omega
             do i = 1, 3
                 !
                 force_contribs(i,nion_2,nion_1) = &
                                        force_contribs(i,nion_2,nion_1) &
                                      + g_vector(i) * common_force_term
             end do
             !
             structure_factor = structure_factor &
                              + double_counting_factor * z1 * z2 &
                              * phase_factor_1 * conjg(phase_factor_2)
           end do
           end do
           !
           common_stress_term = -1.0_q / (latt_cur%omega) &
                              * kernel  * structure_factor &
                              * summand
           !
           do i = 1, 3
           do j = i, 3
              stress_(i,j) = 2._q * g_vector(i) * g_vector(j) &
                            / g_norm**2._q &
                           * ( g_norm**2._q / (4._q * alpha) + 1._q)
              if (i==j) stress_(i,j) = stress_(i,j) - 1._q
           end do
           end do
           !
           stress_(2,1) = stress_(1,2)
           stress_(3,1) = stress_(1,3)
           stress_(3,2) = stress_(2,3)
           !
           stress = stress + stress_ * common_stress_term
        end do
        end do
        !
        energy = energy * factor_grid
        force_contribs = force_contribs * factor_grid
        stress = stress * factor_grid
        !
        end associate
    end subroutine

    subroutine ewald_reciprocal_space(gridc, latt_cur, p, t_info, &
                                ewald_params, coulomb_pot, alpha, &
                                energy, force, stress)
        type(grid_3d), intent(in), target   :: gridc
        type(latt), intent(in)              :: latt_cur
        type (type_info), intent(in)        :: t_info
        type (potcar), intent(in)           :: p(t_info%ntyp)
        type(ewald_parameters), intent(in)  :: ewald_params
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        real(q), intent(out)                :: energy
        real(q), allocatable, intent(out)   :: force(:,:), stress(:,:)
        real(q), allocatable                :: force_contribs(:,:,:)

        PROFILING_START('ewald_reciprocal_space')
        !
        allocate(force_contribs(3, t_info%nions, t_info%nions))
        allocate(force(3, t_info%nions))
        allocate(stress(3,3))
        !
        energy = 0.0_q
        force = 0.0_q
        force_contribs = 0.0_q
        stress = 0.0_q
        !
        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality)
        !
        if (is_truncated .and. dimensionality == 0) then
            !
            ! There is no reciprocal space contribution to coulomb
            ! energy for molecular systems
            !
        else if (is_truncated .and. dimensionality == 2) then
            !
            call reciprocal_space_two_dimensional(gridc, latt_cur, p, &
                                                coulomb_pot, &
                                                ewald_params, &
                                                t_info, alpha, energy, &
                                                force_contribs)
        else
            !
            call reciprocal_space_three_dimensional(gridc, latt_cur, p,&
                                                coulomb_pot, &
                                                ewald_params, &
                                                t_info, alpha, energy, &
                                                force_contribs, stress)
        end if
        !
        CALLMPI( M_sum_d(gridc%comm, energy, 1))
        !
        CALLMPI( M_sum_d(gridc%comm, force_contribs, 3*t_info%nions*t_info%nions))
        call midpoint_forces(t_info%nions, force_contribs, force)
        !
        CALLMPI( M_sum_d(gridc%comm, stress, 9))
        !
        end associate
        PROFILING_STOP('ewald_reciprocal_space')
    end subroutine

    pure subroutine midpoint_forces(nions, force_contribs, force)
        integer, intent(in)     :: nions
        real(q), intent(inout)  :: force_contribs(:,:,:)
        real(q), intent(out)    :: force(:,:)
        integer :: i, nion_1, nion_2

        force = 0._q
        !
        do nion_1 = 1, nions
        do nion_2 = 1, nion_1 - 1
           do i = 1, 3
              force_contribs(i,nion_2,nion_1) = &
                                        -force_contribs(i,nion_1,nion_2)
           end do
        end do
        end do
        !
        do nion_1 = 1, nions
        do nion_2 = 1, nions
           do i = 1, 3
              force(i,nion_1) = force(i,nion_1) &
                              + force_contribs(i,nion_2,nion_1)
           end do
        end do
        end do
    end subroutine

    subroutine ewald_self_correction(p, t_info, coulomb_pot, alpha, &
                                    energy)
        type(type_info), intent(in)         :: t_info
        type(potcar), intent(in)            :: p(t_info%ntyp)
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        integer                             :: n_type
        real(q), intent(out)                :: energy

        PROFILING_START('ewald_self_correction')
        !
        energy = 0._q
        !
        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality)
        !
        if (is_truncated .and. dimensionality == 0) then
            !
        else
            do n_type = 1, t_info%ntyp
               !
               energy = energy + p(n_type)%zvalf**2._q &
                      * t_info%nityp(n_type) * t_info%vca(n_type)
            end do
            !
            energy = -1.0_q * sqrt(alpha / pi) * felect * energy
        end if
        !
        end associate
        !
        PROFILING_STOP('ewald_self_correction')
    end subroutine

    pure function homogeneous_two_dimensional_summand(z_ij, alpha) &
                                                    result(summand)
        real(q), intent(in)    :: z_ij
        real(q), intent(in)    :: alpha
        real(q)                :: summand

        summand = z_ij * erf(alpha**0.5_q * z_ij) &
                + exp(-1._q * alpha * z_ij**2._q) / sqrt(alpha * pi)
    end function

    pure function distance_vector_two_atoms(nion_1, nion_2, positions, &
                                            lattice_vectors, mode, &
                                            aperiodic_dimension) &
                                            result(distance_vector)
        integer, intent(in)           :: nion_1, nion_2
        real(q), intent(in)           :: positions(:,:)
        real(q), intent(in)           :: lattice_vectors(:,:)
        character(len=1), intent(in)  :: mode
        real(q), allocatable          :: distance_vector(:)
        integer, intent(in), optional :: aperiodic_dimension
        real(q), allocatable          :: dposion(:)
        real(q)              :: x1, x2, x3, xdiff, ydiff, zdiff

        allocate(dposion(3))
        allocate(distance_vector(3))
        !
        dposion(:) = positions(:,nion_1) - positions(:,nion_2)
        !
        x1 = modulo(dposion(1)+100.5_q, 1._q)-0.5_q
        x2 = modulo(dposion(2)+100.5_q, 1._q)-0.5_q
        x3 = modulo(dposion(3)+100.5_q, 1._q)-0.5_q
		!
        if (present(aperiodic_dimension)) then
            !
            ! the distance along a non-periodic dimension should not
            ! cross pbc and must be reported as the real space distance
            !
            if (aperiodic_dimension == 1) then
                !
                x1 = dposion(1)
            else if (aperiodic_dimension == 2) then
                !
                x2 = dposion(2)
            else if (aperiodic_dimension == 3) then
                !
                x3 = dposion(3)
            end if
        end if
        !
        if (mode == "d") then
            !
            distance_vector = [x1, x2, x3]
        else
            !
            xdiff = x1 * lattice_vectors(1,1) &
                  + x2 * lattice_vectors(1,2) &
                  + x3 * lattice_vectors(1,3)
            ydiff = x1 * lattice_vectors(2,1) &
                  + x2 * lattice_vectors(2,2) &
                  + x3 * lattice_vectors(2,3)
            zdiff = x1 * lattice_vectors(3,1) &
                  + x2 * lattice_vectors(3,2) &
                  + x3 * lattice_vectors(3,3)
            !
            distance_vector = [xdiff, ydiff, zdiff]
        end if
    end function

    pure function distance_two_atoms_axis(nion_1, nion_2, positions, &
                                          lattice_vectors, axis, mode, &
                                          aperiodic_dimension) &
                                          result(distance)
        integer, intent(in)           :: nion_1, nion_2
        real(q), intent(in)           :: positions(:,:)
        real(q), intent(in)           :: lattice_vectors(:,:)
        integer, intent(in)           :: axis
        integer, intent(in), optional :: aperiodic_dimension
        character(len=1), intent(in)  :: mode
        real(q)                       :: distance
        real(q), allocatable          :: dposion(:)
        real(q)                       :: x1, x2, x3, xdiff, ydiff, zdiff
        real(q), allocatable          :: rdiff(:)

        rdiff = distance_vector_two_atoms(nion_1, nion_2, positions, &
                                            lattice_vectors, mode, &
											aperiodic_dimension)
        !
        distance = rdiff(axis)
    end function

    subroutine homogeneous_two_dimensional(gridc, latt_cur, p, t_info, &
                                           coulomb_pot, &
                                           alpha, energy, force)
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        type(type_info), intent(in)         :: t_info
        type(potcar), intent(in)            :: p(t_info%ntyp)
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        real(q), intent(inout)              :: energy
        real(q), intent(inout)              :: force(:,:)
        real(q), allocatable                :: force_contribs(:,:,:)
        integer                             :: nion_1, nion_2
        integer                             :: type_1, type_2
        real(q)                             :: z1, z2, z_ij
        real(q)                             :: surface_area
        real(q)                             :: summand
        real(q)                             :: double_counting_factor

        associate(dimensionality => coulomb_pot%kernel_truncate%dimensionality, &
                  surface_normal_direction => coulomb_pot%kernel_truncate%surface_normal_direction)
        !
        allocate(force_contribs(3, t_info%nions, t_info%nions))
        !
        force_contribs = 0._q
        !
        do nion_1 = 1, t_info%nions
        do nion_2 = nion_1, t_info%nions
          !
          if (nion_2 == nion_1) then
              double_counting_factor = 0.5_q
          else
              double_counting_factor = 1.0_q
          end if
           !
           type_1 = t_info%ityp(nion_1)
           type_2 = t_info%ityp(nion_2)
           !
           z1 = p(type_1)%ZVALF * t_info%vca(type_1)
           z2 = p(type_2)%ZVALF * t_info%vca(type_2)
           !
           z_ij = distance_two_atoms_axis(nion_1, nion_2, &
                                          t_info%posion, &
                                          latt_cur%A, &
                                          surface_normal_direction, &
                                          mode="c", &
                                          aperiodic_dimension=&
                                          surface_normal_direction)
           !
           summand = homogeneous_two_dimensional_summand(z_ij, alpha)
           !
           energy = energy + double_counting_factor * z1 * z2 * summand
           !
           ! There is a force contribution owing to the z-dependence
           ! of this g=0 term, i.e. d/dz of this energy is non-zero
           !
           force_contribs(surface_normal_direction,nion_2,nion_1) = &
              force_contribs(surface_normal_direction,nion_2,nion_1) &
                + double_counting_factor * z1 * z2 &
                  * erf(alpha**0.5_q * z_ij)
        end do
        end do
        !
        surface_area = latt_cur%omega &
                     / latt_cur%anorm(surface_normal_direction)
        !
        energy = -2._q * pi * energy * felect / surface_area
        !
        call midpoint_forces(t_info%nions, force_contribs, force)
        !
        force = 2._q * pi * force * felect / surface_area
        end associate
    end subroutine

    subroutine homogeneous_three_dimensional(gridc, latt_cur, p, &
                                             t_info, coulomb_pot, &
                                             alpha, energy, stress)
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        type(type_info), intent(in)         :: t_info
        type(potcar), intent(in)            :: p(t_info%ntyp)
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        real(q), intent(out)                :: energy
        real(q), intent(inout)              :: stress(:,:)
        real(q)                             :: sum_z_ion
        integer                             :: nion_1, nion_2
        integer                             :: type_1, type_2
        integer                             :: n_type, i

        sum_z_ion = 0._q
        !
        do n_type = 1, t_info%ntyp
           !
           sum_z_ion = sum_z_ion + p(n_type)%zvalf &
                                 * t_info%nityp(n_type) &
                                 * t_info%vca(n_type)
        end do
        !
        energy = -2._q * pi / latt_cur%omega / (4._q * alpha) &
               * sum_z_ion**2._q
        energy = energy * felect
        !
        do i = 1, 3
           stress(i,i) = energy
        end do
    end subroutine


    subroutine ewald_homogeneous_background(gridc, p, t_info, latt_cur,&
                                            coulomb_pot, alpha, &
                                            energy, force, stress)
        type(grid_3d), intent(in)           :: gridc
        type(type_info), intent(in)         :: t_info
        type(potcar), intent(in)            :: p(t_info%ntyp)
        type(latt), intent(in)              :: latt_cur
        type(latt)                          :: latt_cur_padded
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), intent(in)                 :: alpha
        real(q), intent(out)                :: energy
        real(q), allocatable, intent(out)   :: force(:,:)
        real(q), allocatable, intent(out)   :: stress(:,:)

        PROFILING_START('ewald_homogeneous_background')
        !
        allocate(force(3, t_info%nions))
        allocate(stress(3,3))
        !
        energy = 0._q
        force = 0._q
        stress = 0._q
        !
        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality)
        !
        if (is_truncated .and. dimensionality == 0) then
            !
            ! No homogeneous background charge contribution for coulomb
            ! interaction in molecular systems
            !
        else if (is_truncated .and. dimensionality == 2) then
            !
            call homogeneous_two_dimensional(gridc, latt_cur, p, &
                                             t_info, coulomb_pot, &
                                             alpha, energy, force)
        else
            !
            call homogeneous_three_dimensional(gridc, latt_cur, p, &
                                             t_info, coulomb_pot, &
                                             alpha, energy, stress)
        end if
        !
        end associate
        !
        PROFILING_STOP('ewald_homogeneous_background')
    end subroutine

!
!> old way of calling fewald which was as a library not using and of
!> the derived types used by VASP, this routine is just a wrapper
!
    subroutine fewald_lib(posion,ewifor,a,b,anorm,bnorm, &
                          omega,ewsif,tewen,ntyp,zval,vca, &
                          nions,niond,ityp,nityp,comm)
        real(q), target ::  posion(3,nions)    !< position in direct coordinates
        real(q) ::  ewifor(3,nions)    !< forces in ions
        real(q) ::  a(3,3)             !< direct lattice vectors
        real(q) ::  b(3,3)             !< reciprocal lattice vectors
        real(q) ::  anorm(3)           !< length of direct lattice vectors
        real(q) ::  bnorm(3)           !< length of reciprocal lattice vectors
        real(q) ::  omega
        real(q) ::  ewsif(3,3)         !< stress tensor
        real(q) ::  tewen              !< ewald energy
        integer ::  ntyp, nions, niond !< number of types and ions
        real(q), target ::  zval(ntyp)         !< ionic charge of each species
        real(q), target ::  vca(ntyp)          !< vca weight of each species
        integer, target ::  nityp(ntyp),ityp(nions)
        integer ::  iu
        type(communic), intent(in), target :: comm
        ! derived types used by the current fewald routine
        type(grid_3d) :: gridc
        type(latt)    :: latt_cur
        type (type_info) :: t_info
        type (potcar), allocatable :: p(:)
        type(coulomb_potential) :: coulomb_pot

        latt_cur%a = a
        latt_cur%b = b
        latt_cur%anorm = anorm
        latt_cur%bnorm = bnorm
        latt_cur%omega = omega
        !
        t_info%posion => posion
        t_info%ntyp = ntyp
        t_info%vca => vca
        t_info%nions = nions
        t_info%ityp => ityp
        t_info%nityp => nityp
        !
        allocate(p(t_info%ntyp))
        p(:)%zvalf = zval(:)
        !
        gridc%comm => comm
        !
        call fewald(gridc, latt_cur, p, t_info, coulomb_pot, &
                    tewen, ewifor, ewsif)
        !
        nullify(t_info%vca)
        nullify(t_info%ityp)
        nullify(t_info%nityp)
    end subroutine

!************************ subroutine fewald ****************************
! rewritten by Georg  Kresse (original version Mike Payne?)
! update 14.03.2012: Dario Alfe supplied a version that
!        performs parallelization over ions
!        additional cleanup and implicit removed by gK
! update 02.06.2024:
!        1. Split fewald into subroutines
!        2. Incorporate the coulomb kernel parameters for systems
!           with lower dimensionality
!        3. Takes in derived types instead of reals and complexes
!        4. Exposes the ewald_cutoff_parameter (alpha)
!        5. Encapsulte defaults into ewald_parameters type
!        6. Unit tests in ebs_test.pf
!
!> Calculates the ewald energy due to Coulomb energy
!> between the ions and the neutralising background.
!> Arbitrary cellshape and arbitrary atomic species are allowed.
!> In addition the subroutine calculates forces and stress corresponding
!> to the ewald energy.
!
! MIND: that this routine assumes that all positions are within
! [0,1], if this is not the case the real space summation might not
! be done for sufficient cells
!
!***********************************************************************
      subroutine fewald(gridc, latt_cur, p, t_info, coulomb_pot, &
                        ewald_energy, ewald_force, ewald_stress, &
                        ewald_cutoff_parameter)
        type(grid_3d), intent(in)           :: gridc
        type(latt), intent(in)              :: latt_cur
        type (type_info), intent(in)        :: t_info
        type (potcar), intent(in)           :: p(t_info%ntyp)
        type(coulomb_potential), intent(in) :: coulomb_pot
        real(q), optional                   :: ewald_cutoff_parameter
        real(q)                :: alpha
        real(q)                :: reciprocal_energy
        real(q), allocatable   :: reciprocal_force(:,:)
        real(q), allocatable   :: reciprocal_stress(:,:)
        real(q)                :: real_energy
        real(q), allocatable   :: real_force(:,:)
        real(q), allocatable   :: real_stress(:,:)
        real(q)                :: self_correction_energy
        real(q)                :: homogeneous_background_energy
        real(q), allocatable   :: homogeneous_background_force(:,:)
        real(q), allocatable   :: homogeneous_background_stress(:,:)
        real(q), intent(out)   :: ewald_energy
        real(q), intent(out)   :: ewald_force(:,:)
        real(q), intent(out)   :: ewald_stress(:,:)

        ewald_energy = 0._q
        ewald_force = 0._q
        ewald_stress = 0._q
        !
        if (present(ewald_cutoff_parameter)) then
            !
            alpha = ewald_cutoff_parameter
        else
            !
            alpha = ewald_truncation_parameter(latt_cur)
        end if
        !
        call ewald_reciprocal_space(gridc, latt_cur, p, t_info, &
                                    ewald_params, &
                                    coulomb_pot, &
                                    alpha, reciprocal_energy, &
                                    reciprocal_force, &
                                    reciprocal_stress)
        !
        call ewald_real_space(gridc, latt_cur, p, t_info, ewald_params,&
                              coulomb_pot, &
                              alpha, real_energy, real_force, &
                              real_stress)
        !
        call ewald_self_correction(p, t_info, &
                                   coulomb_pot, &
                                   alpha, self_correction_energy)
        !
        call ewald_homogeneous_background(gridc, p, t_info, latt_cur, &
                                         coulomb_pot, &
                                         alpha, &
                                         homogeneous_background_energy,&
                                         homogeneous_background_force, &
                                         homogeneous_background_stress)
        !
        ewald_energy = reciprocal_energy &
                     + real_energy &
                     + self_correction_energy &
                     + homogeneous_background_energy
        !
        ewald_force = reciprocal_force &
                    + real_force &
                    + homogeneous_background_force
        !
        ewald_stress = reciprocal_stress &
                     + real_stress &
                     + homogeneous_background_stress
      end subroutine

      subroutine fewald_data_writer(gridc, latt_cur, t_info, coulomb_pot, write_unit, &
                                    ewald_cutoff_parameter)
        type(grid_3D), intent(in)           :: gridc
        type(latt), intent(in), target      :: latt_cur
        type(type_info), intent(in)         :: t_info
        type(coulomb_potential), intent(in) :: coulomb_pot
        integer, intent(in)                 :: write_unit
        real(q), optional                   :: ewald_cutoff_parameter
        real(q)                             :: alpha
        integer                             :: rcell1, rcell2, rcell3
        integer                             :: gcell1, gcell2, gcell3
        type(latt), pointer              :: latt_cur_
        type(latt), target               :: latt_cur_padded

        associate(is_truncated => coulomb_pot%kernel_truncate%active, &
                  dimensionality => coulomb_pot%kernel_truncate%dimensionality, &
                  fine_grid => coulomb_pot%kernel_truncate%fine_grid)
        !
        if (present(ewald_cutoff_parameter)) then
            !
            alpha = ewald_cutoff_parameter
        else
            !
            alpha = ewald_truncation_parameter(latt_cur)
        end if
        !
        if (is_truncated) then
            !
            call generate_padded_lattice(coulomb_pot%kernel_truncate, &
                                         latt_cur, &
                                         latt_cur_padded)
            !
            latt_cur_ => latt_cur_padded
        else
            !
            latt_cur_ => latt_cur
        end if
        !
        call maximum_cells_real_space(ewald_params, latt_cur_, &
                                      coulomb_pot, &
                                      alpha, &
                                      rcell1, rcell2, rcell3)
        !
        if (is_truncated) then
            !
            if (dimensionality == 0) then
                !
                gcell1 = 0
                gcell2 = 0
                gcell3 = 0
            else if (dimensionality == 2) then
                !
                call maximum_cells_reciprocal_space(ewald_params, &
                                                    gridc, &
                                                    latt_cur, &
                                                    t_info, &
                                                    coulomb_pot, &
                                                    alpha, &
                                                    gcell1, gcell2, gcell3)
                !
#ifdef FULL_GRID_REC_EWALD_2D
                gcell1 = fine_grid%NGX
                gcell2 = fine_grid%NGY
                gcell3 = fine_grid%NGZ
#endif
            end if
        else
            !
            call maximum_cells_reciprocal_space(ewald_params, &
                                                gridc, &
                                                latt_cur, &
                                                t_info, &
                                                coulomb_pot, &
                                                alpha, &
                                                gcell1, gcell2, gcell3)
        end if
        !
        if (write_unit >= 0) then
            write(write_unit,5) sqrt(alpha), rcell1, rcell2, rcell3, &
                                             gcell1, gcell2, gcell3
        end if
5       FORMAT(' First call to EWALD:  gamma=',F8.3/ &
  &         ' Maximum number of real-space cells',I4,'x',I4,'x',I4/ &
  &         ' Maximum number of reciprocal cells',I4,'x',I4,'x',I4/)
        end associate
      end subroutine
  end module
