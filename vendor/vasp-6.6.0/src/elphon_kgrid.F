#include "symbol.inc"
!> This module implements symmetry handling using electron-phonon matrix elements
!> Formatter is applied:
!> fprettify -i 4 src/elphon_kgrid.F
module elphon_kgrid
    use prec
    use iso_c_binding, only: c_int, c_long, c_double, c_char
    use wave_struct_def, ONLY: wavespin, wavedes
    implicit none
    private
    integer, parameter :: bz_grid_points_estimation = 1000
    public elph_kgrid, elph_triplets
#ifdef ELPH_KSPCLIB
    public elphon_kgrid_get_dataset
    public elphon_kgrid_show_dataset
    public elphon_kgrid_get_points
    public elphon_kgrid_get_fermi
    public elphon_kgrid_get_dos
    public elphon_kgrid_get_dos_full
    public elphon_kgrid_get_spectral
    public elphon_kgrid_get_spectral_full
    public elphon_kgrid_get_triplets_at_q
    public elphon_kgrid_get_triplet_integration_weights
#endif

    !>
    !> @brief Data structure to handle reciprocal generalized regular grid.
    !>
    !> Smith normal form is given by D = P grid_matrix Q.
    !>
    !> In the following documentation, GR grid and BZ gird means generalized
    !> regular grid and Brillouin zone grid, respectively.
    !>
    !> One of BZ grid point index can be retrieved by
    !>
    !>     bz_grid_index = bz_grid_map(gr_grid_index + 1)
    !>
    !> Ir-grid point index is obtained from GR-grid point index by
    !>
    !>     ir_grid_index = ir_grid_indices(grg2irgrg(gr_grid_index + 1))
    !>
    !>
    !> D_diag(3) : c_long
    !>     Diagonal elements of D in Smith normal form, which correspond to mesh
    !>     numbers in Q-transformed reciprocal primitive lattice.
    !> P(3, 3) : c_long
    !>     P matrix in Smith normal form. Transposed in fortran.
    !> Q(3, 3) : c_long
    !>     Q matrix in Smith normal form. Transposed in fortran.
    !> QDinv(3, 3) : c_double
    !>     QD^{-1}. Transposed in fortran.
    !> grid_matrix(3, 3) : c_long
    !>     Grid generating matrix with respect to reciprocal primitive basis
    !>     vectors. Transposed in fortran.
    !> regular_rotations(3, 3, 48) : c_long
    !>     Rotation matrices with respect to reciprocal primitive basis
    !>     vectors. q'=qR in fortran.
    !> rotations(3, 3, 48) : c_long
    !>     QD^-1 transformed rotation matrices with respect to transformed
    !>     reciprocal primitive basis vectors. q'=qR in fortran.
    !> PS(3) : c_long
    !>     Shift in transformed double grid.
    !> sampling_length : c_double
    !>     Like that of KPOINTS auto R_k (Angst).
    !> reciprocal_lattice(3, 3) : c_double
    !>     Basis vectors of reciprocal primitive cell (Angst^-1). Without 2pi,
    !>     row vectors in fortran.
    !> microzone_lattice(3, 3) : c_double
    !>     Basis vectors of microzone cell (Angst^-1). Without 2pi,
    !>     row vectors in fortran.
    !> n_rotations : c_long
    !>     Number of point group operations. May or may not with time reversal
    !>     symmetry.
    !> n_ir_grid_points : c_long
    !>     Number of irreducible grid points.
    !> n_bz_grid_points : c_long
    !>     Number of grid points in first BZ including its surface. This can be
    !>     larger than product(D_diag).
    !> grg2irgrg(product(D_diag)) : c_long
    !>     Mapping table of GR-grid point index to index in ir_grid_indices.
    !>     This array is created in this fortran module.
    !>     This array is used only with ir_grid_indices.
    !>     Index of ir_grid_indices stored starts with one, i.e., fortran style.
    !> ir_grid_indices(size(unique(grg2irgrg))) : c_long
    !>     List of unique GR-grid indices.
    !>     This array is created in this fortran module.
    !>     Note that grid point index stored starts with zero.
    !> ir_grid_weights(size(unique(grg2irgrg))) :  c_long
    !>     Counts of unique grid indices in grg2irgrg. Sum = product(D_diag).
    !>     This array is created in this fortran module.
    !>     This array is used only with ir_grid_indices.
    !> bz_grid_addresses(3, n_bz_grid_points) :  c_long
    !>     Shortest grid points from Gamma point among translationally equivalent
    !>     points in reciprocal lattice. All grid points having distances are
    !>     stored.
    !> bz_grid_map(product(D_diag) + 1) :  c_long
    !>     List of first BZ-grid indices corresponding to GZ-grid points among
    !>     translationally equivalent points.
    !>     Note that grid point index stored starts with zero.
    !> bzg2grg(n_bz_grid_points) :  c_long
    !>     Mapping table of BZ-grid point index to GR-grid point index.
    !>     Note that grid point index stored starts with zero.
    !> relative_grid_addresses(3, 4, 24) :  c_long
    !>     Relative vectors of neighboring grid points defined to use tetrahedron
    !>     method on generalized regular grid.
    !> transformation_matrix(3, 3) :  c_long
    !>     Transpose of transformation matrix (T) from standardized conventional
    !>     unit cell basis vectors to input cell basis vectors.
    !>     It is assumed that the input cell is a primitive cell.
    !>     (a_c, b_c, c_c) T = (a_p, b_p, c_p)
    !> spacegroup_number :  c_long
    !>     Space group number (from 1 to 230)
    !> direct_rotations(3, 3, :) :  c_long
    !>     Rotation matrices in direct space with respect to input cell basis
    !>     vectors. The input cell is expected to be a primitive cell and
    !>     number of rotations is equal to or smaller than 48.
    !> n_direct_rotations :  c_long
    !>     Number of direct rotations.
    !>
    !> by Atsushi Togo
    !>
    type elph_kgrid
        integer(c_long) :: D_diag(3)
        integer(c_long) :: P(3, 3)
        integer(c_long) :: Q(3, 3)
        real(c_double) :: QDinv(3, 3)
        integer(c_long) :: grid_matrix(3, 3)
        integer(c_long) :: regular_rotations(3, 3, 48)
        integer(c_long) :: rotations(3, 3, 48)
        integer(c_long) :: PS(3)
        real(c_double) :: sampling_length
        real(c_double) :: reciprocal_lattice(3, 3)
        real(c_double) :: microzone_lattice(3, 3)
        integer(c_long) :: n_rotations
        integer(c_long) :: n_ir_grid_points
        integer(c_long) :: n_bz_grid_points
        integer(c_long), allocatable :: grg2irgrg(:)
        integer(c_long), allocatable :: ir_grid_indices(:)
        integer(c_long), allocatable :: ir_grid_weights(:)
        integer(c_long), allocatable :: bz_grid_addresses(:, :)
        integer(c_long), allocatable :: bz_grid_map(:)
        integer(c_long), allocatable :: bzg2grg(:)
        integer(c_long) :: relative_grid_addresses(3, 4, 24)
        real(c_double) :: transformation_matrix(3, 3)
        integer(c_long) :: spacegroup_number
        integer(c_long), allocatable :: direct_rotations(:, :, :)
        integer(c_long) :: n_direct_rotations
    end type elph_kgrid

    !>
    !> @brief Data structure for (q, k, k') triplets
    !>
    !> ir_triplets(3, num_ir_triplets) :
    !>     Irreducible triplets at specified q in BZ-grid.
    !> ir_weights(num_ir_triplets) :
    !>     Integer weights of irreducible triplets at specified q.
    !> num_ir_triplets
    !>     Number of irreducible triplets at specified q.
    !> rotation_indices
    !>     Indices of rotations that map k to irreducible k
    !>     The values are in {1..kgrid%n_rotations}.
    !> rotated_triplets(3, num_ir_triplets) :
    !>     {(Rq, Rk, Rk')} where Rk is in the list of kgrid%ir_grid_indices.
    type elph_triplets
        integer(c_long), allocatable :: ir_triplets(:, :)
        integer(c_long), allocatable :: ir_weights(:)
        integer(c_long) :: num_ir_triplets
        integer(c_long), allocatable :: rotation_indices(:)
        integer(c_long), allocatable :: rotated_triplets(:, :)
    end type elph_triplets

#ifdef ELPH_KSPCLIB
contains

    !>
    !> @brief A public function to return kgrid dataset
    !>
    !> Parameters
    !> ----------
    !> kgrid : elph_kgrid (inout)
    !>     kgrid dataset.
    !> lattice : real array
    !>     Basis vectors (Angst.) in direct space given as column vectors.
    !>     shape=(3, 3)
    !> positions : real array
    !>     Fractional coordinates of atomic points.
    !>     shape=(3, natoms)
    !> ntypes : integer
    !>     Number of atomic species.
    !> nitypes : integer array
    !>     Numbers of atoms of respective atomic species.
    !>     shape=(ntypes, )
    !> natoms : interger
    !>     Number of atoms.
    !> kspacing : real
    !>     Distance mearure same as VASP INCAR KSPACING
    !>
    !> Returns
    !> -------
    !> succeeded : logical
    !>     Whether grid generation is succeeded or not.
    !>
    !> by Atsushi Togo
    !>
    function elphon_kgrid_get_dataset(kgrid, lattice, positions, ntypes, &
                                      nitypes, natoms, kspacing) result(succeeded)
        use constant

        type(elph_kgrid), intent(inout) :: kgrid
        integer, intent(in) :: ntypes, natoms
        integer, intent(in) :: nitypes(ntypes)
        real(q), intent(in) :: lattice(3, 3)
        real(q), intent(in) :: positions(3, natoms)
        real(q), intent(in) :: kspacing
        logical :: succeeded

        integer :: int_succeeded

#ifdef use_spglib
        call get_symmetry_info_from_spglib(kgrid, lattice, positions, ntypes, &
                                           nitypes, natoms)
#else
        call get_symmetry_info_from_hdf5(kgrid)
#endif
        kgrid%PS(:) = 0
        kgrid%sampling_length = TPI/kspacing
        int_succeeded = get_grid_matrices(kgrid, lattice)
        if (int_succeeded /= 0) then
            int_succeeded = get_reciprocal_rotations(kgrid)
            if (int_succeeded /= 0) then
                call get_relative_grid_addresses(kgrid)
                call get_bz_grid_addresses(kgrid)
                call get_grg2irgrg(kgrid)
            end if
        end if

        if (int_succeeded /= 0) then
            succeeded = .true.
        else
            succeeded = .false.
        end if

    end function elphon_kgrid_get_dataset

    !>
    !> @brief A public function to return a set of irreducible k-points
    !>
    !> Parameters
    !> ----------
    !> kpoints_inter : kpoints_struct (inout)
    !>     k-points data structure. The following members are modified:
    !>     nkpts : integer
    !>         Number of irreducible k-points
    !>     vkpt : real array
    !>         Irreducible k-points in fractional coordiates with respect to
    !>         reciprocal basis vectors.
    !>         shape=(3, nkpts)
    !>     wtkpt : real array
    !>         Weights of irreducible k-points. The sum is one.
    !>         shape=(nkpts, )
    !>     nkpx, nkpy, nkpz : integer
    !>         nkpx * nkpy * nkpz will give number of grid points.
    !> kgrid : elph_kgrid
    !>     kgrid dataset.
    !>
    !> by Atsushi Togo
    !>
    function elphon_kgrid_get_points(kpoints_inter, kgrid) result(succeeded)
        use mkpoints, only: kpoints_struct, tetirr
        use m_tetra_hash, only: tetirr_hash
        use gridsysf, only: gridsys_get_grid_address_from_index
        use mathtools, only: matrix_inverse_3x3, matrix_inverse_3x3i

        type(kpoints_struct), intent(inout) :: kpoints_inter
        type(elph_kgrid), intent(in) :: kgrid
        logical :: succeeded

        integer :: i, nkpts_full
        integer(c_long) :: bz_adrs(3)
        real(q) :: bk(3, 3), ak(3, 3)
        integer, allocatable:: idtet(:, :)

        succeeded = .false.

        !
        ! Grid and irreducible points
        !
        kpoints_inter%nkpts = kgrid%n_ir_grid_points
        if (associated(kpoints_inter%vkpt)) then
            deallocate (kpoints_inter%vkpt)
            nullify (kpoints_inter%vkpt)
            allocate (kpoints_inter%vkpt(3, kgrid%n_ir_grid_points))
        end if
        if (associated(kpoints_inter%wtkpt)) then
            deallocate (kpoints_inter%wtkpt)
            nullify (kpoints_inter%wtkpt)
            allocate (kpoints_inter%wtkpt(kgrid%n_ir_grid_points))
        end if
        nkpts_full = product(kgrid%D_diag)
        kpoints_inter%nkdim = nkpts_full
        ! ikpt maps GR-grid address to ir-grid-point-index.
        ! Assuming i1 + i2 * n1 + i3 * n1 * n2 + 1 where
        ! i1 in (0..n1), i2 in (0..n), i3 in (0..n3).
        ! gridsys_get_grid_address_from_index returns (i1, i2, i3) for
        ! the input i1 + i2 * n1 + i3 * n1 * n2.
        if (allocated(kpoints_inter%fbz2ibz)) then
            deallocate (kpoints_inter%fbz2ibz)
        end if
        allocate (kpoints_inter%fbz2ibz(2,kgrid%D_diag(1)*kgrid%D_diag(2)*kgrid%D_diag(3)))
        do i = 1, nkpts_full
            kpoints_inter%fbz2ibz(1,i) = kgrid%grg2irgrg(i)
        end do
        ! Reduced coordinates of ir-kpoints are obtained from kgrid%ir_grid_indices,
        ! which is consistent with ikpt.
        do i = 1, kgrid%n_ir_grid_points
            bz_adrs(:) = kgrid%bz_grid_addresses(:, kgrid%bz_grid_map(kgrid%ir_grid_indices(i) + 1) + 1)
            kpoints_inter%vkpt(:, i) = real(matmul(transpose(kgrid%QDinv), bz_adrs), kind=q)
            kpoints_inter%wtkpt(i) = real(kgrid%ir_grid_weights(i), kind=q)
        end do
        kpoints_inter%wtkpt(:) = kpoints_inter%wtkpt(:)/nkpts_full

        !
        ! Basic parameters that are used in various locations
        !
        kpoints_inter%nkpx = kgrid%D_diag(1)
        kpoints_inter%nkpy = kgrid%D_diag(2)
        kpoints_inter%nkpz = kgrid%D_diag(3)
        kpoints_inter%nkfull = product(kgrid%D_diag)
        kpoints_inter%ILSNF(:, :) = int(transpose(kgrid%P(:, :)))
        kpoints_inter%ILSNFINV = matrix_inverse_3x3i(kpoints_inter%ILSNF)
        kpoints_inter%ISNF(:, :) = 0
        do i = 1, 3
            kpoints_inter%ISNF(i, i) = kgrid%D_diag(i)
        end do
        bk(:, :) = real(transpose(kgrid%microzone_lattice(:, :)), kind=q)
        kpoints_inter%B(:, :) = bk(:, :)
        CALL RECIPS(1._q, bk(1, 1), bk(1, 2), bk(1, 3), ak(1, 1), ak(1, 2), ak(1, 3))
        kpoints_inter%A(:, :) = ak(:, :)
        kpoints_inter%shift(:) = 0

        !
        ! Preparation for tetrahedron method
        !
        call tetirr_hash(kpoints_inter, &
                         transpose(real(kgrid%reciprocal_lattice(:, :), kind=q)), &
                         kpoints_inter%ntet, idtet, -1)
        if (associated(kpoints_inter%idtet)) deallocate (kpoints_inter%idtet)
        allocate (kpoints_inter%idtet(0:4, kpoints_inter%ntet))
        kpoints_inter%idtet = idtet(0:4, :kpoints_inter%ntet)
        kpoints_inter%volwgt = 1.0/(6*nkpts_full)
        deallocate (idtet)

        succeeded = .true.

    end function elphon_kgrid_get_points

    !>
    !> @brief A public subroutine to show the contents of kgrid.
    !>
    !> Parameters
    !> ----------
    !> kgrid : elph_kgrid
    !>     kgrid dataset.
    !> iu6 : integer
    !>    Unit for write.
    !>
    !> by Atsushi Togo
    !>
    subroutine elphon_kgrid_show_dataset(kgrid, iu6)
        type(elph_kgrid), intent(in) :: kgrid
        integer, intent(in) :: iu6

        integer :: i, j, bz_gp, ga(3)

        if (iu6<0) return
        write (iu6, '("Sampling length (2pi/KSPACING):", f8.4)') kgrid%sampling_length
        write (iu6, '("Grid matrix:")')
        write (iu6, '(3i4)') kgrid%grid_matrix(:, 1)
        write (iu6, '(3i4)') kgrid%grid_matrix(:, 2)
        write (iu6, '(3i4)') kgrid%grid_matrix(:, 3)
        write (iu6, '("D_diag: ", 3i4)') kgrid%D_diag(:)
        write (iu6, '("P:")')
        write (iu6, '(3i4)') kgrid%P(:, 1)
        write (iu6, '(3i4)') kgrid%P(:, 2)
        write (iu6, '(3i4)') kgrid%P(:, 3)
        write (iu6, '("Q:")')
        write (iu6, '(3i4)') kgrid%Q(:, 1)
        write (iu6, '(3i4)') kgrid%Q(:, 2)
        write (iu6, '(3i4)') kgrid%Q(:, 3)
        write (iu6, '("Microzone basis vectors:")')
        write (iu6, '("a_m*", 3f13.9)') kgrid%microzone_lattice(1, :)
        write (iu6, '("b_m*", 3f13.9)') kgrid%microzone_lattice(2, :)
        write (iu6, '("c_m*", 3f13.9)') kgrid%microzone_lattice(3, :)
        write (iu6, '("Reciprocal regular rotations wrt prim. cell basis vecs")')
        do i = 1, kgrid%n_rotations
            write (iu6, '(i2, ":", 9i4)') i, kgrid%regular_rotations(:, :, i)
        end do
        write (iu6, '("SNF transformed reciprocal rotations")')
        do i = 1, kgrid%n_rotations
            write (iu6, '(i2, ":", 9i4)') i, kgrid%rotations(:, :, i)
        end do
        write (iu6, '("Number of gr-grid-points:", i5)') product(kgrid%D_diag)
        write (iu6, '("Number of ir-grid-points:", i5)') kgrid%n_ir_grid_points
        write (iu6, '("Number of bz-grid-points:", i5)') kgrid%n_bz_grid_points
        write (iu6, '("BZ grid addresses:")')
        write (iu6, '("Index (BZ-mult-index) address k-point sym-map *:ir-grid-point")')
        do i = 1, product(kgrid%D_diag)
            bz_gp = kgrid%bz_grid_map(kgrid%ir_grid_indices(kgrid%grg2irgrg(i)) + 1)
            do j = kgrid%bz_grid_map(i), kgrid%bz_grid_map(i + 1) - 1
                ga(:) = kgrid%bz_grid_addresses(:, j + 1)
                write (iu6, &
                       '(i0, " (", i0, "): QD^-1 ( ", i0, " ", i0, " ", i0, " ) = (", 3f8.4, " ) map to ", i0)', &
                       advance="no") j + 1, j + 1 - kgrid%bz_grid_map(i), ga, &
                    matmul(transpose(kgrid%QDinv), ga), bz_gp + 1
                if (kgrid%bz_grid_map(i) == bz_gp) then
                    write (iu6, '(" *")')
                else
                    write (iu6, *)
                end if
            end do
        end do
        write (iu6, '("")')
        write (iu6, '("SNF transformed relative grid addresses:")')
        do i = 1, 24
            do j = 1, 4
                write (iu6, '(3i3, "  ")', advance="no") kgrid%relative_grid_addresses(:, j, i)
            end do
            write (iu6, *)
        end do
        write (iu6, *)

    end subroutine elphon_kgrid_show_dataset

    !>
    !> @brief A public subroutine to get triplet integration weights at a grid
    !> point.
    !>
    !> Parameters
    !> ----------
    !> kgrid : elph_kgrid
    !>     kgrid dataset.
    !>
    !> Returns
    !> -------
    !> triplets : elph_triplets
    !>        Irreducible triplets.
    !>
    !> by Atsushi Togo
    !>
    function elphon_kgrid_get_triplets_at_q(kgrid, ikp_ibz, iu6) &
        result(triplets)
        use gridsysf, only: gridsys_get_triplets_at_q, &
                            gridsys_get_bz_triplets_at_q

        type(elph_kgrid), intent(in) :: kgrid
        integer, intent(in) :: ikp_ibz
        integer, intent(in) :: iu6

        type(elph_triplets) :: triplets
        integer num_ir_triplets, num_ir_triplets_bz, num_ir_weights, i
        integer(c_long) :: grid_point, is_time_reversal, swappable, num_gps
        ! map_triplets : Indices of triplets is mapped to indices of irreducible
        ! triplets.
        integer(c_long), allocatable :: map_triplets(:)
        ! map_q : Grid indices are mapped to ir-grid-indices stabilized by q.
        integer(c_long), allocatable :: map_q(:)

        integer(c_long), allocatable :: weights(:)

        grid_point = kgrid%ir_grid_indices(ikp_ibz)
        is_time_reversal = 1
        swappable = 0
        num_gps = product(kgrid%D_diag)

        ! Search ir-triplets of q (grid_point) under q+k'+k=G.
        ! The triplets ((q, k', k), ...) are obtained below by
        ! gridsys_get_bz_triplets_at_q.
        triplets%num_ir_triplets = 0
        allocate (map_triplets(num_gps))
        allocate (map_q(num_gps))
        num_ir_triplets = gridsys_get_triplets_at_q(map_triplets, map_q, &
                                                    grid_point, kgrid%D_diag, is_time_reversal, kgrid%n_rotations, &
                                                    kgrid%rotations, swappable)

        ! Compute weights of ir-triplets of q.
        allocate (weights(num_gps))
        weights(:) = 0
        do i = 1, num_gps
            weights(map_triplets(i) + 1) = weights(map_triplets(i) + 1) + 1
        end do
        allocate (triplets%ir_weights(num_ir_triplets))
        num_ir_weights = 0
        do i = 1, num_gps
            if (weights(i) > 0) then
                num_ir_weights = num_ir_weights + 1
                triplets%ir_weights(num_ir_weights) = weights(i)
            end if
        end do

        ! Sanity check
        if (num_ir_triplets /= num_ir_weights) then
            deallocate (triplets%ir_weights)
            triplets%num_ir_triplets = -2
            return
        end if

        ! Generate triplets array ((q, k', k), ...) under q+k'+k=G.
        allocate (triplets%ir_triplets(3, num_ir_triplets))
        num_ir_triplets_bz = gridsys_get_bz_triplets_at_q( &
                             triplets%ir_triplets, kgrid%bz_grid_map(grid_point + 1), kgrid%bz_grid_addresses, &
                             kgrid%bz_grid_map, &
                             map_triplets, num_gps, kgrid%D_diag, kgrid%Q, int(2, c_long))

        ! Sanity check
        if (num_ir_triplets /= num_ir_triplets_bz) then
            deallocate (triplets%ir_triplets)
            deallocate (triplets%ir_weights)
            triplets%num_ir_triplets = -1
            return
        end if

        triplets%num_ir_triplets = num_ir_triplets

        deallocate (weights)
        deallocate (map_q)
        deallocate (map_triplets)

        if (iu6 >= 0) then
            write (iu6, '("Number of irreducible triplets at ", i0, ": ", i0)') &
                ikp_ibz, triplets%num_ir_triplets
        end if

        ! Invert k in triplets ((q, k', k), ...) -> ((q, k', -k), ...)
        call invert_k_in_triplets(triplets, kgrid, iu6)
        ! Rotate triplets ((q, k', k), ...) -> ((Rq, Rk', Rk), ...) to let
        ! Rk becomes one of irreducilbe k-points in KPOINTS_INTER.
        call rotate_triplets(triplets, kgrid, iu6)

    end function elphon_kgrid_get_triplets_at_q

    !>
    !> @brief A public subroutine to get triplet integration weights at a grid
    !> point.
    !>
    !> Parameters
    !> ----------
    !> rec_rotation : integer
    !>     Rotation matrix in reciprocal space with respect to the primitive
    !>     cell. Defined as k'=kR.
    !> kgrid : elph_kgrid
    !>     kgrid dataset.
    !>
    !> by Atsushi Togo
    !>
    subroutine elphon_kgrid_get_triplet_integration_weights(rec_rotation, kgrid)
        use gridsysf, only: gridsys_transform_rotations

        integer, intent(in) :: rec_rotation(3, 3)
        type(elph_kgrid), intent(in) :: kgrid

        integer :: succeeded
        integer(c_long) :: rec_rotation_c_long(3, 3, 1), rotation(3, 3, 1)

        rec_rotation_c_long(:, :, 1) = int(rec_rotation(:, :), c_long)

        succeeded = gridsys_transform_rotations(rotation, rec_rotation_c_long, &
                                                int(1, c_long), kgrid%D_diag, kgrid%Q)

    end subroutine elphon_kgrid_get_triplet_integration_weights

    !>
    !> @brief compute total dos per cell and energy grid
    !>
    !> by Laurent Chaput
    !>
    subroutine elphon_kgrid_get_dos(WDES_DENSE_IBZ, W_DENSE_IBZ, kgrid, energy_points, dos)
        use tutor, only: vtutor
        use string, only: str
        TYPE(wavedes), intent(in) :: WDES_DENSE_IBZ
        TYPE(wavespin), intent(in) :: W_DENSE_IBZ
        type(elph_kgrid), intent(in) :: kgrid
        real(q), allocatable, intent(inout) :: energy_points(:)
        real(q), allocatable, intent(out) :: dos(:, :)
        integer :: nbe
        ! For Atsushi settings
        real(c_double), allocatable :: iw(:, :, :, :)
        integer(c_long) :: ir_gi, is, ib
        integer(c_long):: nbe_c
        real(c_double), allocatable :: energy_points_c(:)

        PROFILING_START('elphon_kgrid_get_dos')
        if (.not. allocated(energy_points)) then
            call vtutor%bug('energy_points is not allocated',__FILE__,__LINE__)
        end if

        if (allocated(dos)) then
            deallocate (dos)
        end if

        nbe = size(energy_points)
        allocate (dos(nbe, 2))
        dos = 0
        allocate (iw(WDES_DENSE_IBZ%NB_TOT, nbe, WDES_DENSE_IBZ%ISPIN, 2))

        ! very bad ...
        nbe_c = int(nbe, c_long)
        allocate (energy_points_c(nbe_c))
        energy_points_c = real(energy_points, c_double)

        do ir_gi = 1, kgrid%n_ir_grid_points
            iw(:, :, :, :) = 0
            call get_integration_weights(iw, kgrid%ir_grid_indices(ir_gi), &
                                         energy_points_c, nbe_c, &
                                         W_DENSE_IBZ%CELTOT, &
                                         WDES_DENSE_IBZ%NB_TOT, &
                                         WDES_DENSE_IBZ%ISPIN, &
                                         kgrid)

            do is = 1, WDES_DENSE_IBZ%ISPIN
                do ib = 1, WDES_DENSE_IBZ%NB_TOT
                    dos(:, :) = dos(:, :) + real(iw(ib, :, is, :)*kgrid%ir_grid_weights(ir_gi), q)
                end do
            end do
        end do

        dos(:, :) = dos(:, :)/product(kgrid%D_diag)
        dos(:, :) = dos(:, :)*WDES_DENSE_IBZ%RSPIN !!! is it correct ? the name would have said NUP_DOWN ....

        deallocate (iw)
        deallocate (energy_points_c)
        PROFILING_STOP('elphon_kgrid_get_dos')
    end subroutine elphon_kgrid_get_dos

    subroutine elphon_kgrid_get_dos_full(WDES_DENSE_IBZ, W_DENSE_IBZ, kgrid, energy_points, dos)
        use tutor, only: vtutor
        use string, only: str
        TYPE(wavedes), intent(in) :: WDES_DENSE_IBZ
        TYPE(wavespin), intent(in) :: W_DENSE_IBZ
        type(elph_kgrid), intent(in) :: kgrid
        real(q), allocatable, intent(inout) :: energy_points(:)
        real(q), allocatable, intent(out) :: dos(:, :)
        integer :: nbe
        ! For Atsushi settings
        real(c_double), allocatable :: iw(:, :, :, :)
        integer(c_long) :: ir_gi, is, ib
        integer(c_long):: nbe_c
        real(c_double), allocatable :: energy_points_c(:)

        PROFILING_START('elphon_kgrid_get_dos_full')
        if (.not. allocated(energy_points)) then
            call vtutor%bug('energy_points is not allocated',__FILE__,__LINE__)
        end if

        if (allocated(dos)) then
            deallocate (dos)
        end if

        nbe = size(energy_points)
        allocate (dos(nbe, 2))
        dos = 0
        allocate (iw(WDES_DENSE_IBZ%NB_TOT, nbe, WDES_DENSE_IBZ%ISPIN, 2))

        ! very bad ...
        nbe_c = int(nbe, c_long)
        allocate (energy_points_c(nbe_c))
        energy_points_c = real(energy_points, c_double)

        do ir_gi = 1, product(kgrid%D_diag)
            iw(:, :, :, :) = 0
            call get_integration_weights(iw, ir_gi - 1, &
                                         energy_points_c, nbe_c, &
                                         W_DENSE_IBZ%CELTOT, &
                                         WDES_DENSE_IBZ%NB_TOT, &
                                         WDES_DENSE_IBZ%ISPIN, &
                                         kgrid)

            do is = 1, WDES_DENSE_IBZ%ISPIN
                do ib = 1, WDES_DENSE_IBZ%NB_TOT
                    dos(:, :) = dos(:, :) + real(iw(ib, :, is, :), q)
                end do
            end do
        end do

        dos(:, :) = dos(:, :)/product(kgrid%D_diag)
        dos(:, :) = dos(:, :)*WDES_DENSE_IBZ%RSPIN !!! is it correct ? the name would have said NUP_DOWN ....

        deallocate (iw)
        deallocate (energy_points_c)
        PROFILING_STOP('elphon_kgrid_get_dos_full')
    end subroutine elphon_kgrid_get_dos_full

    !>
    !> @brief get fermi energy using bisection method.
    !>
    !> by Laurent Chaput
    !>
    subroutine elphon_kgrid_get_fermi(NELEC, WDES_DENSE_IBZ, W_DENSE_IBZ, kgrid, fermi)
        use elphon_base, only: celtot_max
        real(q), intent(in) :: NELEC
        TYPE(wavedes), intent(in) :: WDES_DENSE_IBZ
        TYPE(wavespin), intent(in) :: W_DENSE_IBZ
        type(elph_kgrid), intent(in) :: kgrid
        real(q), intent(out) :: fermi
        real(q), allocatable :: energy_points(:)
        real(q), allocatable :: dos(:, :)
        real(q) E0, N0, N01, E1, E2
        real(q), parameter :: TINY = 1E-10_q ! to be consistent with converged_num_electron value in tet.F

        allocate (energy_points(1))
        allocate (dos(1, 2))
        E1 = minval(real(W_DENSE_IBZ%CELTOT))
        E2 = celtot_max(W_DENSE_IBZ%CELTOT,WDES_DENSE_IBZ%NB_TOTK)
        do
            E0 = 0.5_q*(E1 + E2)
            energy_points(1) = E0
            call elphon_kgrid_get_dos(WDES_DENSE_IBZ, W_DENSE_IBZ, kgrid, energy_points, dos)
            N0 = dos(1, 2)
            if (N0 < NELEC) E1 = E0
            if (N0 > NELEC) E2 = E0
            if (abs(N0 - N01) < TINY) exit
            if (abs(N0 - NELEC) < TINY) exit
            N01 = N0
        end do
        deallocate (energy_points)
        deallocate (dos)
        fermi = E0
    end subroutine elphon_kgrid_get_fermi

    !>
    !> @brief compute spectral functions
    !>
    !> spectral_1(E) = 1/N sum_{n k s} A(n k s) delta(E-E_{n k s})
    !> spectral_2(E) = 1/N sum_{n k s} A(n k s) theta(E-E_{n k s})
    !>
    !> If the energy_points are given, they are used, otherwise they are
    !> computed from the eigenvalue spectrum, +-5%.
    !>
    !> by Laurent Chaput
    !>
    subroutine elphon_kgrid_get_spectral(A, WDES_DENSE_IBZ, W_DENSE_IBZ, kgrid, &
                                         n_energy_points, energy_points, spectral)
        use tutor, only: vtutor
        use string, only: str
        real(q), intent(in) :: A(:, :, :, :, :)
        TYPE(wavedes), intent(in) :: WDES_DENSE_IBZ
        TYPE(wavespin), intent(in) :: W_DENSE_IBZ
        type(elph_kgrid) :: kgrid
        integer, intent(in) :: n_energy_points
        real(q), allocatable, intent(inout) :: energy_points(:)
        real(q), allocatable, intent(out) :: spectral(:, :, :, :)
        integer :: ie
        integer :: nbe, alpha, beta, m_alpha, m_beta
        ! For Atsushi setting
        real(c_double), allocatable :: iw(:, :, :, :)
        integer(c_long) :: ir_gi, is, ib
        integer(c_long):: nbe_c
        real(c_double), allocatable :: energy_points_c(:)

        PROFILING_START('elphon_kgrid_get_spectral')
        m_alpha = size(A, 1)
        m_beta = size(A, 2)

        if (.not. allocated(energy_points)) then
            call vtutor%bug('energy_points is not allocated',__FILE__,__LINE__)
        end if

        if (size(energy_points,1) /= n_energy_points) then
            call vtutor%bug('incorrect dimensions of energy_points expected'//str(n_energy_points)//&
                            'found '//str(size(energy_points,1)),__FILE__,__LINE__)
        endif

        nbe = size(energy_points)
        allocate (spectral(nbe, m_alpha, m_beta, WDES_DENSE_IBZ%ISPIN))
        spectral = 0.0_q
        allocate (iw(WDES_DENSE_IBZ%NB_TOT, nbe, WDES_DENSE_IBZ%ISPIN, 2))

        ! very bad ...
        nbe_c = int(nbe, c_long)
        allocate (energy_points_c(nbe_c))
        energy_points_c = real(energy_points, c_double)

        do ir_gi = 1, kgrid%n_ir_grid_points
            iw(:, :, :, :) = 0
            call get_integration_weights(iw, kgrid%ir_grid_indices(ir_gi), &
                                         energy_points_c, nbe_c, &
                                         W_DENSE_IBZ%CELTOT, &
                                         WDES_DENSE_IBZ%NB_TOT, &
                                         WDES_DENSE_IBZ%ISPIN, &
                                         kgrid)
            do beta = 1, m_beta
                do alpha = 1, m_alpha
                    do is = 1, WDES_DENSE_IBZ%ISPIN
                        do ib = 1, WDES_DENSE_IBZ%NB_TOT
                            do ie = 1, nbe
                                spectral(ie, alpha, beta, is) = spectral(ie, alpha, beta, is) + &
                                                                real(iw(ib, ie, is, 1)*kgrid%ir_grid_weights(ir_gi), q)* &
                                                                A(alpha, beta, ib, ir_gi, is)
                      !! the energy array index has to be checked with Atsushi.
                      !! That for example the order of the irreducible kpoints in kgrid are
                      !! the same as in CELTOT(ib, ir_gi, is)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        spectral = spectral/product(kgrid%D_diag)
        spectral = spectral*WDES_DENSE_IBZ%RSPIN

        PROFILING_STOP('elphon_kgrid_get_spectral')
    end subroutine elphon_kgrid_get_spectral

    subroutine elphon_kgrid_get_spectral_full(A, WDES_DENSE_IBZ, W_DENSE_IBZ, kgrid, &
                                              n_energy_points, energy_points, spectral)
        use tutor, only: vtutor
        use string, only: str
        real(q), intent(in) :: A(:, :, :, :, :)
        TYPE(wavedes), intent(in) :: WDES_DENSE_IBZ
        TYPE(wavespin), intent(in) :: W_DENSE_IBZ
        type(elph_kgrid) :: kgrid
        integer, intent(in) :: n_energy_points
        real(q), allocatable, intent(inout) :: energy_points(:)
        real(q), allocatable, intent(out) :: spectral(:, :, :, :)
        integer :: ie
        integer :: nbe, alpha, beta, m_alpha, m_beta
        ! For Atsushi setting
        real(c_double), allocatable :: iw(:, :, :, :)
        integer(c_long) :: ir_gi, is, ib
        integer(c_long):: nbe_c
        real(c_double), allocatable :: energy_points_c(:)

        m_alpha = size(A, 4)
        m_beta = size(A, 5)

        if (.not. allocated(energy_points)) then
            call vtutor%bug('energy_points is not allocated',__FILE__,__LINE__)
        end if

        if (size(energy_points,1) /= n_energy_points) then
            call vtutor%bug('incorrect dimensions of energy_points expected'//str(n_energy_points)//&
                            'found '//str(size(energy_points,1)),__FILE__,__LINE__)
        endif

        nbe = size(energy_points)
        allocate (spectral(nbe, m_alpha, m_beta, WDES_DENSE_IBZ%ISPIN))
        spectral = 0
        allocate (iw(WDES_DENSE_IBZ%NB_TOT, nbe, WDES_DENSE_IBZ%ISPIN, 2))

        ! very bad ...
        nbe_c = int(nbe, c_long)
        allocate (energy_points_c(nbe_c))
        energy_points_c = real(energy_points, c_double)

        do ir_gi = 1, product(kgrid%D_diag)
            iw(:, :, :, :) = 0
            call get_integration_weights(iw, ir_gi - 1, &
                                         energy_points_c, nbe_c, &
                                         W_DENSE_IBZ%CELTOT, &
                                         WDES_DENSE_IBZ%NB_TOT, &
                                         WDES_DENSE_IBZ%ISPIN, &
                                         kgrid)
            do beta = 1, m_beta
                do alpha = 1, m_alpha
                    do is = 1, WDES_DENSE_IBZ%ISPIN
                        do ib = 1, WDES_DENSE_IBZ%NB_TOT
                            do ie = 1, nbe
                                spectral(ie, alpha, beta, is) = spectral(ie, alpha, beta, is) + real(iw(ib, ie, is, 1), q)* &
                                                                A(alpha, beta, ib, kgrid%grg2irgrg(ir_gi), is)
                      !! the energy array index has to be checked with Atsushi.
                      !! That for example the order of the irreducible kpoints in kgrid are
                      !! the same as in CELTOT(ib, ir_gi, is)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        spectral = spectral/product(kgrid%D_diag)
        spectral = spectral*WDES_DENSE_IBZ%RSPIN

    end subroutine elphon_kgrid_get_spectral_full

    !>
    !> @brief A private subroutine to get grid generating matrices.
    !>
    !> by Atsushi Togo
    !>
    function get_grid_matrices(kgrid, lattice) result(succeeded)
        use gridsysf, only: gridsys_get_snf3x3
        use mathtools, only: matrix_inverse_3x3

        type(elph_kgrid), intent(inout) :: kgrid
        real(q), intent(in)  :: lattice(3, 3)

        real(q) :: rec_lattice(3, 3)
        real(q) :: inv_tmat(3, 3), std_lattice(3, 3), rec_std_lattice(3, 3)
        integer :: i, j, succeeded
        integer :: mesh_numbers(3)

        ! row vectors
        rec_lattice(:, :) = matrix_inverse_3x3(lattice)
        kgrid%reciprocal_lattice(:, :) = real(rec_lattice(:, :), kind=c_double)

        if (kgrid%sampling_length > 0) then
            ! a, b, c are column vectors.
            inv_tmat(:, :) = transpose(matrix_inverse_3x3(real(kgrid%transformation_matrix, kind=q)))
            std_lattice(:, :) = matmul(lattice, inv_tmat)

            ! row vectors
            rec_std_lattice(:, :) = matrix_inverse_3x3(std_lattice)
            mesh_numbers(:) = get_mesh_numbers(kgrid%sampling_length, rec_std_lattice)

            ! tmat of p->c
            do i = 1, 3
                do j = 1, 3
                    kgrid%grid_matrix(i, j) = nint(inv_tmat(i, j))*mesh_numbers(j)
                end do
            end do

            kgrid%microzone_lattice(:, :) = real( &
                                            matmul(matrix_inverse_3x3(real(kgrid%grid_matrix, kind=q)), &
                                                   kgrid%reciprocal_lattice(:, :)), kind=c_double)

            succeeded = gridsys_get_snf3x3(kgrid%D_diag, kgrid%P, kgrid%Q, &
                                           kgrid%grid_matrix)

            do i = 1, 3
                kgrid%QDinv(i, :) = real(kgrid%Q(i, :), kind=q)/kgrid%D_diag(i)
            end do
        else
            mesh_numbers(:) = get_mesh_numbers(-kgrid%sampling_length, rec_lattice)
            kgrid%grid_matrix(:, :) = 0.0
            kgrid%P(:, :) = 0
            kgrid%Q(:, :) = 0
            kgrid%QDinv(:, :) = 0.0
            do i = 1, 3
                kgrid%grid_matrix(i, i) = mesh_numbers(i)
                kgrid%D_diag(i) = mesh_numbers(i)
                kgrid%P(i, i) = 1
                kgrid%Q(i, i) = 1
                kgrid%QDinv(i, i) = 1.0_q/kgrid%D_diag(i)
            end do
            kgrid%microzone_lattice(:, :) = real( &
                                            matmul(matrix_inverse_3x3(real(kgrid%grid_matrix, kind=q)), &
                                                   kgrid%reciprocal_lattice(:, :)), kind=c_double)
        end if

    end function get_grid_matrices

    !>
    !> @brief A private function to get reciprocal sampling mesh numbers
    !> for the standardized basis vectors (not primitive basis vectors.)
    !>
    !> by Atsushi Togo
    !>
    function get_mesh_numbers(sampling_length, rec_lattice) result(mesh_numbers)
        real(c_double), intent(in) :: sampling_length
        real(q), intent(in) :: rec_lattice(3, 3)
        integer :: mesh_numbers(3)

        integer :: i

        do i = 1, 3
            mesh_numbers(i) = nint(norm2(rec_lattice(i, :))*sampling_length)
            if (mesh_numbers(i) < 1) then
                mesh_numbers(i) = 1
            end if
        end do
    end function get_mesh_numbers

    !>
    !> @brief A private function to get reciprocal rotation matrices.
    !>
    !> by Atsushi Togo
    !>
    function get_reciprocal_rotations(kgrid) result(succeeded)
        use gridsysf, only: gridsys_get_reciprocal_point_group, &
                            gridsys_transform_rotations

        type(elph_kgrid), intent(inout) :: kgrid
        integer :: succeeded

        integer(c_long) :: is_time_reversal = 1

        kgrid%regular_rotations(:, :, :) = 0
        kgrid%n_rotations = gridsys_get_reciprocal_point_group(kgrid%regular_rotations, &
                                                               kgrid%direct_rotations, kgrid%n_direct_rotations, is_time_reversal)

        kgrid%rotations(:, :, :) = 0
        succeeded = gridsys_transform_rotations(kgrid%rotations, &
                                                kgrid%regular_rotations, kgrid%n_rotations, kgrid%D_diag, kgrid%Q)
    end function get_reciprocal_rotations

    !>
    !> @brief A private subroutine to get relative grid addresses used by
    !> tetrahedron method.
    !>
    !> by Atsushi Togo
    !>
    subroutine get_relative_grid_addresses(kgrid)
        use gridsysf, only: gridsys_get_thm_relative_grid_address
        use mathtools, only: matrix_inverse_3x3

        type(elph_kgrid), intent(inout) :: kgrid

        integer(c_long) :: relative_grid_addresses(3, 4, 24)
        integer :: i, j, k, l

        call gridsys_get_thm_relative_grid_address(relative_grid_addresses, &
                                                   kgrid%microzone_lattice)

        kgrid%relative_grid_addresses(:, :, :) = 0
        do i = 1, 24
            do j = 1, 4
                do k = 1, 3
                    do l = 1, 3
                        kgrid%relative_grid_addresses(k, j, i) = &
                            kgrid%relative_grid_addresses(k, j, i) + &
                            relative_grid_addresses(l, j, i)*kgrid%P(l, k)
                    end do
                end do
            end do
        end do

    end subroutine get_relative_grid_addresses

    !>
    !> @brief A private subroutine to get irreducible grid point map.
    !> The result is stored in kgrid.
    !>
    !> by Atsushi Togo
    !>
    subroutine get_grg2irgrg(kgrid)
        use gridsysf, only: gridsys_get_ir_grid_map

        type(elph_kgrid), intent(inout) :: kgrid
        integer :: i, ir_index
        integer(c_long), allocatable :: gr2ir(:)
        integer(c_long), allocatable :: weights(:)

        allocate (gr2ir(product(kgrid%D_diag)))

        call gridsys_get_ir_grid_map(gr2ir, kgrid%rotations, &
                                     kgrid%n_rotations, kgrid%D_diag, kgrid%PS)

        allocate (weights(product(kgrid%D_diag)))
        weights(:) = 0
        kgrid%n_ir_grid_points = 0
        do i = 1, product(kgrid%D_diag)
            if (gr2ir(i) == i - 1) then
                kgrid%n_ir_grid_points = kgrid%n_ir_grid_points + 1
            end if
            weights(gr2ir(i) + 1) = weights(gr2ir(i) + 1) + 1
        end do

        if (allocated(kgrid%ir_grid_indices)) then
            deallocate (kgrid%ir_grid_indices)
        end if
        allocate (kgrid%ir_grid_indices(kgrid%n_ir_grid_points))

        if (allocated(kgrid%ir_grid_weights)) then
            deallocate (kgrid%ir_grid_weights)
        end if
        allocate (kgrid%ir_grid_weights(kgrid%n_ir_grid_points))

        if (allocated(kgrid%grg2irgrg)) then
            deallocate (kgrid%grg2irgrg)
        end if
        allocate (kgrid%grg2irgrg(product(kgrid%D_diag)))

        ir_index = 0
        kgrid%grg2irgrg(:) = -1
        do i = 1, product(kgrid%D_diag)
            if (weights(i) /= 0) then
                ir_index = ir_index + 1
                kgrid%ir_grid_indices(ir_index) = i - 1
                kgrid%ir_grid_weights(ir_index) = weights(i)
                kgrid%grg2irgrg(i) = ir_index
            else
                kgrid%grg2irgrg(i) = kgrid%grg2irgrg(gr2ir(i) + 1)
            end if
        end do

        deallocate (gr2ir)
        deallocate (weights)
    end subroutine get_grg2irgrg

#ifdef use_spglib
    !>
    !> @brief A private function to store symmetry information
    !> obtained by spglib.
    !>
    !> by Atsushi Togo
    !>
    subroutine get_symmetry_info_from_spglib(kgrid, lattice, positions, ntypes, &
                                             nitypes, natoms)
        use sym_prec, only: tiny
        use spglib_f08, only: SpglibDataset, spg_get_dataset

        type(elph_kgrid), intent(inout) :: kgrid
        integer, intent(in)  :: ntypes, natoms
        integer, intent(in)  :: nitypes(ntypes)
        real(q), intent(in)  :: lattice(3, 3)
        real(q), intent(in)  :: positions(3, natoms)
        integer :: i, j, k
        real(c_double) :: lattice_c_t(3, 3)
        real(c_double) :: positions_c(3, natoms)
        integer(c_int) :: atom_types(natoms)
        type(SpglibDataset) :: dataset

        k = 1
        do i = 1, ntypes
            do j = 1, nitypes(i)
                atom_types(k) = i
                k = k + 1
            end do
        end do

        ! Need to transpose
        lattice_c_t(:, :) = real(transpose(lattice(:, :)), kind=c_double)
        positions_c(:, :) = real(positions(:, :), kind=c_double)

        dataset = spg_get_dataset(lattice_c_t, positions_c, atom_types, &
                                  int(natoms, kind=c_int), real(tiny, kind=c_double))

        kgrid%spacegroup_number = dataset%spacegroup_number
        kgrid%transformation_matrix(:, :) = dataset%transformation_matrix(:, :)

        kgrid%n_direct_rotations = dataset%n_operations
        if (allocated(kgrid%direct_rotations)) then
            deallocate (kgrid%direct_rotations)
        end if
        allocate (kgrid%direct_rotations(3, 3, dataset%n_operations))
        do i = 1, dataset%n_operations
            kgrid%direct_rotations(:, :, i) = int(dataset%rotations(:, :, i), kind=c_long)
        end do

    end subroutine get_symmetry_info_from_spglib
#endif

    !>
    !> @brief A private subroutine to store symmetry information
    !> obtained from phelel_params.hdf5.
    !>
    !> by Atsushi Togo
    !>
    subroutine get_symmetry_info_from_hdf5(kgrid, filename)
        use vhdf5

        type(elph_kgrid), intent(inout) :: kgrid
        character(len=*), optional, intent(in) :: filename

#ifdef VASP_HDF5
        integer :: spacegroup_number, rank
        integer(hid_t):: fileid
        character(len=100) :: my_filename
        integer, allocatable :: direct_rotations(:, :, :)
        real(q) :: transformation_matrix(3, 3)
        integer(hsize_t) :: dimensions(7)

        my_filename = 'phelel_params.hdf5'; if (present(filename)) my_filename = filename
        VH5_CHECK(vh5_file_open_read(my_filename, fileid))

        VH5_CHECK(vh5_get_dimensions(fileid, 'direct_rotations', dimensions, rank))
        kgrid%n_direct_rotations = dimensions(3)
        allocate (direct_rotations(3, 3, kgrid%n_direct_rotations))

        VH5_CHECK(vh5_read(fileid, 'transformation_matrix', transformation_matrix))
        VH5_CHECK(vh5_read(fileid, 'direct_rotations', direct_rotations))
        VH5_CHECK(vh5_read(fileid, 'spacegroup_number', spacegroup_number))

        VH5_CHECK(vh5_file_close(fileid))

        kgrid%spacegroup_number = spacegroup_number
        kgrid%transformation_matrix(:, :) = transformation_matrix(:, :)

        if (allocated(kgrid%direct_rotations)) then
            deallocate (kgrid%direct_rotations)
        end if
        allocate (kgrid%direct_rotations(3, 3, kgrid%n_direct_rotations))
        kgrid%direct_rotations(:, :, :) = direct_rotations(:, :, :)

        deallocate (direct_rotations)
#endif
    end subroutine get_symmetry_info_from_hdf5

    !>
    !> @brief A private subroutine to get integration weights of tetrahedron
    !> methodat a grid point.
    !>
    !> Parameters
    !> ----------
    !> iw : real(q)
    !>     Integration weights.
    !>     shape=(n_bands, n_energy_points, ispin, 2)
    !> grid_index : integer(c_long)
    !>     Grid point index.
    !> energy_points : real(c_double)
    !>     Energy points where integration weights are computed.
    !>     shape=(n_energy_points, )
    !> n_energy_points : integer
    !>     Number of energy points.
    !> energies : real(q)
    !>     Energies at band indices, irreducible k-points, and spins.
    !>     shape=(n_bands, n_ir_grid_points, ispin)
    !> n_bands : integer
    !>     Number of bands.
    !> ispin : integer
    !>     ispin.
    !> kgrid : elph_kgrid
    !>     kgrid dataset.
    !>
    !> by Atsushi Togo
    !>
    subroutine get_integration_weights(iw, grid_index, &
                                       energy_points, n_energy_points, energies, n_bands, ispin, kgrid)
        use gridsysf, only: gridsys_get_thm_integration_weight, &
                            gridsys_get_double_grid_address, gridsys_get_double_grid_index, &
                            gridsys_get_grid_address_from_index

        integer(c_long), intent(in) :: n_energy_points
        real(c_double), intent(inout) :: iw(n_bands, n_energy_points, ispin, 2)
        integer(c_long), intent(in) :: grid_index
        real(c_double), intent(in) :: energy_points(n_energy_points)
        complex(q), intent(in) :: energies(:, :, :)
        integer, intent(in) :: ispin
        integer, intent(in) :: n_bands
        type(elph_kgrid), intent(in) :: kgrid

        integer :: ij, ie, is, ib, i, j
        integer(c_long) :: gi_neighbor
        integer(c_long) :: ir_gi_neighbor
        integer(c_long) :: address(3)
        integer(c_long) :: adrs_neighbor(3)
        integer(c_long) :: dadrs_neighbor(3)
        real(c_double) :: tetrahedra_energies(4, 24, n_bands, ispin)
        character(kind=c_char) :: function_chars(2)

        call gridsys_get_grid_address_from_index(address, grid_index, kgrid%D_diag)
        do i = 1, 24
            do j = 1, 4
                adrs_neighbor(:) = kgrid%relative_grid_addresses(:, j, i) + &
                                   address(:)
                call gridsys_get_double_grid_address(dadrs_neighbor, adrs_neighbor, &
                                                     kgrid%PS)
                gi_neighbor = gridsys_get_double_grid_index(dadrs_neighbor, &
                                                            kgrid%D_diag, kgrid%PS)
                ir_gi_neighbor = kgrid%grg2irgrg(gi_neighbor + 1)
                do ib = 1, n_bands
                    do is = 1, ispin
                        tetrahedra_energies(j, i, ib, is) = &
                            real(energies(ib, ir_gi_neighbor, is),c_double)
                    end do
                end do
            end do
        end do

        function_chars(:) = ["I", "J"]

        do ij = 1, 2
            do is = 1, ispin
                do ie = 1, n_energy_points
                    do ib = 1, n_bands
                        iw(ib, ie, is, ij) = gridsys_get_thm_integration_weight( &
                                             energy_points(ie), tetrahedra_energies(:, :, ib, is), &
                                             function_chars(ij))
                    end do
                end do
            end do
        end do

    end subroutine get_integration_weights

    !>
    !> @brief A private subroutine to get BZ grid addresses.
    !>
    !>
    !> Parameters
    !> ----------
    !> kgrid : elph_kgrid
    !>     kgrid dataset.
    !>
    !> Returns
    !> -------
    !> kgrid%bz_grid_addresses
    !> kgrid%bz_grid_map
    !> kgrid%n_bz_grid_points
    !>
    !> by Atsushi Togo
    !>
    subroutine get_bz_grid_addresses(kgrid)

        use gridsysf, only: gridsys_get_bz_grid_addresses

        type(elph_kgrid), intent(inout) :: kgrid

        integer num_gps, num_bz_grid_points, num_tmp_bz_gps
        integer(c_long), allocatable :: bzg2grg(:)
        integer(c_long), allocatable :: bz_grid_addresses(:, :)

        num_gps = product(kgrid%D_diag)

        if (num_gps > bz_grid_points_estimation) then
            num_tmp_bz_gps = (kgrid%D_diag(1) + 1)*(kgrid%D_diag(2) + 1) &
                             *(kgrid%D_diag(3) + 1)
        else
            num_tmp_bz_gps = num_gps*8  ! This 8 comes from something to make it enough.
        end if
        allocate (bz_grid_addresses(3, num_tmp_bz_gps))
        allocate (bzg2grg(num_tmp_bz_gps))

        if (allocated(kgrid%bz_grid_map)) then
            deallocate (kgrid%bz_grid_map)
        end if
        allocate (kgrid%bz_grid_map(num_gps + 1))

        num_bz_grid_points = gridsys_get_bz_grid_addresses(bz_grid_addresses, &
                                                           kgrid%bz_grid_map, bzg2grg, kgrid%D_diag, &
                                                           kgrid%Q, kgrid%PS, kgrid%reciprocal_lattice, int(2, c_long))

        kgrid%n_bz_grid_points = num_bz_grid_points
        if (allocated(kgrid%bz_grid_addresses)) then
            deallocate (kgrid%bz_grid_addresses)
        end if
        allocate (kgrid%bz_grid_addresses(3, kgrid%n_bz_grid_points))
        kgrid%bz_grid_addresses(:, :) = &
            bz_grid_addresses(:, 1:kgrid%n_bz_grid_points)
        deallocate (bz_grid_addresses)

        if (allocated(kgrid%bzg2grg)) then
            deallocate (kgrid%bzg2grg)
        end if
        allocate (kgrid%bzg2grg(kgrid%n_bz_grid_points))
        kgrid%bzg2grg(:) = bzg2grg(1:kgrid%n_bz_grid_points)
        deallocate (bzg2grg)

    end subroutine get_bz_grid_addresses

    !>
    !> @brief A private subroutine to rotate triplets
    !>
    !> @details This subroutine rotates triplets so that second point of each
    !> triplet, i.e., k of (q, k, k') is sent to its ir-kpoint. The ir-kpoint
    !> is that in the list of ir_grid_indices.
    !>
    subroutine rotate_triplets(triplets, kgrid, iu6)
        use gridsysf, only: gridsys_rotate_bz_grid_index, &
                            gridsys_rotate_grid_index

        type(elph_triplets), intent(inout) :: triplets
        type(elph_kgrid), intent(in) :: kgrid
        integer, intent(in) :: iu6

        integer i, j, k, rot_gr_gp
        integer(c_long) :: PS(3), adrs(3), rot_adrs(3)
        integer(c_long) :: gr_gp_k, rot_bz_gp, gr_gp_irk
        logical :: is_found
        real(q) :: pos(3), rot_pos(3)

        PS(:) = 0

        if (allocated(triplets%rotated_triplets)) then
            deallocate (triplets%rotated_triplets)
        end if
        allocate (triplets%rotated_triplets(3, triplets%num_ir_triplets))
        if (allocated(triplets%rotation_indices)) then
            deallocate (triplets%rotation_indices)
        end if
        allocate (triplets%rotation_indices(triplets%num_ir_triplets))

        do i = 1, triplets%num_ir_triplets
            ! Search ir-grid point in GR-grid.
            ! triplets%ir_triplets is given in BZ-grid.
            ! In some case, rotation can not bring BZ-grid point to ir-BZ-grid point.
            ! So we forget about translation by G vector.
            gr_gp_k = kgrid%bzg2grg(triplets%ir_triplets(3, i) + 1)
            gr_gp_irk = kgrid%ir_grid_indices(kgrid%grg2irgrg(gr_gp_k + 1))
            is_found = .false.
            do j = 1, kgrid%n_rotations
                rot_gr_gp = gridsys_rotate_grid_index( &
                            gr_gp_k, kgrid%rotations(:, :, j), kgrid%D_diag, kgrid%PS)
                if (rot_gr_gp == gr_gp_irk) then
                    is_found = .true.
                    triplets%rotation_indices(i) = j
                    do k = 1, 3
                        rot_bz_gp = gridsys_rotate_bz_grid_index( &
                                    triplets%ir_triplets(k, i), kgrid%rotations(:, :, j), &
                                    kgrid%bz_grid_addresses, kgrid%bz_grid_map, kgrid%D_diag, &
                                    kgrid%PS, int(2, c_long))
                        triplets%rotated_triplets(k, i) = rot_bz_gp
                    end do
                    exit
                end if
            end do

            if (.not. is_found) then
                write (iu6, '("Rotation of triplets failed.")')
            end if
        end do

        if (iu6 >= 0) then
            write (iu6, '("Rotation indices to map k -> ir-k|q+kp-k(G) -> Rq+Rkp-Rk(rot-G)")')
            do i = 1, triplets%num_ir_triplets
                adrs(:) = 0
                rot_adrs(:) = 0
                do j = 1, 3
                    if (j == 3) then
                        adrs(:) = adrs(:) - kgrid%bz_grid_addresses(:, triplets%ir_triplets(j, i) + 1)
                        rot_adrs(:) = rot_adrs(:) - kgrid%bz_grid_addresses(:, triplets%rotated_triplets(j, i) + 1)
                    else
                        adrs(:) = adrs(:) + kgrid%bz_grid_addresses(:, triplets%ir_triplets(j, i) + 1)
                        rot_adrs(:) = rot_adrs(:) + kgrid%bz_grid_addresses(:, triplets%rotated_triplets(j, i) + 1)
                    end if
                end do
                pos(:) = matmul(transpose(kgrid%reciprocal_lattice), matmul(transpose(kgrid%QDinv), adrs))
                rot_pos(:) = matmul(transpose(kgrid%reciprocal_lattice), matmul(transpose(kgrid%QDinv), rot_adrs))
                write (iu6, '(i0, ": ", i0, " -> ", i0, " |", &
&                    " [", i0, " ", i0, " ", i0, "]", &
&                    "(", f6.3, " ", f6.3, " ", f6.3, ")", " -> ", &
&                    " [", i0, " ", i0, " ", i0, "]", &
&                    "(", f6.3, " ", f6.3, " ", f6.3, ")")') &
                triplets%rotation_indices(i), &
                triplets%ir_triplets(3, i) + 1, triplets%rotated_triplets(3, i) + 1, &
                (adrs(j), j=1, 3), (pos(j), j=1, 3), &
                (rot_adrs(j), j=1, 3), (rot_pos(j), j=1, 3)
            end do
        end if

    end subroutine rotate_triplets

    !>
    !> @brief A private subroutine to invert k in triplets
    !>
    !> @details This subroutine invert k -> -k in triplets.
    !>
    subroutine invert_k_in_triplets(triplets, kgrid, iu6)
        use gridsysf, only: gridsys_rotate_bz_grid_index

        type(elph_triplets), intent(inout) :: triplets
        type(elph_kgrid), intent(in) :: kgrid
        integer, intent(in) :: iu6

        integer i, j
        integer(c_long) :: rot_bz_gp
        integer(c_long) :: inversion(3, 3)

        inversion(:, :) = reshape([-1, 0, 0, 0, -1, 0, 0, 0, -1], [3, 3])

        do i = 1, triplets%num_ir_triplets
            rot_bz_gp = gridsys_rotate_bz_grid_index( &
                        triplets%ir_triplets(3, i), inversion, &
                        kgrid%bz_grid_addresses, kgrid%bz_grid_map, kgrid%D_diag, &
                        kgrid%PS, int(2, c_long))

            do j = 1, 3
                if (kgrid%bz_grid_addresses(j, triplets%ir_triplets(3, i) + 1) /= -kgrid%bz_grid_addresses(j, rot_bz_gp + 1)) then
                    rot_bz_gp = -1
                    write (iu6, '(i0, ": k -> -k failed.")') i
                    return
                end if
            end do

            triplets%ir_triplets(3, i) = rot_bz_gp
        end do

    end subroutine invert_k_in_triplets

#endif

end module elphon_kgrid
