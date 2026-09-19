#include "symbol.inc"

module elphon_superconductivity
    use prec
    use constant
    use elphon_accumulators, only: elph_accumulator_selfen_ph
    implicit none

#ifdef no_elph_superconductivity
#else
    private
    public superconductivity_elphon

    !> Structure to hold superconducting properties in the isotropic eliashberg approximation
    type sc_iso_eliashberg
        integer :: ntemps !< number of temperatures
        integer :: nspin !< number of spin components
        integer :: spindeg !< spin degeneracy (1 or 2)
        integer :: natoms !< number of atoms in the unit cell
        real(q), allocatable :: efermi(:) !< fermi energy in eV
        real(q), allocatable :: tempsk(:) !< temperatures in Kelvin
        integer :: nedos_ph ! number of frequency points for the phonon DOS and a2F
        real(q), allocatable :: w(:) !< phonon frequency grid
        real(q), allocatable :: a2F(:, :, :) !< isotropic eliashberg function (frequency, temperature, spin)
        real(q), allocatable :: phonon_dos(:) !< phonon DOS (frequency,)
        real(q), allocatable :: dosef(:, :) !< dos per spin at the chemical potential (temperature, spin)
    end type sc_iso_eliashberg

contains

    subroutine superconductivity_elphon(elph_set, elph_accumulators, phonon_freqs_ibz, celtot, kpoints_inter, latt_cur, io, comm_global)
        use constant
        use mpimy, only: communic
        use mkpoints, only: kpoints_struct
        use lattice, only: latt
        use elphon_accumulators_high, only: elph_accumulators_selfen, THZTOEV
        use elphon_accumulators_ph, only: elph_accumulator_selfen_ph_get_lifetimes
        use elphon_base, only: elph_settings
        use base, only: in_struct
#ifdef VASP_HDF5
        USE vhdf5, ONLY: IH5OUTFILEID
#endif
        type(elph_settings), intent(inout) :: elph_set
        type(elph_accumulators_selfen), intent(inout) :: elph_accumulators
        type(sc_iso_eliashberg), allocatable :: eliash(:)
        real(q), intent(in) :: phonon_freqs_ibz(:, :) !< phonon frequencies in THz
        complex(q), intent(in) :: celtot(:, :, :) !< electron eigenvalues for all bands in the IBZ ! w_dense_ibz%celtot(:, :, :)
        type(kpoints_struct), intent(in) :: kpoints_inter
        type(latt), intent(in) :: latt_cur
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm_global
        !local variables
        integer :: nphonon_selfen, i
        real(q), allocatable :: phonon_lifetimes_ibz(:, :, :, :) !< phonon lifetimes in THz

        if (elph_set%emax_ph .lt. 0.0_q) then
          elph_set%emax_ph=1.01_q*maxval(phonon_freqs_ibz)*thztoev
        endif

        nphonon_selfen = size(elph_accumulators%selfen_ph)
        allocate(eliash(nphonon_selfen))
        ! compute eliashberg function. DOS at Ef is a byproduct.
        do i=1,nphonon_selfen
            ! get phonon lifetimes in a band,kpoint,spin array
            call elph_accumulator_selfen_ph_get_lifetimes(elph_accumulators%selfen_ph(i), phonon_lifetimes_ibz)
            call eliashberg_superconductivity_init(eliash(i), elph_accumulators%selfen_ph(i)%tempsk, elph_accumulators%selfen_ph(i)%natoms, &
                                                   elph_accumulators%selfen_ph(i)%nspin, elph_accumulators%selfen_ph(i)%spindeg, elph_set%nedos_ph, elph_accumulators%selfen_ph(i)%efermi)
            call elph_isotropic_eliashberg_function( &
                eliash(i), kpoints_inter, latt_cur, phonon_freqs_ibz, celtot, phonon_lifetimes_ibz, elph_set%nedos_ph, elph_set%emax_ph)
            deallocate(phonon_lifetimes_ibz)
        end do

        ! compute critical temperature, isotropic gap and renormalization
        call elph_isotropic_Tc(eliash, elph_set%nedos_ph, elph_set%emax_ph,&
                               elph_set%matsubara_temps, elph_set%matsubara_ntemps , &
                               elph_set%max_matsubara*maxval(phonon_freqs_ibz)*THZTOEV,&
                               elph_set%mu_star,io%iu0,comm_global)

        ! write data to hdf5 file
        do i=1,nphonon_selfen
#ifdef VASP_HDF5
            call elph_write_iso_eliashberg_hdf5(eliash(i), ih5outfileid, i)
#endif
            call eliashberg_superconductivity_free(eliash(i))
        end do
    end subroutine superconductivity_elphon

    !> Initialize sc_iso_eliashberg structure
    subroutine eliashberg_superconductivity_init(self, temps, natoms, nspin, spindeg, nedos_ph, efermi)
        type(sc_iso_eliashberg), intent(inout) :: self
        real(q), intent(in) :: temps(:) ! temperatures in Kelvin
        integer, intent(in) :: nspin
        integer, intent(in) :: spindeg
        integer, intent(in) :: natoms
        integer, intent(in) :: nedos_ph ! number of frequency points for the phonon DOS and a2F
        real(q), intent(in) :: efermi(:) !< fermi energy in eV

        self%nspin = nspin
        self%spindeg = spindeg
        self%natoms = natoms
        self%ntemps = size(temps)
        self%nedos_ph = nedos_ph
        allocate(self%tempsk(self%ntemps))
        self%tempsk = temps
        allocate(self%w(self%nedos_ph)) ! initial allocation of 1000 points
        allocate(self%a2F(self%nedos_ph, self%ntemps, self%nspin)) ! initial allocation of 1000 points, 2 spins
        self%a2F = 0.0_q
        allocate(self%phonon_dos(self%nedos_ph))
        self%phonon_dos = 0.0_q
        allocate(self%dosef(self%ntemps, self%nspin)) ! 2 spins
        self%dosef = 0.0_q
        allocate(self%efermi(size(efermi)))
        self%efermi = efermi
    end subroutine eliashberg_superconductivity_init

    subroutine eliashberg_superconductivity_free(self)
        type(sc_iso_eliashberg), intent(inout) :: self
        deallocate(self%phonon_dos)
        deallocate(self%a2F)
        deallocate(self%w)
        deallocate(self%DOSEF)
    end subroutine eliashberg_superconductivity_free

    subroutine get_lambda0(hbar_omega, a2F, lambda0)
        real(q), intent(in) :: hbar_omega(:)
        real(q), intent(in) :: a2F(:)
        real(q), intent(out) :: lambda0
        integer i

        lambda0 = 0.0_q
        do i = 2, size(hbar_omega) - 1
            lambda0 = lambda0 + (hbar_omega(i + 1) - hbar_omega(i))*(2.0_q/hbar_omega(i))*a2F(i)
        end do
    end subroutine get_lambda0

    subroutine get_omega_log(hbar_omega, a2F, lambda0, omega_log)
        real(q), intent(in) :: hbar_omega(:)
        real(q), intent(in) :: a2F(:)
        real(q), intent(in) :: lambda0
        real(q), intent(out) :: omega_log
        integer i

        omega_log = 0.0_q
        do i = 2, size(hbar_omega) - 1
            omega_log = omega_log + (hbar_omega(i + 1) - hbar_omega(i))*(2.0_q/hbar_omega(i))*log(hbar_omega(i))*a2F(i)
        end do
        omega_log = exp(omega_log/lambda0)
    end subroutine get_omega_log

    subroutine isotropic_eliashberg_scf(temperature, n_max, lambda, hbar_omega_n, mu_star, lambda0, omega_log, omega_max)
        use constant
        real(q), intent(in) :: temperature
        integer, intent(in) :: n_max
        real(q), intent(in) :: lambda(-2*n_max:2*n_max)
        real(q), intent(in) :: hbar_omega_n(-n_max:n_max)
        real(q), intent(in) :: mu_star
        real(q), intent(in) :: lambda0
        real(q), intent(in) :: omega_log
        real(q), intent(in) :: omega_max
        real(q), allocatable :: kernel_p(:, :)
        real(q), allocatable :: kernel_m(:, :)
        real(q), allocatable :: delta(:)
        real(q), allocatable :: Z(:)
        real(q), allocatable :: v_m(:)
        real(q), allocatable :: v_p(:)
        real(q), allocatable :: norm(:)
        real(q), allocatable :: accuracy(:)
        real(q) delta_BCS, Tc_Allen_Dynes
        integer n, np
        integer counter

        allocate (kernel_m(-n_max:n_max, -n_max:n_max))
        allocate (kernel_p(-n_max:n_max, -n_max:n_max))
        call eliashberg_kernel_pm(n_max, lambda, hbar_omega_n, mu_star, kernel_p, kernel_m)

        allocate (delta(-n_max:n_max))
        allocate (Z(-n_max:n_max))
        allocate (v_m(-n_max:n_max))
        allocate (v_p(-n_max:n_max))
        allocate (norm(-n_max:n_max))
        allocate (accuracy(-n_max:n_max))

        ! initialization: (omega_log put to Kelvin)
        Tc_Allen_Dynes = ((omega_log/BOLKEV)/1.2_q)*exp(-1.04*(1.0_q + lambda0) &
                                                        /(lambda0 - mu_star*(1.0_q + 0.62_q*lambda0)))
        delta_BCS = (3.528_q/2.0_q)*BOLKEV*Tc_Allen_Dynes
        delta(:) = 0.0_q
        do n = -n_max, n_max
            if (abs(hbar_omega_n(n)) .lt. 2.0_q*omega_max) then
                delta(n) = delta_BCS
            end if
        end do

        ! convergence loop

        accuracy = 1.0_q
        counter = 0
        do while (((maxval(abs(delta)) .gt. 1e-10) .and. (sum(accuracy) .gt. 1e-4_q)) .and. counter .lt. 1000)

            accuracy = delta

            norm = sqrt(hbar_omega_n**2 + delta**2)
            v_m = delta/norm
            v_p = hbar_omega_n/norm

            Z = matmul(kernel_p, v_p)
            Z = Z*PI*temperature*BOLKEV*(-1.0_q)
            Z = 1.0_q - Z/hbar_omega_n

            delta = matmul(kernel_m, v_m)
            delta = delta*PI*temperature*BOLKEV
            delta = delta/Z

            accuracy = abs(delta - accuracy)/maxval(abs(accuracy))
            counter = counter + 1
        end do

        ! outputs
        if (counter .lt. 1000) then
            write (*, *) 'gap value ', temperature, delta(0), Z(0)
        else
            write (*, *) 'At ', temperature, ' K, isotropic Eliasberg gap equation could not be solved.', delta(0), Z(0)
        end if

        deallocate (delta)
        deallocate (Z)
        deallocate (kernel_p)
        deallocate (kernel_m)
        deallocate (v_m)
        deallocate (v_p)
        deallocate (norm)
        deallocate (accuracy)
    end subroutine isotropic_eliashberg_scf

    subroutine eliashberg_kernel_pm(n_max, lambda, hbar_omega_n, mu_star, kernel_p, kernel_m)
        integer, intent(in) :: n_max
        real(q), intent(in) :: lambda(-2*n_max:2*n_max)
        real(q), intent(in) :: hbar_omega_n(-n_max:n_max)
        real(q), intent(in) :: mu_star
        real(q), intent(out) :: kernel_p(-n_max:n_max, -n_max:n_max)
        real(q), intent(out) :: kernel_m(-n_max:n_max, -n_max:n_max)
        integer i1, i2

        do i1 = -n_max, n_max
            do i2 = -n_max, n_max
                kernel_m(i1, i2) = lambda(i1 - i2) - mu_star
                kernel_p(i1, i2) = lambda(i1 - i2)
            end do
        end do

    end subroutine eliashberg_kernel_pm

    !> compute the eliashberg kernel
    !> K(n,n') = [lambda(omega_n-omega_n')-mu*delta_nn' - sum_l lambda(omega_n-omega_l)sign(omega_l)sign(omega_n)delta_nn']/sqrt(|2n+1||2n'+1|)
    !> where omega_n = pi*(2n+1)*kB*T/hbar are the matsubara frequencies
    !> lambda(omega_n-omega_n') = 2 int d omega (omega/(omega^2+(omega_n-omega_n')^2)) a2F(omega)
    subroutine eliashberg_kernel(n_max, lambda, hbar_omega_n, mu_star, kernel)
        integer, intent(in) :: n_max
        real(q), intent(in) :: lambda(-2*n_max:2*n_max)
        real(q), intent(in) :: hbar_omega_n(-n_max:n_max)
        real(q), intent(in) :: mu_star
        real(q), allocatable, intent(out) :: kernel(:, :)
        integer i1, i2
        real(q) sum_l

        allocate (kernel(-n_max:n_max, -n_max:n_max))

        do i1 = -n_max, n_max
            sum_l = 0.0_q
            do i2 = -n_max, n_max
                kernel(i1, i2) = lambda(i1 - i2) - mu_star
                sum_l = sum_l + lambda(i1 - i2)*sign(1.0_q, hbar_omega_n(i1))*sign(1.0_q, hbar_omega_n(i2))
            end do
            kernel(i1, i1) = kernel(i1, i1) - sum_l
        end do

        do i1 = -n_max, n_max
            do i2 = -n_max, n_max
                kernel(i1, i2) = kernel(i1, i2)/sqrt(abs(2.0_q*i1 + 1.0_q)*abs(2.0_q*i2 + 1.0_q))
            end do
        end do

    end subroutine eliashberg_kernel

    !> compute the largest eigenvalue of a real symmetric matrix
    subroutine compute_eigenvalue(kernel, ncol, eigenvalues, eigenvectors, nb_eigen)
        use tutor, only: vtutor
        real(q), intent(in) :: kernel(ncol, ncol)
        integer, intent(in) :: ncol
        real(q), intent(inout) :: eigenvalues(nb_eigen)
        real(q), intent(inout) :: eigenvectors(ncol, nb_eigen)
        integer, intent(in) :: nb_eigen
        integer M ! number of found eigenvalues
        real(q), allocatable :: eigenvalues_(:)
        real(q), allocatable :: WORK(:)
        integer, allocatable :: IWORK(:)
        integer, allocatable :: IFAIL(:)
        integer VL, VU, IL, IU, LWORK, INFO
        REAL(q) ABSTOL
        integer i1

        ABSTOL = 1E-10_q  ! tolerence on eigenvalues
        VL = 0; VU = 0
        IL = ncol - nb_eigen + 1  ! index of smallest eigenvalue
        IU = ncol    ! index of largest eigenvalue
        LWORK = 8*ncol

        allocate (eigenvalues_(1:ncol))
        allocate (WORK(1:LWORK))
        allocate (IWORK(5*ncol))
        allocate (IFAIL(ncol))
        CALL WDSYEVX('V', 'I', 'U', ncol, kernel, ncol, VL, VU, IL, IU, &
                     ABSTOL, M, eigenvalues_, eigenvectors, ncol, WORK, LWORK, IWORK, IFAIL, INFO)
        eigenvalues(1:nb_eigen) = eigenvalues_(1:nb_eigen)
        if (M .ne. nb_eigen) then
            call vtutor%error('Eliashberg equation: eigenvalues not computed')
        end if
        deallocate (eigenvalues_)
        deallocate (WORK)
        deallocate (IWORK)
        deallocate (IFAIL)
    end subroutine compute_eigenvalue

    !> compute electron-phonon coupling function lambda(omega_n-omega_n')
    !> lambda(omega_n-omega_n') = 2 int d omega (omega/(omega^2+(omega_n-omega_n')^2)) a2F(omega)
    subroutine elph_coupling_function(hbar_omega, a2F, hbar_omega_n, n_max, lambda)
        real(q), intent(in) :: hbar_omega(:)
        real(q), intent(in) :: a2F(:)
        real(q), intent(in) :: hbar_omega_n(-n_max:n_max)
        integer, intent(in) :: n_max
        real(q), allocatable, intent(inout) :: lambda(:)
        integer n, np, p, i
        real(q) hbar_omega_nnp

        do n = -n_max, n_max
            do np = -n_max, n_max
                p = n - np
                hbar_omega_nnp = hbar_omega_n(n) - hbar_omega_n(np)
                lambda(p) = 0.0_q
                do i = 2, size(hbar_omega) - 1
                    lambda(p) = lambda(p) + (hbar_omega(i + 1) - hbar_omega(i)) &
                                *((2*hbar_omega(i))/(hbar_omega(i)**2 + hbar_omega_nnp**2))*a2F(i)
                end do
            end do
        end do
    end subroutine elph_coupling_function

    !> set matsubara frequencies
    !> The frequencies are given by hbar*omega_n = pi*(2n+1)*kB*T
    subroutine set_matsubara_frequencies(temperature, emax, hbar_omega_n, n_max)
        use constant
        real(q), intent(in) :: temperature
        real(q), intent(in) :: emax
        real(q), allocatable, intent(out) :: hbar_omega_n(:)
        integer, intent(out) :: n_max
        integer n

        n_max = int(emax/(2.0_q*pi*temperature*BOLKEV) - 0.5_q) + 1
        allocate (hbar_omega_n(-n_max:n_max))

        ! we consider a symmetric interval
        do n = -n_max, n_max
            hbar_omega_n(n) = pi*(2.0_q*n + 1.0_q)*temperature*BOLKEV
        end do

    end subroutine set_matsubara_frequencies

    !> Compute eliashberg function a2F(omega,temperature) from phonon lifetime
    subroutine elph_isotropic_eliashberg_function(self, kpoints, latt_cur, phonon_freqs_ibz, celtot, phonon_lifetimes_ibz, nedos_ph, emax_ph)
        use constant
        use mkpoints, only: kpoints_struct, get_indexes, kpoint_integers, kpoints_init_tetrahedron
        use tet_macdonald, only: weight_24tetra_delta
        use lattice, only: latt
        use tutor, only: vtutor
        use elphon_accumulators, only: phon_tol, thztoev, get_energies_complex, get_energies
        type(sc_iso_eliashberg), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints
        type(latt), intent(in) :: latt_cur
        real(q), intent(in) :: phonon_freqs_ibz(:, :)
        real(q), intent(in) :: phonon_lifetimes_ibz(:, :, : , :)
        complex(q), intent(in) :: celtot(:, :, :)
        integer, intent(in) :: nedos_ph
        real(q), intent(in) :: emax_ph
        integer i, iq_ibz, imode, itemps, ispin
        integer qint(3)
        integer :: idx_q(4, 24)
        real(q) :: pheig_tetra(4, 24)
        real(q) deltaw(1:nedos_ph)
        integer :: kcut24(3, 4, 24) !< displacements of 24 tetrahedra that contribute to a k-point
        integer :: kcut24_origin(24) !< index of the corner of each of the 24 tetrahedra for which the k-point has 0 displacement

        ! output arrays
        do i = 1, nedos_ph
            self%w(i) = emax_ph*(i - 1)/(nedos_ph - 1)! in eV
        end do

        call kpoints_init_tetrahedron(kpoints,kcut24,kcut24_origin)
        do iq_ibz = 1, kpoints%nkpts
            call KPOINT_INTEGERS(kpoints, latt_cur, kpoints%vkpt(:, iq_ibz), qint)
            call get_indexes(kpoints, qint, kcut24, idx_q)
            do imode = 1, self%natoms*3
                call get_energies(idx_q, imode, phonon_freqs_ibz, pheig_tetra)
                pheig_tetra = pheig_tetra*thztoev
                deltaw = weight_24tetra_delta(kpoints%volwgt, kcut24_origin, nedos_ph, self%w, pheig_tetra)
                do itemps = 1, self%ntemps
                do ispin = 1, self%nspin
                    self%a2F(:, itemps, ispin) = self%a2F(:, itemps, ispin) &
                                                 + deltaw(:)*kpoints%wtkpt(iq_ibz)*kpoints%nkfull &
                                                 *phonon_lifetimes_ibz(itemps, iq_ibz, imode, ispin) &
                                                   /(phonon_freqs_ibz(imode, iq_ibz)*THZTOEV)
                end do
                end do
                self%phonon_dos(:) = self%phonon_dos(:) + deltaw(:)*kpoints%wtkpt(iq_ibz)*kpoints%nkfull
            end do
        end do

        ! dos at chemical potential
        call elph_electron_DOS_at_mu(self, kpoints, latt_cur, celtot)

        ! isotropic eliashberg function
        do itemps = 1, self%ntemps
            do ispin = 1, self%nspin

                if (self%nspin == 1 .and. self%spindeg == 2) then   ! non spin polarized
                    ! 2 times the phonon lifetime, and 2 times the electron dos
                    self%a2F(:, itemps, ispin) = 2.0_q*self%a2F(:, itemps, ispin)/(self%DOSEF(itemps, 1)*2.0_q*2.0_q*pi)
                end if
                if (self%nspin == 2 .and. self%spindeg == 1) then   ! spin polarized
                    ! total eliashberg will be the sum of the 2 spin components
                    self%a2F(:, itemps, ispin) = self%a2F(:, itemps, ispin) &
                                                 /((self%DOSEF(itemps, 1) + self%DOSEF(itemps, 2))*2.0_q*pi)
                end if
                if (self%nspin == 1 .and. self%spindeg == 1) then   ! non collinear
                    self%a2F(:, itemps, ispin) = self%a2F(:, itemps, ispin)/(self%DOSEF(itemps, 1)*2.0_q*pi)
                end if

            end do
        end do

        self%w = self%w/thztoev ! put back in THz
        self%phonon_dos(:) = self%phonon_dos(:)*thztoev ! change to 1/THz
    end subroutine elph_isotropic_eliashberg_function

    !> compute e-dos a thte chemical potential
    subroutine elph_electron_DOS_at_mu(self, kpoints, latt_cur, celtot)
        use constant
        use mkpoints, only: kpoints_struct, get_indexes, kpoint_integers, kpoints_init_tetrahedron
        use tet_macdonald, only: weight_24tetra_delta
        use lattice, only: latt
        use tutor, only: vtutor
        use mathtools, only: det_3x3
        use elphon_accumulators, only: phon_tol, thztoev, get_energies_complex, get_energies
        type(sc_iso_eliashberg), intent(inout) :: self
        type(kpoints_struct), intent(in) :: kpoints
        type(latt), intent(in) :: latt_cur ! elph_mels_pw%latt_cur
        complex(q), intent(in) :: celtot(:, :, :) !< electron eigenvalues for all bands in the IBZ ! w_dense_ibz%celtot(:, :, :)
        integer i, ispin, ik_ibz, ibk, nw
        integer kint(3)
        integer :: idx_k(4, 24)
        real(q) :: eleig_tetra_k(4, 24)
        real(q) w(1:size(self%efermi(:)))
        real(q) deltaw(1:size(self%efermi(:)))
        integer :: kcut24(3, 4, 24) !< displacements of 24 tetrahedra that contribute to a k-point
        integer :: kcut24_origin(24) !< index of the corner of each of the 24 tetrahedra for which the k-point has 0 displacement

        ! output arrays
        nw = size(self%efermi(:))
        w(:) = self%efermi(:)

        call kpoints_init_tetrahedron(kpoints,kcut24,kcut24_origin)
        do ispin = 1, self%nspin
            do ik_ibz = 1, kpoints%nkpts
                call KPOINT_INTEGERS(kpoints, latt_cur, kpoints%vkpt(:, ik_ibz), kint)
                call get_indexes(kpoints, kint, kcut24, idx_k)
                do ibk = 1, size(celtot(:, 1, 1))
                    call get_energies_complex(idx_k, ibk, celtot(:, :, ispin), eleig_tetra_k)
                    deltaw = weight_24tetra_delta(kpoints%volwgt, kcut24_origin, nw, w, eleig_tetra_k)
                    self%DOSEF(:, ispin) = self%DOSEF(:, ispin) + deltaw(:)*kpoints%wtkpt(ik_ibz)*kpoints%nkfull
                end do
            end do
        end do

    end subroutine elph_electron_DOS_at_mu

    subroutine elph_isotropic_Tc(self, nedos_ph, emax_ph, matsubara_temps, matsubara_ntemps, &
                                 emax_matsubara, mu_star, io, comm_global)
        use constant
        use mpimy, only: communic
        use elphon_accumulators, only: thztoev
        type(sc_iso_eliashberg), intent(inout) :: self(:)
        integer, intent(in) :: nedos_ph
        real(q), intent(in) :: emax_ph
        real(q), allocatable, intent(in) :: matsubara_temps(:)
        integer, intent(in) :: matsubara_ntemps
        real(q), intent(in) :: emax_matsubara
        real(q), intent(in) :: mu_star
        integer, intent(in) :: io
        type(communic), intent(in) :: comm_global
        real(q), allocatable :: hbar_omega_n(:) ! matsubara frequencies in eV
        real(q), allocatable :: lambda(:)
        integer i       ! index for accumulators
        integer itemp   ! index for phonon temperatures
        integer itm     ! index ofr Matsubara temperatures
        integer ispin
        integer n_max
        real(q), allocatable :: kernel(:, :)
        real(q), allocatable :: eigenvalues(:)
        real(q), allocatable :: eigenvectors(:, :)
        real(q), allocatable :: lambda_tmp(:)
        real(q) lambda0, omega_log, omega_log_tmp, omega_max
        integer i1, i2, ncol, nb_eigen
        integer global_index_max, j

        i=1
        global_index_max = size(self)*self(i)%ntemps*matsubara_ntemps
#ifdef MPI
        do j = 1 + (comm_global%NODE_ME - 1), global_index_max, comm_global%NCPU
#else
        do j = 1, global_index_max
#endif
            i = (j - 1)/(self(i)%ntemps*matsubara_ntemps) + 1
            itemp = mod(j - 1, self(i)%ntemps*matsubara_ntemps)/matsubara_ntemps + 1
            itm = mod(j - 1, matsubara_ntemps) + 1

            ! log average of the phonon spectrum (in eV)
            omega_log = 1.0_q
            do ispin = 1, self(i)%nspin
                call get_omega_log(self(i)%w(:)*thztoev, self(i)%a2F(:, itemp, ispin), lambda0, omega_log_tmp)
                omega_log = omega_log*omega_log_tmp
            end do

            ! max of the phonon spectrum (i eV)
            omega_max = maxval(self(i)%w(:)*thztoev)

            ! set matsubara frequencies
            call set_matsubara_frequencies(matsubara_temps(itm), emax_matsubara, hbar_omega_n, n_max)

            ! electron-phonon coupling function
            allocate (lambda(-2*n_max:2*n_max))
            allocate (lambda_tmp(-2*n_max:2*n_max))
            lambda = 0.0_q
            do ispin = 1, self(i)%nspin
                lambda_tmp = 0.0_q
                call elph_coupling_function(self(i)%w(:)*thztoev, self(i)%a2F(:, itemp, ispin), &
                                            hbar_omega_n, n_max, lambda_tmp)
                lambda = lambda + lambda_tmp
            end do
            deallocate (lambda_tmp)

            ! electron-phonon coupling constant
            do i1 = -n_max, n_max
            do i2 = -n_max, n_max
                if ((i1 .eq. n_max) .and. (i2 .eq. n_max)) then
                    lambda0 = lambda(i1 - i2)
                    write (*, *) 'electron-phonon coupling constant', lambda0
                end if
            end do
            end do

            ! solve linearized Eliashberg equation
            ncol = 2*n_max + 1
            nb_eigen = 1 ! number of eigenvalues we wish to compute
            allocate (eigenvalues(nb_eigen))
            allocate (eigenvectors(ncol, nb_eigen))
            call eliashberg_kernel(n_max, lambda, hbar_omega_n, mu_star, kernel)
            call compute_eigenvalue(kernel, ncol, eigenvalues, eigenvectors, nb_eigen)
            write (*, *) 'eigenvalues', matsubara_temps(itm), eigenvalues
            !write(*,*) 'eigenvectors',eigenvectors
            deallocate (kernel)
            deallocate (eigenvalues)
            deallocate (eigenvectors)

            ! solve Eliashberg gap equation
            call isotropic_eliashberg_scf(matsubara_temps(itm), n_max, lambda, hbar_omega_n, &
                                            mu_star, lambda0, omega_log, omega_max)

            deallocate (lambda)
            deallocate (hbar_omega_n)
        end do
    end subroutine elph_isotropic_Tc

#ifdef VASP_HDF5
    subroutine elph_write_iso_eliashberg_hdf5(self, fileid, idx, group, subgroup)
        use vhdf5
        use string, only: str
        type(sc_iso_eliashberg) :: self
        integer, intent(in) :: idx
        integer(HID_T), intent(in) :: fileid
        character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
        character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
        ! local variables
        integer(HID_T) :: groupid, subgroupid
        character(len=MAX_LEN_GROUP) :: my_group
        character(len=MAX_LEN_GROUP) :: my_subgroup

        my_group = GRP_RESULTS; if (present(group)) my_group = trim(group)
        my_subgroup = GRP_ELPHON; if (present(subgroup)) my_subgroup = trim(subgroup)

        ! open group
        VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group)//'/'//trim(my_subgroup)//'/superconductivity', groupid))
        VH5_CHECK(vh5_group_open_or_create(groupid, 'isotropic_eliashberg_'//str(idx), subgroupid))

        VH5_CHECK(vh5_write(subgroupid, 'frequency_grid', self%w))
        VH5_CHECK(vh5_write(subgroupid, 'a2F', self%a2F))
        VH5_CHECK(vh5_write(subgroupid, 'phonon_dos', self%phonon_dos))
        if (self%nspin == 1 .and. self%spindeg == 2) then   ! non spin polarized
            VH5_CHECK(vh5_write(subgroupid, 'dos_at_ef', self%DOSEF*2))
        else   ! spin polarized and non collinear
            VH5_CHECK(vh5_write(subgroupid, 'dos_at_ef', self%DOSEF))
        end if

        ! close writting
        VH5_CHECK(vh5_group_close_writing(subgroupid))
        VH5_CHECK(vh5_group_close_writing(groupid))
    end subroutine
#endif

#endif

end module elphon_superconductivity
