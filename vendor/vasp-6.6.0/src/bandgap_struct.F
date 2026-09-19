#include "symbol.inc"
module bandgap_struct

    use prec, only: q

    implicit none

    integer, parameter :: band_not_found = -1

    enum, bind(c)
        enumerator :: use_fermi_weights = 1, use_local_fermi_energy
    end enum

    !> Store the properties of an orbital like the eigenvalue, the indices of band, kpoint, and spin
    type orbital_property
        real(q) :: eigenvalue  !< Kohn-Sham eigenvalue of the orbital
        real(q) :: kvector(3) = 0.0_q  !< k-point coordinates of the orbital
        integer :: band = band_not_found  !< band index of the orbital
        integer :: kpoint = band_not_found  !< kpoint index of the orbital
        integer :: spin = band_not_found  !< spin index of the orbital
    end type orbital_property


    !> Store information about the bandgap of the system
    type bandgap_info
        real(q) fermi_energy
             !< Fermi energy of the system.
        type(orbital_property) fundamental_valence, fundamental_conduction
             !< orbitals that form the fundamental gap
        real(q) fundamental_gap
             !< size of the direct band gap
        type(orbital_property) direct_valence, direct_conduction
             !< orbitals that for the direct gap
        real(q) direct_gap
             !< size of the direct band gap
    end type bandgap_info

end module bandgap_struct

