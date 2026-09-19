#include "symbol.inc"
MODULE locproj_struct
  USE prec

!
!> Data structure that stores the definitions of the local functions
!> on which we will project our Bloch orbitals.
!
  TYPE LPRJ_function
     ! radial part
     INTEGER radial_type    !< 1="optimal" PAW projector, 2=PAW partial wave (PS), 3=Hydrogen-like (=wannier90)
     INTEGER species        !< index of the species at which the local orbital is located
     INTEGER n              !< principal qunatum number
     REAL(q) za             !< real number determining how large the radius of the hydrogenic orbitals should be
     LOGICAL WANNIER90_ORBITAL_DEFINITIONS !< logical governing wether to use wannier90 or LOCPROJ orbital definitions
     ! spherical part
     INTEGER l              !< l angular quantum number
     INTEGER m              !< m magnetic quantum number
     LOGICAL LROTYLM        !< logical to rotate the radial part of the orbitals
     REAL(q) proj_x(3)      !< x axis to use for the rotated projections
     REAL(q) proj_z(3)      !< z axis to use for the rotated projections
     ! site part
     INTEGER poscar_site    !< index of the position in the poscar file
     REAL(q) R(3)           !< position of the localized orbital in reduced coordinates
     ! spin part
     INTEGER spinor         !< spinor index 1 for spin up and 2 for spin down
     REAL(q) spin_qaxis(3)  !< spin quantization axis
     ! related to PS or PR projectors
     INTEGER iproj          !< index of projector or partial wave to be projected on
     INTEGER ibase          !< first pseudopotential channel with angular moment L=l
     INTEGER nproj          !< number of projectors with same L

     !> contains the optimized projectors computed by SPHPRO or SPHPRO_OPTPROJ
     !> and stored in the pseudo datastructure
     REAL(q), POINTER :: optproj(:)

     !> spline representation of the radial part
     REAL(q), ALLOCATABLE :: SPLINE(:,:,:)
  END TYPE LPRJ_function
END MODULE locproj_struct
