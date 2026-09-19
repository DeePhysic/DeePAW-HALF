#include "symbol.inc"
MODULE hamil_struct_def
  USE prec
  USE pot_struct_def, ONLY : coulomb_potential
  USE setexm_struct_def, ONLY : xc_info

  ! handle for orbital magnetization and metaGGAs
  TYPE ham_handle
     TYPE(coulomb_potential) :: COULOMB_POT
     TYPE(xc_info) :: XC
     COMPLEX(q), POINTER, CONTIGUOUS :: AVTOT(:,:)    => NULL() ! local vector magnetization potential
     RGRID     , POINTER, CONTIGUOUS :: AVEC(:,:)     => NULL() ! soft part of vector magnetization potential
     COMPLEX(q), POINTER, CONTIGUOUS :: MUTOT(:,:)    => NULL() ! derivative of energy density with respect to kinetic energy density
     RGRID     , POINTER, CONTIGUOUS :: MU(:,:)       => NULL() ! same as MUTOT, but on GRID instead of GRIDC
     RGRID     , POINTER, CONTIGUOUS :: MUJVEC(:,:,:) => NULL() ! same as MUJVECTOT, but on GRID instead of GRIDC
  END TYPE ham_handle

END MODULE hamil_struct_def
