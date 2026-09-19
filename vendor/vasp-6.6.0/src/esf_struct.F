#include "symbol.inc"

MODULE esf_struct_def

  USE prec
  USE poscar_struct_def, ONLY: latt
  USE mkpoints_struct_def, ONLY: skpoints_trans 
  USE wave_struct_def, ONLY: wavedes
  IMPLICIT NONE

  ! handle for electronic structure factor 
  TYPE elstufac_handle
     !> number of k-points in x-direction
     INTEGER                 :: NKPX = 0
     !> number of k-points in y-direction
     INTEGER                 :: NKPY = 0
     !> number of k-points in z-direction
     INTEGER                 :: NKPZ = 0
     !> wavefunction descriptor of response function
     TYPE (wavedes), POINTER ::  WGW => NULL()
     !> number of irreducible k-points in 1. BZ
     INTEGER                 :: NUMBER_OF_NQ = 0 
     !> is symmetry switched on or off 
     INTEGER                 :: ISYM = -1
     !> kpoint map from IRZ to full BZ
     TYPE( skpoints_trans )  :: KPOINTS_TRANS
     !> number of G-vectors 
     INTEGER, ALLOCATABLE    :: NG(:)         
     !> electronic structure factor 
     !> similar to diagonal of frequency-integrated responsefunction
     !> S(q+G)=\int dw xi(q+G,q+G,w)
     !> 2nd and 3rd index correspond to cutoff and lambda value considered
     !> electronic stucture factor in q+G domain 
     REAL(q),ALLOCATABLE     :: S(:,:)

     !> used q-points, incase NKRED is set
     INTEGER, ALLOCATABLE    :: NQ(:)         

     !> energy cutoff 
     REAL(q)                 :: ENCUT=0._q
     !> soft energy cutoff 
     REAL(q)                 :: ENCUTSOFT =0._q
     !> number of G-vectors in 1st direction  
     INTEGER                 :: N1 = 0         
     !> number of G-vectors in 2nd direction  
     INTEGER                 :: N2 = 0    
     !> number of G-vectors in 3rd direction  
     INTEGER                 :: N3 = 0      
     !> cell volume
     REAL(q)                 :: OMEGA=0._q
     !> reciprocal Bravais matrix 
     REAL(q)                 :: B(3,3)
     !> Bravais matrix 
     REAL(q)                 :: A(3,3)
     !> cutoff 
     !> splined result
     REAL(q)                 :: ESPLINED=0._q
     !> integer error code
     INTEGER                 :: ISTATUS = -1

     TYPE(latt)              :: LATT_CUR       !< stores lattice information

     !> singularity correction for exchange potential requires number of electrons
     REAL(q)                 :: NELECT = -1._q 
 
     !> auxillary arrays that probably should be removed 
     REAL(q),ALLOCATABLE     :: V(:)           !< Coulomb potential in PW basis 
     INTEGER,ALLOCATABLE     :: OLD_INDEX(:)   !< stores old index in PW array that corresponds to NE=1 cutoff  
  END TYPE elstufac_handle 

END MODULE esf_struct_def 
