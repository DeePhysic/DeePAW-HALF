#include "symbol.inc"
MODULE bse_struct
  USE prec
  IMPLICIT NONE


  INTEGER :: BSE_DESC(9)         ! descriptor for the BSE Hamiltonian
  INTEGER :: VEC_DESC(9)         ! descriptor for non-distributed vector
  INTEGER :: DVEC_DESC(9)        ! descriptor for distributed vector

  TYPE banddesc
    INTEGER :: VBMIN, VBMAX      ! band index for the VB minimum and VB maximum
    INTEGER :: CBMIN, CBMAX      ! band index for the CB minimum and CB maximum
    INTEGER :: CBMIN4, CBMAX4    ! band index for the CB minimum and CB maximum for n4
                                 ! this might be reduced for parallelization over bands
    INTEGER :: NGLB              ! block size for conduction bands = total number of conduction bands
    INTEGER :: NGLB4             ! block size for conduction bands for fourth index (N4 usually)
                                 ! [ N4 is distributed over nodes when band parallelization is used
                                 !   in this case NGLB4 is equal NGLB/ number of nodes]
    INTEGER :: NCORE             ! number of core states
    INTEGER :: VBMAXpCORE        ! VBMAX plus NCORE
    LOGICAL, POINTER:: IS_CORE(:) ! given an index between VBMIN and VBMAXpCORE, logical returns .TRUE. if core orbital
  END TYPE banddesc

  TYPE bse_matrix_index
    INTEGER :: NCV                  !< total rank
    INTEGER, ALLOCATABLE ::  N1(:)  !< valence band index
    INTEGER, ALLOCATABLE ::  N3(:)  !< conduction band index
    INTEGER, ALLOCATABLE ::  NK(:)  !< k-point index
    INTEGER, ALLOCATABLE ::  ISP(:) !< spin index
    !< index into Hamilton matrix for given valence, conduction, k-point and spin
    INTEGER, ALLOCATABLE ::  INDEX(:,:,:,:)

    ! used to index full matrix
    ! column indices K1,K2,ISP1,ISP2,NPOS1, and NPOS3
    ! index along one axis columns
    INTEGER :: COLUMNS
    INTEGER, ALLOCATABLE ::  CK1(:)
    INTEGER, ALLOCATABLE ::  CISP1(:)
    INTEGER, ALLOCATABLE ::  CNPOS1(:)
    INTEGER, ALLOCATABLE ::  CNPOS3(:)

    ! index in full matrix
    ! block indices K1,K2,ISP1,ISP2,NPOS1, and NPOS2
    INTEGER :: BLOCKS
    INTEGER, ALLOCATABLE ::  BK1(:)
    INTEGER, ALLOCATABLE ::  BK2(:)
    INTEGER, ALLOCATABLE ::  BISP1(:)
    INTEGER, ALLOCATABLE ::  BISP2(:)
    INTEGER, ALLOCATABLE ::  BNPOS1(:)
    INTEGER, ALLOCATABLE ::  BNPOS2(:)
    INTEGER, ALLOCATABLE ::  BNPOS3(:)
    INTEGER, ALLOCATABLE ::  BNPOS4(:)

    ! map indices from compact triangular to full ractangular matrix
    ! this way we can easily find starting and ending position for mask
    INTEGER, ALLOCATABLE ::  BMAP(:)

    ! mask empty rows/columns
    INTEGER, ALLOCATABLE ::  MASK(:,:)
  END TYPE bse_matrix_index


  TYPE selfenergy_from_bse
     INTEGER :: NOMEGA
     INTEGER :: NBANDS
     INTEGER :: NKPTS
     INTEGER :: ISPIN
     REAL(q), POINTER :: OMEGA(:)  => NULL()
     REAL(q), POINTER :: REAL_PART(:,:,:,:) => NULL() !1st index over frequencies, 2nd index over bands
     REAL(q), POINTER :: IMAG_PART(:,:,:,:) => NULL()!3rd index over kpoints, 4th index over spin
     COMPLEX(q), POINTER :: CELNEW(:,:,:) => NULL()  !bands, kpoints, spin
  END TYPE selfenergy_from_bse

END MODULE bse_struct
