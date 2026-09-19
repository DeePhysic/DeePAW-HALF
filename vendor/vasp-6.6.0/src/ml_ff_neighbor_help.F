#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

      MODULE NEIGHBOR_HELP
        USE ML_FF_CONSTANT

        IMPLICIT NONE

        CONTAINS

        SUBROUTINE UPDATE_SINGLE_NEAREST_NEIGHBOR_ARRAY( NEAREST_NEIGHBOR_LIST, &
                                                         LNIONS_EST, &        
                                                         NNEIB_EST, &         
                                                         NORMALIZATION_XYZ, & 
                                                         XYZ_EST, & 
                                                         XYZ_EST_NORM,&
                                                         INTYP0,&
                                                         INIONS )
           USE NEIGHBOR_LIST_MOD, ONLY : NEIGHBOR_LIST_TYPE
           IMPLICIT NONE
           !> nearest neighbor object for descriptor
           TYPE( NEIGHBOR_LIST_TYPE ), INTENT( INOUT ) :: NEAREST_NEIGHBOR_LIST
           !> (index pointing to POSION, type index, INIONS )
           INTEGER, INTENT( INOUT ) :: LNIONS_EST( :, : )
           !> number nearest neighbours , type index, INIONS )
           INTEGER, INTENT( INOUT ) :: NNEIB_EST( : )
           !> norm of connection vector  ( neighbor number , type index, INIONS )
           REAL( q ), INTENT( INOUT ) :: NORMALIZATION_XYZ( :, : )
           !> connection vector  ( neighbor number , xyz, type index, INIONS )
           REAL( q ), INTENT( INOUT ) :: XYZ_EST( :, : , : )
           !> normalized connection vector  ( neighbor number , xyz, type index, INIONS )
           REAL( q ), INTENT( INOUT ) :: XYZ_EST_NORM( :, : , : )
           !> ion type of central ion
           INTEGER, INTENT( IN ) :: INTYP0
           !> ion number in POSION array
           INTEGER, INTENT( IN ) :: INIONS


           ! index of nearest neighbor
           INTEGER :: INN
           INTEGER :: JNIONS
           INTEGER :: JNTYP0


           LNIONS_EST      = 0
           NNEIB_EST       = 0
           XYZ_EST         = 0.0_q
           XYZ_EST_NORM    = 0.0_q


           LNIONS_EST( 1, INTYP0 )        =  INIONS
           XYZ_EST( 1, 1:3, INTYP0 )      =  0.0_q
           XYZ_EST_NORM( 1, 1:3, INTYP0 ) =  0.0_q

           DO INN =1 , SIZE( NEAREST_NEIGHBOR_LIST%NEIGHBORS )
              JNIONS = NEAREST_NEIGHBOR_LIST%NEIGHBORS( INN )
              JNTYP0 = NEAREST_NEIGHBOR_LIST%NEIGHBORS_TYPES( INN )
              ! update nearest neighbor list 
              NNEIB_EST( JNTYP0 ) = NNEIB_EST( JNTYP0 ) + 1
              LNIONS_EST( NNEIB_EST( JNTYP0 ) + 1, JNTYP0 ) = JNIONS
              ! update norm of connection vector
              NORMALIZATION_XYZ( NNEIB_EST( JNTYP0 ) + 1, JNTYP0 ) = &
                       NEAREST_NEIGHBOR_LIST%NEIGHBOR_DIST( INN )
              
              ! update normed connection vector
              XYZ_EST_NORM(NNEIB_EST( JNTYP0 ) + 1, 1, JNTYP0 )= &
                           -NEAREST_NEIGHBOR_LIST%NEIGHBOR_CONNECT( 1 , INN )
              XYZ_EST_NORM(NNEIB_EST( JNTYP0 ) + 1, 2, JNTYP0 )= &
                           -NEAREST_NEIGHBOR_LIST%NEIGHBOR_CONNECT( 2 , INN )
              XYZ_EST_NORM(NNEIB_EST( JNTYP0 ) + 1, 3, JNTYP0 )= &
                           -NEAREST_NEIGHBOR_LIST%NEIGHBOR_CONNECT( 3 , INN )
              
              ! update connection vector
              XYZ_EST( NNEIB_EST( JNTYP0 ) + 1, 1, JNTYP0) = &
                        -NEAREST_NEIGHBOR_LIST%NEIGHBOR_CONNECT( 1 , INN ) * &
                               NEAREST_NEIGHBOR_LIST%NEIGHBOR_DIST( INN )
              XYZ_EST( NNEIB_EST( JNTYP0 ) + 1, 2, JNTYP0) = &
                        -NEAREST_NEIGHBOR_LIST%NEIGHBOR_CONNECT( 2 , INN ) * &
                               NEAREST_NEIGHBOR_LIST%NEIGHBOR_DIST( INN )
              XYZ_EST( NNEIB_EST( JNTYP0 ) + 1, 3, JNTYP0) = &
                        -NEAREST_NEIGHBOR_LIST%NEIGHBOR_CONNECT( 3 , INN ) * &
                               NEAREST_NEIGHBOR_LIST%NEIGHBOR_DIST( INN )
           END DO

        END SUBROUTINE UPDATE_SINGLE_NEAREST_NEIGHBOR_ARRAY
      END MODULE NEIGHBOR_HELP
#endif
