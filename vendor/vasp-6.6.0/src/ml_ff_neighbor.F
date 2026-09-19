#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

module neighbor_list_mod

   use ML_FF_PREC
   implicit none
   public :: neighbor_list_type


   type neighbor_list_type
      private
      !> neighbor list
      integer, public, allocatable ::neighbors( : )
      !> type of neighbor
      !integer, public, allocatable ::neighbors_types( : )
      integer, public, allocatable ::neighbors_types( : )
      !> distances between neighbor pair
      real( kind = q ), public, allocatable     ::neighbor_dist( : )
      !> connection vector between neighbor pair
      real( kind = q ), public, allocatable     ::neighbor_connect( : , : )
      !> cutoff radius
      real( kind = q )  ::cutoff
      !> cutoff radius suqared
      real( kind = q )  ::cutoff_squared
      !> index of owner of list, central atom
      integer  ::my_index
      !> number of different atom types
      integer :: number_types
      !> maximum number of neighbors per atom type
      integer, public, allocatable :: max_neighbor( : )

      contains
      !> nearest neighbor N^2 algorithm
      procedure, pass( this )  :: init_neighbor_list => neighbor_list_init_n_squared
      !> compute nearest neighbor list with the N^2 algorithm
      procedure, public, pass( this ) :: compute_nn_list => &
                                             create_neighbour_list_single_atom_n_squared
      !> create a cutoff list from a larger one
      procedure, public, pass( this ) :: create_subset_cutoff_list
      !> return cutoff squared
      procedure, public, pass( this ) :: return_cutoff_squared
      !> updating the data arrays of the nearest-neighbor instance
      procedure, public, pass( this ) ::update_data_arrays
      !! *************************************
   end type neighbor_list_type

   contains

   !> initialize neighbor list in the N^2 mode
   !> this does a minimal init where
   !> as input the cutofff radius is needed and
   !> the index of the central atom in the atom
   !> position list
   !>
   subroutine neighbor_list_init_n_squared( this, cutoff, my_index, ntypes )

      class( neighbor_list_type ),intent( inout ) ::this
      real( kind = q ),intent( in ) ::cutoff
      integer,intent( in )  :: my_index
      integer, intent( in ) :: ntypes

      this%cutoff         =  cutoff
      this%cutoff_squared =  cutoff * cutoff
      this%my_index       =  my_index
      this%number_types   =  ntypes
      ! give dummy sizes to arrays
      allocate( this%neighbors( 1:1 ) )
      allocate( this%neighbors_types(1:1) )
      allocate( this%neighbor_dist( 1:1 ) )
      allocate( this%neighbor_connect( 1:1, 1:1 ) )
      allocate( this%max_neighbor( 1:ntypes ) )
   end subroutine neighbor_list_init_n_squared
   

   !> compute nearesr neigbor list with a N^2 algorithm
   !> input:
   !> positions ( 3,N ) real array
   !> lattice   ( 3,3 ) real array
   !> types     ( 1,N ) character array
   subroutine create_neighbour_list_single_atom_n_squared( this ,&
                                  positions, lattice, types, inverse_lattice )

      class( neighbor_list_type ),intent( inout ) :: this
      real( kind = q ),intent( in ) :: positions( : , : )
      real( kind = q ),intent( in ) :: lattice( : , : )
      integer,intent( in )  :: types( : )
      real( kind = q ),intent( in ), optional :: inverse_lattice( : , : )

      !! loop counters
      integer ::nindx  ! neighbor index
      integer ::x      ! x-shift
      integer ::y      ! y-shift
      integer ::z      ! z-shift

      ! positions arrays
      real( kind = q ) ::atom(1:3)
      real( kind = q ) ::neighbor(1:3)
      ! cartesian distance components
      real( kind = q ) ::dx
      real( kind = q ) ::dy
      real( kind = q ) ::dz
      ! squared distance
      real( kind = q ) ::dist

      ! counters
      integer :: n_counter

      ! temprorary lists
      integer,allocatable ::n_list( : )   ! neighbor index
      integer,allocatable ::n_types( : ) ! neighbor atom type
      real( kind = q ),allocatable ::n_dist( : )       ! neighbor distance
      real( kind = q ),allocatable ::n_vectors( :, : ) ! neighbor connection unit vector
      real( kind = q ) ::eps

      ! determines how far atoms are shifted
      integer ::xyz_max( 1:3 )
      integer ::xyz_min( 1:3 )

      Integer :: multi


      atom = positions( : , this%my_index )
      eps = 1e-6_q


      if ( present( inverse_lattice ) ) then
         call set_periodic_size_n_square( this%cutoff, atom, inverse_lattice ,&
                                          xyz_max, xyz_min )
      else
         xyz_min( : )  =  -1
         xyz_max( : )  =   1
      end if

      ! determine maximum number of periodic images within 
      ! cut-off radius plus 1 (current cell)
      multi=(xyz_max(1)-xyz_min(1)+1)*(xyz_max(2)-xyz_min(2)+1)*(xyz_max(3)-xyz_min(3)+1)

      call allocate_working_arrays_for_data_computation(  n_list, n_types, &
                       n_dist, n_vectors,size( positions,2 )*multi )


      n_counter = 0
      this%max_neighbor = 0
      do x = xyz_min( 1 ), xyz_max( 1 )
         do y = xyz_min( 2 ), xyz_max( 2 )
            do z = xyz_min( 3 ), xyz_max( 3 )
               do nindx = 1, Size( positions,2 )
                  neighbor(1)  =  positions( 1 , nindx )
                  neighbor(2)  =  positions( 2 , nindx )
                  neighbor(3)  =  positions( 3 , nindx )
                  neighbor(1)  =  neighbor(1) + real( x, kind=q ) * lattice( 1, 1 ) + &
                                                real( y, kind=q ) * lattice( 1, 2 ) + &
                                                real( z, kind=q ) * lattice( 1, 3 )
                  neighbor(2)  =  neighbor(2) + real( x, kind=q ) * lattice( 2, 1 ) + &
                                                real( y, kind=q ) * lattice( 2, 2 ) + &
                                                real( z, kind=q ) * lattice( 2, 3 )
                  neighbor(3)  =  neighbor(3) + real( x, kind=q ) * lattice( 3, 1 ) + &
                                                real( y, kind=q ) * lattice( 3, 2 ) + &
                                                real( z, kind=q ) * lattice( 3, 3 )
                  dx = atom(1) - neighbor(1)
                  dy = atom(2) - neighbor(2)
                  dz = atom(3) - neighbor(3)
                  dist = dx*dx + dy*dy + dz*dz
                  if ( dist .le. this%cutoff_squared .and. dist .gt. eps ) then
                     n_counter = n_counter + 1
                     n_list( n_counter )  =  nindx
                     n_types( n_counter ) =  types( nindx )
                     n_dist( n_counter )  =  dsqrt( dist )
                     n_vectors( 1 , n_counter ) = dx / n_dist( n_counter )
                     n_vectors( 2 , n_counter ) = dy / n_dist( n_counter )
                     n_vectors( 3 , n_counter ) = dz / n_dist( n_counter )
                     this%max_neighbor( types( nindx ) ) = &
                                   this%max_neighbor( types( nindx ) ) + 1
                  end if
               end do
            end do
         end do
      end do


      call this%update_data_arrays( n_list, n_types, n_dist, n_vectors, n_counter )

   end subroutine create_neighbour_list_single_atom_n_squared




   subroutine update_data_arrays( this, n_list, n_types, n_dist, n_vectors, n_counter )
      
      class( neighbor_list_type ),intent( inout ) :: this
      
      integer,allocatable,intent( inout ) ::n_list( : )   ! neighbor index
      integer,intent(inout),allocatable ::n_types( : ) ! neighbor atom type
      real( kind = q ),intent( inout ),allocatable ::n_dist( : )       ! neighbor distance
      real( kind = q ),intent( inout ),allocatable ::n_vectors( :, : ) ! neighbor connection unit vector
      integer,intent( in ) ::n_counter



      ! organizing array sizes
      if ( n_counter .ne. size( this%neighbors ) ) then
         deallocate( this%neighbors )
         allocate( this%neighbors( 1:n_counter ) )
         ! if size does not fit for one array it won't fit for any
         deallocate( this%neighbors_types )
         allocate( this%neighbors_types( 1:n_counter ) )
         deallocate( this%neighbor_dist )
         allocate( this%neighbor_dist( 1:n_counter ) )
         deallocate( this%neighbor_connect )
         allocate( this%neighbor_connect( 1:3, 1:n_counter ) )
      end if
      ! copy data to public arrays for usage
      this%neighbors        = n_list( 1:n_counter )
      this%neighbors_types  = n_types( 1:n_counter )
      this%neighbor_dist    = n_dist( 1:n_counter )
      this%neighbor_connect = n_vectors( 1:3, 1:n_counter )

      

      ! deallocate working arrays
      deallocate( n_list )
      deallocate( n_types )
      deallocate( n_dist )
      deallocate( n_vectors )

   end subroutine update_data_arrays

      
   subroutine create_subset_cutoff_list( this, other )

      class( neighbor_list_type ),intent( inout ) :: this
      class( neighbor_list_type ),intent( in ) :: other

      integer :: i
      
      ! counters
      integer :: n_counter
      
      integer,allocatable ::n_list( : )   ! neighbor index
      integer,allocatable ::n_types( : ) ! neighbor atom type
      real( kind = q ),allocatable ::n_dist( : )       ! neighbor distance
      real( kind = q ),allocatable ::n_vectors( :, : ) ! neighbor connection unit vector

      allocate( n_list( 1:size( other%neighbors ) ) )
      allocate( n_types( 1:size( other%neighbors_types ) ) )
      allocate( n_dist( 1:size( other%neighbor_dist ) ) )
      allocate( n_vectors( 1:3, 1:size( other%neighbor_connect ) ) )

      n_counter  =  0
      this%max_neighbor = 0
      do i = 1, size( other%neighbors )
         if ( other%neighbor_dist( i ) .le. this%cutoff ) then
            n_counter  =  n_counter + 1
            n_list( n_counter )  =  other%neighbors( i )
            n_types( n_counter ) =  other%neighbors_types( i ) 
            n_dist( n_counter )  =  other%neighbor_dist( i )
            n_vectors( 1, n_counter ) = other%neighbor_connect( 1, i )
            n_vectors( 2, n_counter ) = other%neighbor_connect( 2, i )
            n_vectors( 3, n_counter ) = other%neighbor_connect( 3, i )
            this%max_neighbor( n_types( n_counter ) ) = &
                         this%max_neighbor( n_types( n_counter ) ) + 1
         end if
      end do

      call this%update_data_arrays( n_list, n_types, n_dist, n_vectors, n_counter )

   end subroutine create_subset_cutoff_list
      


   subroutine allocate_working_arrays_for_data_computation( n_list, n_types, n_dist, n_vectors, n_size )

      integer,allocatable,intent( inout ) ::n_list( : )   ! neighbor index
      integer,intent(inout),allocatable ::n_types( : ) ! neighbor atom type
      real( kind = q ),intent( inout ),allocatable ::n_dist( : )       ! neighbor distance
      real( kind = q ),intent( inout ),allocatable ::n_vectors( :, : ) ! neighbor connection unit vector
      integer,intent( in ) ::n_size
      
      allocate( n_list( 1:n_size ) )
      allocate( n_types( 1:n_size ) )
      allocate( n_dist( 1:n_size ) )
      allocate( n_vectors( 1:3, 1:n_size ) )

   end subroutine allocate_working_arrays_for_data_computation


   function return_cutoff_squared( this ) result( res )
      
      class( neighbor_list_type ),intent( inout ) :: this
      real( kind = q ) :: res
      res = this%cutoff_squared

   end function return_cutoff_squared



   subroutine set_periodic_size_n_square( cutoff, pos, inverse_lattice,&
                                          imax, imin )

      integer, intent( inout ) :: imax( 1:3 )
      integer, intent( inout ) :: imin( 1:3 )
      real( kind = q ),intent( in ) :: cutoff
      real( kind = q ),intent( in ) :: pos( 1:3 )
      real( kind = q ),intent( in ) :: inverse_lattice( 1:3, 1:3 )
      
      real( kind = q ) :: rr(1:3)
      real( kind = q ) ::xyz_max( 1:3, 1:3 )
      real( kind = q ) ::xyz_min( 1:3, 1:3 )
      integer :: ixyz
      integer :: jxyz
      real( kind = q ),parameter ::pi = 2.0_q * dacos( 0.0_q )

      !! taken from Ferenc, Ryosuke
      rr(1:3) = dsqrt( inverse_lattice( 1, 1:3 )*inverse_lattice( 1, 1:3 ) + &
                       inverse_lattice( 2, 1:3 )*inverse_lattice( 2, 1:3 )+ &
                       inverse_lattice( 3, 1:3 )*inverse_lattice( 3, 1:3 ) )
      do ixyz = 1, 3
         do jxyz = 1, 3
            xyz_max(jxyz, ixyz) = pos(jxyz) + cutoff * &
                                  inverse_lattice( jxyz,ixyz ) / rr(ixyz)
            xyz_min(jxyz, ixyz) = pos(jxyz) - cutoff * &
                                  inverse_lattice( jxyz,ixyz ) / rr(ixyz)
         end do
      end do

      do ixyz = 1, 3
         imax(ixyz) = ceiling( 0.5_q * ( inverse_lattice( 1, ixyz ) * xyz_max( 1, ixyz ) + &
                                         inverse_lattice( 2, ixyz ) * xyz_max( 2, ixyz ) + &
                                         inverse_lattice( 3, ixyz ) * xyz_max( 3, ixyz ) ) / pi )-1
         imin(ixyz) = ceiling( 0.5_q * ( inverse_lattice( 1, ixyz ) * xyz_min( 1, ixyz ) + &
                                         inverse_lattice( 2, ixyz ) * xyz_min( 2, ixyz ) + &
                                         inverse_lattice( 3, ixyz ) * xyz_min( 3, ixyz ) ) / pi )-1
      end do

   end subroutine set_periodic_size_n_square

end module neighbor_list_mod


!  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
!  "          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
!  "                     ~~~~~~~~~~~~~~~~~~~~~~~~~~"  
!  "                               ~~~~~"
!  "                               ~~~~~"
!  "                     ~~~~~~~~~~~~~~~~~~~~~~~~~~"  
!  "          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
!  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"



module neighbor_cell_list_mod
   
   use ML_FF_PREC
   use neighbor_list_mod, only :neighbor_list_type
   implicit none
   public :: neighbor_cell_list_type

   private



   !> nearest neighbor data between boxes
   type box_neighbor_list
      integer,allocatable  :: nn_box_list( : )
      integer,allocatable  :: box_to_atoms( : )
      integer,allocatable  :: box_to_atoms_sorted( : )
   end type box_neighbor_list
   
   ! Neighbor list algorithm
   enum, bind(c)
      enumerator :: NLA_ENUM    = 0
      enumerator :: NLA_CELL    = 1
      enumerator :: NLA_SQUARED = 2
   end enum


   type neighbor_cell_list_type
      private
      !> enum to determine algorithm type
      integer(kind(NLA_ENUM)) :: algo_type
      !> cutoff radius
      real( kind = q ) :: cutoff
      !> length of box in x,y,z direction
      real( kind = q )  :: box_length( 1:3 )
      !> number of boxes in x,y,z direction
      integer    :: number_boxes( 1:3 )
      !> box to replicate to cover whole supercell
      real( kind = q )  :: unit_box( 1:3,1:3 )
      !> positions of boxes
      real( kind = q ), allocatable  :: box_positions( : , : )
      !> list assigning atoms to boxes
      integer, allocatable  :: atom_to_box( : )
      !> data for box neighbor list
      type( box_neighbor_list ), allocatable :: nn_boxes( : )
      !> nearest neighbor list data for atoms
      type( neighbor_list_type ),public,allocatable :: nn_atoms( : )
      !> atom positions sorted accoring to their box index
      real( kind = q ), allocatable :: positions_sorted( : , : )
      !> storing indices for the sorted atom positions
      integer, allocatable :: sorted_indexes( : )
      !> shifts which have to be included for nearest neighbor boxes
      integer, allocatable :: box_shift_indices( :, : )
      !> integer maximum number of nearest neighbors per atom
      integer :: max_neighbors
      !> number of different atom types
      integer :: number_types
      !> maximum number of neighbors per atom type
      integer, allocatable :: max_neighbor( : )
      integer  :: nboxes_x
      integer  :: nboxes_y
      integer  :: nboxes_z

      contains
      private
      !> create the cell list type
      procedure, public, pass( this ) :: init_nn_cell_list => &
                                         neighbor_list_init_cell_list_algo
      !> create or update the box neighbor list and assign atoms
      procedure, public, pass( this ) :: update_box_list
      !> create or update the nearest neighbor list between atoms
      generic, public :: compute_nn_list => update_nearest_neighbor_list_main_all,&
                                                            create_subset_cutoff_list
      !> get maximal number of neighbors for vasp
      procedure, public, pass( this ) :: get_max_number_neighbors


      !! private procedures
      !> create nearest neighbor list from super set
      procedure, pass( this )  :: create_subset_cutoff_list
      !> update the nearest neighbor list for all atoms
      procedure, pass( this )  :: update_nearest_neighbor_list_main_all
      !> initialzing the box neighbor list
      procedure, pass( this )  :: initialize_box_list
      !> create box parameters
      procedure, pass( this ) :: create_box_parameters
      !> create nearest neighbor list between boxes
      procedure, pass( this ) :: create_box_neighbor_list_main
      !> compute which shifts are needed to compute box neighbors
      procedure, pass( this ) :: compute_nearest_box_indices
      !> allocating neighbor list arrays for boxes
      procedure, pass( this ) :: allocate_box_neighbor_list
      !> allocating neighbor list arrays for boxes
      !> return box number from 3D integer position
      procedure, pass( this ) :: return_box_number
      !> create nearest neighbor list for single box 
      procedure, pass( this ) :: get_nearest_neighbours_single_box
      !> assigning atoms to boxes of cell list algorithm
      procedure, pass( this ) :: assign_atoms_to_box
      !> allocate atoms to box array
      procedure, pass( this ) :: allocate_my_boxes
      !> sorting atoms for nearest neighbor list construction
      procedure, pass( this ) :: sort_atoms
      !> allocate arrays for sorted atoms positions
      procedure, pass( this ) :: allocate_arrays_sorting
      !> compute nearest neighbor list for single atom
      procedure, pass( this ) :: compute_nn_list_single_atom => &
                      create_nearest_neighbor_list_single_atom_cell_list_sorted
   end type neighbor_cell_list_type

   contains
   subroutine neighbor_list_init_cell_list_algo( this, cutoff, &
                                                 lattice, positions, ntypes, ncpu, node_me )

      implicit none
      class( neighbor_cell_list_type ),intent( inout ) ::this
      real( q ), intent( in ) :: cutoff
      real( q ), intent( in ) :: lattice( : , : )
      real( q ), intent( in ) :: positions( : , : )
      integer, intent( in )   :: ntypes
      integer, intent( in ), optional   :: ncpu
      integer, intent( in ), optional   :: node_me

      integer :: i_atom


      this%cutoff   =  cutoff
      ! initialize atom nearest neighbor lists
      this%max_neighbors  = 1
      this%algo_type      = 1
      this%number_types   = ntypes
      allocate( this%nn_atoms( 1:size( positions,2 ) ) )
      do i_atom = 1 , size( positions,2 )
         call this%nn_atoms( i_atom )%init_neighbor_list( this%cutoff, i_atom, ntypes )
      end do
      ! ~~~~~~~~~~~~~~~~~~~~~
      if ( allocated( this%max_neighbor ) ) deallocate( this%max_neighbor )
      allocate( this%max_neighbor( ntypes ) )
      this%max_neighbor = 0

      ! construct boxes for given lattice
      PROFILING_START( "UPDATE_BOX_LIST_INIT" )
      if ( present( ncpu ) .and. present( node_me ) ) then
         call this%initialize_box_list( positions, lattice, ncpu, node_me )
      else
         call this%initialize_box_list( positions, lattice )
      end if
      PROFILING_STOP( "UPDATE_BOX_LIST_INIT" )


   end subroutine neighbor_list_init_cell_list_algo


   !> initializing the the boxes and the nearest neighbor-list
   !> between the boxes
   subroutine initialize_box_list( this, positions, lattice, ncpu, node_me )

      class( neighbor_cell_list_type ), intent( inout ) :: this
      real( kind = q ),                 intent( in )    :: positions( : , : )
      real( kind = q ),                 intent( in )    :: lattice( : , : )
      integer, intent( in ), optional  :: ncpu
      integer, intent( in ), optional  :: node_me
 
      integer :: cpu_node
      integer :: n_step
 
 
      if ( present( node_me ) .and. present( ncpu ) ) then
         n_step   =  ncpu
         cpu_node =  node_me
      else  
         n_step   = 1
         cpu_node = 0
      end if
 
      call this%create_box_parameters( lattice )
      
      if ( this%algo_type .eq. 1 ) then
         call this%create_box_neighbor_list_main
         call this%assign_atoms_to_box( positions, cpu_node, n_step )
         call this%sort_atoms( positions )
      end if

   end subroutine initialize_box_list
 
 

   !> initialize neighbor list in the cell mode
   !> this does a init where
   !> as input the cutofff radius is needed and
   !> the index of the central atom in the atom
   !> position list, and the lattice parameters
   !> to generate the boxes. The boxes will only be created
   !> when the first call to init is done. All the other atoms
   !> will use the same cell list
   subroutine create_box_parameters( this, lattice )
      
      class( neighbor_cell_list_type ), intent( inout ) :: this
      real( kind = q ),                 intent( in )    :: lattice( 1:3, 1:3 )
      
      integer ::i
      real( kind = q )  :: box_len   ! real box len in different dirs
      real( kind = q )  :: super_len ! length of supercell in xyz direction

      box_len = this%cutoff * 0.5_q
      

      this%algo_type = 1
      do i = 1, size( lattice, 2 ) 
         super_len  =  dsqrt( lattice( 1, i )*lattice( 1, i ) +&    
                              lattice( 2, i )*lattice( 2, i ) +& 
                              lattice( 3, i )*lattice( 3, i ) )
         this%number_boxes( i ) =  max( floor( super_len / box_len ), 1 )
         this%box_length( i )   =  super_len / &
                                        real( this%number_boxes( i ), kind = q )
         this%unit_box( :, i )  =  lattice( :, i ) / real( this%number_boxes( i ), kind = q )
      end do


      !! count number of neighbor boxes needed
      call number_of_needed_neighbor_boxes( this%cutoff, this%unit_box( :, : ), &
                                            this%nboxes_x, this%nboxes_y, this%nboxes_z )
      if ( this%nboxes_x * 2 .gt. this%number_boxes( 1 ) .or. &
           this%nboxes_y * 2 .gt. this%number_boxes( 2 ) .or. &
           this%nboxes_z * 2 .gt. this%number_boxes( 3 ) ) then
           this%algo_type = 2
      end if

   end subroutine create_box_parameters


   !> compute nearest neighbor boxes for cell list algorithm
   !> needed in the 3 different coordinate axis x,y,z.
   subroutine number_of_needed_neighbor_boxes( cutoff, unit_box, nx, ny, nz )

      !> cutoff radius
      real( kind = q ),     intent( in )  :: cutoff
      !> size of a box with which the supercell is covered
      real( kind = q ),     intent( in )  :: unit_box( 1:3, 1:3 )
      integer,              intent( out ) :: nx
      integer,              intent( out ) :: ny
      integer,              intent( out ) :: nz

      real( kind = q ) :: unit_vector( 1:3, 1:3 )
      real( kind = q ) :: cross( 1:3 )


      integer :: i

      do i = 1, 3
         unit_vector( :, i )  =  unit_box( :, i ) / dsqrt( unit_box( 1, i ) * unit_box( 1, i ) + &
                                                           unit_box( 2, i ) * unit_box( 2, i ) + &
                                                           unit_box( 3, i ) * unit_box( 3, i ) )
      end do

      nx = ceiling( cutoff / abs( dot_product( unit_box( :, 1 ), &
                         crossprod( 3, unit_vector( :, 2 ), unit_vector( :, 3 ) ) ) ) )
      
      ny = ceiling( cutoff / abs( dot_product( unit_box( :, 2 ), &
                         crossprod( 3, unit_vector( :, 1 ), unit_vector( :, 3 ) ) ) ) )
      
      nz = ceiling( cutoff / abs( dot_product( unit_box( :, 3 ), &
                         crossprod( 3, unit_vector( :, 2 ), unit_vector( :, 1 ) ) ) ) )


   end subroutine number_of_needed_neighbor_boxes


     
   function crossprod( ndim, a, b )

      integer, intent( in ) :: ndim
      real( kind = q ), intent( in )  :: a( ndim )
      real( kind = q ), intent( in )  :: b( ndim )
      real( kind = q )                :: crossprod( ndim )

      crossprod( 1 ) = a( 2 ) * b( 3 ) - b( 2 ) * a( 3 )
      crossprod( 2 ) = a( 3 ) * b( 1 ) - b( 3 ) * a( 1 )
      crossprod( 3 ) = a( 1 ) * b( 2 ) - b( 1 ) * a( 2 )
   end function



   subroutine compute_nearest_box_indices( this )
   
      class( neighbor_cell_list_type ), intent( inout ) :: this

      real( kind = q ) :: a( 1:3 )
      real( kind = q ) :: b( 1:3 )
      real( kind = q ) :: c( 1:3 )
      real( kind = q ) :: vector( 1:3 )
      real( kind = q ) :: box_cutoff
      real( kind = q ) :: norm

      integer  :: nx
      integer  :: ny
      integer  :: nz
      integer  :: counter
      integer  :: ndim

      integer,allocatable  :: index_array( :, : )

      real( kind = q ) :: room_diag( 1:3 )


      integer :: i
      integer :: j
      integer :: k

      nx  =  this%nboxes_x
      ny  =  this%nboxes_y
      nz  =  this%nboxes_z

      ndim = ( 2 * nx + 1 ) * ( 2 * ny + 1 ) * ( 2 * nz + 1 ) 


      a  =  this%unit_box( :, 1 )
      b  =  this%unit_box( :, 2 )
      c  =  this%unit_box( :, 3 )

      
      room_diag( : ) = a + b + c
      
      box_cutoff   = ( this%cutoff + dsqrt( dot_product( room_diag, room_diag ) ) / 2 )**2


      allocate( index_array( 1:3, 1:ndim ) )
      counter = 0

      do i = -nx, nx
         do j = -ny, ny
            do k = -nz, nz
               vector( 1 ) =  real( i, kind = q ) * a( 1 ) + &
                              real( j, kind = q ) * b( 1 ) + &
                              real( k, kind = q ) * c( 1 )
               vector( 2 ) =  real( i, kind = q ) * a( 2 ) + &
                              real( j, kind = q ) * b( 2 ) + &
                              real( k, kind = q ) * c( 2 )
               vector( 3 ) =  real( i, kind = q ) * a( 3 ) + &
                              real( j, kind = q ) * b( 3 ) + &
                              real( k, kind = q ) * c( 3 )
               norm = vector( 1 ) * vector( 1 ) + vector( 2 ) * vector( 2 ) + &
                                    vector( 3 ) * vector( 3 )
               if ( norm .le. box_cutoff ) then
                  counter  =  counter + 1
                  index_array( 1, counter )  =  i
                  index_array( 2, counter )  =  j
                  index_array( 3, counter )  =  k
               end if
            end do
         end do
      end do

      if ( allocated( this%box_shift_indices ) ) deallocate( this%box_shift_indices )
      allocate( this%box_shift_indices( 1:3, 1:counter ) )
      this%box_shift_indices( :, : ) = index_array( 1:3, 1:counter )


   end subroutine compute_nearest_box_indices



   !> create a box neighbor list for the cell
   !> list algorithm with linear scaling by
   !> by exploiting the build up principle
   !> of the boxes
   subroutine create_box_neighbor_list_main( this )

      class( neighbor_cell_list_type ), intent( inout ) :: this

      integer ::x_box
      integer ::y_box
      integer ::z_box

      integer :: nx
      integer :: ny
      integer :: nz

      integer ::box_indx


      PROFILING_START( "create_box_neighbor_list_main" )
      nx = this%number_boxes( 1 )
      ny = this%number_boxes( 2 )
      nz = this%number_boxes( 3 )

      call this%compute_nearest_box_indices()
      call this%allocate_box_neighbor_list( size( this%box_shift_indices, 2 ) )
      do box_indx = 1, nx * ny * nz
         z_box = modulo( box_indx + 2*nz-1, nz ) + 1
         y_box = modulo( ceiling( real( box_indx )/real( nz ) ) + 2 * ny - 1, ny ) + 1
         x_box = modulo( ceiling( real( box_indx )/real( nz * ny ) ) + 2 * nx - 1, nx ) + 1
         call this%get_nearest_neighbours_single_box( x_box , y_box, z_box ,&
                                                          box_indx, nx, ny, nz &
                                                        )
      end do
      PROFILING_STOP( "create_box_neighbor_list_main" )

   end subroutine create_box_neighbor_list_main



   !> return the box number when
   !> the position of the box is supplied on the 
   !> integer grid
   function return_box_number( this, x , y , z ) result( indx )

      class( neighbor_cell_list_type ), intent( inout ) :: this
      integer,                          intent( in ) :: x
      integer,                          intent( in ) :: y
      integer,                          intent( in ) :: z

      integer ::indx

      indx = (x-1) * this%number_boxes(2) * &
                     this%number_boxes(3) + &
             (y-1) * this%number_boxes(3) + z
   end function return_box_number
      

   
   !> create the nearest neighbor list for a single box
   !> this routine uses the building up principle of
   !> the boxes
   subroutine get_nearest_neighbours_single_box( this, x, y, z , box_indx,&
                                                 x_box_max, y_box_max, z_box_max )

      class( neighbor_cell_list_type ), intent( inout ) :: this
      integer,                          intent( in ) :: x
      integer,                          intent( in ) :: y
      integer,                          intent( in ) :: z
      integer,                          intent( in ) :: box_indx
      integer,                          intent( in ) :: x_box_max
      integer,                          intent( in ) :: y_box_max
      integer,                          intent( in ) :: z_box_max


      integer ::nx
      integer ::ny
      integer ::nz


      integer ::n_xx
      integer ::n_yy
      integer ::n_zz

      integer :: x_pos
      integer :: y_pos
      integer :: z_pos
      

      integer ::counter

      real( kind = q ) ::dist

      nx  =  this%nboxes_x
      ny  =  this%nboxes_y
      nz  =  this%nboxes_z

      counter = 0
      ! modulos are slower than the if else caluse.
      ! 
      do counter = 1, size( this%box_shift_indices, 2 )
         n_xx  = this%box_shift_indices( 1, counter )
         n_yy  = this%box_shift_indices( 2, counter )
         n_zz  = this%box_shift_indices( 3, counter )
         x_pos = return_box_position_1d( x , x_box_max , n_xx )
         !x_pos  = modulo( x + n_xx + 2 * x_box_max - 1 , x_box_max ) + 1 
         y_pos = return_box_position_1d( y , y_box_max , n_yy )
         !y_pos  = modulo( y + n_yy + 2 * y_box_max - 1 , y_box_max ) + 1 
         z_pos = return_box_position_1d( z , z_box_max , n_zz )
         !z_pos  = modulo( z + n_zz + 2 * z_box_max - 1 , z_box_max ) + 1 
         !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         this%nn_boxes( box_indx )%nn_box_list( counter ) = &
                                           (x_pos-1) * y_box_max * z_box_max + &
                                           (y_pos-1) * z_box_max + z_pos
      end do

   end subroutine get_nearest_neighbours_single_box
      

   !> return the box position by unsing periodic boundary conditions
   !> within the integer grid of the boxes
   function return_box_position_1d( x, n_max , add ) result( pos )

      integer, intent( in ) ::x
      integer, intent( in ) ::n_max
      integer, intent( in ) ::add

      integer ::pos

      pos = x + add
      if ( pos .lt. 1 ) then
         pos =  n_max + pos
      else if ( pos .gt. n_max ) then
         pos =  pos - n_max
      end if

   end function return_box_position_1d
     
   !> compute the direct lattice B for a given lattice A
   !> on output B contains the direct lattice, scaling
   !> positions into a range [0,1]
   !> and OMEGA contains the volume of lattice A
   subroutine direct_lattice ( a, b, omega )
      implicit none
 ! input variables
      real(q), intent(in) :: a(:,:) !(1:3,1:3)
 ! output variables
      real(q), intent(out) :: b(:,:) !(1:3,1:3)
      real(q), intent(out) :: omega
 ! calculate reciprocal lattice vectors and cell volume
      omega=(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
            +a(2,1)*(a(1,3)*a(3,2)-a(1,2)*a(3,3)) &
            +a(3,1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2)))
      b(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2)) / omega
      b(2,1)=(a(3,2)*a(1,3)-a(3,3)*a(1,2)) / omega
      b(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2)) / omega
      b(1,2)=(a(2,3)*a(3,1)-a(2,1)*a(3,3)) / omega
      b(2,2)=(a(3,3)*a(1,1)-a(3,1)*a(1,3)) / omega
      b(3,2)=(a(1,3)*a(2,1)-a(1,1)*a(2,3)) / omega
      b(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1)) / omega
      b(2,3)=(a(3,1)*a(1,2)-a(3,2)*a(1,1)) / omega
      b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1)) / omega
   end subroutine direct_lattice



   !> assigning atoms to boxes with a histogram method
   !> and using the build up principle of the boxes
   subroutine assign_atoms_to_box( this, positions, cpu_node, n_step )
      !use lattice_ml, only :DIRECT_LATTICE 
      class( neighbor_cell_list_type ), intent( inout ) :: this
      real( kind = q ),                 intent( in )    :: positions( :, : )
      integer, intent( in ) :: cpu_node
      integer, intent( in ) :: n_step

      real( kind = q ) :: unit_vectors( 1:3, 1:3 )

      integer :: i
      integer :: n_atom_tmp
      integer :: n_atom

      real( kind=q ) ::dir_x
      real( kind=q ) ::dir_y
      real( kind=q ) ::dir_z
      
      integer :: indx_x
      integer :: indx_y
      integer :: indx_z

      integer, allocatable :: counter( : )
      integer :: ierror
      integer :: box_indx
      integer :: nx
      integer :: ny
      integer :: nz
      integer :: n2
      real( kind = q ) :: omega


      PROFILING_START( "ASSIGN_ATOMS_TO_BOX" )
      
      call direct_lattice( this%unit_box( 1:, 1: ), &
                           unit_vectors( 1:, 1: ), omega )


      call this%allocate_my_boxes( size( positions, 2 ) )

      allocate( counter( 1:product( this%number_boxes ) ) )
      counter   = 0

      nx = this%number_boxes( 1 )
      ny = this%number_boxes( 2 )
      nz = this%number_boxes( 3 )

      n2 = ny * nz

      do n_atom = 1, size( positions,2 )
         dir_x  =  positions( 1, n_atom ) * unit_vectors( 1, 1 ) + &
                   positions( 2, n_atom ) * unit_vectors( 2, 1 ) + &
                   positions( 3, n_atom ) * unit_vectors( 3, 1 )
         indx_x =  ceiling( dir_x )

         dir_y  =  positions( 1, n_atom ) * unit_vectors( 1, 2 ) + &
                   positions( 2, n_atom ) * unit_vectors( 2, 2 ) + &
                   positions( 3, n_atom ) * unit_vectors( 3, 2 )
         indx_y =  ceiling( dir_y )

         dir_z  =  positions( 1, n_atom ) * unit_vectors( 1, 3 ) + &
                   positions( 2, n_atom ) * unit_vectors( 2, 3 ) + &
                   positions( 3, n_atom ) * unit_vectors( 3, 3 )
         indx_z =  ceiling( dir_z )

         indx_x = check_bin_bound( indx_x , nx )
         indx_y = check_bin_bound( indx_y , ny )
         indx_z = check_bin_bound( indx_z , nz )
         
         box_indx = (indx_x-1) * n2 + & 
                           (indx_y-1)*nz + indx_z
         this%atom_to_box( n_atom ) =  box_indx 

         counter( box_indx ) = counter( box_indx ) + 1

      end do
      
      do i = 1, size( this%nn_boxes, 1 )
         if ( allocated( this%nn_boxes( i )%box_to_atoms ) ) &
                      deallocate( this%nn_boxes( i )%box_to_atoms )
         allocate( this%nn_boxes( i )%box_to_atoms( 1:counter(i) ) )
         counter( i ) = 0
      end do

      do n_atom = 1, size( positions,2 )
         box_indx  =  this%atom_to_box( n_atom )
         counter( box_indx )   =  counter( box_indx ) + 1
         this%nn_boxes( box_indx )%box_to_atoms( counter( box_indx ) )  =  n_atom
      end do
      

      deallocate( counter )
      
      PROFILING_STOP( "ASSIGN_ATOMS_TO_BOX" )

   end subroutine assign_atoms_to_box


   

   !> take the new lattice and update the box paremeters
   !> if number of boxes changed in a direction
   !> recompute neighbor list between boxes
   !> assign atoms and sort them for faster access when
   !> computing neighbor list data
   subroutine update_box_list( this, positions, lattice, &
                                ncpu, node_me )

       class( neighbor_cell_list_type ), intent( inout ) :: this
       real( kind = q ),                 intent( in )    :: positions( : , : )
       real( kind = q ),                 intent( in )    :: lattice( : , : )
       integer, intent( in ), optional  :: ncpu
       integer, intent( in ), optional  :: node_me

       integer :: cpu_node
       integer :: n_step

       integer :: number_boxes_old( 1:3 )
       integer :: algo_old



       if ( present( node_me ) .and. present( ncpu ) ) then
          n_step   =  ncpu
          cpu_node =  node_me
       else  
          n_step   = 1
          cpu_node = 0
       end if


       number_boxes_old( 1 ) =  this%number_boxes( 1 )
       number_boxes_old( 2 ) =  this%number_boxes( 2 )
       number_boxes_old( 3 ) =  this%number_boxes( 3 )
       algo_old              =  this%algo_type



       PROFILING_START( "MAKE_BOX_PARAMS" )
       call this%create_box_parameters( lattice )
       PROFILING_STOP( "MAKE_BOX_PARAMS" )

       
       if ( this%algo_type .eq. 1 ) then
          PROFILING_START( "CREATE_BOX_NEIGH" )
          if ( number_boxes_old( 1 ) .ne. this%number_boxes( 1 ) .or. &
               number_boxes_old( 2 ) .ne. this%number_boxes( 2 ) .or. &
               number_boxes_old( 3 ) .ne. this%number_boxes( 3 ) .or. &
               ( algo_old .ne. this%algo_type .and. this%algo_type .ne. 2 ) ) then
             call this%create_box_neighbor_list_main
          end if
          PROFILING_STOP( "CREATE_BOX_NEIGH" )
          
          call this%assign_atoms_to_box( positions, cpu_node, n_step )

          PROFILING_START( "SORT_ATOMS" )
          call this%sort_atoms( positions )
          PROFILING_STOP( "SORT_ATOMS" )
          
       end if

   end subroutine update_box_list
      

      
   !> main routine to update the nearest neighbor list with
   !> the cell list algorithm
   !> first update cell_list data for current lattice
   !> recompute the nearest neighbor list with either the 
   !> nsquared algorithm or the cell list, depending what suits better
   subroutine update_nearest_neighbor_list_main_all( this, positions, types,&
                                                     lattice, inverse_lattice, &
                                                     node_me, ncpu )
     
      class( neighbor_cell_list_type ),intent( inout ) :: this
      real( kind = q ),intent( in )    :: positions( : , : )
      real( kind = q ),intent( in )    :: lattice( : , : )
      real( kind = q ),intent( in )    :: inverse_lattice( : , : )
      integer,intent( in )             :: types( : )
      integer, intent( in ), optional  :: node_me
      integer, intent( in ), optional  :: ncpu


      integer :: n_atom
      integer :: n_atom_tmp

      integer :: n_step
      integer :: cpu_node

      integer :: itype

      if ( present( node_me ) .and. present( ncpu ) ) then
         n_step   =  ncpu
         cpu_node =  node_me
      else
         n_step   = 1
         cpu_node = 0
      end if


      PROFILING_START( "UPDATE_BOX_LIST" )
      call this%update_box_list( positions, lattice, n_step, cpu_node )
      PROFILING_STOP( "UPDATE_BOX_LIST" )

      this%max_neighbors = 1
      ! select different algorithms
      this%max_neighbor = 0
      select case ( this%algo_type )
         case( NLA_CELL )
            call this%sort_atoms( positions )
            do n_atom_tmp = 1, size( positions, 2 ), n_step
               n_atom  =  n_atom_tmp + cpu_node
               if ( n_atom .le. size( positions , 2 ) ) then

                  call this%compute_nn_list_single_atom( positions, lattice,&
                                                         inverse_lattice,&
                                                         types, n_atom )
                  do itype = 1, this%number_types
                     this%max_neighbor( itype ) = max( this%max_neighbor( itype ), &
                           this%nn_atoms( n_atom )%max_neighbor( itype ) )
                  end do
               end if
            end do
         case ( NLA_SQUARED )
            do n_atom_tmp = 1, size( positions, 2 ), n_step
               n_atom  =  n_atom_tmp + cpu_node
               if ( n_atom .le. size( positions , 2 ) ) then
                  call this%nn_atoms( n_atom )%compute_nn_list( positions, lattice,&
                                                types, inverse_lattice )
                  do itype = 1, this%number_types     
                     this%max_neighbor( itype ) = max( this%max_neighbor( itype ), & 
                     this%nn_atoms( n_atom )%max_neighbor( itype ) )           
                  end do
               end if
            end do
      end select

      this%max_neighbors = maxval( this%max_neighbor )

   end subroutine update_nearest_neighbor_list_main_all


   !> return the maximum number of neighbor atoms
   !> the max is taken between different atom types
   function get_max_number_neighbors( this ) result( nn_neighbors )

      class( neighbor_cell_list_type ),intent( in ) :: this
      integer :: nn_neighbors
      
      nn_neighbors = this%max_neighbors

   end function get_max_number_neighbors


   !> create a nearest neighbor list from a super list with a larger
   !> cutoff radius. Routine assumes that the larger nearest neighbor
   !> list was already computed and copies data from
   !> the super list to construct the sublist
   subroutine create_subset_cutoff_list( this, other, node_me, ncpu )

      class( neighbor_cell_list_type ),intent( inout ) :: this
      class( neighbor_cell_list_type ),intent( in ) :: other
      integer, intent( in ), optional  :: node_me
      integer, intent( in ), optional  :: ncpu

      integer :: n_step
      integer :: cpu_node
      integer :: n_atom_tmp
      integer :: n_atom

      integer :: itype

      if ( size( this%nn_atoms ) .ne. size( other%nn_atoms ) ) then
         write( *, * ) "Error in cell list list algorithm"
         write( *, * ) "atom numbers do not agree"
         return
      end if
      
      if ( present( node_me ) .and. present( ncpu ) ) then
         n_step   =  ncpu
         cpu_node =  node_me
      else
         n_step   = 1
         cpu_node = 0
      end if

      this%max_neighbors = 1
      this%max_neighbor  = 0
      do n_atom_tmp = 1, size( this%nn_atoms ), n_step 
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( this%nn_atoms ) ) then
            call this%nn_atoms( n_atom )%create_subset_cutoff_list( other%nn_atoms( n_atom ) )
            do itype = 1, this%number_types
               this%max_neighbor( itype ) = max( this%max_neighbor( itype ), &
                           this%nn_atoms( n_atom )%max_neighbor( itype ) )
            end do
         end if
      end do

      this%max_neighbors = maxval( this%max_neighbor )

   end subroutine create_subset_cutoff_list



   !> compute the nearest neighbor list for a single atom
   !> with the cell list algorithm. Shows a linear scaling
   subroutine create_nearest_neighbor_list_single_atom_cell_list_sorted( this,&
                                                    positions, lattice, inverse_lattice, types, indx )

      use neighbor_list_mod, only : set_periodic_size_n_square
      class( neighbor_cell_list_type ), intent( inout ) :: this
      real( kind = q ),                 intent( in ) :: positions( :, : )
      real( kind = q ),                 intent( in ) :: lattice( :, : )
      real( kind = q ),                 intent( in ) :: inverse_lattice( :, : )
      integer,                          intent( in ) :: types( : )
      integer,                          intent( in ) :: indx

      integer :: box_indx
      integer :: n_box_indx
      integer :: ibox
      integer :: j
      integer :: n_atom

      real( kind = q )  :: atom( 1:3 )
      real( kind = q )  :: neighbor( 1:3 )
      real( kind = q )  :: cutoff_squared

      integer ::x
      integer ::y
      integer ::z
      integer ::n_counter

      integer :: xyz_max( 1:3 )
      integer :: xyz_min( 1:3 )

      
      ! temprorary lists
      integer,allocatable ::n_list( : )   ! neighbor index
      integer,allocatable ::n_types( : )  ! neighbor atom type
      real( kind = q ),allocatable ::n_dist( : )       ! neighbor distance
      real( kind = q ),allocatable ::n_vectors( :, : ) ! neighbor connection unit vector

      real( kind=q ) :: eps
      real( kind=q ) :: dx
      real( kind=q ) :: dy
      real( kind=q ) :: dz
      real( kind=q ) :: dist




      eps = 1.0e-6_q


      call allocate_working_arrays_for_data_computation( n_list, n_types, &
                       n_dist, n_vectors,size( positions,2 ) )

      box_indx  =  this%atom_to_box( indx )
      atom( 1 ) = positions( 1 , indx )
      atom( 2 ) = positions( 2 , indx )
      atom( 3 ) = positions( 3 , indx )

      call set_periodic_size_n_square( this%cutoff, atom, inverse_lattice ,&
                                       xyz_max, xyz_min )

      cutoff_squared = this%nn_atoms( indx )%return_cutoff_squared()
      n_counter = 0
      !loop over nearest neighbor boxes
      this%max_neighbor = 0
      do x = xyz_min(1), xyz_max(1)
         do y = xyz_min(2), xyz_max(2)
            do z = xyz_min(3), xyz_max(3)
               do ibox = 1 , size( this%nn_boxes( box_indx )%nn_box_list, 1 )
                  ! extract neighbor box index
                  n_box_indx =  this%nn_boxes( box_indx )%nn_box_list( ibox )
                  ! loop over atoms in box
                  do j =1 , size( this%nn_boxes( n_box_indx )%box_to_atoms_sorted, 1 )
                     n_atom = this%nn_boxes( n_box_indx )%box_to_atoms_sorted( j )

                  neighbor(1)  =  this%positions_sorted( 1, n_atom )
                  neighbor(2)  =  this%positions_sorted( 2, n_atom )
                  neighbor(3)  =  this%positions_sorted( 3, n_atom )
                  neighbor(1)  =  neighbor(1) + real( x, kind=q ) * lattice( 1, 1 ) + &
                                                real( y, kind=q ) * lattice( 1, 2 ) + &
                                                real( z, kind=q ) * lattice( 1, 3 )
                  neighbor(2)  =  neighbor(2) + real( x, kind=q ) * lattice( 2, 1 ) + &
                                                real( y, kind=q ) * lattice( 2, 2 ) + &
                                                real( z, kind=q ) * lattice( 2, 3 )
                  neighbor(3)  =  neighbor(3) + real( x, kind=q ) * lattice( 3, 1 ) + &
                                                real( y, kind=q ) * lattice( 3, 2 ) + &
                                                real( z, kind=q ) * lattice( 3, 3 )
                  dx  = atom(1) - neighbor(1)
                  dy  = atom(2) - neighbor(2)
                  dz  = atom(3) - neighbor(3)
                  dist = dx*dx + dy*dy + dz*dz
                  if ( dist .le. cutoff_squared .and. dist .gt. eps ) then
                     n_counter = n_counter + 1
                     n_list( n_counter )  =  this%sorted_indexes( n_atom )
                     n_types( n_counter ) =  types( this%sorted_indexes( n_atom ) )
                     n_dist( n_counter )  =  dsqrt( dist )
                     this%max_neighbor( n_types( n_counter ) ) = &
                          this%max_neighbor( n_types( n_counter ) ) + 1
                     n_vectors( 1 ,n_counter ) = dx / n_dist( n_counter )
                     n_vectors( 2 ,n_counter ) = dy / n_dist( n_counter )
                     n_vectors( 3 ,n_counter ) = dz / n_dist( n_counter )
                  end if
               end do
               end do
            end do
         end do
      end do

      call this%nn_atoms( indx )%update_data_arrays( n_list, n_types, &
                                          n_dist, n_vectors, n_counter )
      
   end subroutine create_nearest_neighbor_list_single_atom_cell_list_sorted
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! helper routines 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      

   !> sorting atoms for faster memory access when computing
   !> nearest neighbors. Only gains speed up for
   !> large cutoff radius
   subroutine sort_atoms( this, positions )

      class( neighbor_cell_list_type ), intent( inout ) :: this
      real( kind = q ), intent( in ) :: positions( : , : )

      integer  :: n_box
      integer  :: n
      integer  :: n_atom
      integer  :: col


      call this%allocate_arrays_sorting( size( positions, 2 ) )



      col = 1
      do n_box = 1, size( this%nn_boxes, 1 )
         if ( allocated( this%nn_boxes( n_box )%box_to_atoms_sorted ) ) &
            deallocate( this%nn_boxes( n_box )%box_to_atoms_sorted )
            allocate( this%nn_boxes( n_box )%box_to_atoms_sorted(&
                        1:size( this%nn_boxes( n_box )%box_to_atoms, 1 ) ) )
         do n = 1, size( this%nn_boxes( n_box )%box_to_atoms, 1 )
            n_atom = this%nn_boxes( n_box )%box_to_atoms( n )
            this%positions_sorted( : , col ) = positions( : , n_atom )
            this%nn_boxes( n_box )%box_to_atoms_sorted( n ) = col
            this%sorted_indexes( col ) = n_atom
            col = col + 1
         end do
      end do

   end subroutine sort_atoms



   !> checking the boundary conditions when assigning the
   !> atoms to boxes by binning
   function check_bin_bound( indx, maxi ) result( index_new )
      
      integer, intent( in )  ::indx
      integer, intent( in )  ::maxi
      
      integer :: index_new

      index_new = indx
      index_new = min( index_new, maxi )
      index_new = max( index_new, 1 )

   end function check_bin_bound
     


   !> allocating box connection list
   subroutine allocate_my_boxes( this, n )

      class( neighbor_cell_list_type ), intent( inout ) :: this
      integer,                          intent( in )    :: n

      if ( allocated( this%atom_to_box ) ) then
         if ( size( this%atom_to_box , 1 ) .ne. N ) then
            deallocate( this%atom_to_box )
            allocate( this%atom_to_box( 1:n ) )
         end if
      else
            allocate( this%atom_to_box( 1:n ) )
      end if

   end subroutine allocate_my_boxes
   


   !> allocating neighbor lists betwen boxes
   subroutine allocate_box_neighbor_list( this, nsize )

      
      class( neighbor_cell_list_type ), intent( inout ) :: this
      integer,   intent( in ) :: nsize
      
      integer :: x_box


      if ( allocated( this%nn_boxes ) ) then
         if ( size( this%nn_boxes, 1 ) .ne. product( this%number_boxes ) ) then
            deallocate( this%nn_boxes )
            allocate( this%nn_boxes( 1:product( this%number_boxes ) ) )
            do x_box = 1, size( this%nn_boxes )
               allocate( this%nn_boxes( x_box )%nn_box_list( 1:nsize ) )
               this%nn_boxes( x_box )%nn_box_list( : ) = 0
            end do
         end if
      else
         allocate( this%nn_boxes( 1:product( this%number_boxes ) ) )
         do x_box = 1, size( this%nn_boxes )
            allocate( this%nn_boxes( x_box )%nn_box_list( 1:nsize ) )
            this%nn_boxes( x_box )%nn_box_list( : ) = 0
         end do
      end if

   end subroutine allocate_box_neighbor_list

      
   !> allocate arrays for the sorted positions
   !> positions are sorted according to their assigned boxes
   subroutine allocate_arrays_sorting( this, n )

      class( neighbor_cell_list_type ), intent( inout ) :: this
      integer,                          intent( in )    :: n

      if ( allocated( this%positions_sorted ) ) then
         if ( size( this%positions_sorted, 2 ) .ne. n ) then
            deallocate( this%positions_sorted )
            allocate( this%positions_sorted( 1:3 , 1:n ) )
            deallocate( this%sorted_indexes )
            allocate( this%sorted_indexes( 1:n ) )
         end if
      else
         allocate( this%positions_sorted( 1:3 , 1:n ) )
         allocate( this%sorted_indexes( 1:n ) )
      end if

   end subroutine allocate_arrays_sorting


   !> allocate temporaray data arrays for nearest neighbor list data
   !> computation
   subroutine allocate_working_arrays_for_data_computation( n_list, n_types, &
                                                            n_dist, n_vectors, n_size )

      integer,allocatable,intent( inout ) :: n_list( : )   ! neighbor index
      integer,intent(inout),allocatable   :: n_types( : ) ! neighbor atom type
      real( kind = q ),intent( inout ),allocatable :: n_dist( : )       ! neighbor distance
      real( kind = q ),intent( inout ),allocatable :: n_vectors( :, : ) ! neighbor connection unit vector
      integer,intent( in ) ::n_size
      
      allocate( n_list( 1:n_size ) )
      allocate( n_types( 1:n_size ) )
      allocate( n_dist( 1:n_size ) )
      allocate( n_vectors( 1:3, 1:n_size ) )

   end subroutine allocate_working_arrays_for_data_computation
         

   !> copy the data from the cell list algorithm from temprorary
   !> arrays for computing to the actual member arrays of the
   !> neighbor_cell_list_type type
   subroutine update_data_arrays( this, n_list, n_types, n_dist, n_vectors, n_counter, n_atom )
      
      class( neighbor_cell_list_type ),intent( inout ) :: this
      
      integer,allocatable,intent( inout ) ::n_list( : )   ! neighbor index
      integer,intent(inout),allocatable ::n_types( : ) ! neighbor atom type
      real( kind = q ),intent( inout ),allocatable ::n_dist( : )       ! neighbor distance
      real( kind = q ),intent( inout ),allocatable ::n_vectors( :, : ) ! neighbor connection unit vector
      integer,intent( in ) ::n_counter
      integer,intent( in ) ::n_atom



      ! organizing array sizes
      if ( n_counter .ne. size( this%nn_atoms( n_atom )%neighbors ) ) then
         deallocate( this%nn_atoms( n_atom )%neighbors )
         allocate( this%nn_atoms( n_atom )%neighbors( 1:n_counter ) )
         ! if size does not fit for one array it won't fit for any
         deallocate( this%nn_atoms( n_atom )%neighbors_types )
         allocate( this%nn_atoms( n_atom )%neighbors_types( 1:n_counter ) )
         deallocate( this%nn_atoms( n_atom )%neighbor_dist )
         allocate( this%nn_atoms( n_atom )%neighbor_dist( 1:n_counter ) )
         deallocate( this%nn_atoms( n_atom )%neighbor_connect )
         allocate( this%nn_atoms( n_atom )%neighbor_connect( 1:3, 1:n_counter ) )
      end if
      ! copy data to public arrays for usage
      this%nn_atoms( n_atom )%neighbors        = n_list( 1:n_counter )
      this%nn_atoms( n_atom )%neighbors_types  = n_types( 1:n_counter )
      this%nn_atoms( n_atom )%neighbor_dist    = n_dist( 1:n_counter )
      this%nn_atoms( n_atom )%neighbor_connect = n_vectors( 1:3, 1:n_counter )

      

      ! deallocate working arrays
      deallocate( n_list )
      deallocate( n_types )
      deallocate( n_dist )
      deallocate( n_vectors )

   end subroutine update_data_arrays

end module neighbor_cell_list_mod





module verlet_cell_list_mod
   

   use ML_FF_PREC
   use MPI_HELP
   use neighbor_cell_list_mod, only : neighbor_cell_list_type
   use neighbor_list_mod, only : neighbor_list_type
   implicit none
   public :: verlet_cell_list_type

   ! Neighbor list algorithm
   enum, bind(c)
      enumerator :: NLA_ENUM     = 0
      enumerator :: NLA_VERLET   = 1 
      enumerator :: NLA_SQUARED  = 2
   end enum

   
   private 
   type verlet_cell_list_type
      private
      !> enum to determine algorithm type
      integer(kind(NLA_ENUM)) :: algo_type
      !> cell list type to update nearest neighbor list
      type( neighbor_cell_list_type ), public :: cell_list
      !> cutoff radius
      real( kind = q ) :: cutoff
      !> cutoff radius squared
      real( kind = q ) :: cutoff_squared
      !> verlet cutoff radius 1.5 cutoff
      real( kind = q ) :: verlet_cutoff
      !> verlet travel distance of particle after which
      !> update has to be done
      real( kind = q ) :: verlet_distance
      real( kind = q ),allocatable :: positions_old( :, : )
      real( kind = q ),allocatable :: travel_dist( :, : )
      real( kind = q ),allocatable :: shift_xyz( :, : )
      type( neighbor_list_type ), public, allocatable :: nn_list( : )
      !> logical if main routine has to be recomputed
      logical  :: recompute
      !> maximum number of neighbors
      integer  :: max_neighbors
      !> number of different atom types
      integer :: number_types
      !> maximum number of neighbors per atom type
      integer, public, allocatable :: max_neighbor( : )

      contains
      private
      procedure, public, pass( this ) :: init_nn_cell_list => verlet_list_algo_init
      generic, public :: compute_nn_list => verlet_list_algo_main,&
                                            create_verlet_subset_cutoff_list
      procedure, public, pass( this ) :: forced_update => verlet_list_algo_main_forced_update
      procedure, public, pass( this ) :: get_max_number_neighbors

      procedure, pass( this ) :: verlet_list_algo_main
      procedure, pass( this ) :: create_verlet_subset_cutoff_list
      procedure, pass( this ) :: compute_verlet_distances
      procedure, pass( this ) :: update_verlet_distances
      procedure, pass( this ) :: update_nn_list_after_cell_recomputation
      procedure, pass( this ) :: recompute_neighbor_data
      procedure, pass( this ) :: update_single_atom
      procedure, pass( this ) :: make_square_algo_update
      procedure, pass( this ) :: rezero_verlet_distances
      procedure, pass( this ) :: verlet_algo_type_test

   end type verlet_cell_list_type


   contains


   !> initlialize the verlet cell list algorithm
   subroutine verlet_list_algo_init( this, cutoff, lattice, positions, ntypes )


      class ( verlet_cell_list_type ),intent( inout ) :: this
      real( kind = q ),intent( in )  :: cutoff
      real( kind = q ), intent( in ) :: lattice( :, : )
      real( kind = q ), intent( in ) :: positions( :, : )
      integer, intent( in ) :: ntypes

      integer :: i_atom


      this%cutoff          =  cutoff
      this%cutoff_squared  =  cutoff * cutoff
      this%verlet_cutoff   =  1.5_q * cutoff
      this%verlet_distance =  this%verlet_cutoff - this%cutoff
      this%verlet_distance =  this%verlet_distance * this%verlet_distance
      this%number_types    =  ntypes
      

      call this%cell_list%init_nn_cell_list( this%verlet_cutoff, &
                                             lattice, positions, ntypes )

      allocate( this%positions_old( 1:3 , 1:size( positions,2 ) ) )
      this%positions_old =  positions
      allocate( this%travel_dist( 1:3, 1:size( positions,2 ) ) )
      this%travel_dist  =  0.0_q
      allocate( this%shift_xyz( 1:3, 1:size( positions,2 ) ) )
      this%recompute     =  .true.
      allocate( this%nn_list( 1:size( positions, 2 ) ) )

      do i_atom = 1, size( positions,2 )
         call this%nn_list( i_atom )%init_neighbor_list( this%cutoff, i_atom, ntypes )
      end do
      ! ~~~~~~~~~~~~~~~~~~~~~

      allocate( this%max_neighbor( 1:ntypes ) )
      this%max_neighbor = 0


   end subroutine verlet_list_algo_init



   !> force a complete update of the nearest neighbor list with the 
   !> verlet algorithm even when the distance criterion is not satiesfied
   !> routine also resets the travel_dist array to zero
   subroutine verlet_list_algo_main_forced_update( this, positions, types, lattice, &
                                                   inverse_lattice, node_me, ncpu )

      class( verlet_cell_list_type ),intent( inout ) :: this
      real( kind = q ),intent( in )    :: positions( : , : )
      real( kind = q ),intent( in )    :: lattice( : , : )
      real( kind = q ),intent( in )    :: inverse_lattice( : , : )
      integer,intent( in )             :: types( : )
      
      integer, intent( in ), optional  :: node_me
      integer, intent( in ), optional  :: ncpu

      integer :: n_step
      integer :: cpu_node


      if ( present( node_me ) .and. present( ncpu ) ) then
         n_step   =  ncpu
         cpu_node =  node_me
      else
         n_step   =  1
         cpu_node =  0
      end if
         
      call this%cell_list%compute_nn_list( positions, types, lattice,&
                                           inverse_lattice, cpu_node, n_step )
      call this%update_nn_list_after_cell_recomputation( cpu_node, n_step )
      this%travel_dist  =  0.0_q

   end subroutine verlet_list_algo_main_forced_update


   !> check the size of the simulation box and compare to
   !> the verlet cutoff radius, if the application
   !> of the verlet algorithm is suitable.
   !> otherwise the routine will decide to use 
   !> the Nsquare algorithm as fallback
   !>
   subroutine verlet_algo_type_test( this, lattice )   
      
      class( verlet_cell_list_type ),intent( inout ) :: this
      real( kind = q ), intent( in ) :: lattice( : , : )

      real( kind = q ) :: norm_square
      real( kind = q ) :: cutoff_squared
      integer :: i

      if ( this%algo_type .eq. 2 ) then

      end if
      this%algo_type = 1

      cutoff_squared = 4.0_q * this%verlet_cutoff * this%verlet_cutoff 

      norm_square = lattice( 1, 1 ) * lattice( 1, 1 ) + &
                    lattice( 2, 1 ) * lattice( 2, 1 ) + &
                    lattice( 3, 1 ) * lattice( 3, 1 )
      if ( norm_square .lt. cutoff_squared ) then
         this%algo_type = 2
         return
      end if

      norm_square = lattice( 1, 2 ) * lattice( 1, 2 ) + &
                    lattice( 2, 2 ) * lattice( 2, 2 ) + &
                    lattice( 3, 2 ) * lattice( 3, 2 )
      if ( norm_square .lt. cutoff_squared ) then
         this%algo_type = 2
         return
      end if

      norm_square = lattice( 1, 3 ) * lattice( 1, 3 ) + &
                    lattice( 2, 3 ) * lattice( 2, 3 ) + &
                    lattice( 3, 3 ) * lattice( 3, 3 )
      if ( norm_square .lt. cutoff_squared ) then
         this%algo_type = 2
         return
      end if

   end subroutine verlet_algo_type_test



   !> main routine of the verlet algorithm. This function
   !> has to be called by the users. Checks if update of
   !> the verlet list is necessary. Otherwise update the
   !> nearest neighbor list from the verlet list
   subroutine verlet_list_algo_main( this, positions, types, lattice, &
                                     inverse_lattice, node_me, ncpu, mpi_communicator )

      class( verlet_cell_list_type ),intent( inout ) :: this
      real( kind = q ), intent( in )    :: positions( :, : )
      integer, intent( in )             :: types( : )
      real( kind = q ), intent( in )    :: lattice( :, : )
      real( kind = q ), intent( in )    :: inverse_lattice( :, : ) 

      integer, intent( in ), optional  :: node_me
      integer, intent( in ), optional  :: ncpu
      integer, intent( in ), optional  :: mpi_communicator

      integer :: n_step
      integer :: cpu_node
      integer :: ierror
      integer :: algo_type_old


      algo_type_old  =  this%algo_type

      if ( present( node_me ) .and. present( ncpu ) ) then
         n_step   =  ncpu
         cpu_node =  node_me
      else
         n_step   = 1
         cpu_node = 0
      end if

      !! determine if verlet nearest-neighbor algorithm can be used
      call this%verlet_algo_type_test( lattice )

      select case ( this%algo_type )
         case( NLA_VERLET )
            if ( algo_type_old .eq. this%algo_type ) then
               this%recompute = .false.
               call this%compute_verlet_distances( positions, lattice, cpu_node, n_step )
               if ( present( mpi_communicator ) .and. ncpu .gt. 1 ) then
                  call MPI_allreduce( MPI_IN_PLACE, this%recompute, 1, MPI_LOGICAL, MPI_LOR, mpi_communicator,ierror )
               end if
            else
               ! force recompute because in last step the nsquared algo was used 
               this%recompute = .True.
            end if

            if ( this%recompute ) then
               call this%cell_list%compute_nn_list( positions, types, lattice,&
                                                    inverse_lattice, cpu_node, n_step )
               call this%update_nn_list_after_cell_recomputation( cpu_node, n_step )
               call this%rezero_verlet_distances( cpu_node, n_step )
            else 
               call this%recompute_neighbor_data( positions, lattice, inverse_lattice,&
                                                  types, cpu_node, n_step )
            end if
         case ( NLA_SQUARED )
           call this%make_square_algo_update( positions, types, lattice, &
                                              inverse_lattice, cpu_node, n_step )
      end select

      this%positions_old  =  positions

   end subroutine verlet_list_algo_main


   !> if the verlet nearest neighbor algorithm is not suitable
   !> to update the nearest neighbor list. This routine
   !> implements the nsquare algorithm as fallback
   subroutine make_square_algo_update( this, positions, types, lattice, &
                                        inverse_lattice, cpu_node, n_step )

      class( verlet_cell_list_type ),intent( inout ) :: this
      real( kind = q ), intent( in )    :: positions( :, : )
      integer, intent( in )             :: types( : )
      real( kind = q ), intent( in )    :: lattice( :, : )
      real( kind = q ), intent( in )    :: inverse_lattice( :, : ) 

      integer :: n_step
      integer :: cpu_node

      integer :: n_atom
      integer :: n_atom_tmp
      integer :: itype
            
            
      this%max_neighbors =  1
      this%max_neighbor  =  0

      do n_atom_tmp = 1, size( this%nn_list ), n_step
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( this%nn_list ) ) then
            call this%nn_list( n_atom )%compute_nn_list( &
                               positions, lattice, types, inverse_lattice )
            do itype = 1, this%number_types
               this%max_neighbor( itype ) = max( &
                            this%nn_list( n_atom )%max_neighbor( itype ),&
                            this%max_neighbor( itype ) )
            end do
         end if
      end do
      this%max_neighbors  =  maxval( this%max_neighbor )


   end subroutine make_square_algo_update


   !> if there are two verlet lists used in the main
   !> code with differnt cutoff radii. This routine 
   !> can be usde to create a subset nearest neighbor
   !> list from the the list with the larger cutoff.
   !> routine assumes that the superlist was already 
   !> computed
   subroutine create_verlet_subset_cutoff_list( this, other, node_me, ncpu )

      class( verlet_cell_list_type ),intent( inout ) :: this
      class( verlet_cell_list_type ),intent( in ) :: other
      integer, intent( in ),optional :: node_me
      integer, intent( in ),optional :: ncpu

      integer :: cpu_node
      integer :: n_step

      integer :: n_atom
      integer :: n_atom_tmp
      integer :: itype

   

      if ( present( node_me ) .and. present( ncpu ) ) then
         n_step   =  ncpu
         cpu_node =  node_me
      else
         n_step   = 1
         cpu_node = 0
      end if


      this%max_neighbors =  1
      this%max_neighbor  =  0
      do n_atom_tmp = 1, size( this%nn_list ), n_step
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( this%nn_list ) ) then
            call this%nn_list( n_atom )%create_subset_cutoff_list( &
                                                   other%nn_list( n_atom ) )
            do itype = 1, this%number_types
                this%max_neighbor( itype ) = max( &
                              this%nn_list( n_atom )%max_neighbor( itype ),&
                              this%max_neighbor( itype ) )
            end do
         end if
      end do
      this%max_neighbors = maxval( this%max_neighbor )


   end subroutine create_verlet_subset_cutoff_list



   !> compute the shortest traveling distances of the atoms under periodic
   !> boundary conditions. Copy the computed data from temporary arrys to
   !> the class member array
   subroutine compute_verlet_distances( this, positions, lattice, cpu_node, n_step )
      

      class( verlet_cell_list_type ),intent( inout ) :: this
      real( kind = q ),intent( in )    :: positions( : , : )
      real( kind = q ),intent( in )    :: lattice( : , : )
      integer, intent( in ) :: cpu_node
      integer, intent( in ) :: n_step

      integer :: n_atom
      integer :: n_atom_tmp

      real( kind = q ),allocatable :: temp_dist( : , : )
      real( kind = q ),allocatable :: dist_square( : )
      real( kind = q ) :: max_distance
      real( kind = q ) :: delta( 1:3 )
      real( kind = q ) :: pos( 1:3 )
      real( kind = q ) :: pos_old( 1:3 )
      real( kind = q ) :: norm

      integer :: xx
      integer :: yy
      integer :: zz
      
      if ( size( this%travel_dist, 2 ) .ne. size( positions, 2 ) ) then
         deallocate( this%travel_dist )
         allocate( this%travel_dist( 1:3, 1:size( positions, 2 ) ) )
         this%travel_dist =  0.0_q
         this%recompute    =  .true.
         deallocate( this%nn_list )
         allocate( this%nn_list( 1:size( positions, 2 ) ) )
         do n_atom = 1, size( positions, 2 )
            call this%nn_list( n_atom )%init_neighbor_list( this%cutoff, n_atom, &
                                                            this%number_types )
         end do
         return
      end if


      allocate( temp_dist( 1:3 , 1:size( positions, 2 ) ) )
      allocate( dist_square( 1:size( positions, 2 ) ) )
      temp_dist    =  huge( temp_dist( 1, 1 ) )
      dist_square  =  huge( dist_square( 1 ) )
      do xx = -1,1
         do yy = -1,1
            do zz = -1,1
               do n_atom_tmp = 1, size( positions, 2 ), n_step
                  n_atom  =  n_atom_tmp + cpu_node
                  if ( n_atom .le. size( positions, 2 ) ) then
                     pos_old( 1 )  =  this%positions_old( 1, n_atom )
                     pos_old( 2 )  =  this%positions_old( 2, n_atom )
                     pos_old( 3 )  =  this%positions_old( 3, n_atom )
                     pos( 1 )      =  positions( 1, n_atom )
                     pos( 2 )      =  positions( 2, n_atom )
                     pos( 3 )      =  positions( 3, n_atom )
                     pos( 1 ) =  pos( 1 ) + real( xx, kind=q ) * lattice( 1, 1 ) + &
                                     real( yy, kind=q ) * lattice( 1, 2 ) + &
                                     real( zz, kind=q ) * lattice( 1, 3 )
                     pos( 2 ) =  pos( 2 ) + real( xx, kind=q ) * lattice( 2, 1 ) + &
                                     real( yy, kind=q ) * lattice( 2, 2 ) + &
                                     real( zz, kind=q ) * lattice( 2, 3 )
                     pos( 3 ) =  pos( 3 ) + real( xx, kind=q ) * lattice( 3, 1 ) + &
                                     real( yy, kind=q ) * lattice( 3, 2 ) + &
                                     real( zz, kind=q ) * lattice( 3, 3 )
                     delta( 1 ) =  pos( 1 ) - pos_old( 1 )
                     delta( 2 ) =  pos( 2 ) - pos_old( 2 )
                     delta( 3 ) =  pos( 3 ) - pos_old( 3 )
                     norm       =  delta( 1 ) * delta( 1 ) + &
                                   delta( 2 ) * delta( 2 ) + &
                                   delta( 3 ) * delta( 3 )
                     if ( norm .lt. dist_square( n_atom ) ) then
                        temp_dist( 1, n_atom )  =  delta( 1 )
                        temp_dist( 2, n_atom )  =  delta( 2 )
                        temp_dist( 3, n_atom )  =  delta( 3 )
                        dist_square( n_atom )   =  norm
                     end if
                  end if
               end do
            end do
         end do
      end do
      call this%update_verlet_distances( temp_dist, cpu_node, n_step )

      deallocate( temp_dist )
      deallocate( dist_square )


   end subroutine compute_verlet_distances


   !> copy shortest travelling distances from atoms during an MD
   !> step to the class member array this%travel_dist
   !> and check if the verlet liust has to be updated
   subroutine update_verlet_distances( this, temp_dist, cpu_node, n_step )

      class( verlet_cell_list_type ), intent( inout ) :: this
      real( kind = q ), intent( in ) :: temp_dist( :, : )
      integer, intent( in ) :: cpu_node
      integer, intent( in ) :: n_step

      integer :: n_atom
      integer :: n_atom_tmp

      real( kind = q ) :: square_dist( 1:3 )


      
      do n_atom_tmp = 1, size( temp_dist, 2 ), n_step
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( temp_dist, 2 ) ) then
            this%travel_dist( 1 , n_atom ) = this%travel_dist( 1, n_atom ) + temp_dist( 1, n_atom )
            this%travel_dist( 2 , n_atom ) = this%travel_dist( 2, n_atom ) + temp_dist( 2, n_atom )
            this%travel_dist( 3 , n_atom ) = this%travel_dist( 3, n_atom ) + temp_dist( 3, n_atom )
            square_dist( 1 )  = this%travel_dist( 1 , n_atom ) * this%travel_dist( 1 , n_atom )
            square_dist( 2 )  = this%travel_dist( 1 , n_atom ) * this%travel_dist( 1 , n_atom )
            square_dist( 3 )  = this%travel_dist( 1 , n_atom ) * this%travel_dist( 1 , n_atom )
            !if ( square_dist( 1 ) .gt. this%verlet_distance .or. & 
            !     square_dist( 2 ) .gt. this%verlet_distance .or. &
            !     square_dist( 3 ) .gt. this%verlet_distance ) then
            if ( square_dist( 1 ) + square_dist( 2 ) + square_dist( 3 ) .gt. this%verlet_distance ) then
                 this%recompute  =  .true.
                 exit
            end if
         end if
      end do

   end subroutine update_verlet_distances


   subroutine rezero_verlet_distances( this, cpu_node, n_step )

      class( verlet_cell_list_type ), intent( inout ) :: this
      integer, intent( in ) :: cpu_node
      integer, intent( in ) :: n_step
      
      integer :: n_atom
      integer :: n_atom_tmp
      
      do n_atom_tmp = 1, size( this%travel_dist, 2 ), n_step
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( this%travel_dist, 2 ) ) then
            this%travel_dist( 1 , n_atom ) = 0.0_q
            this%travel_dist( 2 , n_atom ) = 0.0_q
            this%travel_dist( 3 , n_atom ) = 0.0_q
         end if
      end do


   end subroutine rezero_verlet_distances



   !> copy shortest traveling distances from atoms during an MD  
   !> step to the class member array this%travel_dist            
   !> and check if the verlet liust has to be updated
   subroutine update_nn_list_after_cell_recomputation( this, cpu_node, n_step )
      
      class( verlet_cell_list_type ), intent( inout ) :: this
      integer, intent( in ) :: cpu_node
      integer, intent( in ) :: n_step

      integer :: n_atom
      integer :: n_atom_tmp
      integer :: itype


      this%max_neighbors = 1
      this%max_neighbor  = 0
      do n_atom_tmp = 1, size( this%nn_list ), n_step
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( this%nn_list ) ) then
            call this%nn_list( n_atom )%create_subset_cutoff_list(  &
                                   this%cell_list%nn_atoms( n_atom ) )
            do itype = 1, this%number_types
               this%max_neighbor( itype ) = max( this%max_neighbor( itype ), &
                      this%nn_list( n_atom )%max_neighbor( itype ) )
            end do
         end if
      end do

      this%max_neighbors = maxval( this%max_neighbor )

   end subroutine update_nn_list_after_cell_recomputation



   !> get maximum number of neighbors, whereby
   !> a max over different neighbor types was taken
   function get_max_number_neighbors( this ) result( nn )
      
      class( verlet_cell_list_type ), intent( in ) :: this
      integer :: nn

      nn  =  this%max_neighbors

   end function get_max_number_neighbors



   !> recompute the nearest neighbor list data with the use of
   !> the verlet list, in case the atoms in the verlet list were not updated
   subroutine recompute_neighbor_data( this, positions, lattice, inverse_lattice,&
                                             types, cpu_node, n_step )

      class( verlet_cell_list_type ), intent( inout ) :: this
      real( kind = q ), intent( in ) :: positions( :, : )
      real( kind = q ), intent( in ) :: lattice( :, : )
      real( kind = q ), intent( in ) :: inverse_lattice( :, : )
      integer, intent( in ) :: types( : )
      integer, intent( in ) :: cpu_node
      integer, intent( in ) :: n_step

      integer :: n_atom_tmp
      integer :: n_atom
      integer :: itype
      integer :: type_counter( 1:this%number_types )

      this%max_neighbors  =  1
      this%max_neighbor   =  0
      do n_atom_tmp = 1, size( positions, 2 ), n_step
         n_atom  =  n_atom_tmp + cpu_node
         if ( n_atom .le. size( positions, 2 ) ) then
            call this%update_single_atom( n_atom, positions, lattice, &
                                          inverse_lattice, types, type_counter )
            do itype = 1, this%number_types                                      
               this%max_neighbor( itype ) = max( this%max_neighbor( itype ), &   
                                                 type_counter( itype ) )
            end do
         end if
      end do
      this%max_neighbors  =  maxval( this%max_neighbor )

   end subroutine recompute_neighbor_data 



   !> helper routine for recompute_neighbor_data. Updating the nearest
   !> neighbor data for a single atom with the use of
   !> the verlet list, in case the atoms in the verlet list were not updated
   subroutine update_single_atom( this, my_indx, positions, lattice,&
                                  inverse_lattice, types, type_counter )

      use neighbor_list_mod, only : set_periodic_size_n_square 
      class( verlet_cell_list_type ), intent( inout ) :: this
      integer, intent( in ) :: my_indx
      real( kind = q ), intent( in ) :: positions( :, : )
      real( kind = q ), intent( in ) :: lattice( :, : )
      real( kind = q ), intent( in ) :: inverse_lattice( :, : )
      integer, intent( in ) :: types( : )
      integer, intent( inout ) :: type_counter( 1:this%number_types )

      real( kind = q ) :: pos( 1:3 )
      real( kind = q ) :: pos_neighbor( 1:3 )
      real( kind = q ) :: delta( 1:3 )

      real( kind = q  ) :: dist
      integer :: n_indx
      integer :: xx
      integer :: yy
      integer :: zz
      integer :: i_atom

      integer :: xyz_max( 1:3 )
      integer :: xyz_min( 1:3 )

      integer,allocatable :: n_list( : )
      integer,allocatable :: n_types( : )
      real( kind = q ), allocatable :: n_dist( : )
      real( kind = q ), allocatable :: n_vectors( :, : )
      integer :: n_counter 




      pos( 1:3 ) = positions( 1:3, my_indx )

      allocate( n_list( 1:size( this%cell_list%nn_atoms( my_indx )%neighbors ) ) )
      allocate( n_types( 1:size( this%cell_list%nn_atoms( my_indx )%neighbors ) ) )
      allocate( n_dist( 1:size( this%cell_list%nn_atoms( my_indx )%neighbors ) ) )
      allocate( n_vectors( 1:3, 1:size( this%cell_list%nn_atoms( my_indx )%neighbors ) ) )


      call set_periodic_size_n_square( this%cutoff, pos, inverse_lattice ,&
                                       xyz_max, xyz_min )

      n_counter     =  0
      type_counter  =  0
      do xx = xyz_min(1), xyz_max(1)
         do yy = xyz_min(2), xyz_max(2)
            do zz = xyz_min(3), xyz_max(3)
               do i_atom = 1, size( this%cell_list%nn_atoms( my_indx )%neighbors )
                  n_indx = this%cell_list%nn_atoms( my_indx )%neighbors( i_atom )
                  pos_neighbor( 1 )  =  positions( 1, n_indx )
                  pos_neighbor( 2 )  =  positions( 2, n_indx )
                  pos_neighbor( 3 )  =  positions( 3, n_indx )
                  
                  pos_neighbor( 1 ) =  pos_neighbor( 1 ) + real( xx, kind=q ) * lattice( 1, 1 ) + &
                                                           real( yy, kind=q ) * lattice( 1, 2 ) + &
                                                           real( zz, kind=q ) * lattice( 1, 3 )
                  pos_neighbor( 2 ) =  pos_neighbor( 2 ) + real( xx, kind=q ) * lattice( 2, 1 ) + &
                                                           real( yy, kind=q ) * lattice( 2, 2 ) + &
                                                           real( zz, kind=q ) * lattice( 2, 3 )
                  pos_neighbor( 3 ) =  pos_neighbor( 3 ) + real( xx, kind=q ) * lattice( 3, 1 ) + &
                                                           real( yy, kind=q ) * lattice( 3, 2 ) + &
                                                           real( zz, kind=q ) * lattice( 3, 3 )
                  
                  delta( 1 ) = pos( 1 ) - pos_neighbor( 1 ) 
                  delta( 2 ) = pos( 2 ) - pos_neighbor( 2 )
                  delta( 3 ) = pos( 3 ) - pos_neighbor( 3 )
                  
                  dist  =  delta(1) * delta(1) + delta(2) * delta(2) + delta(3) * delta(3)
                  if ( dist .le. this%cutoff_squared ) then
                     n_counter  =  n_counter + 1
                     n_list( n_counter )  =  n_indx
                     n_types( n_counter ) =  types( n_indx )
                     n_dist( n_counter )  =  dsqrt( dist ) 
                     n_vectors( 1, n_counter )  =  delta( 1 ) / n_dist( n_counter )
                     n_vectors( 2, n_counter )  =  delta( 2 ) / n_dist( n_counter )
                     n_vectors( 3, n_counter )  =  delta( 3 ) / n_dist( n_counter )
                     type_counter( types( n_indx ) )  =  type_counter( types( n_indx ) ) + 1
                  end if
               end do
            end do
         end do
      end do


      call this%nn_list( my_indx )%update_data_arrays( n_list, n_types, n_dist, n_vectors, n_counter )

   end subroutine update_single_atom



end module verlet_cell_list_mod

#else

module neighbor_list_mod
    type neighbor_list_type
    end type neighbor_list_type
end module neighbor_list_mod

module verlet_cell_list_mod
    type verlet_cell_list_type
    end type verlet_cell_list_type
end module verlet_cell_list_mod

#endif
