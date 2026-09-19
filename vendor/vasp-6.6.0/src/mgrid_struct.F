#include "symbol.inc"
MODULE mgrid_struct_def
  USE prec
  USE mpimy

  TYPE layout
     INTEGER NP                              !< local number of grid points
     INTEGER NALLOC                          !< allocation required
     INTEGER NFAST                           !< which index is fast (1-x, 2-y, 3-z)
     INTEGER NCOL                            !< number of columns in the grid
     INTEGER NROW                            !< number of elements in each column
     INTEGER,POINTER :: I2(:) => NULL()      !< y/z/x-index of each column
     INTEGER,POINTER :: I3(:) => NULL()      !< z/x/y-index of each column
     INTEGER,POINTER :: INDEX(:,:) => NULL() !< column index for each yz, zx or xy pair
  END TYPE layout

  TYPE grid_map
     INTEGER,POINTER :: N(:) => NULL()       !< number of elements send by each node
     INTEGER,POINTER :: PTR(:) => NULL()     !< sum_j=1,I N(j)
     INTEGER,POINTER :: TBL(:) => NULL()     !< address of each element
     ! inverse transformation (i.e. receiver information)
     INTEGER,POINTER :: NI(:) => NULL()
     INTEGER,POINTER :: PTRI(:) => NULL()
     INTEGER,POINTER :: TBLI(:) => NULL()
     LOGICAL LOCAL                !< all information is local
     LOGICAL LOCAL_COPY           !< no data redistribution required
  END TYPE grid_map
  
  !> Only  GRID
  !>
  !> @note
  !> REAL2CPLX determines wether a real to complex FFT is use
  !> LREAL     determines wether the data in real space are 
  !>           stored in real a complex array
  !>
  !> @note
  !> The decision whether a serial of parallel FFT is performed
  !> is presently decided by the GRID%RL%NFAST tag
  !> if GRID%RL%NFAST =1 -> serial FFT
  !> if GRID%RL%NFAST =3 -> parallel FFT 
  TYPE grid_3d
     INTEGER NGX          !< number of grid points in x
     INTEGER NGY          !< number of grid points in y
     INTEGER NGZ          !< number of grid points in z
     !> in the complex mode the _rd values
     !> are equal to NGX, NGY, NGZ
     !> if a real to complex FFT is used only half of the data
     !> are stored in one direction and the corresponding NG?_rd is
     !> set to (NG?+1)/2
     INTEGER NGX_rd 
     INTEGER NGY_rd !< see NGX_rd
     INTEGER NGZ_rd !< see NGX_rd
     INTEGER NPLWV                !< total number of grid points NGX*NGY*NGZ
     INTEGER MPLWV                !< allocation in complex words required to do in place FFT's
     INTEGER NGPTAR(3)            !< equivalent to /(NGX,NGY,NGZ/)
     INTEGER,POINTER :: LPCTX(:) => NULL()  !< loop counter in x
     INTEGER,POINTER :: LPCTY(:) => NULL()  !< loop counter in y
     INTEGER,POINTER :: LPCTZ(:) => NULL()  !< loop counter in z
     ! loop counters, in which the unbalanced contribution is zeroed
     INTEGER,POINTER :: LPCTX_(:) => NULL() !< loop counter in x
     INTEGER,POINTER :: LPCTY_(:) => NULL() !< loop counter in y
     INTEGER,POINTER :: LPCTZ_(:) => NULL() !< loop counter in z
     !> reciprocal space layout (x is always  fast index)
     TYPE(layout)    :: RC
     !> intermediate layout (y is always  fast index, used only in parallel version)
     TYPE(layout)    :: IN
     !> real space layout   (x or z is the fast index)
     TYPE(layout)    :: RL
     !> real space layout for FFT  (x or z is the fast index)
     !> this structure is usually equivalent to RL (and hence points to RL)
     !> only if the serial version is used for the FFT of wavefunctions,
     !> the structure differs from RL for the GRID_SOFT structure
     !> in this case, the FFT has z as fast index as required for the parallel FFT
     !> but the RL structure has x as fast index to be compatible to the FFT
     !> of wavefunctions
     TYPE(layout), POINTER  :: RL_FFT
     ! information only required for real space representation
     INTEGER NGZ_complex          !< number of grid points for z fast
     TYPE(grid_map)  :: RC_IN     !<  mapping for parallel version: recip -> intermediate
     TYPE(grid_map)  :: IN_RL     !<  mapping for parallel version: intermediate -> real space
     TYPE(communic), POINTER :: COMM,COMM_KIN,COMM_KINTER
     LOGICAL         :: LREAL     !< are data stored as complex or real numbers in real space
     LOGICAL         :: REAL2CPLX !< real to complex FFT
     REAL(q), POINTER:: FFTSCA(:,:) !< scaling factors for real to complex fft (wavefunction fft) in reciprocal space (only defined for PW in cutoff-sphere)
     !> for real space version, weight of each coefficient in reciprocal space for all coefficients
     !> plane waves with N1=0 (NX=0) within the cutoff sphere correspond to two
     !> G vectors on the 3D FFT grid
     REAL(q), POINTER:: FFTWEIGHT(:)     => NULL() 
     INTEGER, POINTER:: NINDPWCONJG(:)   => NULL() !< position where the conjugated coefficient needs to be stored
     INTEGER, POINTER:: IND_IN_SPHERE(:) => NULL() !< PW indices (within the sphere) with N1=0
  END TYPE grid_3d

  !
  !> Transition table used to go from a large to a small grid
  !> or vice versa
  !
  TYPE transit
     INTEGER,POINTER :: IND1(:) => NULL()   !< fast index transition table
     INTEGER,POINTER :: INDCOL(:) => NULL() !< column to column transition table
  END TYPE transit


  !> If LPLANE_WISE is set, the data are distributed in real and reciprocal
  !> space plane by plane i.e. one processor holds all elements of
  !> a plane with a specific x index
  !> this  reduces the communication in the FFT considerably
  !> the default for LPLANE_WISE can be set in this file (see below),
  !> or using the flag LPLANE in the INCAR reader
#ifdef plane_wise
  LOGICAL :: LPLANE_WISE=.TRUE.
#else
  LOGICAL :: LPLANE_WISE=.FALSE.
#endif

  !> compatibility modus to vasp.4.4
  !> the flag determine among other things whether the charge at unbalanced
  !> lattice vectors NGX/2+1 are zeroed or not
  LOGICAL :: LCOMPAT

  !> determines the effort FFTW will invest in constructing its "wisdom" (calls to FFTMAKEPLAN)
  !> FFTW_PLAN_EFFORT=0 -> FFTW_ESTIMATE, FFTW_PLAN_EFFORT=1 -> FFTW_MEASURE (default).
  INTEGER :: FFTW_PLAN_EFFORT=1

END MODULE
