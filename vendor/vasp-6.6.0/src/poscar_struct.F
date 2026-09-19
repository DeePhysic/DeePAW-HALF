#include "symbol.inc"
   MODULE poscar_struct_def
      USE prec
      IMPLICIT NONE 
!
! poscar description input file
! only included if MODULES are not supported
!
      TYPE type_info
!only T_INFO
        CHARACTER*40 SZNAM2           !< name of poscar file
        INTEGER NTYPD                 !< dimension for types
        INTEGER NTYP                  !< number of types
        INTEGER NTYPPD                !< dimension for types inc. empty spheres
        INTEGER NTYPP                 !< number of types empty spheres
        INTEGER NIOND                 !< dimension for ions
        INTEGER NIONPD                !< dimension for ions inc. empty spheres
        INTEGER NIONS                 !< actual number of ions
        INTEGER NIONP                 !< actual number of ions inc. empty spheres
        LOGICAL LSDYN                 !< selective dynamics (yes/ no)
        LOGICAL LDIRCO                !< positions in direct/recproc. lattice
        REAL(q), POINTER :: POSION(:,:)=> NULL() !< positions usually same as DYN%POSION
        LOGICAL,POINTER ::  LSFOR(:,:) => NULL() !< selective dynamics
        INTEGER, POINTER :: ITYP(:)    => NULL() !< type for each ion
        INTEGER, POINTER :: NITYP(:)   => NULL() !< number of ions for each type
        REAL(q), POINTER :: POMASS(:)  => NULL() !< mass for each ion type
        REAL(q), POINTER :: RWIGS(:)   => NULL() !< wigner seitz radius for each ion type
        REAL(q), POINTER :: ROPT(:)    => NULL() !< optimization radius for each type
        REAL(q), POINTER :: ATOMOM(:)  => NULL() !< initial local spin density for each ion
        REAL(q), POINTER :: DARWIN_R(:)=> NULL() !< parameter for darwin like mass term at each ion
        REAL(q), POINTER :: DARWIN_V(:)=> NULL() !< parameter for darwin like mass term at each ion 
        REAL(q), POINTER :: VCA(:)     => NULL() !< weight of each species for virtual crystal approximation 
        REAL(q), POINTER :: ZCT(:)     => NULL() !< "charge transfer" charges for non-scf calculations
        REAL(q), POINTER :: RGAUS(:)   => NULL() !< widths for Gaussian CT charge distributions
        CHARACTER (LEN=2), POINTER :: TYPE(:) => NULL() !< type information for each ion
#ifdef VASP_HDF5
        CHARACTER (LEN=64), POINTER :: SHA256(:) => NULL() !< POTCAR SHA for each type
        CHARACTER (LEN=20), POINTER :: TYPEF(:)  => NULL() !< full POTCAR name for each type (e.g Li_sv_GW)

        CHARACTER (LEN=64), POINTER :: SHA256_ON_POSCAR(:) => NULL() !< POTCAR SHA for each type read from POSCAR
        CHARACTER (LEN=20), POINTER :: TYPEF_ON_POSCAR(:)  => NULL() !< full POTCAR name for each type read from POSCAR
#endif
      END TYPE


      TYPE dynamics
!only DYN
        REAL(q), POINTER :: POSION(:,:) !< positions
        REAL(q), POINTER :: POSIOC(:,:) !< old positions
        REAL(q), POINTER :: EFOR(:,:)  !< external forces
        REAL(q), POINTER :: VEL(:,:)  !< velocities
        REAL(q), POINTER :: D2(:,:)   !< predictor corrector/coordinates
        REAL(q), POINTER :: D2C(:,:)  !< predictor corrector/coordinates
        REAL(q), POINTER :: D3(:,:)   !< predictor corrector/coordinates
        REAL(q) A(3,3)                !< current lattice (presently unused)
        REAL(q) AC(3,3)               !< old lattice (presently unused)
        REAL(q) :: SNOSE(4) = 0.0_q   !< nose thermostat
        INTEGER IBRION                !< mode for relaxation
        INTEGER ISIF                  !< mode for stress/ ionic relaxation
        REAL(q) POTIM                 !< time step
        REAL(q) EDIFFG                !< accuracy for ionic relaxation
        REAL(q), POINTER :: POMASS(:) !< mass of each ion for dynamics
        REAL(q) SMASS                 !< mass of nose thermostat
        REAL(q) PSTRESS               !< external pressure
        REAL(q) TEBEG, TEEND          !< temperature during run
        REAL(q) TEMP                  !< current temperature
        INTEGER NSW                   !< number of ionic steps
        INTEGER NBLOCK,KBLOCK         !< blocks
        INTEGER INIT                  !< predictore corrector initialized
        INTEGER NFREE                 !< estimated ionic degrees of freedom
      END TYPE

      ! formaly in lattice.inc
      TYPE latt
         REAL(q) :: SCALE             !< universal scaling factor that scales lattice vectors
         REAL(q) :: A(3,3)            !< fractional lattice vectors (already scaled)
         REAL(q) :: B(3,3)            !< reciprocal lattice vectors w/o factor 2pi (use e.g. GZ=(G1*B(3,1)+G2*B(3,2)+G3*B(3,3))*TPI)
         REAL(q) :: ANORM(3)          !< norm of fractional lattice vectors
         REAL(q) :: BNORM(3)          !< norm of reciprocal lattice vectors
         REAL(q) :: OMEGA             !< volume of unit cell
!tb start
         REAL(q) AVEL(3,3)            !< lattice velocities
         INTEGER INITlatv             !< lattice velocities initialized                  !
!tb end

      END TYPE

      TYPE nh_chains
        LOGICAL :: LINIT = .FALSE.
        INTEGER :: nchains = 0
        INTEGER :: nchainsmax = 20
        INTEGER :: period = 40            !c period of thermostat in time steps
        REAL(q) :: G(20) = 0._q
        REAL(q) :: Q(20) = 0._q !c masses of thermnostat DOFs 
        REAL(q) :: P(20) = 0._q !c momenta of thermnostat DOFs 
        REAL(q) :: X(20) = 0._q !c coordinates of thermnostat DOFs 
        INTEGER :: NS = 1       !c number of Yoshida grid points
        REAL(q), POINTER :: WS(:)  => NULL() !c weights of the Yoshida grid points
        INTEGER :: NRESPA = 1 !c number of RESPA steps

      END TYPE nh_chains


   END MODULE poscar_struct_def
