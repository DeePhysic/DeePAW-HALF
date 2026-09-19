#include "symbol.inc"

#ifdef VASP_HDF5

module vhdf5_struct_def
  use hdf5

  integer, parameter :: VH5_MAX_RANK = 10

  character(LEN=*), parameter :: VASPIN = "vaspin.h5"
  character(LEN=*), parameter :: VASPOUT = "vaspout.h5"
  character(LEN=*), parameter :: VASPWPOTFULL = "vaspwpotfull.h5"
  character(LEN=*), parameter :: VASPWPOT = "vaspwpot.h5"

  character(LEN=*), parameter :: GRP_INFOHEADER = "infoheader"
  character(LEN=*), parameter :: GRP_GRID = "grid"
  character(LEN=*), parameter :: GRP_EIGENVAL = "electron_eigenvalues"
  character(LEN=*), parameter :: GRP_ELECTRON = "electron"
  character(LEN=*), parameter :: GRP_PHONONS = "phonons"
  character(LEN=*), parameter :: GRP_ELPHON = 'electron_phonon'
  character(LEN=*), parameter :: GRP_WPOT = "wpot"
  character(LEN=*), parameter :: GRP_DOS = "electron_dos"
  character(LEN=*), parameter :: GRP_NMRCUR = "nmrcur_meta"
  character(LEN=*), parameter :: GRP_PAW = "paw"
  character(LEN=*), parameter :: GRP_POSITIONS = "positions"
  character(LEN=*), parameter :: GRP_HISTORY = "ion_dynamics"
  character(LEN=*), parameter :: GRP_CHARGE = "charge"
  character(LEN=*), parameter :: GRP_TAU = "kinetic_energy_density"
  character(LEN=*), parameter :: GRP_CONDUCTIVITY = "conductivity"
  character(LEN=*), parameter :: GRP_LINEAR_RESPONSE = "linear_response"
  character(LEN=*), parameter :: GRP_DIELECTRIC = "dielectric"
  character(LEN=*), parameter :: GRP_PROJECTORS = "projectors"
  character(LEN=*), parameter :: GRP_LOCPROJ = "locproj"
  character(LEN=*), parameter :: GRP_WANPROJ = "wanproj"
  character(LEN=*), parameter :: GRP_PAIR_CORRELATION = "pair_correlation"
  character(LEN=*), parameter :: GRP_POTENTIAL = "potential"
  character(LEN=*), parameter :: GRP_LOCPOT = "locpot"
  character(LEN=*), parameter :: GRP_LCAO = "lcao"
  character(LEN=*), parameter :: GRP_RHFATOM = "rhfatom"
  character(LEN=*), parameter :: GRP_WAVE = "wave"
  character(LEN=*), parameter :: GRP_WAVEDER = "waveder"
  character(LEN=*), parameter :: GRP_INTERMEDIATE = "intermediate"
  character(LEN=*), parameter :: GRP_RESULTS = "results"
  character(LEN=*), parameter :: GRP_INPUT = "input"
  character(LEN=*), parameter :: GRP_ORIGINAL = "original"
  character(LEN=*), parameter :: GRP_VERSION = "version"
  character(LEN=*), parameter :: GRP_POT = "pot"
  character(LEN=*), parameter :: GRP_OEP = "oep"
  character(LEN=*), parameter :: GRP_SELFENERGY = "selfenergy"
  character(LEN=*), parameter :: SUBGRP_KPOINTS = "kpoints"
  character(LEN=*), parameter :: SUBGRP_KPOINTS_LINE = "kpoints_along_line"
  character(LEN=*), parameter :: SUBGRP_INCAR = "incar"
  character(LEN=*), parameter :: SUBGRP_POSCAR = "poscar"
  character(LEN=*), parameter :: SUBGRP_POTCAR = "potcar"
  character(LEN=*), parameter :: SUBGRP_PARCHG = "partial_charges"
  character(LEN=*), parameter :: SUBGRP_CRPA = "crpa"
  character(LEN=*), parameter :: SUBGRP_GW = "gw"

  integer(HID_T):: ih5intermediategroup_id = -1
  integer(HID_T):: ih5infileid = -1
  integer(HID_T):: ih5outfileid = -1
  integer(HID_T):: ih5wavefileid = -1
  integer(HID_T):: ih5wpotfileid = -1
  integer(HID_T):: ih5wpotfullfileid = -1
  integer(HID_T):: ih5ininputgroup_id = -1
  integer(HID_T):: ih5outinputgroup_id = -1
  integer(HID_T):: ih5outoriginalgroup_id = -1

  logical :: incar_found
  logical :: hdf5_found
  logical :: kpoints_found
  logical :: kpoints_hdf5_found
  logical :: poscar_found
  logical :: potcar_found       ! POTCAR found
  logical :: potcar_h5_found  ! POTCAR.h5 found
  logical :: vaspinh5_potcar_found  ! POTCAR dataset in vaspin.h5 found
  logical :: wavecar_found
  logical :: vaspwave_found = .false.
  logical :: wpot_found = .false.
  logical :: wpotfull_found = .false.

  logical :: lwriteh5 =.FALSE.
  logical :: lsynch5 = .FALSE. !< should the hdf5 file be up to date or is buffering fine
  ! TODO: this variable should be eliminated in favor of passing down IO

  integer,parameter :: vh5_default_chunk_size = 100

  integer(HID_T) :: create_inter_prop

end module vhdf5_struct_def
#else

module vhdf5_struct_def
   contains
   subroutine vhdf5_dummy
   write(*,*)'Im a DEC compiler so I need this line'
   end subroutine vhdf5_dummy
end module vhdf5_struct_def

#endif

