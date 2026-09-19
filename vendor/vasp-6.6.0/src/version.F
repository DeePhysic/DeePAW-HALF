#include "symbol.inc"
!> Stores information about the used VASP version.
!>
!> The information about the version is also written to the HDF5 output and can
!> be used to verify that certain functionality is available.
module version

   implicit none

   integer, parameter :: major = 6  !< Major version number of VASP
   integer, parameter :: minor = 6  !< Minor version number of VASP
   integer, parameter :: patch = 0  !< Patch version number of VASP
   character(len=*), parameter :: date = '06Mar2026'  !< When this version of VASP was released

   type date_t
      integer year, month, day
   end type date_t
   type(date_t), parameter :: release_date = date_t(year = 2026, month = 3, day = 6)

contains

   !> Generate a short string describing the VASP version based on the version
   !> number and the preprocessor options.
   pure character(len=80) function vasp()
      use string, only: str
      character(len=:), allocatable :: label, build, version_
#ifdef VASP_VERSION
      label = VASP_VERSION
#else
      label = 'vasp.' // str(major) // '.' // str(minor) // '.' // str(patch) // ' ' // date
#endif
      build = '(build ' // __DATE__ // ' ' // __TIME__ // ')'
#ifdef gammareal
      version_ = 'gamma-only'
#else
      version_ = 'complex'
#endif
      vasp = label // ' ' // build // ' ' // version_
   end function vasp

end module version
