#include "symbol.inc"

module build_info

   implicit none
   private

   public :: cpp_options, link_line

#include "build_info.inc"

#ifdef BUILD_INFO
   public :: dump_build_info

   contains

   !==============================================================================
   !> Prints build information to the provided file.
   !==============================================================================
   subroutine dump_build_info(io)

      use string, only: str
      use base, only: in_struct
#ifdef MPI
      use mpi_f08
#endif

      !> Output file.
      type(in_struct), intent(in) :: io

      character(len=5)            :: is_defined ! "true" or "false".
      character(len=*), parameter :: pp_fmt = '("   ",A,": ",A)'

      if ( io%iu6 >= 0 ) then
         write(io%iu6, '(A)') &
            ' ============================= BUILD_INFO ============================== '
         ! Verbosity.
         write(io%IU6, '(" ",A,": ",I0)') 'verbosity', BUILD_INFO

         ! Build date and time.
         write(io%IU6, '(" ",A,": ",A)') 'build date', &
            __DATE__ // ' ' // __TIME__

         ! Main contents of makefile.include.
         write(io%IU6, '(" ",A,": ")') 'makefile.include'
         write(io%IU6, '("   ",A,": ",A)') 'CPP_OPTIONS', cpp_options
         write(io%IU6, '("   ",A,": ",A)') 'FC', fc
         write(io%IU6, '("   ",A,": ",A)') 'FCL', fcl
         write(io%IU6, '("   ",A,": ",A)') 'FFLAGS', fflags
         write(io%IU6, '("   ",A,": ",A)') 'LLIBS', llibs
         write(io%IU6, '("   ",A,": ",A)') 'INCS', incs

         ! Common precompiler flags.
         write(io%IU6, '(" ",A,": ")') 'preprocessor'
         write(io%IU6, '("   ",A,": ",A)') 'HOST', HOST
         write(io%IU6, '("   ",A,": ",I0)') 'MPI_BLOCK', MPI_BLOCK
         write(io%IU6, '("   ",A,": ",I0)') 'CACHE_SIZE', CACHE_SIZE
#ifdef MPI
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "MPI", is_defined; is_defined = "false"
#ifdef use_collective
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "use_collective", is_defined; is_defined = "false"
#ifdef scaLAPACK
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "scaLAPACK", is_defined; is_defined = "false"
#ifdef avoidalloc
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "avoidalloc", is_defined; is_defined = "false"
#ifdef vasp6
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "vasp6", is_defined; is_defined = "false"
#ifdef tbdyn
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "tbdyn", is_defined; is_defined = "false"
#ifdef fock_dblbuf
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "fock_dblbuf", is_defined; is_defined = "false"
#ifdef _OPENMP
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "_OPENMP", is_defined; is_defined = "false"
#ifdef ACC_OFFLOAD
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "ACC_OFFLOAD", is_defined; is_defined = "false"
#ifdef USENCCL
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "USENCCL", is_defined; is_defined = "false"
#ifdef qd_emulate
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "qd_emulate", is_defined; is_defined = "false"
#ifdef VASP_HDF5
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "VASP_HDF5", is_defined; is_defined = "false"
#ifdef use_shmem
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "use_shmem", is_defined; is_defined = "false"
#ifdef shmem_bcast_buffer
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "shmem_bcast_buffer", is_defined; is_defined = "false"
#ifdef shmem_rproj
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "shmem_rproj", is_defined; is_defined = "false"
#ifdef sysv
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "sysv", is_defined; is_defined = "false"
#ifdef VASP2WANNIER90
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "VASP2WANNIER90", is_defined; is_defined = "false"
#ifdef libbeef
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "libbeef", is_defined; is_defined = "false"
#ifdef DFTD4
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "DFTD4", is_defined; is_defined = "false"
#ifdef SDFTD3
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "SDFTD3", is_defined; is_defined = "false"
#ifdef LIBMBD
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "LIBMBD", is_defined; is_defined = "false"
#ifdef PROFILING
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "PROFILING", is_defined; is_defined = "false"
#ifdef USELIBXC
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "USELIBXC", is_defined; is_defined = "false"
#ifdef LAPACK36
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "LAPACK36", is_defined; is_defined = "false"
#ifdef CUSOLVERMP
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "CUSOLVERMP", is_defined; is_defined = "false"
#ifdef CUBLASMP
         is_defined = "true"
#endif
         write(io%IU6, pp_fmt) "CUBLASMP", is_defined; is_defined = "false"

#ifdef MPI
         ! MPI information.
         write(io%IU6, '(" ",A,": ")') 'mpi'
         write(io%IU6, '("   ",A,": ",I0)') "MPI_VERSION", MPI_VERSION
         write(io%IU6, '("   ",A,": ",I0)') "MPI_SUBVERSION", MPI_SUBVERSION
#endif

         write(io%IU6,'(A)') &
            ' ============================= BUILD_INFO ============================== '
      endif

   end subroutine dump_build_info
#endif
end module build_info
