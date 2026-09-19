# General

Offloading from a host (CPU) to devices (GPU) will involve the following:
* Data management
* Library calls
* Custom kernels
* Communication
* Flow control

## Data management

Data will need to be transferred to-and-fro the device.
A complication in this context is the extensive use of derived types in VASP: derived types that have members that are again derived types (or pointers to derived types!).
* Deep-copy of derived types

## Library calls

As much as possible/sensible, the heavy lifting is done by libraries: ``FFTW``, ``BLAS``, and ``(Sca)LAPACK``.
The library calls need to be appropriately re-routed for device-sided execution when applicable (see [flow control](#flow-control)).
Note that these libraries will be hardware-dependent.
* NVIDIA: ``cuFFT``, ``cuBLAS``, ``cuSolver``, ...
* AMD: ``rocFFT``, ``rocBLAS``, ...
* Intel: ``MKL``

## Custom kernels

There always remains a part of the computational work that can not be efficiently taken care of by libraries.
The custom compute kernels that are written in this context need to be ported to run device-side.
* Refactored kernels
* Batched versions

## Communication

VASP is parallelized using MPI.
The possibilty to communicate device-present data directly between devices, without first transferring it to the host is crucial for performance.
Library solutions to this effect will be again be hardware dependent.
* NVIDIA: ``CUDA-aware OpenMPI``, ``NCCL``
* AMD: ``GPU-aware MPI (CRAY-MPICH, MVAPICH and OpenMPI support that)``, ``RCCL``

## Flow control

* Runtime control (switch on/off) of device-sided execution
* Asynchrounous execution of offloaded kernels (OpenACC)

# Runtime execution control
* `offload_struct.F`
  ```fortran
  #include "symbol.inc"
  #ifdef OFFLOADING
      MODULE moffload_struct_def
        ...

        !> Flags whether device-sided execution is allowed in principle
        LOGICAL, PUBLIC :: OFFLOAD_ALLOW=.FALSE.
        !> Flags whether device-sided execution is momentarily activated
        LOGICAL, PUBLIC :: OFFLOAD_ON=.FALSE.

      END MODULE moffload_struct_def
  #endif // OFFLOADING
  ```

* `offload.F`
  ```fortran
  #undef PUSH_OFFLOAD_ON
  #undef POP_OFFLOAD_ON
      MODULE moffload_control

        USE mooffload_struct_def, ONLY : OFFLOAD_ON, OFFLOAD_ALLOW

        PRIVATE

        PUBLIC :: PUSH_OFFLOAD_ON, POP_OFFLOAD_ON

        INTEGER, PARAMETER :: MAXLEVEL=20
        INTEGER :: OFFLOAD_ON_LEVEL=0
        LOGICAL :: OFFLOAD_ON_STACK(MAXLEVEL)=.FALSE.

        CONTAINS

  !****************** SUBROUTINE PUSH_OFFLOAD_ON *************************
  !
  !***********************************************************************

        SUBROUTINE PUSH_OFFLOAD_ON(VAR)
        USE tutor, ONLY : vtutor
        LOGICAL :: VAR
  !$ACC WAIT
        IF (OFFLOAD_ON_LEVEL==MAXLEVEL) THEN
           CALL vtutor%error("PUSH_OFFLOAD_ON: ERROR: stack is full")
        ENDIF
        OFFLOAD_ON_LEVEL=OFFLOAD_ON_LEVEL+1
        OFFLOAD_ON_STACK(OFFLOAD_ON_LEVEL)=OFFLOAD_ON
        OFFLOAD_ON=VAR.AND.OFFLOAD_ALLOW
        END SUBROUTINE PUSH_OFFLOAD_ON


  !****************** SUBROUTINE POP_OFFLOAD_ON **************************
  !
  !***********************************************************************

        SUBROUTINE POP_OFFLOAD_ON
        USE tutor, ONLY : vtutor
  !$ACC WAIT
        IF (OFFLOAD_ON_LEVEL==0) THEN
           CALL vtutor%error("POP_OFFLOAD_ON: ERROR: stack is empty")
        ENDIF
        OFFLOAD_ON=OFFLOAD_ON_STACK(OFFLOAD_ON_LEVEL)
        OFFLOAD_ON_LEVEL=OFFLOAD_ON_LEVEL-1
        END SUBROUTINE POP_OFFLOAD_ON

      END MODULE moffload_control


      MODULE moffload
        ...
        USE moffload_control
        ...
      END MODULE mofload
  ```

* Switching offloading on/off at runtime:
  ```fortran
  #ifdef OFFLOADING
    USE moffload
  #endif
    ...

    ! Save current offload_on to stack, and switch offloading "on"
    PUSH_OFFLOAD_ON( .TRUE. )
    ...

    ! Return offload_on to previous state
    POP_OFFLOAD_ON
    ...
  ```
* `symbol.inc`
  ```c
  #ifdef OFFLOADING
  #define PUSH_OFFLOAD_ON(x)  CALL PUSH_VAR_OFFLOAD_ON(x)
  #define POP_OFFLOAD_ON      CALL POP_VAR_OFFLOAD_ON
  #else
  #define PUSH_OFFLOAD_ON
  #define POP_OFFLOAD_ON
  #endif
  ```

# Programming models

For the moment we only consider two directive-based parallel-programming extensions to Fortran:
* OpenACC
* OpenMP

Our efforts with OpenACC have progressed quite far.
Unfortunately, however, OpenACC has not been (and will not be) adopted by all major compiler/hardware vendors, meaning we will either have to switch to OpenMP or support both.
At the moment, however, switching to OpenMP is not an option as there are currently no compilers with a mature OpenMP offloading capability (Intel, AMD, Cray, and GNU are all working on it, though).

The case of OpenMP is complicated:
* Firstly, by the fact that it will not solely be used for offloading, it is already in use for CPU-sided threading.
These two cases have to be distinguished in some manner.
* And furthermore, OpenMP (targeting the CPU) is used in combination with OpenACC offloading.
In the latter case some of the OpenMP directives that target the CPU have to be skipped, because the kernels they affect are offloaded.

Possibly the directives will have to be specialized with respect to the target hardware as well ...

## Directives
In addition to the standard OpenACC and OpenMP Fortran keys `!$ACC`, `!$OMP`, and `!$`, we define the following precompiler macros:

* `CPOMP`, `GPOMP`, `DOOMP`, and `NOOMP`
* `CPACC`, `GPACC`, `DOACC`, and `NOACC`
* `__CPU` and `__GPU`

The definition and meaning of these macros is explained in the [following section](#precompiler-definitions).

## Precompiler definitions
* `symbol.inc`
  ```c
  #if defined(ACC_OFFLOAD) || defined(OMP_OFFLOAD)
  #define OFFLOADING
  #endif

  #ifdef OFFLOADING
  ! Line is included when   OFFLOADING
  #define __GPU
  ! Line is a comment when  OFFLOADING
  #define __CPU          !!
  #else
  ! Line is included when  !OFFLOADING
  #define __CPU
  ! Line is a comment when !OFFLOADING
  #define __GPU          !!
  #endif

  #if defined(_OPENACC) && !defined(ACC_OFFLOAD)
  ! Line is included when   _OPENACC and !ACC_OFFLOAD
  #define CPACC
  #else
  ! Line is a comment when !_OPENACC or   ACC_OFFLOAD
  #define CPACC          !!
  #endif

  #ifdef ACC_OFFLOAD
  ! As of this point _OPENACC will be set when ACC_OFFLOAD
  #define _OPENACC
  ! Line is included when   ACC_OFFLOAD
  #define GPACC
  #else
  ! Line is a comment when !ACC_OFFLOAD
  #define GPACC          !!
  #endif

  #ifdef _OPENACC
  ! Line is included when   _OPENACC
  #define DOACC
  ! Line is a comment when  _OPENACC
  #define NOACC          !!
  #else
  ! Line is included when  !_OPENACC
  #define NOACC
  ! Line is a comment when !_OPENACC
  #define DOACC          !!
  #endif

  #if defined(_OPENMP) && !defined(OMP_OFFLOAD)
  ! Line is included when   _OPENMP and !OMP_OFFLOAD
  #define CPOMP
  #else
  ! Line is a comment when !_OPENMP or   OMP_OFFLOAD
  #define CPOMP          !!
  #endif

  #ifdef OMP_OFFLOAD
  ! As of this point _OPENMP will be set when OMP_OFFLOAD
  #define _OPENMP
  ! Line is included when   OMP_OFFLOAD
  #define GPOMP
  #else
  ! Line is a comment when !OMP_OFFLOAD
  #define GPOMP          !!
  #endif

  #ifdef _OPENMP
  ! Line is included when   _OPENMP
  #define DOOMP
  ! Line is a comment when  _OPENMP
  #define NOOMP          !!
  #else
  ! Line is included when  !_OPENMP
  #define NOOMP
  ! Line is a comment when !_OPENMP
  #define DOOMP          !!
  #endif
  ```

### Meaning:
Lines marked with the following macros are included:

`CPACC`: when using OpenACC to target host only (currently unused)

`GPACC`: when using OpenACC to target device as well

`DOACC`: when using OpenACC

`NOACC`: when not using OpenACC

`CPOMP`: when using OpenMP to target host only

`GPOMP`: when using OpenMP to target device as well

`DOOMP`: when using OpenMP

`NOOMP`: when not using OpenMP

`__CPU`: no offloading at all

`__GPU`: OpenACC and/or OpenMP offloading

Otherwise these macros translate to `!!` (and comment out the line).

**N.B.I**: This list of macros could potentially be expanded to target specific device hardware: *e.g.* introducing `GPOMP_X` and `GPOMP_Y`.
Let's hope this will not be necessary, though.

**N.B.II**: Since we (currently) use OpenACC only for offloading purposes, we could have dropped the `CPACC` and `GPACC` macros: *i.e.*, we do not need `CPACC`, and `DOACC` can be (and is) used instead of `GPACC` if needed.
However, to be consistent with the macro definitions for OpenMP, I chose to include `CPACC` and `GPACC` as well.

### Common use cases:
* Both OpenACC as well as OpenMP offloading directives:
  ```fortran
  GPACC !$ACC .. IF(offload_on)
  GPOMP !$OMP .. IF(offload_on)
    ...
    some kernel to be offloaded
    ...
  ```
  These directives will be active when `-DACC_OFFLOAD` or `-DOMP_OFFLOAD`, respectively.
* OpenACC offloading directives and OpenMP directives that should only be included when compiling *without* OpenACC support:
  ```fortran
  NOACC !$OMP ...
  GPACC !$ACC ... IF(offload_on)
     ...
     some kernel to be offloaded
     ...
  ```
* OpenACC offloading directives and OpenMP directives that should only be included when compiling *with* OpenACC support:
  ```fortran
  DOACC !$OMP ...
     ...
     some code executed in parallel on the host
     ...
  GPACC !$ACC ... IF(offload_on)
     ...
     some kernel to be offloaded
     ...
  ```
  **N.B.I**: This covers for instance the case where we use OpenMP CPU-side to schedule work into different OpenACC asynchronous execution queues (to hide launch latency).

  **N.B.II**: As mentioned [above](#meaning), currently the `GPACC` macro is superfluous, since we use OpenACC only for offloading purposes.
  The following would accomplish the same goal:
  ```fortran
  DOACC !$OMP ...
     ...
     some code executed in parallel on the host
     ...
  !$ACC ... IF(offload_on)
     ...
     some kernel to be offloaded
     ...
  ```
  This is the kind of construct we use in the current OpenACC-only GPU port of VASP.

* OpenACC and OpenMP offloading directives, and OpenMP directives when compiling *without* any offloading support.
  ```fortran
  __CPU !$OMP ...
  GPOMP !$OMP ... IF(offload_on)
  GPACC !$ACC ... IF(offload_on)
    ...
    some kernel to be offloaded or executed in parallel on the host
    ...
  ```

None of the above is particularly clever.
It just replaces a lot of nested `#ifdef` clauses that would seriously mess up the readability of the code.

## Runtime execution control

Device-sided execution of all offloading regions (OpenACC and OpenMP) must be made conditional upon the value of `offload_on`:

```fortran
#ifdef OFFLOADING
  use moffload_struct_def, only :: offload_on
#endif

  ...
!$OMP ... IF(offload_on)
  ...
!$ACC ... IF(offload_on)
  ...
```

# Library calls

## BLAS, LAPACK
### Wrappers
BLAS and LAPACK library calls are routed through wrappers, defined in `blas_wrappers.F` and `lapack_wrappers.F`, and respectively.
The way this works out in practice is shown below using BLAS `ZCOPY` as an example.

* `blas_wrappers.F` contains the wrapper `WZCOPY`:
  ```fortran
  #include "symbol.inc"
        ...

  !************************* SUBROUTINE WZCOPY ***************************
  !
  !>
  !
  !***********************************************************************
  #undef WZCOPY
        SUBROUTINE WZCOPY(N,ZX,INCX,ZY,INCY)

        USE prec
  __GPU USE moffload_struct_def
  __GPU USE moffload_blas, ONLY : ZCOPY_OFFLOAD

        IMPLICIT NONE

        COMPLEX(q) :: ZX(*),ZY(*)
        INTEGER :: N,INCX,INCY

  __GPU IF (OFFLOAD_ON) THEN
  __GPU    CALL ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
  __GPU    RETURN
  __GPU ENDIF

        CALL ZCOPY(N,ZX,INCX,ZY,INCY)

        RETURN
        END SUBROUTINE WZCOPY
        ...
  ```
  where the lines related to *offloading* are protected by the precompiler `__GPU` macro. This means that `WZCOPY` will just call `ZCOPY` when the code is compiled **without** offloading support.

* Where `moffload_blas` is defined in `offload.F`:
  ```fortran
      MODULE moffload_blas
  #ifdef NVCUDA
        USE mcublas, BLAS_OFFLOAD_INIT   => CUBLAS_INIT, &
                     BLAS_OFFLOAD_FINISH => CUBLAS_FINISH
  #endif // NVCUDA
  #ifdef CRAYHIP
        USE mrocblas, BLAS_OFFLOAD_INIT   => ROCBLAS_INIT, &
                      BLAS_OFFLOAD_FINISH => ROCBLAS_FINISH
  #endif // CRAYHIP
       ...

    END MODULE moffload_blas
  ```

* One of the implementations of `ZCOPY_OFFLOAD` is found in the module `mrocblas` in `crayhip.F`:
  ```fortran
  #include "symbol.inc"
  #ifdef CRAYHIP
  #ifdef OFFLOADING
      MODULE mrocblas
        ...

        TYPE(c_ptr) :: ROCBLAS_HANDLE

        INTERFACE
           ...

           SUBROUTINE HIP_ZCOPY(PTR, N, ZX, INCX, ZY, INCY) BIND(C)
              USE iso_c_binding
              TYPE(c_ptr)                      :: PTR
              INTEGER(c_int), VALUE            :: N, INCX, INCY
              TYPE(c_ptr), VALUE               :: ZX, ZY
           END SUBROUTINE HIP_ZCOPY
           ...

        END INTERFACE

      CONTAINS
        ...

        SUBROUTINE ZCOPY_OFFLOAD(N,ZX,INCX,ZY,INCY)
        IMPLICIT NONE
        COMPLEX(q), TARGET :: ZX(*),ZY(*)
        INTEGER :: N,INCX,INCY

  GPOMP !$OMP TARGET DATA USE_DEVICE_PTR(ZX,ZY)
  GPOMP CALL HIP_ZCOPY(ROCBLAS_HANDLE,N,c_loc(ZX),INCX,c_loc(ZY),INCY)
  GPOMP !$OMP END TARGET DATA
        END SUBROUTINE ZCOPY_OFFLOAD
        ...

      END MODULE mrocblas
        ...

  #endif // OFFLOADING
  #endif // CRAYHIP
  ```
  Note that the compiler will see this particular implementation of `ZCOPY_OFFLOAD` when the code is compiled with `-DOMP_OFFLOAD -DCRAYHIP`. Another implementation of the same routine can be found in *e.g.* `nvcuda.F` and will be used when compiling with `-DACC_OFFLOAD -DNVCUDA`.

* As shown above, `HIP_ZCOPY` in `crayhip.F` is an interface to a C subroutine. The C subroutine `HIP_ZCOPY` is defined in `src/HIP/rocblas_interfaces.cpp` and is compiled using `hipcc`. It is just a wrapper around `rocblas_zcopy`:
  ```c
  extern "C" {
      ...

      void hip_zcopy(void *ptr, int N, double _Complex *zx, int incx, double _Complex *zy, int incy) {
          rocblas_handle *handle = (rocblas_handle *) ptr;
          rocblas_double_complex *zx2 = reinterpret_cast<rocblas_double_complex*>(zx);
          rocblas_double_complex *zy2 = reinterpret_cast<rocblas_double_complex*>(zy);
          rocblas_zcopy(*handle, N, zx2, incx, zy2, incy);
      ...

  }
  ```

* Porting a `ZCOPY`call to GPU: just it by a call to the corresponding wrapper `WZCOPY`:
  ```fortran
  #include "symbol.inc"
        ...

        CALL WZCOPY(N,ZX,INCX,ZY,INCY)
        ...
  ```
* The following definition in `symbol.inc`:
  ```c
  ...

  #ifdef OFFLOADING
  ...

  #else
  ...

  #define WZCOPY            ZCOPY
  ...

  #endif
  ...
  ```

  will re-substitute the wrapper calls with the original library call when the code is compiled **without** offloading support. This is not strictly necessary, of course: we could decide to call wrappers per default.o

### Initialization

## FFT
### Planning
### Execution

## MPI

## Building
The code is then build with
```
-DACC_OFFLOAD -DUSE_ACC_LIBx
```
or
```
-DOMP_OFFLOAD -DUSE_OMP_LIBx
```

**N.B.**: There are no additional function calls when neither `ACC_OFFLOAD` nor `OMP_OFFLOAD` are set.

# Custom kernels

Wrap the custom kernel calls like [library calls](#library-calls-1):

Supposing `src_file.F` contains a module `m_some_module` with a subroutine `some_subroutine` that needs to be substantially modified for effective offloading (*e.g* re-ordering loops).

* Add wrapper `wsome_subroutine` and implement offloaded version `some_subroutine_offload` directly in `src_file.F` (preferred solution):
  ```fortran
        SUBROUTINE SOME_SUBROUTINE( ARGS )
        ...
        ...
        END SUBROUTINE SOME_SUBROUTINE

  #ifdef OFFLOADING
        ! Implementation
        SUBROUTINE SOME_SUBROUTINE_OFFLOAD( ARGS )
  #ifdef ACC_OFFLOAD
    ...
    ...
  #endif
  #ifdef OMP_OFFLOAD
    ...
    ...
  #endif
        END SOME_SUBROUTINE_OFFLOAD

        ! Wrapper
        SUBROUTINE WSOME_SUBROUTINE( ARGS )

        USE MOFFLOAD_STRUCT_DEF, ONLY : OFFLOAD_ON

        IF (OFFLOAD_ON) THEN
           CALL SOME_SUBROUTINE_OFFLOAD( ARGS )
           RETURN
        ENDIF

        CALL SOME_SUBROUTINE( ARGS )

        END SUBROUTINE WSOME_SUBROUTINE
  #endif // OFFLOADING
  ```

* Replace the call to `some_subroutine` by a call to the wrapper `wsome_subroutine`:

  ```fortran
  #include "symbol.inc"

        CALL WSOME_SUBROUTINE()
  ```

* And add the following definition to `symbol.inc` to avoid calling the wrapper in case the code is built without offloading support:
  ```c
  #ifdef OFFLOADING
  ...

  #else
  ...

  #define WSOME_SUBROUTINE  SOME_SUBROUTINE
  ...

  #endif
  ```

# Exceptions and special cases
For two libraries we already have wrappers in place (albeit messy and in need of refactoring):
* MPI

# Modules and files
