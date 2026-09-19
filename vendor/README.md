# Private VASP source

`vasp-6.6.0/` contains the VASP source supplied by the repository owner and the
DeePAW-HALF integration changes tested on the Pro 6000 host. This directory is
intended only for the private `vasp-6.6-half-integration` branch and must remain
subject to the owner's VASP license.

Do not build in this source directory. From the repository root run:

```bash
tools/half-cmake vasp all
```

CMake copies the tree under `build/`, generates the machine-specific
`makefile.include`, and writes the executable to
`build/vasp-cuda<CUDA>-cc<CC>-release/vasp-6.6.0/bin/vasp_std`.
