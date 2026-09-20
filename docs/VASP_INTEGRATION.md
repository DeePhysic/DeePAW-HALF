# VASP 6.6.0 direct initialization with DeePAW-HALF

This note describes an in-process integration. With the repository owner's
authorization, this private branch stores the user-supplied
`vasp-edge-release.6.6.0` source. It must not be redistributed from this
private repository to parties without a VASP license. The insertion points
below were verified and tested in that source tree.

The private `vasp-6.6-half-integration` branch contains the complete source and
a one-command build:

```bash
git switch vasp-6.6-half-integration
tools/half-cmake vasp all
```

CMake copies the private source into the build tree, generates a
`makefile.include` matching the detected NVHPC/CUDA/compute capability, builds
the FFTW3 wrapper from oneMKL sources, builds `libhalf.so`, and then links
`vasp_std`. The vendored source tree remains free of object, module, and binary
build artifacts.

## What VASP does at startup

The legacy restart path is split into two phases:

1. `main.F` calls `INWAV_HEAD` before wavefunction allocation. The header may
   change restart metadata and validates `NKPTS`, `NBANDS`, `ENCUT`, and the
   lattice.
2. After `ALLOCW`, `main.F` calls `INWAV_FAST`. For every spin, k point, and
   band, it reads the global serial plane-wave coefficients, remaps them when
   the lattice or cutoff changed, and calls `DIS_PW_BAND` to place them in the
   distributed `W%CPTWFP` array. It also fills `W%CELTOT` and `W%FERTOT`.

The PAW projector coefficients must not be imported from HALF. Later, VASP
calls `PROALL` and `ORTHCH`; those routines regenerate `W%CPROJ` with VASP's
own projectors and orthogonalize the imported orbitals with the PAW overlap.

Charge initialization is independent of wavefunction initialization. With
`ICHARG=2` or `3`, VASP calls its overlapping-atom density path (SAD). Therefore
importing HALF orbitals alone does **not** bypass SAD. A DeepAW density must be
provided as `CHGCAR` and selected with `ICHARG=1` (or `ICHARG=11` for a fixed
density). HALF then diagonalizes the Hamiltonian built from the same density.

## Direct-call contract

Before solving, VASP now sends its parsed host request through memory:

```c
int half_set_request_geometry(
    half_handle handle, const half_request_geometry_v1 *geometry,
    char *error, int error_capacity);
```

The adapter packs `LATT_CUR%A` as row-major lattice vectors,
`T_INFO%POSION` as atom-major fractional coordinates, `T_INFO%ITYP` as
one-based species indices, and uses the dense charge grid `GRIDC%NGPTAR` as
`grid[3]`. HALF deep-copies all fields. They are retained for a future direct
DeepAW API call and intentionally do not alter the current file-backed solve.
`half_get_request_geometry()` provides an auditable query/copy path.

Wavefunction transfer uses the additive ABI-v1 function:

```c
int half_solve_kpoint_mapped(
    half_handle handle, const double kpoint[3], int nbands,
    int64_t host_npw, const int32_t *host_gvec,
    double *eigenvalues_eV, double _Complex *eigenvectors, int64_t ld,
    double *overlap_min, double *overlap_max,
    char *error, int error_capacity);
```

`host_gvec` contains `(Gx,Gy,Gz)` triples in the host's desired serial order.
HALF verifies that both bases contain exactly the same vectors and returns each
eigenvector already permuted to that order. Basis-size mismatch, a missing G
vector, or duplicate host vectors is a hard error; silently truncated initial
states are never returned.

For VASP, construct the serial order from `wavedes1%PL_INDEX`, `%PL_COL`, and
the local `%IGX/%IGY/%IGZ` arrays. Pass the returned column for global band `n`
to VASP's existing `DIS_PW_BAND`. Set the corresponding global eigenvalue in
`W%CELTOT`; let VASP initialize occupations and run `PROALL`/`ORTHCH` normally.

## Recommended VASP control flow

Add a dedicated logical input tag, for example `LHALF_INIT`, rather than
overloading `ISTART`. Overloading `ISTART=1` incorrectly activates unrelated
restart assumptions.

Immediately after `ALLOCW`, use this branch:

```text
if LHALF_INIT:
    reject an unsupported calculation mode
    create HALF context from DeepAW CHGCAR + POTCAR
    transfer lattice + fractional positions + species + dense grid from memory
    for each irreducible k point:
        build VASP serial G-vector list
        call half_solve_kpoint_mapped
        call DIS_PW_BAND for every global band
        fill W%CELTOT
    destroy HALF context
    mark wavefunctions as initialized (skip WFINIT only)
else if ISTART > 0:
    existing WAVECAR path
else:
    existing WFINIT path
```

DeepAW smooth CHGCAR files normally contain the real-space density grid but
not the PAW onsite-occupancy tail written by VASP. Standard VASP 6.6.0 treats
that missing tail as failure of the entire CHGCAR and falls back to SAD. HALF
mode needs one narrowly scoped exception: after the structure and grid have
passed validation, preserve the loaded DeepAW smooth density and retain the
atomic PAW onsite matrix that VASP initialized with `DEPATO/SET_RHO_PAW` before
`READCH`. Do not zero that matrix or relax the check for ordinary `ICHARG=1`
runs.

Run with:

```text
LHALF_INIT = .TRUE.
ISTART     = 0
ICHARG     = 1
INIWAV     = 1
```

The new flag must suppress only `WFINIT`; it must not suppress the `CHGCAR`
read or the later `PROALL` and `ORTHCH` calls.

To let HALF obtain its reconstruction density directly from DeePAW-eSCN,
select the source explicitly in INCAR:

```text
LHALF_INIT = .TRUE.
LHALF_API = .TRUE.
HALF_ESCN_URL = http://127.0.0.1:8265
```

The adapter then passes its existing lattice, fractional positions, one-based
type indices, and dense grid to `half_create_from_escn()`. Atomic numbers come
from the matching POTCAR, so the VASP-to-HALF descriptor does not change. HALF
requests density only, converts the API's C-order grid, and explicitly
normalizes it to the POTCAR valence-electron total. `LHALF_API=.FALSE.` is the
default and always selects the CHGCAR constructor. When API mode is enabled,
the INCAR URL takes precedence; the environment variable of the same name is
used only as a compatibility fallback for an empty INCAR URL. Missing CHGCAR
does not silently enable the network path. This switch changes HALF's
reconstruction density only. VASP still loads CHGCAR for its own `ICHARG=1`
density until a separate host charge-grid injection is implemented.

The standalone adapter that can be copied into a licensed VASP source tree is
[`examples/vasp/half_vasp_init.F`](../examples/vasp/half_vasp_init.F). It is
original HALF bridge code and does not contain or redistribute VASP sources.

## Parallel and GPU rules

- First integration target: collinear full-complex VASP, `ISPIN=1`, `KPAR=1`.
  `NCORE>1` is supported by using `DIS_PW_BAND` instead of writing distributed
  storage directly.
- Call the dense HALF solve on one rank per k-point group. All ranks still call
  VASP's collective distribution routine. Broadcast HALF status and
  eigenvalues before continuing.
- Extend to `KPAR>1` by creating one HALF context per k-point-group root and
  assigning only that group's k points. Never let every MPI rank launch the
  same dense GPU solve.
- Reject `LNONCOLLINEAR`, spinors, `ISPIN=2`, and Gamma-reduced storage until
  their explicit spin/density and conjugate-pair contracts are implemented.
- A CUDA-enabled VASP executable should link `libhalf` with the same NVHPC/CUDA
  major toolchain. Device selection must follow the rank-local device chosen by
  VASP. The current ABI returns the completed coefficient block through host
  memory; this is correct but not yet zero-copy.

## Pro 6000 validation

The integration was run on 2026-09-19 with an RTX PRO 6000 Blackwell, NVHPC
26.5 and its bundled CUDA 13.2, and a full-complex VASP 6.6.0 build. The HfO2
case used 520 eV, a 54-by-54-by-54 density grid, 3407 plane waves, 56 bands,
and the Gamma point.

- One rank: 4.90 s wall time for HALF initialization plus one VASP Davidson
  step; first-step rms was 2.81e-4.
- Two ranks with `NCORE=2`: 4.66 s; energy and rms matched the single-rank run,
  validating the `DIS_PW_BAND` distribution path.
- Native SAD plus random `WFINIT`: 2.63 s, but its first-step rms was 158,
  about 5.6e5 times larger. The HALF overhead is a one-time startup cost, so this one-step run is not an SCF
  time-to-solution benchmark.
- Both HALF runs read 88.0000237 electrons, did not enter the overlapping-atom
  smooth-density path, and completed energy and force evaluation.

See [`docs/validation/vasp_hfo2_pro6000.json`](validation/vasp_hfo2_pro6000.json).
The build used the MKL FFTW3 wrapper from `/opt/intel/oneapi/mkl/2026.0`. The
actual CUDA compile toolkit was NVHPC 26.5's bundled 13.2; a documentation-only
`/usr/local/cuda-13.4` directory was not reported as a working toolkit.

A full SCF comparison used `KSPACING=0.35` (36 irreducible k points),
`LMAXPAW=-1`, `ALGO=All`, `EDIFF=1e-4`, and the 500 eV cutoff selected by
`PREC=High`. HALF converged in 4 electronic iterations versus 16 for SAD;
end-to-end wall times, including HALF initialization at all k points, were
193.44 s and 266.96 s. The detailed Chinese report, complete convergence
traces, checksums, and limitations are in
[`HFO2_VASP_HALF_VS_SAD_KSPACING035.md`](validation/HFO2_VASP_HALF_VS_SAD_KSPACING035.md).
With all `LMAX*` tags omitted, VASP reported automatic-mode `LMAXPAW=-100`
and default `LMAXMIX=2`; HALF and SAD converged in 5 and 16 iterations. The
same report includes this ablation and the cross-run timing caveat.

## Remaining validation

1. Compare the mapped coefficient block with an offline HALF result after
   undoing arbitrary band phases.
2. Verify `C^H S C = I` after VASP's `ORTHCH`.
3. Run the same DeepAW `CHGCAR` with random `WFINIT`, WAVECAR restart, and HALF
   initialization. Compare first-step residuals and converged energies.
4. The no-SAD smooth-density check has passed for `LHALF_INIT + ICHARG=1`.
5. `NCORE=1` and `NCORE=2` have passed; `KPAR>1` remains a separate extension.

The present file-based constructor still reparses `CHGCAR` and `POTCAR`. A
future zero-copy integration should add an in-memory constructor for VASP's
already parsed lattice, density grid, and PAW data; that is an optimization,
not a correctness prerequisite.
