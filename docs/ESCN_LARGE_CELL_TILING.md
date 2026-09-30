# DeePAW-eSCN large-cell tiling design

## Current limit

The v1 client sends one complete periodic structure and requests one complete
`[nx, ny, nz]` density grid.  It rejects grids above 2,000,000 points and
requires the returned shape to equal the requested global shape.  Splitting
only the response buffer cannot remove this limit.

Both layers currently require three-dimensional periodicity.  The HALF request
structure has no PBC field and its transport hard-codes `[true,true,true]`.
Direct requests to the running eSCN service confirmed that `[true,true,false]`
and `[false,false,false]` return HTTP 400 with the server error that the current
model supports only three-dimensional periodic crystals.  A slab must therefore
remain a 3D-periodic cell with explicit vacuum.  Turning PBC off is not an
available shortcut for independent spatial tiles.

Large-cell support must preserve periodicity, translational covariance, local
atomic environments, negative density values, and one global electron-count
normalization.  Independent normalization or clipping of tiles is forbidden.

## Preferred implementation: coefficients first, grid second

When the density head is an atom-centred expansion, expose two service stages:

1. `/v1/predict_coefficients` evaluates the atomic graph and returns the
   learned per-atom density coefficients, model cutoff metadata, and model
   identity/hash.
2. `/v1/render_density_tiles` evaluates those fixed coefficients on requested
   fractional-grid tiles.  Each grid point receives every periodic atomic
   contribution within the density basis cutoff.

The client assembles non-overlapping core tiles into the global C-order array
and performs one POTCAR-valence normalization after all tiles are complete.
Because all tiles use the same atomic coefficients and analytic density basis,
the result should be independent of tile shape, order, and decomposition, up to
floating-point reduction order.

This is preferable to running the neural model independently on cropped
subcells: changing a subcell changes its periodic graph and can change atomic
embeddings even when a geometric halo is added.

## Immediately usable method: shifted sparse probes

The probe formulation also permits a client-side decomposition without
cropping the periodic structure.  Let the desired dense grid be
`N=(Nx,Ny,Nz)` and choose integer strides `S=(Sx,Sy,Sz)` that divide the three
grid dimensions.  For every residue `r` in `[0,S)`, request the sparse grid
`N/S` after shifting every atom by

`-[(rx/Nx) a1 + (ry/Ny) a2 + (rz/Nz) a3]`.

Translational covariance gives `rho_(R-delta)(x)=rho_R(x+delta)`.  The sparse
result therefore fills dense indices `r + S*j`.  All requests retain the full
cell, full atom list, and `[true,true,true]` PBC, so there are no spatial block
boundaries or artificial cropped-cell environments.  Assemble every residue
class first and normalize the global density only once.

For Si, eight shifted `30^3` requests reconstructed one `60^3` prediction with
relative density L1 error `3.42e-4`, mean-density difference `1.04e-7`, and
nearly identical aggregate inference time (`9.52 s` sparse versus `9.57 s`
monolithic).  Peak probe-grid size was reduced by eight.  The result was not
bitwise identical, so this path still requires downstream HALF energy/band and
VASP SCF gates.  Different residues can be dispatched to multiple service
workers or GPUs; caching the unchanged atomic graph would remove repeated graph
inference.

This method requires exact stride divisibility.  If no useful divisor exists,
the service needs explicit probe coordinates or coefficient-first rendering;
silently changing the FFT grid is not allowed.

## Fallback: halo-based model tiling

If coefficients cannot be exposed, add `/v1/predict_tile` with the global cell,
global atom list, global grid shape, integer core origin/shape, and a halo in
angstrom.  The service, not the caller, constructs periodic images and returns
only the core grid.

The required halo is a model property.  It must cover the complete graph
message-passing receptive field plus the density-head support radius.  A single
neighbor cutoff is insufficient when multiple message-passing layers expand
the receptive field.  Models with global pooling or global attention cannot be
made exact by a finite halo.

Overlapping predictions may be combined only with a partition-of-unity window.
The overlap result is still an approximation and must be labelled separately
from coefficient-first exact rendering.

## Periodic geometry

- Tiles are indexed in the global fractional FFT grid, including triclinic
  cells; they are not independent Cartesian boxes.
- Halo membership uses the global cell metric and periodic minimum images.
- Core grid indices are owned exactly once.
- Tile assembly preserves the API's C-order `[nx, ny, nz]` convention before
  HALF converts it to its internal order.
- Density, `nu`, `alpha`, `beta`, and `risk` use identical tile ownership and
  assembly rules.

## Proposed API metadata

Each response records:

- global grid and tile origin/core/halo shapes;
- model hash and density-basis version;
- graph cutoff, message-passing depth, density support radius, and effective
  halo;
- coefficient-first or approximate-halo mode;
- raw global electron integral and final normalization scale;
- per-tile elapsed time, peak memory, retries, and checksum.

Restart manifests mark a tile complete only after shape, model hash, geometry
hash, and payload checksum validation.

## Validation gates

Use existing cells small enough for monolithic inference as numerical oracles.
For every validation cell compare monolithic inference against several tile
shapes, tile orders, and translated tile origins.

Required checks are:

1. one-tile output is byte- or tolerance-equivalent to the existing endpoint;
2. decomposition invariance of the assembled density and uncertainty fields;
3. no elevated error in boundary slabs relative to tile interiors;
4. identical global electron count before and after assembly, followed by one
   recorded valence normalization;
5. low-G Fourier components, HALF energy, sampled bands, and VASP SCF outcome
   remain within the matched numerical tolerances;
6. memory scales with tile volume rather than global grid volume;
7. failed tiles can be retried without recomputing accepted tiles.

The 360-structure density set supplies the first regression suite.  The next
gate is a replicated bulk/surface system whose monolithic result still fits in
memory.  Only after both gates pass should the tiled path be used for the large
NRR interface.

## Scope boundary

Tiling removes the DeePAW inference and density-grid memory bottleneck.  It does
not make the subsequent plane-wave Hamiltonian or VASP diagonalization linear
scaling.  The large-system claim must therefore report DeePAW inference,
HALF initialization, and VASP SCF costs separately.
