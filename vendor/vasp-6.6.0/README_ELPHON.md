# VASP electron-phonon driver

This document contains information regarding the electron-phonon driver implementation in VASP.

# Changelog

## velph-1.4.1 (1/12/2023)
- Added ELPH_MODE metatag with current option "transport". This changes the default
settings for a few other INCAR tags.
- Added test for non-collinear electron-phonon calculations

## velph-1.4.0 (17/11/2023)
- Added support for NCORE/=1 paralelisation
- Improved scalling of the openMP paralelisation

## velph-1.3.0 (11/11/2023)
- Implement band paralelism in the electron-phonon driver
- Average FAN and DW terms over degenerate states
- Added experimental support for spin-polarized and noncollinear spin calculations
- Introduce k-point path selector to compute the self-energy in points from
a regular grid that fall along a path specified by KPOINTS_LINE.
- Introduce transport selector to select for which KS states the lifetimes need
to be computed to perform phonon-limited transport calculations.
- Compute of the transport quantities using simple linear grids and
Simpson integration (ELPH_TRANSPORT_DRIVER=1) or
Gauss-Legendre integration (ELPH_TRANSPORT_DRIVER=1).
- Introduce the elph_ismear tag to choose a different smearing method to
compute the chemical potential in the dense k-point grid
- Refactoring of the output of the electron-phonon driver to the hdf5 file
- Fix bug in the rotation of the WF spinors
- Numerous other bug fixes and improvements

## velph-1.2.1 (10/03/2022)
- Fix computation of the velocity matrix elements
- The computation of the ddijdu was wrongly deactivated

## velph-1.2.0 (17/01/2022)

- Added support for running electron-phonon calculation with different supercell size
for phonons and electron-phonon potential
- Fix in the computation of ewald parameter from encutlr
- Multiple bugfixes, code cleanup and documentation of some routines by Laurent Chaput
- Compute long-range part of the potential in the plane-wave grid using FFT
- Fixed important bug when using kspclib to generate kpoint mesh in
electron-phonon driver
- Use vasp linear optics routines from VASP to compute band-velocities
in electron-phonon driver

## velph-1.1.1 (12/10/2021)

- Added elph_selfen_kpts incar tag and and support for multiple entries in
elph_selfen_ikpt. This allows to compute the self-energy on a list of k-points.

## velph-1.1.0 (23/09/2021)


- Added possibility to specify an array of smearings. Each smearing will
create a different accumulator and the results are written to hdf5 file.
- Implemented tetrahedron method for evaluating the imaginary part of the
electron-phonon self-energy
- Added support for generating k-point meshes and finding irreducible triplets
using kspclib
- Added caching algorithm for the WFs to avoid repeated communication
- Improve memory and timing reports
- Compute fermi level using DENSTA
- Avoid computation of electron-phonon matrix elements that are not required
and add ELPH_NBANDS to choose number of bands treated in the electron-phonon driver
- Treatment of the long-range part of the Dijs
- Refactoring of the main electron-phonon driver (ELECTRON_PHONON routine) by
encapsulating parts into separate routines)
- Added consistency checks for atom positions in phelel and VASP


## velph-1.0.0 (06/05/2021)


- Read `phelel_params.hdf5` produced by phelel. This contains the
  derivative of the potential for all the displacements.
- Use `wave_interpolator` to compute distributed WFs for batches of
  k-points.
- Use `wave_window` to expose the WFs to all the MPI ranks using MPI
  one-sided comunication.
- Use `wave_rotator` to rotate the WFs from the IBZ to the FBZ.
- Backfolding of the electron-phonon potential to the primitive cell
  by means of Fourier interpolation.
- Build the dynamical matrices at arbitrary q-points and diagonalize
  to get phonons eigenvectors and eigenvalues.
- Special treatment of the long-range interaction for the potential
  (dipole and quadrupole) and dynamical matrix (dipole-dipole).
- Computation of the electron-phonon matrix elements in VASP (only
  tested for spin unpolarized calculations).
- Computation of the electron self-energy due to electron-phonon
  interaction using Lorentzian smearing.
- Computation of transport quantities (conductivity, electronic
  thermal conductivity, seebeck, peltier).

## Description of the output of the electron-phonon driver in VASP

A) The lowest order ('Fan') selfenergy is computed (see
   <https://www.cond-mat.de/events/correl17/manuscripts/heid.pdf>)

   $$
   \Sigma^{(r)}_{\mathbf{k}n}(\hbar \omega, T, \mu) = \frac{1}{N}
   \sum_{\mathbf{k}'} \sum_{j=1}^{3 N_a} \sum_{n'=1}^{N_b}
   |g(\mathbf{k}n,\mathbf{k}'n', \mathbf{k}-\mathbf{k}'j)|^{2}
   \Big[
   \frac{1-f^0_{\mathbf{k}'n'}+n^0_{\mathbf{k}-\mathbf{k}'j}}{\hbar
   \omega -\epsilon_{\mathbf{k}'n'}-\hbar
   \omega_{\mathbf{k}-\mathbf{k}'j}+i\eta}+
   \frac{f^0_{\mathbf{k}'n'}+n^0_{\mathbf{k}-\mathbf{k}'j}}{\hbar
   \omega -\epsilon_{\mathbf{k}'n'}+\hbar
   \omega_{\mathbf{k}-\mathbf{k}'j}+i\eta} \Big]
   $$

   with

   $$
   g(\mathbf{k}n,\mathbf{k}'n', \mathbf{q}j)
   =  -i   \frac{\epsilon_{\mathbf{k}n}-\epsilon_{\mathbf{k}'n'}}{
   \hbar }   \sum_{ \tau}  \langle \varphi_{\mathbf{k}n} | - i \hbar
   \frac{\partial}{\partial \mathbf{R}_{\tau}} |
   \varphi_{\mathbf{k}'n'} \rangle
   \cdot \Big( \sqrt{\frac{\hbar }{2M_{\tau} \omega_{\mathbf{q}j}}}
   \mathbf{e}^{\mathbf{q}j}_{\tau}  \Big)
   \Delta(\mathbf{q}+\mathbf{k}'-\mathbf{k}).
   $$

   In this formula, time reversal symmetry is assumed.

   The selfenergy can be found in vaspout.h5, in `/results/self_energy_elph_*`.
   The electron-phonon matrix elements can be found in the file
   `vaspelph.h5`, which is writen if `ELPH_WRITEMELS =.TRUE.` in the
   `INCAR.01`

B) The code internally creates multiple instances of the selfenergies
   based on the values of the following variables:

   - `elph_selfen_carrier_den` OR `elph_selfen_carrier_per_cell` OR
     `elph_selfen_mu` which the doping through the electron
     concentration or the position of the chemical potential.

   - `elph_nbands_sum` which control the upper bound in the summation
     over $n'$.

   For example, in the `INCAR` if we give
   `elph_selfen_carrier_per_cell = 0 1e-12 1e-11 1e-10`
   `elph_nbands_sum = 10 20`,
   then 6 selfenergies with every possible combinations of the values
   will be created.

   In `vaspout.h5` you will find one entry for each of the
   corresponding selfenergies `/results/self_energy_elph` for the
   first one `/results/self_energy_elph_*` for all the other ones.

   Inside each of these entries you will find the values that were
   used to compute the considered self-energy:

   - `/results/self_energy_elph/efermi` contains the chemical
     potential that was computed based on one of
     `elph_selfen_carrier_den`, `elph_selfen_carrier_per_cell` or
     `elph_selfen_mu`.
   - `/results/self_energy_elph/nbands_sum` contains the number of
     bands that were used in the sum based on `elph_nbands_sum`.

   The fan self-energy is written in
   `/results/self_energy_elph/selfen_fan`.

   When reading from the hdf5 file in python you will find the
   dimensions of `/results/self_energy_elph/selfen_fan` (real) to be
   `(nspin,nkpoints,nbands,nw,ntemps,2)`

   - 2 holds the real and imaginary part of the complex number
   - `ntemps`: number of temperatures at which the self-energy is evaluated
   - `nw`: number of frequencies at which the self-energy is evaluated
   - `nbands`: number of bands for which the self-energy is computed
   - `nkpoints`: number of k-points in the irreducible Brillouin zone
   - `nspin`: number of spin channels (1 for non spin polarized or
     noncollinear spin calculations or 2 for spin-polarized
     calculations)

   The frequency grid at which the self-energy is evaluated is built
   around the Kohn-Sham eigenvalue $\epsilon_{\mathbf{k}n}$ using
   the following prescription:

   ```
   energies = [enk + ((iw-1)/nw-0.5)*enwin for iw in range(1,nw+1)]
   ```

   where `enwin` is set, in the `INCAR`, by `elph_selfen_enwin`
   (default: 0.0), and `nw` by `elph_selfen_nw` (default: 1).

## Definition of the inputs in INCAR

   - `elph_run` : to activate the electron-phonon calculation (logical).
   - `elph_selfen_carrier_den` : carrier densities at which to compute
     the self-energy in $\text{cm}^{-3}$ (real array, default=0).
   - `elph_selfen_carrier_per_cell` : doping  electron concentration per
     cell at which to compute the self-energy (real array, default=0).
   - `elph_selfen_mu` : chemical potential (eV) wrt Ef (T=0, no
     doping) (real array, default=0).
   - `elph_selfen_temps` : temperatures (K) at which to compute the
     self-energy (real array, default= 0,100,200,300,400,500).
   - `elph_selfen_delta` : imaginary energy shift to compute the
     self-energy (eV ?) (real, default= 0.1).
   - `elph_nbands_sum`: number of bands to sum over = $N_b$ ( integer,
     default: all bands).
   - `elph_selfen_enwin`: energy window at which to compute the
     self-energy (real, default= 0.0).
   - `elph_selfen_nw`: number of frequency points at which to compute
     the self-energy (integer, default=1).
   - `elph_selfen_fan`: compute Fan selfenergy (logical).
   - `elph_selfen_dw`: compute DW selfenergy (logical).

## How to make a computation: the example of diamond (17/05/2021)


STEP 1: prepare the directory (phelel) to compute potential derivatives

```
mkdir phelel
cd phelel
```

STEP 2: build the `POSCAR` that will be used to generate the supercells

```
C
   1.00000000000000
     3.5670000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    3.5670000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    3.5670000000000000
   C
     8
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.5000000000000000  0.5000000000000000 -0.0000000000000000
  0.5000000000000000  0.0000000000000000  0.5000000000000000
  0.0000000000000000  0.5000000000000000  0.5000000000000000
  0.2500000000000000  0.2500000000000000  0.2500000000000000
  0.7500000000000000  0.7500000000000000  0.2500000000000000
  0.7500000000000000  0.2500000000000000  0.7500000000000000
  0.2500000000000000  0.7500000000000000  0.7500000000000000
```


STEP 3: generate displaced supercells and put them into sub-directories
`phelel/perfect`, `phelel/disp-001`, `phelel/disp-002`, ...

```
phonopy -d --dim 2 2 2 -c POSCAR --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"
mkdir perfect disp-001
mv SPOSCAR perfect/POSCAR
mv POSCAR-001 disp-001/POSCAR
```

STEP 4: prepare the sub-directories to run VASP calculations. The
`INCAR` must contains `LWAP = .TRUE.` , and could be

```
PREC = Accurate
ENCUT = 600
IBRION = -1
NSW = 0
NELMIN = 2
EDIFF = 1.0e-08
IALGO = 38
ISMEAR = 0; SIGMA = 0.2
LREAL = .FALSE.
ADDGRID = .TRUE.
LWAVE = .FALSE.
LCHARG = .FALSE.
LWAP = .TRUE.
```

STEP 5: run the VASP calculations `phelel/perfect`, `phelel/disp-001`,
`phelel/disp-002` ...

STEP 6 (to be removed): determine the size of the FFT mesh along the
axes of the primitive cell (as defined in `phelel/phonopy_disp.yaml`)

method : run separate VASP calculation for the primitive cell, and
read `NGX`, `NGY` and `NGZ` in the `OUTCAR`. They give the fft mesh
size.

STEP 7: compute the potential derivative on the previously determined
FFT mesh

```
phelel -c POSCAR --dim 2 2 2 --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0" --create-derivatives perfect disp-001 --fft-mesh 20 20 20
```

STEP 8: prepare the directory to compute the electron-phonon
interactions

```
cd ..
mkdir ep
cd ep
```

create `POSCAR` from `phelel/phonopy_disp.yaml`,

```
C
1.00000000000000
0.000000000000000     1.783500000000000     1.783500000000000
1.783500000000000     0.000000000000000     1.783500000000000
1.783500000000000     1.783500000000000     0.000000000000000
C
2
Direct
0.000000000000000  0.000000000000000  0.000000000000000
0.250000000000000  0.250000000000000  0.250000000000000
```

STEP 9: run VASP with the following `INCAR`

```
PREC = Accurate
ENCUT = 600
IBRION = -1
NSW = 0
NELMIN = 2
EDIFF = 1.0e-08
IALGO = 38
ISMEAR = 0; SIGMA = 0.2
LREAL = .FALSE.
ADDGRID = .TRUE.
LWAVE = .FALSE.
LCHARG = .FALSE.

ELPH_RUN=.TRUE.
ELPH_SELFEN_TEMPS= 0 1 10 100 1000
ELPH_SELFEN_MU= 0 -0.1 0.1
ELPH_SELFEN_FAN = .TRUE.
ELPH_SELFEN_DELTA = 0.15
ELPH_KSPCLIB=.TRUE.
ELPH_KSPACING=6
```

The last 2 tags can be removed if instead you provided a file
`KPOINTS_DENSE` which will determine the k-points to be used the
compute the self enenrgy. `KPOINTS_DENSE` has the same format that the
usual KPOINTS file. Instead `ELPH_KSPCLIB` and `ELPH_KSPACING` are
used to generate generalize k-points grids.

STEP 10: plot the self energy

For example using the python script `plot_lifetimes.py`.
For example,

```
python plot_lifetimes.py --temp=3 --mu=2
```

will plot the absolute value of the imaginary part of the self energy
for the third temperature given in the `INCAR` (10 K) and the second
chemical potential (Ef-0.1 eV).

## Suggestions

1. at STEP 3 call to phelel rather than phonopy could be nice.
2. at STEP 8 if the POSCAR of the primitive cell was already prepared,
   it would be nice.
3. `ELPH_KSPACING` should be renamed.

## questions

1. definition of 'bandstop' in `plot_lifetimes.py`
2. energies = [enk + ((iw)/(nw+1)-0.5)*enwin for iw in range(1,nw+1)
3. for iw in range(1,nw+1) is it with python ?  or is range [1..nw]?
4. meaning of `elph_selfen_enwin` is unclear. May be change input to `delta_w`?
5. `elph_nbands_sum` default?
