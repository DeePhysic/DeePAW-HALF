# Spin-resolved NRR application plan

## Scientific claim

DeePAW-HALF is used as a universal, spin-channel-aware initializer that makes
large ensembles of fixed-moment electrocatalytic structures affordable.  The
application target is not another single 0 K free-energy diagram.  It is a
spin- and environment-resolved account of NRR/HER competition on M-N4 carbon,
with final chemical energetics verified by matched VASP calculations.

The present `spin-channel` implementation uses a shared total-density Harris
operator.  It supplies controlled NUPDOWN occupations and two-channel VASP
wavefunctions; it is not a replacement for self-consistent spin-polarized VASP
energies.  Mechanistic claims therefore use converged VASP energies and forces.

## Pilot: Fe-N4/graphene

Use one relaxed Fe-N4/graphene slab with at least 15 A vacuum and four adsorbate
states:

1. clean Fe-N4;
2. end-on *N2;
3. *NNH;
4. *H, as the competing HER state.

For each state run `NUPDOWN = 0, 2, 4, 6` when allowed by the electron count.
Every calculation must use the same POSCAR, POTCAR, ENCUT, k mesh, smearing,
dipole correction, and solvation settings in its paired comparison.

For each fixed moment:

1. predict and valence-normalize the DeePAW density;
2. run `half energy --ispin 2 --nupdown M` and export `vaspwave.h5`;
3. start VASP from the exported two-spin wavefunctions with `ISPIN=2` and the
   same `NUPDOWN=M`;
4. run a cold-start VASP control from identical inputs;
5. compare the converged state, energy, local moments, SCF iterations, wall
   time, and failure/restart rate.

## Pilot acceptance gates

- Electron-count and NUPDOWN errors below 1e-8 electrons in every exported HDF5.
- Spin-reversal symmetry checked for isolated numerical controls.
- DeePAW-HALF-started and cold-start VASP converge to the same fixed-moment
  state within 1 meV/atom and 0.02 mu_B on the Fe local moment.
- At least 2x median VASP electronic-iteration or wall-time speedup without an
  increased failure rate.
- The lowest spin state of *N2, *NNH, and *H is reproduced in three independent
  reruns and remains stable after ionic relaxation.

Failure of the state-matching gate is a method diagnostic, not chemistry data.

## Publication-scale extension

After the pilot passes, extend to Fe-, Mo-, and Ru-N4 sites and the alternating
and distal NRR networks:

`*N2 -> *NNH -> *NHNH/*NNH2 -> *NHNH2/*NH2NH2 -> *NH2 -> *NH3`.

For each chemically distinct intermediate, retain all fixed-moment states
within 0.30 eV of the minimum.  Include *H and H2 formation as explicit
competitors.  Compute transition states for N2 adsorption/activation, first
protonation, N-N cleavage where relevant, and the rate-controlling HER step.

Use explicit-water configurations rather than one optimized solvent geometry:

- generate equilibrated interfacial trajectories;
- cluster configurations by N2 orientation, hydrogen-bond topology, and local
  electric field;
- evaluate at least 50 independent snapshots per decisive intermediate;
- refine representative minima and transition states with spin-polarized VASP;
- report distributions and confidence intervals, not only minimum energies.

Potential dependence must be treated with a constant-potential or validated
grand-canonical correction for decisive states.  Feed potential-dependent NRR
and HER barriers into a microkinetic model to predict coverage crossover,
selectivity, and rate as functions of potential and pH.

## Required controls

- VASP cold-start versus DeePAW-HALF-start at identical NUPDOWN.
- Total-density normalization on/off diagnostic, with the normalized path used
  for production.
- At least one higher-level electronic-structure check for spin ordering and
  the first protonation barrier.
- Cell-size, k-point, slab-thickness, vacuum, solvation, and dipole convergence.
- Multiple explicit-solvent seeds and transparent exclusion criteria.
- HER competition and clean-surface controls.

## Go/no-go decision

Proceed to the publication-scale ensemble only if the four-state Fe-N4 pilot
passes every acceptance gate.  The Nature Communications-level contribution is
the combination of universal 87-element validation, verified acceleration, and
a new spin/environment-dependent NRR mechanism.  A single optimized-path free
energy diagram is not sufficient.
