# Review report: single_species_proto (SINGLE atom detection without a species input)

Run of 2026-10-08 on dell, carried out unattended under delegation, against
base revision `92d684c6f`. The plan is
`doc/implementation_notes/planned/species_discovery.md` (the note);
its section 16 (Progress) holds the detail, the logs and every number. This
report is the summary for review. All four phases ended with their exit
criteria met. No maintainer ruling was issued during the run.

Paths below are relative to the run directory
`/home/elmlundho/agent_runs/single_species_proto`.

## What changed, phase by phase

**Phase 0, baseline** (`review/phase_0.diff`, the note only). Debug build of
the base with tests (gfortran 15.2.1), fast gate 15 of 15. `single_atoms_stats`
(the regression test of the present detection: a simulated Pt particle,
`detect_atoms element=Pt`, `atoms_stats`) gave 285 atoms simulated and
detected, recall and precision 1.0000 within 1 A, root-mean-square position
error 0.0148 A and a correlation of 0.9997 between the simulated map from the
detected atoms (`_SIM.mrc`) and the input. The same simulation and detection
were run by hand and kept in `scratch/keep/phase0_ref/` as the reference for
phase 2.

**Phase 1, foundations** (`review/phase_1.diff`).
- Pseudo-atoms `X1`, `X2`, `X3` (reserved atomic numbers 901 to 903, clear of
  the 999 of `CDSE`) in `simple_defs_atoms.f90`. `atoms%convolve` renders them
  as one Gaussian, `q (4 pi / B)^(3/2) exp(-4 pi^2 r^2 / B)`, with `q` from
  the PDB occupancy and `B` from the B column. The element branches are
  unchanged.
- `simulate_nanoparticle pdbfile=` renders pseudo-atom PDBs, and the getter
  `atoms%get_occupancy` was added. The unused helper `egau` was deleted.
- New module `simple_nano_species.f90`:
  - the one-dimensional Gaussian mixture of section 3.6 (`fit_species_mixture`),
    with the Bayesian information criterion (BIC), the admissibility rules
    and classes ordered by decreasing intensity;
  - the class separation `D`, the enclosed Gaussian fraction, and the
    threshold calibration of section 3.2.
- Its tester is the sub-suite `species` of the fast suite `unit_single`; the
  atoms tester gained closed-form pseudo-atom tests.

**Phase 2, level 1 without an element** (`review/phase_2.diff`).
- `element` is optional in `detect_atoms` only.
- Without it the present flow runs with a pseudo-atom template
  (`B_ref` = 13.9 A^2). Its length scales come from the median
  nearest-neighbour distance `d_NN`, measured after the first binarisation:
  atom radius `0.4 d_NN`, neighbour cutoff `1.267 d_NN`, contact-score
  ceiling 12, splitting exclusion `0.21 d_NN`.
- On the phase 0 volume the species-free detection found the same 285 atoms,
  with a root-mean-square position difference of 0.014 A, a `_SIM.mrc`
  correlation of 0.9949 (against 0.9997) and `d_NN` of 2.750 A.
- `element=Pt` reproduced the phase 0 products byte for byte.

**Phase 3, residual recovery and species call** (`review/phase_3.diff`, written
by the driver). New keys on `detect_atoms`: `discover_species` (default `no`),
`nspecies` (default 0, number of classes from the data) and `nn_gate`
(renamed `min_nbrs`, an integer, on the Mac after the run: the number of
found neighbours a recovered atom needs, 0 for none)
(default `yes`). The existing `vol_even` and `vol_odd` serve as the optional
half maps.

With `discover_species=yes`, after the present products are written, the
program runs the pipeline of sections 3.2 to 3.9:
- a stage-1 fit of the level-1 atoms;
- the noise region, and the half-map reference when half maps are given;
- calibrated thresholds;
- stage A, then stage B with the neighbour gate;
- pruning of the recovered atoms on the atom table;
- both fit stages, the aperture intensities and the mixture species call;
- `valid_corr` of the recovered atoms;
- three diagnostic files, `_species.csv`, `_species_radial.csv` and
  `_species.txt`; a fourth, `_species.pdb` (every atom with its class in
  the element column and chain, the class intensity ratio as occupancy and
  the stage-2 B factor), was added on the Mac after the run, with a check in
  the test that it agrees with the table.

The recovered atoms live in their own table, and nothing enters the present
products. A new high-level CTest entry, `species_discovery`, raises
`SIMPLE_CTEST_BUDGET` from 34 to 35.

## Decisions taken under delegation

All are recorded in Progress (section 16) under the phase named.
- Phase 1: the resolution blur that `convolve` adds to every element (`lp`)
  is not applied to pseudo-atoms, whose width is their `B`.
- Phase 1: the test policy's list of fast sub-suites names `species`. In
  phase 3 its CTest budget and table of high-level entries were updated. The
  file was added to the note's file table before each edit.
- Phase 2: with no element, `nanoparticle%new` names the atoms `X1`. The
  routines that only the programs still requiring `element` reach stop with
  a message instead of falling back silently to a default lattice.
- Phase 3 choices:
  - discovery runs after the present output block;
  - the half maps leave the command line in `exec_detect_atoms`, because
    `parameters` refuses `vol1` together with them;
  - `d_NN` is measured on the element path as well when discovery is on;
  - positions are held fixed, as in the emulation;
  - the first stage-1 width prior is centred on `B_ref`;
  - the noise region is built from the same soft mask as the map, with a
    20 000-voxel minimum and 1 000 phantom sites taken at a fixed stride;
  - calibration is done once per map;
  - the gate counts only the atoms of earlier levels;
  - pruning applies the contact-score rule of `discard_atoms` in one pass,
    to recovered atoms only;
  - the predicted detection signal-to-noise scales as `I sigma^(-3/2)`,
    normalised to `s_A` at the template width.

## Deviations from the plan, and why

- Discovery runs after the present products are written, not before as the
  code map of section 9 lists it. This makes the identity of the products
  with and without the key structural.
- Positions are not refined each sweep (section 3.4 suggests it). This is
  as the emulation did, and the floors were met without it.
- Phantom sites are on a fixed stride, not random (section 3.2), so that
  `detect_atoms` stays deterministic.
- The test ran `detect_atoms` on one thread during the run, to work around
  a race in the present code (below), not a defect of the new code. The
  race was fixed on the Mac after the run (see below) and the test now runs
  `detect_atoms` with the test's thread count.
- The first code version defined the predicted signal-to-noise as amplitude
  over `s_A` (`I sigma^(-3)`). The first test run exposed it, and the code
  and test were corrected to section 7's `I sigma^(-3/2)` before any further
  run. The floors were not changed.

## Evidence against the baseline

- `single_atoms_stats` gave the same atom counts, recall, precision,
  root-mean-square error and correlation as phase 0 at the end of every
  phase (285 / 285, 1.0000, 1.0000, 0.0148 A, 0.9997).
- The fast gate passed at the end of every phase. The phase 0 Pt reference
  is reproduced byte for byte with `element=Pt` at phase 2.
- `species_discovery` passed, in 505 s Debug. Every floor of section
  10.4 was set before the first run and met; each of the following held in
  both runs, with the map alone and with half maps:
  - for case 2 and R6, light-class recall where the predicted
    signal-to-noise is 6.5 or more was 1.000, 1.000, 0.983 and 0.955;
  - no false atom anywhere;
  - `K` = 2 with every label right;
  - the intensity ratio was within 0.9% of 0.17;
  - per-shell widths were within 0.5% for the strong class and 3.4% for the
    light class;
  - for case 9 and R3, `K` = 1 with no recovered atom;
  - the five present products were identical with and without the key.

## Left open for the maintainer

1. **A data race in the present detection.** `nanoparticle%find_centers`
   (the same loop at the base revision) accumulates per-component centre
   sums from an OpenMP loop without a reduction.
   - At 24 threads one atom centre moved by 0.29 A between two runs of the
     same command, and results at 8 and 24 threads differ by up to 0.034 A.
   - The fix is a `reduction` or a serial loop. It changes an existing
     routine for existing callers, so it was not made by the run.
   - `single_atoms_stats` (8 threads) never showed it during the run.
   - Fixed after the run, on the Mac (2026-10-08): the loop now carries
     `reduction(+:m,sum_mass)`, the house pattern for per-element sums in
     OpenMP loops, and the map value is read once per voxel. The test's
     one-thread workaround was removed at the same time, so its product
     comparison between the three runs also checks that detection is
     thread-invariant. Both need a build and a run of `single_atoms_stats`
     and `species_discovery` to confirm.
   - Also after the run, on the Mac: the species mixture is no longer fitted
     by a private expectation-maximisation in `simple_nano_species` but by
     the clustering library's extreme deconvolution (`simple_xd_gmm`, one
     dimension, measurement noise `s_I` as the per-point noise variance),
     which is the same likelihood family as the variance-floored fit and
     returns the intrinsic class spread directly. The two deterministic
     starts are kept (`xd_gmm%init` gained an optional `means0`); BIC,
     admissibility and the interface of `fit_species_mixture` are unchanged,
     so the `species` sub-suite and `species_discovery` are the checks.
2. **Recovered atoms are diagnostics only.** Their promotion into the
   products, the independent forward-model validation and the Pt/Ni A/B test
   at matched false-positive rate are section 10.5.
3. **Expected false count at stage B.** The calibrated stage-B count before
   the gate is about 1.1 per particle; the gate's further reduction is not
   modelled. The plan's target of 0.2 applies to stage A (Q7).
4. **Existing records not changed.**
   - The `CdSeW` truncation noted at the end of section 9 is recorded, not
     fixed.
   - `atoms%extract_atom` does not copy the occupancy, so a pseudo-atom
     extracted with it renders with zero intensity.
5. **The open questions Q1 to Q7 of section 13 remain open.** The prototype
   took the recommended defaults: radius from the centre, the close-packed
   cutoff, both widths reported, and the splitting exclusion at `0.21 d_NN`.
