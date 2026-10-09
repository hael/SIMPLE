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

# Second run: species_element_model, phases 4 to 6 (2026-10-08 to 2026-10-09)

Carried out unattended on dell under delegation, starting from the
maintainer's checkout after the first run (base revision `d3480ae79`). The
note's section 16 (Progress) holds the detail and every number. Paths below
are relative to the run directory
`/home/elmlundho/agent_runs/species_element_model`. Six maintainer rulings
were issued during the run (`spec/rulings.md`), each recorded in the note
where it belongs. Phases 4 and 5 ended with their exit criteria met. Phase 6
was added by a ruling and then withdrawn by another; its source changes were
undone.

## What changed, phase by phase

**Phase 4, B factors in the element kernels, the species list, the compound
selectors** (`review/phase_4.diff`).
- `atoms%convolve` takes an optional `bfac_pdb`: an element atom's B column
  is then added to the width of every one of its five Gaussians, beside the
  resolution blur. Every existing call is unchanged.
  `simulate_nanoparticle pdbfile= pdb_bfac=yes` exposes it.
- `element=Pt,Ni` on `detect_atoms`: a comma-separated list of the species,
  brightest first, kept in the order given. It implies
  `discover_species=yes`, fixes the number of classes to its length and
  names the classes in `_species.pdb`. Level 1 runs as with the first symbol
  alone.
  - The parser `parse_element_list` is in `simple_atoms`, splitting with
    the new `list_of_strs2arr`.
  - `parameters%element` is now `STDLEN` long, with `species(:)` and
    `l_species_list`. The command-line rules of the note's truth table are
    in `exec_detect_atoms`, and every other program rejects a list.
  - The class mixture, the posterior columns `POST1` to `POSTK` and the
    chains `A` to `Z` then `0` to `9` are sized by the number of classes.
- The compound selectors `CdSeW`, `CdSeZ`, `CdSeR` are no longer cut to `Cd`
  before the lattice lookups. The nanoparticle keeps the selector for the
  lattice and the symbol for radii and names. Under a ruling, the lattice
  analysis of `atoms_stats` was repaired for the binary crystals:
  - one first-shell cutoff function;
  - a bond-length lattice fit;
  - strain skipped for zincblende and wurtzite.
  `single_atoms_stats` gained entries for the three selectors, which assert
  the interior coordination and the fitted lattice constant (budget 35 to
  38).
- `files_identical` (from the last commit) returned false for a file
  compared with itself, which failed the fast gate at the start of the run.
  It was fixed.

**Phase 5, the species test on element ground truth** (`review/phase_5.diff`).
- `species_discovery` is rebuilt on particles rendered with the
  five-Gaussian element kernels, made to look like a real Pt map supplied
  by the maintainer:
  - per-atom B rising to the surface;
  - the measured signal transfer;
  - noise shaped by the measured background spectrum, at 17 noise standard
    deviations for a core atom;
  - half maps for one case.
- Cases: pure Pt, a random quarter Ni (alloy), a Ni core under a Pt skin
  (core), a random quarter Al (light); 16 `detect_atoms` runs judged against
  the generating models through `test_gate`.
- Discovery now prunes recovered atoms by the present contact-score policy
  (a recovered atom faces the level-1 zone and threshold). `_species.txt`
  flags a class confined to the surface layer. Discovery stops with a clear
  message when level 1 finds fewer than two atoms.

**Phase 6, fitting atoms with the map-filtered kernel, withdrawn**
(`review/phase_6.diff`: the note and this report only).
- The approach: a per-shell transfer estimated against the level-1 atoms, a
  radial kernel table, and every Gaussian atom of the discovery branch
  replaced by the kernel.
- On the fixtures the ruled correlation cut zeroed the transfer from 1.68 A
  on. The alloy's intensity ratio came out 52% low with the cut and 37%
  without (phase 5: 25%). On the real map the class spread grew from 21% to
  77% (42% without the cut).
- The sixth ruling withdrew the phase. The transfer estimate and kernel
  table add parameters and failure modes to correct a secondary output.
- The source was restored to the end of phase 5, and the note records what
  was tried and why it failed (section 10.7).

## The rulings and why

1. **Lattice analysis of the compound crystals repaired in phase 4.** The
   selector repair had exposed `fit_lattice` and the coordination numbers of
   `atoms_stats` to the wurtzite and zincblende branches. Those had only ever
   run with the fcc fallback, and the coordination numbers came out at 24 to
   76 instead of 1 to 4.
2. **Fixtures made like a real Pt map; level 1 not changed.** On the first
   fixtures, white noise on a positive pedestal made the present level-1
   threshold search land on a plateau at some noise seeds (3 of 285 atoms).
   The maintainer measured a real map and had the fixtures take its
   properties. The ground-truth intensity became the isolated atom's
   aperture intensity.
3. **Floors, the ratio reported, recovered surface atoms pruned as level 1
   prunes.**
   - The present path's pruning removes surface Pt atoms next to unseen
     light atoms, so Pt recall is floored only with discovery.
   - The intensity ratio is biased low on filtered kernels, so it is
     reported and the class order is floored.
   - Discovery had reinstated an outer shell of 240 weak surface atoms on
     the real map, which `detect_atoms` removes on purpose.
4. **The shell case replaced by a core case.** A light species confined to
   the surface layer cannot be told from partial occupancy, and the pruning
   removes it by design.
5. **The core case's label floor moved to phase 6.** One Ni atom at the
   centre was called Pt, from the radial bias of the aperture intensity.
6. **Phase 6 withdrawn.** It is described above.

## Decisions taken under delegation

All are recorded in Progress under the phase named.
- Phase 4:
  - A compound selector is only one of `CdSeR`, `CdSeW`, `CdSeZ`.
    Four-character pairs such as `CdSe` stay on the present path.
  - The `nspecies` key keeps its bound of 3 for the blind mode; a list sets
    `nspecies` without a bound.
  - A selector names its species only under `discover_species=yes`.
  - The atom name of a `CdSeW` run becomes `CD` (it was `CDSE`).
  - For the binary crystals, the displacement fits of `atoms_stats` take
    half the bond as their radius.
- Phase 5:
  - The Fourier filters of the fixtures act on the continuous spatial
    frequency.
  - The isolated atoms of the ground truth are rendered at the core B.
  - The eligibility under the pruning policy uses the mean Pt position and
    the fcc cutoff.
  - The test reads every table by its header.
- Phase 6: the review report is this run's section of the first run's
  report.

## Deviations from the plan, and why

- The fixtures of 10.6 as first written (B rising to 10 A^2, white noise) and
  several of its floors were replaced by rulings 2 to 5.
- The intensity-ratio floor is reported, not floored, and the core case's
  label check is reported (one wrong label at the centre). Both follow the
  rulings; phase 6, which was to restore them, was withdrawn.
- The shell case is replaced by the core case (ruling 4).
- Phase 6 was added and withdrawn. No source change of it remains.

## Evidence against the baseline

- The phase 0 reference was regenerated at phase 4 with a build of the
  unchanged checkout and with the new build: all 13 files byte-identical
  (285 atoms, 0.0148 A, correlation 0.99969).
- On a noisy copy, `element=Pt`, the key, `Pt,Ni` and `Pt,Ni,Al,Si` wrote the
  five present products byte-identical.
- `single_atoms_stats` was unchanged at the end of every phase: 285 / 285,
  recall and precision 1.0000, 0.0148 A, 0.9997, 19.527 A.
- The three compound entries pass:
  - `CdSeW`: 147 atoms, interior coordination 4, fitted `a` 4.279 against
    4.2985 A;
  - `CdSeZ`: 154, coordination 4, 6.049 against 6.077;
  - `CdSeR`: 193 of 199, coordination 6, 5.468 against 5.49.
- `species_discovery` passes in about 300 s (Debug), with identical metrics
  in two runs:
  - no false atom;
  - Pt recall 1.000; Ni recall 1.000; Al 0.952, where the `element=Pt` run
    finds none of the Al;
  - two-class fits admissible, classes in the right order;
  - half-map agreement 1.00;
  - B rise 1.28 to 1.42 of the generated rise;
  - present products identical across runs.
  The fitted ratios (reported) are 25% (alloy), 13% (core) and 38% (light)
  low.
- The real map (the maintainer's `Pt_example`, by hand):
  - `element=Pt` gives the supplied 476 atoms (474 within 1 A, 0.016 A);
  - with discovery, 267 atoms are recovered and 249 pruned; one species;
  - `Pt,Ni` gives an inadmissible two-class fit.

## Left open for the maintainer

1. **The intensity bias** of the aperture intensity on filtered element
   kernels: a known limitation (note, section 7), and the one central Ni
   atom of the core case. Phase 6's attempt, and why it failed, is in 10.7.
2. **`single_workflow_wurtzite` fails at the base revision too.** It
   crashes in `refine3D` (`simple_strategy3D_srch.f90:211`, `subspace_inds`
   past its bound), and its runs are not reproducible from one run to the
   next. Recorded, not fixed.
3. **The supplied real-map half maps** (`_even`, `_odd`) correlate 0.49 with
   the map and 0.96 with its `_SIM.mrc`. They are not half reconstructions
   of it.
4. **The fixtures match part of the real map.** They match its
   signal-to-noise and transfer, but not all of its surface broadening
   (surface over core peak 0.91 against 0.77) nor its negative wells between
   atoms.
5. **Promotion of recovered atoms**, the validation beyond the element model
   and the A/B on real Pt/Ni data remain as section 10.8 describes.
6. **Found after the merge, on the maintainer's Mac (2026-10-09).**
   - The B-rise floor of `species_discovery` was withdrawn: on filtered maps
     the fitted rise carries a coordination artefact (note, section 7), and
     the floor failed at 1.60 in the core case. The rise is reported.
   - Run by hand without `SIMPLE_SEED`, the test is not reproducible: every
     commander it runs in process reseeds from `/dev/urandom`. CTest sets
     `SIMPLE_SEED=20260923`, the setting the floors were validated with; by
     hand, export it first.
   - With such an unseeded draw, `element=Pt,Ni` on the pure Pt fixture gave
     an admissible two-class split, the second class a block of weaker
     atoms at one end of the radius (most likely surface atoms). The pure
     case's "inadmissible" floor therefore holds for the validated draw, not
     for every draw. How often an asserted second species is found in a
     single-species particle is for the validation on real maps to measure;
     the surface-confined flag is the check meant to catch it.
