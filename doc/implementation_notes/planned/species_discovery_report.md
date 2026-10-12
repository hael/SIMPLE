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
   - Both are settled since. Phase 4 of the second run repaired the
     compound selectors, and on 2026-10-12 `extract_atom` was made to copy
     every per-atom field (note, Progress, before phase 9).
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
   The central atom is labelled right since the common width of item 7.
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
   - Run by hand, the test was not reproducible: every commander it runs in
     process reseeds from `SIMPLE_SEED`, which only CTest set (20260923, the
     setting the floors were validated with), and from `/dev/urandom`
     otherwise. Since 2026-10-11 `simple_test_exec` sets that value for its
     own process, so a run by hand of this or any other test draws what a
     CTest run draws.
   - With such an unseeded draw, `element=Pt,Ni` on the pure Pt fixture gave
     an admissible two-class split. The pure case's "inadmissible" floor
     held for the validated draw, not for every draw.
7. **The cause of that split, and the fix (2026-10-11).** A numpy emulation
   of the pipeline on the pure fixture, over 35 noise draws, reproduced the
   split in 9 of them.
   - It is not the surface splitting off as a weaker class, as first
     guessed. The brighter class was the low-coordination surface atoms, and
     the class called Ni was the bulk of the particle. The Mac run's numbers
     agree: the Pt class's inner shell moved out (B 28.7 against about 25
     for all atoms in the other runs) while its outer shell stayed.
   - The cause was the aperture intensity's division by the enclosed
     fraction at each atom's own fitted width. On filtered maps that width
     is narrowed by the neighbours' halos, so the correction lowered the
     intensity of well-coordinated atoms. Noise-free, the intensity fell by
     19% from 5 to 12 neighbours, enough for a two-class fit at `D` of about
     3, the admissibility threshold. The stage-2 tie fed the class back into
     the width and pushed 6 of the 9 draws over the threshold.
   - The fix divides every aperture sum by the enclosed fraction at one
     common width, the median fitted B. In the emulation the forced split
     was then admissible in 0 of 35 draws (`D` at most 2.54). On the alloy,
     core and light fixtures every atom was labelled right in 7 draws of
     each, including the Ni atom at the centre of the core case that phase
     5 had called Pt; `D` rose by 1 to 3.
   - Note sections 3.5, 3.7 and 7 record it. The pure case's report now
     shows the forced classes' `D`, the floor's margin. The maintainer's
     build and test run passed, and the fix was committed as `0029ed77e`.
     On 2026-10-12 a build that also carried the coordination term labelled
     every atom of the core case right, the centre Ni atom included, and
     that case's label floor was restored.

# Third run: species_validation_r2, phases 7 and 8 (2026-10-11)

Carried out unattended on dell under delegation, on the maintainer's
checkout at `c28bfc101`. That revision carries the fix of the automatic mask
of `detect_atoms` and takes discovery's noise region from the map as read
(note, sections 7 and 3.2). It is a measurement run: no source file changed.
The note's section 16 (Progress) holds the detail. Paths below are relative
to the run directory `/home/elmlundho/agent_runs/species_validation_r2`. The
tables, the scripts that made them and every `detect_atoms` log are in
`scratch/keep/validation/`; this section does not reproduce them. No ruling
was issued. Both phases ended with their exit criteria met. Of 10.9's six
expectations two held (one class on 95% of the maps, the control never a
class) and four did not (every map reproduced, `Pt,Ni` inadmissible on 95%,
the planted class found, the diameter within 15% everywhere), as reported
below.

A first attempt at this validation (run `species_validation`, the same day)
found that the mask estimate covered the whole box on these unmasked maps.
The present detection then reproduced only 33 of the 66 maps, and on 11 it
returned background. Its numbers are compared below wherever the fix
changed them.

Terms. `detect_atoms` first runs the present detection, called level 1:
template correlation, a threshold, connected components, splitting, and
pruning by contact score. It writes five products, the atoms
(`_ATMS.pdb`) and four maps. With `discover_species=yes` it then adds the
recovered atoms: weaker atoms found in the residual map and pruned by the
same contact-score policy. It measures each atom's aperture intensity (the
density summed in a sphere of half the nearest-neighbour distance) and fits
one to three intensity classes. A two-class fit is admissible when both
classes are large enough and separated by `D` of 3 or more (`D` is the
distance of the class means over their pooled spread). The five runs per map
are run 1 `element=Pt` (the present detection), run 2 `element=Pt
discover_species=yes`, run 3 `element=Pt,Ni` (two classes asserted), run 4
no element (level 1 with a Gaussian template), and run 5 no element with
`discover_species=yes`.

## What changed, phase by phase

**Phase 7, five runs of `detect_atoms` on the 66 maps**
(`review/phase_7.diff`: the note's Progress only).
- Release build with tests, the fast gate (the quick unit tests every build
  runs) and the five CTest entries, all passing before any map. The fast
  gate included the two new mask checks.
- The 66 maps of `/data/NanoX/zenodo` with their stored coordinates, five
  runs each with the program's own mask, 330 runs.
- Per-map tables `validation.tsv` and `environment.tsv`.

**Phase 8, planted classes and the comparison across maps**
(`review/phase_8.diff`: the note's Progress and this section).
- 20% of the level-1 atoms inside the pruning zone of each map were dimmed
  to a fraction `r` of their density, with the map's own background noise
  put back. The levels were `r` = 0.6, 0.4 and 0.25, plus a control that
  adds the noise alone. That is 264 planted maps and 3792 planted atoms.
- Runs 2, 3 and 5 on each planted map, 792 runs, in `planted.tsv`.

## Decisions taken under delegation

All are recorded in Progress under the phase named.
- Phase 7:
  - The first attempt's scripts were reused, with two log lines added: the
    measured particle diameter and its ratio to the coordinate extent. The
    two runs therefore measure the same way.
  - The kept recovered atoms were also judged against the stored
    coordinates, which the first attempt did not do.
- Phase 8: the planting is the first attempt's `plant.py`, unchanged but for
  paths. That covers four choices:
  - the zone centre is the mean of the level-1 atoms;
  - the random seed is `crc32` of the map name;
  - the background patch comes from beyond the outermost atom plus two
    nearest-neighbour distances;
  - the header is copied byte for byte.
  The planted atoms are therefore the same wherever the two runs' level-1
  atoms agree.

## Deviations from the plan, and why

None in what was run. 10.9 gives the stored atom count's median as 315; it is
313 (the two middle values are 311 and 315), as the first attempt also
found.

## Evidence against the baseline

- `single_atoms_stats` shows exactly the numbers of phase 0 of the first run:
  - 285 of 285 atoms;
  - recall and precision 1.0000;
  - 0.0148 A root-mean-square error;
  - correlation 0.9997.
  Its three CdSe entries show the numbers of phase 4.
- `species_discovery` passes with no failed check. Every metric that moved
  against the first attempt is named in Progress. They moved by at most a
  few percent, from the new mask on the simulated box and the new noise
  region, for example:
  - the alloy's intensity ratio went from 0.4350 to 0.4324;
  - the present detection found one more atom on the core particle.
- The five present products are byte-identical between runs 1, 2 and 3 and
  between runs 4 and 5 on all 66 maps.

## Across maps, against the expectations of 10.9

| Item | Expected | This run | First attempt |
| --- | --- | --- | --- |
| present path, maps reproduced (98% of stored atoms within 1 A) | every map | 56 of 66 (run 1), 47 (run 4) | 33, 44 |
| one class found | 95% of maps | 95.5% (run 2), 97.0% (run 5) | 93.9%, 93.9% |
| `Pt,Ni` two-class fit inadmissible | 95% | 93.9% (62 of 66) | 92.4% |
| planted class found at `r` 0.6 and 0.4 (recall 0.9, at most one false label), best of runs 2 and 5 | 90% of maps | 0 and 1 of 66 | 0 and 2 |
| control reported as a second class | never | never | never |
| measured particle diameter within 15% of the coordinate extent | every map | 64 of 66 (ratio 0.91 to 1.27) | none (its mask was 2.9 to 4.0 times the extent) |

- **Present path.** The fix brought every map's particle back: no map
  matches below 84%, where 11 were below 50%. The ten maps short of 98% in
  run 1 are ESC1_time28s 0.849, FCHCS1_time23s 0.923 (0.996 in the first
  attempt), ESC4_time04s 0.947, EFC2_time39s 0.956, FC1_time21s 0.960,
  EFC1_time29s 0.961, ESIP2_time27s 0.966, ECT2_time15s 0.973, ECT2_time22s
  0.977 and ESIP1_time03s 0.979. The remaining misses are
  surface atoms in seven of the ten maps that fall short, and interior atoms
  on three of them (ESIP2_time27s, FCHCS1_time23s, EFC1_time29s). The mask
  cuts no stored atom on any of them. ESC1_time28s is the worst (0.849), and
  it is sensitive to the mask:
  - run 4 matches 0.989 on the same map;
  - a mask of 37.8 A instead of the automatic 40.1 A matched 568 of 569 in
    the first attempt.
  Eight maps the first attempt reproduced without an element now fall just
  short (0.960 to 0.979).
- **False positives.** The three maps with a second class are:
  - EFC4_time03s and ESC4_time01s: 22 and 11 atoms at about 0.4 to 0.5 of
    the Pt intensity, every one in the outer zone, so the surface-confined
    flag fires;
  - ESIP2_time05s: 12 atoms brighter than the rest.
  The forced `Pt,Ni` separation `D` runs from 0.01 to 7.72 (median 1.92).
  The surface-confined flag fires on 2, 11 and 1 maps in runs 2, 3 and 5.
  Kept recovered atoms number 0 to 117 per map (median 12.5) in run 2.
  Of the 1534 recovered atoms kept in run 2, 187 sit on stored atoms that
  level 1 missed. The other 1347 sit one bond length (median 2.57 A) from
  the nearest stored atom, at sites the stored model does not hold. Their
  detection z is 4.8 to 9.5 for 80% of them (z is the residual peak in noise
  standard deviations). They never form a class of their own.
- **Environment.** These numbers changed by at most 0.1 A^2 or a few
  thousandths with the fix:
  - the aperture intensity rises with coordination, by +3.8% per neighbour
    at the median map, and is 26% higher at 12 neighbours than at 5;
  - it falls from 1.14 to 0.87 of the map mean between the inner and the
    outer fifth of the atoms;
  - the fitted width is flat: B 28.9 A^2 at 6 neighbours, 28.3 at 12, and
    flat with radius (median slope -0.45 A^2 against `(r / r_max)^2`).
  On the simulated particles of the note's section 7, the intensity fell
  with coordination and the width fell by 4 A^2. The real maps show neither
  artefact, and their intensity trend has the opposite sign. Partial
  occupancy and disorder of the surface are a likely reason, not tested
  here; the B-factor work needs to know which it is.
- **Time series.** There is one class in runs 2 and 5 at every time point
  of 12 of the 15 particles. The exceptions are single time points (EFC4 at
  3 s in run 2 only, ESC4 at 1 s, ESIP2 at 5 s). The recovered counts do not
  stay the same within a particle (EFC2 1 to 56, ESC1 24 to 111).
- **Sensitivity.**
  - At `r` = 0.6 the planted atoms are detected (99% at the median map) but
    form no class. Their measured intensity, 0.73 of the rest, lies inside
    the spread of real Pt intensities, which is about twice the measurement
    noise.
  - At 0.4 and 0.25 a second class appears on about a third of the maps
    (21 and 25 of 66), clustered by particle: every time point of ESC4 and
    FCHCS1, and five of EFC4's six. There its recall is 0.83 and 0.72,
    its precision about 0.8, and it carries about ten false labels.
  - The fitted class ratio is above `r` by +0.13 at 0.4 and +0.20 at 0.25
    (medians; range +0.06 to +0.33), and equals the measured intensity of
    the planted atoms on those maps (0.53 against 0.54 at 0.4, 0.45 against
    0.45 at 0.25). At 0.6 it is 0.05 below `r` on the four maps with two
    classes, whose class 2 is only about half planted atoms (precision 0.47).
  - On the other maps one class is reported.
  - `Pt,Ni` forces two classes everywhere, with median precision 0.26,
    0.41 and 0.33 at the three levels and the fitted ratio +0.13 to +0.40
    above `r`.
  - In 10.9's sense the class is found on at most one map at any level,
    so there is no `r` below which it is found.
  - The first attempt measured the same picture, so the fix did not change
    the sensitivity.

## Left open for the maintainer

1. **Sensitivity on real maps is far below the expectation.** A class at
   0.6 of the Pt intensity is not separated, and at 0.4 it is separated
   only on a third of the maps. The limit is the spread of the real
   aperture intensities, about twice their noise, with the coordination
   and radius trends above. It is not detection. Whether to model that
   spread is open: per-shell or per-coordination normalisation, or the
   density likelihood of note section 5. That decision belongs before any
   promotion (10.8). Ruled on 2026-10-12: the class fit gets one intensity
   term in the coordination number, shared by the classes (note, 10.10).
   It was implemented by hand the same day and passed the maintainer's
   tests. A fourth run, phase 9 (10.11), measures it on the 66 maps.
2. **The recovered atoms off the stored model.** These are about 20 per
   map, one bond from a stored atom, mostly in the outer half of the
   particle (median radius 0.87 of the outermost atom) but often with many
   found neighbours (median 9). Whether they are partly occupied
   sites or noise cannot be decided from the stored coordinates, and half
   maps would tell. Recall also that the expected number of noise maxima
   above the stage-B threshold is 2 to 4 per map on these boxes, above the
   target of 0.2 that the stage-A threshold meets.
3. **The present path on ten maps** matches 85% to 98% of the stored atoms.
   On ESC1_time28s it depends on the mask diameter (40.1 A gives 0.849,
   37.8 A gave 0.998). The 1.2 factor between the particle diameter and the
   mask may deserve a look on non-spherical particles. Two maps' diameters
   are 19% and 27% above the coordinate extent (EFC2_time18s, EFC2_time39s).
4. **`Pt,Ni` admissible on four maps**, one more than 10.9's 5% allows:
   EFC4_time03s, ESC4_time01s, ESIP1_time08s and ESIP2_time05s.
5. Figures are left for the Mac, from `validation.tsv`, `environment.tsv`
   and `planted.tsv`.
