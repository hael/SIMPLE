# The atomic density model: report of the unattended run

Report of the run "atomic_density_model" on dell, 2026-10-25, for the maintainer's review of the
complete raw diff. It reads on its own; the plan,
`doc/implementation_notes/planned/atomic_density_model.md`, is the contract, and its Progress
section (section 12) holds every number and log path. Each phase's diff is in
`/home/elmlundho/agent_runs/atomic_density_model/review/phase_N.diff` (N = 0 to 4). The tables,
scripts and logs of the validation are in
`/home/elmlundho/agent_runs/atomic_density_model/scratch/keep/validation/` (and the earlier phases'
in `scratch/keep/baseline`, `phase2`, `phase3`); this report does not reproduce them.

Terms. SINGLE is SIMPLE's nanoparticle workflow; `detect_atoms` finds the atoms of a 3D map. The
density model (plan section 2) renders every atom as a Gaussian of its own integrated intensity `I`
and width `B`, seen through the reconstruction's halo: the transfer `1 - a exp(-B_h s^2/4)`, two
numbers per map. Profiles describe how `ln I` and `ln B` change with an atom's number of neighbours
(coordination, 12 when fully surrounded): `beta` per neighbour for the intensity, `gamma` for the
width, `B_p` the width of a fully coordinated atom. Strong atoms are those the first detection finds
(level 1); weak atoms are those the residual recovery adds. `d_NN` is the nearest-neighbour distance.
The residual over noise is the root-mean-square of map minus model inside the particle over the
noise standard deviation (1 would be a perfect fit).

## What changed, phase by phase

**Phase 0, the baseline** (`review/phase_0.diff`: the plan's Progress entry only). The base build at ca33f927d, its
tests, and `detect_atoms` on the 66 maps without an element and with `element=Pt`, plus the mask-factor scan. Today's
detection reproduces the stored coordinates on 35 and 33 maps, not the species validation's 47 and 56. One earlier
commit explains the whole change: b57538b7b made the size rule of `discard_atoms` act, and it now removes about five
surface atoms per map. One map's result depends on the thread count, so every later run used `nthr=24`.

**Phase 1, the model as a module** (`review/phase_1.diff`). New `simple_nano_density.f90` with its tester (the
`density` sub-suite of `unit_single`): the kernel and the atom seen through the halo, rendering, the per-atom fit
with measurement variances, the position update, the profile regressions, the transfer search, the noise estimate
and the background. `atoms%convolve` gained the halo term. `gauss_filter3D` and `robust_spread` moved there from
the species module. Particles rendered with known parameters at the real maps' noise come back within every floor.

**Phase 2, the fit on the present path** (`review/phase_2.diff`). `fit_density_model` runs on every `detect_atoms`
after `discard_atoms`, and writes the products of plan section 5 (`_SIM`, `_FIT`, `_RES`, `_atoms.csv`,
`_model.txt`; `_ATMS.pdb` with occupancy `I/S` and B column `B_i`). The strong atoms and the binary products are
identical to phase 0's on all 66 maps. The model beats the equal-atom render on every map. The phase stopped at the
fade item: the particle's edge was not reproduced. The first ruling accepted that as a finding and moved the item to
phase 3.

**Phase 3, the weak atoms** (`review/phase_3.diff`). The species branch's residual recovery was moved onto the
present path, between the provisional and the full fit: a calibrated threshold, two stages, pruning by contacts
with strong atoms only. Labels go in every product (`_ATMS.pdb` chain B, `_atoms.csv`). `discover_species` was
reduced to the class call on the fitted intensities. Map correlations (third ruling) and residual peak and well
counts (second ruling) go in `_model.txt`, positions in `_resid_extrema.csv`. `species_discovery` keeps every
floor, and the pure synthetic particle gets no weak atom.

**Phase 4, the real maps and this report** (`review/phase_4.diff`). The validation of section 6 on the 66 maps and on FePt, and the
detection-gap table. The second detection pass with the fitted kernel was written and measured, then removed
(below), so phase 4 leaves no code change: `simple_nanoparticle.f90` is in its phase-3 state and
`simple_nanoparticle_utils.f90` is at HEAD. The diff of the phase is documents only: the two species notes'
pointer to the moved recovery, the plan's Progress, and this report.

## Decisions taken under delegation, and where they are recorded

Every decision is recorded with its measurement in the phase's journal (`state/journal_phase_N.md`) and in the
plan's Progress entry for the phase (section 12). The ones that shape the result:

- **`nthr=24` on every run** (phase 0, then the first ruling). `nthr` is a required key, and one map's atoms depend
  on the thread count through a near tie in the binarisation threshold.
- **The transfer is fitted between the atoms, without the background** (phase 1). Over the atom cores, `a` and the
  intensities trade off, and the plan's two fits stop at 0.54 for a true 0.6. A background smoothed at `d_NN`
  absorbs the halos' flat plateau, so with it any `a` fits.
- **Strong atoms for the profiles are those above ten times the fixed-width intensity noise** (phase 1), as the
  species branch counted. Against the joint-fit variance, no atom of the correlated-noise particle qualified.
- **The no-halo check uses the species measurement's estimator** (phase 1), chosen after the production variant
  gave 2.13 Å² against the floor's 2.2. Both numbers are in Progress for that reason.
- **The products** (phase 2). `_ATMS.pdb` keeps the detection positions, so they compare with phase 0, and the
  fitted positions are in `_atoms.csv`. `_SIM` is the module's render, which makes `single_atoms_stats`' density
  correlation 0.994 instead of 0.9997 on its Pt particle. `detect_atoms` stops when level 1 finds fewer than two
  atoms.
- **The moved recovery** (phase 3). It runs between the provisional and the full fit, calibrated on the map's own
  empty region or the residual between the atoms; half maps no longer feed it. Distances are measured from the
  provisionally fitted positions. Weak atoms are pruned by contacts with strong atoms only.
- **One `species_discovery` control changed form** (phase 3). "With `element=Pt` the present path finds fewer than
  0.90 of the aluminium" assumed a path without recovery. It now counts the strong atoms of that run, the same
  quantity as before. Its floor and name are unchanged.
- **The third ruling's "fraction below the pruning floor"** (phase 3). The present code has no `valid_corr` floor,
  so the fractions below 0.5 and 0.7 are reported instead.
- **The second detection pass is removed** (phase 4). By the second ruling's real-map measures it helps. It fails
  the synthetic floors, and no floor is relaxed. Recorded in the phase-4 journal and Progress.
- **The mask factor stays at 1.2** (phase 4). 1.3 is better by the same measures but by little, and the constant
  also sets the element path's mask, which was not measured. Left as a recommendation.
- **How the phase-4 table was measured.** Through two temporary environment switches in one build, both removed;
  four maps at a time, each at `nthr=24`, so its wall times are not single-run timings.

## Deviations from the plan, and why

- **Phase 0's ten short maps became 31.** The plan's ten came from older code. Today's code is short on 31 maps
  without an element, and the scan ran on all 31. The hand-set `mskdiam` is rounded to an even pixel count, so the
  hand factors are effective 0.97-1.34. Phase 4's table set the factor exactly.
- **The 47 and 56 of the species validation were not reproduced.** Commit b57538b7b (the size rule of
  `discard_atoms`, 2026-10-15) explains the whole difference: two bisect builds, and today's atoms are an exact
  subset of the older ones. The second ruling then made the stored coordinates a reference only.
- **Section 4's text was changed in three places in phase 1** (the transfer's residual region, the background
  left out of it, the definition of strong atoms), each because a measurement showed the text's version missing
  the floors. See the decisions above.
- **The `density` sub-suite takes 1.9 s on one thread**, more than the test policy's "well under a second" per
  sub-suite. The plan puts two particle fits in it, and they are the cost. The fast gate's wall time is unchanged.
- **Phase 2 stopped at its fade item.** The first ruling accepted the stop and moved the item to phase 3. It is
  still not met at the particle's edge (below).
- **The ruling-2 yardsticks replaced the matched fraction from phase 3 on**, including section 6 item 3's "the
  matched fraction per map decides".
- **FePt's `beta` is not compared with their 5.2 percent per neighbour.** The species-free detection finds only
  60 percent of their platinum and none of their iron, so the profile comes from an incomplete, single-element
  set.
- **Phase 4 has no code change.** The plan foresaw the second pass as one "if the table says it helps". On the
  real maps it did; on the synthetic particles it failed.

## The evidence against the baseline

The tables are in `scratch/keep/baseline/` (phase 0), `phase2/`, `phase3/` and `validation/` (phase 4); each
has a README or its script list, and the Progress entries carry the numbers.

- **Level 1 is unchanged.** The strong atoms equal phase 0's in count and coordinates on all 66 maps, and `_BIN`,
  `_CC` and `_MSK` are byte-identical. Every earlier difference from the species validation's numbers is the
  size-rule commit named above.
- **The model fits every map better than today's render.** The residual over noise inside the envelope has a
  median of 1.54 (1.13 to 2.96), against 2.11 for the equal-atom render at 13.9 Å² and 1.92 at the best single
  width. It is below the plan's 1.3 on 7 maps. Of the map correlations (third ruling), `_FIT` gives 0.783
  against 0.681 for the identical atoms, and the per-atom widths help in every resolution band.
- **The interior is reproduced, the edge is not.** The radial profile follows the map in every shell inside the
  particle on 35 maps (30 before the weak atoms). At the edge no map passes. There the map lies above the model
  (median of the largest deviations 14 standard errors), and 4679 of the 5041 residual peaks above the threshold
  lie in that band: atoms still missing.
- **The weak atoms are real-looking, and the pruning removes the rest.** Synthetic particles of one species get no
  weak atom against an expected false count of 4.6 to 6.1. The real maps get 1379, a median of 11 per map,
  almost all one bond from a stored atom, weaker (occupancy 0.57 against 0.89) and as sharp (36 against 37 Å²).
  287 of the 305 surface atoms the size rule removes sit on a residual peak, but the strong-only contact pruning
  removes them again; 11 survive.
- **The wells between the atoms are too shallow in the model.** The mean level there is right (within 0.25 noise
  units on 58 of 66 maps), but where the model has its wells the map is a median 1.26 times deeper (58 of 64
  maps).
- **Across time points the atoms agree and the parameters do not.** Consecutive time points share 0.953 of their
  atoms, and ln `I` of matched atoms correlates beyond noise on 41 of 51 pairs. `beta`, `a` and `B_h` vary within
  a particle by about their own size; only `B_p` holds (range 6.5 against 37.8 Å²).
- **The detection-gap table.** The second pass lowers the residual peaks and wells by a median of 10.5 per map and
  the residual by 0.018, and fails `species_discovery`'s synthetic floors (`beta`, width rise, `B_h`, false weak
  atoms). The mask at 1.3 lowers them by 3 at the same residual.
- **FePt.** The model fits a map from outside SINGLE without an empty region: residual-based noise, the transfer
  fitted, 55 weak atoms, all on their platinum and none on their 618 removed peaks. The species-free detection
  finds 9909 of their 16 627 platinum atoms and none of their 6569 iron, missing about 30 percent at every radius.
- **Tests at the end.** Fast gate 15 of 15; the five CTest entries pass with `species_discovery`'s floors as in
  phase 0 (its metrics listed beside phase 0's in `scratch/keep/phase3/species_discovery_beside_phase0.tsv`);
  `unit_single` passes in Release and in Debug with bounds checking; products are identical at 11 and 24 threads.
- **Time.** `detect_atoms` takes a median 6.1 s per map at 24 threads (1.1 s in phase 0), the fit and the recovery
  being largely serial.

## What is left open for the maintainer

- **The size rule and the pruning of weak atoms.** The recovery finds 287 of the 305 surface atoms the size rule
  removes, and the strong-only contact score prunes them again. Counting weak neighbours, or letting the residual
  peak vouch for a surface atom, would change that. It is the first ruling's open decision.
- **The edge.** Atoms are still missing there (4679 residual peaks), and the outside band is negative: halos the
  model does not match.
- **The transfer.** `a` reaches its cap of 0.99 on 15 maps, `B_h` ranges 21 to 682 Å² and varies within a
  particle by its own size, and the model's wells are a quarter too shallow. The phase-1 probe found `B_h` and the
  intensity scale on a flat valley. A transfer shared by the time points of a particle, or fitted with the edge's
  atoms in, are the obvious next measurements.
- **The measured spreads** (`tau_B`, the intensity spreads) come out zero on many maps, because the correlated
  noise inflates the measurement variances beyond the observed scatter.
- **The mask factor.** 1.3 is recommended for the species-free path. Measure it with `element=Pt` before changing
  `MSKDIAM_OVER_DIAM`.
- **A second pass that keeps level 1's atoms.** The tested pass replaced level 1 and merged neighbours. A variant
  that uses the fitted kernel only to propose candidates for the recovery was not tried.
- **The neighbour gate** (`min_nbrs` = 3, plan section 9). It admitted 12 percent of the weak atoms and no false
  atom on the synthetic particles. It is left on.
- **`image%win2arr_rad` returns one pixel too many** (`src/main/image/simple_image_geom.f90`). Every `valid_corr`
  includes one leftover pixel of the previous window. It was found in phase 3 and not changed, being outside the
  run's files.
- **FePt** needs element-aware detection for its iron, and its noise has no empty region to come from. Its numbers
  are a first look, not a validation.
- **`autorefine3D_nano`** now aligns against the model's `_SIM` (phase 2), which section 5 intends. It was not
  exercised in this run.
- **The `density` sub-suite** (1.9 s) is above the test policy's guideline per sub-suite.

# The second run: the measured transfer, a species-free `detect_atoms`, `detect_atomic_species`

Report of the run "atomic_density_model_2" on dell, 2026-10-26 to 2026-10-27, phases 5 to 9 of the plan (its
section 12). It reads on its own. The plan's Progress (section 13) holds every number and log path, and the
maintainer's rulings 4 to 7 are recorded in section 12, each where it applies. The tables are under
`/home/elmlundho/agent_runs/atomic_density_model_2/scratch/keep/`: `validation/` holds the phase-5 baseline, and its
`README_second_run.txt` and `../README.txt` name the later phases' tables beside it. The diffs are in
`/home/elmlundho/agent_runs/atomic_density_model_2/review/`.

Terms. The transfer `T(s)` is what the reconstruction does to an atom at each spatial frequency `s`. The first run
modelled it as a notch with two numbers; this run measures it from each map, shell by shell, as the projection of the
map onto its fitted atoms rendered without any transfer. A round is one measurement of that curve followed by a full
fit of the atoms with the curve held fixed. The gauge is the convention that fixes what a map alone cannot: the
curve's scale and its Gaussian tilt are set to unity over 3.1-2.0 Å, so the atoms' absolute brightness `S` and width
`B_p` are numbers in that gauge. A class call sorts the atoms into intensity classes (candidate species). A posterior
is an atom's probability of belonging to a class. The forced two-class fit is the call with a two-element list;
admissible means its classes are separated by at least 3 noise widths (D ≥ 3) and each holds enough atoms.

## What changed, phase by phase

The driver writes one diff per phase. Phases 6 and 8 were reopened after later phases had run, and the fence's
snapshot was then reset, so their diffs are incomplete:
- `phase_6.diff` holds only the last session: the seventh ruling's form of the rounds and the gauge fix of the intensity noise.
- `phase_8.diff` is empty: its work was already in the snapshot when the phase was resumed.

The complete change of the run is `full_after_8.diff` plus `phase_9.diff` (documents only).

**Phase 5, the baseline** (`phase_5.diff`).
- **The tree as handed over** built on dell and passed the fast gate and the five CTest entries.
- **`single_workflow_fcc` failed:** the final map's Fourier shell correlation (FSC) with the truth fell to 0.143 at 5.73 Å, against a floor of 5.00 Å; HEAD gives 4.41 Å. The fourth ruling gave the nanoparticle refinement loop's atomicity prior free amplitudes per atom (its positions and widths from the model). That was applied in `simple_rec3D_pcg_strategy.f90` and `simple_pcg_atomic_prior.f90`, but the test still failed, and the ruling moved the judgement to phase 9.
- **The baseline on the 66 Pt maps and on FePt** used the first run's scripts:
  - the model's residual over the noise: median 1.537;
  - 59 residual peaks per map;
  - the wells within 10 % of the model on 9 maps;
  - the map correlation of the full model (`_FIT`): 0.783;
  - 7.0 s per map.
- **The maintainer's pooled transfer** from the Mac was reproduced on dell (median difference 0.004).

**Phase 6, the measured transfer** (its whole code is in `full_after_6.diff`, its last session in `phase_6.diff`). The transfer is measured per map and held fixed during each fit.
- The measurement is the projection per shell of 1/50 Å⁻¹, smoothed by a running median of three shells, and held at the last reliable value where the atoms carry less than 1 % of their peak power.
- The halo the curve gives each atom is rendered by one FFT per colour class of the parallel fit.
- `_transfer.tsv` holds the curve, and the key `transfer_tab` hands a known curve to the fit; only the unit tests use it.

It took three stops and three rulings to reach the final form.
- **The first form alternated measurement and fit.** It has an exact degeneracy and ran away on two maps and on FePt.
- **The fifth ruling** fixed the gauge and allowed one measurement. That left the synthetic particles' scale 0.22 off.
- **The sixth ruling** asked for rounds stopped by the residual. On the white-noise particle the residual is flat at its noise, so it could not stop them.
- **The seventh ruling** stops the rounds when the curve changes by 0.01 or less, with at most six rounds. The last round is kept unless its residual is more than 0.1 % above the best round's, and the measurement's regressor is the model as fitted.

**Phase 7, `element` out of `detect_atoms`** (`phase_7.diff`).
- **Species-free detection:** detection takes no element on any path. `element`, `discover_species`, `nspecies`, `vol_even` and `vol_odd` are refused with a message that names `atoms_stats` and `detect_atomic_species`.
- **Contact ceiling:** the ceiling of `discard_atoms` is now taken from the data, as the median contact score of the inner 30 % of the atoms. It is 12 for Pt; the CdSe wurtzite, zincblende and rocksalt test particles get 4, 4 and 6.
- **`single_atoms_stats`** runs `detect_atoms` without an element and `atoms_stats` with it.

**Phase 8, `detect_atomic_species`** (`full_after_8.diff`). A new program of `single_exec`.
- **What it reads:** the products of a `detect_atoms` run (`pdbfile=` its `_ATMS.pdb`, with `_atoms.csv` and `_model.txt` beside it) and the map. It takes the species keys `element` (a list such as `Pt,Ni`, brightest first, whose length is the class count), `nspecies`, `vol_even` and `vol_odd`.
- **The class fit** is `ln I = ln S_k + beta (c - 12)`, with one `beta` per particle (`c` is the coordination).
- **Each atom's posterior** carries its measurement variance and the coordination spread `tau_I(c)`.
- **What it writes:** the `_species` files as before, with the posterior's standard error and the coordination spread as two more columns.
- **In `detect_atoms`:** `discover_species` and its code are gone.
- **The test:** `species_discovery` is now the test of the pair. The phase stopped on two of its checks, and the sixth ruling reframed them: no confident call (largest posterior 0.9 or more) may be wrong, and the softness check's floor scales with the class gap.

**Phase 9, the real maps and this section** (`phase_9.diff`). The final code was run on the 66 maps and on FePt, with `detect_atomic_species` on every map; `single_workflow` was judged again. This section was written, and the species status and unfinished notes now point to the new program and the moved keys.

## Decisions taken under delegation, and where they are recorded

Each is in the phase's journal (`state/journal_phase_N.md`) and in its Progress entry; those that became rulings are in section 12.
- **The build environment** (phase 5). The run uses its own build, `repo/build`, in Release with tests, through `scratch/env.sh`. Phase 5 also built HEAD once to diagnose the failure of `single_workflow_fcc`.
- **The transfer's representation** (phase 6):
  - shells of 1/50 Å⁻¹;
  - a running median of three shells and the hold below 1 % of the atoms' power;
  - the halo by FFT per colour class, with the FFT threads capped at 16, because the system FFTW is slow above that on this machine;
  - `transfer_tab` registered as a string key, because the file-key registry refuses `.tsv`.
- **The provisional fit stays plain** (phase 6). A transfer measured on the strong atoms alone took up the missing atoms' density and lost the aluminium recovery of the light test particle.
- **The gauge band of 3.1-2.0 Å** (phase 6). It became the fifth ruling.
- **The intensity noise takes the gauge's scale** (phase 6, reopened). The strong-atom test compared gauged intensities with an ungauged noise. On FePt and ECT2_time15s it kept the profiles of the provisional fit through every round.
- **The data-derived contact ceiling and its definition** (phase 7): the median contact score of the inner 30 % of the atoms, rounded.
- **The design of `detect_atomic_species`** (phase 8):
  - a new domain module, `simple_atomic_species.f90`;
  - `discover_species` removed outright: an unknown key is refused, as release 4 requires;
  - the half maps' aperture intensities scaled to the fitted intensities by the map's own apertures;
  - the posterior's standard error by the delta method.
- **`tau_I` is measured within the classes** (phase 8). The `tau_I` that `detect_atoms` writes is fitted under one species. On a two-species particle it holds the species' separation (0.32 on the Pt/Ni alloy), and every two-class call was then inadmissible. The sixth ruling kept the decision.
- **The test of the pair** (phase 8):
  - one `detect_atoms` run per case;
  - `pt_disc` dropped: with no element in detection it is the same run as `noel_disc`;
  - the old products check became "the species call leaves the detection products unchanged".
- **The extra measurements of phase 9:** `element=Pt,Ni` on every map, the check the species work's open note asks for; and FePt's classes against Miao's model by coordination.

## Deviations from the plan, and why

- **Phase 6 is not "measured once and fixed".** The plan's single measurement left the synthetic particles' scale and slope outside their floors. The rulings replaced it with rounds, each holding the curve fixed during its fit.
- **Phase 7's "the 66 maps equal phase 6" did not hold on 14 maps.** The data-derived ceiling applies to every run, not only to runs that gave an element. On all 14 maps the residual is lower, and they have fewer residual peaks.
- **Phase 8's coordination spread is not the one `detect_atoms` writes.** The reason is above; the sixth ruling confirmed it.
- **Phase 8's element check became "no confident call is wrong"** (the sixth ruling), and the softness check's floor became the softness gate times the class gap.
- **The species work's planted-class runs were not repeated.** Phase 9's scope is the species call on every map; the unfinished note says so.

## The evidence against the baseline

- **The final detection is the phase-6 run.** On all 66 maps and on FePt, the final code's `_atoms.csv`, `_model.txt` and `_transfer.tsv` are byte-identical to phase 6's final run (`validation_p9/runs.tsv`).
- **The fit is better than phase 5 on every map, and better than the single measured curve on 63.**
  - the residual over noise: median 1.284 against 1.537 in phase 5, lower on 66 of 66 maps, and below 1.3 on 37;
  - 16 residual peaks per map against 59;
  - the wells within 10 % on 39 maps against 9;
  - `_FIT` correlation 0.840 against 0.783.

  Against the single measured curve of the fifth ruling, the residual is lower on 63 maps and more than 1 % higher on none (`phase6/validation_r7/summary_r7.txt`).
- **The synthetic particles meet the fifth ruling's floors at the round the rule keeps.**
  - The uniform particle settles in 5 rounds: `ln S` -0.018, `beta` 0.0509 against 0.05, curve error 0.050.
  - The gradient particle reaches the 6-round limit: `ln S` 0.042, `beta` -0.0246 against -0.03, curve error 0.050.
- **`species_discovery`'s density cases pass.** The `beta` error is 0.0003 and 0.0038, against 0.012 and 0.022 with the single curve.
- **The measured curves look like the reconstruction's band-pass.**
  - Medians by band: 0.22 at 50-6.2 Å, 0.20 at 6.2-3.3 Å, 0.905 at 3.3-2.0 Å, 0.30 at 1.9-1.6 Å.
  - The edge, where the curve crosses one half, has a median of 3.33 Å and is near 5 Å on four particles. Its spread within a particle is 0.94 Å (median over the particles), and the band 3.3-2.0 Å varies within a particle by only 0.07.
- **One species on every real map.**
  - With the class count free, `detect_atomic_species` finds one class on all 66 maps.
  - With `element=Pt,Ni` the forced two-class fit is inadmissible on all 66. The species work's open note asked for at least 95 %, and the four maps that accepted `Pt,Ni` in its phase 7 now refuse it.
  - Most forced fits leave the second class empty. Nine put atoms in both classes, and there many calls are uncertain: EFC4_time18s has 266 of 304 atoms below a posterior of 0.9, FCHCS1_time02s 185 of 246, ECT2_time15s 27 of 340 (`validation_p9/summary_p9.txt`).
- **FePt** (`validation_p9/fept/`):
  - With the class count free: one class.
  - With `element=Pt,Fe` the two-class fit is inadmissible (D 1.38). Detection finds none of Miao's 6569 iron atoms: all 9916 matched atoms are platinum in his model. So the forced split divides platinum, and this is not a test of the class call.
  - The residual over noise is 0.955, against 0.963 in phase 5; the `_FIT` correlation 0.739, against 0.709.
- **The tests:**
  - **`single_workflow_fcc` passes:** FSC 0.143 against the truth at 3.58 Å (floor 5.00 Å; phase 5 5.73 Å; HEAD 4.41 Å), map correlation 0.984.
  - **`single_workflow_wurtzite` passes:** 2.49 Å (phase 5 1.59 Å), within its profile (median 1.59 Å, spread 0.97 Å, worst 3.0 Å), correlation 0.978.
  - The fast gate passes 15 of 15. The five CTest entries pass, the `density` and `species` sub-suites at 1, 11 and 24 threads with identical prints.
- **Time.** `detect_atoms` takes a median of 8.8 s per map at 24 threads when the machine is quiet (7.0 s in phase 5, 3.8 s with one curve). The rounds are the cost.

## What is left open for the maintainer

- **The curve oscillates on many real maps.** The rounds reach their limit on 31 of 66 maps and on FePt, and the guard keeps an earlier round on 40, so a measurement from the fitted model is not a fixed point on many real maps. A damped update, a curve shared by a particle's time points, or the product of the first two rounds would be the next things to measure; the plan's per-map rule rules out the second without a new ruling.
- **EFC2_time10s keeps the profiles of its provisional fit.** After the first measurement its atoms are about five intensity-noise units bright, below the strong-atom threshold of ten.
- **The posterior's standard error is linearised.** It exceeds 1 on two forced fits (ESC5_time08s 1.08, FCHCS1_time06s 1.05). That happens where the class gap is several noise widths and the posterior is near one half. A bounded form (in the log-odds) would read better.
- **An inadmissible forced fit still writes its labels.** On FePt, `element=Pt,Fe` names 8730 platinum atoms iron, 3404 of them with a posterior of at least 0.9, while `_species.txt` says the fit is inadmissible. Whether the program should refuse to label, or fall back to one class, is a product decision.
- **FePt's iron is not detected.** Species-free detection finds 60 % of Miao's platinum and none of his iron, so species on FePt is first a detection question.
- **`S` and `B_p` are gauge-dependent.** `B_p` ranges from 1.1 to 93 Å² across the maps. `B` by coordination now falls from about 39 Å² at 4-5 neighbours to 31.9 Å² at 12, where phase 5 was flat.
- **The species work's planted-class validation** (the open note's expectations 2 to 4) has not been run with the new program.
- **The first run's open items** that this run did not touch:
  - the size rule and the strong-only pruning of weak surface atoms;
  - the missing atoms at the particle's edge;
  - the mask factor.
