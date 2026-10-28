# A nano context for model_cavgs_rejection: class-average rejection for nanoparticle time series

Proposal, 2026-10-28, written for one autonomous run. Written against the working tree of that day, which
holds, uncommitted, the three things this plan builds on: the rank-based AUC in
`src/main/cavg_quality/simple_cavg_quality_stats.f90`, the bank-free logistic core
`src/main/quality/simple_quality_logistic.f90` that the class-average learner now calls (with its tester and
`unit_numerics` sub-suite "quality logistic"), and the lattice measurements of
`src/main/frame_quality/single_frame_quality_{types,feats,floor}.f90`, which this plan re-homes and the run
deletes from where they are.

    simple_exec prg=model_cavgs_rejection quality_mode=analyze quality_context=nano projfile=... mskdiam=...
    simple_exec prg=model_cavgs_rejection quality_mode=learn filetab=... fname=...
    simple_exec prg=model_cavgs_rejection quality_mode=apply quality_model=nano projfile=... mskdiam=... prune=yes

## 1. What this is

The nanoparticle workflow (`tseries_track_particles`, `trajectory_denoise`, `trajectory_import`,
`analysis2D_nano`, `autorefine3D_nano`) selects frames by hand on the class averages of `analysis2D_nano`:
classes whose aligned average shows atomic contrast are kept, their frames go to 3D. This plan makes that
selection a model, the same way the cryo-EM path does it: `model_cavgs_rejection` gains a `quality_context`
`nano`, with its own hard gates, a feature bank extended by a lattice family, a normalization against the
project's floor of lattice-free classes, and a learned preset `nano` trained from the maintainer's class
labels and promoted when it reproduces them and the 3D reconstruction from its selection matches the one
from the hand selection. Nothing of the per-frame approach survives except the measurements it developed;
section 9 records what it established.

## 2. Why the class average is the unit

The maintainer's judgement works on aligned averages, the labels he makes are class states, and the
reconstruction he validates runs from class selections. The per-frame experiment of 2026-10-27 showed the
other side: frames below single-frame detection that a class selection keeps cannot be told apart from the
ones it rejects by anything measured on a frame, raw or denoised (section 9). An aligned average of twenty
frames has the signal a frame lacks, and the cryo-EM machinery — `analyze` from manual `cls2D` states, the
learner, `evaluate`, `promote`, `apply` mapping class states to particles and pruning — already does the whole
loop. The development is a context and a bank, not a program.

## 3. The nano context

`quality_context=nano` is the fourth context beside `chunk`, `pool` and `sieve`
(`simple_cavg_quality_types.f90`, `validate_quality_context`, `cavg_hard_reject_for_context` in
`simple_cavg_quality_feats.f90`; the UI choices of `quality_context` and `quality_model` in
`simple_ui_params_common.f90`; the per-mode checks in `simple_parameters_phases.f90`). It differs from the
cryo-EM contexts in three places.

**Hard gates.** The population gate of the cryo-EM contexts would reject the good classes: on the delivered
labels the selected classes are the smaller ones (the lattice-free frames pile into large junk classes;
section 5), and the FRC "resolution" of a nanoparticle class is quantized at a few shells and carries
nothing. `nano` keeps the validity gates only: invalid pixels; an empty
class (`pop <= 0`); both local variances at or below `EPS`. The foreground-geometry gates (Otsu after a
30 A low-pass, centroid outside the mask, component outside the mask) are not applied: at 0.358 A per
pixel in a 57 A box a 30 A low-pass leaves two shells, and the particle fills the mask.

**Feature bank.** The bank is per context. `chunk`, `pool` and `sieve` keep their fourteen features,
unchanged in number, order, extraction constants and presets. `nano` has those fourteen plus a `lattice`
family of six, measured on the class average by the code of `single_frame_quality_feats.f90`
re-homed as `src/main/cavg_quality/simple_cavg_quality_lattice.f90` (a shared module, as the backend it
serves; section 10 has the naming question): `lattice_peak` (log10 of the largest 3 by 3 smoothed power
over its ring mean, over 3 A to Nyquist with the graphene shells excluded), `spot_count` (local maxima of
that ratio above the floor's threshold), `spot_arc` (the strongest spot's extent along the ring over
across it; straight fringes give a round spot), `spot_diam` (the projected-sphere diameter fitted to the
strongest spot: the particle's size per class), `spot_snr` (the fit's amplitude over background),
`band_presence` (band-passed energy in the central disc over the annulus, unmasked). The fourteen
cryo-EM features are computed for `nano` with their resolution constants scaled to the specimen:
`FOREGROUND_SEG_LP`, `SIGNAL_METRIC_LP` and the overfit band are the cryo-EM values divided by a context
factor so that the foreground segmentation sees the particle and the band-pass sees the lattice — the
factor and the values are fixed by the study of section 6, not guessed here; the fallback is to drop the
features that stay degenerate.

Bank size: `CAVG_QUALITY_NFEATS` becomes the maximum (20); the active feature list is a property of the
context (`cavg_quality_nfeats(context)`, the names, families and directions by context), and the model
spec, the training-table writer and reader (columns keyed by name, the `quality_context` line selecting
the bank), the normalizer, the feature-policy masks and the model-file reader (`feature_index_from_token`)
all follow the context's list. The cryo-EM presets' coefficient arrays are written for fourteen
(`spec%linear_coefficients(1:14) = [...]`, `set_pairwise_interactions_for_feature_count(spec, 14, ...)`),
and the gate of section 7 holds them byte-identical.

**Normalization.** The lattice family is normalized against the project's floor: a two-component mixture
(`xd_gmm`, the two starts of `fit_k`, the code of `single_frame_quality_floor.f90`) over the classes'
`lattice_peak`; the lower component is the floor, the classes without a lattice; `lattice_peak` enters as
distance from the floor's mean in floor sd, the other lattice features by median and MAD of the floor
classes with a minimum scale (one count for `spot_count`, a tenth of the interquartile range otherwise);
the guard (floor weight at least 5% of the classes, separation at least 4 floor sd) falls back to the
cryo-EM median/MAD normalization and records it in the table header. The fourteen cryo-EM features keep
their median/MAD normalization. The clip is the bank's `CLIP_Z`.

**The relational feature** (`simple_cavg_quality_relations.f90`, each average correlated against its k
nearest neighbours after in-plane alignment) is kept: averages of one particle correlate, junk does not.
Its correlation band is a property of the model spec (`relational_corr_hp`, `_lp`, `_trs`), so the nano
preset carries nano values (high-pass about 15 A, low-pass at the lattice, 1.5 A, shift range a few
pixels), fixed by the study.

**Policies.** The learner's feature policies for `nano` are `lattice`, `lattice_plus_signal`,
`lattice_plus_microchunk`, `all`; the cryo-EM policies are unchanged.

**Presets.** `nano`, the learned model, added by `promote` when section 8 is met; until then `analyze` and
`evaluate` run with the context and no preset (the hard-gate-only decision), as the program already allows.

## 4. Pipeline position

`tseries_track_particles` → `trajectory_denoise` → `trajectory_import` → `analysis2D_nano` →
`model_cavgs_rejection quality_mode=apply quality_model=nano prune=yes` → `simulate_nanoparticle` →
`autorefine3D_nano`. `apply` annotates `cls2D` with `quality`, `quality_cluster` and `accept`, maps the
class states to `ptcl2D` and `ptcl3D`, prunes, and writes the project once; the frames keep their original
index in `indstk` through the prune, and any time-windowed consumer reads it (the project's row order is
no longer time). The particle diameter per class (`spot_diam`) is written to the class rows and its median
over the accepted classes to the table header and `os_out`, beside the mask diameter, which the 3D stage
takes from there instead of from `estimate_diam`.

## 5. Inputs and labels

A labeled dataset is one `analysis2D_nano` run with the maintainer's selection. As delivered (2026-10-28,
`~/Processing/NanoX/NanoX_trajs/<dataset>/<trajectory>/`), each holds `trajectory_denoised.mrcs` (the stack
the 2D ran on; its header pixel size is wrong, 0.0358 A, so every import passes `smpd=0.358`),
`subtr_frames.mrc` (the raw frames, not needed here), `trajectory_denoise_flex.txt`, and the two project
fields written as `oris` text, one row per entry: `cls2D_field.txt` (`class pop res corr state`; the
maintainer's selection is the `state`) and `ptcl2D_field.txt` (the full particle record per frame: `class`,
the in-plane alignment `e3 x y`, `state`, `indstk`, `eo`, `corr`). No project file and no class-average
stack are delivered; the run rebuilds both: `import_particles stk=trajectory_denoised.mrcs smpd=0.358
ctf=no oritab=ptcl2D_field.txt` (the `oris` reader takes the text table), `make_cavgs ncls=<classes>
refs=cavgs.mrc` from the stored assignment and alignment (which fills `cls2D` with population,
resolution and correlation and registers the stack in `os_out`), then the class states from
`cls2D_field.txt` through `selection oritype=cls2D infile=<the 0/1 column>`, checked against the
delivered populations class by class. The regenerated averages are the originals up to the even/odd
weighting of `make_cavgs`; if the original `cavgs_iter*.mrc` stacks become available they are used
instead and the regenerated ones compared to them.

Checked on delivery (2026-10-28): all three are complete and consistent: `71-17/traj04` (649 classes,
185 selected, 12,971 frames), `71-6/traj01` (684 classes, 425 selected, 13,662 frames) and `71-6/traj04`
(684 classes, 264 selected): every class's population equals its frame count in `ptcl2D`, every frame's
state equals its class's state, the alignment is present on every frame, `indstk` is the row,
populations run from 0 to 50 (two empty classes in each `71-6` set). The classes are not time chunks:
only 0.5% of consecutive frames share a class, the frames of one class are a median of 270 to 300 frames
apart and span most of the trajectory. That is not what the code reads as doing: `refine2D_nano` sets
`tseries=yes` and `nptcls_per_cls=20`, the class assignment starts as consecutive chunks
(`init_tseries_refs`), and with `nptcls_per_cls` defined the greedy search is `strategy2D_tseries`, which
lets a frame move only to the nearest populated class before or after its own per iteration, so after
the 15 iterations of the run (`updatecnt`) a class should hold frames from within a few hundred frames of
its chunk, not from everywhere. Either the labeled runs were made with other settings or an earlier
build, or the constraint does not hold in practice; Phase 0 settles it with one `analysis2D_nano` on a
hold-out trajectory and the same statistic. For this plan the scattered assignment is the easier case:
a class is then a group of alike frames from the whole trajectory, as a cryo-EM class is, and its average
is a cleaner object than a 20-frame time chunk. Three labeled trajectories are the minimum for
leave-one-out; more make the held-out numbers mean more.

What the delivered `cls2D` fields already say, before any image is read (class-level AUC against the
selection, empty classes excluded): the stored class correlation `corr` (the frames' mean correlation to
their average, what the cryo-EM bank reads as `corr_frc_proxy`) 0.965, 0.907, 0.908 on the three sets;
population 0.23, 0.44, 0.32, that is inverted: the selected classes are the smaller ones (median 15 to 19
frames against 20 to 22), because the lattice-free frames pile into large junk classes; the FRC
resolution 0.47 to 0.62, nothing (its values are quantized at a few shells). Two consequences are already
fixed by this: the population hard gate of the cryo-EM contexts would reject the good classes and is off
for `nano` (section 3), and the stored correlation alone is the baseline every model has to beat
(section 6). Labeled datasets train and are evaluated leave-one-out; every other trajectory of
`~/Processing/NanoX/NanoX_trajs` is a hold-out whose `analysis2D_nano` output the run makes itself and
scores for review.

## 6. The study, before the Fortran

On the labeled datasets, in order:

1. `model_cavgs_rejection quality_mode=analyze quality_context=chunk` as it is (the population gate off
   through the existing `threshold_pop` switch, `mskdiam` from the project), to see what the fourteen
   cryo-EM features do on nanoparticle averages: the training table's per-feature AUC against the labels,
   and which features are degenerate at these resolutions;
2. the same with the resolution constants scaled (a numpy reproduction of the six affected features at
   two or three candidate factors, on the class-average stacks), to fix the factor of section 3;
3. the lattice family in numpy on the class-average stacks (`single_frame_quality_feats.f90`'s definitions,
   the prototype scripts of the previous run in `~/agent_runs/runs/frame_rejection_mac` and the
   trajectory folder's `trajectory_focus_prototype/` are the references), the floor over classes, per-feature
   AUC against the labels, and the floor's separation per dataset;
4. a numpy pairwise-logistic fit mirroring the learner, leave-one-dataset-out, with and without the lattice
   family, against three baselines: the stored class correlation alone (AUC 0.91 to 0.97 on the delivered
   labels, section 5: the bar), the lattice floor alone (`lattice_peak` of the class at 4 floor sd) and the
   relational feature alone. If the stored correlation alone matches the fit, the honest product is a
   per-project threshold on that one number with the floor as its reference, and the bank is diagnostic.

Decision rules, written in the journal: a feature family enters the nano bank when its best member's
leave-one-out AUC against the labels is at least 0.6 and the fit's class-level balanced accuracy drops by at
least 0.01 without it; a cryo-EM feature that is degenerate on nano averages at the chosen factor (robust
spread zero on more than half the datasets) is dropped from the nano list.

## 7. Tests and gates

- The cryo-EM path is byte-identical: `model_cavgs_rejection` `learn` (model file and report) and `evaluate`
  on saved class-average training tables, and `analyze` on a class-average project, before and after, all
  outputs compared byte for byte; the saved tables are made in Phase 0 from the repository's test fixtures,
  failing that from synthetic tables as the previous run did (`~/frame_rejection_mac/scripts/cavg_gate.sh`);
  `unit_numerics`' three sub-suites ("cavg quality relations", "cavg quality selection", "quality logistic").
- A `nano` sub-suite of `unit_numerics` (no new CTest entry; `scripts/check_test_registry.py`): the
  lattice family on synthetic class averages (a lattice at a known spacing, the same low-passed, noise at
  the SNR of a 20-frame nanoparticle average): band and spot shell, separation, the floor and its guard
  (an all-lattice project and an all-empty one reject nothing), the diameter against disc and sphere
  objects (the calibration of the previous run as assertions), `spot_arc` on straight against bent
  fringes, `band_presence` mask-free, the per-context bank sizes and name lists, the table round trip for
  both bank sizes, a model file round trip for a 20-feature spec, and the gates of the context.
- The commander: `analyze` on a labeled nano project writes the table with the context's columns and the
  floor in its header; `apply` with a nano model maps states, prunes, keeps `indstk`, writes once.

## 8. Acceptance and promotion

The learned `nano` model is promoted when, on every held-out labeled dataset: its class-level agreement with
the maintainer's labels (balanced accuracy at the model's threshold, and AUC) is not below the best of the
three baselines of section 6, step 4 (the stored correlation alone being the one to beat); and, on at least one held-out trajectory, `autorefine3D_nano` from its
selection gives a map and atomic model that match the hand selection's (the FSC bands outside the lattice
band within the run-to-run spread, atoms within 1 A of the hand selection's model at the same fraction).
A model that matches the baselines but does not beat them is promoted only if it beats them on the 3D
criterion; one that fails either is reported, not promoted, and the floor baseline is recorded as the gate
the maintainer can use by hand.

## 9. What the per-frame experiment established (2026-10-26 and 27)

Kept here because it fixes several design choices above; the full record is in
`~/agent_runs/runs/frame_rejection_mac/` and in the prototype folder of `71-17/traj01`.

- The lattice measurement: the maximum of the 3 by 3 smoothed power over its ring mean, over a
  species-agnostic band with the graphene shells excluded, is the right statistic; the smoothing matters
  (a Bragg spot of a 30 A particle spans two Fourier pixels at this sampling); the plain band maximum (the
  collaborators' quantity) does as well on raw frames and is dose-dependent.
- The floor: the classes or frames without a lattice form one tight population; a two-component mixture
  finds it, Otsu's midpoint does not (it splits the broad visible population), a white-noise null sits
  below the real floor, and the floor's left side underestimates its spread. Normalize against the fitted
  floor; guard on its weight and separation.
- On single raw frames the floor is flat along a trajectory; on flex-denoised frames it drifts with each
  frame's pooled population and no threshold holds, while the rank still does. On class averages the
  question does not arise: the floor is over classes of one project.
- Lattice visibility flickers frame to frame as the particle tumbles (runs of consecutive strong frames
  have a median length of one); a running median alone rejects the flashes. Class averages of time chunks
  absorb this by construction.
- Frames below single-frame detection that a class selection keeps are indistinguishable from the ones it
  rejects by every per-frame feature tried, on raw and denoised frames, and by their correlation to the 3D
  map (AUC 0.46 to 0.60); the 3D from the floor rule's half of the frames matched the hand selection's
  outside the lattice band. The class average is where the judgement can be made.
- The Bragg-spot diameter (projected-sphere model fitted to the spot) is unbiased on discs and 12% low on
  spheres at a single frame's SNR; on an average it is the better size estimate and replaces
  `estimate_diam`.
- Defects found and left outside the per-frame work: `make_cavgs` in time-series mode re-splits `ptcl2D`
  halves by row parity even when halves exist (`simple_make_cavgs_strategy.f90`), so `center2D_nano`
  cannot follow a prune; the MRC writer sets a stack's z cell before its frame count (cosmetic); the
  particle record has no slot for a per-row probability or diameter (class rows do, in `cls2D`, which
  this plan uses).

## 10. Open rulings

- The home and name of the lattice module: `simple_cavg_quality_lattice.f90` in `src/main/cavg_quality/`
  (shared, as the backend it serves and the crystalline specimens it measures are not only nanoparticles),
  or `single_cavg_quality_lattice.f90` under `src/main/nano/` by the SINGLE naming rule.
- Whether `estimate_diam` is retired from `analysis2D_nano` once the class diameter exists.
- The scaling factor of the cryo-EM resolution constants for `nano` (fixed by the study; the ruling is
  whether to keep those features at all if only the lattice family carries weight).

## 11. Phases

Every phase builds in the run's out-of-tree Release build and the worktree's Debug build with the fast
gate when sources changed, runs the touched sub-suites, writes the journal, leaves no job running, writes
its status file; evidence under `scratch/keep/phase_<n>/`.

- **Phase 0, audit.** Builds; the labeled datasets found, read and tabulated (classes, populations, label
  counts); the cryo-EM gate's inputs saved (section 7). Exit: both builds, the fast gate, the audit table,
  the gate baselines.
- **Phase 1, the study** (section 6) in numpy and with the existing program; the bank decision and the
  constants in the journal and in section 14 of this note. Exit: the AUC tables, the leave-one-out fits,
  the decision.
- **Phase 2, the context and the bank.** Per-context bank size and lists; the lattice module re-homed
  from `src/main/frame_quality/` (then that directory deleted); the nano gates, constants, normalization
  with the floor, policies, relational parameters; the training-table and model-file formats by context;
  the tester; the cryo-EM byte-identity gate. Exit: tests and gate pass, fast gate passes.
- **Phase 3, learn and evaluate.** `analyze` on every labeled dataset with `quality_context=nano`; `learn`
  leave-one-out; `evaluate` on the held-out dataset against the two baselines; the hold-out trajectories'
  `analysis2D_nano` outputs made where absent and scored for the maintainer's review (tables, the
  selected and rejected class-average stacks the program writes). Exit: the reports and tables.
- **Phase 4, the 3D criterion.** On one held-out trajectory: `apply` with the learned model (an external
  model file, since nothing is promoted yet), `simulate_nanoparticle` from the class diameter,
  `autorefine3D_nano`; the same from the hand selection and from the lattice floor baseline; compared as
  section 8 says. The promotion decision; if promoted, the `nano` preset through `promote` and the UI
  choice. Exit: the comparison table, the decision.
- **Phase 5, close.** `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` gains the nano context
  (gates, bank, normalization, presets, the labels, the commands); the skills `simple-cluster-cavgs-quality`
  and `simple-main-nano` updated; the report
  (`doc/implementation_notes/planned/nano_cavgs_rejection_report.md`: the learning first, the AUC and
  3D tables, what was built and run, what is left); `full.diff`.

## 12. File table

| file | phases |
| --- | --- |
| `doc/implementation_notes/planned/nano_cavgs_rejection.md` (this note: the journal in section 14, nothing else) | each |
| `doc/implementation_notes/planned/nano_cavgs_rejection_report.md` (new) | 5 |
| `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` | 5 |
| `.github/skills/simple-cluster-cavgs-quality/SKILL.md`, `.github/skills/simple-main-nano/SKILL.md` | 5 |
| `src/main/cavg_quality/simple_cavg_quality_{types,feats,model,analysis,learn,selection,relations}.f90` | 2, 4 |
| `src/main/cavg_quality/simple_cavg_quality_lattice.f90` (new, from `src/main/frame_quality/`), its tester | 2 |
| `src/main/frame_quality/` (deleted) | 2 |
| `src/main/commanders/simple/simple_commanders_cavgs.f90` | 2 |
| `src/main/ui/simple_ui_params_common.f90` (the `quality_context` and `quality_model` choices) | 2, 4 |
| `src/main/params/simple_parameters_phases.f90` (the per-mode checks of `quality_context`) | 2 |
| `src/main/commanders/test/simple_commanders_test_class.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90` (the sub-suite) | 2 |
| `src/main/commanders/single/single_commanders_nano2D.f90` (only if `estimate_diam` is retired by ruling) | 4 |

No other repository file is edited. Scratch holds builds, projects, the study's scripts and tables.

## 13. The clock

Twenty-four hours from launch. The order of value if time runs short: the Phase 1 study (it answers the
question without the Fortran), then Phases 2 and 3, then Phase 4; Phase 5 always runs, with what exists.
A phase that would start less than ninety minutes before the deadline is skipped for the close; a phase
that overruns is cut by the loop and the close reports where it stood.

## 14. Journal of decisions (written by the run)

**Phase 0 (2026-10-28).**

- The commit the worktree was made from (8e0ea79a9) does not build alone: it uses
  `src/main/nano/simple_atomic_species.f90`, `simple_nano_density.f90` and `simple_nano_density_tester.f90`,
  which exist only untracked in the maintainer's checkout. They were copied in verbatim, are not part of this
  work, are never edited and are left out of `full.diff`. With them the Debug fast gate fails one entry,
  `unit_single`, on a bounds error in `simple_nano_density.f90` (`profile`, line 373: `prof` 0:447 against a
  1001-row table, in `test_transfer_kernels`); the Release fast gate passes 15 of 15. That failure is the
  baseline the later phases are held to.
- Ruled for the run, as the brief takes them: the lattice module is
  `src/main/cavg_quality/simple_cavg_quality_lattice.f90`; `estimate_diam` is not retired; the scaled cryo-EM
  features stay or go by the rule of section 6.
- The three labeled datasets rebuild exactly by the recipe of section 5 (with `kv`, `cs`, `fraca` given to
  `import_particles`, which requires them even with `ctf=no`, and `nparts=1` to `make_cavgs`): populations and
  states equal class by class, the frames' class, alignment, state, `indstk` and `eo` equal, and even the
  stored `corr` and `res` come back equal.
- The scattered classes (section 5) have a cause in the code. On a fresh 2D start the search sets each frame's
  previous class to zero and then to a random populated class (`prep4srch` in `simple_strategy2D_srch.f90`,
  through `is_fresh_2D_start`), so the time-chunk assignment `init_tseries_refs` writes is never read, and
  `strategy2D_tseries` lets each frame start from a random class and move only to that class's neighbours. The
  previous run's `refine2D_nano` with `nptcls_per_cls=20` on the same trajectories scatters exactly as the
  delivered runs do, and so does this build's `analysis2D_nano` (defaults, `element=Pt`) on the hold-out
  71-6/traj03: 661 non-empty classes of 684, 0.40% of consecutive frames share a class, the frames of a class
  are a median 329 frames apart and span 0.87 of the trajectory. The file is outside the table; this is
  reported, not fixed. The plan goes on with scattered classes, the case it calls the easier one.
- Project-backed `evaluate` fails at HEAD for every relational model ("model requires relational feature"):
  `evaluate_cavg_quality_result` builds its dataset without the relational feature. To be fixed in Phase 2;
  until then the byte gate leaves that run out.
- The byte gate (`~/nano_cavgs_mac/scripts/cavg_gate.sh`): `learn`, `promote` and `evaluate` (the three presets
  and the learned file) on five synthetic tables, since the repository has no class-average fixtures; `analyze`
  in the chunk, pool and sieve contexts on a cryo-EM project with a manual selection (20 classes, 3 deselected)
  and on the rebuilt nano 71-6/traj01. 49 files, reproducible byte for byte. The FFTW wisdom cache is left out
  of the comparison. Every run starts from a pristine copy of the project, because each program appends its
  command line to the `jobproc` segment, and checks that `cls2D`, `ptcl2D` and `out` come back unchanged.
- In the chunk, pool and sieve contexts every nano class is hard rejected, as section 3 expects: the 100-40 A
  band is empty in a 43 A box, and the geometry features are constant.

**Phase 1 (2026-10-28): the study and the bank decision.** The record is `~/nano_cavgs_mac/study/DECISION.md`.

- How it was measured. Step 2 was not a numpy reproduction. It ran the program itself on projects rebuilt at
  `smpd = 0.358 f`: every extraction constant is in A and the mask in pixels, so that is the program with the
  constants divided by f, exactly (pixels verified identical). The relational feature and the full bank of a
  validity-only context came from a scratch copy of the source (outside the repository). There, the `sieve` gates
  were replaced by the nano validity gates, because the cryo-EM contexts hard-reject every nano class, the
  relational feature is computed only on surviving classes, and the population gate has no command-line switch.
  The lattice family and the leave-one-dataset-out fits were done in numpy.
- Lattice family: `lattice_peak` reaches AUC 0.990/0.975/0.983 and `spot_snr` 0.987/0.972/0.974, against the
  stored correlation's 0.965/0.907/0.908.
- The floor guard is changed from 4 to 3 floor sd for class averages. On 71-6/traj01 the deselected classes form
  two modes, so the floor is wide and the separation is 3.7 sd. Under the 4-sd guard that dataset falls back to
  median/MAD normalization on a different scale, and every lattice model fails there (balanced accuracy 0.50
  at AUC 0.975). The weight guard stays at 5%. A separation of 6.7 sd was seen on the hold-out 71-6/traj03.
- Scaled cryo-EM constants: at f = 8 or 16 the local variance and center/edge features become strong (AUC
  0.95-0.97). At f = 8, `centered` and `cc_area_frac` are degenerate; at f = 16 nothing is.
- The relational feature at a nano band (15 A to 1.5 A, 3 px) carries nothing (leave-one-out AUC 0.44).
- Held-out balanced accuracy / AUC of the prototype fit: lattice family + relational 0.978/0.992, 0.933/0.971,
  0.944/0.983. Against it:
  - stored correlation alone at a learned threshold: 0.890/0.965, 0.756/0.907, 0.850/0.908;
  - the stored correlation's best threshold on each project itself: 0.937, 0.880, 0.860;
  - `lattice_peak` alone at a learned threshold: 0.960, 0.915, 0.940;
  - the floor rule at 4 sd: 0.984, 0.640, 0.873.
- Decision by the rule of section 6. Without the lattice family the balanced accuracy drops by 0.08-0.11; without
  any other family (microchunk, signal, overfit at f = 4, 8, 16, the score, the relational feature) it changes
  by less than 0.01.
  - The nano bank is the six lattice features, and the policy is `lattice`.
  - No scaling factor is fixed, because no scaled feature enters.
  - The relational term stays, because the model schema requires it. It uses the nano band and a learned
    coefficient.
  - The other constants of section 3 are kept: the band, the mask at box/2 - 2, the floor mixture, the spot
    threshold at 4 floor sd, the normalization and its minimum scales, and `CLIP_Z`.
- Two defects found, both to be fixed in Phase 2:
  - A model file that names a `preset` without restating its interaction terms is refused (`read_model`,
    "interaction_coefficients count must match"), contrary to the policy document.
  - Project-backed `evaluate` passes no relational feature (Phase 0).

**Phase 2 (2026-10-28): the context and the bank.**

- Files changed: the plan's table for Phase 2 (`simple_cavg_quality_{types,feats,model,analysis,learn}.f90`, the
  commander, `simple_ui_params_common.f90`, `simple_parameters_phases.f90`, the test class and its UI). New:
  `src/main/cavg_quality/simple_cavg_quality_lattice.f90` and `simple_cavg_quality_lattice_tester.f90`.
  Deleted: `src/main/frame_quality/`.
- The bank is indexed by position in the context's own list. The nano list (6) is shorter than the cryo-EM list
  (14), so `CAVG_QUALITY_NFEATS`, the maximum, stays 14. `cavg_quality_bank_size(context)` and the names and
  families by context drive everything that reads the bank: extraction, normalization, both training-table
  writers and the reader, the policies, the model file (feature weights, coefficients, interaction terms by index
  or name), scoring and the promotion snippet. The chunk, pool and sieve code paths, arrays and presets are
  unchanged; the byte gate is identical on all 49 files.
- The lattice module is a class (`cavg_quality_lattice`: `new`, `extract`, `fit_floor`, `kill`) plus the floor
  normalization. Dropped as dead: the band maximum, the running median, the white-noise reference and the
  per-frame table. On the labeled class averages its raw values equal the Phase 1 numpy reference within 1e-4.
- Guard: a third condition. The upper component must lie at least 0.2 above the floor in log10 ratio. On a
  64-class all-lattice synthetic project the mixture carved out a narrow lower component (weight 0.09, sd 0.004)
  that passed both other conditions (gap 0.02). The real projects have gaps of 0.46-0.66. A guard on the larger
  component sd was rejected, because it fails the hold-out 71-6/traj03, whose lattice population is broad.
- Nano `analyze` without a nano model is the hard-gate-only decision. A coefficient-free nano model names the
  context, the bank and the relational band (15 A to 1.5 A, 3 px) in the table. The table header records the
  floor and the median `spot_diam` over the accepted classes. `apply` writes `diam` on every scored class row
  and its median on the class-average entry of `os_out`.
- Defects fixed:
  - project-backed `evaluate` now passes the relational feature;
  - `read_model` keeps a preset's interactions when the file restates none;
  - `evaluate` refuses training files and a model of different banks.
  - My own port bug, found by the tester before any use: an auto-allocated half-plane mask lost its lower
    bound.
- Commander test: a model learned on 71-17/traj04 and 71-6/traj04 was applied with `prune=yes` to a copy of the
  71-6/traj01 project. It accepted 443 of 684 classes. Balanced accuracy against the labels was 0.903 (TPR 0.953,
  TNR 0.852), before Phase 3's leave-one-out runs. The prune kept exactly the selected frames of the accepted
  classes with `indstk` and order intact. `map_cavgs_selection` only deselects, so frames Hans had already
  deselected stayed out.
- Left: the `quality_context` comment in `simple_parameters.f90` still lists three contexts. That file is outside
  the table.

**Phase 3 (2026-10-28): learn and evaluate.** Tables: `~/nano_cavgs_mac/keep/phase_3/`. No source changed.

- Leave one labeled dataset out: `learn` runs on two nano tables and project-backed `evaluate` on the third
  project. The learned models use the `lattice` policy and pick thresholds of 0.20, 0.15 and 0.10.
- Held-out agreement, AUC / balanced accuracy, the model at its threshold against the best baseline:

  | held out | model | best baseline |
  | --- | --- | --- |
  | 71-17/traj04 | 0.991 / 0.916 | floor rule at 4 sd: 0.990 / 0.984 |
  | 71-6/traj01 | 0.973 / 0.903 | `lattice_peak` at a learned threshold: 0.975 / 0.915 |
  | 71-6/traj04 | 0.983 / 0.850 | `lattice_peak` at a learned threshold: 0.983 / 0.940 |

  - The stored correlation at a learned threshold reaches 0.890, 0.756 and 0.850; with its best threshold on
    each project, 0.937, 0.880 and 0.860. The relational feature alone stays near 0.5.
- Decision: the model does not meet the agreement criterion of section 8 on any held-out dataset, so it is not
  promoted on agreement. It ranks the classes exactly as `lattice_peak` alone does (AUC within 0.002). Its
  thresholds come from the learner's objective, specificity under a 2% false-negative tolerance, which was built
  for cryo-EM chunks; they keep recall at 0.95-0.99 and cost specificity (0.71-0.85). At its best threshold on each
  project it reaches 0.985, 0.940 and 0.950, the same as `lattice_peak`. The gate by hand: `lattice_peak` at
  least about 2 floor sd above the floor mean (2.09 pooled over the three sets, in sample, balanced accuracy
  0.941). The labels sit well below the 4-sd floor rule.
- Hold-outs (the model learned on all three, `apply` on copies with `prune=no`, stacks for review). The floor
  guard passes on four of seven. It fails on 71-17/traj01 (2.6 sd), on 71-6/traj02 (2.9 sd), and on 71-17/traj03,
  whose classes show no lattice population (gap 0.24). Where the guard fails, the lattice features fall back to
  median/MAD and the model's decision is not calibrated. Where it passes, the model accepts every class at or above
  the 2.09-sd gate, plus more. The class diameters (median over accepted) run from 20 to 32 A.
- Phase 4 note: the hold-outs have no hand selection, so the 3D comparison of section 8 (model, hand selection,
  floor) is run on a labeled dataset held out of the model's training: 71-6/traj01 with its leave-one-out model.

**Phase 4 (2026-10-28): cancelled by ruling 1.** No 3D experiment. Of the five planned refinements, two (the hand
selection and the leave-one-out model) had started and were stopped; none finished; their partial outputs were
deleted. No `nano` preset is promoted.

**Phase 5 (2026-10-28): close.** `doc/policies/sieving_and_rejection/model_cavgs_rejection.md` gains the nano
context and the gate by hand (`analyze` in context nano, a filter of the table at 2.09 floor sd, `selection` of
`cls2D` then of `ptcl2D state=1 prune=yes`; checked end to end on a copy of 71-6/traj01: 7244 frames, the frames
of the 378 gated classes; `selection oritype=cls2D` alone does not prune). The skills `simple-cluster-cavgs-quality`
and `simple-main-nano` updated. The report is `nano_cavgs_rejection_report.md`; `full.diff` leaves out the three
copied files.
