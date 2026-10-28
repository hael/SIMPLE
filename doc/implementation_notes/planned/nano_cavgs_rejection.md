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

**Hard gates.** The class averages of `refine2D_nano` start as time chunks of `nptcls_per_cls` frames, so
population carries no quality information, and the FRC "resolution" of a nanoparticle class is a statement
about the lattice band, not a quality gate. `nano` keeps the validity gates only: invalid pixels; an empty
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

A labeled dataset is one `analysis2D_nano` run with the maintainer's selection: the project file (its
`cls2D` with population, resolution, correlation, the class averages registered in `os_out`, the frame to
class map in `ptcl2D`, the alignment) and the class-average stack, with the labels either as the `cls2D`
`state` of the project or as a text file of `class, state` rows that the run applies to a copy of the
project before `analyze` (the `oris` setters; the project is never modified in place). The maintainer is
writing out the class labels for the trajectories he selected (2026-10-28); their location is given at
launch in the run's prompt. Labeled datasets train and are evaluated leave-one-out; every other trajectory
of `~/Processing/NanoX/NanoX_trajs` is a hold-out whose `analysis2D_nano` output the run makes itself if
absent and scores for review.

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
   family, against two baselines: the lattice floor alone (`lattice_peak` of the class at 4 floor sd) and
   the relational feature alone.

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
the maintainer's labels (balanced accuracy at the model's threshold, and AUC) is not below the better of the
two baselines of section 6, step 4; and, on at least one held-out trajectory, `autorefine3D_nano` from its
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

(empty)
