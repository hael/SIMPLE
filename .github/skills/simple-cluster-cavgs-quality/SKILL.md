---
name: simple-cluster-cavgs-quality
description: Use when working on SIMPLE's model_cavgs_rejection application, the src/main/cavg_quality class-average quality backend (features, hard gates, relational feature, pairwise-logistic models, learn/evaluate/promote), quality_mode/quality_context/quality_model/score_states, cavgs_quality_training.txt analysis outputs, chunk vs pool vs sieve vs nano contexts, the bank-free logistic core in src/main/quality, or the model-based class-average scoring used by the stream particle sieve, stream initial analysis, solve3D state assignment and the lattice gates of trajectory_solve3D.
---

# SIMPLE `model_cavgs_rejection`

Use this skill for the class-average quality application, its backend, and for
analysis of `quality_mode=analyze` outputs.

## Mental Model

`model_cavgs_rejection` is a `simple_exec` class-average processing program. It is not routed through `simple_exec_refine2D`.

Flow:

`simple_ui_cavgproc -> simple_exec_cavgproc -> commander_model_cavgs_rejection -> simple_cavg_quality_*`

The application reads class averages from an input SIMPLE project, applies the
hard gates of a context (`chunk`, `pool`, `sieve` or `nano`), scores the survivors with
a pairwise-logistic model on normalized features plus one relational feature,
writes selected/rejected class-average stacks, and either updates the project
(`quality_mode=apply`) or compares against the current manual `cls2D` state
without modifying the project (`quality_mode=analyze`, `evaluate`).

## Rules That Fail Silently

- Scores are relative to the scored set: every feature is centred and scaled by the median and MAD of
  the set's hard-gate survivors (`normalize_robust`; with one survivor or none every feature is 0), and
  the relational feature is a kNN over the same set. Scoring a subset or a merged set changes every score.
- Model files and built-in presets are positional: `feature_weights` and `linear_coefficients` follow
  the bank order, and a model file must give exactly one value per feature of the bank (`read_model`
  stops otherwise). Training tables (`z_<name>` columns)
  and `interaction_terms` (`a:b`) go by name. Reordering or inserting a feature re-maps the weights of
  every model file and preset without an error; a feature appended at the end gets weight 0 in them.
- `score_project_cavgs` needs `smpd` (the class averages' pixel size) when it scores with a model: the
  relational feature's mask radius comes from it, and without it the default pixel size is used. It
  builds its own `parameters` (`relation_params%new`) at the caller's team and puts back the caller's
  random stream, `nthr_glob`, team, `CWD_GLOB` and worker identity (`part_glob`, `l_distr_worker_glob`, the
  worker's recorded project) after (`get_rnd_state`/`put_rnd_state`).
- `nano` apply writes `diam` and `lattice_peak` into `cls2D` for every class that passes the hard gates
  (and the median diameter beside the class-average entry of `out`). `trajectory_solve3D`
  (`commanders/single/single_commanders_trajectory_solve3D.f90`) runs `quality_model=nano` at
  `thres=LAT_PEAK_MIN_RAW`, reads `lattice_peak` back (0 when absent) and cuts with `model_gate_states`
  and `K_MODEL` of `simple_cavg_quality_lattice.f90`: change these together.
- A standing rule: every change to code the contexts share (the logistic core, the learner, the table
  and model readers), and any `nano` work, leaves the cryo-EM `learn`, `evaluate` and `analyze` outputs
  byte-identical. No test checks this: compare the outputs before and after by hand.
- Do not extend the `nano` context here: `doc/implementation_notes/planned/trajectory_cavgs_rejection.md`
  moves it into its own program.

## Key Files

- `src/main/ui/simple/simple_ui_cavgproc.f90`: UI registration for `model_cavgs_rejection`.
- `src/main/ui/simple_ui_params_common.f90`: `quality_mode`, `quality_context`, `quality_model`, `score_states`.
- `src/main/params/simple_parameters_phases.f90`: per-mode validation of those keys.
- `src/main/exec/simple_exec_cavgproc.f90`: routes `model_cavgs_rejection` to the commander.
- `src/main/commanders/simple/simple_commanders_cavgs.f90`: `exec_model_cavgs_rejection`, output writing, project annotation, `score_states`.
- `src/main/cavg_quality/simple_cavg_quality_feats.f90`: feature bank, extraction, hard gates, normalization.
- `src/main/cavg_quality/simple_cavg_quality_relations.f90`: the relational (CC-neighbour) feature.
- `src/main/cavg_quality/simple_cavg_quality_model.f90`: built-in presets, model file I/O, scoring, promotion snippets.
- `src/main/cavg_quality/simple_cavg_quality_analysis.f90`: evaluation, hard-gate-only evaluation, training and feature tables.
- `src/main/cavg_quality/simple_cavg_quality_learn.f90`: learn/evaluate model fitting and reports.
- `src/main/quality/simple_quality_logistic.f90`: the bank-free core the learner calls for every context (the logistic problem and its L-BFGS-B fit, dataset roles and class weights, the lambda and threshold grids, the robust macro score); `simple_cavg_quality_learn.f90` is its only caller. Order: `new(feature indices, nobs, has_extra, lambda)` -> `set_row` for every row -> `finalize` -> `fit` -> `kill`. A row never set keeps weight 0 and drops out without an error; `fit` without `finalize` stops. It knows no feature names: the indices point into the caller's bank (at most `QUALITY_MAX_FEATS`). Tests: `unit_numerics suite=quality_logistic`.
- `src/main/cavg_quality/simple_cavg_quality_selection.f90`: `score_project_cavgs` and `write_cavg_selection_stacks`, shared by every caller.
- `src/main/cavg_quality/simple_cavg_quality_stats.f90`, `_helpers.f90`, `_types.f90`: metrics, helpers, shared types and constants.
- `src/main/cavg_quality/simple_cavg_quality_lattice.f90`: the `nano` context's lattice bank (Fourier spot
  measurements of a class average) and the lattice grade `LAT_PEAK_MIN`, the preset's threshold and `spot_count`'s bar;
  `LAT_PEAK_MIN_RAW` (0.52, class averages of raw frames) and the model gate (`model_gate_states`, `K_MODEL`) of
  `trajectory_solve3D`. `doc/implementation_notes/planned/trajectory_cavgs_rejection.md` plans to move the `nano`
  context out of `model_cavgs_rejection` into its own program.
- `doc/policies/sieving_and_rejection/model_cavgs_rejection.md`: current modes, contexts, hard gates, feature bank, built-ins, learning and promotion.
- `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`: how the stream particle sieve calls this backend.

## Parameters

- `quality_mode=apply` (default): writes `cavgs_quality_features.txt`, the stacks below, maps selection into particles, annotates `cls2D` (and `cls3D` when its count matches) with `quality`, `quality_cluster`, `accept`, optionally prunes (`prune=yes`), and writes the project.
- `quality_mode=analyze`: requires `quality_context`, treats the existing `cls2D` state as the manual reference, writes `cavgs_quality_training.txt` and the stacks, and leaves the project unchanged.
- `quality_mode=learn`: fits a model from the `cavgs_quality_training.txt` files listed in `filetab` (all of one context); writes the model (`fname`) and `cavgs_quality_learn_report.txt`.
- `quality_mode=evaluate`: scores a fixed model on `filetab` training files, or on one project against its manual `cls2D` state; writes `cavgs_quality_evaluate_report.txt`.
- `quality_mode=promote`: emits a Fortran snippet for adding a model file (`infile`) as a built-in preset.
- `quality_context=chunk|pool|sieve|nano`: the hard gates and the feature bank; accepted only by `analyze`. Every other mode takes the context from the model artifact.

## The bank by context

`chunk`, `pool` and `sieve` share the fourteen-feature bank. `nano` has its own six lattice features
(`lattice_peak`, `spot_count`, `spot_arc`, `spot_diam`, `spot_snr`, `band_presence`) and the validity gates only;
`lattice_peak` enters raw (a log10 peak-over-ring ratio with an absolute threshold, `LAT_PEAK_MIN`), the rest by
median/MAD like the cryo-EM bank. Every bank array is indexed by position in the context's list: use
`cavg_quality_bank_size(context)` and `cavg_quality_feature_name(i, context)`, never `CAVG_QUALITY_NFEATS` (the largest
bank) as a loop bound.
Without a nano model, `analyze quality_context=nano` is the hard-gate-only decision.
- `quality_model`: built-ins `chunk100mics` (default), `sieve`, `pool`, `nano` (`CAVG_QUALITY_MODEL_*_DEFAULT` in `simple_cavg_quality_model.f90`; `nano` is the lattice gate `lattice_peak >= LAT_PEAK_MIN` (0.65, set for recall by ruling), `thres=` moves it for a run; `quality_model=nano` is how a nanoparticle time series is rejected). `infile` is a complete external model and wins; `fname` is output-only.
- `score_states=yes`: diagnostic, `apply` only, no overrides or pruning: scores every built-in model by `cls3D` state into `cavgs_quality_state_scores.txt` and changes nothing.

## Analyzing a Run

The `quality_mode=analyze` outputs, the analysis script (`scripts/analyze_quality_run.py`, Python 3.7
or later) and the interpretation checklist are in `references/analysis.md`.

## Other Callers

The library callers score through `score_project_cavgs` and keep their own policy (the commander calls
`evaluate_cavg_quality` with its own `params`):

- the stream particle sieve (`src/main/sieve/simple_ptcl_sieve.f90`): the
  `sieve` preset when `use_model` is on, else the `sieve` hard gates only;
  then its coarse/fine `class_compatibility` models. Preserve the sieve chunk
  lifecycle and sentinels described in
  `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` and the skill
  `simple-main-sieve`;
- the stream's initial analysis: the `chunk100mics` preset;
- `score_cavgs_with_model` in `src/main/solve/simple_state_assignment.f90`,
  used by `solve3D_pyramid` and the stream's solve3D stage: the `pool` preset;
- `trajectory_solve3D` (SINGLE) runs the program with `quality_model=nano`, then its own model gate
  (see the rules above).

Prefer small local changes:

- share or extract project annotation/state-mapping code before changing rejection policy;
- avoid exposing calibration margins as casual CLI knobs;
- a new built-in needs a learned artifact, `quality_mode=promote`, and the UI choices of `quality_model`.

## Tests

- `build/bin/simple_test_exec test=unit_numerics` (about 1 s; `suite=` `cavg_quality_relations`,
  `cavg_quality_selection`, `cavg_quality_nano` (team against one-thread identity included), `quality_logistic`).
- No test runs the commander: check a behaviour change with an `analyze` run on a manually selected project,
  and the byte-identity rule above with saved outputs.
- Callers: `test=stream_chains suite=sieve_to_3d` (the sieve), `single_workflow_fcc` (the nano gates of
  `trajectory_solve3D`; its gates may lose at most 5% of the frames).
