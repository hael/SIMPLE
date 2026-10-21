---
name: simple-cluster-cavgs-quality
description: Use when working on SIMPLE's model_cavgs_rejection application, the src/main/cavg_quality class-average quality backend (features, hard gates, relational feature, pairwise-logistic models, learn/evaluate/promote), quality_mode/quality_context/quality_model/score_states, cavgs_quality_training.txt analysis outputs, chunk vs pool vs sieve contexts, or the model-based class-average scoring used by the stream particle sieve, stream initial analysis and solve3D state assignment.
---

# SIMPLE `model_cavgs_rejection`

Use this skill for the class-average quality application, its backend, and for
analysis of `quality_mode=analyze` outputs.

## Mental Model

`model_cavgs_rejection` is a `simple_exec` class-average processing program. It is not routed through `simple_exec_refine2D`.

Flow:

`simple_ui_cavgproc -> simple_exec_cavgproc -> commander_model_cavgs_rejection -> simple_cavg_quality_*`

The application reads class averages from an input SIMPLE project, applies the
hard gates of a context (`chunk`, `pool` or `sieve`), scores the survivors with
a pairwise-logistic model on normalized features plus one relational feature,
writes selected/rejected class-average stacks, and either updates the project
(`quality_mode=apply`) or compares against the current manual `cls2D` state
without modifying the project (`quality_mode=analyze`, `evaluate`).

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
- `src/main/cavg_quality/simple_cavg_quality_selection.f90`: `score_project_cavgs` and `write_cavg_selection_stacks`, shared by every caller.
- `src/main/cavg_quality/simple_cavg_quality_stats.f90`, `_helpers.f90`, `_types.f90`: metrics, helpers, shared types and constants.
- `doc/policies/sieving_and_rejection/model_cavgs_rejection.md`: current modes, contexts, hard gates, feature bank, built-ins, learning and promotion.
- `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`: how the stream particle sieve calls this backend.

## Parameters

- `quality_mode=apply` (default): writes `cavgs_quality_features.txt`, the stacks below, maps selection into particles, annotates `cls2D` (and `cls3D` when its count matches) with `quality`, `quality_cluster`, `accept`, optionally prunes (`prune=yes`), and writes the project.
- `quality_mode=analyze`: requires `quality_context`, treats the existing `cls2D` state as the manual reference, writes `cavgs_quality_training.txt` and the stacks, and leaves the project unchanged.
- `quality_mode=learn`: fits a model from the `cavgs_quality_training.txt` files listed in `filetab` (all of one context); writes the model (`fname`) and `cavgs_quality_learn_report.txt`.
- `quality_mode=evaluate`: scores a fixed model on `filetab` training files, or on one project against its manual `cls2D` state; writes `cavgs_quality_evaluate_report.txt`.
- `quality_mode=promote`: emits a Fortran snippet for adding a model file (`infile`) as a built-in preset.
- `quality_context=chunk|pool|sieve`: the hard gates; accepted only by `analyze`. Every other mode takes the context from the model artifact.
- `quality_model`: built-ins `chunk100mics` (default), `sieve`, `pool` (`CAVG_QUALITY_MODEL_*_DEFAULT` in `simple_cavg_quality_model.f90`). `infile` is a complete external model and wins; `fname` is output-only.
- `score_states=yes`: diagnostic, `apply` only, no overrides or pruning: scores every built-in model by `cls3D` state into `cavgs_quality_state_scores.txt` and changes nothing.

## Analyze Outputs

Files from `quality_mode=analyze`:

- `cavgs_quality_training.txt`: the canonical learner input. `# key=value`
  comments (version, dataset, model, `quality_context`, relational schema,
  counts, `soft_decision`, confusion metrics, AUC, current and best
  thresholds, model fields), `# feature_summary,` rows (per-feature AUC,
  medians, robust separation, current and suggested weight), `# threshold_scan,`
  rows, then one CSV row per class (state, hard reject, score, manual state,
  raw and `z_` features, the relational feature).
- `quality_selected_cavgs.mrc`, `quality_rejected_cavgs.mrc`, `hard_gate_rejections.mrc`, `quality_ranked_cavgs.mrc` and `quality_ranked_cavgs.txt`: stacks for visual inspection.

When analyzing a run directory, first run:

```bash
python3 .github/skills/simple-cluster-cavgs-quality/scripts/analyze_quality_run.py /path/to/run
```

Use `--recursive` on a parent directory to summarize many run folders:

```bash
python3 .github/skills/simple-cluster-cavgs-quality/scripts/analyze_quality_run.py --recursive /path/to/test-data
```

Then inspect the reported false positives, false negatives, best thresholds, weak/inverted features, and hard-rejected manual-good classes. For chunk tuning, pay special attention to whether manual-good classes are lost through the context's hard gates or through the model's probability threshold.

## Interpretation Checklist

1. Confirm `quality_mode=analyze` was used and which `quality_context` it ran. In analyze mode, project selection should be unchanged.
2. Compare current threshold metrics with `best_balacc_threshold` and `best_f1_threshold`.
3. Separate false negatives caused by `hard_reject=T` from false negatives caused by score thresholding.
4. Rank features by AUC and robust separation. A high current weight with poor or inverted AUC is a tuning warning.
5. Check whether false positives are accepted by population/resolution alone or by local signal and histogram support.
6. Compare contexts (`chunk`, `pool`, `sieve`) only on the same project and manual states.
7. Use the MRC stacks when numeric diagnostics disagree with visual quality.

## Other Callers

All callers score through `score_project_cavgs` and keep their own policy:

- the stream particle sieve (`src/main/sieve/simple_ptcl_sieve.f90`): the
  `sieve` preset when `use_model` is on, else the `sieve` hard gates only;
  then its coarse/fine `class_compatibility` models. Preserve the sieve chunk
  lifecycle and sentinels described in
  `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` and the Particle
  Sieving section of `.github/skills/simple-main-stream/SKILL.md`;
- the stream's initial analysis: the `chunk100mics` preset;
- `score_cavgs_with_model` in `src/main/solve/simple_state_assignment.f90`,
  used by `solve3D_pyramid` and the stream's solve3D stage: the `pool` preset.

Prefer small local changes:

- share or extract project annotation/state-mapping code before changing rejection policy;
- avoid exposing calibration margins as casual CLI knobs;
- a new built-in needs a learned artifact, `quality_mode=promote`, and the UI choices of `quality_model`.

## Tests

`simple_test_exec test=unit_numerics` (sub-suite `cavg quality relations`);
for the commander, an `analyze` run on a manually selected project.
