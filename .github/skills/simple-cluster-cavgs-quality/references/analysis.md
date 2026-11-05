# Analyzing a `quality_mode=analyze` run

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
python3 .github/skills/simple-cluster-cavgs-quality/scripts/analyze_quality_run.py /path/to/run  # Python 3.7 or later
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
