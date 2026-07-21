# Memory model maintenance

This directory contains the periodic maintenance tools for the SIMPLE memory
estimator. Normal users only need `scripts/memory_estimator.py`; the
active runtime coefficients are in `scripts/memory_estimator_models.json`.

## Collect measurements

The benchmark programs create isolated inputs, run the selected commander with
native SIMPLE memory telemetry, and write `results.csv` plus metadata and logs.

```bash
python3 scripts/memory/benchmark_motion_correct.py --help
python3 scripts/memory/benchmark_abinitio2d.py --help
python3 scripts/memory/benchmark_abinitio3d.py --help
```

Keep completed benchmark output under `output/`. Dataset paths are passed
explicitly to the fitter so a periodic refit cannot accidentally consume a
stale or partial run.

## Fit and review

```bash
python3 scripts/memory/fit_models.py \
  --motion-csv output/motion_grid/results.csv output/motion_followup/results.csv \
  --abinitio2d-csv output/abinitio2d_grid/results.csv \
  --abinitio3d-csv output/abinitio3d_grid/results.csv
```

This does not change the active model. It writes:

- `output/memory_estimator_models_fitted.json`
- `output/memory_estimator_validation.txt`

Review the validation report, calibration ranges, allocation coverage, and
outliers. Test the candidate explicitly with:

```bash
python3 scripts/memory_estimator.py \
  --models output/memory_estimator_models_fitted.json \
  motion_correct --xdim 4096 --ydim 4096 --frames 16 --smpd 1.3 --threads 4
```

## Promote a reviewed model

After review, refit and install the result used by the runtime estimator:

```bash
python3 scripts/memory/fit_models.py \
  --motion-csv output/motion_grid/results.csv output/motion_followup/results.csv \
  --abinitio2d-csv output/abinitio2d_grid/results.csv \
  --abinitio3d-csv output/abinitio3d_grid/results.csv \
  --install
```

After installation, smoke-test all three runtime entry points:

```bash
python3 scripts/memory_estimator.py motion_correct --help
python3 scripts/memory_estimator.py abinitio2D --help
python3 scripts/memory_estimator.py abinitio3D --help
```

Installing a fit changes only the standalone Python estimator. It does not
alter SIMPLE queue requests or require recompiling SIMPLE.
