# SIMPLE memory estimator and benchmark guide

The user-facing estimator is `scripts/memory_estimator.py`; its active
coefficients are stored beside it in `scripts/memory_estimator_models.json`.
Periodic data collection, fitting, and reporting tools live separately under
`scripts/memory/`. The calibrated estimator covers `motion_correct`,
`abinitio2D`, and `abinitio3D`. It reports both a fitted peak-RSS value and a
larger recommended allocation.

## Estimate memory from the command line

Run the estimator from the repository root:

```bash
python3 scripts/memory_estimator.py motion_correct \
  --xdim 4096 --ydim 4096 --frames 16 --smpd 1.3 --threads 4

python3 scripts/memory_estimator.py abinitio2D \
  --particles 2000 --box 192 --threads 4 --references 32

python3 scripts/memory_estimator.py abinitio3D \
  --particles 500 --box 128 --smpd 1.3 --mask-diameter 60 \
  --threads 2 --partitions 4 --states 1
```

Put `--json` before the commander name for machine-readable output. A
`sampled-particles` value of zero means all particles. Inputs outside the
calibration range are accepted but produce explicit extrapolation warnings.

### Model targets

- `motion_correct` estimates the peak RSS of one distributed worker. This is
  the same boundary used by `job_memory_per_task`.
- `abinitio2D` estimates the complete commander peak, not one `cluster2D`
  partition.
- `abinitio3D` estimates the conservative process-tree bound used by the
  benchmark: parent peak plus the largest `nparts` child peaks. It is not an
  assertion that all those peaks occurred at exactly the same instant.

The estimator is intentionally kept outside the SIMPLE executables. It can be
updated and evaluated without rebuilding SIMPLE or silently changing existing
queue requests.

### Calibration and validation

The current models use:

- 200 successful motion-correction runs, including the original 1K-10K movie
  grid and a 36-run follow-up varying 4/16 frames and 1/4 threads;
- 30 successful ab-initio 2D screening runs;
- 36 successful ab-initio 3D screening runs.

Run the reproducible fitter with a Python runtime that provides NumPy:

```bash
python3 scripts/memory/fit_models.py \
  --motion-csv output/motion_grid/results.csv output/motion_followup/results.csv \
  --abinitio2d-csv output/abinitio2d_grid/results.csv \
  --abinitio3d-csv output/abinitio3d_grid/results.csv
```

It writes `output/memory_estimator_models_fitted.json` and
`output/memory_estimator_validation.txt` without changing the active model.
After reviewing those files, rerun the same command with `--install` to promote
the new model. The fitted equations are ordinary least squares over physically
meaningful features. The saved JSON model is the Python estimator's source of
truth.

Current raw-model validation:

| Commander | Rows | R² | Leave-one-out R² | Conservative coverage |
|---|---:|---:|---:|---:|
| motion_correct | 200 | 0.970 | 0.965 | 100% |
| abinitio2D | 30 | 0.946 | 0.908 | 100% |
| abinitio3D | 36 | 0.976 | 0.922 | 100% |

Coverage is in-sample after safety margins and rounding; it is not a guarantee
for a different SIMPLE build, operating system, allocator, algorithm, or input
domain. The 2D and 3D models should be recalibrated as more screening and
interaction runs become available.

## Collect benchmark data

Each benchmark runs the measured workflow in a fresh process and saves the
inputs, peak-memory response, metadata, logs, and native SIMPLE telemetry needed
to refit the models. Use repeated measurements and representative production
data before treating a result as a capacity-planning limit.

### `motion_correct`

`scripts/memory/benchmark_motion_correct.py` measures peak resident memory
across movie dimensions and sampling distances (`smpd`). The current baseline
fixes `algorithm=iso`, `mcpatch=no`, and `downscale=yes`; these settings are
written to `metadata.json` with every run.

Each grid point runs in a fresh `simple_private_exec` process because the
operating system's peak-RSS counter cannot be reset within a process. Synthetic
movie creation and project import happen before the measured process starts.

Run a grid from the repository root after building SIMPLE:

```bash
python3 scripts/memory/benchmark_motion_correct.py \
  --dimensions 1024 2048 4096x3072 \
  --smpds 0.8 1.2 2.0 \
  --frames 20 \
  --threads 4 \
  --repeats 3 \
  --output-dir motion_correct_memory_results
```

The harness finds executables in `PATH` first, then in `build/bin` or
`build/production`. Select nonstandard builds explicitly:

```bash
python3 scripts/memory/benchmark_motion_correct.py \
  --simple-exec /path/to/simple_exec \
  --private-exec /path/to/simple_private_exec \
  --output-dir motion_correct_memory_results
```

The output directory must not already exist. It contains:

- `results.csv`: one summary row per `(dimensions, smpd, repeat)` case;
- `metadata.json`: the complete grid, executable paths, platform, and
  measurement definition;
- `cases/<case>/prepare.log`: synthetic-project setup output;
- `cases/<case>/measure/motion_correct.log`: commander output;
- `cases/<case>/measure/memory_usage_<pid>.csv`: native periodic SIMPLE
  telemetry.

Generated movies and corrected-image products can be large. Add
`--discard-case-data` to retain the CSV, metadata, logs, and native telemetry
while removing those inputs and products after each completed case.

`peak_delta_rss_bytes` is the main response variable: process peak RSS minus
current RSS at telemetry start. `peak_rss_bytes` is also retained for capacity
planning. Failed cases remain in the CSV and make the harness exit nonzero, so
a partially completed grid is not lost.

By default, `motion_correct` downsamples movies finer than 1.3 A/pixel. The
harness fixes and records that target with `--smpd-downscale 1.3`. Consequently,
the CSV includes both input dimensions and `effective_xdim`/`effective_ydim`;
the latter often explain memory use more accurately.

Motion correction currently produces 512-pixel GUI power spectra. The harness
rejects combinations whose effective width or height would fall below 512,
because those inputs cannot complete the commander. Increase the input movie
dimensions or use a smaller `--smpd-downscale` target. For rectangular studies
that intentionally cross this boundary, add `--skip-invalid` to record those
combinations as `unsupported_effective_dimension` while running valid cases.

For production estimates, fit peak or delta RSS against
`effective_pixels_per_frame` while holding frame and thread counts fixed. Both
counts are recorded so they can also be modeled explicitly.

### `abinitio2D`

`scripts/memory/benchmark_abinitio2d.py` measures peak resident memory for
isolated `abinitio2D` runs using deterministic synthetic particle stacks.

The default screening design varies:

- particle count: 100, 250, 500, 1000, and 2000;
- particle box size: 64, 96, 128, 160, and 192 pixels;
- OpenMP thread count: 1, 2, 4, and 8;
- class/reference count (`ncls`): 2, 4, 8, 16, and 32;
- automatic downscaling in selected comparisons with `autoscale=yes`;
- jointly varied interaction cases covering the usable range.

The controlled baseline is 500 particles, box 96, 4 threads, 8 references,
`autoscale=no`, and all particles sampled. Box-size sweeps use `autoscale=no` so
downscaling does not hide the requested box's allocation effect. Separate
autoscaling cases record the effective box parsed from the SIMPLE log.

The screening run uses one bounded initialization stage:

```text
nstages=1 nits_per_stage=1 extr_lim=4 eo_stage=no
refine=prob_snhc sigma_est=global ctf=no smpd=1.3 mskdiam=60
```

SIMPLE stage policy executes the first stage with sampled SNHC even though the
requested overall refinement mode is `prob_snhc`. The short workflow compares
peak allocations; it does not measure scientific convergence or production
runtime.

`nsample` is fixed to all particles in this bounded design. A sampled one-stage
run invokes the full-assignment coverage guard and is not a bounded proxy for a
normal multi-stage sampled workflow. The harness exposes `--include-nsample`,
but requires at least three stages when it is selected.

Run the screening design with:

```bash
python3 scripts/memory/benchmark_abinitio2d.py \
  --discard-case-data \
  --discard-inputs \
  --output-dir /path/to/abinitio2d-memory-results
```

Use `--dry-run` to print the case matrix. Use `--resume` with the same output
directory to skip case IDs already present in `results.csv`. The explicit
`--design factorial` mode creates the full Cartesian product of the four core
variable lists. `--max-cases` prevents accidental submission of an unexpectedly
large grid.

Every case runs in a fresh `simple_exec` process with native SIMPLE memory
telemetry enabled. The primary response is:

```text
peak_delta_rss_mib = peak_rss - current_rss_at_telemetry_start
```

Project creation and particle import happen before the measured process. The
CSV also records total peak RSS, elapsed time, effective box, all input
variables, normalized bytes per input pixel, status, logs, and telemetry paths.

Synthetic data and one repeat are suitable for allocation screening. For
production capacity planning, repeat selected configurations with real data,
the intended CTF mode, the full stage schedule, and the deployed thread count.

### `abinitio3D`

`scripts/memory/benchmark_abinitio3d.py` measures memory across a screened set
of data and workflow variables:

- particle count, raw box size, pixel size, and mask diameter;
- OpenMP threads and local job partitions;
- state count and independent multi-volume mode;
- sampled particles, symmetry, and starting symmetry;
- low-pass limit, projection reconstruction, and initialization route;
- the observed number of refinement iterations;
- four combined interaction cases.

The default design has 33 cases. It is a one-factor screening design with a
small interaction set, not a full Cartesian product.

`abinitio3D` launches `simple_private_exec` workers. Every measured process uses
SIMPLE native memory telemetry. The harness reports:

- peak RSS of the parent `simple_exec` process;
- a conservative worker component calculated as the largest `nparts` worker
  peaks;
- a process-tree planning bound equal to their sum.

The planning bound is suitable for capacity planning but may exceed exact
simultaneous RSS because the component peaks can occur at different times.
Parent and worker components are retained separately in the raw data.

Input generation, import, and the required one-iteration `abinitio2D` class
preparation occur outside the measured interval. Synthetic particles are
realistic projections of a deterministic asymmetric 3D volume with CTF
disabled and SNR 0.1.

The benchmark does not add or require a diagnostic SIMPLE input. Stage 1 uses
the unmodified built-in schedule of up to 20 iterations, and the CSV records
the number of iterations actually observed before completion or early stopping.

Run it from a built SIMPLE checkout:

```bash
python3 scripts/memory/benchmark_abinitio3d.py \
  --simple-exec build/production/simple_exec \
  --output-dir output/abinitio3d_memory_screening_raw \
  --discard-case-data
```

The harness supports `--resume`, `--repeats`, `--dry-run`, and input/case-data
retention controls. Run `--help` for the complete interface.

Generate reports with a Python runtime that provides the required document
libraries:

```bash
/path/to/python3 scripts/memory/report_abinitio3d.py \
  output/abinitio3d_memory_screening_raw/results.csv \
  --metadata output/abinitio3d_memory_screening_raw/metadata.json
```

The report generator writes the complete CSV, a tab-separated TXT file for
future modeling, and a PDF report to `output/pdf/`.
