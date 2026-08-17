# Continuous 3-D PCG pose production validation

This directory contains the complete Oracle Linux validation runner and its
checker. The runner does not change the source checkpoint. It copies the
checkpoint into five audit arms and writes all evidence under one timestamped
directory.

## Input checkpoint

Use a completed, single-state beta-gal project. It must contain active `ptcl3D`
particles, fixed even/odd ownership, particle images, and CTF metadata. No
refinement/search iteration is run. The checker supports `pgrp=c1` and
`pgrp=d2`.

For the frozen simulated beta-gal project, use the PCG reconstruction directory
that contains `betagal_frozen.simple`. Each audit arm reconstructs its own base
PCG half-maps from the frozen project before the optional polish.

## Run

Build SIMPLE first. Confirm that `simple_exec`, `simple_private_exec`, and
`simple_test_continuous_3D_pcg_refinement` are on `PATH`. The runner loads
`miniconda/25.9.1` automatically if the current `python3` is too old and the
Oracle module command is available.

```bash
cd ~/Projects/hael_SIMPLE-rsync-test

bash production/tests/continuous_3D_pcg_pose_validation/run_oracle_validation.sh \
  --projfile "$HOME/Projects/frozen_noisy_betagal_20260813_153030/reconstruction_comparison/pcg/betagal_frozen.simple" \
  --box 144 \
  --smpd 1.3 \
  --mskdiam 120 \
  --pgrp c1 \
  --objfun cc \
  --nthr 8 \
  --distributed-parts 2 \
  --output-parent "$HOME/Projects"
```

The command creates one directory named
`continuous_3D_pose_validation_YYYYMMDD_HHMMSS`. Download only that directory.
After the initial argument, executable, and checkpoint preflight, the runner
does not stop at the first failed test. It records the failure and continues
with every independent policy case and beta-gal arm. If one arm cannot finish,
its dependent exports are skipped, but the remaining arms still run.

## Method

The runner first executes `case=pose_contract` and the complete mother suite.
It then confirms that four invalid command combinations stop for the expected
reason without changing their copied projects. Finally, it creates one default
audit arm and four matched A/B arms:

- `shared_default`: `nparts=1`, option omitted;
- `shared_off`: `nparts=1`, `pcg_pose_polish=no`;
- `shared_on`: `nparts=1`, `pcg_pose_polish=yes`;
- `distributed_off`: `nparts>1`, `pcg_pose_polish=no`;
- `distributed_on`: `nparts>1`, `pcg_pose_polish=yes`.

Within each pair, `pcg_pose_polish` is the only planned difference. Every arm
runs deterministic fixed-pose `reconstruct3D` with `rec_backend=pcg`,
`combine_eo=no`, `ml_reg=no`, two fixed PCG iterations, and no orientation
search. An enabled arm then polishes against those base half-maps and runs one
final fixed PCG reconstruction. Use `--objfun cc` for the frozen synthetic
project. The runner also accepts `--objfun euclid`; the selected value is
recorded in `run_config.env`.

## Acceptance conditions

- All focused and mother tests pass.
- Every invalid route fails for its specified reason.
- Every invalid route leaves its copied project byte-for-byte unchanged.
- The omitted option produces the same poses, metadata, and FSC as explicit
  `pcg_pose_polish=no`.
- Disabled arms contain no pose-polish marker.
- Enabled arms contain one complete terminal summary and exactly one final PCG
  reconstruction marker.
- Terminal counts balance and invalid numerical outcomes are zero.
- At least one particle has an accepted, measurable pose change.
- Accepted poses change at least one final even/odd half-map artifact.
- State, half ownership, defocus, and astigmatism metadata remain identical
  between each matched A/B pair.
- Angular and shift changes remain inside the cumulative eight-step LM bounds.
- Enabled FSC area is not lower than the matched control by more than `0.01`.
- Shared and distributed pose RMS differences remain below `0.002` radians and
  `0.02` pixels, and their FSC areas differ by no more than `0.01`.

These tolerances are fixed before the run. Do not change them after seeing the
result. A failure can identify an implementation defect, a non-matched A/B
checkpoint, or a scientific limitation. Inspect the evidence before deciding.

## Output

Read these files first:

- `STATUS.txt`;
- `validation.log`, the complete chronological console log;
- `analysis/summary.md`;
- `analysis/summary.json`;
- `analysis/analyzer.log`;
- `failures.txt` when any scheduled check fails;
- `MANIFEST.sha256`.

The `arms/` directory contains complete arm projects, logs, metadata exports,
FSC tables, and volumes. The `policy/` and `unit/` directories contain the
focused contract evidence.

## Truth-diagnostic matrix after an A/B failure

Use `run_truth_diagnostic.sh` when the production A/B runner reports that the
enabled reconstruction has worse FSC. This test does not repeat the
shared/distributed contract. It asks why an accepted pose update can reduce the
internal particle objective but reduce truth agreement or half-map agreement.

The runner uses the frozen truth volume, truth orientation/CTF table, and noisy
stack. It also generates a no-noise stack with SIMPLE's standard projector and
a deterministic perturbed orientation table. It runs matched polish-off and
polish-on arms for:

- clean exact poses and clean perturbed poses at the full band;
- noisy exact poses and noisy perturbed poses at the full band;
- noisy exact and perturbed poses at the FSC=0.5 cutoff;
- noisy exact and perturbed poses at the FSC=0.143 cutoff.

The clean cases test operator consistency and controlled local recovery. The
noisy exact cases test overfitting. The noisy perturbed cases test whether the
method recovers useful pose signal. The two low-pass settings test whether the
problem is caused by fitting frequencies that do not have reliable signal.

Use the FSC resolutions measured in the frozen disabled control. For the
downloaded beta-gal result discussed in the plan, run:

```bash
cd ~/Projects/hael_SIMPLE-rsync-test

bash production/tests/continuous_3D_pcg_pose_validation/run_truth_diagnostic.sh \
  --truth-volume "$HOME/Projects/frozen_noisy_betagal_20260813_153030/frozen_fixture/truth_1JYX_box144.mrc" \
  --truth-oris "$HOME/Projects/frozen_noisy_betagal_20260813_153030/frozen_fixture/betagal_truth_oris.txt" \
  --noisy-stack "$HOME/Projects/frozen_noisy_betagal_20260813_153030/frozen_fixture/betagal_noisy.mrcs" \
  --box 144 \
  --smpd 1.3 \
  --mskdiam 120 \
  --pgrp c1 \
  --kv 300 \
  --cs 2.7 \
  --fraca 0.1 \
  --lp-fsc05 4.255 \
  --lp-fsc0143 3.671 \
  --nthr 8 \
  --output-parent "$HOME/Projects"
```

The command continues through all 16 reconstruction arms after preflight and
creates one directory named
`continuous_3D_pose_truth_diagnostic_YYYYMMDD_HHMMSS`. Download only that
directory. Read `analysis/truth_matrix.md` first.

The microscope values must match the truth orientation table. The runner now
checks that match, supplies the values required by `import_particles`, and
requires the imported project to contain all particles in two nonempty fixed
half-sets before it schedules a reconstruction arm. If the truth table has no
split, it assigns the same deterministic alternating `eo=0/1` ownership used
by `reconstruct3D` and preserves that ownership in every exact/perturbed and
off/on arm.

The feature gate requires both clean cases to pass and at least one noisy band
to pass both exact-pose stability and perturbed-pose recovery. Exact poses can
increase by at most `0.001745` rad rotation RMS and `0.05` pixel shift RMS.
Clean perturbed errors must decrease by at least 50 percent; noisy perturbed
errors must decrease by at least 10 percent. No evaluated FSC area can decline
by more than `0.005`. These limits are fixed before the run.

A clean-exact failure points to forward-operator or objective inconsistency. A
clean-recovery failure points to the derivative, pose convention, scaling, or
solver. Stable clean cases with a noisy-exact failure point to noise
overfitting. If a low-pass pair passes, a data-supported polishing band is a
candidate correction. Do not change the tolerances after the result.
