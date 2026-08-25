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
  --nthr 64 \
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
  --nthr 64 \
  --output-parent "$HOME/Projects"
```

The command continues through all 16 reconstruction arms after preflight and
creates one directory named
`continuous_3D_pose_truth_diagnostic_YYYYMMDD_HHMMSS`. Download only that
directory. Read `analysis/truth_matrix.md` first.

With the particle-parallel polisher, each enabled arm must report
`REQUESTED 64 ACTIVE 64 BATCH 500` before pose work. The requested thread
count controls independent particle LM workers; it does not create 500
threads. The 500-particle batch matches reconstruct3D's bounded work unit and
keeps memory independent of the total project size. Progress is reported every
five complete batches and at the end of each half. Reconstruction setup,
particle I/O, summary reduction, and small final batches can still use fewer
cores. A run started before rebuilding cannot acquire a later parallel or batch
policy while it is running.

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

A clean-exact failure shows that the truth generator and fixed-map objective do
not have the same stationary pose. The matrix alone cannot tell whether the
cause is the forward operator or bias in the reconstructed reference map. A
clean-recovery failure can also involve the derivative, pose convention,
scaling, or solver. Stable clean cases with a noisy-exact failure point to noise
overfitting. If a low-pass pair passes, a data-supported polishing band is a
candidate correction. Do not change the tolerances after the result.

## Fixed-reference isolation diagnostic

Run this focused diagnostic after the truth matrix reports clean exact-pose
drift:

```bash
simple_test_continuous_3D_pcg_refinement case=fixed_reference \
  2>&1 | tee continuous_3D_pcg_fixed_reference.log
```

The test runs eleven collect-all arms. Three observation generators separate the
numerical conventions:

1. the earlier native-box projector helper;
2. the actual noiseless `simulate_particles` geometry: padded volume, padded
   projection, inverse FFT, real-space clip, and native FFT;
3. an independently evaluated padded exact-KB Fourier slice without clipping.

Those observations are compared with raw truth, inverse-envelope-corrected
truth, globally amplitude-fitted raw and corrected truth, a PCG-matched
inverse-crime control, and both raw and corrected two-iteration PCG
reconstructions. The amplitude-fitted arms change only one real volume scale;
they do not fit a pose, shell, or per-particle gain.

Before LM, every arm reports the exact-pose objective, raw and optimally scaled
relative residual, amplitude scale, rotation/shift gradient RMS, and fitted
low/mid/high-shell residuals and correlations. It then reports final pose RMS,
objective change, and all terminal outcomes. The test prints a diagnosis ladder
in this order: LM/Jacobian, exact-versus-fast padded gather, raw versus
inverse-envelope simulator model, simulator padding/clipping, then
reconstructed-reference bias.

Only numerical integrity, complete terminal accounting, and stationarity of the
PCG-matched control are hard gates. Open scientific hypotheses do not stop the
matrix early. A normal exit means `EVIDENCE COMPLETE`, not that production pose
polishing is scientifically accepted.

The shared half-set fixture now uses the same padded-and-clipped geometry. After
the focused matrix, run the mother suite once to recheck half-set FSC,
reconstruction trajectories, noise isolation, and all derivative/recovery
contracts with the corrected observations:

```bash
simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_fixed_reference.log
```

This diagnostic is available by explicit `case=` only. It is not in the mother
suite while the scientific failure is under investigation.

## Forward-path boundary diagnostic

After the eleven-arm matrix identifies a simulator-versus-PCG mismatch, run:

```bash
simple_test_continuous_3D_pcg_refinement case=forward_path \
  2>&1 | tee continuous_3D_pcg_forward_path.log
```

This diagnostic does not run LM. It holds the raw truth volume, 48 exact
orientations, shell range, and PCG workspace fixed while capturing three
versions of each noiseless observation:

1. native frequencies read directly from the padded projected Fourier plane;
2. the same frequencies after `simimg`'s inverse/forward Fourier round trip;
3. the native FFT after the padded real-space image is clipped to the particle
   box.

Every stage is compared with the raw-truth PCG gather. The two adjacent stage
transitions are also compared directly after fitting one global real amplitude.
The report includes total and low/mid/high-shell residuals and correlations,
amplitude scale, and exact-pose rotation/shift gradient RMS.

The predeclared diagnosis tolerance is a fitted relative residual of `1e-3`.
The diagnosis order is padded projector versus PCG gather, Fourier round trip,
then real-space clip/native FFT. Only finite evidence is a hard gate. A normal
exit means `EVIDENCE COMPLETE`, not that a production correction is approved.

The case is diagnostic-only and is not scheduled by the mother suite. After it
completes, run the mother suite once to verify that adding the selector did not
change the ten established cases:

```bash
simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_forward_path.log
```

## Matched-window diagnostic

After `case=forward_path` identifies `REALSPACE_CLIP_NATIVE_FFT`, run:

```bash
simple_test_continuous_3D_pcg_refinement case=matched_window \
  2>&1 | tee continuous_3D_pcg_matched_window.log

simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_matched_window.log
```

The focused case is an explicit correctness prototype. It samples the PCG
volume and its five pose derivatives on the padded 2-D grid. It then applies
the simulator's Fourier round trip, real-space native-box clip, and native FFT
to the model and every Jacobian column. Production code is not changed.

The test holds the raw truth volume, 48 exact noiseless orientations, shell
range, and one global amplitude fixed. It reports:

- fitted residual and amplitude scale for all exact observations;
- exact-pose rotation and shift gradient RMS;
- three centered finite-difference errors for the transformed Jacobian;
- exact-pose drift under a bounded test-only LM solve;
- recovery from one known joint rotation-and-shift perturbation.

Finite evidence and a minimum finite-difference error below `5e-3` are hard
numerical gates. The `0.001745` rad and `0.05` pixel exact-pose limits are
scientific diagnosis thresholds. The recovery arm must reduce both rotation
and shift error. Scientific outcomes are printed without stopping the remaining
measurements.

`MATCHED_WINDOW_SUPPORTED` means the explicit clip model is a viable numerical
correction. It does not approve its expensive FFT-per-particle implementation
for production and does not clear reconstructed-reference bias. The mother
suite remains ten cases because this prototype is available only through its
explicit selector.

## Reconstructed-reference bias diagnostic

The truth-volume matched-window control passed on Oracle Linux. The remaining
`case=reference_bias` matrix keeps eight clean exact probe particles fixed and
changes only the reference construction or objective shell range. Its 22 arms
cover:

- the exact truth-volume control;
- own-half and opposite-half PCG references after 2, 4, and 8 iterations;
- a two-iteration reference that excludes all eight probe particles;
- low (`2:4`), middle (`2:8`), and full (`2:12`) shell ranges;
- masked and unmasked two-iteration references;
- two-iteration references with lambda `1`, `100`, and `2000`, compared with
  the normal `1e-3` reference.

Every arm uses the explicit matched-window value and five-column Jacobian. It
fits one global real amplitude at the exact poses and holds that nuisance scale
fixed during both pose batches. It reports fitted residual, exact-pose
gradients, raw pose drift, drift
after removal of the best common rigid gauge, controlled perturbed-pose
recovery, accepted steps, and objective changes. The component flags report
whether ownership, gauge, frequency, iteration count, regularization, or the
support mask reduces the normalized exact-pose error by at least 25 percent.

Only finite accounting and the previously validated truth-volume stationarity
and recovery controls are hard gates. Open reference-bias hypotheses are
reported without stopping later arms. The fitted amplitude is an isolation
device, not an approved production scaling policy. The bounded test-only solve
measures the matched objective landscape; it does not replace the production
LM gain-ratio contract. `EVIDENCE COMPLETE` is not production approval.

Run the new focused case and the mother suite directly:

```bash
simple_test_continuous_3D_pcg_refinement case=reference_bias \
  2>&1 | tee continuous_3D_pcg_reference_bias.log

simple_test_continuous_3D_pcg_refinement \
  2>&1 | tee -a continuous_3D_pcg_reference_bias.log
```

For one unattended package containing every operator-isolation diagnostic and
the mother suite, run:

```bash
bash production/tests/continuous_3D_pcg_pose_validation/run_oracle_operator_diagnostics.sh \
  --output-parent "$HOME/Projects"
```

The runner continues after an independent case failure and writes `summary.md`,
individual logs, a combined log, status, and source/checksum manifests in one
timestamped `continuous_3D_operator_diagnostics_*` directory.

## Complete operator-contract matrix

`case=operator_contract` is the pre-production decision gate added after the
forward-path review. It does not change the production operator. In one run it
separates four prediction models:

1. bare PCG gather, $G(V)$;
2. inverse-envelope gather, $G(E^{-1}V)$;
3. finite-box gather, $W G(V)$;
4. finite-box inverse-envelope gather, $W G(E^{-1}V)$.

Here $W$ is executed explicitly as padded inverse FFT, simulator round trip,
central native clip, optional production edge taper, native FFT, and native
plane extraction. The same transform is applied to each trial prediction.
The test uses an independent padded SIMPLE projector for observations; it does
not generate observations with the PCG gather.

The matrix repeats the same asymmetric object in boxes 24, 32, and 48 while
holding the support radius at 10 pixels. This tests whether the previously
observed hard-clip mismatch is mainly a two-pixel-margin fixture effect. For
each margin it measures clean exact truth with and without the production edge
taper, then reconstructs tapered observations with 2, 4, and 8 matrix-free PCG
iterations. Every truth or reconstructed reference is evaluated with all four
models. Each arm reports fitted residual, fixed amplitude, normalized
objective, and finite-difference rotation and shift gradient RMS.

The independent SIMPLE projector uses the same KB interpolation family and
therefore carries the same gather envelope. Its exact finite-box match cannot
test deapodization. A separate envelope-free control follows the established
`test=pcg_recon` Stage 8 construction: it generates $G(E^{-1}V)$ observations
and compares the bare and inverse-envelope pose models. This control must match
only the inverse-envelope model and must distinguish the bare gather. It
isolates the PCG envelope contract; it does not claim independence for the KB
value or derivative implementation.

Observed-particle variance normalization is deliberately not applied to this
clean deterministic matrix: a noiseless solvent variance is not a valid
normalization control. The fitted real amplitude removes one global gain and
shells 0 and 1 are excluded. A later noisy production test must freeze the
normalization statistics measured from each observed particle and reuse them
for its prediction and Jacobian, as required by the SPEC.

Run the focused diagnostic:

```bash
simple_test_continuous_3D_pcg_refinement case=operator_contract \
  2>&1 | tee continuous_3D_pcg_operator_contract.log
```

Or run the unattended operator package shown above. The focused case hard-gates
finite accounting, the independent truth controls for $WG$ with and without
taper, and the envelope-free PCG control for $G(E^{-1}V)$. Reconstructed-
reference, iteration, and margin comparisons are collect-all scientific
evidence. `EVIDENCE COMPLETE` does not authorize a finite-box production
correction.

## Required diagnostic sequence

Use these gates in order. Do not interpret later workflow output as a substitute
for an earlier numerical contract.

1. **Operator and reference isolation:** run
   `run_oracle_operator_diagnostics.sh`. This includes the fixed-reference,
   forward-path, matched-window, reference-bias, complete operator-contract,
   and mother-suite cases.
2. **Truth-controlled production workflow:** after a correction is selected,
   rerun `run_truth_diagnostic.sh` for clean/noisy, exact/perturbed, and
   full/FSC-limited production arms.
3. **Scientific handoff:** retain the complete truth-matrix directory and use
   `doc/implementation_notes/completed/continuous_3D_pose_end_polishing_history_and_handoff.md` to
   review the result with the refinement owners.
4. **Deferred component comparisons:** `run_oracle_validation.sh` and
   `.codex/hpc_betagal_prepare.sh` end in an explicit top-level
   `prg=reconstruct3D` A/B. They do not prove integration with standard
   `abinitio3D` or `refine3D_auto`, which do not currently propagate
   `pcg_pose_polish`.

Gate 1 passed in the retained Oracle operator package. Gate 2 is the only
remaining runtime gate for the current investigation. Gate 3 is a review and
product-decision boundary. The explicit-command comparisons in Gate 4 are
deferred and must not be reported as normal workflow validation. Apoferritin
also remains deferred until both the workflow contract and octahedral symmetry
analysis are defined.
