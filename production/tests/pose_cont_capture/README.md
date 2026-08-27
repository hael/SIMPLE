# Cartesian pose-capture diagnostic

This package runs one isolated numerical experiment. It does not activate pose
polishing in `refine3D`, `reconstruct3D`, or another production workflow.

## Question

The test asks how far the existing five-parameter Jacobian and bounded
Levenberg--Marquardt solver can return to a known pose under its clean, matched
Cartesian Fourier forward model.

The test keeps the following items fixed:

- one asymmetric 24-cubed volume;
- one known rotation and image shift;
- Fourier shells 2 through native Nyquist;
- one clean observation made by the tested operator;
- one LM policy and iteration limit.

The separate sweeps change only one starting coordinate. They cover rotations
through 15 degrees and shifts through 15 pixels on each axis and sign. Shift
cases at 5.5, 6, 6.5, and 7 pixels bracket the previously observed boundary.
The full LM system can still update all three rotation coordinates and both
shifts.

The joint matrix contains these cases:

- two small controls at 2 degrees/0.5 pixel and 5 degrees/1 pixel;
- a complete 3 by 3 boundary grid at 10, 12.5, and 15 degrees crossed with 3,
  4, and 5 pixels;
- three mixed-sign cases;
- three simultaneous multi-axis rotation and two-axis shift cases.

The final matrix has 138 trials. It brackets the clean matched-operator capture
boundary. It does not define a production operating range.

The same package also runs a mechanism diagnostic on three deterministic
asymmetric volumes. It uses five boundary seeds and compares:

- direct joint LM;
- joint LM with a cumulative 15-degree/5-pixel trust region around the input
  seed;
- shift-only LM followed by joint LM;
- rotation-only LM followed by joint LM.

The cumulative guard is optional and does not change existing callers. It
limits total movement from the discrete input pose. It does not claim that a
15-degree or 5-pixel update is scientifically correct.

## Oracle Linux run

Compile SIMPLE first. Then load the matching compiler runtime and set the
SIMPLE runtime path:

```bash
module load gcc/15.2.0
module load miniconda/25.9.1
export SIMPLE_EMAIL="my.name@uni.edu"
export SIMPLE_QSYS="local"
export SIMPLE_PATH=/home/hossainm7/Projects/hael_SIMPLE-rsync-test/build
export PATH=${SIMPLE_PATH}/scripts:${SIMPLE_PATH}/bin:${PATH}
python3 --version
```

The controlled Oracle environment is Miniconda 25.9.1 with Python 3.13. The
analyzer also avoids unnecessary new typing syntax so an older system Python
cannot invalidate an otherwise complete numerical run.

Run the single diagnostic:

```bash
cd ~/Projects/hael_SIMPLE-rsync-test

bash production/tests/pose_cont_capture/run_oracle_pose_capture.sh \
  --output-parent "$HOME/Projects"
```

The script creates one directory named
`continuous_3D_pose_capture_YYYYMMDD_HHMMSS`. Download that complete directory.
Review `analysis/summary.md` first.

The runner also executes `case=pose_recovery` with no new optional arguments.
This focused regression checks that the trajectory, mask, and cumulative-guard
extensions do not change the established default LM path.

## Evidence

The result directory contains:

- `run.log` and `logs/pose_capture.log`;
- `artifacts/pose_capture.csv`;
- `artifacts/pose_mechanism_summary.csv`;
- `artifacts/pose_mechanism_trajectory.csv`, with all five pose coordinates
  after every accepted step;
- `artifacts/pose_objective_paths.csv`, with sampled paths from each seed to
  truth and to the direct-joint endpoint;
- `artifacts/truth_volume.mrc` and `truth_observation.mrc`;
- the three mechanism-test volumes and their observations;
- initial and final predictions and residuals for each trial;
- `analysis/summary.md` and `summary.json`;
- the exact source inputs and SHA-256 manifests.

`PASS` means that the package is complete and its numerical integrity checks
passed. It does not mean that every perturbation recovered the known pose. The
summary reports each separate axis, sign, and magnitude and each joint vector.
A final-to-initial truth-error ratio below 1 means that the trial moved closer
to truth.

The analyzer reports scientific recovery warnings separately. A warning occurs
when a nonzero rotation or shift error retains more than 10 percent of its
initial value. A warning also occurs when an initially zero coupled error grows
above 0.1 degree or 0.1 pixel. These limits identify gross capture loss. They
are not production acceptance criteria and do not change package `PASS`.

The test does not prove behavior with noise, CTF variation, an independent
forward generator, real particles, or production `refine3D` state.

## Reanalyze an existing completed run

The analyzer supports Python 3.6 and later. If an earlier analyzer failed after
the Fortran log printed `CONTINUOUS_3D_POSE_CAPTURE: EVIDENCE COMPLETE`, do not
rerun the numerical test. Synchronize the corrected scripts and run:

```bash
cd ~/Projects/hael_SIMPLE-rsync-test

bash production/tests/pose_cont_capture/reanalyze_oracle_pose_capture.sh \
  "$HOME/Projects/continuous_3D_pose_capture_YYYYMMDD_HHMMSS"
```

The recovery script verifies the original completion marker, regenerates the
analysis, updates `STATUS.txt`, and refreshes `MANIFEST.sha256`.
