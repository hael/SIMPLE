# PfCRT Nov12 regression investigation

## Contract and scope

Preserve the current sampling policy, size-planned ladder, gold-standard switch,
and the maintainer's uncommitted removal of stage omission with a three-iteration
minimum. Fix premature convergence within that design; do not roll back workflows.
The convergence object owns progress history. No new CLI keys or global state.

## Evidence and limitations

The supplied `PfCRT_regression/Nov09_2026` logs identify `3fb0b00b`;
`Nov12_2026` identifies `f2e907c2`. All ten runs in each set finish normally.
The maintainer confirms all Nov09 maps are acceptable and Nov12's reported
gold-standard resolutions reflect the map errors.

| Measurement | Nov09 | Nov12 |
|---|---|---|
| Active particles | 16,155 | 16,771 |
| Sampling | 100%, no trailing | 100%, no trailing |
| Stage 3 iterations | 17 | 3–5 |
| Stage 4 iterations | 15 | omitted in every run |
| Stage 5 iterations | 12 | 12 |
| Stage 8 iterations | 12–18 | 8–17 |
| Final FSC=0.143 resolution | 4.03–4.09 A | 4.25–6.90 A |

Nov12 stage 5's first mean angular update is 48.531–55.038 degrees, compared
with Nov09's 20.220–21.340 after stage 4. Stage 5 finishes at FSC=0.143
4.777 A in every Nov09 restart and 4.852–6.210 A on Nov12: the deficit is
already visible before the late backend and gold-standard transition.
Nov12 restart 9 exits stage 5
while FSC=0.143 improves 6.75 -> 6.47 -> 6.21 A, because its FSC=0.5 is flat.
That exit is at its twelfth iteration, also its budget: fixing that guard
alone would not extend this recorded stage. The actual shortening is in
stage 3 and the omitted stages; stage 5 demonstrates the missing gain guard.
At the start of stage 8, several runs anneal after three checks even while
the orientation overlap is rising. Stage 3 exits also print 13.500 A as
strictly past a 13.500 A boundary: the comparison has no rounding tolerance.
Eight of ten Nov12 restarts take that equal-shell stage-3 resolution exit;
restart 2 takes the plateau exit and restart 8 crosses to a genuinely finer shell.
Single-precision arithmetic with the logged native box 300, pixel size 1.035 A
and crop 108 reproduces it: crop pixel size 2.874999762 A gives shell 23 at
13.49999905 A; the native-grid shell is 13.5 A. They represent the same shell.

The late backend differs (Nov09 PCG, Nov12 gridding), and half-specific
matching became the default between these versions. These and the changed
particle selection prevent attribution of the entire resolution difference
to one change. The supplied directories contain logs, not the corresponding
particle projects/stacks for a matched real-data rerun.

## Implementation and validation plan

1. Regress the production convergence path for a flat FSC=0.5 with improving
   FSC=0.143; accumulated small overlap gains; annealing during those gains;
   and equality at a shell boundary within floating-point roundoff.
2. Make resolution exits observe both resolution measures when available.
   Share stage-local overlap progress between resolution exits, the plateau
   rule and single-state posterior annealing. Count a cumulative new overlap
   high of 0.02, retaining the existing threshold; a fall and rebound is not
   progress. Preserve the existing one/two-check stall windows and budget tail.
3. Give the strict resolution-boundary comparison a small relative floating
   tolerance. Keep the minimum iteration count and all maximum budgets.
4. Build incrementally in the existing Debug tree, run the convergence and
   stage-rule suites, the fast gate, and simulated_workflow_1jxy. The latter
   explicitly defers gold standard and cannot validate PfCRT recovery.
5. A matched PfCRT ten-restart run remains the acceptance test for recovery:
   all stages execute, progress survives flat FSC shells, and final maps and
   gold-standard resolution recover without disabling the new settings.
   The maintainer will rerun the Nov12 inputs and settings with the combined
   diff; inspect both map quality and the across-restart resolution spread.

## Validation record

- Added thirteen assertions through the existing convergence suite. Before the
  production fix, seven fail, including improving FSC=0.143, cumulative overlap,
  annealing and boundary equality; the executable returns failure.
- After the fix: all 161 convergence assertions pass in the existing macOS
  gfortran-16 Debug build. Both incremental build invocations pass without
  compiler warnings (the first also brought the pre-existing tree up to date).
- Registry, descriptor/source-unit and FLEX DAG checks pass. Full fast gate:
  15/15 pass, 147.2 s in Debug (warning-only 60 s budget overrun). The touched
  `unit_cart_align3D` area also passes in reverse suite order: 1,077 assertions.
  Stage rules pass in `unit_reconstruction`: 258 assertions in that sub-suite.
- Installation passes. `simulated_workflow_1jxy` passes in 4,563.02 s with
  all eight solve3D stages executed in both shared and two-part distributed
  arms. This fixture explicitly sets `gold_std_stage=9` for solve3D, so it
  does not exercise that transition. A snapshot of all 874 Fortran source
  hashes records the compiled source. The final audit finds one subsequent
  comment-only cleanup in `simple_solve3D_controller`; executable source is
  unchanged. The maintainer's other initial changes remain preserved.
- The shared-memory arm reports truth correlation 0.888023, truth
  FSC=0.143 at 5.20 A; expected grid and wrong-truth control pass. Pose
  recovery is a profile miss: median pair error 16.43979 degrees, fraction
  within 5 degrees 0.4975657 (per-run criteria <=10 degrees and >=0.5).
  The 0.80 profile floor is the fraction of successful runs over seeds,
  not the within-5-degree threshold. A single miss is report-only by policy;
  `ea23942a5` already documents a seed-specific miss. No matched pre-fix
  workflow run was performed, so this is not attributed to the patch.
- The two-part distributed arm reports truth correlation 0.8898241, truth
  FSC=0.143 at 4.378947 A; expected grid, wrong-truth control and two-part
  execution checks pass. Pose recovery succeeds: median pair error
  0.9907107 degrees and fraction within 5 degrees 0.9500162. These are
  single-seed observations, not estimates of the across-seed success rate.
- Logs and parsed real-data trajectories:
  `/private/tmp/pfcrt_nov12_review/`. No real-data rerun or Linux build yet.
