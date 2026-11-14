# Frozen bootstrap noise spectra for 2D matching

Status: implemented; workflow validation in progress, 2026-11-14.

## Contract

`sigma2_frozen=yes` is an opt-in experiment on refine2D, solve2D and
solve2D_chunks; the default is `no`. Save each half's native-grid bootstrap
spectrum once, smooth its logarithm with a seven-point quadratic
Savitzky–Golay filter, and exponentiate. The first/last seven samples supply
the endpoint polynomial fits. The spectrum must have at least seven shells.

Matching, continuous in-plane refinement and insertion-posterior evaluation
use these fixed spectra through the terminal greedy pass. Residual estimation,
canonical sigma commits and reconstruction retain measured spectra. There is
no annealing, posterior-weighted residual update, new temperature or floor.
Selected poses and insertion posterior weights can change the class sums
indirectly. The separation preserves measured reconstruction noise weights,
not identical class averages between the two modes.

Use the same fresh solve2D command with `sigma2_frozen=yes`. Keep
`sigma2_state.bin.frozen` beside its canonical store when continuing that run.
Leaving the switch out, or setting `sigma2_frozen=no`, retains the existing path.

## Implementation and review

The refinement owner prepares a small `.frozen` companion to the canonical
sigma store before workers start. It is generated only from a particle-power
bootstrap, atomically published, and checked against the store's grid and
particle layout. A rebuilt bootstrap replaces the companion. A continuation
without its companion cannot substitute a residual spectrum. Workers only
read it. Each chunk has its own bootstrap and companion.

The sigma domain owns smoothing and artifact validation. The existing sigma
object keeps separate measured and scoring arrays; only the matcher attaches
the frozen array to the polar calculator, with its usual cache generation.
The existing general SG routine uses a quartic fit and asymmetric endpoint
handling; it does not implement this contract and is left unchanged.

The matcher checks the frozen shell count before attachment and refreshes its
weights once. Each owner logs both halves' unweighted mean natural-log ratio
`log(measured/anchor)`, the log ratio at the low-pass edge, the matching shell
bounds, beta, and the number of shells above the anchor. These describe the
incoming measured spectrum used by that iteration, before its residual commit.
They appear beside the existing class and in-plane posterior n_eff readouts;
the greedy terminal pass has no posterior draw. The shared matcher supplies its
actual band. The distributed owner obtains the workers' searched low-pass limit
from their merged particle rows, before restoration, rather than a stage-start
estimate. Both workers' benchmark bands matched all 30 diagnostic iterations.

Measured power need not lie below the anchor. The exact residual includes the
reference power and data/reference cross term; pose and model errors change
both. SG smoothing alone can put the anchor below the raw bootstrap at some
shells. Above-anchor counts are diagnostic, not a validity failure. For a fixed
residual, the shell's likelihood coefficient relative to measured-sigma scoring
is `measured/anchor`, on top of the common beta; this need not always flatten it.

## Validation

Check SG against independent polynomial and impulse expectations, positivity,
half separation, round-trip persistence, stale/truncated/nonfinite inputs,
refusal to manufacture a bootstrap from residuals, and preservation through
later measured-spectrum updates. Check that polar scoring and posterior
temperature use the frozen spectrum while reconstruction and residuals use
the measured one. Verify the switch's parser, default, UI and child forwarding.

Build incrementally in the existing Debug build, run the affected unit areas,
the fast gate, and a simulated workflow. Exercise enabled shared-memory and
distributed solve2D through the terminal pass; verify the off path against a
fixed-seed baseline. Real SLC quality remains the maintainer's experiment.

Verified locally so far: Debug incremental build and complete UI validation;
registry, source descriptors, FLEX DAG, CLI chain and tester linkage; all 100
focused `unit_pftc_align2D3D suite=sigma2_frozen` assertions. The tests cover
the analytic SG impulse response and quadratic endpoint preservation, both
halves, residual-generation persistence, rebuilt bootstraps, invalid artifacts,
two search bands, scores and posterior temperature, unweighted residual records,
and unchanged reconstruction inverse weights.

All 15 fast suites pass, including 6033 assertions in the polar-alignment
area. The sandbox denied the gate's CPU query, so the final combined-tree gate ran serially in 146.0 s;
the 60 s budget overrun is a warning. Installation passed.

The 79-particle solve2D replays pass with the switch enabled in shared memory
and with two local partitions. Additional `nsample=40` runs pass in both modes,
including the conditional terminal greedy pass. The saved spectrum retains
bootstrap generation 1 while the canonical store advances through residual
commits. The off replay matches all 84 map/stack/class-state artifacts from the
saved executable, all 79 particle records, all four class records, and canonical
sigma values. Compare numerical project fields after removing JSON thumbnail
indices; ignore the canonical store's run-specific layout digest. The saved
executable writes its project and canonical store beside the input, whereas
the rebuilt executable places them in the execution directory.

These replays use a stable snapshot of the simulated project's completed
2D result, rerun with fresh bootstrap initialization and identical seeds.
The full `simulated_workflow suite=1jxy` regression is still running.

## Fixed-beta truth comparison

The preserved [validation data](../data/frozen_sigma2_2D_validation.json) contain
the commands, seeds, source/input hashes, quality measurements, and all iteration
readouts; the temporary logs are not the only record. The fixture has 300
particles (50 per known projection), six 1JXY projection templates downsampled
to 64 pixels at 3.9 A/pixel, uniform random in-plane rotations, and independent
Gaussian noise at masked signal/noise power ratio 0.15. It has no CTF or shifts.
Both arms use the same input and `SIMPLE_SEED=20261114`, `nsample=150`,
`prob_beta=0.3 prob_beta_end=0`; every observed draw remains at beta 0.30.

| Shared-memory arm | Class accuracy | Median in-plane error | 90th-percentile error |
| --- | ---: | ---: | ---: |
| Default measured sigma | 49.7% | 3.04 degrees | 8.54 degrees |
| Frozen SG bootstrap | 33.0% | 6.02 degrees | 79.56 degrees |

Class labels are matched one-to-one by Hungarian assignment. Pose errors use
correctly assigned particles only (149 control, 99 frozen), remove one circular
reference-angle offset per class, and use the same rotation convention selected
from the control. These are final results after the terminal greedy pass.
The independent two-part frozen run also completes, with 33.0% class accuracy
and 3.41-degree median error; shared and distributed trajectories need not match.

This single difficult fixture favors the default path. It is not evidence that
freezing cures the SLC spokes. In the frozen shared run, iteration 20's even-half
mean log ratio is -0.0888 (geometric mean ratio about 0.915), while the class
n_eff is 1.20 and angular n_eff 1.9. The control's corresponding readouts are
1.09 and 1.8; later trajectories and stage durations diverge. Keep the switch
opt-in and distinguish implementation validation from real-data qualification.

After the concurrent terminology/API renames finished, the combined tree built
and installed without warnings and passed all 15 fast suites. All three truth
runs were repeated on that build: their metrics are identical, and respectively
79, 82 and 94 scientific artifacts match the earlier control, frozen and
distributed runs (canonical layout digests excluded). The validation JSON
preserves the combined source fingerprint and these comparisons.

Logs and replay commands are under `/private/tmp/simple-frozen-sigma/`.
Navigation indexes were regenerated; they also refresh pre-existing stale entries.
No Linux/BOX or real-data quality validation is claimed.
