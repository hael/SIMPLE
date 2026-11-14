# Bootstrap noise annealing for 2D matching

Status: two-stage annealing implemented and locally validated, 2026-11-14.

## Contract

`sigma2_anneal=yes` is the default of solve2D and solve2D_chunks (since
2026-11-14; it was opt-in before); standalone refine2D defaults to `no`. There
is no separate anchor-only mode. On solve2D and solve2D_chunks,
keep the smoothed bootstrap through the early stages, then increase the
measured-spectrum weight from 0 to 0.5 in the penultimate stage and from
0.5 to 1 in the final stage. The terminal greedy pass uses measured sigmas.
Standalone refine2D ramps from 0 to 1 over its iteration budget.

For each half and shell, scoring uses
`exp((1-w)*log(anchor) + w*log(measured))`. Zero weight copies the anchor;
unit weight attaches the measured array directly. Reconstruction and
canonical residual updates always retain measured spectra. Poses and
insertion-posterior weights can change the class sums indirectly.

Each half's native bootstrap is smoothed in log variance with a seven-point
quadratic Savitzky–Golay filter. The first/last seven samples supply the
endpoint polynomial fits. The validated `.anchor` companion holds the saved
anchor; it is reused throughout annealing and on continuation.
A missing anchor cannot be manufactured from residual spectra.

Each ramp spans its stage's actual iteration budget, including the first and
last iterations. A one-iteration stage uses its endpoint. Coverage iterations
beyond the budget hold that endpoint. The two annealing stages cannot be
omitted; an early convergence verdict cannot stop an unfinished sigma ramp.
The beta schedule and its convergence state are unchanged. In particular,
the sigma endpoint does not depend on completion of the beta tail. The seed
pass of a seeded solve2D run stays at zero weight.

Internal typed keys `sigma2_weight_start`, `sigma2_weight_end` and
`sigma2_weight` carry the stage bounds and current weight. The refine2D owner
sets the current weight before each iteration and explicitly forwards it to
partition workers. The readout includes weight, beta, actual matching band,
and the measured/anchor diagnostics. All weights must be finite, in [0,1],
with ordered stage endpoints. No new tuning parameters are exposed in the UI.

A fresh solve2D run anneals without the switch. Keep
`sigma2_state.bin.anchor` beside its canonical store when continuing that run.
`sigma2_anneal=no` selects measured scoring (and is needed to continue a run
whose store has no anchor).

## Two-stage annealing implementation and validation

The SLC run supplied by the maintainer nearly eliminated the spokes with the
anchor alone (weight 0 throughout) but weakened the best classes. Its iteration-30 geometric
mean measured/anchor ratio is 0.303 in both halves, with beta 0.7401; its
five-step beta tail stops at step 4 because the final stage exhausts its
budget. The saved spectra after terminal iteration 31 have mean ratios
0.256/0.254 over shells 2:54. These motivate releasing the anchor over the
last two stages, independently of that beta endpoint. They do not establish
that annealing improves real-data quality.

Validation for this extension: exact blend endpoints and independent square/
fourth-root expectations, scoring/self-term/posterior temperature with both
halves and changing bands, unchanged reconstruction weights and fixed-pose
residuals, stage schedules and omission guards, parser and worker forwarding.
Incremental Debug build, one fast gate and small shared/distributed 2D replays;
no repeat of the long 3D regression.

Verified: macOS GNU Fortran 16 Debug incremental build, complete UI validation
and installation pass without compiler warnings. Focused stage rules pass
159 assertions; anchor/annealed-scoring checks pass 152. All 15 fast
suites pass, including 6214 polar-area assertions. The serial fast gate takes
134.9 seconds (warning-only budget overrun; the sandbox denies its CPU query).

Both 79-particle replays, sampling 40 particles, complete through terminal
iteration 31. Their stages 5 and 6 emit exactly
`[0, .125, .25, .375, .5]` and `[.5, .625, .75, .875, 1]`; terminal weight is 1.
The final-stage beta sequence remains `[.3, .3, .4054, .5477, .7401]`, so the
sigma endpoint is reached despite the truncated beta tail. Earlier weights
are zero. Commands, source-diff/input hashes and the asserted schedules are
preserved under `annealing_extension_verification` in the linked validation
JSON. No new SLC quality result or Linux/BOX validation is claimed.

## Original anchor-only implementation and review

The refinement owner prepares a small `.anchor` companion to the canonical
sigma store before workers start. It is generated only from a particle-power
bootstrap, atomically published, and checked against the store's grid and
particle layout. A rebuilt bootstrap replaces the companion. A continuation
without its companion cannot substitute a residual spectrum. Workers only
read it. Each chunk has its own bootstrap and companion.

The sigma domain owns smoothing and artifact validation. The existing sigma
object keeps separate measured and scoring arrays; only the matcher attaches
the anchor array to the polar calculator, with its usual cache generation.
The existing general SG routine uses a quartic fit and asymmetric endpoint
handling; it does not implement this contract and is left unchanged.

The matcher checks the anchor's shell count before attachment and refreshes its
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
temperature use the anchor while reconstruction and residuals use
the measured one. Verify the switch's parser, default, UI and child forwarding.

Build incrementally in the existing Debug build, run the affected unit areas,
the fast gate, and a simulated workflow. Exercise enabled shared-memory and
distributed solve2D through the terminal pass; verify the off path against a
fixed-seed baseline. Real SLC quality remains the maintainer's experiment.

Verified locally so far: Debug incremental build and complete UI validation;
registry, source descriptors, FLEX DAG, CLI chain and tester linkage; all 100
focused `unit_pftc_align2D3D suite=sigma2_anchor` assertions. The tests cover
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
The long `simulated_workflow suite=1jxy` regressions were stopped at the
maintainer's request after their shared-memory truth checks passed; no complete
full-workflow pass is claimed.

## Fixed-beta truth comparison

The preserved [validation data](../data/sigma2_anchor_anneal_2D_validation.json) contain
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
| Anchor-only SG bootstrap | 33.0% | 6.02 degrees | 79.56 degrees |

Class labels are matched one-to-one by Hungarian assignment. Pose errors use
correctly assigned particles only (149 control, 99 anchor-only), remove one circular
reference-angle offset per class, and use the same rotation convention selected
from the control. These are final results after the terminal greedy pass.
The independent two-part anchor-only run also completes, with 33.0% class accuracy
and 3.41-degree median error; shared and distributed trajectories need not match.

This single difficult fixture favors the default path. It is not evidence that
the anchor alone cures the SLC spokes. In the anchor-only shared run, iteration 20's even-half
mean log ratio is -0.0888 (geometric mean ratio about 0.915), while the class
n_eff is 1.20 and angular n_eff 1.9. The control's corresponding readouts are
1.09 and 1.8; later trajectories and stage durations diverge. Keep the switch
opt-in and distinguish implementation validation from real-data qualification.

After the concurrent terminology/API renames finished, the combined tree built
and installed without warnings and passed all 15 fast suites. All three truth
runs were repeated on that build: their metrics are identical, and respectively
79, 82 and 94 scientific artifacts match the earlier control, anchor-only and
distributed runs (canonical layout digests excluded). The validation JSON
preserves the combined source fingerprint and these comparisons.

The validation JSON records the commands as run, under the switch's name at the
time (`sigma2_frozen`); its logs were kept in a temporary directory.
Navigation indexes were regenerated; they also refresh pre-existing stale entries.
No Linux/BOX or real-data quality validation is claimed.
