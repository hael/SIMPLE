# Making Cartesian continuous refinement as good as the polar search, without losing its speed

Updated 2026-11-11. Active work plan. Delivered implementation and historical
evidence are in [the done record](../completed/refine3d_local_cartesian_priors_done.md).

## 1. The objective

`refine3D_local` can refine a particle's 3D pose two ways:

- **the polar route** (`refine=prob_neigh`), which scores a grid of candidate
  poses and picks one from them;
- **the Cartesian route** (`refine=cont`), which slides continuously downhill
  from the pose we already have.

We have measured that **the polar route produces better maps, and runs about
five times slower.**

**The goal of this project is to make the Cartesian route produce maps as good
as the polar route, while staying substantially faster than it.**

The budget: the Cartesian route may become up to about **twice** as expensive as
it is today. That still leaves it roughly 2.5 times faster than polar. Anything
that buys quality inside that budget is on the table; anything that spends more
needs an explicit decision.

Everything else in this document serves that objective. If a task does not help
close the quality gap or protect the speed advantage, it is not a priority.

### A note on what came before

Most of the recent work on this code added an optional penalty called
`cont_prior`, which discourages a particle from moving far during one
refinement. That work is complete and tested, and section 7 lists the loose ends.

But be clear about its status: **the penalty is not the reason the polar route
wins, and tuning it is unlikely to close the gap.** It restrains motion; the
evidence in section 3 says the Cartesian route's problem is that it cannot move
to the right answer in the first place. Treat `cont_prior` as finished
infrastructure that we may use later, not as the current line of attack.

## 2. Vocabulary

The rest of this document uses these terms. If they are already familiar, skip
to section 3.

**Pose.** The five numbers that place a particle image against the 3D reference:
three for the rotation, two for the shift in the image plane.

**View direction and in-plane angle.** A 3D rotation splits naturally into two
parts: which direction the molecule is facing (the *view direction*, two
numbers), and how much the resulting 2D projection is spun around in the image
plane (the *in-plane angle*, one number). This distinction matters a great deal
below.

**Loss.** A number saying how badly a pose fits the image. Lower is better. Both
routes use the same kind of loss, but they normalize it differently, which is
why section 5 insists on checking they agree before comparing them.

**Posterior.** The probability distribution over poses, given this image. A
sharp posterior means the image pins the pose down; a broad one means many poses
fit nearly as well. A **multimodal** posterior means there are several separated
poses that each fit well — the image is genuinely ambiguous.

**Basin.** A valley in the loss landscape. Downhill search from inside one basin
can only ever reach the bottom of that basin, never a better one on the far side
of a ridge.

**MAP versus draw.** The *MAP* (maximum a posteriori) pose is the single best
one. A *draw* is a random pick, where better poses are more likely but worse
ones are still possible. Committing a draw rather than the MAP lets a particle
escape a wrong answer on a later iteration. The polar route draws; the Cartesian
route takes the best.

**Tempering.** Flattening the posterior before drawing from it, so that worse
poses get picked more often than they strictly deserve. Controlled by `prob_beta`
(default 0.3; below 1 means flatter, i.e. more exploratory). It is a way of
keeping options open early in a refinement.

**Calibration.** Whether a stated uncertainty is honest: if we say we are 90%
sure the pose is within some cone, is the truth inside that cone 90% of the time?
This matters in section 8, not before.

## 3. Why the polar route wins

Three differences. They are listed in the order I believe matters, but
**section 4 exists precisely because that order is a hypothesis and we can test
it cheaply before committing to any of it.**

### Difference 1: the Cartesian route has a wall in the in-plane angle

Both routes limit how far a particle may move from its stored pose, using the
same parameter, `athres_cont` (10 degrees). But they apply it to different
things.

The polar route limits the **view direction** to 10 degrees, and then considers
**every in-plane angle, all 360 degrees of it**
(`simple_eul_prob_tab_neigh.f90:764`, and the in-plane scan in
`simple_polarft_corr.f90:697`).

The Cartesian route limits the **whole rotation** to 10 degrees
(`simple_cartft_pose_opt.f90:745`). An in-plane rotation counts against that
same budget. So if a particle needs to turn 8 degrees in view direction, it has
only about 6 degrees of in-plane left — and it can never fix an in-plane error
larger than 10 degrees at all. Trial poses outside the bound are rejected before
they are even evaluated.

This is an artificial wall, it is the largest difference in what the two routes
are allowed to do, and **removing it costs nothing**. It does not by itself make
the Cartesian search find a distant in-plane angle — downhill search still has to
get there — but at present it is forbidden even from trying.

A smaller related point: the polar bound folds over the molecule's symmetry and
the Cartesian bound does not, so for symmetric molecules the Cartesian support is
narrower than the number suggests.

### Difference 2: one downhill run versus a hundred and fifty candidates plus a lottery

The Cartesian route performs **exactly one** downhill optimization per particle
per iteration, starting from the stored pose, and keeps the result **only if it
is better** (`simple_strategy3D_cont.f90:193`, `simple_cartft_pose_opt.f90:423`).
It is deterministic and strictly monotone.

Two consequences, both bad. It cannot cross a ridge, so a particle that starts in
the wrong basin stays there — not just this iteration, but **every** iteration,
forever. And when the image is genuinely ambiguous, it silently commits to one
interpretation with no record that there were others.

The polar route scores roughly 150 candidate directions inside the cone, each at
its own in-plane angle, and then **draws** the committed pose from a tempered
posterior over all of them (`simple_eul_prob_tab_neigh.f90:961`,
`simple_eul_prob_tab_utils.f90:1633`). It can commit a pose that is worse than
the current one, which is exactly how a particle climbs out of a false basin.
The code even counts how often the draw differs from the best candidate, because
it usually does.

### Difference 3: the polar route hedges in the reconstruction; the Cartesian route does not

When a particle's pose is ambiguous, the polar route does not force a choice. It
inserts the particle into the 3D reconstruction at **up to three poses at once**,
weighted by posterior probability (`topk`, default 3;
`simple_eul_prob_tab_utils.f90:1813` and `simple_matcher_3Drec.f90:314`).
Ambiguous particles therefore contribute a blurred, honest average instead of a
confident wrong answer.

**The Cartesian route has no equivalent of this at all.** Every particle enters
the reconstruction at exactly one pose.

This is a strong suspect for the quality gap, because its effect is on the map
directly rather than on the poses, and because a pose-error metric would not
reveal it.

### What this adds up to

The Cartesian route is a good *optimizer* and a poor *estimator*. It finds the
bottom of its valley efficiently and then behaves as though that valley were
certainly the right one. The polar route is a worse optimizer — its answers sit
on a grid — but it represents the uncertainty, and that is what is reaching the
map.

The encouraging part: of the three differences, one is free to remove, one is
nearly free, and the expensive one (scoring many candidates) is where our
evaluation budget actually is. See section 6.

## 4. Step 0: find out which difference actually matters — four runs, no code

Before writing any code, attribute the gap. All three differences above are
active in the polar route simultaneously, and we can switch two of them off from
the command line:

- committing the **best** candidate instead of a draw;
- `topk=1`, which inserts each particle at one pose instead of three (documented
  as reproducing the single-pose path exactly).

Run the polar route four times on the same project, same seeds, same references:

| Run | Commit rule | Reconstruction |
| --- | --- | --- |
| 1 | draw | `topk=3` (this is the current polar route) |
| 2 | draw | `topk=1` |
| 3 | best | `topk=3` |
| 4 | best | `topk=1` |

Run 4 is the polar route stripped of both mechanisms: a wide discrete search that
commits its best candidate at one pose. If run 4 is **still** much better than
the Cartesian route, the gap is mostly the search support and we should attack
difference 1. If run 4 falls back toward Cartesian quality, the gap is the
estimator and we should attack differences 2 and 3 — and the table tells us which
of the two is carrying it.

This costs four polar runs and no development time, and it decides the order of
everything that follows. **Do this first.**

Record FSC, cFAR and, on synthetic data, pose error against truth. Keep the
particle identities so the comparison can be paired (section 7).

## 5. The one check that must pass before cross-route comparison means anything

The two routes compute their loss with different normalizations: the Cartesian
route scales by whitened particle power over four, the polar route by a different
expression involving the polar ring size. Until we have confirmed that both
report the **same number for the same physical pose on the same particle**, any
statement of the form "polar scored this better" is confounded with a units bug.

So: freeze references, noise, mask, CTF, crop, band and pose; take poses found by
the polar route; rescore them with the Cartesian calculator; compare. Check the
special low-frequency bins, the interpolation error and the normalization
separately. Set the agreement tolerance from interpolation accuracy, not from the
observed difference.

This is a one-off diagnostic, not a feature. It is cheap, it protects every
comparison downstream, and it is likely to surface real bugs. Do it alongside
step 0.

## 6. The work, in order

Costs are in **model evaluations per particle**. The Cartesian route currently
averages **4.7** per particle (measured, beta-galactosidase). The 2x budget
therefore allows about **5 extra evaluations**. For scale, the polar route uses
roughly 300 to 450 evaluations per particle, but each is individually much
cheaper, which is how a seventy-fold difference in evaluation count becomes only
a five-fold difference in wall-clock time.

That per-evaluation ratio is inferred, not measured. Measure it early: one
benchmarked run of each route on the same project gives nanoseconds per polar
row from the fill profile and nanoseconds per Cartesian pass from the existing
per-pixel timing. Every cost estimate below rests on it.

### 6.1 Remove the in-plane wall — cost: zero

Split the single rotation bound into a view-direction bound and a separate
in-plane bound, so the two no longer compete for one budget. Mechanically:
decompose the trial rotation relative to the seed into view change and in-plane
change and test them separately, instead of the single combined test at
`simple_cartft_pose_opt.f90:745`. Fold the view bound over the point group, as
the polar route does.

Done when: a particle whose in-plane angle is wrong by 30 degrees can be
corrected, with the view-direction bound still enforced at 10 degrees; prior-off
results are otherwise unchanged.

### 6.2 Let the commit be a draw rather than strictly the best — cost: zero

Replace "keep the result only if it beats the seed" with a draw over
{seed, result} at the particle's temperature, times `prob_beta`. The machinery
exists and is representation-agnostic: `draw_posterior_row` and
`posterior_temperature` take a plain array of losses and a temperature. Supply
the Cartesian temperature, which is the same quantity the Cartesian loss is
already normalized by.

This is the cheapest available escape mechanism from a false basin, and it costs
no extra evaluations.

Done when: a fixture with a known false basin escapes it over several iterations
with the draw enabled and never escapes with it disabled; the result is
reproducible under the keyed random stream across thread counts and partitions.

### 6.3 Weight the reconstruction over several poses — cost: zero to small

Give the Cartesian route the `topk` behaviour. The insertion machinery
(`prob_insert_pose`, `prob_insert_table`, `prob_insert_poses`, `prob_insert_cone`)
is geometric and has no polar dependency; the consumer in
`simple_matcher_3Drec.f90:314` already handles a particle with no extras as the
plain single-pose path.

The question this raises is *which* several poses, since the Cartesian route
currently produces one. Either this follows 6.4, or an interim version inserts
the single converged pose with a weight reflecting how well determined it is.
Decide after step 0 says how much `topk` is worth.

### 6.4 Several starting points, then weight the results — cost: about 4.7 per extra start

Run the downhill search from two or three starting points instead of one and keep
all converged results, then draw the commit from them and insert the best few
into the reconstruction with posterior weights. `refine_pose` already takes the
pose as an in/out argument, so a multi-start loop needs no new optimizer, only a
per-start record and a selection rule.

Two starts fits the budget. Where the extra starts come from matters more than
how many there are: candidates from a coarse scan (6.5) are worth more than
random perturbations.

### 6.5 A cheap coarse scan to choose the starting points — cost: to be measured

The principled version of 6.4: before optimizing, evaluate the loss on a coarse
grid — in-plane angles at the seed direction are the cheap and probably highest
value choice — and start the optimizer from the best two or three distinct
basins. This is the Cartesian analogue of what the polar route does, at a tiny
fraction of the candidate count.

**Prerequisite:** the calculator has no loss-only entry point. Every evaluation
currently computes the gradient and all fifteen Hessian entries, which a scan
does not need. Add a loss-only path first; the saving is substantial and it is
useful independently.

Budget honestly before building: 36 in-plane angles at 10-degree spacing, even at
a third of the cost of a full evaluation, is roughly 12 full-evaluation
equivalents — over budget on its own. A coarser scan, or scanning only particles
flagged as uncertain, is the way to fit. Measure the loss-only cost first, then
choose the grid.

### 6.6 Only then, the fuller posterior treatment

Sections 8 and 9 describe credible regions, posterior-derived search bounds and
fitted population widths. These are worth doing but they are not the shortest
path to the objective, and they depend on calibration work that depends in turn
on the section 5 gate. Do not start them before 6.1 to 6.4 have been measured.

## 7. How to avoid fooling ourselves

Three recurring traps, all of which have already cost us time.

**Measure whether a change does anything before measuring whether it helps.** The
`cont_prior` penalty was compared across a sixteen-fold range of widths and every
arm returned the same FSC to three digits. That is the signature of a term that
is being added to something much larger. Before any comparison campaign, confirm
the mechanism under test can change the answer at the settings used: report the
penalty actually applied, the width actually used and the particle power
alongside the change in data loss, and include one deliberately extreme setting
as a positive control. A build that compiles and tests that pass prove neither.

**Do not let a change be scored on the quantity it suppresses.** Our convergence
rule declares a particle stable when it stops moving, and several of the
mechanisms here change how far particles move. A change that makes everything
stop early will look like fast convergence and may be a worse map. Always compare
at a matched iteration budget as well as at natural stopping, and report the
iteration count at convergence alongside the quality.

**Use the pairing.** Our arms share particles, seeds, references and starting
poses. Comparing one FSC number per arm throws away thousands of paired
observations and gives no uncertainty at all — which is why five arms differing
by a factor of sixteen all reported 3.33 Å and we could not tell whether that
meant "no effect" or "no measurement". Compare per-particle differences between
arms with a paired bootstrap or signed-rank analysis. That removes
between-particle variance and produces a confidence interval, which is what turns
"no difference" into "no difference larger than X" — the statement needed to
close a question. Report shell-wise FSC differences, not only the crossing, and
on synthetic data report the upper quantiles of pose error, since the tails are
where the failures live.

Fix calibration and equivalence margins before the runs that judge them, not
after.

## 8. Correctness items to fix along the way

These came out of a code review on 2026-11-11. None blocks section 6; fold each
into whatever change next touches that code.

**The rotation penalty's curvature is documented wrongly.** For a proper
rotation, the Gauss-Newton block of the Frobenius-residual rotation penalty is
exactly `2*w*I` everywhere, not only near the seed. The comparison table in the
done record rejects an alternative for having constant curvature, which is the
property the chosen form also has. Derive and test the penalty, its gradient, its
exact Hessian and the Gauss-Newton block separately, then correct the
description. Before replacing the computed loop with the constant, check the real
input contract: stored matrices come from single-precision Euler angles, so they
are not exactly orthogonal, and the identity holds exactly only for proper
rotations.

**The identifiability check does not check identifiability.** A successful
Cholesky factorization with a relative pivot floor establishes that the system is
numerically solvable, not that the data determine the pose. A badly conditioned
or uniformly weak system passes, and the prior then quietly rescues a direction
the data do not constrain, while the code reports a data-identifiable result.
Decide whether to narrow the claim to "numerically solvable" or add a genuine
weak-data rejection. Note that a condition number alone will not catch a
well-conditioned but uniformly uninformative system; both failures need covering.

**The stationarity test is an absolute threshold on a scale-dependent quantity.**
Convergence is declared at a gradient norm below `1e-8`, but that gradient's
scale depends on band, particle power and the rotation scaling. Build fixtures
that change the scaling while preserving the physical problem, and define a
relative test. This matters more once the search domain widens.

**Self-conjugate frequency bins.** The documentation discusses the zero-frequency
bin alone. The general case is bins where `k = -k (mod N)` in both coordinates;
enumerate them, confirm how they are stored, counted and included in the band,
and state the convention. Entire Nyquist rows are not individually
self-conjugate — only the fixed points are. This is a documentation correction
unless the audit finds the band actually includes one.

**Minor.** `prepare_pose_data_system` performs a factorization whose result is
discarded, using it only as a guard flag before the real combined solve; it reads
like a bug and needs a comment. Dummy argument alignment regressed in the two
prior setters, and three renamed declarations carry stray spacing
(`simple_parameters::l_cont_followup`, `simple_convergence::calc_cont_pose_stats`,
`simple_commanders_refine3D::exec_refine3D_workflow`).

## 9. Later: honest uncertainty, and using it

Once the Cartesian route produces several candidate poses with weights
(section 6.4), it becomes possible to ask what its uncertainty actually is and to
use that to make the search cheaper. That is the long-term prize: a search that
spends effort where the pose is uncertain and skips particles that are already
pinned down. The design below is carried forward and is **not implemented**.

**One explicit density.** Work in physical negative log likelihood,
`E(R,s) = T*L(R,s) + U(R,s)`, with `T` the conversion from normalized loss and
`U` the penalty. The two routes must agree on `E` at identical poses — this is
the section 5 gate, and it blocks everything here.

**Credible regions, not fixed multiples of a width.** For a candidate set with
weights `q`, the smallest cone retaining `1-epsilon` of the mass is the honest
bound. Three cautions. The weights must account for the volume each candidate
represents, so non-uniform grids cannot be counted equally. Several separate 99%
bounds on rotation, in-plane angle and shift do not give 99% joint coverage;
budget the lost mass jointly. And a quantile of a posterior that has already been
cropped to a search domain only describes mass inside that crop — it cannot
justify cropping further, so start broad and check the boundary.

**Beware tails estimated from few candidates.** A 99% bound is set by the least
likely candidates, exactly where a small weighted sample is least reliable. Gate
it on how many candidates effectively carry the weight, and fall back to a broad
search when that number is too small; otherwise the tightest bounds get assigned
to the worst-estimated particles, which loses them irrecoverably. There is an
existing entropy-based readout (`prob_neff_readout`) to build on, though it
measures concentration rather than sampling adequacy and the distinction needs
care.

**Uncertainty from curvature, if the posterior is a single narrow peak.** The
local information matrix gives an approximate covariance. Two requirements:
build it from likelihood and prior only, never from the optimizer's damping or
its empirical curvature scaling; and marginalize properly — rotation uncertainty
with shifts integrated out is the Schur complement
`(H_RR - H_Rs H_ss^-1 H_sR)^-1`, not `H_RR^-1`. Radii scale with the square root
of a chi-square quantile whose degrees of freedom differ by coordinate: two for a
direction cone, three for a full rotation, two for shifts. A fixed multiple of
the width is not the same retained mass in those different spaces.

**Fitted population widths, carefully.** Widths can be fitted across particles
rather than set by hand, using completed-pass statistics within each half map.
Two hazards. Reusing a particle's own posterior as its next prior and then
multiplying by the same image again counts the same evidence twice; refinement
iterations are passes of one optimization, not independent measurements. And the
fit has a degenerate direction: a tight width shrinks the posterior moments,
which shrinks the next width. Record the fitted width per pass per half and stop
on repeated decrease without held-out improvement.

**Symmetry.** When distances are minimized over the point group, the reference
distribution is not the unrestricted one. State which space every calibration
claim is measured in and validate against enumeration for both C1 and a
non-trivial group.

**Calibration testing.** Coverage at two or three stated levels is a blunt
instrument. On synthetic data, record the rank of the true pose within the
posterior sample: uniform ranks mean honest uncertainty, a U shape means too
narrow, a dome means too wide, a tilt means bias. Verify the expected rank
distribution on an exact toy model before trusting it on real output.

See also the [polar Laplace-width proposal](polar_shift_laplace_width.md), which
separates a search domain from a probabilistic width on the polar side. Its
`shsig` is an RMS radius in native pixels; the Cartesian
`cont_prior_shift_sigma` is a per-axis standard deviation in Angstroms. Convert
before comparing.

## 10. Current behaviour worth knowing before you touch the code

- `cont_prior=none|shift|pose` applies to single-state Euclidean joint Cartesian
  refinement; `none` is the default. The penalty is re-centred on the stored pose
  at the start of **each** optimizer transaction, so it restrains the size of one
  update and does not accumulate across passes. It is a step restraint, not a
  belief about where the pose is.
- A positive `cont_prior_shift_sigma` is a per-axis width in **Angstroms**,
  converted internally to native pixels by `smpd`. Zero selects an estimator
  built from all active particles' stored shift increments, with a floor.
  `cont_prior_rot_sigma` is a per-axis small-angle width in degrees.
- That estimator is driven by the increments the penalty suppresses, so it can
  tighten itself. The floor is what stops it. Treat the floor as load-bearing.
- Acceptance and optimizer stopping use the penalized objective; published scores
  and noise residuals remain data-only.
- The Cartesian solve returns a local penalized optimum, not a posterior.
- Cartesian bounds (full rotation, shift disk) and polar bounds (view direction,
  shift box per axis) must be matched explicitly before any comparison.

## 11. Ownership

| Source | Responsibility in this plan |
| --- | --- |
| `src/main/cftc/simple_cartft_pose_opt.f90`, `simple_cartft_pose_opt_system.f90`, sibling tester | 6.1 bounds split, 6.2 draw, 6.4 multi-start; section 8 curvature, identifiability, stationarity |
| `src/main/cftc/simple_cartft_calc.f90`, sibling tester | 6.5 loss-only evaluation path; section 5 NLL gate; self-conjugate bins; particle power reporting |
| `src/main/strategies/search/simple_strategy3D_cont.f90`, `simple_strategy3D_alloc.f90`, testers | per-particle configuration, width lifecycle, candidate handling, per-pass diagnostics |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | reused posterior, temperature, top-k and readout helpers; keep Cartesian policy out of the polar table module as it grows |
| `src/main/strategies/search/simple_matcher_3Drec.f90` | 6.3 weighted multi-pose insertion |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, matcher | pass lifecycle, part-to-master handoff, deterministic merge |
| `src/main/simple_convergence.f90`, tester | stopping rule, coverage, natural-convergence diagnostics |
| `src/main/params`, `src/main/ui/simple/simple_ui_refine3D.f90`, policies | controls, units, help text |

Keep numerical terms in `cftc`, configuration and batching in the search
strategies, stopping in its owners, and merging on the master. Workers consume
the master's frozen per-half model; they never estimate it from partition-local
or partially updated data. No project writes in workers, no test-only public
helpers, no per-file compile flags.

Build incrementally, run the fast gate and the touched suites, then the affected
workflow entries. Numerical changes need independently derived expected values.
Anything touching parallel code needs a test at the machine's full thread count,
repeated — a single-thread comparison never demonstrates thread safety.

## 12. Still open, not scheduled

- Whether all-active increment widths or a completed-sample estimator is the
  better width policy; decide from multi-dataset held-out evidence, after
  section 6.
- Single-precision boundaries recorded in the done archive: whitening and
  storage, the interpolation gather, phase constants, pose publication, sigma and
  score narrowing. Measure range and sensitivity separately.
- `objfun=cc` support, which needs a CC-compatible objective and prior scale
  derived first; Euclidean power scaling does not carry over.
- Retiring the internal `cont_route=shift_then_joint` route.
- Qualification on Linux and BOX, Release build, and the intended scheduler, when
  implementation resumes. Historical passing output in the done archive is not
  evidence for new code.
