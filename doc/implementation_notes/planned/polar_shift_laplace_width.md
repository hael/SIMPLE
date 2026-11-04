# The Laplace width: the dynamic centring's prior width from the full band

Planned, not built (Hans's call). Separated on 2026-11-03 from section 2.6 of the dynamic-centring proposal
(`../completed/polar_shift_search_dynamic_centring.md`), when the polar shift search closed
(`../completed/polar_shift_search.md`, ruling 19 and the progress entries of 2026-11-02 and 2026-11-03). The
review's qualifications are from `../completed/polar_shift_search_review.md`, section 14.

## Why

The dynamic centring keeps two numbers per particle: `shincr`, the last move, and `shsig`, the prior width of
its next search. A particle takes the disk when `shsig` is the sentinel or when three widths exceed the window
(`shift_needs_disk`); a move that reached the box's bound inflates the width to the radius searched
(`next_shift_width`). The width is measured on the coarse map, at the physical temperature, over the cells the
screen evaluated (`shiftmap_draw` with `temp_phys`).

On data that width reads 0.00 for nearly every particle from the second iteration on (beta-gal: 90th
percentiles of 0.27 px early, 0.00 later). It is a resolution limit of the measurement, not of the particle:
the coarse band ends at shell `kscreen`, so the map's radial step is `box / (4 kscreen)`, 3.6 px at beta-gal's
stage 4 against a window of 4.4 px, and at the physical temperature the loss differences between cells are
tens of nats, so all the mass sits on one cell. A read of 0.00 means "surer than one cell". The rule "three
widths beyond the window" needs widths resolved below about 1.5 px, which the coarse map cannot give, so today
the width is a two-state flag (multi-modal at the coarse band or not) and the re-centring is carried by the
move at the bound. Two changes within the map do not help: the draw's temperature (0.3 times the physical)
inflates the widths with the exploration's flatness, not the particle's uncertainty; a floor at the cell's size
makes every particle wide against a window of one cell.

## The proposal

The width at the polished optimum on the full band, from the curvature of the loss there (the Laplace
approximation): the standard deviation of a Gaussian fitted to the posterior at the committed pose, resolved
to a fraction of a pixel. A well aligned particle reads a few tenths of a pixel, a weak or wrongly classed one
pixels, and the rule then separates those rather than the two states the coarse map can tell apart.

- Measured by the consumer that commits the pose, after the continuous in-plane polish, for the committed pose
  alone: the Hessian H of T L (L the loss, T the particle's physical temperature) from a stencil of full-band
  evaluations about the optimum.
- The width: the square root of the larger eigenvalue of H^-1 (the conservative choice for a box), or of the
  trace of H^-1 / 2.
- The stored prior width is the Laplace width, with the move's inflation as now. The coarse map's width stays
  beside it as the multi-modality flag: a particle with mass on two coarse cells takes the disk whatever its
  local curvature says.
- The width never enters the likelihood. The posterior re-used as a prior on the same image would count the data
  twice; the map's prior (on the increment from the stored shift, ruling 22) stays the only prior (in the first
  centring's coarse draw and in the insertion's loss),
  and the width chooses the domain of the search, not the weight of its cells.

## Qualifications (the review, section 14)

1. **The angle.** The committed pose after the continuous in-plane polish has a fractional angle, and
   `gen_corr_for_rot_8` takes an integer one. Evaluate at the committed fractional angle; the calculator has
   fractional-angle evaluators.
2. **Conditional or marginal.** A 2 x 2 shift Hessian at a fixed angle gives the conditional covariance, while
   the coarse moments marginalise over the in-plane rotation. With appreciable shift-angle coupling these are
   different criteria: for a joint Hessian partitioned into shift and angle blocks the marginal shift
   covariance is `(Hss - Hsa Haa^-1 Has)^-1`. For the positive-definite `H = [[1, 0, .99], [0, 1, 0],
   [.99, 0, 1]]` the largest conditional standard deviation is 1 and the marginal one 7.09; a window of 5 is
   chosen by one rule and refused by the other. Either marginalise the coupled local model (the 3 x 3 Hessian
   with the angle) or adopt the conditional rule and calibrate it.
3. **Flat or indefinite curvature.** A Hessian that is not positive definite (a saddle, a flat direction, a
   polish that stopped at the box's bound) has no Laplace width; treat it conservatively, as the radius
   searched.
4. **Where the numerics live.** The stencil and the Hessian belong to the calculator (`src/main/pftc`) or the
   search domain; the strategies orchestrate.

## The width's convention

The coarse width `shiftmap_draw` measures is the RMS radius of the shift about its mean, the square root of the
covariance's trace: for an isotropic Gaussian of per-axis sigma s it reads sqrt(2) s. The rule "three widths
beyond the window" (`shift_needs_disk`: 3 shsig > trs) is therefore 3 sqrt(2) = 4.24 per-axis sigmas, more
conservative than the disk's three per-axis sigmas of the prior (98.89 % of an isotropic Gaussian's mass). The
Laplace width as proposed, the square root of the larger eigenvalue of H^-1, is a per-axis sigma along the
principal axis: adopting it unchanged makes the same rule 3 per-axis sigmas, sqrt(2) less conservative than
today. Either keep the RMS convention (the square root of the trace of H^-1) or state the factor with the rule.

## Decisions (Hans)

1. Conditional (2 x 2) or marginal (3 x 3 with the angle) width.
2. The larger eigenvalue or the trace (the trace keeps today's RMS convention, above).
3. What an indefinite Hessian means for the next box (the radius searched is the proposal).
4. The stencil's step (a fraction of the map step) and its cost budget per particle.

## Where it plugs in (as of 2026-11-03)

- The committing consumers, after the continuous route: `strategy2D_prob` and `strategy3D_prob` (the tables),
  the classic 2D `assign_ori` and the classic 3D `record_shift` (from `extract_peak_ori`).
- `record_shift_state` takes the Laplace width as `sdev` and a multi-modality flag from the coarse map.
- The evaluator in the calculator's correlation submodule, beside the fractional-angle objectives.
- Nothing in the tables' fills changes.

## Validation

- The widths' distribution on data against the movers (the disk set the drivers print).
- The fixtures with sampling (solve2D and solve3D a quarter sampled per iteration), where particles enter late.
- The dynamic-state regressions the review asked for (mixed sentinel, moving and settled particles across the
  stages, sampling and injection), deferred when the polar shift search closed; they belong with this change.
