# Pose marginalization in the probabilistic 3D search

Implementation note, 2026-10-21. **Design, with the first diagnostic
implemented** (`prob_post_diag=yes`, this batch; nothing in the reconstruction
changes). Written after the aggressive-frequency-marching experiment
(`freq_march_aggressive`, `solve3D_policy.md` section 4): HR-HAIR's
high-resolution start works inside cryoSPARC's SGD with pose marginalization
over mini-batches, and PfCRT `latest10` showed it does not transfer to a
hard-assignment stochastic hill climb from random orientations. This note
asks what marginalization itself would look like in SIMPLE, and how to find
out whether it is worth building.

## 1. The question

For a small membrane protein the pose likelihood is expected to be broad or
multimodal: several views through the micelle fit nearly equally. A hard
assignment then commits each particle to one of them and the reconstruction
averages the commitments' errors; marginalization lets each particle enter
the reconstruction at every plausible pose with its posterior weight, so a
particle that cannot decide does not vote for one wrong answer. cryoSPARC
offers it as an option ("adaptive marginalization") and reports it as a
secondary help for small, low-SNR targets, behind non-uniform regularization
and masking. Whether it helps *here* is an empirical question with a cheap
answer: measure how broad the posterior actually is on PfCRT.

## 2. What exists

- **The table holds the posterior's support.** In the `prob` stages (3-5)
  each table worker fills, for its particles, a full row
  `dists(iptcl, 1:nrefs)` over every projection direction (and state) at a
  likelihood-sampled in-plane rotation, then shift-refines the `projs_ns`
  best into the per-particle `refined` list, whose refined distances
  overwrite the row (`fill_tab_range`, `refine_best_refs`,
  `simple_eul_prob_tab.f90`). The rows are written to the part file and read
  by the master into `dists(nptcls, nrefs)`. The neighbourhood tables of the
  `prob_neigh` stages (6-8) hold the same over a neighbourhood.
- **The assignment is not a per-particle draw.** `ref_assign` is a
  reference-driven balanced matching: every reference's particles sorted by
  distance, the frontier heads drawn by `draw_likelihood`, one
  (reference, particle) pair committed per step. `strategy3D_prob` carries
  that single pose. The draw's weights, `exp(-(L - L_min))` on the normalized
  loss, are nearly flat (section 3), so the draw is in effect uniform over
  the `prob_athres` window of frontier heads: an exploration device, not a
  posterior.
- **The reconstruction already takes weighted multi-insertion.**
  `calc_3Drec(..., wset, wtab)` lets a particle enter several states with
  weights; `update_state_half_rec` inserts `fplanes(i)` through
  `insert_plane_oversamp(..., w=)` per entry (`simple_matcher_3Drec.f90`).
  The PCG backend takes records (`take_records`, `get_rec_observation`).
- **sigma2 is estimated at the committed pose.** `calc_sigma2_pftc` writes
  the residual at the assigned orientation into the particle's record
  (`gen_sigma_contrib`: the mean squared residual over two per ring, the
  variance per real or imaginary part).

Marginalization is therefore not a new search. It is: keep the top of the
row with weights instead of drawing one entry, and hand the weighted list to
the two consumers that should see it.

## 3. The posterior and its temperature

The row holds the normalized Euclidean loss (`gen_raw_euclid_vals_impl`)

    L_r = sum_k (k/sigma2_k) sum_p |X_kp - CTF_kp R_kp,r|^2 / wsqsum(X),
    wsqsum(X) = sum_k (k/sigma2_k) sum_p |X_kp|^2,

summed over the `pftsz` polar samples of each ring k in the band. The
Gaussian negative log-likelihood runs over the Cartesian half-plane
components, pi*k of them in ring k, so each polar sample stands for
pi*k/pftsz of them; with sigma2 the variance per real or imaginary part,

    NLL_r = sum_k (pi k/pftsz) sum_p |r_kp|^2 / (2 sigma2_k)
          = pi wsqsum(X) L_r / (2 pftsz).

The posterior over the row is p_r proportional to exp(-T (L_r - L_min)) with

    T = pi wsqsum(X) / (2 pftsz)        (posterior_temperature).

L lives in [0, ~2] and wsqsum is of order 2 pftsz sum_k k (1 + SNR_k), so T
is of order pi sum_k k, about 1.5 K^2 for K shells in the band: roughly a
thousand at the first prob stage. A loss difference of 0.01 between two
directions is then ten nats; the posterior is sharp whenever the data can
tell the directions apart and broad only where they fit within about 0.001
in L. Two consequences:

1. The current draw at T = 1 is not a posterior and should not be read as
   one. This is a finding about the existing sampler, not a defect: it is
   doing exploration within an angular window.
2. The posterior self-adapts: concentrated for a high-SNR particle, broad
   for a weak one, which is the behaviour marginalization needs.

The constant pi/(2 pftsz) assumes the pft is taken at integer ring radius on
the cropped grid (one Cartesian component per unit area) and that the
particle-group sigma2 is per real or imaginary part. Both can be checked on
the diagnostic output: at the mode, T L_min is the NLL, which for a calibrated
sigma2 should be close to the number of components in the band, about
pi (k2^2 - k1^2) / 2; and on bgal at the fine stages the effective number of
references should sit near one.

## 4. The design

**M-step marginalization, direction only, `prob` stages, gridding backend
first.** The pose carried forward stays the assignment's; the weights change
only what the reconstruction and sigma2 see.

1. *Weights.* On the worker, where the full row and `wsqsum` are, the
   posterior over the row at temperature T; the candidates holding 99% of
   the mass, capped at `projs_ns` (they are the refined ones), with their
   normalized weights. The `refined` entries already carry (reference,
   shift) per particle; they gain the weight and the in-plane angle, and the
   assignment file carries the list beside the assigned pose.
2. *Reconstruction.* `prep_imgs4rec` prepares the particle's plane once;
   `update_state_half_rec` inserts it at each candidate (direction,
   in-plane, shift) with its weight instead of once. The shift is a phase
   ramp applied in `gen_fplane4rec`, so candidates with different refined
   shifts need the ramp re-applied per candidate (cheap) rather than a
   shared shift (wrong exactly for the far-apart views the micelle
   ambiguity produces). Gridding: k insertions per particle, a small cost
   beside the search. PCG: k records per particle and every CG iteration
   visits every record, so the solve costs k times as much; start on
   gridding, or cap k at 2-4 on PCG.
3. *sigma2.* The posterior-expected residual, sum_r p_r residual_r (k calls
   to `gen_sigma_contrib`, accumulated), in place of the residual at the
   committed pose. Without this the estimate stays biased low (the winner's
   residual is the smallest), the next posterior is too sharp, and the
   scheme collapses back to hard assignment within a few iterations.
4. *Unchanged.* The balanced assignment and the pose it commits, the
   neighbourhood tables, the overlap statistic and convergence, in-plane
   and shift at their maxima (the ambiguity that hurts a small membrane
   protein is between directions; in-plane is sharp), the `prob_neigh`
   stages (their neighbourhood tables give a truncated posterior, fine
   later), multi-state (weights within the particle's assigned state first).

## 5. The diagnostic (implemented)

`prob_post_diag=yes` (internal key, default `no`, `objfun=euclid` only) makes
every dense table worker compute, per particle, once its row is complete:

| column | meaning |
| --- | --- |
| `iptcl` | particle index |
| `istate`, `iproj` | the mode (lowest distance) |
| `nvalid` | valid distances in the row (0: no mode, zeros follow) |
| `lmin` | the mode's normalized loss |
| `temp` | T = pi wsqsum / (2 pftsz) |
| `pmax` | the mode's posterior mass |
| `neff` | exp(entropy): the effective number of references |
| `k99` | references holding 99% of the mass |
| `spread_deg` | posterior-weighted mean angular distance of those references from the mode, in the mode's state, across the symmetry operations (`sym_dists`) |

and, since 2026-10-22, the same row scored against the *other* half's
reference (`gen_best_euclid_val_cross`: the other half's `pfts_refs`, the
particle's own noise weights and self-term, the seed shift, best in-plane
rotation, no refinement), which the particle did not enter:

| column | meaning |
| --- | --- |
| `iproj_x`, `lmin_x` | the cross-half mode and its loss |
| `pmax_x`, `neff_x`, `k99_x`, `spread_x_deg` | as above, for the cross-half row |
| `dist_x_deg` | angular distance between the own and the cross-half modes (symmetry-aware) |
| `rank_x` | the own mode's rank in the cross-half row (1: the halves agree) |
| `bonus_nats` | T x (loss at the own mode, its refined shift, against the other half minus against its own): the self-reinforcement bonus |

written as `prob_posterior_iterNNN_partNN.txt` in the run directory, one
file per part and iteration, a header line first (`write_posterior_diag`).
The cross-half row doubles the scoring cost of the table in diagnostic
mode. `test_cross_half_euclid` (suite "polar correlations") pins that the
switch flips the reference half alone.
`posterior_temperature` and `posterior_row_stats` live in
`simple_eul_prob_tab_utils.f90` and are pinned by `test_posterior_row_stats`
(`unit_...` suite "3D probability table").

The cost is one exponential and a sort per row, negligible beside the
scoring; the files are one line per sampled particle and iteration.

**How to read it.** `neff` and `k99` say how many references the data
cannot tell apart; `spread_deg` says whether they are neighbours (a smooth
peak sampled by several grid points: benign, the posterior is a kernel
around one view) or far apart (true multimodality: the case marginalization
is for). Per stage, on PfCRT against bgal:

- `neff` near 1 and `spread` small: the data decide the pose; marginalization
  changes nothing, and the small-particle problem is elsewhere (the
  reference, the band, the shifts).
- `neff` large with `spread` within a few degrees: finite grid sampling of
  a broad but unimodal likelihood; marginalization smooths, and a finer
  `nspace` or a continuous polish (`refine=cont`) does the same job.
- `neff` large with `spread` of tens of degrees in a sizeable fraction of
  the particles: multimodal; build section 4.

Also worth reading off: `temp * lmin` against the component count, the
calibration check of section 3, and the trajectory of `neff` over the
iterations of a stage (does the posterior sharpen as the reference improves,
or does the alignment churn?).

## 6. What has to be true first

- The temperature constant and the sigma2 convention must be confirmed on
  the calibration check before any weight enters a reconstruction; a wrong
  constant makes everything look ambiguous (too cold) or decided (too hot).
- The balanced assignment and a per-particle posterior are different
  objects; the balancing exists to keep the reference population from
  collapsing early. Marginalizing in the M-step keeps the balancing for the
  carried pose and only redistributes the reconstruction's votes; whether
  the assigned candidate needs a floor weight is an open question.
- Marginalization needs a calibrated noise model; cryoSPARC leaves it out of
  local refinement for exactly that reason (a wrong sigma2 under micelles
  and crowding can make it worse). The sigma2 expectation in section 4 is
  the minimum.

## 7. Test plan

1. `prob_post_diag=yes` on the ten-restart PfCRT regression and on bgal, at
   the default schedule; per stage and iteration, the distributions of
   `neff`, `k99`, `spread_deg`, and the calibration ratio. This decides
   whether section 4 is built.
2. If built: section 4 on the gridding backend, `prob` stages, scored as
   the frequency-marching experiment (stage-5 FSC, iterations to the stage-6
   Nyquist, the gold-standard bootstrap of the final poses), and the
   per-iteration `neff` trajectory, which should sharpen as the map
   improves.

## 8. Results: PfCRT, 2026-10-22

Run `PfCRT_posterior/4_solve3D` (default schedule, 16,832 particles, full
sampling, PCG, `prob_post_diag=yes`), the dense stages: iterations 41-57 are
stage 3 (18.3 A, nspace 1000), 58-72 stage 4 (12.9 A, 1000), 73-78 stage 5
(NU, 2500; early exit); stage 6 omitted, 7-8 `prob_neigh` (no table). The
run's `spread_deg` is 57.2958x too large (a `rad2deg` applied to degrees;
fixed in the code the same day) and was divided out in the analysis.

**Calibration.** `temp * lmin / (pi k2^2 / 2)` = 0.99 through stages 3-4
(k2 from the stage crop and limit). The constant pi/(2 pftsz) and the
sigma2 convention of section 3 are right; the posterior is in nats. (Stage
5 reads 1.8-3.4 only because the NU band is finer than the planned 11.1 A
used for k2.)

**Given the reference of the iteration, the posterior is sharp.** Median
n_eff 2.9 at the first prob iteration, 1.3 by the fifth, 1.1 through stage
4; median p_max 0.93-0.99; k99 median 2. Broad particles (n_eff > 2):
22-30% in stage 3, 13-20% in stage 4, 18-38% in stage 5 (the band jump at
the NU handoff re-broadens it for three iterations). Among the broad, the
spread is mostly one grid step (median 4-9 A against 7.1 deg sampling);
the far tail (> 20 deg) is 5-8% of all particles. Broadness is not a
particle property: the overlap of the broad sets at a stage's first and
last iteration is at chance, and the per-particle mean n_eff correlates
0.2-0.3 between halves of a stage. It is weakly tied to signal (broad
fraction 31% in the lowest wsqsum quartile, 26% in the highest).

**But the mode is not stable across iterations.** Between consecutive
iterations 25-36% of all particles move their mode by more than 20 deg,
and 16-32% of the particles whose posterior put at least 90% of its mass
on one direction do; the 90th percentile of their move is about 100 deg.
Over stage 3, 76% of modes end more than 20 deg from where they started
(median 68 deg); stage 4, 58%; stage 5, 43%. The far moves have the angle
distribution of random direction pairs (no pseudo-symmetry flip), and the
movers' best loss is worse after the move by 2-3 nats.

**Reading.** A posterior that is sharp against the current reference and
flips by 100 deg between iterations is not measuring the data's
information about the pose; it is measuring the reference's memory of
where the particle was last put. With about 8 particles per direction per
half-map at nspace 1000, a particle's own noise in the reference at its
assigned direction is a large part of that projection, worth tens of nats
at T of order 1,000, against the few nats that separate the data's true
alternatives. The balanced stochastic assignment (uniform over the top
~140 frontier heads, section 2) then moves the particle, the reference
follows, and the next posterior points at the new place. The stages behave
as a very hot sampler over a very cold, self-reinforced likelihood; neither
the hard assignment nor a marginalization with these weights measures what
the data say.

**Consequence for section 4.** Weights from the own-half row would
entrench the previous assignment (p_max ~0.95 on the self-reinforced
direction) and help nothing. The posterior has to be taken against a
reference the particle did not enter. The decisive measurement is the
cross-half row (section 5, implemented 2026-10-22, to be re-run): the
cross-half n_eff, the own mode's rank in the cross row, the distance
between the two modes, and the bonus. Expected readings:

- `bonus_nats` large (tens) and `f(modes agree)` low, `neff_x` broad:
  self-reinforcement confirmed and quantified; the cross-half posterior
  is the honest one and the one to marginalize with (for the search it
  breaks the half independence, so it is a solve3D-stage device, not a
  gold-standard one), and the own-half posterior is not a measurement.
- `bonus_nats` small (a few) and the halves' modes agree while the mode
  still jumps between iterations: the posterior is honest and the jumps
  come from the assignment alone (the frontier draw at T = 1); the remedy
  is in `ref_assign`, and marginalization with these weights is
  well-founded but will not change much (n_eff near 1).
- `neff_x` broad with `dist_x_deg` far for a stable subset of particles:
  real multimodality after all, section 4 as designed.
