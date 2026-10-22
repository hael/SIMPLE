# Pose marginalization in the probabilistic 3D search

Implementation note, 2026-10-21, revised 2026-10-22 twice: after the
measurements (sections 8-9) and after the review
(`pose_marginalization_review.md`, section 11). **Design, not yet
justified; the diagnostic is implemented (`prob_post_diag=yes`), the
design is not, and the review gate of section 11 comes before it.**
Sections 2-3 describe what exists and the posterior; section 4 is the
design; sections 5-7 are the diagnostic, the preconditions and the test
plan; sections 8-9 the measurements, with section 9's inference
withdrawn in section 11. Written after the
aggressive-frequency-marching experiment (`freq_march_aggressive`,
`solve3D_policy.md` section 4): HR-HAIR's high-resolution start works
inside cryoSPARC's SGD with pose marginalization over mini-batches, and
PfCRT `latest10` showed it does not transfer to a hard-assignment
stochastic hill climb from random orientations. This note asks what
marginalization itself would look like in SIMPLE, and whether it is worth
building.

**The answer so far (sections 8-9, as corrected in 11).** Given the
reference of the iteration, the posterior over the direction grid is
sharp (n_eff near 1 for the median particle), so there is no broad
posterior to spread; whether that sharpness is calibrated against the
data (and not an artefact of the Gaussian model, the correlated polar
samples or the nuisance search) is open until the known-pose test of
section 11. The assignment's draw ignores that posterior: at T = 1 on the
normalized loss it is near-uniform over the ~140 best frontier heads (a
code fact, section 2), and the production row's mode is unstable between
iterations where a deterministic re-scoring of the same reference is
stable (section 9, corrected). The proposed fix, drawing the committed
pose from the row posterior at its temperature (section 4), is therefore
a hypothesis with a mechanism behind it, not a demonstrated remedy: the
"cross-half" rows of section 9 scored the same merged reference, so the
self-reinforcement chain written there is unproven.

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

## 4. The design (revised 2026-10-22, Monte Carlo form)

**Draw the committed pose from the posterior; let the update scheme do the
averaging.** SIMPLE is already a Monte Carlo marginalizer by
architecture: one pose per particle per iteration, drawn stochastically,
averaged over particles within an iteration and, when the fractional
update is on, over iterations by the trailing chain. Inserting k planes
with posterior weights and inserting one plane drawn from the posterior
give the same reconstruction *accumulators* in expectation, conditional on
fixed weights; the restored map is a ratio of accumulators, so its
expectation is not exactly the weighted map (review finding 4), and the
weighted form removes the within-iteration variance of the draw. With
n_eff 1.05-1.3 and 95% of the mass on the mode for most particles
(sections 8-9) that variance is small. Two qualifications: under forced
full sampling (`nsample` at or above 90% of the particles, the PfCRT
benchmark's mode) `set_refine3D_trailrec_policy` turns the trailing chain
off, so there is no averaging across iterations there; and the weights
change with the reference between iterations, so the chain does not
integrate one fixed posterior. The draw is also the cheaper estimator for
sigma2: the residual at a pose drawn from the posterior is unbiased over
the posterior, where the argmax residual is biased low.

What is wrong today, in these terms, is the draw's distribution, not the
number of planes: `draw_likelihood` draws with exp(-(L - L_min)) on the
normalized loss, where the posterior is exp(-T (L - L_min)) with T of
order 500-4000 (section 3; the constant checks to 1% on the mean residual
scale, section 8, which is not yet a validation of the likelihood ratios,
section 11). At T = 1 the draw is near-uniform over the ~140 best frontier
heads, so a sizeable fraction of particles is committed far from the row's
mode each iteration. That this harms the reconstruction is the hypothesis
of this note, not a measurement (section 11): what the diagnostic showed
is that the production row's mode is itself unstable between iterations.
At the posterior temperature the far-draw rate would fall to the
posterior's own far mass, about 5% on PfCRT, held by the particles whose
row leaves two views open.

1. *The draw.* On the worker, where the row and `wsqsum` are: the posterior
   over the row at the temperature of section 3, the candidates holding
   99% of the mass (capped at `projs_ns`, the refined ones, each with its
   refined shift and in-plane), and one draw from them. The draw is
   written as the particle's assignment; `strategy3D_prob` commits it as
   today, the reconstruction inserts it once as today. The carried pose is
   the draw, which is the mode for 95% of particles and a grid neighbour
   for most of the rest, so the `prob_neigh` stages inherit the data's
   answer. Nothing changes in the tables, the reconstruction, the
   fractional update or the trailing chain.
2. *sigma2.* Unchanged in code: the residual at the committed pose, now a
   posterior draw, is the right estimator in expectation.
3. *The balancing.* `ref_assign` is reference-driven: every reference draws
   its best available particle, which fills every direction. The
   per-particle draw drops that guarantee. In 3D with a thousand
   directions an empty direction is arguably the right answer when no
   particle's data put it there (the reconstruction handles missing views
   through neighbouring slices, and filling the direction with rejected
   particles is exactly the wrong vote of section 9), but whether the
   stage-3 consolidation (overlap 0.23 to 0.50) needs the balancing is a
   hypothesis the test must settle. The middle paths, if it does: keep
   the frontier structure and temper its draw, or cap the per-reference
   population and draw within the cap.
4. *Tempering.* The exact posterior temperature may be too cold in the
   first iterations of stage 3 (n_eff 2.9 at the first, 1.3 by the fifth):
   a factor on T rising to 1 over the first iterations of a stage is the
   annealing knob, in one constant.
5. *Unchanged.* In-plane and shift as the production row leaves them
   (sampled for the unrefined entries, continuously refined for the top
   `projs_ns`: the draw is an approximate direction sampler on that row,
   not a joint pose posterior; section 11), the `prob_neigh` stages at
   first (their neighbourhood draw can follow the same rule as a second
   step), multi-state (the draw within the particle's assigned state).

*The weighted form* (k planes per particle with their posterior weights,
the sigma2 expectation computed explicitly, k insertions on gridding and k
records on PCG) stays as the variance-reduction option, worth building
only if the per-iteration noise of the draw shows in stages 3-4, where
there is no trailing chain. *A colder frontier draw* (section 4.3's first
middle path) is the smaller code change but keeps the reference-driven
assignment; the per-particle draw is the one that follows from the data.

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

and, in the v3 of 2026-10-22 (after the review, section 11), the
controls on the same particle, each answering one question:

| column | meaning |
| --- | --- |
| `null_a_nats` | the identical-pose null test: the kernel re-scored at the production mode's exact (reference, shift, grid rotation) minus the stored loss, in nats, with the stored shift rotated back into the kernel's frame (`shift_ref2ptcl`; the search rotates it to the frame of reference before storing). Expected zero up to precision. |
| `null_b_nats` | the same with the stored shift taken as it is: the alternative convention, expected nonzero whenever the shift is |
| `iproj_d`, `lmin_d`, `pmax_d`, `neff_d`, `k99_d`, `spread_d_deg` | the deterministic row on the own slot: best grid rotation at the seed shift for every direction (`gen_best_objfun_val`), no sampling, no refinement, and its posterior |
| `dist_d_deg`, `rank_d` | the search test: the angle between the production mode and the deterministic mode, and the production mode's rank in the deterministic row (1: the same mode) |
| `iproj_h` ... `spread_h_deg` | the same deterministic row on the particle's **own real half** (the diagnostic model below; `gen_best_objfun_val` on the reference at `iref + nrefs`), and its posterior |
| `dist_h_deg`, `rank_h` | the angle between the matcher-row deterministic mode and the own-half mode, and the matcher mode's rank in the own-half row: what the merged reference and the own half disagree on |
| `iproj_x` ... `spread_x_deg` | the same deterministic row on the **other real half** (`gen_best_euclid_val_cross` at `iref + nrefs`), and its posterior |
| `dist_x_deg`, `rank_x` | the reference test: the angle between the own-half and other-half deterministic modes, and the own-half mode's rank in the other-half row (1: the two halves agree on the mode) |
| `iproj_p` | the pose the particle went into its half with: its previous direction on the current grid (`find_closest_proj`; exact within a stage), 0 when there is none (a virgin field, an invalid state) |
| `self_nats` | the self-reinforcement test: T (L_other - L_own) at that pose (previous direction, previous in-plane rotation, zero residual shift: the image is pre-shifted by its stored shift). Positive when the own half, which holds the particle, favours the pose it was inserted at. The pose was chosen on the previous merged reference, to which both halves contributed alike, so the selection bias is the same for both halves and the mean over particles isolates the particle's own contribution. |
| `ctrl_nats` | the same difference at the same direction a half turn away in-plane, where the particle's own contribution does not correlate with it (a half turn is the Friedel mate on the polar grid: the cross term is Re X^2, zero in expectation): the noise scale of `self_nats`, expected zero on average |

Two more header lines say which case the file is in: `# slots_identical
yes|no` (`polarft_calc%ref_slots_identical`, the matcher's two slots;
solve3D's LP-set stages give `yes`) and `# halves yes|no` (the real halves
loaded beside the matcher's references; `no` leaves the `_h`, `_x` and
self-test columns zero). The columns are written as
`prob_posterior_iterNNN_partNN.txt` in the run directory, one file per
part and iteration (`write_posterior_diag`).

**The diagnostic model of the real halves.** The matcher's references in
solve3D's LP-set stages are one merged map in both slots (section 11,
finding 1), so the master prepares a second polar model for the
diagnostic: after `materialize_reprojection_model_from_volumes` has written
`reprojection_model_{even,odd}.bin` it calls
`read_mask_filter_reproject_refvols(..., map_shift=.false.,
diag_halves=.true.)`, which runs `read_mask_filter_refvols` with the
independent `vol_even`/`vol_odd` maps of every state (the `_nu_filt` halves
under NU filtering, the raw halves otherwise) and the same masking,
centering and filtering as the matcher's references but with the merged-map
copies bypassed (the LP-set average, the merged NU reference, the ICM
average), and writes `reprojection_model_diag_{even,odd}.bin`. No shift
mapping, no band-limit estimation, no Cartesian volumes on that pass; the
files are removed with the model's by `remove_ref_section_files`. The dense
table workers (`exec_prob_tab`, `exec_prob_tab_neigh`) call
`read_reprojection_model(..., diag_halves=params%l_prob_post_diag)`: when
the pair is on disk the polar object is made with `2*nrefs` references, the
matcher's model loaded at 1 and the diagnostic pair at `nrefs+1`
(`read_ref_pfts(..., iref_from=)`, the header's reference count checked
against the range), and `eul_prob_tab` keeps `nrefs = nspace*nstates` for
the table, so the production search never sees the extra references. The
reference memory of a worker doubles in diagnostic mode (`mirr_proj=yes`
is excluded: its memoization assumes `nrefs/2`). The deterministic rows cost
three more scoring passes per particle, the self test four at-pose
evaluations.

`test_cross_half_euclid` (suite "polar correlations") pins that the slot
switch flips the reference half alone, that `gen_euclid_val_at`
reproduces the best value at its rotation, and that the slot detector
tells distinct sections from a copy; `test_ref_pfts_offset_io` (same suite)
that a model written by a two-reference object lands unchanged at
`iref_from` 1 and 3 of a four-reference one, both halves, and that the
plain read is the `iref_from` 1 case.
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

- The temperature constant and the sigma2 convention: the mean residual
  scale checks to 1% on PfCRT (section 8), so the algebra of section 3 is
  right; the likelihood ratios between poses are not validated by that
  (review finding 3), and the known-pose simulation of section 11 is the
  test.
- Reference independence of the diagnostic: not established in solve3D,
  whose staged children score one merged map in both slots (section 11).
  The posterior's honesty against an independent reference is unmeasured.
- Whether the stage-3 consolidation needs the balanced draw's exploration
  is not known. The balancing exists to keep the reference population from
  collapsing early; carrying the mode removes it from the pose. The test
  plan watches the stage-3 overlap trajectory and the cFAR for a
  preferred-orientation collapse, with the stage-3-only fallback of
  section 4.3.
- A calibrated noise model: the sigma2 expectation of section 4.4 is the
  minimum. cryoSPARC leaves marginalization out of local refinement
  because a wrong sigma2 under micelles and crowding makes it worse; here
  the weights are near-hard for most particles, so the exposure is small.

## 7. Test plan

1. Done: `prob_post_diag=yes` on PfCRT at the default schedule, the
   production row (section 8) and the deterministic re-scoring (section 9,
   as corrected). Next: diagnostic v3 and the review gate of section 11,
   then bgal as the control.
2. Build section 4.1 behind one key (`prob_assign=balanced|posterior`,
   default `balanced`: the per-particle posterior draw on the worker,
   written as the assignment), `prob` stages, with the diagnostic on. No
   change to the reconstruction.
3. The ten-restart PfCRT regression, `posterior` against `balanced`. The
   diagnostic itself is the first readout: for particles decided by both
   halves, the own-mode far-move rate should fall from 10-20% toward the
   cross-half's 0-1%, and the bonus tail (8-18% above 20 nats) should
   shrink, iteration by iteration within a stage. Then the stage-3 overlap
   trajectory and cFAR (the collapse check of section 4.3), the stage-5
   FSC, the iterations to the stage-6 Nyquist, the finals, and the
   gold-standard bootstrap of the final poses, scored as the
   frequency-marching experiment.
4. If the stage-3 consolidation suffers: the tempered frontier draw or the
   per-reference cap of section 4.3; if the per-iteration draw noise shows
   in stages 3-4, the weighted form.

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

**Reading at the time (superseded by section 9).** The first reading was
that the sharpness itself was self-reinforcement: the reference's memory
of where the particle was last put, worth tens of nats against the few
that separate the data's alternatives. The cross-half row tested that
reading directly and found the mechanism real but the magnitude wrong:
the bonus at the mode is a fraction of a nat for most particles, and the
data do decide the direction. What the first run had seen is the
assignment moving the mode, section 9.

The cross-half row was added to decide between three readings: a large
bonus with disagreeing halves (self-reinforcement dominant), a small bonus
with agreeing halves while the mode still jumps (the assignment moving the
mode), or a broad cross-half posterior with a stable far subset (real
multimodality). Section 9 gives the second.

## 9. Results of the cross-half re-run: PfCRT, 2026-10-22

Same data and schedule, iterations 41-77 (stage 5 through its fifth
iteration), the 19-column files. Analysis in `PfCRT_posterior/analysis/`
(`cross.py`, `cross2.py`, outputs `run2_*.txt`).

**The cross-half posterior is sharp too, and the halves mostly agree.**
Against the reference the particle did not enter: median n_eff_x 1.05
(stage 3), 1.00 (stage 4), 1.04 (stage 5); median p_max_x 0.99-1.00; the
two halves' modes coincide (rank_x = 1) for 61-67% of particles in stage 3,
48-59% in stage 4, and lie more than 20 deg apart for 15-16% (stage 3-4),
9-12% (stage 5). The self-reinforcement bonus at the own mode is small for
most particles, median 0.3 nats (stage 3), 1 (stage 4), 2 (stage 5), with
a tail: 8-18% above 20 nats, p90 15-38 nats, growing with the band.

**The cross-half mode is stable; the own-half mode is not.** Between
consecutive iterations the cross-half mode moves more than 20 deg for
18-21% of particles (own: 25-30%), and its far-move rate falls with its
confidence (45% below p_max_x 0.5, 11-19% above 0.999) as a data effect
should. The own-half far-move rate is flat across cross-half confidence
(22-31% in every p_max_x bin). The decisive subset: particles decided by
both halves on the same direction (p_max >= 0.99 own and cross, rank_x = 1;
16-27% of all). Their cross-half mode stays within 20 deg next iteration
99-100% of the time; their own-half mode jumps more than 20 deg for 10%
(stage 3), 16-18% (stage 4), 20-21% (stage 5).

**Reading as written on 2026-10-22, withdrawn the same day (section 11).**
The reading was that the own-half mode follows the assignment through the
particle's own noise in its half's reference, while the other half's
reference, which never saw the particle, keeps the data's answer. The
premise is false for this run: the staged `refine3D` children of solve3D
always receive `lp`, `l_lpset` is then true, and `read_mask_filter_refvols`
copies one map into both reference slots (the average in stages 3-4, the
merged NU map in stage 5). The "cross-half" row scored the *same*
reference as the own row. What the two rows differ in is the search: the
own row is the production row (likelihood-sampled in-plane rotation for
every direction, continuous shift-and-rotation refinement of the best
`projs_ns`), the "cross" row is deterministic (best grid rotation at the
seed shift, no refinement).

**What the data of this section do show.** With the reference held fixed
between the two rows, the production row's mode moves more than 20 deg
between iterations for 10-21% of the particles both rows call decided,
while the deterministic row's mode moves for 0-1% of them; and the
production far-move rate is flat across the deterministic row's
confidence. The instability of the carried mode therefore comes from the
search procedure itself, from the in-plane sampling of the unrefined
entries and the refinement of the top set, before any assignment draw is
made. The "bonus" column, which with identical references compares a
continuously refined (shift, rotation) against the best grid rotation at
the same shift, should be small and non-negative; its median is, its tail
(8-18% above 20 nats) is not explained and needs the identical-pose null
test of section 11 before either path's number is trusted. Nothing in
this section measures self-reinforcement or reference independence.

## 10. Implementation: placement, boundaries, decisions (2026-10-22)

### 10.1 Where the draw lives

Two placements fit section 4.1. **On the master, in `ref_assign`** (chosen):
the master already owns the assignment (`prob_align_exec` reads the parts
into `dists(nptcls, nrefs)`, calls `ref_assign`, writes
`assignment_part.dat`); the per-particle draw is a second branch of that
routine over the same rows, and the committed candidate comes through the
same `dense_candidate` / `assign_candidate` path, which pulls the refined
shift and in-plane of the chosen reference from the `refined` list. The
master lacks one number per particle, the temperature T_i (it has no
`wsqsum`), so the worker writes T_i beside the seed shifts in the dense
part file header (one real per particle; `write_seed_shift_table` /
`read_seed_shift_table` grow by one array, the format word of the dense
part file is bumped). The random stream stays on the master, as the
frontier draw's is, so shared-memory and distributed runs draw alike. The
alternative, drawing on the worker at table fill and passing the choice
through the part file, moves assignment ownership into the worker and
duplicates the candidate path; rejected.

### 10.2 Boundaries

| Layer | Change | Owner file(s) |
| --- | --- | --- |
| params | `prob_assign='balanced'` (`balanced|posterior`), `l_prob_assign_posterior`; validation; `posterior` needs `objfun=euclid` | `simple_parameters*.f90` |
| UI | the key on `refine3D` and `solve3D` (advanced); solve3D forwards its command line to refine3D and the table workers unchanged | `simple_ui_refine3D.f90`, `simple_ui_solve3D.f90` |
| table utils | `draw_posterior_row(n, dists, temp, ...)`: the weights of `posterior_row_stats`, one draw by the cdf (the `draw_likelihood` pattern, `ran3`); pure apart from the draw, unit-tested | `simple_eul_prob_tab_utils.f90` |
| table | the dense part file carries T_i per particle; `ref_assign` gains the posterior branch: per particle, draw a reference from its row (within the argmin state for `nstates > 1`), `assign_candidate(i, dense_candidate(ri, i))`, `record_seed_shift`, `frac = 100` as today; the frontier code is untouched | `simple_eul_prob_tab.f90` |
| pre-alignment | none: `prob_align_exec` calls `ref_assign` as today | `simple_prob_align_strategy.f90` |
| search, matcher | none: `strategy3D_prob` commits the assignment as today | |
| reconstruction, assembly, sigma2 | none | |
| `prob_neigh` tables | none in this step (section 4.5) | `simple_eul_prob_tab_neigh.f90` |
| controller | none in code; the convergence question of 10.3 may need a stage-control change | `simple_solve3D_controller.f90` |
| tests | the draw (T to infinity gives the argmin, T = 0 uniform, seeded reproducibility); the dense part round trip with the temperature array; `ref_assign` posterior branch on the synthetic tables of `simple_eul_prob_tab_tester` (every particle gets a reference; with a huge T every particle gets its argmin); the fast gate | `simple_eul_prob_tab_utils`/`_tab` testers |
| docs | `refine3D_policy.md` (the assignment section), `solve3D_policy.md` section 3 if the stage controls change, this note | |

The sample-once-and-reproduce contract between `prob_align`, `prob_tab`
and `refine3D_exec`, the particle read contract, the artifact handoffs
(`assignment_part.dat`, the part files) and the shared-memory/distributed
parity all stand: the only format change is one array in the dense part
header, read and written by the same two routines.

### 10.3 Decisions to take before coding

1. **Convergence under a cold draw.** `check_conv3D` ends a stage on
   `frac_srch > fracsrch .and. mi_proj > overlap`, with `overlap` 0.9 in
   stages 4-6 and 0.95 in 7-8 (`solve3D_controller`). `mi_proj` is the
   fraction of particles keeping their direction between iterations; under
   the frontier draw it sits at 0.2-0.5 through the climb and the stages
   run to their budgets or to `lpexit`. Under the posterior draw it will
   sit near the posterior's own stability, 0.9+ from early in a stage, and
   the overlap rule would end stages at once. `solve3D_policy.md` already
   records that the early stop cannot tell a frozen assignment from
   convergence. Raising `overlap` does not help: `mi_proj > 0.99` is still
   met by a frozen assignment at `mi_proj = 1` (review finding 5). Under
   `prob_assign=posterior` the pose-overlap exit must be bypassed
   explicitly (a `check_conv3D` branch that converges on `lpexit`, the
   res05 gain guard, or `MAXITS` only), with a unit test of frozen
   assignments and a still-improving map. A map-based rule (res05
   stalled) for every mode is the larger change and the better one.
2. **Multi-state.** Keep the state label deterministic (argmin over all
   states, as `assign_state_labels` does) and draw the direction within
   the state; the state-level sampling of `refine3D_states` is a separate
   contract.
3. **Tempering.** Start at the exact posterior temperature (section 4.4's
   factor = 1); add the schedule only if stage 3's first iterations
   collapse.

### 10.4 Not decided by this note

Whether the balancing is needed for the stage-3 consolidation (section
4.3): the regression decides. The `prob_neigh` stages: a later step, same
rule over the neighbourhood row. The weighted form: only if the draw noise
shows.

## 11. Review of 2026-10-22 and what it changes

`pose_marginalization_review.md` found three P1 and two P2 issues; all
five hold against the source, and the first withdraws section 9's
inference.

1. **The "cross-half" row scored the same reference (P1).** `l_lpset =
   cline%defined('lp')` (`simple_parameters_phases.f90`); the controller
   sets `lp` on every stage (`emit_refine3D_stage_cfg`); in
   `read_mask_filter_refvols` the LP-set branch copies the average into
   the odd slot, and `l_use_merged_nu_ref` copies the merged NU map into
   both. `gen_best_euclid_val_cross` flips `ieo` correctly (the unit test
   proves it on distinct slots) but in solve3D's stages the slots are
   identical; the policy says so itself (no gold-standard stage). Section
   9 is corrected above.
2. **The two rows used different searches (P1).** Own: sampled in-plane,
   selective continuous refinement; "cross": best grid rotation at the
   seed shift. Sharpness, mode moves and the bonus confound reference and
   search. Section 4.5's "in-plane and shift at their maxima" was wrong:
   the row holds a mixture of sampled and selectively refined nuisance
   coordinates, and exponentiating it is an approximate direction sampler,
   not a joint pose posterior.
3. **The calibration ratio checks the mean residual scale only (P1).**
   `temp * lmin` over the component count near 1 does not validate the
   likelihood ratios between poses: correlated polar samples, the mask and
   model error can keep the mean and still make the posterior
   overconfident. The temperature's algebra is pinned by tests; its
   coverage is not.
4. **Monte Carlo equivalence (P2).** Holds for the accumulators
   conditional on fixed weights, not for the restored map; and the PfCRT
   benchmark ran under forced full sampling, where the trailing chain is
   off. Section 4 corrected.
5. **`overlap=0.99` does not bypass a frozen exit (P2).** Section 10.3
   corrected.

**Diagnostic v3 (implemented 2026-10-22, not yet run).** The same key,
re-cut so each column answers one question (the columns are in section 5):

- *Identical-slot detection.* Done: the two slots of the first reference
  compared at the table fill, `# slots_identical yes|no` in the header;
  when identical the other-slot row is skipped and its columns copy the
  own-slot row, and the summarizer says the reference test is unavailable.
- *Identical-pose null test.* Done (`null_a`, `null_b`). The search's
  stored value comes from the same FFT kernel as the table
  (`gen_objfun_vals` -> `gen_raw_euclid_vals_impl`), but the stored shift
  is rotated to the frame of reference first (`grad_shsrch_minimize`,
  `matmul(xy, rotmat)`); `null_a` undoes that with `shift_ref2ptcl` and is
  expected to be zero, `null_b` does not and is the convention the v2
  "bonus" column used, which is the likely origin of its tail. A nonzero
  `null_a` would mean the search's objective and the table's scoring
  disagree, a production finding.
- *Three rows, one search.* Done: the production row, the deterministic
  row on the own slot (the search test, `dist_d`, `rank_d`), the
  deterministic row on the other slot (the reference test, `dist_x`,
  `rank_x`). Section 9's data were the search test.
- *A real cross-half in solve3D.* Done (2026-10-22, the first option): the
  master writes a second polar model from the independent even/odd maps
  (`reprojection_model_diag_{even,odd}.bin`), the dense table workers load
  it after the matcher's references, and the own-half, other-half and
  self-reinforcement columns score it (section 5). The matcher's references
  are untouched; the search never sees the extra references.
- *The self-reinforcement bonus, done right.* `self_nats` is T (L_other -
  L_own) at the pose the particle was inserted at, on the real halves, with
  `ctrl_nats` a half turn away in-plane as the noise scale. The v2 "bonus"
  (same reference, rotated-frame shift) is superseded.
- *Known-pose simulation.* Particles projected from a known volume with
  known noise, scored through the production path: the fraction of
  particles whose true direction lies within the 99% set (coverage), the
  stability of n_eff under a finer grid, and a Cartesian evaluation of the
  same likelihood for a handful of particles against the polar one. This
  is the only validation of the temperature as a likelihood scale.

**Review gate for section 4.** (1) reference independence established:
the `_h`/`_x` columns score the real halves (built; to be read off the
re-run); (2) own and cross scoring comparable (both deterministic rows use
the same procedure; built) and the null test passed (`null_a` at zero on
the re-run); (3) coverage validated on known poses; (4) the convergence
bypass in place and unit-tested, full-sampling and fractional modes both
covered; (5) best-pose, single-draw and weighted updates compared on
identical candidate sets and references. The design of section 4 is not
built before (1)-(4); (5) is the experiment itself.

**The one run that settles (1) and (2).** `simple_exec prg=solve3D ...
nstages=3 prob_post_diag=yes` on PfCRT (the stages whose references and
posteriors sections 8 and 9 describe): `null_a` at zero and `null_b` not
settles the null test; `rank_d`/`dist_d` the search test on the
deterministic row; `rank_x`/`dist_x` the reference test on real halves;
`self_nats` against `ctrl_nats` the self-reinforcement test, read against
the own-half posterior's sharpness (`neff_h`, `pmax_h`): a mean `self_nats`
of the order of a nat or more where `neff_h` is near 1 would mean the
inserted pose is favoured by its own insertion by as much as the data
separate neighbouring directions.

