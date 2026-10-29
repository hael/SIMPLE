# Pose marginalization in the probabilistic 3D search

Implementation note, 2026-10-21, revised 2026-10-22 twice: after the
measurements (sections 8-9) and after the review
(`pose_marginalization_review.md`, section 11), and kept as the execution
record since (sections 16-20). **Status, 2026-10-26:** built and in
production on every probability table (sections 16 and 19): one draw per
particle from its tempered row posterior, no balanced frontier, no argmin,
`prob_beta` 0.3 flat with an annealing tail to the physical posterior on the
last stage (19.8); streptavidin 10/10 at the frontier's resolution (19.7),
PfCRT's ten restarts still to be made on the build that regularizes the
gridding map beyond the cutoff (19.9), with the 3.6 A bar set aside (21.1).
Section 20 is the next design step: one exploration level and a particle
count that follow the resolution ladder, in place of the flat `prob_beta`
and the flat `nsample`; section 21 the calibration of the temperature on
real data, which is the review gate's last open item. The sentence below is
the note's original frame; what it asked for is what sections 16-19 built.
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
| `null_a_nats` | the identical-pose null test: the production mode's stored loss re-computed at its exact (reference, shift, grid rotation) with the kernel the stored value came from, minus the stored loss, in nats; expected zero up to precision. `mode_src` (last column) says which case: **1** the mode was shift-refined, the stored shift is the search's particle-frame value (`shift_ref2ptcl` of the kernel-frame minimizer) and `null_a` feeds `shift_ptcl2ref` of it to the shifted-reference kernel (`gen_euclid_val_at`); **2** the mode is an unrefined entry scored through the seed bank, `null_a` is the bank kernel at the stored rotation (`gen_euclid_val_at_bank`); **3** unshifted, the kernel at zero shift |
| `null_b_nats` | the contrast: case 1, the stored shift fed to the kernel as it is (the wrong frame; nonzero whenever the shift is); case 2, the shifted-reference kernel at the frame-converted seed, which differs from the bank kernel by the angular Nyquist bin the former drops (the size of that bin on real data); case 3, equal to `null_a` |
| `iproj_d`, `lmin_d`, `pmax_d`, `neff_d`, `k99_d`, `spread_d_deg` | the deterministic row on the matcher's references: best grid rotation for every direction at the particle-frame seed the production search scored at, through the seed bank when the particle has one (`gen_best_euclid_val_bank`; a reference shift is one frame per rotation, so the seed cannot be fed to the shifted-reference kernel for all rotations at once), unshifted otherwise (`gen_best_objfun_val`); no sampling, no refinement; and its posterior |
| `dist_d_deg`, `rank_d` | the search test: the angle between the production mode and the deterministic mode, and the production mode's rank in the deterministic row (1: the same mode) |
| `iproj_h` ... `spread_h_deg` | the same deterministic row on the particle's **own real half** (the diagnostic model below, the reference at `iref + nrefs`), and its posterior |
| `dist_h_deg`, `rank_h` | the angle between the matcher-row deterministic mode and the own-half mode, and the matcher mode's rank in the own-half row: what the merged reference and the own half disagree on |
| `iproj_x` ... `spread_x_deg` | the same deterministic row on the **other real half** (`cross=.true.`, the bank or unshifted kernel at `iref + nrefs`), and its posterior |
| `dist_x_deg`, `rank_x` | the reference test: the angle between the own-half and other-half deterministic modes, and the own-half mode's rank in the other-half row (1: the two halves agree on the mode) |
| `iproj_p` | the pose the particle went into its half with: its previous direction on the current grid (`find_closest_proj`; exact within a stage), 0 when there is none (a virgin field, an invalid state) |
| `self_nats` | the self-reinforcement test: T (L_other - L_own) at that pose (previous direction, previous in-plane rotation, zero residual shift: the image is pre-shifted by its stored shift). Positive when the own half, which holds the particle, favours the pose it was inserted at. The pose was chosen on the previous merged reference, to which both halves contributed alike, so the selection bias is the same for both halves and the mean over particles isolates the particle's own contribution. |
| `ctrl_nats` | the same difference at the same direction a half turn away in-plane, where the particle's own contribution does not correlate with it (a half turn is the Friedel mate on the polar grid: the cross term is Re X^2, zero in expectation): the noise scale of `self_nats`, expected zero on average |
| `mode_src` | the null test's case, 1-3 as above (0: no valid mode) |

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
plain read is the `iref_from` 1 case; `test_bank_and_at_pose_kernels`
(same suite) that the at-pose loss at a nonzero reference shift is the
shift search's objective (`gen_objfun_vals`) at that pose, that a seed-bank
entry scores at its best rotation as the reference shifted by
`shift_ptcl2ref(u, rot)` up to the Nyquist bin (and not as the unconverted
seed), that the bank's best and at-pose losses agree, and that `cross`
swaps the half for the bank kernels too. The frame conversions themselves
are pinned by `test_shift_frames` (suite "in-plane search").
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
- *Two frame errors in the v3 diagnostic, found by the first full run
  (section 12) and fixed 2026-10-22 evening.* (i) The null test converted
  the stored shift with `shift_ref2ptcl`, which is the conversion the search
  applies when it *stores* the shift (`matmul(x, rotmat)` is `shift_ref2ptcl`
  of the kernel-frame minimizer, `test_shift_frames`); undoing it needs
  `shift_ptcl2ref`. Both v3 columns were therefore wrong-frame, which is why
  `null_a` and `null_b` came out alike and nonzero (section 12). (ii) The
  deterministic rows fed the particle-frame seed to the shifted-reference
  kernel for every rotation; a reference shift is one frame per rotation,
  so the production scores a seeded particle through the seed bank
  (`memoize_shifted_ptcl`, `gen_likelihood_val_bank`), and the rows now do
  the same (`gen_best_euclid_val_bank`, with `cross`). The rows of the first
  run are therefore at a shift that is wrong by up to twice the seed, the
  same for all three rows, so their comparisons with each other (the
  reference test) stand and their comparison with the production row (the
  search test) is indicative only. The null test now also records which
  kernel the stored value came from (`mode_src`), because the production
  row mixes two: the seed-bank kernel sums the angular Nyquist bin, the
  shifted-reference kernels (the unseeded scores and the shift search's
  objective) drop it, and the refined entries of a row are on the latter
  while its unrefined seeded entries are on the former. `null_b` of the
  bank case measures that gap on real data.
- *Known-pose simulation.* Particles projected from a known volume with
  known noise, scored through the production path: the fraction of
  particles whose true direction lies within the 99% set (coverage), the
  stability of n_eff under a finer grid, and a Cartesian evaluation of the
  same likelihood for a handful of particles against the polar one. This
  is the only validation of the temperature as a likelihood scale.

**Review gate for section 4.** (1) reference independence established:
the `_h`/`_x` columns score the real halves (done, sections 12 and 14);
(2) own and cross scoring comparable (both deterministic rows use the same
procedure) and the null test passed (section 14: `null_a` zero to three
decimals for every mode); (3) coverage validated on known poses; (4) the convergence
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
separate neighbouring directions. The first such run (a full one) is
section 12; it settled the reference and self tests and exposed the frame
errors above, so the null and search tests wait for the re-run with the
40-column diagnostic.

## 12. Results of the real-halves run: PfCRT, 2026-10-22 (5_solve3D)

A full solve3D (8 stages, 92 iterations, final 4.0 A at FSC 0.143 /
6.9 A at 0.5) with `prob_post_diag=yes`, 39-column files for iterations
41-79, the stages that fill the dense table (stage 3: iterations 41-57,
box 108, lp 18.3 A; stage 4: 58-72, box 120, 12.9 A; stage 5: 73-79, box
130, 11.1 A; stages 1-2 and 8 do not run the dense fill). `# halves yes`
throughout: the diagnostic model of the real halves loaded on every worker.
Analysis in `PfCRT_posterior/analysis/` (`post5.py` -> `run5_summary.txt`,
`self_lock.py` -> `run5_self_lock.txt`). The null and search tests of this
run are invalid (the frame errors of section 11); the reference and self
tests are at a common, slightly wrong shift for both halves and at zero
shift respectively, and stand.

**Self-reinforcement is real and grows with the band.** `self_nats`, the
own half's advantage at the pose the particle was inserted at, averages
+0.81 nats (stage 3), +2.30 (stage 4) and +4.71 (stage 5), rising
monotonically within stage 5 from +1.6 at its first iteration to +6.8 at
its last; the control a half turn away averages 0.00 +- 0.01 to 0.03 at
every stage (|mean| <= 0.13 in every iteration). The per-particle spread is
large (f(ctrl > 1 nat) 77-87%), so the bonus is a mean effect of a few
nats on a noise of several, not a per-particle certainty.

**It tracks the admitted band, not the stage's low-pass.** Within stages 3
and 4 (uniform low-pass at the stage `lp`) the bonus is flat across
iterations (0.77-0.86 nats; 1.8 then 2.3-2.4). Stage 5 is the first NU
stage (`NU_FILTER_STAGE = 5`): the nominal `lp` is 11.1 A throughout, but
the NU-filtered references admit the shells where the halves' FSC says
there is signal, and the half-map FSC=0.143 resolution runs 9.7, 8.2, 7.8,
7.1, 6.6, 5.9, 5.2 A over iterations 73-79 while the bonus runs 1.6, 3.3,
4.3, 5.1, 5.7, 6.3, 6.8 nats, monotone in the admitted band. The particle's
own noise is a larger share of its half's section at the shells where the
half is noise-dominated, so the band lead the NU filter buys is paid for in
self-reinforcement, at about one nat per Angstrom of resolution gained
here. On the merged reference the matcher scores against, the share is a
half: 0.8 to 3.4 nats across stage 5.

**It is only weakly a share-of-the-direction effect.** By the population
of the inserted direction at that iteration (both halves pooled; stage 3
median 25 of 1000 directions, stage 5 median 10 of 2500), the bonus falls
from 1.04 (lowest quartile) to 0.67 nats (highest) in stage 3, 2.5 to 2.1
in stage 4, 5.5 to 3.3 in stage 5; `self * pop` grows four- to threefold
across the quartiles where a 1/N law would keep it flat. Most of the bonus
is therefore not the particle's share of one section but its contribution
spread by the gridding kernel and the filter over the neighbourhood, and
the sparsely populated views, the ones a small membrane protein with a
preferred orientation is short of, carry the largest bonus. In loss units
the bonus is nearly constant across the particle's signal (self/T 1.3-1.5e-3
in stage 3, 2.1-2.2e-3 in stage 4, 0.9-1.2e-3 in stage 5), so in nats it
scales with the particle's weighted power: the strongest particles, whose
posterior gaps are also the largest, carry 2.2 (lowest quartile) to 6.8
nats (highest) in stage 5. The sharper the own-half posterior, the larger
the bonus (5.8 nats at n_eff_h <= 1.05, 1.2 at n_eff_h > 3 in stage 5),
which is cause as much as consequence: the bonus is in that row.

**It is of the size the independent half separates the mode by.** The
other half's posterior gap, ln(p_max_x / (1 - p_max_x)), has median 3.9
(stage 3), 6.4 (stage 4), 3.6 nats (stage 5). `self_nats` exceeds that gap
for 21%, 26% and 51% of particles, against 18%, 16% and 28% for the
control: beyond noise, the own-half bonus alone outweighs the other half's
separation of the mode for 3%, 10% and 23% of particles. On the merged
reference the matcher uses the particle's share is half, so 1-3 nats in
stage 5 against gaps of a few nats.

**The own half keeps the inserted direction more often than the other
half.** The deterministic own-half mode is the inserted direction for 78%,
84%, 71% of particles (stages 3-5), the other half's for 76%, 79%, 62%, the
matcher's merged reference for 78%, 84%, 69%. Where the two halves
disagree on the mode (12%, 14%, 24% of particles) the own half keeps the
inserted direction 44%, 64%, 59% of the time, the other half 30%, 30%, 23%.
The halves' modes coincide (rank_x = 1) for 88%, 86%, 76% of particles,
median dist_x 0 deg, more than 20 deg apart for 3% at every stage: the two
halves agree far more than the production search agrees with either (next
paragraph). Posteriors on either real half are sharp: median n_eff_h
1.10/1.01/1.07, n_eff_x 1.11/1.01/1.14, f(n_eff_h > 2) 14%/4%/11%.

**The production row's mode is the inserted direction far less often than
any deterministic row's: 46%, 39%, 21%.** The production row (likelihood-
sampled in-plane for every direction, the top `projs_ns` refined) agrees
with the deterministic row on the same references for 61%, 54%, 46% of
particles (`rank_d` = 1), while two deterministic rows on *different data*
(the two halves) agree for 88%, 86%, 76%. With the caveat that the
deterministic rows of this run are at a slightly wrong shift, the ordering
is not in doubt: the sampled in-plane rotation stores, for most
directions, a loss above the direction's best, so the row's minimum and
its exponential are dominated by which directions happened to draw their
best rotation. This is section 9's search effect measured against
independent data, and it says the production row as filled is not the data
posterior over directions; a posterior for section 4 has to be built from
deterministic best-in-plane values (the refined set, or a row scored at the
best rotation), not by exponentiating the table as it stands.

**Null test (invalid in this run, kept for the record).** `null_a` median
+0.5, +1.0, +2.0 nats with p90 22-32 and `null_b` the same to within 5%:
both columns used the wrong frame (section 11), and the re-run with
`mode_src` decides it. A median of a nat with a tail to tens is the size a
frame error on sub-pixel refined shifts gives at these bands (a 0.2 px
error is a third of a radian at the band edge), so the result is consistent
with the diagnostic's error and says nothing yet about the search's
bookkeeping.

**Consequences for the design.** (1) Self-reinforcement at 1-7 nats on the
own half, 0.5-3.5 on the merged reference, against mode gaps of a few nats,
is a locking force from the NU stage on; marginalization over the row
posterior does not remove it (the bonus is in the row). Scoring against
references the particle did not enter removes only the direct term and
couples the halves (section 13, parked); a correction for the particle's
own contribution, or gold-standard matching once the NU machinery is on
(Hans: not in the ab initio stages, and only with extensive validation),
are the open options. (2) The posterior must
come from deterministic best-in-plane values. (3) The halves' disagreement
on the mode (12-24%, mostly neighbouring directions) is the honest
uncertainty the data leave at these bands; it is what a draw from the
posterior should reproduce, and it is small compared with the production
search's own scatter, which should be removed first.

Analysis scripts: `analysis/self_lock.py` (`run5_self_lock.txt`),
`analysis/self_mechanism.py` (`run5_self_mechanism.txt`).

## 13. Cross-half matching: an experimental search policy, parked (2026-10-23)

The first draft of this section (2026-10-23 morning) called scoring each
particle against the half it did not enter "gold-standard" and claimed
"zero self-reinforcement at every shell". Both claims were wrong, as the
review of the same day (`pose_marginalization_review`, finding 1) says.
Gold-standard refinement is the *opposite* arrangement: each half-set is
refined against its own map, so the particle's own insertion is fully in
its reference, and what the scheme buys is the independence of the two
refinement histories, hence an honest FSC. Cross-half scoring removes the
particle's *direct* term from its reference but couples the two histories:
even poses depend on the odd map, which depends on odd poses, which depend
on the even map the particle is in. The indirect self-influence remains,
one loop later, and the FSC between the halves is no longer that of two
independent refinements. The three arrangements are distinct: the merged
reference of solve3D's LP-set stages (half the direct term, coupled
histories), gold standard (the full direct term, independent histories),
cross-half (no direct term, coupled histories).

Hans's constraint, 2026-10-23: full gold-standard matching cannot be used
in solve3D's ab initio stages, where it does not converge on a good model
from the start it has; re-introducing it when the NU machinery is on is
conceivable but would need extensive validation of its own. Cross-half
scoring, if ever tried, is an experimental search policy with its own
validation, not a remedy the design relies on. Nothing in this note's
plan depends on it; section 16 proceeds without it.

What remains useful from the draft: the mechanics are a slot assignment
(the kernels take the reference slot from the particle's half, so loading
the halves swapped makes every score cross-half with no kernel change;
the GPU `gen_many2many_*` kernels read `REF_EVEN` alone and would need
the same rule), and the diagnostic's `_x` columns already measure what
such a reference sees. The direct self term is 0.8-6.8 nats on the own
half at stages 3-5 (sections 12, 14), half of that on the merged
reference; what the coupled loop adds back is not measured.

## 14. Results of the frame-correct re-run: PfCRT, 2026-10-22 (4_solve3D, stages 3-5)

The same schedule with the 40-column diagnostic (section 11, the two
frame errors fixed), read at iteration 68 (stages 3 and part of 4) and
again at iteration 87 (stages 3-5 complete: 41-57, 58-72, 73-80; stage 8
from 81 at box 150, no dense fill). Analysis `post5.py` ->
`run4b_summary.txt`, `search_test.py` -> `run4b_search_test.txt`,
`self_lock.py` -> `run4b_self_lock.txt`. Stage 5 confirms the three
findings at the NU band: `null_a` zero to three decimals for the refined
modes (77%; `null_b` median 2.9 nats, p90 44) and for the bank modes;
production mode equal to the deterministic mode for 49%, worse than the
deterministic best by more than a nat for 27% of all particles, equal to
the other half's mode for 50% against the deterministic rows' 88%;
`self_nats` +4.0 nats on average, rising 1.5 -> 5.4 over iterations
73-80 as the half-map FSC=0.143 resolution runs 10.4 -> 7.6 A (the
schedule of this run promoted later than section 12's), control -0.05.

**The null test passes.** For the shift-refined modes (83% and 80% of the
particles), the shifted-reference kernel at the frame-converted stored
shift reproduces the stored loss to three decimals in nats
(`null_a` median and p90 +0.000, f(|null_a| > 0.1 nat) = 0.000 over
239k and 143k particles); the unconverted shift (`null_b`) is off by a
median 0.8 and 2.0 nats with p90 26 and 40, which is what the first run
had reported as `null_a` and `null_b` alike. For the modes the refinement
left at their seed value (17-20%; `mode_src` 2, the top set whose shift
search found nothing better) the bank kernel reproduces the stored loss
exactly too, and the shifted-reference kernel at the frame-converted seed
agrees with the bank kernel to better than 0.001 nats: on real data the
angular Nyquist bin the two kernels differ by is below 1e-6 of the loss,
so the row's mixture of kernels (section 11) is harmless. The production
search's bookkeeping is clean: every stored value is the kernel's loss at
the stored pose.

**The search test stands, and it is the in-plane sampling.** With the
deterministic rows now at the production's own seed (through the bank),
the production row's mode equals the deterministic mode for 63% (stage 3)
and 54% (stage 4) of particles, as in the first run, and the agreement
falls along each stage as the reference sharpens (0.69 -> 0.60 over stage
3, 0.57 -> 0.50 over stage 4; modes more than 20 deg apart 12% -> 19%).
Where the two modes differ the production's refined loss is *higher* than
the deterministic best in 56% and 63% of cases (median by 0.6 and 1.5
nats, p10 by 6 and 9); the production committed to a pose more than a nat
worse than the best grid pose on the same reference for 17% (stage 3) and
25% (stage 4) of all particles, a quarter of them more than 20 deg away.
Against the independent half's mode the deterministic row on the matcher's
references agrees 94% and 93%, the production row 62% and 54%: the
production search scatters the pose far more than the data's own
ambiguity (the two halves disagree for 12-13%).

The mechanism is in `score_ref`: every direction is scored at an in-plane
rotation *drawn* from the `nsample` lowest losses
(`inpl_likelihood_nsample`, `prob_athres` 10 deg capped by the mean
`dist_inpl`, so 5-20 rotations at these boxes), and the draw
(`draw_likelihood`) weights by exp(-(L - L_min)) on the *normalized* loss,
where neighbouring rotations differ by 1e-3: the draw is uniform over the
candidates, the stored loss is the loss at a random rotation within the
threshold, and the direction whose draw happened to land near its best
wins the row. The refinement of the top `projs_ns` cannot repair this,
because the unlucky direction is not in the top set. This is the
assignment's T = 1 draw (section 2) one level down, applied to every
entry of every row.

**Self and reference tests reproduced.** `self_nats` +0.86 (stage 3) and
+2.05 (stage 4) with the control at -0.02 and -0.07; own-half vs
other-half mode agreement 88% and 87%; the own half keeps the inserted
direction for 79% and 84%, the other for 77% and 79%.

**What follows.** Before any posterior is drawn from the table, the table
has to hold the best grid rotation per direction, or a rotation drawn at
the particle's temperature (`posterior_temperature`, available on the
worker), so that the draw is a posterior sample and not a uniform one.
The cheap experiment was built as a switch (`prob_inpl_draw=argmin`:
every table entry at its best grid rotation, in the 3D dense fill, the 3D
neighbourhood fills including the stochastic styles' extremal count and
the assignment's in-plane fallback, and the 2D table, where
`neighfrac2nsmpl` never goes below 2 and the entries otherwise keep a
one-step dither to the end of the extremal decay), run ten times, and
removed again the same day: it traps the search (section 15). The noise is
the exploration. What remains of it in the code is the comment on
`inpl_likelihood_nsample`. The principled version tempers the direction
draw with T on an explicit schedule (section 15); it is the first item of
section 4's design and no longer needs the search to be trusted first,
since the null test holds.

## 15. The argmin runs: PfCRT, 2026-10-23

The schedule of sections 12 and 14 with `prob_inpl_draw=argmin
prob_post_diag=yes`. The baselines are `prev/5_solve3D` (section 12,
finished, 4.03 / 6.90 A) and `prev2/4_solve3D` (section 14, stopped at
iteration 87 in stage 8, 7.22 / 8.62 A at box 150). Analysis
`run6_summary.txt`, `run6_search_test.txt`, `run6_self_lock.txt`.

**The search noise is gone.** Production mode equal to the deterministic
mode: 89%, 90%, 85% (stages 3-5; before 63%, 54%, 49%). Where they differ
(11%, 10%, 15%) the production's refined loss is lower in 99.9-100% of
cases, by a median 1.1-1.6 nats: the remaining differences are the
continuous refinement finding a better pose at a neighbouring direction,
which is what it is for. No particle commits to a pose more than a nat
worse than the best grid pose (before 17%, 25%, 27%). Production mode
equal to the independent half's mode: 88%, 88%, 81% (before 62%, 54%,
50%; the deterministic rows 95%, 94%, 89%). Production mode equal to the
inserted direction: 82%, 87%, 70% (before 46%, 39%, 26%), 93% at the last
stage-5 iteration read.

**The halves agree more, the bonus is unchanged.** Own-half vs
other-half deterministic modes disagree for 9%, 9%, 17% (before 12%, 14%,
20%): less search noise goes into the two reconstructions. `self_nats`
+0.81, +2.19, +3.5 (rising 1.5 -> 4.9 over the first five stage-5
iterations), control -0.03 to -0.08: the switch did not touch the
reference, as it should not have. The null test passes as before; the
refined modes are now 18%, 13%, 31% of all (before 77-83%), because the
seed shift is already the optimum at the previous pose and the refinement
of a mode that stays there has nothing to improve.

**Convergence and resolution (interim).** Stages 1 and 4 converged in 12
and 9 iterations (before 20 and 15 in both baselines); stages 2 and 3 as
before (20, 17). The stage-end resolutions are the ladder's and identical
across the three runs. Within stage 5, the FSC=0.143 resolution over the
first five iterations ran 10.35, 9.41, 8.39, 7.96, 7.57 A, against 9.70,
8.17, 7.76, 7.06, 6.61 (`prev/5_solve3D`) and 10.35, 9.41, 8.87, 8.62,
8.39 (`prev2/4_solve3D`): the two baselines differ from each other by
more than either differs from the argmin run, so the stage-5 trajectory
is not the arbiter; the final map is.

**Final maps: argmin is worse.** Ten restarts with `prob_inpl_draw=argmin`
(`PfCRT_posterior/greedy_inpl`, 75-82 iterations each, stages 1-5 and 8):
FSC 0.5 / 0.143 of 7.76/4.63, 8.39/6.90, 7.06/4.03, 8.87/7.06, 8.39/6.75,
7.39/4.44, 8.62/7.06, 7.39/4.14, 7.76/4.57, 9.41/7.57 A. Five of ten end
near 7 A, five at 4.0-4.6 A, one at the default's 4.03 A; the default
reaches ~4 A reliably (Hans's runs; `prev/5_solve3D` 4.03). The sharper,
less noisy table of the first paragraph produces a worse and bimodal
outcome: the run converges faster (stages 1 and 4 in 12 and 9 iterations)
to a self-consistent solution that half the time is the wrong one.

**Reading.** The regression is evidence about the complete ablated
workflow: the switch took the in-plane draw out of the dense fill *and* the
neighbourhood fills (stages 1-2 and 6-8), and stages 1 and 4 ended earlier
(12 and 9 iterations), so the worse maps do not separate the loss of
exploration in the dense stages from the early neighbourhood search or
from shorter stages (review finding 4). A dense-stage-only ablation with
matched stage budgets would separate them; it was not run, and the
mechanism below is a reading of the numbers, not a measurement. In that
reading, the in-plane draw's randomness is not a defect to remove but
the search's exploration, and the two measurements of this note are the
two sides of one mechanism. Self-reinforcement (sections 12 and 14) pulls
every particle back to the pose it was inserted at, by 1-7 nats on its
own half; the in-plane noise kicks each entry by a comparable amount
(median 0.6-1.5 nats, p10 6-9), so the mode of a row keeps moving among
near-equivalent directions, the assignment stays diverse, and a wrong pose
can leave. With `argmin` the kick is gone: 93% of particles return to
their inserted direction each stage-5 iteration (section 15, first
paragraph; 46-82% with the draw), the reference confirms them, and the
run locks, right or wrong. `athres` shrinking with the measured in-plane
change is an annealing schedule: the noise is largest while the alignment
moves and dies as it settles. Hans's reading, that the effect is indirect
and acts through the diversity of the direction sampling while the final
in-plane angle is the continuous refinement's anyway, is what the numbers
say.

**Consequence for the design.** Section 4's draw from the row posterior
must carry the exploration the in-plane noise supplies today, or it may
trap as `argmin` did: a draw at the data temperature T over a posterior
with n_eff ~ 1 (sections 8, 12) on deterministic rows would be `argmin`
with extra steps. The design therefore keeps an explicit tempering
factor, section 4.4's: the weights are exp(-T_eff (L - L_min)) with
T_eff = beta T, beta < 1 *hotter* than the data (the first draft of this
paragraph wrote T / beta, which heats in the wrong direction; review
finding 3; the implemented `PROB_POSTERIOR_TEMP_FACTOR` multiplies T and
`test_posterior_row_stats` pins that a factor below one raises the
entropy), annealed towards 1 on the measured pose change if a schedule is
needed, and the noise placed where it belongs, on the draws, not on the
loss values the assignment compares. The schedule's values are to be
calibrated on the posterior runs themselves (section 16), not prescribed
from the ablation; section 14's numbers (the production mode's loss
excess over the deterministic best, a few nats at p50 and 6-9 at p10)
say what today's exploration amounts to and are the first value to try if
the data temperature traps. The `prob_inpl_draw` key was removed after
the ten runs (dead code is confusing code); the ablation is reproducible
by setting `nsample` to 1 in `inpl_likelihood_nsample` and `ninpl_smpl`/
`smpl_ninpl` to 1 in the neighbourhood and 2D fills.

## 16. Execution of sections 4 and 10 (2026-10-23)

The design is sections 4 and 10 as written; this section records what the
measurements of sections 12-15 add to it and the acceptance bar, and
nothing else. The gate of section 11 is open for items (1) and (2): the
reference columns score real halves and the null test passes (sections 12,
14). Item (4), the convergence bypass, is built with the key (10.3, 1).
Item (3), coverage on known poses, is not done and stays on the list.

- **Acceptance.** PfCRT only, until the new path reproduces and then beats
  the production path: 10 of 10 restarts, at least 3.6 A at FSC 0.143
  (the baseline is 4.0). Nothing is tested on other data before that.
- **Production untouched.** `prob_assign=balanced` is the default and
  bit-identical; `posterior` is the second branch of `ref_assign` (10.1),
  the frontier code stays.
- **Argmin is not a step.** Section 15 showed that taking the randomness
  out of the in-plane draw traps the run; the per-particle argmin would do
  the same one level up. The posterior draw of section 4.1 is built
  directly, at the temperature of section 3, with the tempering factor of
  4.4 as the one constant (10.3, 3: start at 1; T_eff = factor x T, a
  factor below one is hotter). Section 14's calibration says today's
  exploration corresponds to a factor of order 0.05-0.1, so if stage 3
  collapses at 1 the first value to try is known; the schedule (rising to 1
  along a stage on the measured pose change) is the follow-up the design
  already names.
- **Weighted form.** As in section 4: only if the per-iteration noise of
  the draw shows in stages 3-4. Not before.
- **2D.** Added 2026-10-23: the 2D table gets the same per-particle rule
  under the same key as a later step (a particle's class drawn from its
  row posterior at its temperature; the frontier default untouched), after
  the 3D path holds on PfCRT. `prob_neigh` as in 4.5, after that.
- **Diagnostic.** `prob_post_diag=yes` on every PfCRT run of the series, so
  the self, reference and search columns say what each step changed.

Order: (a) `draw_posterior_row` with its tests; (b) T_i in the dense part
file, round trip tested; (c) the `ref_assign` posterior branch with the
state rule of 10.3 (2), tested on the synthetic tables (huge T gives the
argmin); (d) the convergence bypass with its unit test; (e) the key, the UI
entries, the policy docs; (f) the PfCRT 10/10.

**Built 2026-10-23 (a)-(e), awaiting (f).** `draw_posterior_row`
(`simple_eul_prob_tab_utils`; tests: huge T a unique argmin every time and
either of two equal minima with equal probability (a flat row is a uniform
draw at any temperature, the posterior's answer rather than the balanced
rule's lowest-index convention), T = 0 uniform over the valid distances and
never an invalid one, T = 1 on [0, 1, 2] within 0.02 of [0.665, 0.245,
0.090] over 20000 draws, seeded reproducibility, an invalid row draws 0). The
dense part file carries `temps(nptcls)` after the seed table (format word
-4; `write_dense_head`/`read_dense_part`; the fill sets
`temps(i) = posterior_temperature(wsqsum, pftsz)` for every particle; the
round-trip test checks them). `ref_assign` branches to
`ref_assign_posterior` under `params%l_prob_assign_posterior`: the state
label as `assign_state_labels` (blocks of 256 particles, strict `<`), one
`draw_posterior_row` per particle within the state's block at
`PROB_POSTERIOR_TEMP_FACTOR * temps(i)` (the constant is 1.0; section 4.4's
knob), `assign_candidate(dense_candidate)`, `record_seed_shift`, `frac`
100; the log line reports the fraction of particles drawn at their row
mode; a particle without a temperature (none expected) takes its argmin
with a warning; the invalid/flat check is shared with the balanced branch
(`warn_invalid_and_flat`). Tests on the synthetic tables of
`compare_with_reference`: huge T gives every particle with a unique row
minimum its argmin, the flat particle a direction of state 1, the invalid
one reference 1, and the state ties their lower state; at T = 2 on
distances in [0, 5] some particles leave the mode and none leaves the
argmin's state; the same seed gives the same assignment, another seed
another. `check_conv3D` bypasses the pose-overlap
rule under the key for `refine=prob` (the dense table; the neighbourhood
stages keep the frontier and their rule in this step, section 4.5; a
Cartesian pass has its own rule) and logs it;
`test_posterior_assignment_bypass`: a frozen assignment converges the
default rule and not the key's on the dense table, a `prob_neigh` pass keeps
its rule, and the `lpexit` exit still ends the stage. Key `prob_assign=balanced|posterior` (`posterior` requires
`objfun=euclid`), on `refine3D` and `solve3D` as advanced search inputs;
`refine3D_policy.md` and `sampling_and_fractional_updates.md` updated.

**Review of 2026-10-23 (`pose_marginalization_review`, four findings),
all accepted.** (1) Cross-half is not gold standard: section 13
rewritten and parked. (2) The diagnostic's second reference block was
offset by `nrefs`, which counts the active states only, while the
calculator holds the full `nspace*nstates` production model and the
halves' model after it: wrong block, and a false "halves present" when
states are inactive. Fixed: `post_ref_offset = nspace*nstates`,
`post_has_halves = get_nrefs() == 2*post_ref_offset`. No effect on the
single-state PfCRT runs; the multi-state diagnostic has no unit coverage
(the fill needs a calculator and images), recorded here as a gap. (3)
The tempering convention: `T_eff = beta T`, beta below one hotter; the
text corrected in sections 15 and 16, the factor's direction pinned by a
test. (4) The argmin ablation does not isolate the dense search; section
15 now says so and the schedule is calibrated on the posterior runs, not
prescribed from the ablation. The review's standing limitation holds and
is section 16's own statement: the direction draw exponentiates a row with
sampled in-plane scores and selective refinement, an experimental
direction sampler until the in-plane draw is moved under the same rule
and the known-pose calibration (gate item 3) is done; final maps decide.

**First posterior runs, 2026-10-23 (f), interim.** PfCRT with
`prob_assign=posterior prob_post_diag=yes`, the schedule of section 12.
Restart 1: 94 iterations (stages 20/20/17/15/12/10), 4.03 / 6.21 A (FSC
0.143 / 0.5) — the same 0.143 as the baseline, the best 0.5 of the series
(6.90 and 7.06 before); stage 5 ran its twelve iterations without the
overlap exit and ended at 4.78 A against 5.17 (`prev/5_solve3D`) and 7.57
(`prev2/4_solve3D`). Restart 2: 104 iterations (stage 5 ended on the
resolution exit after 7, stage 8 ran 25), 3.98 / 7.06 A. Restart 3 in
stage 4 at the time of reading. Hans: the maps are the best-looking of
the series. The fraction of particles drawn at their row mode runs 0.63 at
a stage's first iteration to 0.90 at its last, the 0.8-0.9 the posterior's
width predicted (section 16's read-off), and no pose-overlap bypass ever
fired: the stages never met the overlap rule, so the bypass changed
nothing here. The diagnostic columns of restart 1 are the baselines' to
within noise (search test 64/53/44%, halves' agreement 87/86/78%,
`self_nats` +0.8/+2.3/+6.1, the last larger because stage 5 went deeper;
`run7_*.txt`), as they must be: the fill is unchanged, only the assignment
draws differently, and what the assignment changed shows in the maps, not
in the rows. The 3.6 A bar is not met by the FSC=0.143 of the first two
restarts; the ten-restart tally decides, and the in-plane draw (item 5)
is the next step either way.

**The in-plane draw and the 2D route, 2026-10-23 (built, awaiting the
runs).** Agreed with Hans after the first posterior runs: section 4.5's
"unchanged" in-plane is revised. The in-plane draw goes under the same key
by the same rule, and so does the 2D table; the balanced path stays
bit-identical. What the in-plane draw was (section 14, `inpl_likelihood_nsample`):
a uniform draw over the rotations within the angular threshold, the
search's exploration; removing it traps the run (section 15). What it is
under `prob_assign=posterior`: one draw from exp(-T_i (L_r - L_min)) over
every grid rotation of the reference against the seeded particle, T_i the
particle's temperature under the one tempering constant, no angular window
and no truncation to the best few (`gen_posterior_val`,
`gen_posterior_val_bank` in `simple_polarft_corr`: the Euclidean crvec
kernels, `draw_posterior_row` on the thread's rotation buffer, a row without
a valid rotation returns rotation 1 at the invalid marker as the likelihood
draw does; `objfun=euclid` only, THROW otherwise). The sites: the dense 3D
fill (`score_ref`, `score_ref_bank`), the neighbourhood tables
(`score_seeded_ref`, the three fills behind it, and the assignment's
fallback), the 2D tables (`score_seeded_class`, the shift-free direct call).
The refinement of the top set (alternating shift and rotation, greedy in
rotation) and the continuous in-plane polish are unchanged; so the
exploration the draw provides is now at the data's width, the particle
leaving a reference's remembered rotation as often as its posterior lets
it, and the direction draw exponentiates a row whose in-plane entries are
posterior draws, not uniform ones within a window (the review's standing
limitation, one step closer). Tests (`simple_polarft_corr_tester`,
`test_posterior_inpl_draw`): a huge temperature is the argmin of the same
kernel with its loss, in both forms; temperature 0 spreads over the grid
and returns the drawn rotation's at-pose loss; at the particle's own
temperature the draw's distribution is the posterior computed from the
row (`gen_objfun_vals` gives every rotation's loss): the mode's frequency
and the mass within two grid steps of it over 200 draws agree with the
computed masses to 0.1, in two regimes. The first version of that check
asked for the mode at least 90% of the time and failed, and the printed
temperature said why: at unit sigma2 the fixture's temperature is 0.005
(sigma2 = 1 is enormous against the image's power in the calculator's
units; the loss is normalized, the temperature is not), the posterior is
uniform (mode mass 1/nrots = 0.007, drawn 0.005) and the draw was right.
The test now checks that regime and a second one, sigma2 set at the
image's own power scale so that T = 100 (a constant sigma2 cancels in the
normalized loss and scales T as its inverse, so no FFT normalization has
to be known), where the posterior and the draws concentrate about the
mode. The lesson stands for the real runs: the draw's width is whatever
sigma2 says it is, and section 3's T of 500-4000 is sigma2 estimated from
the data; a sigma2 off by a constant factor moves every draw's
temperature by that factor.

The 2D route (section 16's 2D item, moved up with Hans's agreement): the
2D part files carry T_i per particle after the seed table (formats -6
stream, -7 dense), both fills write it, the master reads it
(`test_parallel_read` checks the placement). Under the key the assignment
is one draw per particle: over the active classes of its dense row
(`ref_assign_likelihood`, `assign_posterior`: the drawn class's
`dense_candidate`, `frac` 100, the first active class when no distance is
valid) or over its valid evaluated candidates
(`ref_assign_sparse_likelihood`, `assign_particles_posterior`: an active
class with an in-plane index, the drawn candidate with its shift, the peak
count and fraction as the frontier records them, the best candidate when
no distance is valid); a particle without a temperature takes its argmin
with a warning; the draws run in particle order on the master's stream; the
log line reports the fraction drawn at the row mode. Tests
(`simple_eul_prob_tab2D_tester`): on the synthetic tables of the reference
tests, a huge temperature gives every particle with a unique minimum its
argmin (distance and in-plane carried), a class tie one of the tied
classes, the invalid particle class 1 (dense) or a candidate at the
smallest valid distance with the frontier's peak count and fraction
(sparse), no temperature the argmin, T = 2 on distances in [0, 4] moves
some particles and none to a smaller distance than their best (sparse),
the same seed the same assignment and another seed another. The key is an
advanced search input of `refine2D`, `solve2D` and `solve2D_chunks`;
`solve2D_policy.md` has the record. The 2D convergence rules are untouched
in this step; whether the class overlap needs the bypass the 3D dense stage
has is for the 2D regression to say.

`PROB_POSTERIOR_TEMP_FACTOR` moved to `simple_eul_prob_tab_utils` so the
three tables share it (it became the run-time `prob_beta` after section
17); the in-plane and the reference/class draws are at
the same T_eff. Still pending from the plan after this step: the beta
schedule if stage 3 needs it; the known-pose calibration (gate item 3).
Next: PfCRT x10 with the in-plane draw on, and a 2D regression of Hans's
choosing for the class route.

**The strategy-level samplers and the `prob_neigh` assignment, 2026-10-23
(built, awaiting the runs).** The temperature every draw runs at is now
one function, `pftc%ptcl_posterior_temp(iptcl)` (the factor on
`posterior_temperature` of the memoized self-term), used by the three
tables and the strategies alike; `draw_posterior_scores` (utils) is
`draw_posterior_row` over a buffer of exp(-L) scores with the strategies'
sentinel (a score at or below zero is unscored), tested beside it.

*4.5, the `prob_neigh` assignment.* The neighbourhood part file carries
T_i after its seed table (both fills write it, `read_sparse_tab_to_glob`
reads it). `ref_assign_neigh` branches under the key after the seeding of
empty particles: the particle's state is the state of its closest valid
candidate over all states (an active reference with an in-plane index,
ties to the first in store order; the dense table's rule on the store),
and within that state one draw from exp(-T_i (L - L_min)) over its valid
candidates; the drawn candidate committed with its shift, the seed shift
recorded, `frac` as the frontier records it; a particle without a valid
candidate keeps its orientation (the frontier's leftover rule), one
without a temperature its argmin with a warning; the log line reports the
fraction at the mode. The pose-overlap bypass now covers `prob_neigh`
(`check_conv3D`; the tester asserts it, and that a `shc` pass keeps its
rule). Test (`test_ref_assign_neigh_posterior`): a synthetic store of 150
particles over 3 states x 20 directions, 1 to 12 candidates each, some
without an in-plane index; a huge temperature gives every particle its
closest valid candidate with its in-plane index, distance and evaluated
fraction; T = 2 moves some and none out of its state; seeds.

*The strategies.* Each keeps its structure and gets the posterior rule
where it drew: `snhc_smpl` and `snhc_smpl_many` draw each class's
rotation from the posterior at the seed shift (`gen_posterior_val`; the
batched scores through `draw_posterior_scores`) and the class from the
posterior over the whole evaluated row, the best `snhc_smpl_ncls` entries
shift-refined and the rest left as scored (`inpl_srch_peaks` with
`l_keep_unrefined`: the balanced path invalidates the unrefined entries so
its power sampler draws among the refined few, the posterior path draws
over the mixed row as the tables do); the drawn class's shift refined as
today. `inpl_smpl` draws the rotation from the posterior at the seed shift
and refines the shift at it (draw, then polish; one shift search instead
of `snhc_smpl_ninpl`). `greedy_smpl` draws each class's rotation from the
posterior, its class selection greedy as named. 3D `refine=neigh`
(`strategy3D_greedy_sub`) with `prob_inpl=yes` draws the neighbourhood's
rotations from the posterior in place of the bounded sampler within the
angular threshold; its coarse pass and peak selection stay greedy, as
named. The old samplers (`gen_likelihood_val`, `gen_prob_objfun_val`,
`gen_prob_power_objfun_val`, `power_sampling`, `squared_sampling`,
`angle_sampling`) remain the default's code; they go once every caller has
moved and the regressions have passed. No unit test covers a strategy
(they need a calculator and images); the regressions do.

The key is inert under `objfun=cc`: the parser sets
`l_prob_assign_posterior` only with `objfun=euclid`, and no program checks
it. There is nothing to draw at without the Euclidean loss (no
temperature, no sigma2), so the balanced path is the only one that can
run there, and an error served nobody: `solve3D_cavgs
prob_assign=posterior` died in its random-start reconstruction, which
inherits the workflow's command line with `objfun=cc` and makes no draw.
The two pftc primitives keep their contract guard.

## 17. The posterior everywhere: PfCRT, 2026-10-23 (4_solve3D, run 8)

The first run with the in-plane draw and the `prob_neigh` assignment under
the key, `prob_assign=posterior prob_post_diag=yes`, restart 1: 116
iterations (stages 20/20/17/15/7/-/12/25), 5.35 / 7.96 A (FSC 0.143 /
0.5), against 4.03 and 3.98 A for the two restarts of section 16's first
posterior runs (the key on the dense assignment only: likelihood in-plane
everywhere, the balanced frontier in the neighbourhood stages). The 2D
route is fine on Hans's regression; `solve3D_cavgs` on streptavidin fails
2 of 10, which is unusual. The trajectory and the comparison with
`prev/5_solve3D` (a 4.0 A run) are in `analysis/run8_trajectory.txt`.

**Where it is lost.** Stages 1-4 end where the earlier runs end (11.1 A
at the end of stage 4 against 10.7-11.1). Stage 5 is slower: 10.7 to 7.6 A
over its seven iterations where the 4.0 A run went 9.4 to 5.1 over the
same seven, and it ended on the resolution exit (FSC 0.5 at 9.70 twice,
the shell value; FSC 0.143 still gaining 0.2-0.4 A per iteration; the
exit's "no longer gaining" is quantised to a Fourier shell, and section
16's restart 2 met the same exit). The neighbourhood stages do not
recover: stage 7 8.4 to 6.8 A in 12 iterations, stage 8 6.8 to 5.3 in 25
and still gaining, where the balanced neighbourhood stages took the 4.0 A
run from 8.0 to 4.0 in eighteen.

**What the draws did.** The fraction drawn at the row mode is 0.93 at the
first iteration of stage 1, 0.80 at the first of stage 5 (0.63 before,
section 16) and 0.85-0.95 through stages 7-8; the pose overlap reaches
0.86 at the second iteration of stage 7 and 0.93 by the end; the
posterior's median mode mass is 0.97 by the tenth iteration of stage 3
and 0.99 in stages 4-5 (n_eff 1.02-1.05). The self test is as before:
+0.8 nats in stage 3, +2.3 in stage 4, +1.6 to +5.2 across stage 5's
seven iterations (p90 13 nats), growing with the NU band; the search and
reference tests are the baselines'.

**Reading (corrected after Hans's objection, 2026-10-23).** The first
reading blamed the self-term: the particle's own contribution to the
reference makes its previous pose look 2 to 5 nats better than the data
say, more in the deep band. The diagnostic does not support that as the
first-order cause. Against the other real half, where there is no
self-term, the row posterior is as narrow as the production row's: in the
run7 reference test n_eff_x is 1.20 against 1.18 for the production row by
the fourth iteration of stage 3, p_max 0.97. The self-term moves the mode
for a minority of particles; it does not make the posterior sharp. At the
physical temperature the posterior over directions is a delta with or
without the bias, the draw is the argmin, and the map-pose alternation is
a coordinate descent that stops at the nearest fixed point. That is what
this run shows: section 16's first posterior runs reached 4.0 A because
the uniform in-plane draw and the balanced frontier were still supplying
the exploration the direction draw no longer had, and this step removed
both. The old samplers were an ad hoc annealing, which is why their
removal (section 15's argmin in-plane, this run) traps the search every
time. The streptavidin failures read the same way: class averages have
high power, T is large, every draw is the argmin from the first iteration.

**What follows.** The design has this case in 4.4: the posterior at T is
the right weight for the reconstruction, but the search from a random
start needs a temperature above the physical one, falling to it as the
stage settles; one constant. Section 14 calibrated the old exploration at
beta 0.05-0.1. The next step is the constant as a run-time parameter
applied to every draw (the in-plane, direction, class and neighbourhood
draws share T_eff already), two restarts each at beta 0.1 and 0.3 with
the diagnostic on, read off on the fraction at the row mode at a stage's
first iteration (0.63 in the working regime, 0.80 here) and on stage 5's
end resolution (5.1 A in the 4.0 A run, 7.6 here); then, if a constant
works, the schedule rising to 1 along each stage on the measured pose
change, so the last iterations weight at the physical posterior. The
self-term stays on the list as a second-order correction (cross-half
scoring from the NU stage on, section 13), after the annealing is in. The
known-good fallback is section 16's configuration, which is not the
design.

**Built 2026-10-23: `prob_beta`.** The tempering constant is the run-time
parameter `prob_beta` (default 1.0, the physical posterior; below one
hotter; an advanced search input of `refine3D`, `solve3D`, `refine2D`,
`solve2D`, `solve2D_chunks`), applied in `pftc%ptcl_posterior_temp` for
every draw at the fills and the strategies and on the stored T_i at the
three assignments, so every draw sits at the same T_eff; the constant is
gone. The 3D tester pins the semantics (temps 4 at beta 0.5 draws as
temps 2 at beta 1, draw for draw). The schedule (beta rising to 1 along a
stage) waits until a constant is calibrated on PfCRT: beta 0.1 and 0.3,
two restarts each, diagnostic on.

**beta = 0.1: PfCRT, 2026-10-23 (`prob_beta0_1`, run 9).** Three
restarts finished, a fourth in stage 5 on the same track: 4.03 / 7.22,
4.03 / 7.22, 4.03 / 7.22 A (FSC 0.143 / 0.5), 105-107 iterations each
(stages 20/20/17/15/8-10/-/-/25). Stage 5 ends at 5.0-5.8 A on the
resolution exit (7.6 A frozen at beta 1, 5.1 in the 4.0 A run), stages
6-7 are omitted because FSC 0.5 is already past their limits, and stage 8
climbs 5.9 to 4.03 A in about seven iterations and sits there for the
remaining eighteen. The search is hot: the fraction drawn at the row
mode is 0.01-0.07 in stage 3, 0.08-0.13 in stage 4, 0.04-0.14 in stage 5
and 0.15-0.29 in stage 8; the overlap never passes 0.16; the committed
pose moves 90 degrees per iteration at the start of stage 3 and still
7.5 degrees at the end of stage 8, about one grid step. The
reconstruction is an average over that neighbourhood and does not care
at 0.143, so three restarts give the same number; it does care at 0.5
(7.22 against 6.21 and 7.06 for the two beta = 1 restarts of section 16
and 6.90 for the baseline). Trajectory: `analysis/run9_beta0.1_trajectory.txt`.

*The ceiling.* 4.03 A is shell 77 (box x smpd = 310.5 A, so shells are
physical at every crop); stage 8 runs at crop 160 of 300, Nyquist 3.88 A
(shell 80), and the NU ladder's finest rung is 3.98 A (shell 78). Every
run of the series that "reaches 4.0" is at this ceiling; section 16's bar
of 3.6 A cannot be met or measured inside this ladder, and FSC 0.5 is
the number with headroom. The bar needs restating against what stage 8
can align, or the measurement moves to the refinement that follows.

*Reading.* beta 0.1 is the regime that searches: three of three at the
ceiling where beta 1 froze at 5.35. It is also hotter than the posterior
the data support in the late stages (the pose wanders a grid step at 4 A),
and the 0.5 FSC pays for it. That is the case for the schedule the plan
names: beta of order 0.1 while the stage is settling, rising to 1 over
its last iterations so the map is weighted at the physical posterior
once the search has found its basin; the read-off is FSC 0.5 and the
fraction at the mode in the last iterations. beta 0.3 is still worth its
pair for the bracket.

*Terminology (Hans's question, 2026-10-23).* At beta below one the draw
is from the tempered posterior, p(theta|X)^beta normalized, the power
posterior of annealing and parallel tempering; it is the posterior only
at beta 1. Two readings of beta < 1 are both exact: annealing, the same
function used as search dynamics rather than as weights; or the posterior
under a noise model whose variance is inflated by 1/beta, since beta T is
what a sigma2 larger by that factor gives. The second reading is why the
known-pose calibration (gate item 3) still matters: if the physical T is
overestimated by some factor, that part of beta is a correction and the
ramp's end should be that factor, not 1. The key `prob_assign=posterior`
names the assignment rule (one draw per particle from its own row at a
data-derived temperature, no reference forced to be populated) against
the reference-driven frontier; the form, the per-particle normalization
and the temperature's origin are the posterior's at any beta, the
exponent is the knob, and the ramp's last iterations are the posterior
proper. The policies and this note say "tempered posterior" wherever a
beta below one is in play.

**Built 2026-10-23: the default and the ramp.** `prob_beta` defaults to
0.1 (Hans, on run 9). `prob_beta_end` (default 0: a flat run at
`prob_beta`) is the factor a refine3D run reaches at its last iteration by
a geometric ramp from `prob_beta` (`prob_beta_at` in utils: beta_k =
beta0 (beta1/beta0)^((k-1)/(n-1)) over the run's n iterations, so every
iteration cools by the same factor, x1.1 per iteration for 0.1 to 1 over
25; a one-iteration run is at beta1; tested); `exec_refine3D` sets the
iteration's factor on the master's parameters and on the command line,
the distributed job description carries it to the workers, one log line
per iteration. solve3D sets `prob_beta_end=1` on its last stage under
the key unless the user set it; the other stages run flat, since they end
on their budgets or on the resolution exit at an unpredictable iteration
and run 9 did not say their ends need cooling. The last stage has no exit
but its budget (the overlap bypass holds), so its ramp always reaches 1
and the final iterations weight the map at the physical posterior. 2D
does not anneal. Run: PfCRT at the defaults, three restarts; the
read-off is FSC 0.5 at the end (7.22 flat at 0.1; 6.21 and 6.90 the
numbers to beat) and the fraction at the mode over stage 8's last
iterations, which should climb from about 0.25 toward the 0.9 of the
beta = 1 regime without the resolution going backwards. If the last
iterations re-freeze and 0.5 does not improve, the cold end is meeting
the self-term, and cross-half scoring (section 13) earns its turn.

## 18. After the posterior path is stable: one temperature, no thresholds (2026-10-23)

Hans's standing instruction: once `prob_assign=posterior` with a
calibrated `prob_beta` reproduces and beats production on PfCRT (section
16's bar) and the 2D regression, the hand-set neighbourhood sizes,
angular thresholds and sampler constants go. One temperature decides
every draw; nothing else should. The inventory, in the order to retire
it, each step with its own regression and nothing deleted until every
caller has moved:

1. **Already replaced by T_eff in the tables**: the in-plane angular
   window of the fills (`prob_athres` and `dist_inpl` through
   `inpl_likelihood_nsample`, `calc_athres` for the in-plane,
   `likelihood_nsample`), the likelihood samplers of the tables
   (`gen_likelihood_val` and its bank form, `sample_likelihood_dist(s)`),
   the balanced frontier and its draw width (the frontier's `prob_athres`),
   the reference check for flat rows that only the frontier needs. Deleted
   from every table on 2026-10-24 (section 19, items 5-6;
   `sample_likelihood_index`, the function-form `sample_likelihood_dist`
   and the neighbourhood fills' likelihood in-plane draw with them); what
   remains of the list is `gen_likelihood_val` and its bank form in the
   pftc, with bench and tester callers only, and `prob_athres` in the
   strategies and the dense fill's refinement count (item 2). The 2D
   `nhood_sz`/`NHOOD_FRAC` now size only the shift refinement (item 2). The
   SHC strategies' samplers (`power_sampling`,
   `squared_sampling`, `gen_prob_objfun_val`, `angle_sampling`,
   `gen_prob_power_objfun_val`, `sample_bounded_dists`,
   `sample_power_dists`) stay: the posterior is the tables' rule
   (section 19, item 3).
2. **Refinement counts, driven by the posterior mass**: tried and taken
   back on 2026-10-25 (section 19, item 10). The refined set is a search
   budget: the row it would be sized from is conditional on the seed shift
   and the previous band, and the shift refinement is what re-ranks it.
   `projs_ns` from `calc_num2sample(prob_athres, dist)`, `npeaks_inpl`,
   `nhood_sz` and `smpl_ncls` stay. If a posterior-based rule is tried
   again it has to size the budget from the pose uncertainty (the shift and
   direction change the previous iteration measured), not from the
   sharpness of the unrefined row.
3. **Neighbourhood construction**, the cost control that decides what is
   scored before anything is scored: the stochastic bound
   (`extremal_decay`, `neigh_frac`, `snhc_nrefs_bound`, `MIN_STOCH_FINE_REFS`),
   the subspace cells and peak counts of `prob_neigh_mode=state|geom`
   (`npeaks`, the pooled capacity), the 2D `snhc_nrefs_bound`. The honest
   replacement is the previous iteration's posterior mass deciding which
   cells to re-score; a larger change, last.
4. **Then the keys themselves**: `prob_assign` went on 2026-10-24 (one
   draw on every table; section 19, items 5-6); `prob_inpl`, and the
   refine modes that only existed to pick a sampler
   (`snhc_smpl`, `snhc_smpl_many`, `inpl_smpl`, `greedy_smpl`, and
   `refine=neigh`'s sampled in-plane), collapsed onto the posterior path.
   What stays is what is not a sampler: the resolution ladder, the shift
   bounds (`trs`), the continuous polish and its bounds, `nspace`.

Not before: the calibration of `prob_beta` (section 17), the schedule if
a constant is not enough, the PfCRT 10/10 at the bar, the 2D regression,
and the known-pose calibration of the temperature (gate item 3), which
is what makes the one number mean what it says.

## 19. Decisions of 2026-10-24

Hans, with the code "a little messy" and several failure modes under test
(streptavidin `solve3D_cavgs` fails 1 of 10 on the current code; its best
results ever came from commit 43de9b34, 2026-10-18, 74 commits back):

1. **`prob_beta` defaults to 0.3.** The runs optimise for small membrane
   proteins on the assumption that the rest follows; 0.1 searched PfCRT
   to the ceiling three of three (section 17) but is the hotter end.
2. **The balanced reference assignment goes.** Methodologically a
   no-brainer: particles are assigned to references, references are not
   forced to be assigned particles, on every 2D and 3D table route.
   The frontier was first kept reachable as `prob_assign=balanced` for one
   more round of testing; Hans: "delete it now" (the same day), and item 5
   says what was deleted and what stays.
3. **The posterior sampling is the tables'.** It is restricted to the
   prob modes of the later stages of abinitio2D and abinitio3D
   (`refine=prob|prob_snhc` in 2D, `prob|prob_neigh` in 3D). The
   migration of the SHC strategies (`snhc_smpl`, `snhc_smpl_many`,
   `inpl_smpl`, `greedy_smpl`, `refine=neigh`) onto the posterior
   primitive, built on 2026-10-23 (section 16), is taken back: the six
   strategy files are restored to their state before it, `inpl_srch_peaks`
   loses `l_keep_unrefined`, and `draw_posterior_scores` goes with its
   only callers (its test too). Section 18's inventory is read with this
   in mind: the strategies' samplers are not in scope; the tables' are.
   solve3D's stages 1-2 are `prob_neigh` too (mode `shc`, from random
   orientations), so they drew from the tempered posterior like the later
   stages. Hans: "This might be the problem I have been seeing. Can we
   put them back to the original sampling approach for now?" Items 5-6.
4. **The streptavidin regression** is not attributed. Between 43de9b34
   and the posterior work sit changes to the stage control that a bisect
   should visit first, each with 10 restarts of the (fast) cavgs run:
   24d50aabf (the lpexit gain guard made able to fire), 0c9542e50 (the
   modified stage progression, nu_stage), a29f0e00c and fd9cfde6f
   (freq_march_aggressive), f078cf683 (2D tables, 3D table, NU label
   optimisation), af13bc828 (the search's crop reused in the 3D
   reconstruction), f4d6accff (the frequency marching with vol1), then
   the posterior commits 0b723362c onwards. The posterior path is off by
   default before this commit, so a failure at balanced points at the
   stage control, not at the draw.
5. **What was deleted, what stays (the same day).** The class-balanced
   frontier is deleted from every table: the dense 3D table
   (`simple_eul_prob_tab.f90`: `ref_assign` is the posterior draw), the 2D
   tables (`simple_eul_prob_tab2D.f90`: `ref_assign`, `ref_assign_sparse`;
   the `ninpl_smpl` plumbing, the in-plane work arrays of the sparse
   workspace and the likelihood branches of `score_seeded_class` with it)
   and the neighbourhood table (`simple_eul_prob_tab_neigh.f90`:
   `ref_assign_neigh` re-seeds a particle without a candidate and draws;
   the graph, the frontier, the greedy state labels and the leftover rule
   are gone). With them: the reference tests
   (`test_ref_assign_against_reference`, `reference_assign`,
   `test_dense/sparse_assign_against_reference`, `reference_assign2D`,
   `reference_sparse_assign2D`), `sample_likelihood_index` and the
   function-form `sample_likelihood_dist` (no caller left; the array form
   stays for the pftc samplers), the 2D bench's frontier timings
   (`prob_assign_profile2D`: reads and draws), and the `prob_assign` key
   itself (one assignment per table mode, no key; the tables' constructors
   require `objfun=euclid`, the loss the posterior is defined by). The
   convergence bypass covers `refine=prob` and `prob_neigh` over the
   subspace neighbourhoods (item 6).
   Hans, on a first version that kept the neighbourhood table's frontier for
   solve3D's stages 1-2: "I do not want reference-particle association to be
   balanced but be done in the same style as the posterior route even if the
   sampling of the orientations is not." Item 6 says what the stages 1-2
   question (item 3) came to.
6. **The early stages are a hill climb; the draw over its proposal set is
   not a search (streptavidin, ten restarts at `prob_beta=0.3`, the
   posterior everywhere).** Stages 1-2 did nothing in any restart: pose
   overlap 0.000-0.002 over their 40 iterations, the pose moving 85-90
   degrees per iteration (what random pairs give), the stage-1 score flat
   at e^-1 and the stage-1 map a round blob. Stage 3 started from random
   poses (at its first iteration the pose stage 2 handed over was the
   data's mode for 1% of the particles) and solved the structure from
   scratch, nine times to the Nyquist limit of the 88-pixel crop (3.12 A,
   a ceiling) and once to 3.27 A. Two causes combine.
   - The shc fill scores references in random order and stops at the first
     that beats the current pose, the current pose scored last. That is a
     hill climb's proposal set: taking the improvement is the climb's move,
     a hot posterior draw over it is a move to a random proposal.
   - The temperature is small where the band is small. T is about the
     band's Fourier component count over L_min (`T x L_min` 356 and 406
     against 347 and 396 components in stages 3 and 4 of that run), so
     stage 1 (rings 2-5 on a 137 A box) sits near T = 33, T_eff = 10 at
     beta 0.3, where even the final map separates a typical pair of
     directions by a third of a nat.
   PfCRT at beta 1 (section 17, run 8) kept stages 1-4 on track because its
   draw was cold enough to pick the improvement; run 9's 90-degree pose
   moves at the start of stage 3 at beta 0.1 are streptavidin's effect on a
   particle whose stage 3 recovers.
   Two answers were tried the same afternoon. The first kept the
   neighbourhood frontier for stages 1-2; Hans: no balanced frontier in any
   assignment ("if all the images represent one projection direction, they
   should be assigned to that projection direction, nothing else makes
   sense"). The second had the direct stochastic neighbourhoods take their
   best candidate, the hill climb's move; Hans: "we can never do argmin
   here ever. It simply doesn't work." So every table draws, and the two
   causes are met where they are:
   - **The current pose is a candidate.** `evaluate_direct_stochastic_refs`
     scores the previous reference when the stop rule left it out, so the
     draw over the store is an acceptance step over the proposals and the
     incumbent, not a forced move.
   - **The temperature of the early stages is not decided.** A first
     version set `prob_beta=1` on the stages before `PROB_REFINE_STAGE`
     (beta 1 read as the working point of PfCRT run 8's early stages);
     Hans: a run at 1, then 0.3, then back to 1 on the last stage is not
     elegant, and "I am not even sure we need prob_beta=1 for the first two
     stages." It was an inference, not a measurement (the diagnosis run
     had both causes at once), so it is out: every stage runs at
     `prob_beta`, and the streptavidin tally with the incumbent in the set
     is the experiment that separates the two causes. If the early stages
     stay dead, the principled next step is not a stage schedule but the
     draw sharpness as the one control (section 18, item 2, made literal):
     per particle, the temperature at which the row posterior has a fixed
     effective number of candidates n_eff = 1/sum p^2, colder than physical
     where the data are flat (stage 1) and hotter where they are a delta
     (stage 5), the final map weighted at the physical posterior. That
     replaces `prob_beta` and `prob_beta_end` with one number.
7. **The streptavidin tally (2026-10-24 evening; `for_claude/streptav_posterior`,
   `analysis/frontier_vs_posterior_beta0_3.txt`).** Ten restarts of the
   posterior-only build (811934ca5: every table draws, the incumbent in the
   shc candidate set, beta 0.3 flat, the last stage ramping to 1) against
   ten of the frontier (the 7dfea61a5 build at its default, which turned
   out to be the pre-posterior production path on every stage; a baseline,
   kept as `frontier_baseline`). 2000 particles, 88-pixel crop, ladder
   8.0 to 4.5 A, stages 6-7 skipped on the resolution exit.
   - Both 10/10. Final full-box maps 3.71-3.81 / 3.19-3.35 A (FSC 0.5 /
     0.143); the posterior at 3.71 in 6 of 10, the frontier in 3 of 10: one
     shell, not a difference. Stage by stage the same resolutions (stage 3
     6.86 A both, the posterior's 0.143 at 4.57 against 5.35).
   - Stages 1-2 are alive with the incumbent in the set (pose overlap 0.19
     and 0.12 against 0.000-0.002 in the diagnosis run): the forced move
     was the dominant cause. But hotter than the frontier: mean pose change
     65 and 69 degrees per iteration against 50 and 27, row-mode fraction
     0.50 and 0.38; stage 2 barely improves on stage 1 and hands stage 3 a
     nearly random start (overlap 0.04 against 0.21), which stage 3 (T in
     the hundreds) recovers in full. The early stages work and contribute
     little at beta 0.3: T about 33 there. The row-mode fraction at each
     stage's first iteration, 0.22 (3), 0.70 (4), 0.09 (5), 0.25 (8), is
     the exploration beta 0.3 buys at each band: hot at the band jumps,
     cold inside stage 4. The n_eff control of item 6, not a stage
     schedule, is the answer to both.
   - Cost: 426 against 323 s per restart (+32%), all of it stage 8: 26
     iterations against 7. The resolution (3.81/3.35 A) is reached at the
     stage's second iteration and does not move for 24 more; the overlap
     climbs with the ramp and crosses 0.95 at beta 0.90, iteration 22 of
     25. The ramp spread over the stage's budget and the pose-overlap
     bypass together make the last stage always run its budget, for
     nothing here.
8. **The annealing tail (the same evening).** A run with `prob_beta_end`
   searches flat at `prob_beta` until neither FSC 0.5 nor FSC 0.143 gained
   on the previous iteration (`check_conv3D`, the gain test of the
   resolution exit on the previous iteration's STATS_FILE record; the
   tail's first iteration is the convergence object's, since the in-memory
   strategy rebuilds `params` from the command line every iteration, which
   a first version keeping it in `params` learnt on a cavgs run), then
   `exec_refine3D` ramps the draws geometrically to `prob_beta_end` over
   `PROB_BETA_TAIL` = 5 iterations and the run ends when the tail is
   complete (`converged`, so the combined even/odd final iteration
   follows as usual, at `prob_beta_end`); a budget too short for a stall
   anneals over its last five iterations regardless. Flat runs
   (`prob_beta_end` = `prob_beta`, the stages before the last) are
   untouched: the resolution exit and the budget end them. On the
   streptavidin trace the last stage would end at its 8th iteration
   (stall at the 3rd, tail 4-8) plus the combined iteration, against 26.
   `prob_beta_end` keeps its meaning; what it took the run to get there is
   no longer the budget.
9. **The "over-sharpened" final maps (2026-10-25).** Hans: the
   `rec_final*pproc*` maps look over-sharpened, the B-factor policy upweights
   the high resolutions too much. Measured on the shell spectra of the
   unfiltered pair, the regularized final map and its postprocessed version
   (`streptav_posterior/analysis/pproc_spectra.py`, `pproc_leak.py`):
   - Streptavidin (box 128, cutoff 3.27 A): the regularized map's amplitude
     is the unfiltered pair's times the FSC, shell for shell (0.57 at FSC
     0.57, 0.20 at 0.20), a stronger attenuation than RELION's
     sqrt(2FSC/(1+FSC)); the windowed B of -47 leaves the postprocessed map
     a Guinier slope of +70 to +80 (RELION flattens to 0); beyond the cutoff
     the postprocessed map holds 0.1% of the signal band's power. Not
     over-sharpened by any spectral measure.
   - PfCRT (box 300, crop 160 for the search, cutoff 3.98 A): the
     regularized final map holds 6.7 times more power beyond the cutoff
     than in its 10 A to cutoff band. Shells of FSC at or below zero pass
     the unfiltered noise nearly unattenuated (amplitude ratio 0.94-0.99),
     shells of FSC 0.007 are shrunk (0.2-0.5): the opposite of the Wiener
     order, the signature of `add_invtausq2rho`'s fallback. tau2 = 0.001
     sig2 with sig2 = 1/mean(rho) falls under the absolute TINY (1e-10) once
     rho passes 1e7 (the Euclidean 1/sigma2 weights, many particles), and
     the fallback ridge min(1e3, 1e3 rho) = 1e3 is nothing against such a
     rho: no regularization. Streptavidin's 2000 particles keep rho below
     the threshold. The postprocess then sharpens those shells by the capped
     gain and the Butterworth lets 6% of the signal band's power through as
     pure noise over the many voxels beyond the cutoff: the grain read as
     over-sharpening. The fix is scale-free (`wiener_ridge`: 1/(tau tau2)
     whenever tau2 > 0, a thousandfold the voxel's own density otherwise),
     with a unit test from rho 1 to 1e9. The B-factor policy stays; the
     references of every Euclidean ML-regularized refinement with a large
     rho carried this noise beyond the FSC cutoff too, within the matching
     low-pass of their stage.
   With that the neighbourhood table's likelihood in-plane draw went too
   (`score_seeded_ref` posterior only; `calc_athres` for `dist_inpl` and
   the `nsample` plumbing out of the neighbourhood fills; `prob_assign` and
   a short-lived `prob_inpl_draw` gone). `gen_likelihood_val` and its bank
   form keep their bench and tester callers for now. No key decides any of
   this: the stage plan does.
10. **The refined set by posterior mass: tried, broke PorF, taken back
   (2026-10-25/26).** Hans, on `calc_num2sample(prob_athres, dist)`: could
   it be its upper bound, or something more elegant? The bound (nspace/18
   for every particle, every iteration) would have refined three to ten
   times more late in a stage than the pose-change rule and kept
   `prob_athres`. The rule tried instead (`posterior_refine_set`, commit
   9ada6a44): take the row's entries in ascending loss while each holds at
   least 1% of the row's tempered posterior mass, until 99% is covered,
   the mode always; per state block in the dense table, per state in the
   neighbourhood table, over the classes in 2D. The reasoning was that
   refinement should be spent where the draw can land, so a sharp row
   refines one entry and a flat row its mode.
   PorF (class-average start, 37,502 particles, box 180, ten restarts on
   aa98b09f against two on 9ada6a44): the posterior-mass build ends at
   8.43/7.09 A (FSC 0.5/0.143) against 5.6-7.1/4.6-5.1 A for all ten
   before. Stages 1-2 are the same (the shc neighbourhoods at 39 and 22
   A). Stage 3 (the first dense stage, 16.4 A) ends a little worse (11.1
   against 9.8-11.6 A) with the pose overlap climbing faster (0.54 against
   0.40 at its tenth iteration): the draws agree with each other sooner
   and on a worse answer. Stage 4 sits at the worse end (10.4-10.8 against
   8.9-10.4). Stage 5, where the band jumps from 8 to 4.8 A, is where it
   breaks: 7.8-8.0 A against 4.9-6.8 A in 8-12 iterations.
   The learning, plainly: the unrefined row's posterior is not the refined
   row's. Its losses are at one seed shift (estimated at the previous
   pose) in the previous band's frame; after a band jump the posterior is
   sharp (T has grown) around poses fitted at the coarser band, and only
   refining many candidates lets the shift search re-rank them. The
   refinement is the search, and the sharpness of a row that has not been
   refined is not evidence of where refinement would move its mass; a flat
   unrefined row most often means the seed shift is wrong, which is the case
   for refining more, not less. The pose-change rule measures the right
   thing (how uncertain the poses still are) and sizes the budget from it.
   So the sizing is back as it was on every table (`calc_num2sample`,
   `npeaks_inpl`, `nhood_sz`, `smpl_ncls`), `posterior_refine_set` and its
   test are deleted, and the diagnostic's k99 (median 2, measured on the
   refined row) is understood as a statement about the finished row only.
   What stays from the day: the likelihood in-plane chain of item 9's last
   paragraph is gone (`gen_likelihood_val`, its bank form and samplers; the
   bench and testers read the posterior kernels' rows through
   `gen_euclid_dists`), and the `refine3D_states` key `local_inpl_bound`,
   which only set `prob_athres` for a `prob_neigh` run and so had been
   dead since item 5.
11. **The refined count is a false-positive control, and the shift prior
   experiment (2026-10-26/27).** Hans ran the cap up: `prob_athres=20`
   (three PorF restarts: 7.09 A failing from stage 3, 4.95, 4.66) and 40
   (two: 5.03, 4.80, stage 5 stalling two or three iterations early at
   5.7-6.0 A, the dense fill 1.8x the wall time), against ten at 10
   (4.59-5.11 A). The bench histogram of where the entry lowest after
   refinement sat in the unrefined ranking (`prob_fill_profile%rank_counts`,
   the `PROB_TAB` bench files) says why: at stage 3's first iteration with
   the cap at 222 of 1000, 31% of the particles keep the unrefined argmin,
   77% land within the top 22 ranks and about 86% within the old cap of 55;
   the rest spread over ranks 56-222 at one or two percent per decile, and
   at stage 5's first iteration (cap 168 of 2500) the tail is flat at
   2.1-2.7% per decile. A flat tail is winners by chance: where the
   unrefined loss carries no information, a win after refinement is a wrong
   reference whose shift and in-plane search found a noise optimum, and
   the largest such gain over K refined wrong candidates grows like log K,
   at these temperatures one or two nats, the size of the gap between the
   true pose and the best wrong one. So the count is a false-positive
   control, not a search budget: 55 catches the real wins at ranks 2-50,
   222 or 555 hand a few percent of the particles to a lucky wrong
   reference each iteration. Ten degrees stays.
   What would let more candidates be refined without the lottery is to
   take the noise gain away from the wrong references, which have no reason
   to sit near the particle's shift: the shift prior, `shift_prior=yes`,
   an experiment (default `no`, the proven path untouched). The tables'
   shift searches and seed estimates carry a Gaussian on the shift's
   distance from the search's start (the seed converted to the reference's
   frame; zero for the seed estimate), folded into the objective as the
   factor exp(-|s - s0|^2 / (2 sigma^2 T)) on the correlation exp(-L), so
   the search minimises and the row holds the MAP loss
   L + |s - s0|^2 / (2 sigma^2 T) (`pftc_shsrch_grad%set_shift_prior`,
   applied in the cost, its gradient, the all-rotation objective and the
   coarse grid). The polish of the committed pose under `inpl_cont=yes`,
   the last optimizer to touch the shift, is the joint route over
   (sx, sy, rotind) within the full `trs`: it carries the same term on its
   raw loss, centred at its own start, the committed shift, armed per
   particle by the prob strategies' `prep4prob`/`prep4srch` with the width
   `prep_strategy3D`/`prep_strategy2D_glob` measured once per pass;
   without that the polish would have undone the tables' regularisation
   at the last step. T is the particle's
   untempered temperature, since the prior is part of the loss the tempering
   acts on; sigma^2 is the variance per axis of the shift increments the
   sampled particles made in their last search (`shincarg`, half its mean
   square; a particle never searched counts as zero, which a particle
   record cannot tell from a search that kept its shift;
   `shift_prior_sigma2`), floored at a pixel or a quarter of the matching
   wavelength, whichever is wider: below a pixel the prior would only fight
   the sub-pixel moves, and since increments measured under the prior can
   never exceed it, a floor that ignored the band could pin the width at a
   coarse stage for good. Both are in the original sampling's pixels, the
   unit of the search shifts and of `shincarg` (the polar kernels shift in
   the full box's frequencies whatever the crop). Until 2026-10-27 the
   wavelength was taken in cropped pixels, `lp / smpd_crop`, which solve3D's
   lp-tied crop holds near 3, so the band term never beat the one-pixel
   guard and sigma sat at one original pixel from the first shift stage on
   (stage 3, after the two shift-free stages left zero increments); that
   cost 3D, where stage 3 is the first correction of the shifts the class
   averages handed over, and passed in 2D, which re-centres its classes.
   The same day the multi-state table's quadratic grid step
   (`try_grid_fit`) began re-scoring its fitted shift with the prior: it
   had used the raw value, so a kept fit was accepted against a penalised
   grid minimum and entered the row unpenalised beside penalised entries
   (`test_shift_prior_grid_fit`). At stage 5's start on PorF the increments average 2 px, so a
   wrong reference wandering 2-3 sigma pays 2-4.5 nats, the size of the
   gain it was getting for free, and the true reference near its seed pays
   nothing. Hans's reading, and Cyril's on Tat and SLC: the shifts are the
   problem for the small particles.
   The readout, PorF first: three restarts at `shift_prior=yes` with the
   cap at 10 against the ten at 10; then three at `prob_athres=20
   shift_prior=yes`, which failed one of three without the prior. The
   histogram's tail should stop being flat before the resolution shows
   anything; if it stays flat the prior is not doing what it is for and
   comes out again. `SHIFT INCR ARG` in the logs is sigma's own trace.
   **The frame (Hans, 2026-10-28):** with a nonzero seed, changing rotation
   could move the stored particle-frame shift without the intended
   displacement penalty. Right: the optimizer's shift at rotation irot is
   the stored shift rotated into the reference's frame, x = R(irot) u, and
   the prior was centred at the start x itself, so at the start's rotation
   the penalty was the stored displacement but at any other rotation the
   same x is a different stored shift, moved by 2 |u0| sin(dtheta/2) for
   nothing, while the x that kept the stored shift paid that amount; the
   particle-frame lattice compared its stored-frame points against the
   reference-frame centre outright. All of it vanishes at u0 = 0, which is
   why the seed-estimate step and the re-centred 2D route never showed it
   and the class-average seeds of stage 3 onward did. Fixed the same day:
   the centre is the stored shift u0 = R(irot_in)^T x_start
   (`set_prior_centre`) and the penalty at rotation irot is
   |x - R(irot) u0|^2 = |u - u0|^2, per rotation in `objfun_vals` and the
   coarse grids (a cos/sin table per object), at `cur_inpl_idx` in the cost
   and gradient, on the stored shift directly for the particle-frame lattice
   and its fitted point, and at the continuous angle on the joint route,
   whose rotation gradient gains the term -2w (x - c0) . dc0/drho. Readouts
   `shift_prior_penalty(x, irot)` and `shift_prior_penalty_ptcl(u)`;
   `test_shift_prior_frame` checks the four cases (seed, same stored shift
   at another rotation, same x at another rotation, zero seed) and that a
   narrow prior keeps the stored shift at a nonzero seed whatever rotation
   the search ends on; `test_shift_prior`'s factor identity became the
   argmax of the raw objective times the rotation's factor. Every PorF and
   exp_gate number above taken at `shift_prior=yes` predates the fix.
   **Default (Hans, 2026-10-29): `shift_prior=yes`**, after the controlled
   comparison of 20.2 ran with it on (PorF 10/10, exp_gate 9/9 at
   prob_neff=0); `shift_prior=no` turns it off. The same day the tempering
   keys left the UI: `prob_beta` and `prob_beta_end` are set by the
   workflows (the solve controllers per stage, the refine commanders per
   iteration) and `prob_beta_ini` was removed.

## 20. Annealing coupled to importance sampling (proposal, 2026-10-26)

Hans, after the streptavidin tally and the PfCRT fix: the implementation
is crude. cryoSPARC's deterministic annealing is coupled to its importance
sampling, so that more and more particles enter as the search proceeds and
the rapid early gains come cheap from few particles; here the sampling
strategy changed but the run keeps one flat `prob_beta` (0.3) until the
tail of the last stage and one flat `nsample`, and the two early stages,
on stochastic neighbourhoods, need a policy for particles that enter with
`updatecnt` 0, or they inject random orientations.

### 20.1 What the three settings are standing in for

The resolution ladder is itself an annealing schedule: T_i is about the
number of Fourier components in the band over L_min (section 17), so the
physical posterior over directions is flat at 25 A (T about 33 at
streptavidin's stage 1, a third of a nat between neighbouring directions)
and a delta at 4 A (hundreds of nats, n_eff 1.0-1.2 in run 8). `prob_beta`
below one exists only to undo the over-cooling late; at stage 1 it makes a
flat posterior flatter, which is why the early stages are alive but
contribute little at 0.3 (19.7: pose change 65 and 69 degrees per
iteration against the frontier's 50 and 27, row-mode fraction 0.50 and
0.38, stage 3 handed a nearly random start). The tail (19.8) exists
because the final map must be weighted at the physical posterior, not at
the search's. And `nsample` is flat because nothing ties the particle
count to what a stage can use. cryoSPARC rides its noise annealing and its
batch growth on its resolution schedule; the analogue here is to ride the
draw's sharpness and the particle count on the stage ladder, with nothing
set per stage by hand.

### 20.2 The draw's sharpness: n_eff, not beta

The quantity to hold fixed across the ladder is the exploration of the
draw, measured in units that do not depend on the band or the particle
count: the effective number of candidates a particle draws among,
n_eff = exp(H) of its row posterior at the draw's temperature (the
diagnostic's `neff`; the entropy of a tempered row falls monotonically with
beta, dH/dbeta = -beta Var(T d), so the solve is a bisection). The
control: one global beta per iteration, chosen so that the mean n_eff over
the particles equals a target N* (2 to 4; the row-mode fraction the logs
carry since 2026-10-24 is the readout, 0.5 at stage 1 and 0.2-0.7 at the
later stages' first iterations under 0.3, 19.7). Global, not per particle:
a well-determined particle keeps a sharp row and an ambiguous one a broad
row, only the overall level is pinned; a per-particle n_eff would force a
flat row to concentrate on its noise-ordered top, which is the hill climb
by another name.

What it does on the ladder: at stage 1 the physical posterior is flat, so
holding N* = 3 comes out colder than physical (beta well above 1) and the
particle goes to one of its best few references: the climbing the early
stages lost without being an argmin. At stage 5 it comes out hotter than
physical, which is what 0.3 was doing, and by exactly the amount the
row's sharpness there asks for. As the particle count grows and the map
sharpens (20.3) the control cools on its own.

Implementation: the master holds the rows at assignment time (`ref_assign`
of the dense table, `ref_assign_neigh` over the candidate store, the 2D
tables), so beta is solved by bisection on the mean n_eff over a subsample
of rows (a thousand rows, twenty passes) before the draws; the 2D tables
take the same rule. The ramp of `prob_beta_at` goes; the tail (19.8) stays
and becomes clean: once the resolution stalls, N* ramps to the physical
posterior's own n_eff over `PROB_BETA_TAIL` iterations, i.e. beta ramps to
1 from wherever the control sits, and the run ends. `prob_beta` goes;
`prob_neff` (N*) replaces it; `prob_beta_end` loses its number (the end is
physical) and becomes the switch that the last stage already sets.

**Built 2026-10-27, behind `prob_neff` (0, the default, is the flat
`prob_beta` path, bit for bit).** `prob_neff_control` in the tables' utils:
rows offered at a fixed stride through the particles, up to 1000 kept, each
with at least 2 N* valid candidates (a row of two candidates lands on one
of its best two at any beta and would pull the level of the rows where the
choice is open: the neighbourhood tables' stop-rule rows), packed above
their minimum; the mean n_eff at a beta in one pass per row (H = log Z +
beta T <d>); 40 bisections in log beta over [1e-3, 1e3]; fewer than 50 kept
rows fall back to `prob_beta` for the iteration. Every table (dense 3D,
neighbourhood, dense and sparse 2D) runs it before its draws and logs
`POSTERIOR DRAW CONTROL` with the solved beta, the mean n_eff at it and at
beta 1 (the physical posterior's own n_eff, which is section 21's first
reading for free). The solved beta goes back to the commander
(`prob_beta_control_last`) and out as `prob_beta` to the next iteration's
fills, so the in-plane draws run one iteration behind the direction draws
and `prob_beta` sets the first iteration only; `prob_beta_ini` is then
moot. The tail: `prob_beta_end` is the switch (any positive value; solve3D
sets 1 on the last stage), `exec_refine3D` passes the tail's progress as
`prob_anneal_frac` and the table ramps the solved beta geometrically to 1.
Unit tests: `test_neff_control` (synthetic geometric rows against numpy,
`neff_control_ref.py`: the solved beta for N* = 3, monotonicity, the
temperature scaling, the stride, the row minimum, the fallback, the ramp)
and the dense table's assignment under the control. What is not yet done
from this subsection: `prob_beta` is still the key the flat path reads, the
20.5 retirements wait for the tally.

Hans, the same day: solve2D and solve3D are to be consistent. The 2D route
had the control but no tail, so a 2D run ended with its draws at the
control's level (or at 0.3 before) and its last class averages were, in
expectation, a soft membership among the best few classes. Built: the tail
in `check_conv2D`, keyed on the class overlap converging past `minits`
(2D has no per-iteration resolution) or on the budget's end, the same
`PROB_BETA_TAIL` ramp in `exec_refine2D`, the tail's state in the 2D
strategies' convergence object (the shared-memory strategy kept a fresh one
every iteration; it keeps one now, as the 3D one), `prob_beta_end` on the
2D inputs, and solve2D setting `prob_beta_end=1` on its last stage unless
the user set it, before its terminal greedy pass. `test_posterior_anneal_tail2D`
covers the start past minits, the short budget completing at `maxits`, the
switch under the control and the flat run.

**The shift prior (19.11) under the schedule** (Hans, 2026-10-27: check
carefully how the two marry). They do, by construction, on three points:

- *What beta tempers.* The row holds the MAP loss L + |s - s0|^2 /
  (2 sigma^2 T) and the draw is p_c proportional to exp(-beta T (d_c -
  d_min)) = exp(-beta [T L_c + |s_c - s0|^2 / 2 sigma^2] + const): the
  tempering acts on likelihood and prior together, the tempered posterior
  (likelihood x prior)^beta. So beta never changes the balance between the
  two, at any value the control picks (a cold beta above 1 at stage 3 does
  not make the prior stiffer relative to the data, it sharpens the whole
  row); the balance is set once, in the fill, where the MAP refinement runs
  at the untempered T whatever the draw's beta (`set_shift_prior` takes
  `posterior_temperature`, not `prob_beta` times it). The alternative, a
  prior left out of the tempering, would make the prior's weight relative
  to the data a function of beta, which is what must not happen; the
  present form is the right one and the control needs nothing from it.
- *The feedback through sigma.* The prior's width is the previous
  iteration's increments, and the draw's temperature sets the increments:
  colder draws commit smaller moves, so sigma tightens as the schedule
  cools (the control's cold early betas, the tail) and loosens where it is
  hot. The direction is the right one, and where it could bite the floor
  decouples the two: in the coarse bands the band term (a quarter
  wavelength in original pixels, 4-5 px at 16-20 A) is above any measured
  increment, so the early cold draws of the control cannot pin the prior;
  in the fine bands the floor is one pixel and the prior can tighten only
  into the sub-pixel regime, where a half-pixel move costs an eighth of a
  nat and the prior constrains nothing the refinement would want. In the
  tail (five iterations, beta to 1) the increments shrink toward the floor
  and the prior with them, one iteration behind: harmless at the last
  band for the same reason. solve3D's stages 1-2 search no shifts, so
  sigma sits at the floor through stage 3's entry under either schedule.
- *The control's readout.* The prior sharpens the rows (far-shift
  candidates pay nats), so holding N* comes out at a hotter beta than
  without it, and the control does not undo what the prior did: the prior
  decides which candidates are credible, the control how many of them the
  draw spreads over. The `n_eff at beta 1` of the control's log line is
  the physical posterior's width with and without the prior, a companion
  to item 11's rank histogram: the prior is working if that number falls
  and the histogram's tail stops being flat.

The weighted insertion of 20.8 takes the MAP losses as the draw does,
prior included, at the draw's beta: nothing to add. The 2D tables carry the
prior the same way, so the 2D tail anneals the same row.

**The regression of 2026-10-28 (`for_claude/prob_neff_regression`: ten
restarts each of exp_gate, FlhB, 132k particles, and PorF, 37.5k, at
`prob_neff=3 shift_prior=yes`, commit ad443249).** Streptavidin passed;
exp_gate ended 6 of 10 at 4.01 A (cFAR 0.62-0.68) and 4 at 4.15-4.30 A
(cFAR 0.40-0.49), PorF 3 of 10 at 4.5-4.9 A and 7 at 5.8-7.1 A, where ten of
ten had been 4.6-5.1 A at `prob_beta=0.3`. Hans: the search is running too
cold. The logs, stage by stage:

- The control is not the stochastic element: the solved beta is the same
  in the good and the bad restarts at every stage, to the second decimal.
- The fate is sealed early. PorF separates during stage 2 and is decided by
  the end of stage 3 (8.4-9.2 against 9.5-10.8 A, score 0.410 against
  0.406-0.408, pose change 10.6-11.4 against 13-15 degrees), stages 4-8
  never recover it; exp_gate splits at the stage 4-5 band jump (the good
  runs 7.0 to 4.47 A in three iterations at the NU entry, the bad ones at
  6-7 A for the stage). The tail works: the good runs still gain 4.38 to
  4.01 A in it.
- Against 0.3 the control is 2-5x colder exactly where 0.3 was hottest:
  the shc stages (PorF stage 1 climbs 0.5 to 1.4 and stage 2 sits at
  1.1-1.5; exp_gate stage 2 at 0.6-1.3) and the first iterations of every
  new band (0.65-0.8 at the stage-3 entry, 0.52-0.69 at the stage-5 entry,
  both data sets), relaxing to about 0.3 only by a stage's end. The
  row-mode fraction sits at 0.70 throughout, where 19.7 had 0.22, 0.09 and
  0.25 at the stage 3, 5 and 8 entries. The one place the control is hotter
  than 0.3 is exp_gate's stage 1 (0.07-0.25: its physical n_eff is 1.08
  already at 20 A, a large particle), and exp_gate's split is not there.
  The within-stage trend is the control warming (beta falling) as the rows
  sharpen with the map, and cooling again at each band jump, where the
  rows are sharp for the wrong reason: T grew with the band and the new
  shells are noise-ordered because the reference has no signal in them yet.
  That is 20.6's risk as written, a cold draw on noise-ordered rows, at the
  argmin's ratio (5 of 10 PfCRT restarts near 7 A in section 15).
- Beside the point but in the logs: stages 1 and 2 reach their plateau
  after about 8 of their 20 iterations (score flat, pose change flat at 33
  degrees, the rest spent hopping at equilibrium); exp_gate restart 1
  reported cFAR 0.0002 for nine consecutive stage-5 iterations, that run
  only, a cFAR artefact to look at rather than a mechanism; stages 1-4 run
  on the same best 10,000 particles by rank (`greedy_sampling=yes`,
  `frac_best` 1.0, then 0.5 at stages 5-6 and 0.85 at 7-8), by design.

Hans, the same day: all of this is too complicated; fix the in-plane
`prob_beta` at 0.3 with a suitably placed constant while the annealing
schedule of the projection directions is worked out. Built: `PROB_BETA_INPL`
(0.3, `simple_defs`, beside `PROB_BETA_TAIL`, which moved there from the
parameters module the same day: constants belong in a defs module) is the
factor of every
in-plane draw of the fills (`ptcl_posterior_temp`), on every path; neither
`prob_beta`, the tail nor the control touches the in-plane temperature any
more, and the handoff of the solved factor to the next iteration's
`prob_beta` only feeds the control's fallback. The direction draws keep
their schedule as the thing under study: `prob_beta` flat with the tail
(`prob_neff=0`) or the control (`prob_neff`). A first version of the day
had put the in-plane choice behind a key (`prob_neff_inpl`, with a tail
ramp of its own) and added a solve3D `prob_neff_stage` for Hans's separate
test of the established draws through stages 1-2; both were taken out the
same day, the constant alone stays (the stages 1-2 test, if wanted, is a
run with the control off on those stages, no key needed yet). Offered, not built: the sweep of N* on PorF (6
and 10, with `prob_neff=0` as the control that pins the regression on the
control rather than on the shift prior); and the statistic the control
pins, the arithmetic mean of exp(H), which the broad rows dominate (70%
delta rows and 30% broad average to 3 with the broad rows at about 8 and
the deltas never moving), against the mean entropy (the geometric-mean
n_eff), which the delta rows pull down and which therefore solves a hotter
beta exactly where the rows are heterogeneous, the band entries, without a
stage schedule; one line in `neff_control_mean_neff`.

**The controlled comparison (2026-10-28 evening, commit 62ffd592: in-plane
at `PROB_BETA_INPL`, the prior's frame fixed, `shift_prior=yes`; logs in
`for_claude/prob_neff_regression/{PorF,exp_gate}/prob_neff{0,3}`).**
`prob_neff=0`: PorF 10 of 10 (one weaker, 6.6 A, the rest 4.5-4.9),
exp_gate 9 of 9 finished at 4.01 A, cFAR 0.63-0.70, in 76-99 iterations,
Hans: the best maps seen for that data set at this stage (restart 10 died
on a project-file I/O error, not the algorithm). `prob_neff=3`: PorF 3 of
10, exp_gate 8 of 10. Read side by side in the row-mode fraction, the one
readout both print: under 0.3 the shc stages are hot (PorF 0.39-0.52,
pose change 53-78 degrees, overlap 0.10-0.14; exp_gate's stage 1 is cold,
0.92 to 0.75, and its stage 2 hot at 0.38) and establish the shape at low
resolution (Hans, 2026-10-29: after stage 2 the shape of PfCRT, exp_gate and
PorF is essentially fully established; the pose statistics churn at the
grid level, not the shape), stage 3 enters at 0.14-0.23 (pose change
54-75 degrees) and re-derives the poses against that shape on the denser
grid with shifts, cooling to 0.71-0.73 by its end as the model sharpens;
the later entries are moderately hot (0.46-0.48 at stage 5, 0.47-0.52 at
stage 7) and cool to 0.6-0.7. The within-stage annealing is real and
data-driven: the entry heat is set by how poor the model is relative to the
new band, the cooling by its improvement. Under N* = 3 the mode fraction
is 0.68-0.73 in every stage and at every entry; stages 1-2 reach a sharper
map (PorF stage 2 ends at 11.1-12.0 A against 13.0-14.2) with the poses in
consensus (overlap 0.33, pose change 17 degrees) and stage 3 enters at
0.68 instead of 0.23 and cannot re-search them; whether that consensus is
a wrong shape some of the time or the right shape with poses the shc grid
cannot resolve, the logs cannot say (the stage-2 volumes of the failed
control restarts can). Either way the thermostat removes the signal that
carries the ladder's own annealing (20.1): the width of the physical
posterior at a stage entry is the model's estimate of how far to re-search.

**Rulings of 2026-10-29 (Hans).** The control is retired; `prob_beta=0.3`
flat with the tail is the schedule. What N* = 3 was attractive for, a
handle on how many times a particle enters the reconstruction with weights
(20.8) without overcooking the compute, is a cap, not a thermostat: it
survives as `topk` on the insertion and never touches the draw. The n_eff
machinery stays as a readout: every table logs `POSTERIOR DRAW N_EFF`, the
mean n_eff at the draw's factor and at beta 1 over the sampled rows
(`prob_neff_readout`), which is the band-independent trajectory of every
stage that the first batch lacked. The `prob_neff` key, `prob_anneal_frac`,
the solve, the handoff and the control's tail are deleted.

**One factor again (Hans, later the same day: tidy the in-plane constant
against the directions' schedule, keeping the last stage's annealing).**
With the control gone the in-plane constant had nothing left to be
separate from, so `PROB_BETA_INPL` is gone too: `PROB_BETA` (0.3,
`simple_defs`) is the one tempering factor of every draw of the tables,
in-plane, direction, class and neighbourhood; `parameters%prob_beta`
defaults to it and is the iteration's value the workflows carry (the refine
commanders put it on the workers' lines during the tail), and
`ptcl_posterior_temp` is `prob_beta` times the particle temperature again,
as before 2026-10-28. The one behavioural difference from the comparison
runs of 20.2: in the tail the fills' in-plane draws cool to 1 with the
direction draws instead of staying at 0.3; since a refined entry's rotation
is the shift search's argmax, this touches the unrefined entries and the
refined set's membership only. `prob_beta`, `prob_beta_end` left the UI the
same day (set by the workflows) and `prob_beta_ini` was deleted.

### 20.3 The particle count: the batch and the map's memory (rewritten 2026-10-27)

**What `nsample` is.** Not the particles in the map. `nsample` is the
batch: the particles searched and inserted fresh in one iteration,
f = nsample / N of the active set (`update_frac`). The trailing
reconstruction keeps the map at the mass of the whole updated population
by the population rule, A_new = (u/f) A_current + (1 - u) A_previous with
u = f, so after one sweep (N / nsample iterations, the count-tiered
sampler visiting every particle once) the map represents all N particles
and its half-map SSNR is theirs. The first draft of this subsection sized
the stage's particle count by Rosenthal-Henderson from the map's FSC and
called it `nsample`; Hans: that is excessive, `nsample` only controls how
many particles are updated. Right. What the batch governs with u = f is
the cost of an iteration and the map's memory: the map is an exponential
average over the last ~1/f iterations of partials, each with the poses as
they were then, so it lags the poses by N / nsample iterations.

**The actual cost of a flat small batch**, then, is at the fine bands and
it is lag, not resolution: N = 200k, `nsample` 10,000, a 12-iteration
stage at 8 A ends with (1 - 0.05)^12 = 54% of its map's mass carrying the
poses of the previous band. The stage's references cannot be better than
the poses they average, more than half of which were found at 12 A, and
the stage grinds its budget on a map that follows the poses at 5% per
iteration. The same arithmetic is why the early stages cannot simply be
made cheap by shrinking the batch under u = f: a batch of 2,000 moves the
map at 1% per iteration, a hundred-iteration horizon, and the poses run
away from it. cryoSPARC's early iterations are cheap *and* fresh because
its momentum is independent of the batch: the map is an average over the
last ~10 batches whatever their size, so early it represents a few
thousand fresh particles, which is all 25 A needs, and the iterations cost
a few thousand evaluations. SIMPLE has that knob already: `ufrac_trec`,
the single-state override u > f (`importance_sampling_fractional_update_policy.md`
section 6). With u > f the chain still carries mass N, but its information
is the last ~1/u iterations' batches: a map of about nsample / u fresh
particles up-weighted to the population, whose half-map FSC is that of
nsample / u particles.

**The schedule, two numbers per stage.** The batch nsample_k and the
update weight u_k, tied by two conditions:

1. *fresh*: u_k >= about 2.3 / J_k, J_k the stage's iteration allowance,
   so that by the stage's exit at least 90% of the map's mass carries the
   stage's own poses (1 - (1 - u)^J >= 0.9);
2. *resolved*: nsample_k / u_k >= N_RH(lp_k), the particle count at which
   the half maps put FSC 0.5 at the stage's band. This is where
   Rosenthal-Henderson belongs: N_RH is read from the trailing map's FSC at
   the next rung and the population it represents, N_RH(lp_next) =
   M_k (1 - FSC) / FSC, floored (a few thousand) and capped at N.

Early, N_RH is small and the two are met with a small batch at a large u:
the cryoSPARC regime, cheap and fresh at low resolution. As the band
tightens N_RH grows until nsample / u reaches N, where u = f and the
population rule holds again: the map is the whole set, the batch is what
freshness asks for, nsample_k = N u_k, and the full set once that passes
the 90% switch. The cross-over is decided by the data (N_RH against N),
not by a stage number. A stage that leaves on the resolution exit has FSC
0.5 past the next rung, so N_RH there is below its own count and nothing
grows: the exit and the growth are the same signal read two ways. Pose
error early reads as a low FSC and a high N_RH, the safe direction.

**A plausible run** (N = 200k, J = 20, 20, 17, 15, 12 ...; the ladder 25,
20, 16, 12, 8 A, NU): stage 1 on a batch of 2,000 (a floor for the
class-balanced units) at u 0.25, a map of ~8,000 fresh particles, 40,000
evaluations for its 20 iterations against 200,000 today; stages 2-3 on
2,000-3,000 at u 0.15-0.2; stage 4 (12 A, N_RH ~40,000) on ~6,000; stage 5
(8 A, N_RH ~ N) is the cross-over: u = f >= 2.3/12 = 0.2, a batch of
40,000; then the full set. Cheaper than today through the dense stages,
dearer at the NU stage, and the NU stage's map is its own. Coverage: with
batches of a few thousand through stages 1-3 most particles reach stage 4
never updated and enter under 20.4 (the dense first scan and the draw),
as cryoSPARC's do.

**The budget** is then evaluations, as before: one constant per stage of
about today's 200,000, iterations = B / nsample clamped by `minits` and a
cap, with the resolution exit and the tail still ending stages early.

**The ruling this needs.** The early regime deliberately breaks the
map-equals-population property (u > f) that the trailing chain was built
to keep, in the stages where the band cannot use the population anyway.
`ufrac_trec` is today a manual single-state override; the schedule would
set it per stage. Whether that is acceptable as the automatic early mode
is the design question, ahead of any constant. Behind a key
(`nsample_auto`), with the flat 10,000 as the control, the readout per
stage: the batch and u the rule chose, N_RH as read, the fraction of the
map carrying the stage's poses at its exit, evaluations and wall time,
the stage-3 handover and the finals.

### 20.4 Particles entering with updatecnt 0

The concern: a particle sampled for the first time in stage 1 or 2 has a
random current pose, the `shc` fill's stop rule ends at the first
reference that beats it, so its candidate set is a few random references
and the draw among them is a random orientation into the reconstruction.
The scan is already handled: `evaluate_direct_stochastic_refs` forces a
dense first pass on a particle that has never been searched
(`l_greedy_first`, keyed on `updatecnt`), so the candidate set of an
entering particle is every reference; `snhc` should do the same (its
random subset is no better for an entering particle than the stop rule).
What makes the entering particle inject noise is the draw, not the scan:
at beta 0.3 and T about 33 a draw over 500 candidates is a random pick.
Under 20.2 it is not: the particle lands on one of its best few, and the
incumbent-in-the-set rule (19.6) keeps a searched particle from being
forced to move. So the policy for entering particles is the dense first
scan plus the N* draw; the alternative of scanning the average fraction
the searched particles scan is weaker than the full scan and does not
address the draw.

### 20.5 What the existing settings become

- `prob_beta`: a consequence of N* and the row; gone as an input.
- `prob_beta_end`: the end is the physical posterior; the key becomes the
  last stage's switch, which solve3D sets.
- `nsample`: the batch, a consequence of the stage's update weight and
  the FSC at the band (20.3: nsample_k = u_k N_RH early, u_k N once the
  band needs the population); the user's `nsample` remains a cap.
- `ufrac_trec`: set per stage by the same rule (u above f in the early
  regime, f once the band needs the population), no longer a manual
  single-state override only.
- `MAXITS` per stage: a budget in particle evaluations, the resolution
  exit and the tail end the stages before it.
- `prob_beta_ini` (Cyril, 2026-10-26, commit f18ebb0f: a separate beta for
  particle stages 1-2, capped at `prob_beta`, so an equal or hotter early
  draw): the opposite sign of what 20.2 predicts for those stages (holding
  N* at stage 1 comes out colder than 0.3, not hotter), so its runs are a
  direct test of that prediction before 20.2 is built: if a hotter stage
  1-2 helps, 20.2's N* is wrong there, or the neighbourhood tables' rows
  are not the dense rows the argument was made on. Under 20.2 the key goes
  with `prob_beta`.
- `prob_athres`, `nhood_sz`, `npeaks_inpl` and the rest of section 18:
  unchanged by this step. Item 2 of section 18 was tried and taken back the
  same day (19.10, PorF): the unrefined row's posterior is not the refined
  row's, so n_eff is not a basis for the refinement budget. The n_eff of
  20.2 is measured on the finished row, the row the draw is from, which is
  the quantity 19.10 says is sound.

### 20.6 The risk, and its readout

N* small is a cold draw, and a cold draw on noise-ordered rows can trap:
the `K_inpl = 1` argmin of section 15 left 5 of 10 PfCRT restarts near
7 A. N* = 2-3 is a draw among the best few, not the argmin, but the
calibration is the same tally as always, and the readouts per stage are
in the logs already: the mean n_eff and the row-mode fraction (20.2
makes the first a reported quantity), the pose change per iteration, and
the stage-3 handover overlap, which 19.7 put at 0.04 against the
frontier's 0.21.

### 20.7 Order of work

1. 20.2 alone, everything else as it is: streptavidin ten restarts against
   the 19.7 tally (stages 1-2 alive and tightening, handover overlap, final
   numbers, wall time), then PfCRT ten against the acceptance bar.
2. 20.3, behind `nsample_auto`, after the ruling on u above f in the
   early stages; with the budget in evaluations.
3. 20.8, the weighted insertion, behind a key: the one-stage comparison,
   then the iteration counts in a full run.
4. The `snhc` first scan (20.4), one line.
5. 20.8's second experiment, the two temperatures, one number.
6. Then section 18 in its order.

### 20.8 Weighted insertion: the draw's expectation without its noise (2026-10-27)

Hans, after 20.2 was built: the next natural step is to put the n_eff
candidates into the reconstruction with weights.

**The argument.** One draw from p^beta per particle, inserted with weight
one, gives a map whose expectation over the draws is the weighted
insertion of the row at p^beta. Inserting the row's posterior mass with
its weights therefore changes nothing about what the map converges to and
removes the per-iteration Monte Carlo noise of one pose per particle
(Rao-Blackwellization). The next iteration's references are then the
expectation the present scheme only approximates, and the gain should
show as fewer iterations to the same resolution: the draw's variance is
what the n_eff control pays for.

That settles two things:

- *The weights are at the draw's beta, not the physical one.* Inserting at
  beta 1 while drawing at the control's beta would sharpen the map
  relative to the search, a different scheme (see the second experiment
  below). The tail anneals the draw and the weights to 1 together, so the
  final map is at the physical posterior either way, and section 21's
  calibration is not a prerequisite: a miscalibrated temperature affects
  the draw and the weighted insertion identically.
- *The committed pose stays the draw.* The search state for the next
  iteration (the seed shift, the neighbourhood, the incumbent in the
  candidate set of 19.6) is one pose, and a draw is the right one-pose
  summary of a broad row. Only the reconstruction takes the expectation,
  which needs no summary.

**Where it lands.** At assignment the master holds the row. Beside the one
`ptcl_ref` of `ASSIGNMENT_FBODY` it writes the entries holding 99% of the
mass at the draw's beta, capped at a few (the diagnostic's k99: median 3-5
at a stage's start, 1-2 as it settles), each with its in-plane, its shift
(the refined shift of a refined entry, the seed shift otherwise: what the
draw commits today) and its normalized weight. The reconstruction inserts
K slices per particle into the particle's half, each at p_k times the
particle weight; the insertion already takes a per-particle weight, so
this is a loop over the entries. The cost is the insertion time times the
mean K, about 1.5 early and 1 late, against a search that dominates. The
trailing blend and `update_frac` are untouched: K poses for the sampled
particles, the momentum for the rest. The 2D tables can take the same
step later (soft class membership in the class averages).

**The readout.** On one stage, the FSC of the weighted map against the
single-draw map at the same iteration from the same assignment (the same
rows: a controlled comparison), then the iterations to the stage exit in a
full run. Built after the `prob_neff` tally: the weighted insertion is the
variance-reduced form of whatever the control's beta does and inherits its
tuning.

**What it separates, and the second experiment.** The weighted insertion
makes the draw's reconstruction role redundant but not its search role,
and that splits the two halves of what `prob_beta` was doing. After this
step the annealing means two things that share a beta: how far the search
state may wander, and how blurred the map is. cryoSPARC couples them
through the noise model because it never commits a pose; here they are
coupled by choice. One knob, and the tail ends both at the physical
posterior, is probably right for now, but it is also the point where they
can legitimately separate: a search held at N* = 3 for exploration while
the map is inserted at a sharper beta gives the references the next
iteration scores against less blur without taking the exploration away.
Whether that is better or a way to lock in early is what section 15's
argmin runs warn about, so it is not to be guessed: once the weighted
insertion exists it is one number, a second key for the insertion's beta
(or its own N*), the same tally.

**Built 2026-10-29, behind `topk` (default 1: the committed pose alone,
bit for bit the single-draw run, its random stream included).** Hans's
rulings on the form: 3D only (2D out, "that I am sure of"); refined rows
only; the extra poses SAMPLED WITHOUT REPLACEMENT from the tempered row at
the draw's temperature rather than the top K, which would be the argmin's
cousin; weights the tempered posterior renormalised over the sampled set;
an extra under 1% of its particle's weight not inserted
(`PROB_INSERT_WMIN`); the committed pose always inserted. Performance
rules: no second read or transform of a particle, the extra poses are
inserted from the committed plane re-shifted by the pose difference (a
phase ramp over the stored lattice, `fplane_shift`), so an extra costs one
gridding insertion in a step that is a tenth of the iteration; the fill is
untouched; the master writes the extras beside the assignment
(`INSERT_FBODY`, `prob_insert_table`: per particle the committed weight and
the extras as reference, in-plane index, refined shift and weight) and logs
`PROB_TAB_INSERT` (insertions per particle, mean committed weight); the
worker reads the table after `read_assignment` into the part's layout
(`prob_insert_poses`, passed to the search through `strategy3D_spec`), the
prob strategy resolves a particle's extras in its `srch` just before it
commits the draw (`resolve_extra_poses`), through the one conversion rule
the committed pose uses (`candidate_pose` in `simple_strategy3D_utils`:
reference Euler pair, in-plane index's angle, stored shift = particle's +
increment), and `calc_3Drec` inserts (`update_state_half_rec`,
`insert_extra_poses`). The companion files of the alignment have one
deleter, `del_prob_align_files`. The shc/snhc neighbourhoods
of stages 1-2 insert the committed pose alone (their candidates are the
hill climb's unrefined proposals), so stages 1-2 are unaffected; gridding
backend only (a warning under `rec_backend=pcg`); the final and bootstrap
reconstructions use the committed poses. Unit tests: the sampler
(`test_draw_posterior_extras`), the table's round trip, the plane re-shift
against a plane generated at the new shift with the insertion
(`test_fplane_shift`), and a two-state table (`test_insert_extras_multistate`:
every extra in the committed state's block, refined, not the committed
reference, weights summing to 1, the file round trip, `topk=1` recording
nothing). **Multi-state (Hans, 2026-10-29): state assignment stays hard.**
The extras are poses within the committed state, drawn from its block
alone and inserted into its half-map; a particle never enters another
state's map through the insertion. Soft state membership is `flex`'s
(`m_estimator=flex`, the state weight set), whose weight multiplies the
pose weights as it multiplies the committed pose's. The readout: PorF ten
restarts at `topk=1` against `topk=3`, iterations to resolution and the
REC3D bench's gridding time; then exp_gate.

### 20.9 Review of the `prob_neff=3 shift_prior=yes` regression (2026-10-28)

**Evidence.** Reviewed the implementations at `04f6f0cb7` and the twenty
supplied Linux restart logs (PorF and exp_gate, ten each, built at
`ad443249`, which includes the original-pixel prior-floor correction).
Every run ends normally. Every control line keeps 1000 rows and reaches
the target; no fallback, bracket saturation or invalid/flat-row warning
appears. This is a reconstruction-quality regression, not a failed entropy
solve. The logs alone do not establish the correctness of a map or isolate
the two switches, which were enabled together.

The maintainer subsequently reports that shift_prior alone works, with
only 1-2 failed runs. This makes the added n_eff policy the stronger
regression suspect. That comparison is user-reported; prior-only logs
were not supplied in this review. The prior-frame finding below remains
a separate correctness issue and is not established as the cause of the
new failure rate.

| Dataset | Final FSC=0.143 resolutions by restart, A |
| --- | --- |
| PorF | 6.37, 5.78, 4.46, 6.50, 4.88, 7.09, 7.09, 6.12, 6.50, 4.80 |
| exp_gate | 4.23, 4.01, 4.01, 4.30, 4.01, 4.01, 4.23, 4.15, 4.01, 4.01 |

The PorF split precedes the final tail. At the end of stage 4, restarts
3/5/10 are at 8.21/8.67/8.67 A, and the others at 8.91-10.40 A. Final
tails finish in every log. Neither a missing tail nor the too-few-rows
branch explains these runs.

**Search-policy concern, supported by the logs but not isolated causally.**
The premise of 20.2 needs calibration on 3D. Holding mean direction
entropy at three permits a substantially colder early search than the
previous beta 0.3. PorF's median last beta is 1.4075 in stage 1 and
1.4063 in stage 2; the mean last row-mode fractions are 0.704 and 0.668.
The preceding direction beta also controls the next fill's in-plane
draws (`polarft_calc%ptcl_posterior_temp`), although the control measures
direction entropy, not in-plane or joint-pose entropy. Thus this changes
both searches, and the direction target does not constrain the latter's
breadth. exp_gate's median last stage-1 beta is only 0.2500: the same
target does not imply the same exploration relative to the old setting.
Early trapping is a leading hypothesis, not a demonstrated sole cause.

**Implementation findings.**

1. `prob_beta_last_solved` is process-wide, with no reset or run/stage
   identity. Both refine commanders read it before their first iteration,
   overriding a newly configured starting beta. This leaks between stages
   and between independent 2D/3D invocations in one process. At PorF's
   stage-3 entry the previous beta is 1.343-1.501, while the newly solved
   direction beta is 0.629-0.693. The first fill therefore uses a coefficient
   from the preceding band and candidate policy. Own the feedback explicitly
   in the invocation and define/reset the stage-entry value. The documented
   one-iteration lag also leaves the last fill at the preceding tail beta
   (0.753-0.788 on PorF), although the final direction assignment uses 1.
   A physical endpoint for both draws requires handling the fill too.

2. The prior is a Gaussian on reference-coordinate optimizer displacement,
   not consistently on the stored particle-frame shift displacement when
   the seed is nonzero and the in-plane angle changes. The dense fill
   converts the seed at the sampled angle; `minimize` saves that reference
   vector as `prior_xy0`, then rescans all angles with one scalar prior
   factor. Returning the shift rotates it at the new angle. For a
   particle-frame seed (3,0), holding the reference vector (3,0) while
   changing angle from 0 to 180 degrees returns (-3,0) with zero current
   penalty; a Gaussian about the particle-frame seed with sigma=1 charges
   18 nats. Conversely, preserving the physical shift can incur a penalty.
   This is an analytic frame counterexample, not a measured contribution
   to the supplied failures. Zero-centered searches are unaffected by this
   discrepancy. The existing prior tests verify its factor in optimizer
   coordinates and do not test this invariant. Correcting it requires
   carrying the seed in a common frame through angle scoring and gradients.

3. The too-few-rows branch of `neff_control_apply` skips the annealing ramp:
   fewer than 50 rows, beta_flat=0.3 and anneal_frac=1 return 0.3, not 1.
   This is a separate correctness defect, absent from these twenty logs.

The additive loss, its shift derivative, the untempered T in the prior
weight, and the corrected `lp/smpd` floor are internally consistent for
the stated optimizer-coordinate prior. The entropy bisection is monotone
on fixed valid rows. A separate robustness edge remains: 995 resolved rows
and five flat rows of 500 candidates have a minimum mean n_eff of 3.495;
target 3 drives beta to 1000 and makes the resolved rows greedy. It did
not occur in the supplied logs.

**Next comparison.** Retain the working prior-only configuration as the
control, on the same input, build and matched initializations, judging
the maps as well as the stage trajectories. First isolate direction-to-
in-plane feedback: keep the fills at the established beta 0.3 while
solving n_eff=3 for direction assignment only. Separately test whether
retaining the established early-stage draws and enabling the control only
after stages 1-2 restores robustness. Neither separation is an existing
CLI option. Do not change the prior and the exploration policy in the
same attribution experiment. Add regressions for independent
invocations/stage entry, the final fill's beta, the nonzero-seed angle
invariant, and fallback annealing before accepting fixes.

**Validation observed here.** Mac Release `make -C build -j8`, the complete
fast gate (15/15, 47.80 s, including polar alignment, entropy-control,
prior and convergence checks), and installation passed. No production
source was changed. The twenty Linux logs were inspected; no new Linux
or real-particle run was executed locally.

*2026-10-29: superseded.* The control was retired after the controlled
comparison of 20.2 (prob_neff=0 against 3 on PorF and exp_gate, built at
62ffd592 with the in-plane factor fixed and the prior frame corrected):
`prob_beta=0.3` flat with the tail is the schedule and n_eff is a readout,
so findings 1 and 3 fall with the control. Finding 2, the prior's frame, was
confirmed and fixed in 62ffd592 (19.11): the penalty follows the stored
particle-frame shift at every rotation.

## 21. Rulings of 2026-10-26, and the calibration on real data

### 21.1 Rulings

- **The 3.6 A bar is set aside.** PfCRT at 3.6 A is probably not reachable
  without a solvent prior or a mask, so the acceptance run of section 7,
  item 3, is ten restarts that converge, judged on their maps and on the
  tally of 19.7 (overlap trajectories, handover, iterations, wall time),
  not on a resolution number. The bar returns with the solvent prior.
- **The streptavidin regression (19.4) is closed.** There is none to
  attribute: the posterior build and the frontier build are both 10 of 10
  on streptavidin (19.7), and the cavgs-route failure that opened the item
  predates the posterior path. The bisect is dropped.
- **The stage exits of `solve3D_cavgs` were not built** when 19.8's tail
  state fix was committed under that title (a01ac3226 carries the tail
  fix, not the exits). Built today: the `lpexit` emission and the
  `solve3D_stage_behind_fsc05` omission on the class-average route, same
  eligibility as the particle route, res05 read from the temporary project
  (`solve3D_cavgs_policy.md` section 5).

### 21.2 Why known poses are not needed

Gate item (3) of section 11 asked for coverage on known poses because a
coverage test needs, per particle, a pose the scoring did not choose. A
simulation gives that trivially, and nothing else: the Gaussian model is
true on simulated data by construction, so a flat coverage histogram there
says the machinery computes the posterior it claims to (the temperature,
the kernel's bin sum, the polar sampling), and says nothing about whether
the data's posterior is that posterior. The question the temperature
raises is the second one. On real data the pose the scoring did not choose
exists too: the particle's **converged pose from the finished run, scored
against the stage's other-half reference**. The other half did not hold
the particle (section 12: the cross-half rows are reference-independent,
`null_a` zero), so that row is a prediction of where the particle is from
the other half's data alone, and the converged pose is where the whole run,
many bands later, put it. Whether the converged pose is the truth is the
one assumption, and it is weak where it matters: at a stage's band the
posterior is broad against the final angular precision (the stage-3 row's
99% set spans tens of degrees where the final pose is good to a degree or
two), so the proxy's error is a small fraction of the width being tested.
At the last stage it is not, and there the gold-standard bootstrap of the
final poses (the two halves' independent poses of the same particle) is
the proxy: the other half's pose, not the merged one.

### 21.3 The test

Per particle and stage k, with p_x the tempered cross-half posterior over
the stage's direction grid at beta and j* the grid direction nearest the
converged pose:

- the **posterior quantile** of the converged pose, `pit_x` = sum of
  p_x(j) over the directions j with p_x(j) > p_x(j*): the mass the
  posterior places above the truth. A calibrated posterior makes `pit_x`
  uniform on [0, 1] over particles (the probability-integral transform; for
  a 1 - alpha highest-density set, coverage is the fraction with `pit_x`
  below 1 - alpha). Too sharp, the histogram piles at 1: the truth falls
  outside the credible set more often than alpha. Too broad, it piles at 0.
- the **truth's depth**, `nats_x` = T (L_x(j*) - L_x,min) in nats: the
  same information as a distance rather than a quantile, comparable across
  stages.

The readout is beta* per stage: the beta at which the `pit_x` histogram is
flat (equivalently, coverage holds at alpha of 0.5, 0.1, 0.01 at once).
Three outcomes, each decisive for something already in the code:

- beta* near 1 at every stage: the temperature is the data's likelihood
  scale, section 3's claim; `prob_beta` 0.3 is exploration and nothing
  else, section 20's n_eff control is its right frame, and the annealing
  tail's endpoint of 1 is the physical posterior.
- beta* near 0.3 (or any constant below 1): the physical posterior is
  broader than the model says, by that factor (coloured noise, CTF error,
  the micelle); the tail then over-sharpens at its end and `prob_beta_end`
  should be beta*, and any weighted reconstruction (section 18's deferred
  item) must use beta*.
- beta* varying with the stage, or with the particle (the decided-by-both
  subset of section 9 against the rest): the mismatch is not one factor,
  and the n_eff control of 20.2 is the only policy that does not care.

The junk fraction sets a floor: a particle whose converged pose is wrong
sits outside every credible set at every beta, so the histogram's pile at 1
has a floor the test cannot remove, read off the decided-by-both subset
against the rest (the rest's pile is the floor).

### 21.4 Where it runs

Nothing new in the search. The diagnostic of section 5 already scores the
deterministic cross-half row (`_x`) and records the particle's stored
direction on the current grid (`iproj_p`). Run on a copy of the finished
project (the converged poses are then the stored poses, so `iproj_p` is j*)
with the stage's half maps as the references and the stage's `lp`, one
iteration with `prob_post_diag=yes`, it needs two more columns: `nats_x` at
`iproj_p`, a few lines beside `self_nats`, which already evaluates the row
at that direction, and `pit_x` at a fixed set of betas (0.1, 0.3, 0.5, 1;
the row itself is not written, and the posterior at a beta is one
renormalization of the distances the worker holds), so one pass per stage
gives the histogram at every beta. Simulation as the control of the machinery, cheap and worth having
once: the final map projected at the converged poses with the particles'
own CTFs and Gaussian noise at the measured per-shell sigma2
(`simulate_particles`), scored the same way; the histogram there must be
flat at beta 1, and if it is not the fault is in the temperature's
computation (section 3's kernel terms), not in the model. Real minus
simulated is the model mismatch.

This runs after a finished run and blocks nothing in 20.7; its answer sets
the tail's endpoint and decides whether section 18's weighted
reconstruction can use the posterior as it stands.
