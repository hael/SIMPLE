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

1. **Already replaced by T_eff under the key** (delete once the default
   flips to `posterior`): the in-plane angular window (`prob_athres` and
   `dist_inpl` through `inpl_likelihood_nsample`, `calc_athres` for the
   in-plane, `likelihood_nsample`), the bounded sampler of `refine=neigh`
   (`gen_prob_objfun_val`, `angle_sampling`, `smpl_inpl_athres`), the power
   and squared samplers and their constants (`power_sampling`,
   `squared_sampling`, `EXTR_POWER`, `POST_EXTR_POWER`, `snhc_smpl_ninpl`,
   `smpl_ninpl`), the likelihood samplers of the tables (`gen_likelihood_val`
   and its bank form, `sample_likelihood_dist(s)`, `sample_likelihood_index`,
   `sample_bounded_dists`, `sample_power_dists`), the balanced frontier of
   the three tables and its draw width (`nhood_sz`, `NHOOD_FRAC`,
   `npeaks_inpl` as a width, the frontier's `prob_athres`), the reference
   check for flat rows that only the frontier needs.
2. **Refinement counts, to be driven by the posterior mass**: `projs_ns`
   from `calc_num2sample(prob_athres, dist)` in the dense 3D fill,
   `nhood_sz`, `snhc_smpl_ncls` and `smpl_ncls` in 2D, `npeaks_inpl` in 3D.
   The rule: shift-refine the entries holding a fixed fraction of the
   row's posterior mass at T_eff (the diagnostic's k99: median 3-5 at a
   stage's start, 1-2 as it settles), so the budget follows the data's
   width and shrinks as the posterior sharpens. This retires
   `prob_athres` from the fills.
3. **Neighbourhood construction**, the cost control that decides what is
   scored before anything is scored: the stochastic bound
   (`extremal_decay`, `neigh_frac`, `snhc_nrefs_bound`, `MIN_STOCH_FINE_REFS`),
   the subspace cells and peak counts of `prob_neigh_mode=state|geom`
   (`npeaks`, the pooled capacity), the 2D `snhc_nrefs_bound`. The honest
   replacement is the previous iteration's posterior mass deciding which
   cells to re-score; a larger change, last.
4. **Then the keys themselves**: `prob_assign` (posterior the only path),
   `prob_inpl`, and the refine modes that only existed to pick a sampler
   (`snhc_smpl`, `snhc_smpl_many`, `inpl_smpl`, `greedy_smpl`, and
   `refine=neigh`'s sampled in-plane), collapsed onto the posterior path.
   What stays is what is not a sampler: the resolution ladder, the shift
   bounds (`trs`), the continuous polish and its bounds, `nspace`.

Not before: the calibration of `prob_beta` (section 17), the schedule if
a constant is not enough, the PfCRT 10/10 at the bar, the 2D regression,
and the known-pose calibration of the temperature (gate item 3), which
is what makes the one number mean what it says.

