# cls_expansion: method, benchmarks, and why it works where the diffusion-map split did not

Status: implemented on branch `flex-cls-split` (2026-09-30 to 2026-10-01). Run with

    simple_exec prg=cls_expansion projfile=... ncls=3 nparts=10 nthr=8

Since 2026-10-01 the flex model is the only mode of `cls_expansion`; the diffusion-map and kPCA split
were removed after the measurements in section 4 (git history before that date keeps them).

Unit suite (34 checks, about 5 s): `simple_test_exec prg=unit_heterogeneity suite=flex_cls_expansion`.

## 1. Method

Per 2D class, with every member brought into the class frame by its own in-plane alignment:

1. **Model.** `y_i = c_i (mu + U z_i + N nu_i) + noise` on the Fourier half-plane up to 30 A
   (`lp=` overrides): `c_i` the member's CTF, `mu` the class mean, `U` a rank-8 basis (`neigs=`),
   `N` four nuisance columns (the pose tangents dx, dy, dtheta of the class mean, and the mean
   itself for contrast), `z_i` the member's latent. Noise weights are the member's CTF-squared over
   the sigma2 shell spectrum (canonical sigma2 if present, else the class residual spectrum).
   Fitted by an ALS probe (3 random starts, best weighted residual) then PPCA EM.
2. **Cross-fitted embedding.** Members are dealt into 5 folds; the basis of the other folds embeds
   each fold, and each held-out latent is carried into the full fit's frame through the bases
   (the fold reconstruction projected onto the full basis). The in-sample posterior is
   overconfident (a column fitted on a member's own noise resolves that noise), the cross-fitted
   one is honest; its scatter over its posterior width per component is printed as the
   calibration (1.0 = pure noise).
3. **Placement.** The placement latent is `[z, log residual]`, the residual being the member's
   cross-fitted weighted fit residual (1 = noise level), detrended against the member's CTF power.
   Divisive bisection into exactly `ncls` leaves: candidates are k-means on the standardised
   latent and on each single coordinate (the residual competes only on its own axis), scored by
   Ashman's D along the cut's own centre axis; the leaf whose best cut separates most is cut next.
   Junk the basis cannot describe is far from everything along the residual axis and forms its
   own leaf. A tied-covariance GMM seeded by the leaves gives the labels.
4. **One greedy in-plane round of cluster2D.** After the split the commander runs one stock
   cluster2D iteration on the split project with the sub-averages as references and the classes
   fixed (`refine=inpl`, iteration 3 so the previous alignment is used and the shift search is on,
   `trs=3`, no joint continuous optimiser), then restores the sub-averages from the refined
   alignment (`SIMPLE_FLEXCLS_LABELS=project` path: subclasses are the project's classes, parents
   its clusters). The commander estimates the sigma2 state with calc_pspec first when the project
   has none that resolves (a project that never went through abinitio2D, or one copied out of its
   abinitio2D directory, which registers the state by a bare file name). A 10-degree rotation
   window around each member's previous rotation was tried and measured against the stock sweep
   on the same 10180 split: 0.877/0.623 against 0.870/0.631 in the two structural bands, so the
   window is not kept (the stock sweep rotates 9 % of the members by more than 90 degrees without
   any effect on reproducibility). The round's gain is in the shifts, not the rotations: median
   rotation change 2-3 degrees, median shift change 7 px, every parent's reproducibility up
   (0.76-0.93 -> 0.83-0.97 on the first 12). This replaced the module's own brute-force pass
   (0.76/0.43 in the two structural bands against 0.88/0.62 now).
5. **Weights and averages.** `w_is = exp(-1/2 (z_i - c_s)^T A_i (z_i - c_s))` with `A_i` the
   member's own posterior precision (a member whose latent is within its uncertainty of a centre
   is pooled, no bandwidth, no population target); labelled members always count fully. The
   delivered sub-averages are CTF-corrected weighted averages regularised like the class averager:
   the CTF-squared sum plus the inverse signal power per shell, the signal power from the
   subclass's own even/odd FRC (so a few-member subclass is damped where it has no signal rather
   than amplifying noise on its CTF-zero rings), written as `cls_expansion_cavgs.mrc` plus the even
   and odd member versions.
6. **Readouts delivered per subclass** (`os_cls2D`, class map): `pop`, `neff` (Kish effective
   population), `sep` (Ashman's D of the parent's first cut) and `repro`, the label-free
   cross-half reproducibility: the correlation, over the fit band and sigma2-weighted, of the
   even-member difference image between the subclass and its most different sibling with the
   odd-member difference image. A real structural difference reproduces (towards 1), a noise-driven
   split does not (towards 0).
7. **Distribution.** The standard cls_expansion master/worker parts, the class being the scheduling
   unit; part files are merged by the master. 46k particles at box 320 in 2 min on 80 cores,
   132k in 3 min.

## 2. Benchmarks against the diffusion-map cls_expansion

Synthetic (planted truth, abinitio2D classes of a mixed stack, SNR 0.05):

| test | flex | diffusion maps |
|---|---|---|
| EMD-8440/8445 50/50 mix, 60 classes, 2-way: pure subclasses / majority purity / NMI | 143 of 180 / 0.90 / 0.24 | 35 of 120 / 0.68 / 0.01 |
| five-state ladder 40/30/15/10/5, 2-way: purity / NMI | 0.62 / 0.18 | 0.42 / 0.01 |

Unit suite, one class each (flex only; the control has no equivalent): noiseless exactness,
two-state recovery 1.0 at rank 1, 2 and 4, 80/20 at 1.0; junk isolation at 15 % and 30 % with
recall 1.0 and a 100 % junk leaf; continuous motion ladder: latent tracks the motion at 0.95,
five rungs ordered with a span of 1.3 of 2.0, sub-averages at 0.90 of the parent's correlation to
the truth; null split reproducibility 0.015 vs 0.32-0.66 for the rungs; labels invariant to
defocus (0.51) and to a global sigma2 rescale; sub-pixel residual shifts recovered at 0.998 with
the pose tangents as nuisance.

Real data: labels from a 3D classification (flex_pca's states for 10180, the deposited classes
for 10076) are not ground truth. 10180: parents 0.615 majority purity, flex 0.648, diffusion maps
0.616. 10076 (converged abinitio2D of the full 132k, 200 classes; the earlier 10076 project had
zero shifts and was unconverged, see below): parents 0.357 (chance 0.347), flex 0.397 with 41
subclasses above 70 % one state, diffusion maps 0.363 with 4. The label-free comparison:

Best-pair cross-half reproducibility of the subclass differences per resolution band, 3-way,
median over parents:

| | > 60 A | 60-30 A | 30-15 A |
|---|---|---|---|
| 10180 flex | 0.93 | 0.88 | 0.62 |
| 10180 diffusion maps | 0.90 | 0.27 | 0.07 |
| 10076 flex | 0.86 | 0.81 | 0.79 |
| 10076 diffusion maps | 0.99 | 0.43 | 0.24 |

Seven real sets, 3-way, same settings, no per-dataset knobs (the 12 most populated parents, the
ones in the montages). Reproducibility as above; effect size = cross-half power of the best
sibling difference over the parent's cross-half power at 60-30 A (how big the difference is, where
reproducibility only says whether it is real); largest share = the biggest subclass's fraction of
the parent (a balanced split of conformers against a minority peeled off a dominant class):

| dataset | particles | 60-30 A | 30-15 A | effect size | largest share | wall (80 cores) |
|---|---|---|---|---|---|---|
| ribosome 10076 | 132k | 0.81 | 0.79 | 2.87 | 0.64 | 6 min |
| spliceosome 10180 | 46k | 0.87 | 0.63 | 1.06 | 0.45 | 4 min |
| "not" (unpublished) | 805k | 0.88 | 0.82 | 0.57 | 0.81 | 56 min |
| trpm4 (sieved) | 576k | 0.78 | 0.70 | 0.40 | 0.81 | 13 min |
| EMPIAR-13553 | 1.34M | 0.98 | 0.93 | 0.33 | 0.53 | 4 h 30 |
| 10028 ribosome ref | 105k | 0.68 | 0.56 | 0.33 | 0.47 | 41 min |
| flipqr (sieved) | 504k | 0.73 | 0.57 | 0.26 | 0.72 | 9 min |

Every set reproduces in the 30-15 A band. The effect size, not the reproducibility, ranks the sets
the way the montages read: the spliceosome shows conformers (balanced split, effect 1.06), the
ribosome ratchets (2.87), 13553 and "not" split reproducibly into differences a third of the
spliceosome's with one dominant subclass. Sibling centring contributes to the effect size (the
spliceosome parent with the largest value also has a 14 px sibling shift), and over all 300
parents a few tiny junk subclasses blow the ratio up (trpm4, flipqr reach 20-30), so a per-parent
readout needs a population floor. Scorers: `repro_bands.py`, `split_amp.py`, `split_amp2.py`.

Neither method's differences are in-plane pose (aligning the sub-averages changes nothing),
contrast (|corr(difference, parent)| 0.3-0.46), or defocus (between-subclass defocus spread over
the within spread 0.14-0.19 for both after the flex fix below).

Ablation, model vs restoration (10180, 3-way; a labels-file mode, since removed, ran external labels
one-hot through the flex restoration):

| labels / restoration | > 60 A | 60-30 A | 30-15 A |
|---|---|---|---|
| flex / flex | 0.90 | 0.74 | 0.36 |
| diffusion maps / flex | 0.48 | 0.45 | 0.38 |
| diffusion maps / make_cavgs | 0.90 | 0.27 | 0.07 |

The 60-30 A reproducibility is the model's (0.74 vs 0.45 under the same restoration); the 30-15 A
value is the restoration's (the same labels reach 0.38); the control's 0.90 above 60 A is its
restoration's background handling, not its labels (0.48 under ours).

Wall time (80 cores): 10180 46k, split 71 s + greedy round 29 s + restoration 36 s, about 4 min
in one call with the sigma2 estimate; the diffusion-map split took 2 min. Per 1000 particles the
split costs 0.4 (box 180) to 1.8 s (box 420), about twice a stock cluster2D pass, the restoration
about 1.5 times; five of the seven sets sit on that curve. The two off it: 10028 is slow in every
phase including the stock cluster2D round (9 s per 1000, its 1083 stacks under the EMPIAR tree),
so its data path, not the method; 13553 (one 550 GB stack) is normal in the cluster2D round and
10 times slow in the split and restoration, which fetch members class by class, a random walk
through the file where cluster2D streams particles in index order through its cache. Reading each
worker's particle range once in index order and serving classes from memory is the fix; not done.

Build note: the git commit hash now lives in `simple_gitinfo` (SimpleGitHash.h) instead of the
root header, so a commit recompiles one object; before, every commit rebuilt the tree.

Input caveat found on the way: the 10076 project used for all earlier numbers
(`10076_flex2D/a200`, a flex2D-branch abinitio2D snapshot at iteration 5 of an unconverged run)
had every particle shift at 0.0; its parents were rotation-only blobs at chance against the
labels and both splits scored at chance on it. The numbers above are from a converged
abinitio2D of the pristine import (`flexcls/10076/a2d_new`, 27 iterations).

## 3. Why it works

- **The metric is the likelihood's.** Members are compared in the CTF-corrected, noise-whitened
  Fourier band where the structural signal is, with the per-member CTF in the model rather than
  in the distance. Defocus, contrast and residual pose are modelled (nuisance columns, weights)
  instead of being left to dominate a distance.
- **The per-particle uncertainty is honest.** Cross-fitting removes the self-noise that makes an
  in-sample posterior overconfident; the posterior-precision kernel then pools members that the
  data cannot tell apart and keeps apart those it can, which is what the delivered weights are.
- **Junk has its own axis.** Junk is not a cluster in a low-rank latent (recall 0.52 with the
  latent alone); the fit residual separates it cleanly (junk 1.09-1.27 vs signal 0.91-1.08 at the
  noise level 1.0), and detrending it against CTF power keeps it from becoming a defocus axis.
- **Bimodality, not variance, picks the cut.** The largest-variance direction of a 2D class is a
  pose residual and the most bimodal one is junk vs signal; divisive bisection by Ashman's D lets
  the state surface at the second cut.
- **Every step was gated by a synthetic suite** that plants states, junk, a motion continuum,
  defocus-correlated noise and sub-pixel shifts. Three ideas that looked right were measured and
  removed because the suite or the real-data readouts said no: a full-band likelihood EM of the
  weights (collapses on a continuum; balanced, it only shuffles), same-view pooling of parents
  (+0.034 vs +0.031), and a Procrustes alignment of fold latents (smears rank-4 columns).

## 4. Why the diffusion-map split does not

- **Its distance is Euclidean between raw masked pixels.** At cryo-EM SNR the squared distance
  between two particle images is a constant noise term plus a tiny structural term; the
  nearest-neighbour graph is a noise graph, and diffusion coordinates on it are noise. On a
  planted, separable 50/50 mix it scores NMI 0.01, chance.
- **What it does split on is the very-low-resolution content**: overall density, background
  and blob size. Its subclass differences reproduce across halves above 60 A (0.90, 0.98) and
  collapse in the structural band (0.27, 0.34). That is why its montage rows look like one average
  at three brightness levels.
- **No CTF, no whitening, no nuisance, no uncertainty**: defocus, contrast and residual pose are
  not in the model, so there is nothing to stop them from entering the distance, and nothing to
  say which members the data cannot separate; its averages are plain unweighted means, so small
  noise-driven subclasses become empty tiles (down to 6 members on 10180, 4 on 10076).

## 5. What limits both on real data

The per-particle state signal within a 2D view class. The flex calibration on real classes is
1.1-1.8 times the noise per component (2.0-2.4 in-sample): enough to order members weakly and to
produce reproducible sub-average differences, not to label individual particles. Doubling the
members per view (same-view pooling) did not change it. Raising it means pooling members across
views, which is the 3D flex path.

---

## Appendix: the plan as written on 2026-09-30

Kept as written (the program was then called `cls_split`); the sections above record what was built and measured.

## flex_pca-style class splitting (`cls_expansion pca_mode=flex`): plan

Status: plan, 2026-09-30 (revised the same day: no even/odd machinery; fixed `ncls` subclasses and
fixed `neigs` rank per parent class; EM fits the basis, the subclasses are placed once afterwards).
Target branch: master. Nothing implemented yet.

### 0. What cls_expansion does today, and the answer to "does it weight?"

`cls_expansion` (`src/main/strategies/parallelization/simple_cls_expansion_strategy.f90`, commander in
`simple_commanders_denoise.f90`) runs per parent 2D class:

1. `transform_ptcls` (`simple_classaverager_restore.f90:780`): reads the class members, optionally
   phase-flips, shifts and rotates each particle into the class frame by KB gridding, returns
   real-space images and their plain average.
2. automask from the parent average, noise normalisation, soft mask.
3. `make_pcavecs` -> diffusion-map or kPCA embedding (`make_split_embedding`), ICM rank selection.
4. Euclidean distance matrix in the embedding -> k-medoids, k chosen by silhouette over
   `nsubcls_min..nsubcls_max` (or fixed `ncls`). Every class is split into at least 2.
5. Hard labels written to `ptcl2D/class` (global subclass) and `ptcl2D/cluster` (parent);
   `os_cls2D` rebuilt with `cluster`/`pop`.
6. A trailing distributed `make_cavgs` regenerates `cls_expansion_cavgs.mrc` from the hard labels.

**It does not build weighted class averages.** `cavger_update_sums` calls `accumulate_fplane`
without the optional `weight` argument (`simple_classaverager_restore.f90:342,345`), so every member
enters its subclass at unit weight. The embedding also has no noise or CTF model: phase-flipped,
masked real-space pixels with a Euclidean/kNN affinity. Sub-class membership is hard and
silhouette-driven, and a homogeneous class cannot answer "do not split".

### 1. Goal

A per-class **projection-aware low-rank covariance fit** (the flex_pca model, in 2D with the
in-plane pose fixed) that yields, per parent class:

- a latent embedding `z_i` with posterior precision, from a CTF- and sigma2-weighted likelihood;
- a fixed rank `neigs` (existing cls_expansion key); no rank gate;
- exactly `ncls` subclasses per parent (existing cls_expansion key): k-means centres on the final
  embedding, soft weights `w_is` from a tied-covariance GMM with those `ncls` components
  (`gmm_state_weights`, `nstates=ncls`); no auto placement, pruning or merging;
- **weighted sub-class averages**: `A_s = sum_i w_is C_i^* y_i / (sum_i w_is |C_i|^2 + eps)`; the 2D
  analogue of `reconstruct_flex_weighted_states`. No even/odd halves: the classes are small and the
  halves only served the rank gate and a per-average `res`;
- hard labels (argmax) kept for the existing project contract, and the weight table delivered
  through the `flex_weights` store so downstream tools can use the soft assignment.

No new user-facing keys: `ncls` and `neigs` are the two inputs, both already on cls_expansion. Poses are not refined
(the low-rank family is a poor aligner; flex2D branch, proposal section 15.1).

### 2. Model (per parent class, Cartesian half-plane, pose fixed)

For member `i` in the class frame (rotated by `-e3`, shifted by `-shift`, as `transform_ptcls`
does, but **without** phase flipping and kept in Fourier space):

```
y_i(h,k) = C_i(h,k) [ mu(h,k) + sum_q U_q(h,k) z_iq ] + n_i(h,k),   n_i ~ N(0, sigma2_i(|k|))
```

`C_i` is the particle's CTF evaluated on the rotated lattice: same defocus, astigmatism angle
rotated by `e3` (sign to be pinned by a one-particle test against explicit rotation). Weights
`W_i(h,k) = |C_i|^2 / sigma2_i(k)`.

E-step (flex_pca's MAP embedding, `_em_embed`): `A_i = I + U^H W_i U`, `b_i = U^H W_i (y_i/C_i - mu)`
(formally: `b_i = U^H C_i^* (y_i - C_i mu)/sigma2`), `m_i = A_i^-1 b_i`, `S_i = A_i^-1`.

M-step (flex2D branch, proposal section 6, transplanted from polar to Cartesian): for every
coefficient `(h,k)` a `(K+1)x(K+1)` weighted least-squares system

```
G(h,k) = sum_i W_i [1 m_i^T; m_i  S_i + m_i m_i^T],   rhs(h,k) = sum_i W_i (y_i/C_i) [1; m_i]
[mu; U] = G^-1 rhs        (Tikhonov ridge on the basis block, relative to G(0,0))
```

over all members, `K = neigs` fixed.

The fit is expectation maximization (PPCA EM, flex_pca's probe loop): `mu` starts from the plain
weighted average and `U` from random projections; prior-free (ALS) probe iterations first, then the
PPCA E-step/M-step alternation with rescaling (flex2D findings: random-projection init is 1000x too small; ALS is
scale-invariant so `U` is rescaled by the latent second moment). Iteration count: flex_pca's
`n_probe_iters` upper bound with the convergence stop.

Fit band: coarse. State differences are low-resolution features and the fit degrades as the band
marches (flex2D 15.11: labels held at 20-30 A beat 8 A). Default `lp_fit` derived, not typed:
the coarser of flex_pca's derived band rule and the parent's own `res` from `os_cls2D`; the
restoration uses the full band. Fit at `box_crop` matched to `lp_fit`, restore at `box`.

Sigma2: the project's canonical sigma2 state (what `make_cavgs ml_reg=yes` reads via
`cavger_read_euclid_sigma2`); if absent, run the same `ensure_canonical_sigma_state` path flex_pca
uses (move it out of `simple_commanders_flex_pca.f90` into a shared helper) rather than a silent
white fallback (memory: the white-noise fallback invalidated real-data runs).

### 3. States, weights, averages

- Placement, once, after the basis EM has converged: `kmeans_latent_targets` (`simple_flex_pca_targets.f90`)
  with `ncls` centres on the MAP embedding, then `gmm_state_weights` (`simple_flex_pca_gmm.f90`)
  with `nstates=ncls` and those centres for the responsibilities. The subclasses are not part of the
  EM (no mixture prior on `z`); `nsubcls_min`/`nsubcls_max` are ignored on this path and `ncls>=2`
  is required.
- Weighted averages: per state, numerator `sum_i w_is C_i^* y_i`, denominator `sum_i w_is |C_i|^2`
  (+ Wiener eps from sigma2 when ml_reg), on the full-box lattice, in the class frame (no
  re-gridding: the members were already rotated once). `res` inherits the parent's `os_cls2D` value.
  Members below the responsibility floor (`RESP_FLOOR`, as flex_pca) contribute nothing.
- Delivery: `cls_expansion_cavgs.mrc` in global subclass order, registered with `add_cavgs2os_out` like
  today; `os_cls2D` rows with `cluster`, `pop` (hard), `neff` (soft), `res` (parent's), `state=1`; `ptcl2D/class` = argmax subclass, `ptcl2D/cluster` = parent (unchanged contract).
  Weights: one `flex_weights_state_NNN.bin` per global subclass via `flex_weights_deliver`
  with `os_ptcl2D` (the store is keyed on box/smpd/layout digest, so a 2D field works; check that
  `flex_weights_consumable` does not assume a `vol_flex` sibling; if it does, add a `which_imgkind`
  argument rather than a second store).

### 4. Architecture and ownership (master)

The flow is the repo's: `ui -> exec -> commander -> strategy -> domain`. Nothing new is added
above the strategy; the new code is one domain module plus one extension of the class subsystem.

| layer | file | owns | must not |
|---|---|---|---|
| ui | `src/main/ui/simple/simple_ui_cluster2D.f90` (`new_cls_expansion`) | `pca_mode` gains `flex`; help text says `ncls` and `neigs` are required on that path | anything else |
| exec | `src/main/exec/simple_exec_denoise.f90` | unchanged: `cls_expansion` -> `commander_cls_expansion` | |
| commander | `src/main/commanders/simple/simple_commanders_denoise.f90` (`exec_cls_expansion`) | cmdline normalisation and defaults before `params%new`; the canonical sigma2 state on the master only (workers consume it); the four lifecycle calls; skipping the trailing `make_cavgs` when `pca_mode=flex` | reading part files, numerics |
| strategy | `src/main/strategies/parallelization/simple_cls_expansion_strategy.f90` | roles (shmem/master/worker), class partitions, the part-file protocol, the merge, all project writes (`ptcl2D`, `os_cls2D`, out-segment registrations of the cavg stack and the weight files). `split_one_parent_class` becomes a dispatcher: the existing body moves to `split_class_diffmap`, the new `split_class_flex` prepares the inputs (planes, CTFs, sigma2 rows via the class and sigma2 subsystems) and calls the domain | the model, the EM, the averaging |
| domain | `src/main/flex/simple_flex_cls_expansion.f90` (new) | `flex_cls_model` and its EM (`fit`), `embed`, `place_states` (k-means + fixed-K GMM through `simple_flex_pca_targets` / `simple_flex_pca_gmm`), `restore_states` (weighted averages). Inputs are arrays: complex planes, CTF planes, sigma2 rows, `ncls`, `neigs`. Outputs: `z`, `weights(ncls,N)`, `labels(N)`, `cavgs(ncls)` as Fourier images | importing `parameters`, `builder`, `sp_project`, `qsys_env`, or any part I/O (the 3D driver `simple_flex_pca_model` does; the 2D module stays pure so the unit test needs no project) |
| class | `src/main/class/simple_classaverager_restore.f90` (`transform_ptcls`) | the one gridding kernel that puts a member into the class frame; gains an optional Fourier-plane output (unflipped) and the rotated CTF parameters per member. The real-space path is unchanged | knowing about flex |
| sigma2 | `src/main/sigma2/simple_sigma2_bootstrap.f90` | `ensure_canonical_sigma_state` lifted here from the flex_pca commander. Today three private copies exist (flex_pca commander, rec3D strategy, cluster2D strategy `prepare_canonical_sigma_update`); cls_expansion calls the lifted one, the other copies can converge later, out of scope | |
| weights store | `src/main/flex/simple_flex_weights_state.f90` + `sp_project%add_flex_weights2os_out` | reused as is with `os_ptcl2D`; if `flex_weights_consumable` turns out to assume a `vol_flex` sibling, the fix is a `which_imgkind` argument in the flex module, not a project change | a second store |
| test | `production/tests/simple_test_flex_cls_expansion.f90` (new) | exercises the domain module alone on synthetic planes (auto-globbed, no CMake edit) | project fixtures |
| docs | `src/main/flex/README.md`, this note | | |

Reuse from the 3D flex code, exactly this and no more: `kmeans_latent_targets`
(`simple_flex_pca_targets.f90`), `gmm_state_weights` (`simple_flex_pca_gmm.f90`), the small
Hermitian solver `solve_real_spd_complex` (`simple_flex_reconstructor_latent_ops.f90`), and the
`flex_weights` store. Not reused: the EM driver, E-step, M-step and their submodules (all bound to
volume projection and the polar bank), `simple_flex_pca_rec3D`, `simple_flex_pca_rounds` and the
flex_pca strategy (its stage/round scheme and producer-owned part I/O are what this design avoids).
The per-coefficient EM in 2D is ~300 lines of new code; the flex2D branch's
`simple_flex2D_model.f90` is the reference for it and is not imported.

Reused unchanged elsewhere: `discrete_read_imgbatch`/`prepimgbatch`, the KB gridding loop,
`ctf%eval`, `load_sigma2_groups`, `add_cavgs2os_out`, the cls_expansion master/worker scheduling and
merge.

### 5. Parts: the standard SIMPLE distribution, not flex_pca's rounds

The existing cls_expansion strategy already follows the usual pattern and is kept as is:

1. The master (`nparts>1`, no `part=`) partitions the parent classes into `nparts` lists
   balanced by population (`prepare_class_partitions`), writes `cls_expansion_classes_partNN.txt`,
   and submits `simple_exec prg=cls_expansion ... part=NN class_assignment=<list>` through `qsys_env`.
2. Every worker runs the same shared-memory code (`run_local_split`) on its classes, in the
   master's directory, and writes only part files: today `cls_expansion_class_map_partNN.txt` and
   `cls_expansion_assignments_partNN.txt`; the flex path adds `cls_expansion_weights_partNN.bin`
   (pind, parent, local sub, `ncls` weights) and `cls_expansion_cavgs_partNN.mrc` (local sub order).
   A worker never touches the project.
3. The master waits on the job-finished flags, merges the part files once
   (`merge_worker_outputs`: global subclass numbering, then labels, weights and the concatenated
   cavg stack in that order), writes the project, registers the stack and the weight files, and
   cleans up. One submission, one merge, no rounds, no stages, no master-side numerics.
4. `nparts=1` runs the same `run_local_split` in process and skips the part files.

The strategy owns every part file, read and write, as it does today. Cost per class per EM
iteration is `N * ncoeff * (K+1)^2` (N = 1000, box_crop 128, K = 8: ~7e8 flops), negligible
against the particle read, so a class is never split across workers.

### 6. Validation (batch2, never the Mac)

1. Unit test (item 7), no data.
2. One-particle convention test: rotated-CTF sign, plane conventions against `transform_ptcls`'
   real-space output (inverse FFT of the returned plane must equal the returned image).
3. 10028 8000-particle reference run (`/mnt/beegfs/elmlund/afan/10028/flex2D_tests/ref.simple`):
   `ncls=2 neigs=4` on every class; do the two sub-averages differ visibly, and how do their
   weights split (a 50/50 split on a homogeneous class is the null result to expect).
4. Synthetic 50/50 mix (EMD-8440 vs 8445, `t4_mix5050`): abinitio2D's 60 classes are all ~50/50;
   score per-subclass purity after `cls_expansion pca_mode=flex` against the existing
   `pca_mode=diffusion_maps` on the same run. This is the decision test.
5. EMPIAR-10076 with the published per-particle labels (`labels/real_recall.py`): purity/NMI per
   subclass vs the diffusion-map path.

### 7. Risks

- **View vs state.** Within a parent class the residual out-of-plane spread can be the first
  dominant component, and a fixed `ncls` split will then partition on view. The current cls_expansion
  has the same confound. Mitigation: report the component's dominant resolution band and leave the
  split, note the confound in the manifest.
- **Small classes.** With `ncls` and `neigs` fixed, a 40-member class is fitted and split whatever
  its content; classes below `3 x (neigs + 1)` members are skipped (as cls_expansion skips pop <= 2).
- **Rotated-CTF sign and the gridding conventions**: covered by validation step 2 before anything
  else is run.
- **Weight store consumers** were written for `vol_flex` siblings; verify item 3's delivery does not
  break `flex_weights_consumable` for 3D runs.

### 8. Open decisions

- `cls_expansion pca_mode=flex` (recommended: reuses scheduling, merge, and the project contract) versus a
  separate `flex_cls_expansion` program (cleaner ownership, duplicated master/worker plumbing).
- Whether the weighted averages should also be re-centred per subclass (cluster2D centres its
  averages every iteration; cls_expansion does not).
- Whether to expose the soft weights to `map2ptcls`/selection downstream, or leave the hard labels
  as the only consumed output for now.

