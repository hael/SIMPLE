# PfCRT audit: what changed between a6ae317ae and HEAD (2026-09-02)

Reference: the best PfCRT reconstruction so far came from `a6ae317ae`
(2026-07-01, gridding backend, `filt_mode=nonuniform` post-hoc filtering).
HEAD (`c9d12d486` + the 2026-09-02 review fixes) converges on PfCRT again
but does not beat that map. This note lists the behavioural differences that
plausibly cost resolution, ranked, with the commits that introduced them and
the cheapest A/B for each. Everything below was verified by reading both
revisions; refactor-only changes are omitted.

## 1. Probabilistic sampling became near-uniform over the top-K candidates

Commit `843c2966c` (2026-07-08, one week after the reference run) replaced
the bounded sampler with a likelihood softmax and made it the only 3D path.

- Old: `sample_bounded_dist` returned the best of the top-K with probability
  about 0.9 and the K-th otherwise (`a6ae317ae:src/main/simple_eul_prob_tab_utils.f90:214-258`).
- New: `sample_likelihood_dist` draws from `exp(-(d - d_min))` over the top-K
  (`src/main/simple_eul_prob_tab_utils.f90:474-530`), used for the in-plane
  angle (`simple_eul_prob_tab.f90:345`), the reference frontier (`:560`), the
  per-state frontier and, via `prob_neigh`, the coarse subspace representative
  (`simple_eul_prob_tab_neigh.f90:485-490, 601-606`).
- The distance is not an absolute negative log-likelihood. `euclid_dist_from_crvec`
  returns `1 + crvec/norm` with `norm = wsqsum_ptcl_now * 2*nrots`
  (`src/main/pftc/simple_polarft_corr.f90:480, 518-529`), where `wsqsum` is
  the whitened particle power `sum |P_k|^2 / sigma2_k`. So `d` is a relative
  residual and the gaps between the best candidates are of order
  `delta(NLL) / (whitened power)`, i.e. much smaller than 1 for cryo-EM
  particles. `exp(-delta d)` is then close to 1 for every candidate and the
  draw is essentially uniform over `K = athres * n / 180` candidates
  (`likelihood_nsample`). The design note on the commit ("distances are
  already noise normalized") is true of the whitening but not of the
  division by the particle power.
- The same commit removed the per-particle min/max normalisation of the
  assignment frontier score (`simple_eul_prob_tab.f90:618-629` vs
  `a6ae317ae:585-595`); the frontier now ranks raw absolute loss, so
  high-contrast particles are always served first. `prepare_ref_score_vectors`
  still computes `score_min/score_spread` that nothing reads.
- Affects both backends and every abinitio3D stage (all stages use
  `refine=prob|prob_neigh`). No CLI knob selects the old sampler.

Suggested fix: scale the exponent by the whitened particle power
(`exp(-(d - d_min) * wsqsum_ptcl_now)` or an equivalent absolute NLL gap)
so the softmax has its physical temperature, or restore the bounded sampler
for an A/B.

## 2. The NU evidence cutoff now caps the matching band in plain `nonuniform` mode and latches

Commit `2b4e5b0b6` (2026-08-05) changed the gate in
`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90:144` from
`if( .not. (l_nu_refine .or. l_nonuniform_lpset) ) return` to
`if( .not. l_nonuniform ) return`.

- abinitio3D stages 6-8 run `filt_mode=nonuniform` with `nu_refine=no`. At
  `a6ae317ae` those stages set `kfromto(2)` from the FSC at `lplim_crit=0.143`
  (plus `incrreslim`); at HEAD they take the project `lp` written by the
  assembly (gridding: finest selected NU label; PCG: `pcg_nu_matching_lowpass`,
  which with `nu_refine=no` is the raw finest supported evidence cutoff of a
  static eight-candidate bank).
- `set_bp_range3D` writes the result back to the project and
  `adopt_reprojection_model_range` re-stamps it, so a cutoff that is coarser
  than the FSC limit never recovers within the run.
- On PCG the evidence state can additionally ride frozen for up to
  `NU_EVIDENCE_REBUILD_MAX_LAG = 5` iterations
  (`simple_rec3D_pcg_strategy.f90:100, 467-499`); the binding-band guard only
  fires at the cap.

A/B: one-line revert of the gate (require `l_nu_refine .or. l_nonuniform_lpset`
again) on an otherwise identical run. This applies to the gridding control
run too.

## 3. PCG references carry no Wiener shrinkage, no band limit, and a hard envelope

What the matcher saw at `a6ae317ae`: every reference voxel was either a hard
low-pass of the raw half at its NU-selected cutoff or, in the auxiliary
replacement label (`use_static_nu_aux_replacement = l_ml_reg .and. .not. l_nu_refine`,
true on the abinitio3D default), the FSC/SSNR-Wiener-shrunk ML half map
(`add_invtausq2rho`, `simple_reconstructor.f90`).

What it sees at HEAD on the PCG backend:

- `set_ml_prior` (P_tau, the FSC/SSNR shell diagonal) and `set_nu_prior`
  (Q_NU) are mode-exclusive (`simple_reconstructor_pcg.f90:617`), so with
  Q_NU active there is no shell-wise SSNR term anywhere in the PCG pipeline,
  and the matcher applies no filter at all to Q_NU references
  (`simple_matcher_refvol_utils.f90:339-342`).
- Q_NU is calibrated to about 15 % amplitude suppression of prior energy
  (`simple_parameters_phases.f90:925`), orders of magnitude gentler than a
  Wiener filter that drives unsupported shells toward zero.
- The regularized replay always solves under the conservative density
  envelope and returns `x = P u`, so everything outside
  `automask3D(envmsklp=20 A)` is exactly zero in the shipped maps. The
  matcher-side rule that references must never hard-remove density present
  in the images is enforced at the masking step but the identical operation
  now happens inside the estimator.
- The base solve is two CG iterations (`maxits_pcg=2`, `rtol=0`) on an
  approximated Toeplitz operator; the regularized replay skips the
  closed-form FSC-shrinkage initial guess in NU mode
  (`simple_rec3D_pcg_strategy.f90:2627-2628`) and warm-starts from the
  previous ML half, so without a usable previous half the shipped map is
  essentially the unregularized solution.
- `RHO_FLOOR_FRAC = 1e-2` of the shell mean damps the least-sampled modes;
  gridding's exact `cmat/rho` has no counterpart.
- The auto-lambda / auto-target loop is closed on the shipped-pair FSC
  crossing, which the code itself labels "never a resolution claim"; nothing
  in the loop measures how well the reference explains the images.

A/B: `rec_backend=gridding` at HEAD isolates the backend. There is no
`pcg_nu_lambda_rel=0` route with NU any more (`validate_nu_replay_request`
rejects it), so a P_tau-plus-post-hoc-NU control on PCG needs code.

## 4. Stage policy changes in abinitio3D

- Stage 1 `prob_neigh_mode` `snhc` -> `shc` and `NSPACE(1)` 1000 -> 500
  (`33c9d3663`, 2026-08-10, "100 % success rate on PfCRT"). This improved
  robustness; not a regression candidate, listed for completeness.
- Hard 20 A stage-1 low-pass floor (`221b8388c`, 2026-08-26,
  `simple_estimate_ssnr.f90:23, 305-338`), not overridable by `lpstart`.
- The rule "stage 1 borrows stage 2's crop" was deleted in `33c9d3663`; the
  `box_crop`/`smpd_crop`/`trslim` ladder for stages 2-7 is now anchored on
  the coarser stage-1 box. Diff the printed
  `lpset lp box_crop smpd_crop trslim` ladders of the two runs before
  assuming this is neutral (the `minbox=88` clamp may hide it).
- `rec_backend=pcg` switches backend mid-run at stage 3 (`96e372e92`,
  2026-08-16); Q_NU is active only in stages 6-8 and the final rec.
- The final map is Q_NU-regularized in-solve instead of classically
  postprocessed.
- Stage-boundary and symmetric reconstructions are now built from
  `cline_reconstruct3D` with an explicit allow-list plus `trail_seed=yes`
  (`simple_abinitio_utils.f90:253-291, 656-665`).
- `objfun_den`/`objfun_den_w` are dropped for the NU stages 6-8 at HEAD
  (`simple_abinitio_controller.f90:566-572`); the July "10/10" note used
  `objfun_den=yes objfun_den_w=0.3`.

## 5. Joint continuous in-plane refinement on by default

`inpl_cont=yes` (`b0fb84470`, 2026-08-14). A joint (sx, sy, theta) L-BFGS
polish runs on every particle after selection, writes fractional Euler
angles, and on no-improvement re-scores the seed with the joint objective.
Absent at `a6ae317ae`. A/B: `inpl_cont=no`.

## 6. Amplitude convention and the gridding control

- `9b82e8760` (2026-08-22) retired the `/box` (reconstructor) and `*box`
  (`expand_cmat`) pair. Designed net-neutral, but `autorescale_old_convention`
  (`simple_matcher_refvol_utils.f90:394-406`) multiplies a reference by
  `params%box` when its foreground sigma is below `10/box`, evaluated on a
  `box_crop` volume; a genuinely low-contrast map could trip it. Watch the
  log for `AUTO-RESCALED OLD-CONVENTION`.
- Gridding at HEAD is not gridding at `a6ae317ae`: shipped halves are now
  deapodized (`118388291`), the FSC spherical mask radius changed from
  `box/2 - COSMSKHALFWIDTH - 1` to `msk_crop` (tighter mask, higher FSC,
  weaker `add_invtausq2rho` shrinkage), and `envfsc=yes` now adds phase
  randomization. A `rec_backend=gridding` control therefore bounds the PCG
  effect but does not reproduce the July map.

## Verified unchanged

3D convergence thresholds; `lplim_crit=0.143`, `incrreslim`, `prob_athres=10`,
`ml_reg` per stage; particle image preparation (`prepimg4align`); reference
masking (spherical soft mask only, no envelope, no B-factor, both versions);
`update_frac`/`nsample`/`fillin`/`trail_rec` gating; `MAXITS` per stage;
symmetry-stage policy; `add_invtausq2rho` itself (byte-identical).

## Recommended order of experiments

1. `inpl_cont=no` on the current PCG route (free).
2. `rec_backend=gridding` at HEAD (isolates section 3).
3. One-line gate revert from section 2 on the better of 1-2.
4. Sampler temperature fix from section 1 (code change; largest expected
   effect, both backends).
5. Compare the printed crop ladders (section 4) against the July log.
