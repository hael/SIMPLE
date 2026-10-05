# Comment novels: replacement drafts

Companion to `comment_novel_inventory_2026-09-30.md`. It gives one draft for each audited block
(the 78 blocks of 150 words or more, plus 7 procedure-level records). Line numbers refer to
commit 714a15a99.

The audit checked each draft against the code. No maintainer has reviewed them yet. The drafts
keep the local comment style. Where a file has a `!@descr:` line, it stays as line 1; "(keep)"
means the existing line is unchanged.

Each entry ends with a **Move to** line, which says where the rest of the old comment should go:

- **drop**: the content is in git history, or nowhere worth keeping.
- **covered**: an existing doc already says it, so nothing needs moving.
- **new**: a doc that doesn't exist yet and would have to be written.

## flex_pca EM (`src/main/flex/simple_flex_pca_em*.f90`)

**`simple_flex_pca_em.f90:129-144`** (OK)
```fortran
!> Per-fit EM state, owned by the driver so one loop can advance two resident fits.
!! Lifecycle: mix_* and their work arrays (rhs0th, mkth, lwth, rkth, mxa_*) live across
!! iterations; only kill_probe_fit or the rank-change resize in fit_iter_begin may free them.
```
**Move to:** drop the crash story. The hoist rationale is covered in finding 5 of
`flex_pca_architecture_audit_and_refactoring_plan_2026_09_17.md`.

**`simple_flex_pca_em_compose.f90:1-20`** (STALE)
```fortran
!@descr: flex_pca EM: multi-band basis composition from finished runs (SIMPLE_COV_COMPOSE)
!! SIMPLE_COV_COMPOSE=<dir>[,<dir>...]: finished runs; per run polished > merged > plain namespace.
!! Columns are Fourier-padded to box_crop (zero beyond their band), Gram-Schmidt'ed coarse box first;
!! residuals below COMPOSE_R2_FLOOR drop, prior variances follow the rescaling. compose_cut_reembed
!! then cuts the union to its signal subspace (SIMPLE_COV_COMPOSE_CUT=0 skips it).
```
**Move to:**
- The "why not march one basis" text goes in a new compose section of
  `doc/algorithms/heterogeneity_analysis/flex_pca.md`.
- The rejected variants of 2026-09-10 go in a flex_pca decision log (new).

**`simple_flex_pca_em_crossfsc.f90:118-131`** (STALE)
```fortran
!> Per-component, per-shell invtau2 from one paired record in this fit's indices: arm 1 = matched
!! cross-fit FSC (unmatched -> kill branch), arm 2 = min(internal, cross) (unmatched -> internal).
!! Always this fit's own H; components beyond the record's rank get no ridge.
```
**Move to:** a paired-engine note (new), or drop if arm 2 is retired. Arm 2 is hard-wired off.

**`simple_flex_pca_em_crossfsc.f90:209-222`** (STALE)
```fortran
!> Append one paired=1 crossfsc record per iteration, after both fits' tails: per-fit internal
!! FSC, Gamma and own e+o H, plus greedy signed |cos| matching on the delivered bases (prev_real)
!! and each matched pair's cross-fit FSC. Greedy, not varimax+Hungarian.
```
**Move to:** drop the `~/ribo_local` path.

**`simple_flex_pca_em_estep.f90:933-947`** (WRONG)

Delete the block together with the uncalled `box_*` routines. The live rationale is in the
"index-list packing" block (1008-1017). That block's "doubled PCG lattice" sentence is also stale:
the kernels ship packed (1182-1206). Fix the "band-boxed" note at `simple_flex_pca_rounds.f90:94`
as well.

**`simple_flex_pca_em_fit.f90:34-53`, and the header at :12** (WRONG)
```fortran
!> Initial EM basis: delegates to init_basis_datafree (master/shared-memory only).
```
Drop the in-body block. The routine looks unreachable: every non-compose master run goes to
`run_flex_pca_paired` (`simple_flex_pca_model.f90:228`), and workers never call it. Whether to
remove the routine is the maintainer's call.

**`simple_flex_pca_em_fit.f90:126-140`** (STALE)
```fortran
!> Data-free EM start: lowest-|k| band lattice points (col_sep apart) realised as masked cos/sin
!! pairs and orthonormalised; deterministic. One capped data pass calibrates sig2 and Gamma^0;
!! the master also writes the it000 copies the paired merge deflates against.
```
**Move to:** covered in `doc/algorithms/heterogeneity_analysis/flex_pca.md:43-51`.

**`simple_flex_pca_em_iter.f90:132-154`** (STALE)
```fortran
        ! Optional probe-stage cap (SIMPLE_COV_PROBE_MAX, default off), applied per halfset; embedding
        ! still uses every particle. The cap is a cross-process total: only a worker divides it by nparts.
```
**Move to:** drop. `cov_stage_subsample` (`simple_flex_pca_em_env.f90:81-85`) already carries the
rationale.

The same procedure has orphan blocks at 48-88, 106-107 and 317-320 with no code behind them:
- a "default 0.97" threshold, where the code uses 0.999999;
- a hybrid split of "0.55*band", where the code uses 0.72;
- line 630 says `SIMPLE_COV_PAIRED=1`, but the paired engine runs unconditionally.

**`simple_flex_pca_em_iter.f90:494-506`, and `fit_stage_config` as a whole** (WRONG)
```fortran
!> Per-fit stage defaults, all hard-wired: polar E-step, MCFA prior (K=COV_EM_MIX), mean-shaped
!! deflation, fixed per-particle contrast (no ECM, no a-scaled M-step), band and convergence settings.
```
In the body, use one line per setting, e.g.
`! MCFA: one basis, K latent Gaussians (Baek et al. TPAMI 2010); K=1 is plain PPCA EM.`

**Move to:** the 10076 measurements go in the flex_pca decision log (new).

**`simple_flex_pca_em_pairmerge.f90:1-34`** (STALE)
```fortran
!@descr: flex_pca paired engine final stage: merge the two fits' statistics, don't refit
!! Each fit stashes its last-iteration raw M-step statistics and entry-frame basis. B is rotated into
!! A's frame (R = polar factor of the entry cross-Gram; Y.R, R^T rho R), the four quarter-sets are
!! summed, the cross-fit-FSC ridge is added with the summed H, and one joint coupled solve runs;
!! deflation, orthonormalisation and a gauge fix to A's frame follow.
```
**Move to:**
- The derivation goes in the paired-engine note (new).
- The measured percentages go in the decision log.
- The mode-2 text goes, together with the dead branch.

**`simple_flex_pca_em_pairmerge.f90:108-123`** (WRONG)
```fortran
!> Cross-half match after projecting the shared it000 init out of both bases (fits start identical).
!! Optionally pa_cos (principal cosines of the deflated spans) and sub_cos (each A component's
!! cosine with the span where pa_cos >= thr), which sets the axis weights. ok=.false. without stamps.
```
Move lines 108-109 to `probe_paired_merge` at :247, which has no header. The same stale "rank
gate" wording appears at 357-358 and 394-396.

**`simple_flex_pca_em_polar.f90:21-35`** (WRONG)
```fortran
!> Bank directions: ~40 particles per direction, clamped to [1000,4000] and even (build_refspiral).
!! The bank holds all ndir directions, so memory scales with ndir.
```
**Move to:** drop. The measurements come from the removed moment estimator.

**`simple_flex_pca_em_polar.f90:142-154`** (OK)
```fortran
!> project_fplane(apply_ctf_amp=.true.) for the mean, cmplx_plane only: identical interpolation, but
!! the plane is zeroed only at (re)allocation and no ctfsq/transfer copies are made. Assumes a fixed
!! disc (frlims/nyq) per process, as subtract_mean_banded does.
```

**`simple_flex_pca_em_solve.f90:191-204`** (WRONG: placement)
```fortran
!> One particle's MAP solve; nml>0 adds ECM contrast updates against the current basis,
!!   a <- (m'y + b'z) / (||m||^2 + 2c'z + z'Gz + tr(G A^-1)),  clamped to [0.1, 5].
!! tr(G A^-1) is the posterior variance; dropping it biases a high.
```
Above `spd_inv_dp` (:479), add: `!> SPD inverse by Cholesky, same rescaling and ridge escalation
as spd_solve_dp; zeros if all attempts fail.`

**`simple_flex_pca_em_mstep.f90` `fit_iter_finish` (28-485)** (WRONG)
```fortran
!> Master tail of one EM iteration for one fit: Gamma and likelihood, optional cross-FSC ridge,
!! coupled per-half solve, FSC-Wiener merge, mean-shaped deflation, re-orthonormalisation and
!! basis swap. Frees per-iteration fields only; mix_* persist (see probe_fit_t).
```
Fixes inside the procedure:
- **91-94:** `SIMPLE_COV_XFSC_REG` has no reader, and the routine is `xfsc_prep_iter`.
- **127-132 and 159-162:** the Wiener filter is hard-wired on, so the off branch is dead.
- **259-264:** `SIMPLE_COV_EM_DEFLATE` has no reader.
- **279-283:** the loop runs once.
- **392-394:** the even/odd stop is disabled (:470).
- **428-429:** "z is zero-initialised" is false (see inventory H.1).

## flex_pca: other files

**`simple_flex_pca_deconv.f90:1-17`** (OK)
```fortran
!@descr: flex_pca latent deconvolution: calibrated per-particle noise + an empirical-Bayes mixture prior fitted through it
!! MAP latents z_i = A_i^-1(D_i x_i + e_i), D_i = A_i - P: E[z|x] = R_i x, R_i = A_i^-1 D_i, noise a*A_i^-1 D_i A_i^-1,
!! with the scalar a calibrated from even/odd half solutions. The prior is a K-Gaussian mixture fitted by extreme
!! deconvolution (Bovy, Hogg & Roweis 2011), K by particle-half held-out log-likelihood; z/precision become
!! posterior means/precisions under it.
```
**Move to:** a "Latent deconvolution" paragraph in `heterogeneity_analysis/flex_pca.md` (new).

**`simple_flex_pca_merge.f90:1-22`** (STALE)
```fortran
!@descr: Two-gate agglomerative merge of over-provisioned flex_pca states.
!        Gate 1 (orientation): a state whose viewing-axis distribution stands out from its peers is a view
!        cluster, folded into a sufficiently similar map. Gate 2 (volume): pairs whose deviation maps agree within
!        their own half-map reproducibility fuse under complete linkage. Latent distance never merges.
!        On when SIMPLE_COV_MERGE is non-zero or preimage_auto=yes; SIMPLE_COV_MERGE=0 always wins.
```

**`simple_flex_pca_pcg.f90:1-19`** (WRONG)
```fortran
!@descr: flex_pca coupled M-step on the PCG operator (rec_backend=pcg); design: flex_pca_envelope_support.md 3.4.
!  CG on the native-lattice basis u: b = S^H y and T = S^H S from 2x-lattice KB deposits (scale OSMPL_PAD_FAC**3),
!  Nyquist-ball band limit, hard support P, floored per-voxel coupled divide as preconditioner. maxits<=0 ships
!  solve_coupled_basis_exp; otherwise CG warm-starts from the masked, LS-scaled gridding solution. put_back writes E*u.
```
**Move to:** an "as implemented" paragraph in §3.4 of `flex_pca_envelope_support.md`. That note's
Status line still reads "Planning only".

**`simple_flex_pca_pcg.f90` `finalize` / `solve`** (WRONG)
```fortran
! finalize: One thread per pair writes a contiguous slab; a blocked transpose builds pair-leading khat. Not counted in solve()'s seconds.
!> M-step solve on one half. maxits<=0 ships the gridding solution (solve_coupled_basis_exp); otherwise CG on the 2x
!! right-hand sides starts from it, masked to the support (cg_core rescales it). Returns E*u on the expanded lattice.
! wthreads: Y(q) transforms run one at a time below: use unthreaded FFTW plans (threaded ones stalled in os_sem_down).
```
Further fixes:
- Delete lines 1883-1884 and 1898-1906, which describe a cold start and an opt-in switch. The code
  hard-wires `iwarm=1` (:1907).
- Line 1848 cites `simple_flex_pca_em_fit.f90:641`; that file has 529 lines.

**`simple_flex_pca_plane_cache.f90:1-19`** (WRONG)
```fortran
!@descr: flex_pca plane cache: the full-box prep's padded transform, restricted to the box_crop grid, kept on disk per particle
!! cache=yes and box_crop<box only. Per project row, the |h|<=box_croppd/2 block of the full-box padded transform
!! (== a box_croppd cmat; self-checked on particle 1), injected by plane_cache_fill ahead of gen_fplane4rec.
!! Direct access: record 1 header, record p+1 row p. Built by the non-worker process; adopted on matching
!! magic/version/rows/box/box_crop/smpd/size (the selection is not checked).
```

**`simple_flex_pca_planes.f90:1-17`** (WRONG)
```fortran
!@descr: flex_pca resident planes: prepped particle Fourier planes kept in memory across E-step passes
!! Shared-memory runs with the plane cache in use only. planes_batch_load serves a fully held batch by copy
!! (passes modify planes in place); otherwise it reads+preps and stores rows until 0.25*MemAvailable is used.
```
Also fix the log message at :62, which names `SIMPLE_COV_RESIDENT_GB`.

**`simple_flex_pca_polar.f90:1-41`** (WRONG)
```fortran
!@descr: polar-Fourier shared-direction basis bank for flex_pca
! Mean + basis sections are projected once per shared direction (cov_polar_ndir) on polar rings; particles are
! sampled at their continuous relative in-plane angle and ring-mean |T|^2 factorises G_qr = sum_k w_i(k) C_qr(k).
! Approximations: direction snap and radial |T|^2 (tazim). Quadrature measure, KB weights and CTF adjoint must
! match the Cartesian cov_herm_inner path.
```

**`simple_flex_pca_tester.f90:1-11`** (OK)

Keep `!@descr:`, then add
`! Fixed seeds. The fast suite deconvolves 4000 particles; the library suite 20000 at realistic noise.`

**`simple_flex_reconstructor_latent_ops.f90:799-816`** (STALE)
```fortran
    !> Flex analog of reconstructor::add_invtausq2rho: adds the per-component, per-shell invtau2(q,sh)
    !! (crossfsc_to_invtau2) to the diagonal rows pair_index(q,q) of the packed coupled density.
    !! Call after any distributed reduction and before solve_coupled_basis_exp.
```

**`simple_flex_weights_state.f90:1-18`** (OK)
```fortran
!@descr: flex per-state weight files: identity, science validation, transactions, delivery and loading
!! Builder-free layer over simple_flex_weights_file; row identity = canonical sigma2 layout digest of the same field.
!! Producer: flex_weights_deliver (validate the set, then publish). Consumers: flex_weights_consumable +
!! flex_weights_load_state / flex_weights_load_all. The range merge has no producer yet.
```

**`fileio/simple_flex_weights_file.f90:1-17`** (OK)
```fortran
!@descr: versioned binary persistence and transaction primitives for the flex per-state weight files
!! One committed flex_weights_state_NNN.bin per state: its weight over the full particle layout (0 outside the
!! selection) and a hard-label flag. Format/transactions follow simple_sigma2_state_file; header and scalars carry checksums.
!! Layout: header(512) | scalars(real64, NFIELDS+ncomp) | weights(real32, nptcls) | flags(int32, nptcls)
```

**`commanders/simple/simple_commanders_flex_pca.f90:10-29`** (WRONG; needs a ruling on UI vs command-line defaults)
```fortran
real, parameter :: COV_LP_DEFAULT = 16.0 !< default lp (A); smpd_target = lp / COV_LP_OVER_NYQUIST
!> rec_backend=pcg defaults: fixed 4-iteration budget of the warm-started basis M-step (rtol=0 disarms the rtol
!! and FLEX_PCG_XTOL stops). NOTE: differs from the UI defaults (gridding, 20, 1e-3). State solves use FLEX_PCG_STATE_MAXITS.
```
Also fix:
- line 110, which says "every positive budget from a ZERO start";
- line 282, which cites a non-existent `COV_SMPD_TARGET_DEFAULT`.

**`strategies/parallelization/simple_flex_pca_strategy.f90:221-234`** (WRONG)
```fortran
        ! Master thread boost: master and worker phases never overlap (the master blocks on each qsys round), so
        ! master-only stages use nparts*nthr threads, capped at owned cores (SLURM_CPUS_PER_TASK, else
        ! omp_get_num_procs). params%nthr follows because builder/reconstructor scratch is sized from it.
```

**`simple_flex_pca_model.f90` `run_flex_pca`** (WRONG)

Replace four blocks with one line each:
- **:138-142:** `! npreimages = state ceiling (>= MIN_NSTATES); preimage_auto=yes raises it to AUTO_NSTATES unless given and forces the merge.`
- **:222-228:** `! Paired fit: two mod-4-half fits advanced together, merged, polished once, embedded on all N; workers skip it.`
- **:242-249:** `! Axis weight from cross-half match cosine c: prior variance x 2c/(1+c), ~0 below 0.143.`
- **:472-478:** `! Drop states below min_neff before reconstruction (they gave artefact maps).`

Further fixes:
- `if( .true. )` at :285 is a dead conditional.
- :838 names `SIMPLE_COV_DECONV`, which has no reader.

## NU filtering

**`nu_filt/simple_nu_filter.f90:1-26`** (STALE)
```fortran
!@descr: volume-domain nonuniform filtering of even/odd volumes
! Sequence: setup_nu_dmats -> optimize_nu_cutoff_finds -> nu_filter_vols -> cleanup_nu_filter.
! Bank: static ladder lowpass_limits, cut at fsc_res/NU_BANK_FSC_HEADROOM when given (>= 2 rungs);
! an auxiliary (ML) pair at or beyond the finest retained rung is appended at that rung's Potts coordinate.
! The last label is the finest member and the matching low-pass handoff.
```
**Move to:** covered in `doc/policies/NU/nonuniform_filtering_policy.md` §8.

**`simple_nu_filter.f90:60-78`** (WRONG)
```fortran
! Static-bank cap: with fsc_res, keep rungs at or coarser than fsc_res/NU_BANK_FSC_HEADROOM (>= 2);
! without it (nu_filt3D, flex_pca) the bank is uncapped. The cap bounds the handoff's lead over the
! FSC from above; rationale in nonuniform_filtering_policy.md sections 8 and 12.
```
The policy repeats the wrong "1.25-1.5x" sentence at :532.

**`simple_nu_filter.f90:179-190`** (WRONG)
```fortran
! Diagnostic only: finest label whose cumulative population reaches this % of the signal voxels
! (mask minus the background clamp), quoted on the NU MATCHING LOW-PASS HANDOFF line.
! The handoff itself is the finest bank member (get_nu_filter_bank_finest_lp).
```

**`simple_nu_filter.f90:261-277`** (OK)
```fortran
! Evidence-envelope null. Spherical base pair: margin median/MAD over the observed support (no shell).
! Envelope-constrained pair (set_nu_evidence_null_shell): null on the density envelope's dilation
! ring at full base-support weight; labels free only on the observed density envelope.
```

**`simple_nu_filter_bank.f90:65-82`** (OK)
```fortran
! Append the ML pair as the last label when its Fourier index is at or beyond the finest retained rung
! (in practice only past 4 A: the fsc/1.5 cut keeps a rung finer than the FSC). It shares that rung's
! Potts coordinate, so the unary alone decides; its filtered pair is never cached.
```

**`simple_nu_filter_bank.f90:295-308`** (WRONG)
```fortran
! Coarse-to-fine like-for-like selection: label l replaces the incumbent only if strictly cheaper with
! both re-smoothed at l's radius from raw_dmats_mask; equal unaries keep the coarser label.
! Adds n(n+1)/2-1 smoothing passes; dmats_mask (own radius) still feeds the Potts prior.
```

**`simple_nu_filter_envmask.f90:1-27`** (WRONG)
```fortran
!@descr: NU-evidence-driven envelope masking for volume-domain nonuniform filtering
! Segments margin = max(0, raw C_coarsest - min_c raw C_c) (nu_ev_base/nu_ev_best), smoothed once at
! lp_smooth; optional scale-free base/best - 1 with a robust floor. Binary MRF by 8-colour ICM (beta =
! boundary area only); topology in simple_image_msk; write_nu_evidence_envmask is the single producer.
! Derivation: doc/algorithms/nu_evidence_envelope_mask.md.
```

**`simple_nu_filter_evidence.f90:104-121`** (OK)
```fortran
! Null candidate: zero cross-half prediction, smoothed at label 1's scale. Offset = lower quartile of
! C_zero - min(C_signal) over observed voxels (null-component centre, not a detection threshold);
! median/MAD are kept as diagnostics.
```
**Move to:** nothing to move. The policy is the stale side here: §8:343 still says
"median-plus-three-MAD".

**`simple_nu_filter_sharpen.f90:1-39`** (OK)
```fortran
!@descr: NU-evidence nonuniform postprocessing, v2 (classical pipeline, local)
! One Guinier B from the evidence-pair average (HPLIM_GUINIER to the finest evidenced cutoff, only if
! finer than NU_SHARP_BFAC_FINEST_A); per voxel sqrt(2FSC/(1+FSC)) stretched so FSC=0.143 sits at its
! evidenced cutoff, then Butterworth there; null/outside-support voxels take the map mean.
! One merged display map; design and v1 record: nu_evidence_local_sharpening.md section 3c.
```

**`volume/simple_nu_state_filter.f90:1-20`** (STALE)
```fortran
!@descr: assembly-owned nonuniform (NU) filtering of one state's half-map pair
! Shared by gridding volassemble and the PCG master: bank from the base pair with the ML pair as optional
! auxiliary member (setup_nu_dmats), the automsk background (density, or valid NU evidence for automsk=nu),
! _nu_filt/_nu_locres products, and the handoff (finest bank member).
! Contract: doc/policies/NU/nonuniform_filtering_policy.md.
```

**`simple_nu_state_filter.f90:52-73`** (OK)
```fortran
!> NU competition for one state on the base pair (consumed); the ML pair (consumed) is the auxiliary
!! member when l_use_aux. res0143 caps the bank and is the auxiliary's resolution; align_lp = finest
!! bank member. An optional apply pair (not consumed) receives the label field instead of the base pair.
```

**`commanders/simple/simple_commanders_postprocess_nu.f90:1-21`** (WRONG)
```fortran
!@descr: NU-evidence nonuniform postprocessing (isolated from the standard postprocess path)
! Evidence pair: the state's _even_unfil/_odd_unfil. With _even/_odd present the refinement's competition
! is rerun for <vol>_locres_nu. nu_evidence_sharpen_vol then sharpens the _solvent pair if present, else
! the unfil pair, to <vol>_pproc_nu (+_mirr). Display maps only, never FSC/resolution inputs.
! Design: doc/implementation_notes/completed/nu_evidence_local_sharpening.md.
```

**`commanders/simple/simple_commanders_volops.f90` `postprocess_volume_from_files`** (OK)
```fortran
!> Isotropic postprocess: cutoff = FSC=0.143 (FSC file, else the _unfil pair); Guinier B (lp < 5 A) from
!! the _unfil pair average, capped at the cutoff shell; sqrt(2FSC/(1+FSC)) unless provenance says ML-shrunk;
!! Butterworth at the cutoff; no mask if the support is already in the map. pair_stem: _unfil pair stem.
```
Replace the blocks at 392, 419, 434, 446, 460 and 489 with one line each. The measurements are
covered in `refine3D_policy.md:665-713` and `pcg_decision_log.md:911-940`.

Also, outside the audited blocks, stale NU comments describe:
- the deleted shell walk: `simple_nu_filter.f90:102-111`, `simple_nu_filter_bank.f90:204-222` and `:535-539`, `simple_nu_filter_potts.f90:38-40`, `simple_nu_filter_evidence.f90:80-82`;
- adaptive bands that cannot fire;
- the garbled sentence at `simple_nu_filter.f90:171-177`;
- product names: `postprocess_nu.f90:133` should say `_locres_nu`, and `sharpen.f90:200` should say `_pproc_nu`.

## Volume, abinitio and high-level tests

**`volume/simple_pcg_solvent_sidecar.f90:1-44`** (WRONG)
```fortran
!@descr: opt-in soft solvent prior of the PCG base solve (pcg_solvent=yes)
!  Per half: w(r) in [0,1] from its own prior-free map (Otsu + logistic on |x| smoothed at 2x FSC=0.143),
!  then a cold re-solve with ridge lambda_s(1-w). Prior-free pair = base pair (FSC, NU, evidence, _unfil);
!  prior'd pair = replay base and NU apply target. lambda_rel: closed-form cross-half CV unless given.
!  Rationale/history: doc/implementation_notes/completed/pcg_decision_log.md (2026-09-18..22).
```

**`volume/simple_frozen_accum.f90:1-25`** (OK)
```fortran
! Raw frozen-particle accumulators (gridding S/rho, PCG B/D per half), one set per box, added
! with coefficient one before any restoration or prior; nothing here scales or restores.
! Sets are bound to the run context (run id, backend, states, row counts, weighting) and
! validated against it and the consumer grid before any payload is read; defects are fatal.
```

**`volume/simple_frozen_accum_tester.f90:1-14`** (OK)
```fortran
! Frozen set + cohort equals the direct union (gridding sums/rho, PCG B/D, restored and solved
! maps); zero cohort gives F exactly. Plus context round trip and refusals. Box 16, seeded noise.
```

**`volume/simple_reconstructor.f90:698-715`** (OK)
```fortran
!>  Floor rho at its shell mean / frac (default 1000, RELION-style; clamped >= 1) before
!!  sampl_dens_correct, which divides unfloored. Needed for kernel-weighted (flex) state maps,
!!  whose rho is small and noisy in low-occupancy regions.
```

**`project/simple_project_superset.f90:1-23`** (OK)
```fortran
!@descr: (keep)
! Current and frozen projects share row indices; shared rows must name the same image (and CTF/optics
! for frozen rows), appended rows come from new stacks. frozen = frozen ptcl3D state>0 & updatecnt>0;
! cohort = current ptcl2D state>0 & not frozen. mask zeroes frozen rows; restore brings them back.
```

**`project/simple_project_superset_tester.f90:1-15`** (OK)
```fortran
! Current project 20 rows / 2 stacks, frozen project 14 rows (its first stack): identity refusals
! naming the first offender, membership and floors, mask and restore.
```

**`abinitio/simple_abinitio3D_addon_report.f90:1-14`** (OK)
```fortran
!@descr: (keep)
! Per state: union vs base FSC (verdict on the FSC=0.143 shell move beyond SHELL_TOL), correlation in
! the mskdiam sphere up to the base resolution, docking below DOCK_CORR_FLOOR, optional cohort-only
! check (addon_diag), stage limits; kind=stage|state|cohort text records that read restores.
```

**`abinitio/simple_abinitio3D_manifest_tester.f90:1-16`** (OK)
```fortran
! Manifest round trip and refusals, registration, >256-char paths, replay onto a cmdline and
! frozen-project validation; in-memory projects, one 8^3 map, a sigma2 stand-in with blanks.
```

**`abinitio/simple_abinitio_controller.f90:549-563`** (OK)
```fortran
! Stage-boundary FSC=0.5 promotion past FSC05_PROMOTE_MIN_STAGE, never finer than lp_cap: once per
! boundary, since without gold-standard halves only an FSC beyond the previous band is clean.
! Add-on promotes from the union FSC (abinitio3D_policy.md sec. 4; abinitio3D_addon_policy.md sec. 5).
```
The procedure's other blocks (536-545, 564-575, 682-690) repeat `abinitio3D_policy.md` §4.

**`commanders/test/simple_commanders_test_highlevel.f90:1159-1202`** (STALE)
```fortran
!>  \brief  Fail-fast gate of the PCG operator/solver, 14 stages (doc/policies/3D/reconstruct3D_pcg_policy.md sec. 9).
!  Stages 1-7 build data with forward_plane (inverse crime: operator algebra only);
!  stage 8 (deapodization, envelope-free data) is the only envelope check. No reconstructor/volassemble.
```
Also cut the history at 1245-1252 and the lab record at 3542-3552.

**`simple_commanders_test_highlevel.f90:2817-2831`** (STALE)
```fortran
!> Gridding vs PCG reconstruct3D on the same project/poses/sigma2 (numbered exec dir unless mkdir=no).
!! Hard gates: agreement band, in-band amplitude ratio and FSC, radial-ratio range; with truth vol1
!! also truth FSC (LS flatness only if ml_reg=no). Thresholds: implementation_notes/completed/drop_legacy_box_division.md
```

**`simple_commanders_test_highlevel.f90:3482-3500`** (OK)
```fortran
!> abinitio3D_addon gate on symmetry-broken 6VXX particles: abinitio3D on a seeded selection of the
!! first NBASE rows, then the add-on on all NPTCLS rows (same project basename). Gates frozen-input
!! integrity, manifests, cohort coverage and poses, union map vs truth and base; metrics.tsv.
```

**`commanders/simple/simple_commanders_refine3D.f90` `exec_refine3D_auto`** (WRONG)

For :64-74:
```fortran
        ! Registration pass: one refine=greedy iteration of all particles against the masked startup
        ! references, banded via lpstop (never lp: l_lpset would break gold standard).
```
For :211-232:
```fortran
        ! Startup: particle-power sigmas, then ONE regularized reconstruct3D with the refinement's filtering,
        ! so iteration 1 matches the masked, NU-filtered references every later iteration uses.
```
Also fix `refine=prob` at :285. `exec_bootstrap_rec3D` is accurate. The doc of
`prepare_bootstrap_rec_cline` (2198-2207) leaves out that `filt_mode` is kept under `automsk=nu`.

**`commanders/simple/simple_commanders_rec_distr.f90` `restore_state_from_parts` / `blend_trailing_accumulators`** (OK)
```fortran
! Population rule (population_blend_weights; partials note sec. 4.1): current *= u/f,
! chain *= (1-u)N/M, so the blended mass is N and the current-map coefficient exactly u.
```
```fortran
! No blend (no chain yet, or u ~ 1): persist the chain at full mass (partials x 1/f, M = N);
! a fractional-mass seed would over-weight the next update.
```

## Search, pftc, ori, class, params, image, misc

**`strategies/search/simple_ptcl_cache.f90:1-53`** (WRONG)
```fortran
!@descr: on-disk cache of noise-normalized, Fourier-cropped particles for the 2D matcher workflows
! An entry is the iteration-independent prefix of prepimg4align (noise norm vs lmsk, FFT, clip to box_crop),
! stored in real space (exact round trip). Covers every active particle. Read by prob_tab2D and cluster2D_exec;
! restoration then runs at box_crop with no second normalization (cavger_init_online(cropped_ptcls=.true.)).
! 3D workflows reject cache=yes. Contract: doc/policies/2D/particle_cache_policy.md
```
Fix the policy's 3D leftovers first (inventory G). Lines 440-457 (`ptcl3D` in
`oritype_cacheable`) and 751-754 (the `lmsk` fallback) are leftovers from 3D support.

**`simple_ptcl_cache.f90:322-339`** (OK)
```fortran
!>  Source fingerprint: each stack (name, range, nptcls_stk, size, mtime) plus every particle's stkind/indstk,
!!  so a stack replaced in place or a remapped project cannot validate a stale cache. One stat per
!!  stack + O(nptcls) folds; only the rank that decides whether to rebuild computes it.
```
`ptcl_cache_ensure` header:
`!> Master-side, before workers: build if missing/stale or adopt a valid leftover (the builder/adopter owns the files); any fallback sets cache=no on cline so all ranks agree.`

**`pftc/simple_pftc_shsrch_grad.f90:416-436`** (WRONG)
```fortran
!> Polish the caller's selected pose over (sx,sy,rotind_frac) within +/-2 cells of irot_in (euclid, cc, hybrid).
!! irot_in is authoritative: no all-angle rescan. An invalid or non-material solve returns the seed
!! cell re-scored at xy_in (irot never 0); evaluation_valid/improved are diagnostics only.
```
In the body:
- Lines 563-564 say "irot still holds the global scan result", which is stale.
- Each tolerance constant gets one line:
  - `JOINT_NEG_COST_TOL ! loss below -tol is a series artifact: invalid`
  - `JOINT_IMPROVE_REL_TOL ! gain must be material, not solver noise`
  - `JOINT_BOUND_TOL ! bound-pinned solutions never displace the seed`

**`pftc/simple_polarft_corr.f90:1548-1566`** (STALE)
```fortran
! cc counterpart of gen_raw_euclid_grad_at_angle in loss orientation (f=-cc, grad=-dcc); on integer
! grid indices it reproduces gen_corrs (unweighted shells, Nyquist included).
```
Above `gen_normalized_corr_grad_at_angle`, add:
`! cc = N/sqrt(D*C); D is shift-independent, so the quotient rule enters only d/dtheta`

**`ori/simple_oris_sampling.f90:78-95`** (WRONG: edge case)
```fortran
!> Population-rule weights of one group: new = s*current + w*previous; f = n/N, u = ufrac or f, s = u/f,
!! w = (1-u)*N/M (0 if M = 0); mnew = s*n + w*M (= N when M > 0). n = 0 keeps previous at mass N; w > 1 if rows return.
!! nrep = N active updated rows, nsmp = n sampled, mrep = M stored population.
```

**`class/simple_cavg_sums.f90:1-15`** (OK)
```fortran
!@descr: unregularized 2D class Fourier sums on disk: the partless carry-over set and per-worker contributions
! Per class: even/odd numerators and CTF^2 sums, captured before restoration. STATE: carried set, owner-written,
! records M(c). CONTRIBUTION: one worker's from-zero sums + centering offsets, e/o pops, l_frac.
! Written via temp name + rename, closed by a payload byte count. Ownership: doc/policies/2D/abinitio2D_policy.md
```

**`params/simple_parameters.f90:1-31`** (OK)
```fortran
!@descr: public parameters type and interfaces for parameter parsing and derivation phases
! New parameter: declare here (scripts/simple_args_generator.pl derives the accepted-argument list from these
! declarations), register in simple_parameters_parse.f90, derive/validate in simple_parameters_phases.f90
! (phase order: new()), expose in src/main/ui. Dynamic type(string) defaults: init_dynamic_defaults (core).
```
Three inline descriptors are stale:
- :222: both backends are valid;
- :303: the default disagrees with the UI;
- :356: the `refine` value list is wrong.

**`image/simple_image_tester.f90:690-702`** (OK)
```fortran
! rmat's first dim is the in-place FFTW buffer (n1+2 even / n1+1 odd); in real space the extra rows must be
! exactly zero, so every tolerance is 0. fft/ifft need even boxes; fft_noshift takes any.
```

**`nano/simple_atoms.f90:12-36`** (WRONG)
```fortran
! PDB ATOM/HETATM fixed columns per the wwPDB v3.3 coordinate section; segID (73-76) is skipped.
! Extended (>99999 atoms / >9999 residues): 'ATOM'/'HETA' + I7 serial (5-11), I5 resSeq (23-27), no iCode.
```
Drop the commented-out legacy formats at 35-36.

**`commanders/simple/simple_commanders_imgops.f90:505-521`** (OK)
```fortran
! pca_mode=ppca_kpca_resid: PPCA gives x_ppca, kPCA models only r = x - x_ppca; output is
! avg + x_ppca + alpha*r_kpca (alpha = ppca_kpca_resid_alpha). transp_pca=no only.
```

**`utils/filter/simple_lpstages_tester.f90:1-10`** (OK)
```fortran
!@descr: unit tests for the low-pass/crop stage schedules, FSC weightings, B-factor cap and Butterworth filter
! Stage constants come from an out-of-tree double-precision emulation (lpstages_ref.py, not in the repo);
! FSC weightings are closed forms.
```

**`image/simple_image_calc.f90` `nu_objective`** (OK)
- Header: `! Cross-half Huber unary, whitened by the radial E/O noise profile (interpolated between shell centres)`
- `HUBER_DELTA ! L2 at noise scale, L1 for outliers`
- `HUBER_LOSS_CAP ! saturate in dp; keeps the smoother stable beside hard-support zeros`

## Stream, GUI, persistent worker

**`stream/simple_stream_p03_initial_analysis.f90:1-35`** (WRONG)
```fortran
!@descr: stream stage 3: two-cycle opening analysis that generates picking references
! Cycle 1 (NMICS_PLAN(1) mics): segdiam pick -> extract -> abinitio2D -> cavg quality selection.
! Cycle 2: pick all mics with cycle-1 bins -> ptcl_sieve -> abinitio2D -> quality -> balance_classes
!          -> abinitio3D_cavgs -> reproject. Exits when done or on a GUI pickrefs selection.
! GUI: progress on ipc_pipe_initial_analysis_in, selections on ipc_pipe_initial_analysis_out.
```

**`stream/simple_stream_p04_refpick_extract_new.f90:1-27`** (WRONG)
```fortran
!@descr: stream stage 4: reference-based picking and extraction
! Waits (<=24 h) for pickrefs, builds templates with make_pickrefs on the first preprocessed
! project (smpd known only then), then queues one pick_extract job per new project and merges
! results into the stage project. ctfres/icefrac/astig thresholds apply only if reject_mics=yes.
```

**`stream/simple_stream_p05_sieve_cavgs_new.f90:1-40`** (STALE)
```fortran
!@descr: stream stage 5: continuous particle sieving via ptcl_sieve (coarse, optional fine, cavg rejection)
! Imports completed refpick projects (<=MAX_MOVIE_IMPORT per loop) and calls sieve%cycle each loop.
! The sieve is created on first import (mskdiam from pickrefs). Final ingestion is set after
! FINAL_INGESTION_IDLE_TIME without imports and cleared on the next import.
```

**`stream/simple_stream_p06_pool2D_new.f90:1-41`** (STALE)
```fortran
!@descr: stream stage 6: global 2D classification of the sieved-particle pool
! Imports sieved sets when the pool is idle (stepwise=yes: stop once the pool reaches nptcls_threshold),
! iterates with an adaptive pause policy, and from EXPORT_3D_START_ITERATION exports new particles
! to DIR_STREAM_COMPLETED for stage 7. GUI in: mskdiam2D, sieverefs, snapshot2D.
```

**`stream/simple_stream_p07_abinitio3D_multistate.f90:1-36`** (STALE)
```fortran
!@descr: stream stage 7: multistate 3D from pool2D exports (abinitio3D, then abinitio3D_addon)
! Imports stage-6 exports (each filtered by model_cavgs_rejection), runs one NSTATES3D abinitio3D,
! then an abinitio3D_addon pass whenever the pool has grown. Ingestion pauses while a job runs.
! After abinitio3D, sends per-state gui_metadata_vol3D and reprojection tiles to the GUI.
```

**`sieve/simple_ptcl_sieve.f90:1-37`** (STALE)
```fortran
!@descr: multi-tier particle sieve with coarse/fine 2D chunking and rejection
! cycle() = collect_and_reject -> generate_chunks_coarse -> generate_chunks_fine (unless single_pass)
! -> submit (fine first). Restart state comes from per-chunk sentinels
! (ABINITIO2D_FINISHED, REJECTION_FINISHED, COMPLETE). Contract: ptcl_sieve_policy.md
```
Further fixes:
- The `new()` comment (:221) says "four output directories"; the code creates two.
- `SIMPLE_CHUNK_PARTITION` (sieve) and `SIMPLE_STREAM_CHUNK_PARTITION` (p05) are two names for the
  same thing.

**`utils/gui/simple_gui_assembler.f90:1-42`** (STALE)
```fortran
!@descr: builds the GUI JSON document (stream and batch) from gui_metadata objects
! Each assemble_* replaces its section of json_root. Content sections are dropped when their FNV-1a
! hash matches the last one sent (heartbeats always go). Call clear_hashes() after a failed send.
```

**`utils/gui/metadata/simple_gui_metadata_vol3D.f90:1-38`** (STALE)
```fortran
!@descr: GUI metadata type for a single 3D volume entry (product paths + stats).
! Optional fields (res0143/res05/cfar/pop, FSC, 72x36 oridist, oridistpath, per-kind MRC min/max,
! reprojtiles) are emitted only when set; i/i_max route IPC batches.
! reprojtiles is allocatable: never set it on objects sent via raw serialise() (review A6).
```

**`utils/gui/metadata/stream/simple_gui_metadata_stream_update.f90:1-36`** (WRONG)
```fortran
!@descr: GUI metadata for a stream quality update — thresholds and user selections broadcast from the GUI
! The master copies GUI response fields here and writes the whole object to every stage pipe.
! Unset values are 0; receivers ignore 0 and unchanged values. Readers: p01 thresholds,
! p03 pickrefs_selection+cycle, p06 mskdiam2D/sieverefs/snapshot2D.
```
The routine comments at :185 and :194 call `pickrefs_cycle` "an integer array"; it is a scalar.

**`utils/persistent_worker/message/simple_persistent_worker_message_base.f90:1-34`** (WRONG)
```fortran
!@descr: polymorphic base type of the persistent-worker wire messages
! msg_type leads every message; 0 = uninitialised. Subtypes set it in new().
! Messages travel as their raw memory image. Each subtype overrides serialise() with its own
! copy of the transfer body; do not send extended types through this routine.
```
The claimed `sizeof` truncation is untested: the tester never sends an extended type through the
base routine. Fix `persistent_worker_policy.md` §6 to match.

**`..._message_heartbeat.f90`, `_status.f90`, `_task.f90`, `_terminate.f90` (headers)**
(heartbeat WRONG, status STALE, task STALE, terminate WRONG)

The pattern for all four is the same. Heartbeat:
```fortran
!@descr: heartbeat message of the persistent-worker protocol: a worker reports liveness and thread load
! Worker -> server each loop: worker_id, worker_uid, heartbeat_time, nthr_used/total (fd set server-side).
! Reply: task, STATUS idle, or TERMINATE. serialise() override: see ..._message_base.
```
Corrections for the other three:
- **Status:** it is also the `queue_task` acknowledgement (idle = queued, error = queue full).
- **Task:** dispatch goes out as `WORKER_NEW_TASK_MSG`, and the timing fields are never set.
- **Terminate:** the worker *cancels* running tasks. Terminate is also sent on scale-down and on
  ID/UID mismatch.

In all four, replace the copied DESIGN CONTRACT paragraph with a one-line pointer to the base.

**`utils/persistent_worker/simple_persistent_worker_server.f90:1-39`** (STALE)
```fortran
!@descr: TCP poll-loop server that hands tasks to persistent SIMPLE worker processes in reply to their heartbeats
! Per heartbeat the listener replies TASK, STATUS idle, or TERMINATE (shutdown, UID clash, scale-down).
! queue_task() submits WORKER_NEW_TASK_MSG over a self-connection.
! Invariants: server owns the mutex; queues/job_count are listener-local; no socket I/O under the
! mutex; port is zeroed only after the join.
```

**`production/simple_persistent_worker.f90:1-30`** (WRONG)
```fortran
!@descr: SIMPLE persistent worker: heartbeats a worker server and runs dispatched bash scripts in pthread slots
! Args: server=<ips> port=<n> (required), worker_id=<slot from qsys_env>, nthr=<slots, default 1>.
! TASK -> free slot runs `bash script`; STATUS ignored; TERMINATE or lost server -> cleanup,
! which cancels unfinished scripts. Script-path validation is currently bypassed.
```

**`utils/simple_forked_process.f90:1-30`** (STALE)
```fortran
!@descr: POSIX fork-based child-process manager with timestamps, auto-restart, and status polling
! Extend and override execute(cline); start() forks and runs it in the child (exit 0).
! status() polls waitpid(WNOHANG) and, if restart=.true., re-forks a failed child (side effect).
! A non-zero wait status (including signals) is FAILED. terminate()=SIGTERM, kill()=SIGKILL.
```
