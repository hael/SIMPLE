---
name: simple-main-sigma2
description: Use when working in SIMPLE's src/main/sigma2 subsystem: the euclid_sigma2 noise model of objfun=euclid (build%esig: two global spectra, one per even/odd half, and a part's per-particle residual records), the canonical sigma2 store (sigma2_state.bin, generations, candidates, local range files, the streaming commit and publication, concatenation, solve3D_addon extend and union means), the calc_pspec power-spectrum bootstrap, ensure_sigma2_for_iteration, the residual refine=sigma pass, the 2D scoring-noise annealing (sigma2_anneal: the smoothed bootstrap anchor and its log blend with the measured spectra), and who refreshes the calculators' sigma2 weights.
---

# `src/main/sigma2`: the Euclidean noise model and its store

Owns the noise power of `objfun=euclid`: `euclid_sigma2` (`build%esig`: the two global spectra, column
`sigma2_half(eo)`, and a part's per-particle records), the store's transaction (`simple_sigma2_state.f90`),
the power-spectrum bootstrap's part sums (`simple_sigma2_pspec_sums.f90`) and the one rebuild entry
(`simple_sigma2_bootstrap.f90`, a controller that runs the `calc_pspec` commander), and the 2D scoring
anchor (`simple_sigma2_anchor.f90`). The file format is
`src/fileio/simple_sigma2_state_file.f90`, the consumability check `src/fileio/simple_sigma2_files.f90`, the
commit step `commit_sigma2_exec` (`strategies/parallelization/simple_strategy_sigma2.f90`); the scoring that
weights by the spectra is `pftc/` and `cartft_calc`.

## Rules

Silent ones first.

- The commit takes the particles' states and halves from the caller's particle field (merge the workers'
  documents first): each half's spectrum is the mean of its active records whose every shell is finite and
  positive. Invalid records are skipped with one warning; an empty half stops.
- After every load, refresh what caches the spectra: `prep_sigmas_objfun` calls `pftc%memoize_sigma2_weights`
  after `read_global` (a Euclidean score on spectra newer than the weights stops: `new_cart`, `read_global`
  and `kill` bump `spectra_gen`); the Cartesian calculator takes them per batch (`cftc%set_sigma2`); `read_global` drops
  the plane grid that `prep_fplane_grid` rebuilds for class sums and reconstruction.
- `esig%new(params, pftc, ...)` gives the calculator a pointer to `sigma2_noise`: make it after `pftc%new`,
  which drops the pointer. `fromp`/`top` are copied from `params` at `new`.
- One maker per generation, verifying consumers. Workers never write the store: a matcher reads its rows
  (`read_part`, `fromp:top`) and the spectra (`read_global`), overwrites the rows of the particles it aligned
  over its band (`calc_sigma2`) and writes its whole range for the next generation (`write_sigma2`). The master
  opens the transaction before the workers (`prepare_canonical_sigma_update`) and commits after them: the
  ranges must belong to that generation and layout, not overlap and cover every row, or the commit stops. So an
  empty partition writes its unchanged range too.
- The store spans shells `1..fdim(box)-1` of the native box and the project's particle layout (box, smpd,
  particle count, layout digest); a pass updates only its band. Its path is relative to the project's
  directory (`update_projinfo`).
- A follow-up pass opens no transaction: `read_pending_range` starts from the discrete pass's range and
  `write_sigma2` replaces it atomically.
- A workflow that consumes or updates the store first calls `ensure_sigma2_for_iteration`: nothing when the
  store is consumable, otherwise `calc_pspec` seeds and publishes it and registers its path, after which the
  caller re-reads `projinfo`. Where no Euclidean pass follows, `prepare_residual_sigma2_pass_cline`
  (`refine=sigma`) upgrades the seed.
- `sigma2_anneal=yes` (2D only, `objfun=euclid`): matching scores against `exp((1-w) log anchor + w log measured)`;
  restoration, reconstruction and residual records always use the measured spectra (`sigma2_noise`); only the
  scores, self-terms and posterior temperatures (so the insertion weights) see the blend. The anchor is
  each half's bootstrap, log-smoothed by a seven-point quadratic Savitzky-Golay filter, saved once by the refine2D
  owner (`ensure_canonical_sigma2` -> `prepare_sigma2_anchor`) as the companion `<store>.anchor`, only while the
  store still holds the bootstrap; a missing companion stops a continuation and is never made from residuals.
  Workers only read it (`read_sigma2_anchor`, checked against the store's grid and layout).
- `anneal_for_scoring` (called by `prep_sigmas_objfun` instead of `memoize_sigma2_weights`) points the polar
  calculator at the blend in `sigma2_scoring`, or at `sigma2_noise` itself when `w = 1`, and bumps `spectra_gen`.
  The weight `sigma2_weight` is set per iteration by `exec_refine2D` from `sigma2_stage_weight` and the internal
  stage bounds `sigma2_weight_start|end`, and goes to workers on the job line.
- Halves: `sigma2_half(eo) = eo + 1`; an eo outside 0..1 stops, in a worker and in the commit.
  `cc_emit_sigma=yes` (cc scoring that writes residuals; an `internal:` key, ruled to be retired) needs a
  seeded store and refuses `euclid`.

## Where to start

- `simple_euclid_sigma2.f90` (a pass's view), `simple_sigma2_state.f90` (transaction names, assembly,
  concatenation for merged projects and chunks, `sigma2_state_extend`/`union_means` for `solve3D_addon`),
  `simple_sigma2_bootstrap.f90`, `simple_sigma2_anchor.f90` (anchor companion, `blend_sigma2`,
  `sigma2_stage_weight`, the `>>> SIGMA2 ANNEAL` readout's band statistics); callers that show the order: `prep_sigmas_objfun`
  (`strategies/search/simple_matcher_ptcl_batch.f90`) and the masters' `execute_iteration`.

## Tests

- `build/bin/simple_test_exec test=unit_pftc_align2D3D suite=sigma2_state` (transaction, coverage, recovery,
  concatenation, add-on copy-on-write, invalid records; under a second) and `suite=sigma2_bootstrap` (part
  sums, command lines, the rebuild entry) and `suite=sigma2_anchor` (smoothing, companion persistence and
  refusals, the log blend, scoring versus reconstruction separation).
- Workflow: `simulated_workflow` (bootstrap and commits in solve2D and solve3D), `solve3D_addon`.

## See also

`doc/policies/3D/refine3D_policy.md` (the sigma2 bootstrap rule), `doc/policies/3D/solve3D_addon_policy.md`,
`doc/policies/layering_and_project_ownership_policy.md` section 5.7, `doc/policies/2D/solve2D_policy.md`
(scoring-noise annealing), `doc/implementation_notes/planned/sigma2_anchor_anneal_2D.md`; skill `simple-main-pftc`.
