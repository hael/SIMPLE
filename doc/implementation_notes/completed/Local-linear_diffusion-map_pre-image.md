# Local-linear diffusion-map pre-image — validation plan

## Current implementation

`flex_analysis` now defaults to `preimage_mode=linear`.  It retains the
local-constant (Nadaraya–Watson) path as `preimage_mode=constant` for
comparison.  The local-linear implementation is local at each medoid and at
each Fourier voxel; it does not fit a global linear subspace.

## Implemented numerical contract
Per state `s` with target `z_c`, per Fourier voxel, fit weighted least squares of the voxel value on the **centered** latent coordinate `u_i = z_i − z_c`, weight `w_i = w_i(z_c)`:
- design row per particle: `[ 1, u_i(1..d) ]`  (length `d+1`)
- normal equations: `A = Σ w_i · gᵢ gᵢᵀ · CTF_i²`, `b = Σ w_i · gᵢ · D_i`
- solve `A x = b`; **pre-image = intercept `x(1)`** (constant term). Basis slopes `x(2:)` are local tangents and are not written as state volumes.

The active path folds `w_i` into `data_scales = w_i·g_i` and
`density_scales = w_i·(g_i g_iᵀ)` for the existing coupled batched
accumulator.  `A` is the packed second-moment block and `b` the RHS block.  A
relative ridge and an observability floor keep rank-deficient cells finite.

---

## Completed checks

- T1 is implemented with a numerical tolerance: an all-zero local coordinate
  recovers the local-constant weighted average.
- T2 uses a one-sided affine response and a curved response.  The local-linear
  intercept recovers the affine target and materially improves on the
  local-constant curved estimate.  It does not claim exact recovery of a
  quadratic response.
- T5 checks finite fallback on a rank-deficient synthetic design.
- T6 checks that a dense state's bandwidth is independent of an unrelated
  sparse state's population.

These checks run from `production/tests/simple_test_flex_projected_latent_model.f90`.

## Remaining tests

### T1 — Reduction: local-linear ⊇ local-constant
Completed; retain as a regression test.  The current tolerance is numerical,
not a bit-for-bit guarantee through the ridge fallback.

### T2 — Bias removal on a known curved ground truth (the key test)
Completed in the current synthetic form.  Extend it only with analytically
justified tolerances; local linear removes first-order design bias but is not
an exact quadratic interpolator.

### T3 — Boundary / endpoint behavior
Place `z_c` at the extreme of the sampled range (one-sided neighborhood). N-W bias is worst here; local-linear must still recover the linear trend to tol. Explicitly assert local-linear endpoint error ≤ interior error × small constant. Endpoints are your rare conformations — this test protects them.

### T4 — Adjoint / accumulation contract (extend, don't replace)
You already have `test_cartesian_projection_contract` and `test_coupled_batch_accumulation`. Add a case where `data_scales`/`density_scales` carry the `w_i·g_i` / `w_i·g_iᵀ` local-linear payload, and assert the single-record and batch accumulation paths still agree to 2e-5. Reuses existing gather/splat adjoint machinery.

### T5 — Per-voxel solve regularization (high-freq sparse cells)

The finite rank-deficient check is implemented.  A future test should measure
the current relative-ridge bias on a well-conditioned cell and define any
change to the ridge policy.  The implementation does not presently make
`lambda` a function of local `neff`.

### T6 — Per-state bandwidth, not global-min-gated

Implemented.  The regression builds two problems with the same dense cluster
and different sparse populations; the dense-state bandwidth must be unchanged.

### T7 — Frequency-dependent bandwidth (if implemented)
If bandwidth is made shell-dependent: assert low-freq shells use a narrower effective kernel (fewer particles, sharper) and high-freq shells wider, and that a flat-bandwidth run is recovered as the degenerate case. Skip if not yet implementing per-shell.

### T8 — Statistical validation harness (not a unit test — a routine)
Implement as a callable diagnostic, assert on synthetic data:
- **Half-set FSC:** split each state's particles in two, reconstruct both, inter-state FSC must exceed intra-state noise floor on curved synthetic data.
- **Permutation null:** randomly relabel particles across states; assert real between-state volume variance > shuffled variance by a set margin. If it fails on real data, the states are consensus + noise.

---

## Notes for the implementer
- `d` (latent dim) is small; `A` is `(d+1)×(d+1)` per voxel — direct symmetric solve (Cholesky with ridge fallback), no iterative solver.
- Keep everything in `real(dp)` for the moment accumulation as the existing code does; the `A` solve especially.
- Do NOT introduce a global mean/basis volume — that reintroduces the global-linear assumption you explicitly rejected. The linear fit is **local per target**.
- Tolerances above assume noise-free synthetic inputs; loosen only for the explicitly noisy T8 harness.
