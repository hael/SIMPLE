# Local-linear diffusion-map pre-image — implementation & test spec

## Goal
Replace the local-**constant** (Nadaraya–Watson) centroid pre-image with a local-**linear** one that removes the O(h²) curvature bias (the "collapse to consensus" failure), while reusing the existing per-voxel moment accumulation. Manifold curvature is preserved — do NOT fit a single global linear subspace.

## Core change (context for the tests)
Per state `s` with target `z_c`, per Fourier voxel, fit weighted least squares of the voxel value on the **centered** latent coordinate `u_i = z_i − z_c`, weight `w_i = w_i(z_c)`:
- design row per particle: `[ 1, u_i(1..d) ]`  (length `d+1`)
- normal equations: `A = Σ w_i · gᵢ gᵢᵀ · CTF_i²`, `b = Σ w_i · gᵢ · D_i`
- solve `A x = b`; **pre-image = intercept `x(1)`** (constant term). Basis slopes `x(2:)` are the local tangent.

Reuse `accumulate_plane_oversamp_coupled_stats`: fold `w_i` into the scales and pass `data_scales = w_i·g_i`, `density_scales = w_i·(g_i g_iᵀ)`. `A` = the packed second-moment block; `b` = the RHS block.

---

## Tests to implement (in priority order)

### T1 — Reduction: local-linear ⊇ local-constant
Set all `u_i = 0` (i.e. `z_i = z_c` for every particle). The local-linear intercept MUST equal the existing Nadaraya–Watson weighted average bit-for-bit (tol 1e-5). Guards that the new path collapses to the old one when there is no local variation.

### T2 — Bias removal on a known curved ground truth (the key test)
Synthesize `d=1` toy: voxel value `D(z) = a + b·z + c·z²` sampled at known `z_i`, no noise.
- N-W estimator at `z_c` must show error ≈ `c·h²·(kernel 2nd moment)` — nonzero, sign of `c`.
- Local-linear estimator at `z_c` must recover `a + b·z_c` to tol ~1e-4 (residual O(h⁴), not O(h²)).
Assert local-linear error < N-W error by the expected large factor. **This is the test that proves the method does what it's for.**

### T3 — Boundary / endpoint behavior
Place `z_c` at the extreme of the sampled range (one-sided neighborhood). N-W bias is worst here; local-linear must still recover the linear trend to tol. Explicitly assert local-linear endpoint error ≤ interior error × small constant. Endpoints are your rare conformations — this test protects them.

### T4 — Adjoint / accumulation contract (extend, don't replace)
You already have `test_cartesian_projection_contract` and `test_coupled_batch_accumulation`. Add a case where `data_scales`/`density_scales` carry the `w_i·g_i` / `w_i·g_iᵀ` local-linear payload, and assert the single-record and batch accumulation paths still agree to 2e-5. Reuses existing gather/splat adjoint machinery.

### T5 — Per-voxel solve regularization (high-freq sparse cells)
In high-resolution shells `A` is near-singular (few effective particles). Test that:
- ridge term `A + λI` with `λ` tied to local `neff` keeps the solve finite (`ieee_is_finite` on all outputs) for a deliberately rank-deficient `A`;
- as `λ → 0` on a well-conditioned cell, result → unregularized solution (tol 1e-5).
Assert graceful degradation to the **intercept-only** (= N-W) answer when `A` is fully degenerate, never NaN.

### T6 — Per-state bandwidth, not global-min-gated
Regression guard on the bug in `build_flex_preimage_kernel_weights`: construct 2 states, one sparse. Assert the dense state's bandwidth is NOT inflated by the sparse state (each state's `bw_scale` is independent, or capped). Currently `bw_scale` is global and gated on `minval(neff)` — this test should FAIL against the old code and PASS after the fix.

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