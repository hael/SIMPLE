# Continuous in-plane rotation in the polar representation (`objfun=euclid`)

> Historical mathematical design note. Its discussion of adding a continuous
> angle callback predates the current policy. Production exposes only
> `inpl_cont=no|yes`: `no` uses the legacy discrete-angle callback, while
> eligible `yes` uses the fused joint `(sx,sy,rotind_frac)` evaluator at every
> former callback site, carries rounded angles through probabilistic artifacts,
> and persists a fractional angle only after final assignment. Each joint solve
> first discards prior in-plane metadata and performs one callback-equivalent
> all-angle selection at the supplied shift; it does not use the legacy 5-by-5
> shift initializer. The mathematics
> below remains the basis of the joint evaluator; callback-oriented sequencing
> is not current policy.

Implementation proposal. The in-plane angle is the only alignment parameter
SIMPLE never refines below its search grid. In the polar representation a
continuous rotation is a *phase factor* — exact, analytically differentiable,
requiring no interpolation — and for the Euclidean objective the angular
dependence turns out to be a **single trigonometric polynomial**, which makes
this simpler than it would be for `cc`.

Scope: `objfun=euclid` throughout. `cc` and the hybrid/denoised variants are
follow-ons, §10.

## 0. The claim in one line

Shifts and in-plane rotations have the same structure in the polar
representation:

```text
shift    by s:  multiply by exp(-2πi·(h,k)·s / box)   ->  ∂/∂s = -2πi·(h,k)/box · (·)
rotation by θ:  multiply by exp(-i·m·θ)               ->  ∂/∂θ = -i·m · (·)
```

SIMPLE already exploits the first. The second is available at low cost, because
the angular Fourier coefficients are **already computed and memoized**.

## 1. Current state and the size of the gap

| Parameter | How it is searched | Where |
| --- | --- | --- |
| shift `(sx, sy)` | continuous, analytic gradients, L-BFGS-B | `simple_pftc_shsrch_grad.f90` |
| in-plane angle | **discrete index only**, `maxloc` over `nrots` | `grad_shsrch_optimize_angle`, `simple_pftc_shsrch_grad.f90:165` |

`grad_shsrch_minimize` alternates: L-BFGS-B on `(sx,sy)` at frozen integer
`cur_inpl_idx`, with `opt_callback => grad_shsrch_optimize_angle_wrapper`
re-snapping the angle to the best grid index between shift solves. `irot` is an
integer through the whole call chain.

`magic_pftsz_1` (`simple_magic_boxes.f90:34`) sets `pftsz ≈ π·msk`, and
`nrots = 2·pftsz` (`simple_polarft_core.f90:61`), so `nrots ≈ 2π·msk` — the arc
step between angular samples is **one pixel at the mask radius**. Worst-case
residual after snapping is half a step, so density at radius `r` is misplaced by
up to `r·π/nrots`: ~0.5 px at `r = msk`, proportionally less toward the centre.

Two consequences. It is the same order as the shift precision the code already
works hard for, so refining shifts sub-pixel while leaving the angle on a grid is
unbalanced. And because the error scales with `r`, it is a **radially increasing
blur** — worst exactly at the periphery, where the high-resolution information is.

## 2. Why polar rather than Cartesian

A rotation maps every ring onto itself, so it is a pure translation along the
angular axis: a phase in angular-frequency space, no interpolation and no
interpolation *derivative*. A Cartesian formulation (the `dgather_window`
machinery of `continuous_3D_refinement_on_pcg_operator.md`) would pay a KB gather
derivative — and its non-C¹ behaviour at the window edge — for a degree of
freedom polar gets exactly.

**Corollary, so nobody is surprised later:** this does **not** de-risk the 3D
continuous-refinement design. In-plane rotation is the one rotational DOF that is
a phase; out-of-plane rotation moves between central sections and is a genuine
interpolation on the volume. That is where `dgather_window` and the non-C¹
problem live, and none of it appears here. Two separate projects.

## 3. The mathematics, in SIMPLE's conventions

### 3.1 Storage

- `pfts_*(pftsz, nk, n)` — `pftsz` angular samples spanning `[0, π)`, `nk` shells.
- The full `nrots` ring is recovered by Friedel extension, as `gen_euclids` does
  it (`simple_polarft_corr.f90:308`): `c(pftsz+1:nrots) = conjg(c(1:pftsz))`.
- `dang = 2π/nrots`, `angtab(irot) = (irot-1)·dang` (`simple_polarft_core.f90:65`).
- `polar(1,k,irot) = k·sin(φ)` (h), `polar(2,k,irot) = -k·cos(φ)` (kc)
  (`:70-71`); `argtransf` is these scaled by `2π/ldim` (`:91-92`).
- Rotation by `irot` is the index map `rp = p - irot + 1` with Friedel
  conjugation on wraparound (`simple_polarft_corr.f90:1523,1536`), i.e.
  `REF(φ_p − θ)` with `θ = (irot−1)·dang`. **Note the minus sign**; it propagates
  into every derivative below.

`polar` is a table over integer `irot`, but at continuous `θ` the same quantities
are closed-form (`k·sin(φ_p−θ)`, `−k·cos(φ_p−θ)`), so **the shift factor is exact
at any angle**.

### 3.2 The euclid residual is a single trig polynomial

This is the structural reason euclid is the easier case. `gen_euclids`
(`simple_polarft_corr.f90:282`) builds one half-spectrum and inverse-transforms
it once:

```text
crvec1%c  =  Σ_k  w_k · [ FT(CTF²)·FT(REF²)  −  2·FT(X·CTF)·conj(FT(S·REF)) ]
w_k       =  k / σ²(k, iptcl)                                    ! :323
r(grid)   =  c2r( crvec1%c )                                     ! :367
euclids   =  exp( −1 − r / A ),   A = wsqsums_ptcls(i)·2·nrots    ! :369-370
```

So the weighted squared residual over rotation is **one** trig polynomial, not a
ratio. Evaluating it off-grid is exact:

```text
ω     = 2πθ / nrots                                   (θ in grid units)
r(θ)  = R_0 + 2·Σ_{m=1}^{pftsz-1} Re[ R_m e^{imω} ] + Re[ R_pftsz e^{i·pftsz·ω} ]
r'(θ) = (2π/nrots)·( 2·Σ_{m=1}^{pftsz-1} Re[ i·m·R_m e^{imω} ] + Re[ i·pftsz·R_pftsz e^{i·pftsz·ω} ] )
```

`O(pftsz)` per evaluation, derivative the same loop with one extra factor, **no
new FFT**.

For optimisation, do **not** work with `exp(−1 − r/A)`. It is a strictly
decreasing function of `r`, so minimising `r(θ)` directly is equivalent and
removes the exponential from the inner loop entirely. The 1D refinement becomes a
root-find on `r'(θ) = 0` — a real trig polynomial with an analytic derivative.
Convert to the objective only when reporting.

Contrast with `cc`, where `gen_corrs` produces a *ratio* `N(θ)/sqrt(D(θ)·…)` and
every derivative needs the quotient rule. Euclid avoids that.

### 3.3 The two routes are algebraically identical — use it as an anchor

`gen_euclid_residual_grad` (`:1499`) computes the same quantity by direct sum:
`f = Σ_k w_k Σ_p |Rot(Shift(REF))·CTF − PTCL|²`, then
`gen_euclid_grad_for_rot_8` (`:1583`) forms `exp(−f/denom)` with
`denom = wsqsums_ptcls(i)`.

Expanding the square gives `f = r + Σ_k w_k|PTCL|² = r + denom`, hence
`exp(−f/denom) = exp(−1 − r/denom)` — which is exactly the `−1` in
`gen_euclids:370`, with `A` carrying FFTW's extra `2·nrots`. The two routes agree
by construction. **That identity is the strongest validation anchor available**
(§6, test 1) and it costs nothing to check.

### 3.4 Nyquist — and an open question in the current code

`nrots ≈ 2π·msk` means the angular representation is **critically sampled**: the
`m = pftsz = nrots/2` coefficient sits exactly on Nyquist. Three points:

1. **Weight.** In the `c2r` half-spectrum convention the Nyquist term carries
   weight **1**, not 2 like the interior terms, and interpolates as
   `R_pftsz·cos(πθ)`. Giving it weight 2 is the obvious slip; test 1 catches it.

2. **`gen_euclids` appears to drop it already.** The accumulation loops run
   `do p = 1, self%pftsz` (`:325, :335, :349, :358`) while `crvec1%c` is
   dimensioned `pftsz+1`, so index `pftsz+1` — the Nyquist bin — is left at the
   zero it was initialised to. `gen_corrs` by contrast uses whole-array
   expressions (`:258, :273`) and *does* include it. **Resolve this before
   building on it.** If the omission is deliberate (zeroing Nyquist is a
   legitimate way to remove the ambiguity, and `pftsz` is the vectorisation-
   friendly length) then hazard 1 disappears for euclid and the interpolant is
   unambiguous — good news, worth writing down. If it is an oversight, the
   grid values themselves carry a small systematic, and that is worth knowing
   independently of this project. Either way: **do not "fix" it silently**, and
   do not let the answer change under you mid-implementation.

3. **Aliasing — the error mode no grid test can catch.** The trig-polynomial
   interpolant is the true residual only if the signal is genuinely band-limited
   to `msk`. A hard mask edge, density beyond the radius, or noise puts angular
   content above `m = nrots/2`, which aliases. On-grid values stay exactly right;
   **off-grid values are wrong and no comparison against `gen_euclids` will ever
   reveal it.** Test 3 exists solely for this.

If aliasing proves material, the fix is to raise `nrots` above `2π·msk` for the
refinement stage only — a parameter change, not a redesign.

### 3.5 The three gradients

Write `crefctf(θ) = REF(φ_p−θ, k)·SH(φ_p−θ, k)·CTF(φ_p, k)` and
`cdiff = crefctf − PTCL`. Then `∂f/∂x = Σ_k w_k Σ_p 2·Re[ (∂crefctf/∂x)·conj(cdiff) ]`
for each parameter `x`.

**Shifts** — unchanged from today (`:1529-1533`):
`∂crefctf/∂s_x = i·argtransf(rp,k)·crefctf`.

**Rotation** — chain rule through `φ_p − θ`, so a leading minus:

```text
∂crefctf/∂θ = −[ ∂_φREF(φ_p−θ,k)·SH  +  REF·∂_φSH(φ_p−θ,k) ] · CTF
```

Both pieces are cheap:

- `∂_φREF` is the **band-limited angular derivative of the reference**,
  `IFFT[ i·m · ft_ref ]`. It depends only on the reference, so it is memoized
  **once** into a `dref_dphi(pftsz, nk, nrefs)` array alongside `pfts_refs_*` —
  one extra copy of the reference stack, built in `memoize_refs`
  (`simple_polarft_memo.f90:126`).
- `∂_φSH` is **closed form**, and expressible in the table that already exists.
  With `h = k·sin φ` and `kc = −k·cos φ`: `∂_φh = −kc` and `∂_φkc = h`, so the
  derivative of the phase argument is the *perpendicular* combination of the same
  polar coordinates:

```text
∂_φ arg = −argtransf(pftsz+rp, k)·s_x  +  argtransf(rp, k)·s_y
∂_φ SH  = i · (∂_φ arg) · SH
```

(assumes `A(1) = A(2)`, i.e. a square box — true in practice, but assert it.)

That the angular derivative is `perp` of the position vector is the same fact as
`∂loc/∂ω = [loc]×` in the 3D note, specialised to one rotational DOF. It is worth
recognising, because it is also a free sanity check on the sign convention.

**So the θ-gradient is one more accumulator in the existing loop, with one extra
memoized array — the same cost and shape as the shift gradient already there.**

## 4. Staged implementation

### Stage 0 — the baseline you must beat

Implement **parabolic interpolation of the residual minimum**: fit a parabola to
`r(j−1), r(j), r(j+1)` around the `minloc` index. Three lines, no new machinery,
and it recovers a large fraction of the pure sub-grid gain.

Everything below is measured against *this*, not against the raw grid. It also
separates the two effects cleanly: the parabola should capture the continuity
gain and **none** of the rotation/shift coupling gain (§4.2).

### Stage 1 — continuous angular minimum

Replace `maxloc` in `grad_shsrch_optimize_angle` with continuous minimisation of
the trig polynomial `r(θ)`.

**New in `simple_polarft_corr.f90`:**

```fortran
! gen_euclids with the coefficients retained instead of discarded
module subroutine gen_euclid_angular_coeffs( self, iref, iptcl, shift, rcoeffs, jmin )
! r(theta), r'(theta), r''(theta) at continuous theta in grid units
module subroutine eval_euclid_resid_at_angle( self, rcoeffs, theta, r, dr, d2r )
```

`gen_euclid_angular_coeffs` is `gen_euclids` with the `c2r` retained (for `jmin`)
and `crvec1%c` handed back.

**In `simple_pftc_shsrch_grad.f90`:**

- add `real :: cur_inpl_ang` alongside `cur_inpl_idx`;
- `grad_shsrch_optimize_angle` (`:165`): get `jmin`, then 2–3 Newton steps on
  `dr = 0` from `θ = jmin`, clamped to `|θ − jmin| ≤ 0.5` grid steps, with a
  fallback to `jmin` if `d2r ≤ 0`;
- keep `cur_inpl_idx = nint(cur_inpl_ang)` so every existing consumer still works.

**API ripple.** `grad_shsrch_minimize` returns integer `irot`; callers convert to
`e3` via `get_rot`. Add an **optional** `theta` output — callers that ignore it
behave exactly as today. The shift rotate-back at
`simple_pftc_shsrch_grad.f90:249` (`rotmat2d(get_rot(irot))`) should use the
continuous angle when available; the shift lives in the rotated frame, so this
one matters.

### Stage 2 — joint 3-parameter refinement

Stage 1 still alternates. **Rotation and shift are strongly correlated** — for
density at radius `r`, a small rotation is nearly a perpendicular shift of
`r·δθ` — so coordinate descent crawls along that valley and can stall short of
the joint optimum. Separate gain from sub-grid continuity, and probably larger.

Extend L-BFGS-B from `ndim=2` to `ndim=3` over `(s_x, s_y, θ)`:

- `simple_pftc_shsrch_grad.f90:81`: `specify('lbfgsb', 3, ...)`, third limit a
  bounded window (say ±2 grid steps) around the Stage-1 angle — this is a local
  polisher, not a searcher;
- drop `opt_callback`; θ is now a real variable;
- the three wrappers (`:113-163`) pass a 3-vector and receive a 3-gradient.

Two routes for the objective at continuous θ, both needing measurement:

- **2a — coefficient route.** Per shift value, three forward angular FFTs
  (`S·REF`, `argtransf_x·S·REF`, `argtransf_y·S·REF`) give three half-spectra;
  then `r`, `∂r/∂θ`, `∂r/∂s_x`, `∂r/∂s_y` all evaluate at continuous θ in
  `O(pftsz)`. **The `FT(CTF²)·FT(REF²)` term is shift-independent, so it is
  hoisted entirely out of the search loop** — a saving specific to this
  structure. The `argtransf` multiplication must happen on the polar samples
  *before* the angular FFT (post-FFT pointwise multiplication is not
  equivalent — it would be a convolution in angular frequency), and the
  shift-phase arguments flip sign across the Friedel mate, so the derivative
  products extend to the second half-circle as `−conjg(...)`, not `conjg(...)`.
- **2b — direct-sum route.** Rotate reference and `dref_dphi` to continuous θ
  (`IFFT[e^{-imθ}·ft_ref]`, likewise for the derivative) once per evaluation,
  then reuse `gen_euclid_residual_grad`'s loop with `∂_φSH` from §3.5 added.

Rough flop counts at `nrots = 256`: the existing direct sum is ~2.5k flops per
shell, an FFT of length 256 ~10k. So 2a is ~3 FFTs and 2b ~2, both roughly an
order of magnitude per evaluation. **But the comparison is less bad than it
looks:** the current alternating scheme already pays a full `gen_euclids` FFT set
at every `opt_callback`, so the right metric is FFTs per *outer* iteration, not
per objective evaluation. Measure that, not the microbenchmark.

### Stage 3 — wiring

Only after 1–2 are validated in isolation. Continuous angles must survive into
`e3`, into class-average assembly (`rotmat2d` in
`simple_classaverager_core.f90:174` already takes a real angle — that side is
clean), and into project metadata. Expect the integer-`irot` assumption in more
call sites than a grep suggests.

## 5. Validation

Each test gates the next. All cheap, no project I/O.

1. **Route identity (§3.3).** For every grid `irot`, `exp(−1 − r(irot)/A)` from
   the coefficient route must equal `gen_euclid_grad_for_rot_8`'s `f` to
   single-precision round-off. Catches the Nyquist weight, every index and
   conjugation slip, and the `w_k`/`denom` normalisation. **If this fails,
   nothing downstream means anything.**

2. **Analytic vs numerical gradient.** Central differences in θ at several
   off-grid angles and several shifts; `O(h²)` agreement. Repeat for the shift
   components through the new 3-vector interface, to confirm no regression in the
   existing shift gradient.

3. **Ground truth vs real-space rotation — the aliasing test.** Rotate an image
   in real space by a known non-grid angle, build its PFT, confirm the continuous
   residual minimises there. **Run twice: once with a reference low-passed well
   below the ring Nyquist, once hard-edged/full-band.** The gap between them *is*
   the aliasing error of §3.4, and it sets the credible precision of the method.
   Report it as a number, not pass/fail.

4. **Recovery on synthetic particles.** Known class average, known continuous
   `(θ, s_x, s_y)`, realistic CTF and per-particle `sigma2`; recover. RMS error
   per parameter for: grid-only (today), parabola (Stage 0), Stage 1, Stage 2.
   **This four-column table is the deliverable** — it is what says whether
   Stage 2 was worth building.

5. **Monotonicity.** Stage 2's residual at its solution must be `≤` Stage 1's at
   its solution, per particle, to L-BFGS-B tolerance. A violation is a sign or
   frame-convention error, not a tuning issue.

6. **Zero-width regression.** With the Stage-2 angular window clamped to zero,
   results must reproduce current behaviour exactly.

## 6. What to measure on real data

The angular error scales with radius, so **the gain is concentrated at the
periphery**. A global FRC will dilute it, possibly to invisibility. Use a
radially resolved or ring-wise correlation between independently refined
half-sets and report improvement as a function of radius. Curves separating at
large `r` and coinciding at small `r` confirms the mechanism; uniform separation
means something else is happening and is worth understanding before claiming the
win.

Secondary: cluster2D iterations to convergence. If rotation/shift coupling was
genuinely slowing the alternating scheme, Stage 2 should show up as faster
convergence, not only better averages.

## 7. Risks

- **Aliasing at critical sampling** (§3.4) — the main technical risk and the one
  the obvious tests miss. Test 3 exists for it.
- **The `pftsz` vs `pftsz+1` question** (§3.4.2) — resolve before building.
- **Overfitting.** Three continuous parameters per particle against a reference
  that particle helped build. The reference a particle is refined against must
  not contain it. Thinner margin in 2D than 3D, where a class may hold only a few
  hundred particles.
- **Scope creep into the global search.** Polar's structural advantage is that
  rotation is an index shift, so all `nrots` residuals come from one FFT. Nothing
  here competes with that and nothing here should try — this is a **local polish
  on the final assignment**, a handful of L-BFGS-B iterations per particle.
- **Frame convention on the shift.** `shmat` is indexed at the *rotated* angular
  index (`:1537`) and the caller rotates the result back (`:249`). Getting this
  wrong yields a subtly wrong answer that still converges. Test 3 catches it; no
  self-consistency test will.

## 8. Sequencing

Stage 0 + Stage 1 with tests 1–4 is a self-contained first deliverable worth
doing regardless of what follows: small, exact, cheap, and it produces the
four-column table that decides whether Stage 2 is justified. **Do not start
Stage 2 before that table exists.**

## 9. Objective-function coverage

Do `OBJFUN_EUCLID` first and completely.

- **`cc`** (`gen_corrs`, `:206`) has the same angular-FFT structure but produces
  a *ratio* `N(θ)/sqrt(D(θ)·…)`, so every derivative needs the quotient rule.
  Mechanical once euclid is proven. Note `D(θ)` is shift-independent, the same
  hoist as §4.2.

  *Status:* implemented as `gen_corr_grad_at_angle`
  (`simple_polarft_corr.f90`), matching the `gen_corrs` selection score
  (unweighted shells, Nyquist bin included — deliberately NOT the k-weighted
  correlation of the legacy discrete shift gradient
  `gen_corr_cc_grad_for_rot_8`). The `D(θ)` series is assembled from the
  memoized `FT(CTF²)·FT(REF²)` products without an extra FFT and recomputed
  per evaluation; hoisting it across L-BFGS-B iterations remains available
  if profiling justifies the added state. Validated by
  `simple_test_continuous_inplane_cc_grad` (on-grid identity vs `gen_corrs`
  at 2.3e-7, three-variable FD gradient checks at 1.1e-4). Wired into the
  joint continuous route: `new_joint`/`minimize_joint` dispatch on the
  objective (`is_joint_grad_objfun`), with cc-aware score mapping
  (`score = -loss` clamped to `[-1,1]` instead of `exp(-loss)`) and a
  `|cc| > 1` invalidity guard replacing the Euclidean negative-loss guard.
  The hybrid/denoised blend remains excluded everywhere.
- **Hybrid / denoised** (`gen_hybrid_grad_for_rot_8_local`, `:1352`;
  `gen_denoised_corrs`, `:502`) is a *linear* blend of raw euclid and a
  cc-shaped denoised term with weight `objfun_den_w`. Once each component has a
  continuous-θ form the blend adds nothing new — but it does mean the denoised
  term needs the `cc` work above, so it lands last.

Doing these in parallel with euclid will make any failure impossible to
localise.
> Historical design note. The experimental standalone angular-coefficient API
> described below was removed with the continuous callback route. Production
> `inpl_cont=yes` uses the fused joint coefficient evaluator documented in
> `doc/algorithms/continuous_inplane_refinement_abinitio2D.md`.
