# `reconstruct3D_pcg` workflow integration: implementation note

## Status

**Specification. No code written yet.**

Scope: what must be built so `reconstruct3D_pcg` can be validated inside
`abinitio3D`, `refine3D_auto` and `refine3D_multi`. Three mechanisms are
specified in full -- symmetry (section 2), fractional update / trailing
(section 3), memory and the single-object scheme (section 4). Remaining gaps
are listed in section 6 without design.

Companion documents: `doc/policies/reconstruct3D_pcg_policy.md` (contract of
the code that runs today; supersedes the deleted
`ctf_sigma_weighted_pcg_reconstruction.md`) and
`continuous_3D_refinement_on_pcg_operator.md` (pose refinement, and its section
5.3 for the NU-driven resolution limit, which is assumed here).

**Validated against the tree at commit `872db5bc5`.** The `rec_pcg refactoring`
commit `5858d023f` renamed the module and commander and replaced the operator
note with the policy document; all symbol and path references below were
re-checked after that.

## 1. Sufficient statistics

Everything in this note rests on one fact about the current implementation.

After `accumulate_absT2`, the entire particle-dependent state of the operator is
**two arrays**:

| Symbol | Content | Layout | Size at box 256 |
| --- | --- | --- | --- |
| `acc` | `sum_i G_i^dagger \|T_i\|^2` | full-range complex, `lims3` | 1.08 GB |
| `b` | `sum_i G_i^dagger conjg(T_i) y_i / sigma^2` | real, `box^3` | 67 MB |
| `acc_folded` | real part of `acc` on packed `cmat` layout | real, `cdim = (257,512,512)` | 269 MB |

`precond` and `Khat` are both **derived** from `acc_folded`
(`precond_from_accum`, `fold_accum_to_khat` + `finalize_khat`) and carry no
information not already in it. `finalize_khat` is FFT work only -- no particle
pass -- with one exception noted immediately below.

Persisted artifact per `(state, eo)` is therefore `(acc_folded, b)` =
**337 MB at box 256**, which is the same order as production's
`(halfmap, rho)` pair. Do not persist `precond` or `Khat`.

### 1.1 Required change: remove `calibrate_kernel`'s particle pass

`calibrate_kernel` calls `apply_normal_matrixfree(probe)` -- a full particle
pass -- purely to fit one scalar. Measured values:

| Dataset | Factor |
| --- | --- |
| synthetic hetero-CTF fixture | 6.43e1 |
| bgal, 5000 ptcls, box 256 | 6.398380e1 |

`padsc**2 = (padf**3)**2 = 64`. Two independent datasets land within 0.5% of
it. **The calibration is measuring an analytic constant.**

Action: set `Khat = Khat * self%padsc**2` directly; keep the measured fit only
inside `test=pcg_recon_kernel` as an assertion that it stays within tolerance
of 64.

Two reasons this is required, not optional:

- It removes roughly 10 s of the measured 26.9 s setup.
- **The trailing scheme in section 3 cannot calibrate at all**, because it
  derives `Khat` from an accumulator with no particles resident.

## 2. Symmetry

### 2.1 Settled: constraint and replication are the same estimator

Replication gives `H_sym x = Pi b` with
`H_sym = (1/M) sum_g S_g^dagger H S_g`; for symmetric `x`,
`H_sym x = Pi H Pi x`. Same normal equations, same solution.

Consequences, stated once and not revisited:

- **No SNR difference between the two.** The `M` replicas are one measurement
  written `M` times; no information is created. Symmetry's gain comes from
  cutting the number of unknowns by `~M` (amplitude SNR `~sqrt(M)`), and both
  formulations get it identically.
- Fewer unknowns means better conditioning, so expect fewer PCG iterations
  under symmetry, not more.

### 2.2 Decision: coordinate replication at scatter time

The policy note's section 10 leaves this open -- symmetry "must be handled
either by expanding the particle set over the asymmetric unit or by symmetrising
inside the operator; the choice interacts with the kernelized operator's
shift-invariance assumption and needs its own design pass". This section is that
design pass, and it picks the first option.

Scatter each plane pixel at all `M` orientations `S_g R_i` inside
`accumulate_absT2`, and symmetrize `b` the same way in `apply_adjoint_all`.
Rotation matrices from `sym%get_sym_rmat`; `sym%get_nsym` for `M`.

This reverses the earlier draft of this note, which recommended rotating
`Khat`/`rho` once at setup. Reasons for the reversal:

| | replicate coordinates | rotate accumulator |
| --- | --- | --- |
| extra memory | none | +1.08 GB (second full-range array) |
| cost, full batch (5000 ptcls, M=60) | ~5 min | ~1-2 min |
| cost, `update_frac=0.1` (500 ptcls) | **~30 s** | ~1-2 min, unchanged |
| extra interpolation error | none | yes, compounds with kernel's ~3% |
| matches production | yes | no |

The fractional-update regime of section 3 is the one that matters, and there
replication wins on time as well as on accuracy and memory, because its cost
scales with the subset while rotation is a fixed cost. `absT2_plane` is
evaluated once per particle regardless of `M`, so only the KB weight
computation and the 27-tap scatter are replicated.

### 2.3 Implementation contract

- `pcg_reconstruction` gains a `sym` handle set at `prep_particles` time, or a
  cached `symmats(3,3,M)` array alongside `rotmats`.
- `accumulate_absT2`: innermost loop over `g`, `loc = padf * [h,k,0] . (S_g R_i)`.
  The h-strided OpenMP colouring is computed from `h` only, so it remains valid
  -- but see 2.5.
- `apply_adjoint_all`: same replication for `b`.
- `apply_normal_matrixfree` must replicate identically, or it stops being the
  reference for the kernel. This makes it `M` times slower and effectively
  unusable for `pgrp /= c1` at large `M`; that is acceptable, but it means
  **symmetry makes `pcgop=kernel` mandatory in practice**.
- Nothing else changes. `precond_from_accum`, `finalize_khat` and `solve` are
  untouched: `H_sym` preserves the symmetric subspace by construction, so no
  per-iteration projection.

### 2.4 Do not

- Do not apply `Pi` to the iterate each PCG iteration. `M` volume rotations per
  iteration is unusable.
- Do not blend symmetric and asymmetric accumulators.
- Do not add soft symmetry (`lambda_sym * (I - Pi)`) before the hard constraint
  validates against production. It is a cheap follow-on -- `(I - Pi)` is a
  projector so the operator stays SPD -- but it is not part of this work.

### 2.5 Two things to verify

**Scatter race under replication.** The colouring stride
`ceiling(sqrt(3)*wdim)` guarantees that two `h`-lines in the same colour map
`>= padf*stride` apart *for a single rotation*. With `M` rotations scattered
inside the same worksharing iteration, two threads holding different `h` now
write `M` footprints each, and separation must hold for every pair `(g, g')`,
not just `g = g'`. It does not in general. Options: hoist the `g` loop outside
the `!$omp do` (one colour sweep per symmetry element), or serialize `g`.
**Hoisting the `g` loop outside the colour sweep is the safe default.**

**Axial views.** Particles whose central section is invariant under a subgroup
produce coincident replicas, which `rho` counts as independent. This affects
production too and is not introduced by this work, but it is worth a synthetic
test with an axis-concentrated orientation distribution. Open, section 9.

## 3. Fractional update and trailing reconstruction

### 3.1 Production contract, for reference

- Selection: `sample4update_class` / `sample4update_cnt` / `sample4update_all` /
  `sample4update_fillin` in `simple_matcher_smpl_and_lplims.f90`; probabilistic
  path samples once in `prob_align` and reproduces with `sample4update_reprod`.
  `sampled` marks the current round, `updatecnt` the persistent history.
- Realized fraction: `oris%get_update_frac` =
  `count(sampled == sampled_max .and. state > 0) / count(updatecnt > 0 .and. state > 0)`;
  per-state via `get_state_update_fracs`. Overridden by `ufrac_trec` when
  defined (single-state only).
- Blend, in `simple_commanders_rec_distr.f90::trail_restored_halves_if_needed`,
  applied to **restored** halfmaps after `sampl_dens_correct_eos`:
  `vol = uf*vol_current + (1-uf)*vol_prev`. Skipped when `uf >= 0.99`.
- Previous-artifact compatibility: `reconstructor_eo%read_eos_parallel_io` pads
  smaller previous dims when `l_update_frac` is active, rejects larger.

**The selection policy is not this note's concern and must not be
reimplemented.** Consume `pinds` from the existing samplers and the realized
fraction from `get_update_frac` / `get_state_update_fracs`.

### 3.2 Decision: blend the accumulators, not the volumes

```text
acc_t = uf * acc_new + (1 - uf) * acc_{t-1}
b_t   = uf * b_new   + (1 - uf) * b_{t-1}
x_t   = solve( acc_t, b_t )          starting from x_{t-1}
```

where `acc_new`, `b_new` are accumulated over the sampled subset only.

This is recursive weighted least squares with exponential forgetting: the
system solved at step `t` is the exact normal equations of a problem in which
each particle's contribution decays geometrically with the age of its last
update. Effective sample size is `~ subset/uf`, i.e. the full set.

Do **not** mirror production by blending solved volumes. Two reasons:

- A solve from `uf = 0.1` of the particles is badly underdetermined -- large
  null space, preconditioner floor doing most of the work. Blending two
  solutions of two different systems is not the solution of any system.
  Blending accumulators never constructs an underdetermined system in the first
  place.
- `(acc, b)` is the **complete sufficient statistic** (section 1). Production
  blends volumes because a gridding halfmap is the only summary it has; PCG has
  the actual statistic and should carry that instead. This is the one place the
  iterative formulation is structurally better suited to the existing workflow
  policy rather than merely equivalent to it.

Note this differs from production in *what* is blended and *when*: production
blends after density correction and filtering, this blends before the solve.
Consequence: PCG trailing is linear and exact, production's is an
approximation. That is a divergence to sign off deliberately, not a bug.

### 3.3 Warm start

`solve` already accepts nonzero `x` (`if( all(x == 0.0) )` branch). Pass
`x_{t-1}`. Between outer iterations the system changes by `uf` in the
accumulators, so the previous solution is close and 2-3 PCG iterations should
suffice against 10-15 cold.

This is the mechanism that makes the per-iteration cost acceptable. Without it
the solve is a fixed ~25 s regardless of `update_frac`, which would defeat the
purpose of fractional updates entirely.

Estimated per outer iteration at `update_frac = 0.1`, 5000 particles, box 256:

| Phase | Cold, full batch | Fractional + warm start |
| --- | --- | --- |
| I/O + prep | ~8 s | ~1 s |
| accumulate (`acc`, `b`) | ~7 s | ~1 s |
| `finalize_khat` (FFTs) | ~2 s | ~2 s |
| `calibrate_kernel` | ~10 s | **removed** (section 1.1) |
| blend | -- | <1 s |
| solve | ~25 s (15 it) | ~5 s (3 it) |
| **total** | **~52 s** | **~10 s** |

Estimates, not measurements. `finalize_khat` is the only phase that does not
shrink with `update_frac` and will dominate at small subsets.

### 3.4 Artifact contract

- Previous artifact per `(state, eo)`: `acc_folded` (real, `cdim`) and `b`
  (real, `box^3`). Written by the producing stage, required when
  `trail_rec=yes`, exactly as production requires previous halfmaps.
- **Do not derive the previous artifact from a volume.** The contract is
  accumulator-to-accumulator.
- Dimension compatibility: mirror `read_eos_parallel_io`'s rule -- pad smaller
  previous dims with zero when `l_update_frac` is active, reject larger. Note
  `acc_folded` is on the **padded** lattice (`cdim` from `boxpd`), so a box
  change rescales both dimensions and content; padding a Fourier-domain
  accumulator across a box change is not the same operation as padding a
  real-space halfmap and needs its own test.
- `x_{t-1}` must also persist for the warm start. It is the output halfmap, so
  no new artifact -- but it must be the **unfiltered** solve output, not a
  postprocessed map.

### 3.5 Interaction with the resolution limit

`acc` depends on `sqlp` through `absT2_plane`'s masking. When the NU-driven
limit advances between outer iterations (companion note section 5.3), the
previous `acc` was accumulated under a *lower* limit than the current one.
Blending them mixes two band limits.

This is benign in one direction only: the previous accumulator is zero above
its limit, so the blend under-weights the newly admitted shells by `(1-uf)`
until they refresh. It self-corrects over `~1/uf` iterations. Production has the
same property. **Record it, do not attempt to correct it**, but note it as a
reason to prefer a limit that advances slowly relative to `1/uf`.

### 3.6 Gates

- Synthetic: `uf = 1` must reproduce the full-batch solve bit-for-bit.
- Synthetic: a stationary particle set under `uf < 1` must converge to the same
  answer as full batch as `t` grows.
- Real data: FSC against a full-batch reference at matched total compute.

## 4. Memory

### 4.1 Production scheme

`simple_matcher_3Drec.f90::calc_3Drec` holds **one** `type(reconstructor)`.
Particles are pre-grouped by `(state, eo)` into a permutation `grouped_pinds`
with CSR offsets `state_eo_offsets` (`group_pinds_by_state_eo`,
`state_eo_group = 2*(state-1) + eo + 1`). The loop is:

```text
do state; do eo
   init_state_half_rec(recvol)
   batches -> update_state_half_rec
   write_state_half_partial(recvol, state, eo)
   kill_state_half_rec(recvol)
```

Peak memory is one half of one state, independent of `nstates`. Note it uses
plain `reconstructor`, not `reconstructor_eo` -- the eo split comes from the
grouping and the loop, not from holding both halves.

### 4.2 PCG mirror

Same structure, same grouping helpers, one `type(pcg_reconstruction)` object:

```text
do state; do eo
   reset accumulators (NOT kill -- keep wimg and its FFTW plans)
   batches -> accumulate_absT2 into acc, apply_adjoint_all into b
   write (acc_folded, b) partial for (state, eo)
   free acc
```

Then, separately and serialized over `(state, eo)`:

```text
do state; do eo
   read + reduce partials -> acc_folded, b
   blend with previous artifact (section 3.2)
   precond_from_accum; fold_accum_to_khat; finalize_khat
   solve (warm start from previous x)
   write halfmap; free precond, Khat
```

**Rule: `kill` the `pcg_reconstruction` object never, `reset` it per
`(state, eo)`.** `wimg` is state-independent (same `boxpd`), and killing it
destroys and rebuilds four FFTW plans per group.

Splitting accumulate from solve in time is what keeps peak memory at one
`(state, eo)` and is also the natural distribution boundary (4.4).

### 4.3 Budget at box 256

| Array | Size | Phase |
| --- | --- | --- |
| `acc` full-range complex | 1.08 GB | accumulate only |
| `wimg` padded image | ~0.54 GB | both |
| `Khat` real packed | 269 MB | solve |
| `precond` real packed | 269 MB | solve |
| `get_cmat()` copy | 539 MB | transient, per call |
| CG vectors `b,r,p,hp,z,x,rtr` | 7 x 67 MB = 470 MB | solve |

Accumulate peak ~2.7 GB, solve peak ~2.1 GB. Matches the measured ~3 GB.

Two reduction opportunities, both optional:

- **`acc` as packed real instead of full-range complex** would save ~800 MB and
  is the single largest win. **A previous attempt at a packed accumulator was
  implemented and reverted**: production's `expand_cmat` materializes the
  periodic extension up front so its gather and scatter are consistent by
  construction, whereas this operator wraps on the fly, and a non-wrapping
  scatter breaks the exact self-adjointness CG requires. Any retry must
  reproduce the adjoint dot-product test first.
- `get_cmat()` returns a copy; `image` has `get_rmat_ptr` but no cmat
  equivalent. Adding one removes 539 MB of transient and measurable memcpy
  from both `apply_normal` and `apply_precond`.

### 4.4 Distributed split

`acc` and `b` are sums over particles, so they partition across parts and
reduce exactly like gridding's partial reconstructions. The solve does not
partition -- and does not need to: once reduced, it is particle-independent and
runs on one node with a serial tail that does not grow with `nptcls`.

Per-part memory is the accumulate peak (~2.7 GB at box 256), one `(state, eo)`
at a time, matching production's per-part footprint.

## 5. Particle I/O and preprocessing

### 5.1 `y_planes` must go -- this is the blocking defect

`simple_commanders_reconstruct3D_pcg.f90` allocates
`y_planes(lims2, lims2, nptcls)` and holds every particle's Fourier plane
resident for the whole run:

| `nptcls` | `y_planes` |
| --- | --- |
| 5 000 | 2.5 GB |
| 50 000 | 24.6 GB |
| 500 000 | 246 GB |

Two facts make the fix structural rather than a compromise:

- **`y_planes` is read in exactly one place**: `apply_adjoint_all`, called once
  at the top of `solve` to build `b`. The CG loop never touches it again.
- **`accumulate_absT2` needs no image data at all** -- only `rotmats`,
  `ctfparms` and `sig2`, all cached scalars.

So the images are needed for one thing only, and for one pass only. Restructure
to accumulate `b` and `acc` inside the batch loop and discard each batch:

```text
prepimgbatch(params, build, MAXIMGBATCHSZ)          ! MAXIMGBATCHSZ = 500
do ibatch
   discrete_read_imgbatch(...)
   per particle: norm_noise(lmsk); taper_edges_particle; fft; extract_native_plane
   accumulate this batch into b        (needs the plane)
   accumulate this batch into acc      (needs only cached scalars)
end do
killimgbatch(build)
solve(b, x, ...)                                     ! no y_planes argument
```

`solve` loses its `y_planes` argument and takes a prebuilt `b`. Particle-image
residency drops from `nptcls` planes to `MAXIMGBATCHSZ` planes -- 500 x 528 kB
= 264 MB at box 256, independent of dataset size.

This is a prerequisite for everything else in this note at real dataset sizes,
and it is the natural place to introduce the accumulate/solve split that
sections 3 and 4 already require.

### 5.2 What stays resident

| Item | Per particle | 500 k particles |
| --- | --- | --- |
| `rotmats(3,3)` | 36 B | 18 MB |
| `ctfparms` | ~40 B | 20 MB |
| `shifts(2)` | 8 B | 4 MB |
| `sig2(0:R)` | 516 B | 258 MB |
| **total** | **~600 B** | **~286 MB** |

`sig2` dominates and is the only one worth revisiting: it is stored per
particle at full shell resolution, but `euclid_sigma2` groups particles, so
storing a group index plus a per-group spectrum would cut it by the group
count. Not urgent; note it before scaling past ~10^5 particles.

### 5.3 Disk caching of preprocessed planes: do not

**The reason is not that gridding does not cache.** Gridding inserts a particle
once per iteration and is done; PCG iterates, so the analogy needs justifying
rather than asserting. The justification is specific to this formulation.

**PCG on the normal equations never compares against the particle data.** `y`
enters exactly once, when `b = A^dagger y` is formed. In `solve`, `y_planes` is
read on one line and never again:

```fortran
b   = self%apply_adjoint_all(y_planes)   ! the only read of y
hp  = self%apply_normal(p)               ! rotmats, absT2, Khat -- no y
z   = self%apply_precond(r)              ! precond -- no y
rtr = b - self%apply_normal(x)           ! residual replacement: b, not y
```

`H = A^dagger A` depends on the **geometry and weights** -- orientations, CTF,
sigma^2 -- not on the measured intensities. Forming the normal equations is the
step that discards the data, keeping only its projection `b`. This is why
residual replacement (companion note, `solve`) recomputes `b - Hx` rather than
`A^dagger(Ax - y)`.

Consequence: per outer iteration the I/O is *read the sampled subset once, form
`acc` and `b`, discard* -- the same particle traffic as gridding, not
`n_iterations` times more. There is nothing for a cache to amortize.

This property is also load-bearing for section 3: the trailing scheme blends
`(acc, b)` and solves without revisiting old particles, which is only possible
because the solve is data-free.

Two supporting points:

- **Cache entries are not smaller than what they replace.** At box 256 a raw
  image is 256 kB and a packed preprocessed plane is 259 kB, so a cache trades
  disk footprint for CPU (`norm_noise` + `taper_edges_particle` + `fft`), not
  for I/O.
- **Invalidation has no owner.** `box_crop`, `mskdiam` (via `lmsk`),
  `sigma_est`, `ctfflag` and the particle selection all change what a cached
  plane should contain, several of them per iteration.

### 5.3.1 Where this argument stops holding

The above is a property of the normal-equations formulation, not of iterative
reconstruction in general. Two cases invert it, and both should be revisited
deliberately rather than discovered:

**Continuous pose refinement.** The pose objective is
`r = C S V_hat(loc) - y` evaluated per Levenberg--Marquardt step, so `y` is
compared repeatedly. Per outer iteration rather than per CG iteration, but a
genuine repeated comparison. See
`continuous_3D_refinement_on_pcg_operator.md`; caching should be re-evaluated
there on its own terms.

**A move off the normal equations.** CGLS or LSQR applied to `A` directly needs
`A p` and `A^dagger r` every iteration:

| | normal equations (current) | CGLS / LSQR |
| --- | --- | --- |
| conditioning | `cond(A)^2` | `cond(A)` |
| data touched | once, forming `b` | every iteration |
| cost per iteration | 1.5 s (kernel), independent of `nptcls` | ~2 particle passes, scales with `nptcls` |
| trailing accumulators (section 3) | works | impossible |
| plane caching | pointless | essential |

Better conditioning is a real attraction given CTF zeros, but it costs roughly
10x per iteration, removes the kernel's particle-count independence, and breaks
the trailing design. **Trigger for revisiting:** the residual stalling at a floor
traceable to squared conditioning rather than to the operator or the
preconditioner. Nothing observed so far suggests this.

**The one case worth measuring before dismissing:** when
`box_crop << box`, a plane cached at the cropped size is much smaller than the
raw read -- at `box_crop = 64` against `box = 256` it is roughly 4.5 kB versus
256 kB. Early workflow stages run both the smallest boxes and the most
iterations. If profiling shows read I/O dominating those stages, a cropped-plane
cache is worth revisiting; do not build it speculatively. Open, section 9.

### 5.4 The real preprocessing cost is `absT2`, and the fix is compute

`absT2_plane` evaluates, per pixel per particle:

```fortran
spafreqsq = (real(h)/real(self%box))**2 + (real(k)/real(self%box))**2
ang       = atan2(real(k), real(h))
```

Both depend only on `(h,k)` -- **not on the particle** -- yet both are
recomputed inside the per-particle loop. That is `3.3e8` `atan2` plus `sqrt`
calls per accumulate pass at 5 000 particles, and `3.3e10` at 500 000.

Fix: hoist them into two `(h,k)`-indexed lookup tables shared by all particles.
258 kB each, allocated once, **not** per particle. Additionally adopt
production's kernel form: `gen_fplane4rec` precomputes `shell_lut` and
`sum_df / diff_df / angast / amp_contr_const / wl / half_wl2_cs` per particle,
then evaluates `ft_map_ctf_kernel` with no per-pixel transcendentals.

Note what this settles: the CTF part `C_i^2` is genuinely invariant across
outer iterations (defocus does not change during refinement), so it looks like a
caching candidate. It is not -- storing it costs 258 kB per particle, the same
problem as caching planes. A shared LUT plus a cheap kernel is strictly better
than caching a per-particle array. `1/sigma^2` is a per-shell scalar and does
change per iteration, so it must stay outside anything cached anyway.

### 5.5 Rules

- One I/O pass per outer iteration. Never two, and never one per PCG iteration
  -- the solve is data-free (5.3).
- Never hold more than `MAXIMGBATCHSZ` particle images or planes.
- No new on-disk artifact per particle. The only per-particle persistence is
  what the project already carries.
- Match `calc_3Drec`'s read -> prep -> consume -> discard shape, including its
  `(state, eo)` grouping (section 4.2), so the batch loop serves both.

## 6. Remaining gaps

| Gap | Status | Blocks |
| --- | --- | --- |
| all particles resident (`y_planes`) | see 5.1 | any real dataset size |
| even/odd half sets + FSC | absent | everything; no measurement possible |
| symmetry | **hard-errored**, `pgrp=c1` enforced | `refine3D` on symmetric data |
| distributed | **hard-errored**, `nparts=1` enforced | real dataset sizes |
| `box_crop` / `smpd_crop` | uses `params%box`/`smpd` | like-for-like comparison |
| driven resolution limit | `sqlp` fixed at native Nyquist | workflow integration |
| particle weights | `prep_particles` ignores `w` | frac-update workflows |
| multi-state | single state | `refine3D_multi` |
| project write-back | never writes | stage handoff; `set_bp_range3D` reads FSC files |
| postprocess handoff | none | double filtering with NU path |

None of these need design work; they are ports of existing contracts.

## 7. Sequencing

0. Remove `y_planes` and fuse the accumulate into the batch loop (5.1).
   Blocking defect; nothing runs at real dataset size until this is done.
1. Remove `calibrate_kernel`'s particle pass (1.1), and hoist the `absT2`
   `(h,k)` lookup tables (5.4). Prerequisite for 3, and both pay for themselves
   immediately.
2. Even/odd + FSC. Nothing is measurable before this.
3. `box_crop` + driven `sqlp`.
4. Symmetry (2). First item with scientific content, and now measurable.
5. Particle weights.
6. Fractional update + trailing (3), including the accumulate/solve split.
7. Multi-state + single-object memory scheme (4.2), distributed reduction (4.4).
8. Write-back and postprocess handoff.

Steps 1-4 are validatable in isolation with `reconstruct3D_pcg`. Steps 5-8 are
what make it substitutable inside `refine3D`.

## 8. Acceptance gates

- `uf = 1` reproduces full batch bit-for-bit.
- `pgrp` symmetric run agrees with production `reconstruct3D` at the same
  `pgrp` to within the kernel's measured error.
- Half-set FSC no worse than gridding at matched compute.
- Peak RSS per part no worse than production at the same `nstates`.
- At least one regime measurably better. Current candidates: heterogeneous CTF,
  low particle counts, axial-view-heavy symmetric data.

If only equivalence is achieved, the correct conclusion is that this is a
better-founded implementation of the same estimator and the value is downstream
(priors, soft symmetry, continuous refinement), not in substitution.

## 9. Open questions

1. Does the axial-view multiplicity (2.5) matter on real symmetric data?
2. Does warm-started PCG actually converge in 2-3 iterations under trailing, or
   does the changing operator break the warm start's usefulness? The whole cost
   argument in 3.3 depends on this and it is unmeasured.
3. Is the kernel adequate at the small `box_crop` values early workflow stages
   use? Its ~3% error was measured at box 256 only.
4. How should `lambda` and `RHO_FLOOR_FRAC` scale with `uf` and with `M`? Both
   are calibrated against the data term's absolute scale, which both mechanisms
   change.
5. Padding a Fourier-domain accumulator across a `box_crop` change (3.4) -- is
   zero-padding the right operation, or does it need rescaling?
6. Does read I/O dominate the small-`box_crop` early stages? Only that would
   justify revisiting the cropped-plane cache dismissed in 5.3. Profile before
   building anything.
7. `sig2` at 258 MB per 500 k particles (5.2) -- worth replacing with a group
   index plus per-group spectra, or is it below the noise floor?
8. Is squared conditioning (5.3.1) ever the binding constraint? If preconditioned
   CG stalls at a floor not explained by the operator or preconditioner, the
   CGLS/LSQR fork is on the table -- and it inverts every I/O conclusion in
   section 5.

## References

- `src/main/strategies/search/simple_matcher_3Drec.f90` -- `calc_3Drec`,
  `group_pinds_by_state_eo`, single-reconstructor loop.
- `src/main/commanders/simple/simple_commanders_rec_distr.f90` --
  `trail_restored_halves_if_needed`, `determine_trailing_update_fraction`.
- `src/main/ori/simple_oris_getters.f90` -- `get_update_frac`,
  `get_state_update_fracs`.
- `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` -- samplers,
  `set_bp_range3D`.
- `src/main/volume/simple_reconstructor.f90` -- `insert_plane_oversamp`
  symmetry replication, `sampl_dens_correct`.
- `src/main/volume/simple_reconstructor_eo.f90` -- `read_eos_parallel_io`,
  `set_sh_lim`.
- `src/main/simple_sym.f90` -- `get_nsym`, `get_sym_rmat`.
- `.claude/skills/simple-frac-update-trailing/` -- contract this note conforms
  to.
