# `reconstruct3D_pcg` workflow integration: implementation note

## Status

**Specification. No code written yet.**

Scope: what must be built so `reconstruct3D_pcg` can be validated inside
`refine3D_auto` and `refine3D_multi`.

**`abinitio3D` and the early workflow stages are out of scope and stay on
gridding.** While sampling there is sticky -- `sample_ptcls4update` reproduces
the initial subset rather than redrawing it while `startit == 1` -- the
accumulator carries no new particles between iterations, so the trailing
formulation of section 3 has nothing to amortize and gridding is the right tool.
This is a deliberate boundary, not a deferral: it removes the smallest
`box_crop` values and the highest iteration counts from the target regime. Three mechanisms are
specified in full -- symmetry (section 2), fractional update / trailing
(section 3), memory and the single-object scheme (section 4). Remaining gaps
are listed in section 7 without design.

Companion documents: `doc/policies/reconstruct3D_pcg_policy.md` (contract of
the code that runs today; supersedes the deleted
`ctf_sigma_weighted_pcg_reconstruction.md`) and
`continuous_3D_refinement_on_pcg_operator.md` (pose refinement, and its section
5.3 for the NU-driven resolution limit, which is assumed here).

**Validated against the tree at commit `9839d8639` (master).** The `rec_pcg refactoring`
commit `5858d023f` renamed the module and commander and replaced the operator
note with the policy document; all symbol and path references below were
re-checked after that.

## Where this stands (end of session, uncommitted on top of `9839d8639`)

Phase A item 1 is **built and validated on real data**. Nothing is committed;
eight files are modified in the working tree.

### Done and verified

| | |
| --- | --- |
| streaming accumulation (5.1) | `begin_accum` / `accumulate_batch` / `end_accum` / `solve_accum`; `y_planes` gone from the command |
| `calibrate_kernel`'s particle pass removed (1.1) | now `Khat * padsc**2`; the fit survives as `measure_kernel_scale` for the test only |
| rim pre-filter | `sq_rim` rejects on `h^2+k^2` before computing anything |
| four tests merged into `test=pcg_recon` | eight fail-fast stages, including two that did not exist before |

Measured on 5000 particles, box 256, `mskdiam=220`, `nthr=10`, `pcgop=kernel`:

| | before | after |
| --- | --- | --- |
| setup | 38.0 s | **24.7 s** |
| solve | 13.4 s | 13.4 s (1.6 s/iteration) |
| total | 51.4 s | **38.4 s** (~1.8x gridding's ~21 s) |
| iterations | 8, `rtol=1e-3` **reached** | same |
| peak memory | ~5.5 GB | ~3.0 GB, and now constant in `nptcls` |

Two consecutive runs gave identical results. That was the actual point of the
race fix below, and it is now confirmed on real data rather than only at box 24.

### The significant find: a race in the h-strided colouring

**Pre-existing, not introduced by this work, and it affected every result this
operator has ever produced.** The colouring guarantees two same-colour `h`-lines
map at least `padf*stride` apart -- but that separation is computed on
UNWRAPPED coordinates. The accumulator is periodic, so a coordinate past
`wlims(2)` folds to `wlims(1)` and two points most of a period apart can land
within a voxel of each other. Two threads in the same colour sweep then write
the same voxel.

At box 24, `|loc| <= 24` while `wlims = [-24,23]`, and `h = -12` and `h = +12`
differ by `4*stride` -- same colour, overlapping windows after folding.

It was present at **three** scatter sites: `scatter_plane` (RHS),
`accumulate_absT2` (preconditioner and kernel), and `apply_normal_matrixfree`
(the operator, every CG iteration). All three now split the plane: interior in
the parallel coloured pass, wrapping rim in a serial pass, gated by `win_wraps`.
Serial rather than `!$omp atomic` deliberately -- an atomic removes the race but
leaves summation order varying, so results still would not reproduce.

Consequence worth remembering: **maps produced before this fix were not
reproducible at the Nyquist rim.** The effect is small and confined to the
outermost shell, but it was there.

It was found because the merged test compares a batched and a monolithic
accumulation of identical data. Nothing before it had ever compared two runs
whose thread scheduling differs. The single-thread pass was the diagnostic that
identified it as a race in one run.

### Not yet verified -- do this first next session

`win_wraps` was made module-level (matching `gather_window` / `scatter_window`,
which were moved out of type-binding for exactly this reason) on the hypothesis
that ~2.6e8 dispatched calls per accumulation account for the ~8 s of the rim
regression that the `sq_rim` pre-filter did not recover. **This is unbuilt and
untested.**

Expected: the residual trace must stay bitwise identical
(`1.0862E-01 ... 9.6795E-04`) since only the call form changed; setup should
fall from 24.7 s toward ~17 s. If setup barely moves, the hypothesis is wrong
and setup should be instrumented per phase rather than guessed at a third time.

### Coverage gaps to keep in mind

- The command's own I/O loop is not covered by any test: `discrete_read_imgbatch`,
  `norm_noise`, `taper_edges_particle`, `extract_native_plane`, and the
  `y_batch` indexing. Only a real project exercises them.
- In particular the **short final batch** never fires at 5000 particles, since
  `MAXIMGBATCHSZ = 500` divides it. A particle count that is not a multiple of
  500 would test it.
- A uniform scaling of `H` leaves the CG residual sequence exactly invariant
  (`alpha -> alpha/c` cancels in `r <- r - alpha*Hp`), so an identical residual
  trace does **not** validate a change to the kernel scale -- only the map
  amplitude moves. Do not use the trace as evidence for that class of change.

### Next

Re-run box-256 setup timing to quantify the `ft_map_ctf_kernel` form now adopted
in `absT2_plane` / `build_transfer` (5.4). Per-phase profiling (added to the
commander) already localized the cost: on a box-256, 3000-particle `objfun=cc`
run, setup 23.6 s splits as accumulation 17.4 s (of which **`scatter |T|^2`
15.1 s**, read I/O 0.17 s, prep 1.8 s), `end_accum` 4.2 s, operator alloc 1.7 s
-- so OQ2 is answered: the scatter is the whole story. The post-CTF-form re-run
is now in: `scatter |T|^2` = 15.9 s vs the 15.1 s pre-form baseline (setup
25.1 vs 23.6 s), i.e. **no measurable win** -- within run-to-run noise. The CTF
form removes call/arithmetic overhead but that was never the bottleneck: the
15.9 s is the KB-window lattice scatter (`apod_mat_3d_fast` + 27-tap
`scatter_window`), which 2b does not touch. That is now firmly the next target.
Phase A item 2 is complete. Then symmetry (section 2), phase A item 3 and the
first with scientific content.

First bit-identical pass at the scatter is in: the interior (non-wrapping)
colour sweeps of `accumulate_absT2` and `scatter_plane` now call
`scatter_window_nowrap` / `scatter_window_cmplx_nowrap`, which index the
accumulator directly instead of through `self%wrap`. `win_wraps` has already
proven the window cannot wrap there, so `self%wrap` is the identity over its
whole span -- the direct form is bit-identical and makes the inner run
contiguous (stride-1, vectorizable); the `|T|^2` variant also applies its real
weight at half the FMAs of the complex path. The wrapping rim keeps the
`self%wrap` scatter, and the matrix-free reference path is untouched. Measured:
**no win.** `scatter |T|^2` = 15.1 s, identical to the 15.1 s baseline (the
three runs -- pre-CTF, post-CTF, post-nowrap -- are all 15.1-15.9 s, one noise
band). The cost is cache-miss traffic writing into the ~1 GB accumulator at
rotated locations, not index arithmetic; removing the `wrap` indirection cannot
move a memory-bandwidth wall. Reducing it would need a structural change
(cache-blocked scatter, a smaller working set, or a different accumulation
order) -- out of scope for phase A. The nowrap scatter stays because it is
bit-identical (`test=pcg_recon` green, deapod-OFF corr `0.97155` unchanged) and
strictly cheaper in instruction count; the scatter cost is "it is what it is".

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

### 1.1 Required change: remove `calibrate_kernel`'s particle pass -- DONE

`calibrate_kernel` calls `apply_normal_matrixfree(probe)` -- a full particle
pass -- purely to fit one scalar. Measured values:

| Dataset | Factor |
| --- | --- |
| synthetic hetero-CTF fixture | 6.43e1 |
| bgal, 5000 ptcls, box 256 | 6.398380e1 |

`padsc**2 = (padf**3)**2 = 64`. Two independent datasets land within 0.5% of
it. **The calibration is measuring an analytic constant.**

Action: set `Khat = Khat * self%padsc**2` directly; keep the measured fit only
inside `test=pcg_recon` stage 5 as an assertion that it stays within tolerance
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

### 2.2 Symmetry of the output, and lattice-exact point groups

Under coordinate replication the output is symmetric **by construction**, not
approximately: `acc` and `b` are `G`-invariant, `H_sym` preserves the symmetric
subspace, and CG started from `x = 0` never leaves it. Same mechanism and same
guarantee as gridding. (The "constraint does not guarantee symmetry" concern is
real, but it applies to the *soft* penalty of 2.6, not to the hard constraint.)

The residual question is discretization, and the answer splits by point group.
A symmetry operator maps the sampling lattice onto itself exactly **iff its
matrix is a signed permutation matrix** -- entries in `{0, +-1}`, one non-zero
per row and column. Then `(h,k,m)` maps to `(+-h', +-k', +-m')`, an exact
relabelling of lattice sites; and because the KB stencil is a separable product
of three identical 1-D windows, permuting axes maps the stencil onto itself.
Replication is then exact.

- **`d2` qualifies.** Its operators are `{I, 180deg about x, y, z}`, e.g.
  `180deg_z: (h,k,m) -> (-h,-k,m)`. A `d2` map from this operator is
  **perfectly** symmetric, up to floating-point summation order -- not
  interpolation-limited.
- `c2` about a Cartesian axis, `c4`, `d4`, `t` and `o` also qualify: all their
  operators are signed permutations, including the 3-folds about body diagonals,
  `(h,k,m) -> (m,h,k)`.
- `c3` about z, `c5`, `c6`, `d3`, `d5`, `d6` and icosahedral do **not**. Those
  operators need interpolation, and the output is symmetric only to
  interpolation accuracy.

Gridding replicates the same way and inherits exactly the same split, so neither
method is better here. The difference is that it is now written down rather than
assumed.

**Implementation consequence.** For a lattice-exact group, replication can be
done as an index permutation of `acc` -- one pass per operator, no particle
loop, no KB weights, and exact:

```text
acc_sym = sum_g  acc . S_g^-1        (exact index remap)
```

At `d2` with 5000 particles that is three passes over `acc` (order seconds)
rather than 4x the scatter (order 20 s), and it *removes* interpolation error
instead of adding it. Detect the case at runtime by testing each
`sym%get_sym_rmat` result for the signed-permutation property; do not
special-case `pgrp` strings.

**Gate:** on a `d2` fixture the permutation path must agree with coordinate
replication to round-off. If it does not, the reasoning above is wrong somewhere
and coordinate replication is the fallback.

### 2.3 Decision: coordinate replication at scatter time

**DONE.** Coordinate replication is implemented in `reconstructor_pcg`
(`symmats(3,3,nsym)` cached at `set_sym` time; `g`-loop replication in
`accumulate_absT2`, `scatter_plane`/`apply_adjoint_all`, and
`apply_normal_matrixfree`) and wired into `reconstruct3D_pcg` (the `pgrp/=c1`
hard-error is lifted; `set_sym(build%pgrpsyms)` is called after `pcgop%new`).
`nsym=1` (c1) is bit-identical to the pre-symmetry path because `symmats(:,:,1)`
is exactly the identity. Validated by `test=pcg_recon` stage 9: an in-operator
c2 build produces the same normal system -- kernel `Khat` and RHS `b` -- as a c1
build of the c2-expanded particle set, to accumulator round-off
(`rel_err(Khat) = 1.3e-6`, `rel_err(b) = 1.3e-7`). (Compared at the operator
level, not through a solve: the kernelized operator is only approximately SPD,
so a zero-tolerance solve driven to convergence loses positive-definiteness --
that is stages 5-7's concern, not this stage's.) The lattice-exact permutation
optimization (2.2) and the asymmetry diagnostic (2.4) remain future work.

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

### 2.4 Search-space coupling: the reconstruction is mode-agnostic

`sym%build_refspiral` gates projection directions on
`within_asymunit(o, incl_mirr=.false.)`, so for `pgrp /= c1` production's search
grid covers only the asymmetric unit. That policy assumes the reference is
symmetric **and** aligned to the point group's principal axis.

**Neither assumption is introduced by the constraint formulation, nor weakened
by it.** `get_sym_rmat` returns operators in the point group's fixed frame, so
replication -- in gridding or here -- assumes the map is aligned to that frame;
misalignment destroys signal in both.

What does change is that the reconstruction no longer cares which mode the
search ran in. If the volume is `G`-symmetric, a particle at `R` and the same
particle at `S_g R` produce identical projections, so whichever ASU
representative the search assigns, replicating over all `M` operators yields the
*same* accumulator:

| search | cost | assumes symmetric | assumes aligned | reconstruction |
| --- | --- | --- | --- | --- |
| ASU-restricted (production) | 1x | yes | yes | identical |
| full sphere | `M`x | no | no | identical |

A full-sphere mode is therefore a **search** change, not a reconstruction
change, and belongs in the `refine3D` integration rather than in
`reconstruct3D_pcg`. It is affordable at low symmetry (`d2`: 4x) and not at high
(`i`: 60x).

**Build the diagnostic instead.** Solve once without the constraint and report

```text
asym = || x - Pi x || / || x ||
```

This measures directly whether the map is symmetric and aligned -- the
assumption the ASU restriction makes silently. It costs one unconstrained solve,
reuses the `Pi` the constraint work already provides, and gridding cannot
produce it at all. Report it per shell as well as globally: mismatch
concentrated at high resolution indicates flexibility or genuine symmetry
breaking, whereas a low-resolution component indicates misalignment to the
principal axis.

### 2.5 Implementation contract

- `reconstructor_pcg` gains a `sym` handle set at `prep_particles` time, or a
  cached `symmats(3,3,M)` array alongside `rotmats`.
- `accumulate_absT2`: innermost loop over `g`, `loc = padf * [h,k,0] . (S_g R_i)`.
  The h-strided OpenMP colouring is computed from `h` only, so it remains valid
  -- but see 2.7.
- `apply_adjoint_all`: same replication for `b`.
- `apply_normal_matrixfree` must replicate identically, or it stops being the
  reference for the kernel. This makes it `M` times slower and effectively
  unusable for `pgrp /= c1` at large `M`; that is acceptable, but it means
  **symmetry makes `pcgop=kernel` mandatory in practice**.
- Nothing else changes. `precond_from_accum`, `finalize_khat` and `solve` are
  untouched: `H_sym` preserves the symmetric subspace by construction, so no
  per-iteration projection.

### 2.6 Do not

- Do not apply `Pi` to the iterate each PCG iteration. `M` volume rotations per
  iteration is unusable.
- Do not blend symmetric and asymmetric accumulators.
- Do not act on the 2.4 asymmetry diagnostic automatically -- report it, do not
  let it retune `pgrp`, the search space, or a symmetry weight.
- Do not add soft symmetry (`lambda_sym * (I - Pi)`) before the hard constraint
  validates against production. It is a cheap follow-on -- `(I - Pi)` is a
  projector so the operator stays SPD -- but it is not part of this work.

### 2.7 Two things to verify

**Scatter race under replication.** The colouring stride
`ceiling(sqrt(3)*wdim)` guarantees that two `h`-lines in the same colour map
`>= padf*stride` apart *for a single rotation*. With `M` rotations scattered
inside the same worksharing iteration, two threads holding different `h` now
write `M` footprints each, and separation must hold for every pair `(g, g')`,
not just `g = g'`. It does not in general. Options: hoist the `g` loop outside
the `!$omp do` (one colour sweep per symmetry element), or serialize `g`.
**Hoisting the `g` loop outside the colour sweep is the safe default.**
**DONE** -- the `g` loop is hoisted outside the colour sweep in
`accumulate_absT2` and `scatter_plane`.

**Axial views.** Particles whose central section is invariant under a subgroup
produce coincident replicas, which `rho` counts as independent. This affects
production too and is not introduced by this work, but it is worth a synthetic
test with an axis-concentrated orientation distribution. Open, section 10.

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

Estimates, not measurements, and the cold column predates the section 1.1
kernel-calibration change: the measured cold solve is now 13.4 s / 8 iterations
(see Status), not the ~25 s / 15 iterations estimated here. `finalize_khat` is
the only phase that does not shrink with `update_frac` and will dominate at
small subsets.

### 3.4 Artifact contract

- Previous artifact per `(state, eo)`: `acc_folded` (real, `cdim`) and
  `b_raw` (real, `box^3`).
- **`b` is persisted UN-deapodized.** `apply_adjoint_all` currently returns
  `E^-1 b_raw`; the envelope is box-dependent, so carrying the deapodized form
  across a box change would require undoing and redoing it. Persist `b_raw` and
  apply `deapod_mul` at solve time with the current box's envelope. This is
  free: `E^-1` is linear and fixed within an iteration, so blending raw and then
  deapodizing equals deapodizing and then blending.
- **Do not derive the previous artifact from a volume.** The contract is
  accumulator-to-accumulator.
- Dimension compatibility: see 3.5.
- `x_{t-1}` must also persist for the warm start. It is the output halfmap, so
  no new artifact -- but it must be the **unfiltered** solve output, not a
  postprocessed map.

### 3.5 Box changes: exact zero-pad, no rescaling

`autoscale` sets `scale = box_new/box_in` and `smpd_new = smpd_in/scale`
(`simple_magic_boxes.f90`), so `box * smpd` -- the **physical extent** -- is
exactly preserved across a downsampling change.

That single fact makes the whole problem trivial, and it is worth stating why.
The padded lattice has frequency spacing `1/(padf * box_crop * smpd_crop)`
`= 1/(padf * extent)`, which is **independent of `box_crop`**. So a padded
Fourier index `(h,k,m)` denotes the *same physical frequency* at every
`box_crop`; a smaller box is a centred **crop** of the same Fourier volume, not
a rescaled one.

Therefore:

```text
acc_new(h,k,m) = acc_old(h,k,m)   for |index| within the old lattice
               = 0                otherwise
```

**Index-for-index copy, zeros outside, and no scale factor.** This resolves what
was open question 8. There is no rescaling because the plane sample spacing in
*frequency* is also `1/extent` in both boxes -- the larger box samples the same
frequencies plus additional higher ones -- and `padf`, `padsc` and the KB
stencil are all box-independent.

Consistency with the resolution mask: at `box_crop`, `sqlp = (box_crop/2)^2`, so
`absT2_plane` already zeroed everything beyond the old ball. The zero-pad
therefore reproduces exactly what was accumulated, rather than inventing
content.

`b_raw` and the warm-start `x` are real-space and are carried by the same
transform: pad to the old padded lattice, FFT, zero-pad the spectrum to the new
padded lattice, inverse FFT, crop to `box_new`. Values are preserved pointwise
because the added coefficients are zero. Two FFTs each, once per box change.

Rules:

- **Growing the box is supported; shrinking is rejected**, mirroring
  `read_eos_parallel_io`. A shrink would discard accumulated high-frequency
  evidence with no way to recover it, and no workflow needs it.
- Rebuild `precond` and `Khat` from the padded `acc` at the new box. Both are
  derived (section 1), so nothing needs migrating. The kernel calibration
  constant is `padsc**2` and is box-independent (1.1); `env`/`invenv` and the
  deposition envelope are rebuilt by `new` and `finalize_khat` at the new box.
- The newly admitted shells are zero in the carried accumulator, so they are
  under-weighted by `(1-uf)` until they refresh -- the same self-correcting
  behaviour as 3.6, now with an explicit mechanism.

**Gate:** accumulate a fixture at `box_crop = 64`, grow to 128, and confirm the
padded accumulator plus a fresh full-batch accumulation at 128 give the same
solve as a full-batch 128 accumulation alone, to within the `(1-uf)` shell
deficit predicted above.

### 3.6 Interaction with the resolution limit

`acc` depends on `sqlp` through `absT2_plane`'s masking. When the NU-driven
limit advances between outer iterations (companion note section 5.3), the
previous `acc` was accumulated under a *lower* limit than the current one.
Blending them mixes two band limits.

This is benign in one direction only: the previous accumulator is zero above
its limit, so the blend under-weights the newly admitted shells by `(1-uf)`
until they refresh. It self-corrects over `~1/uf` iterations. Production has the
same property. **Record it, do not attempt to correct it**, but note it as a
reason to prefer a limit that advances slowly relative to `1/uf`.

### 3.7 Gates

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

Same structure, same grouping helpers, one `type(reconstructor_pcg)` object:

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

**Rule: `kill` the `reconstructor_pcg` object never, `reset` it per
`(state, eo)`.** `wimg` is state-independent (same `boxpd`), and killing it
destroys and rebuilds four FFTW plans per group.

Splitting accumulate from solve in time is what keeps peak memory at one
`(state, eo)` and is also the natural distribution boundary (4.4).

### 4.3 Budget at box 256

| Array | Size | Phase |
| --- | --- | --- |
| `acc` full-range complex | 1.08 GB | accumulate only |
| `b_accum` full-range complex | 1.08 GB | accumulate only |
| `wimg` padded image | ~0.54 GB | both |
| `Khat` real packed | 269 MB | solve |
| `precond` real packed | 269 MB | solve |
| `get_cmat()` copy | 539 MB | transient, per call |
| CG vectors `b,r,p,hp,z,x` | 6 x 67 MB = 402 MB | solve |

Accumulate peak ~3.0 GB (two full-range complex accumulators, see 5.1), solve
peak ~2.1 GB. Matches the measured ~3 GB.

Two reduction opportunities, both optional:

- **`acc` and `b_accum` as packed real instead of full-range complex** would
  save ~1.6 GB and is the single largest win. **A previous attempt at a packed accumulator was
  implemented and reverted**: production's `expand_cmat` materializes the
  periodic extension up front so its gather and scatter are consistent by
  construction, whereas this operator wraps on the fly, and a non-wrapping
  scatter breaks the exact self-adjointness CG requires. Any retry must
  reproduce the adjoint dot-product test first.
- `get_cmat()` returns a copy; `image` has `get_rmat_ptr` but no cmat
  equivalent. Adding one removes 539 MB of transient and measurable memcpy from
  both `apply_normal` and `apply_precond`.

  **This violates `image`'s encapsulation and the justification has to be
  written at the call site, not assumed.** `get_rmat_ptr` already carries the
  comment `VIOLATES ENCAPSULATION`, so the precedent exists and the cost is
  understood: a caller holding a raw pointer can outlive or bypass the owning
  image's state machine, in particular its FT flag, and no compiler check will
  catch it. The justification here is specific -- the copy is 539 MB per call in
  the two hottest routines, it is pure memcpy doing no arithmetic, and the
  pointer's lifetime is a single procedure body with no reallocation and no
  FT-state transition in between. If a future change makes any of those three
  false, the pointer must go back to being a copy. Any `get_cmat_ptr` added must
  carry that same `VIOLATES ENCAPSULATION` marker and this rationale.

### 4.4 Distributed split

`acc` and `b` are sums over particles, so they partition across parts and
reduce exactly like gridding's partial reconstructions. The solve does not
partition -- and does not need to: once reduced, it is particle-independent and
runs on one node with a serial tail that does not grow with `nptcls`.

Per-part memory is the accumulate peak (~3.0 GB at box 256), one `(state, eo)`
at a time, matching production's per-part footprint.

## 5. Particle I/O and preprocessing

### 5.1 `y_planes` must go -- this is the blocking defect -- DONE

(Historical; the defect below has been fixed -- see the DONE marker above and
the Status section. The commander now streams a bounded `y_batch` through the
accumulate/solve split.) The original commander allocated
`y_planes(lims2, lims2, nptcls)` and held every particle's Fourier plane
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

**What this does and does not buy, stated exactly.** The figure above is
*particle-image* residency, and that is the part which becomes independent of
dataset size. Total peak is not 264 MB, because accumulating `b` across batches
means its Fourier accumulator is live alongside `acc` for the whole batch loop:

| Live during the batch loop | box 256 |
| --- | --- |
| `acc` (full-range complex) | 1.08 GB |
| `b_accum` (full-range complex) | 1.08 GB |
| `wimg` padded work image | ~0.54 GB |
| batch images + planes | ~0.26 GB |
| **accumulate peak** | **~3.0 GB** |

So the two Fourier accumulators are the floor, and they do not shrink with
`MAXIMGBATCHSZ`. The win is that this floor is **constant in `nptcls`**, where
today's `y_planes` adds 2.5 GB at 5 000 particles and 246 GB at 500 000 on top
of it. At 5 000 particles the change is roughly 5.5 GB -> 3.0 GB; at 500 000 it
is the difference between running and not running.

`b_accum` must be accumulated in Fourier and folded **once** after the batch
loop, not folded per batch: `fold_and_ifft` is an inverse FFT plus a crop, and
doing it per batch would be both wrong (the crop is not linear over partial
sums in the deapodized domain) and needlessly expensive.

If the 3.0 GB floor becomes binding, the packed-real accumulator in 4.3 is the
lever -- it would remove ~1.6 GB of the 2.16 GB -- but note the prior failed
attempt recorded there.

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

The one case that would have been worth measuring -- `box_crop << box`, where a
cached cropped plane is far smaller than the raw read (~4.5 kB versus 256 kB at
`box_crop = 64`) -- **no longer arises**. That regime is the early workflow
stages, which stay on gridding (see Status). In the target regime the boxes are
large and the cache-entry-size argument above applies in full, so the case is
closed rather than deferred.

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

**DONE.** The LUT hoist (`build_hk_luts`: `spafreqsq_lut`, `ang_lut`,
`shell_lut`) and the kernel form are both in the tree. `absT2_plane` and
`build_transfer` now hoist the per-particle CTF constants via `get_ctfvars` and
evaluate the flat `ft_map_ctf_kernel` expression (one `cos` + one `sin`, no
3-deep type-bound dispatch, no per-pixel `canonical_phshift`). `ft_map_ctf_kernel`
is **inlined, not called**: its memoized `(h,k)` maps (`memoize_ft_maps`) span the
`h >= 0` non-redundant half only, but `lims2` is a full both-sign-h disk, so the
library routine would index its maps out of bounds for every `h < 0`; the LUTs
are ranged over `lims2` and carry no such restriction. The kernel form is
numerically-equivalent-not-bit-identical (one FP reassociation in
`half_wl2_cs = 0.5*wl*wl*Cs`, factored out of the per-pixel product), so it is
validated by the operator stages, **not** the residual trace: `test=pcg_recon`
stage-4 corr `0.97559` vs prior `0.97560`, stage-5 interior rel_err `2.10e-2`,
stage-6 shift-invariance `0.0`, all eight stages PASS. Transcendental *count* is
unchanged, so the win is call/arithmetic overhead only -- and the setup re-run
confirms that overhead was never the bottleneck: `scatter |T|^2` = 15.9 s
post-form vs 15.1 s pre-form (within noise). The 15.9 s is the KB-window lattice
scatter (`apod_mat_3d_fast` + 27-tap `scatter_window`), the next target.

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

## 6. Standalone trailing emulator

Fractional update is required for `refine3D_auto` / `refine3D_multi` integration but
**not** for establishing that the PCG reconstructor works. The two should not be
entangled: a trailing bug and a reconstruction bug look identical in an FSC
curve, and debugging them inside a refine3D loop means paying a full workflow
iteration per experiment.

So the recursion of section 3 is developed and validated in a **separate module
driven by a test commander**, taking particles with assigned orientations and
emulating the fractional-update schedule directly. No workflow integration, no
project write-back, no dependency on the samplers.

**Module** `src/main/volume/simple_pcg_trailing.f90`, type `pcg_trailing`,
owning the carried state and nothing else:

| Component | Content |
| --- | --- |
| `acc_folded` | carried accumulator, packed real |
| `b_raw` | carried RHS, un-deapodized (3.4) |
| `x` | previous solution, warm start (3.3) |
| `box_crop`, `smpd_crop` | current lattice |

| Procedure | Contract |
| --- | --- |
| `new(box, smpd, ...)` | wraps one `reconstructor_pcg` |
| `set_box(box_new, smpd_new)` | zero-pads the carried state (3.5); rejects shrink |
| `update(pinds, uf)` | accumulate subset -> blend -> warm-started solve |
| `get_vol(x)` / `write_vol` | current estimate |
| `read_artifact` / `write_artifact` | persistence, for the restart test |

`update` performs exactly the recursion of 3.2 and nothing else. The subset is
**supplied by the caller**, not chosen: selection policy stays in
`simple_matcher_smpl_and_lplims.f90` and is out of scope here (3.1).

**Test commander** `test=pcg_trailing`: reads a project with assigned
orientations, partitions particles into `R` rounds of size `nptcls*uf`, runs the
recursion, and compares against the full-batch solve. Assertions:

1. `uf = 1` reproduces the full-batch solve bit-for-bit.
2. A stationary particle set under `uf < 1` converges to the full-batch answer
   as `R` grows -- report `||x_R - x_full|| / ||x_full||` against `R`.
3. **Warm-start convergence**: iterations to `rtol` per round, which is the
   unmeasured assumption the entire cost argument of 3.3 rests on. If this does
   not drop to 2-3, the fractional-update design does not pay and the note needs
   revisiting.
4. **Box growth**: run rounds at `box_crop = 64`, switch to 128 mid-schedule,
   confirm the 3.5 gate.
5. Round ordering does not change the converged answer beyond the expected
   forgetting-window effect.
6. Write/read artifact round-trip is exact, and a restart from a persisted
   artifact matches an uninterrupted run.

Emulating rather than integrating also makes the awkward experiments cheap:
adversarial round orderings, `uf` schedules, deliberately stale subsets, and box
changes at arbitrary points -- none of which a real workflow would let you
schedule.

## 7. Remaining gaps

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

## 8. Sequencing

Three phases. Phase A is everything needed to test the reconstructor on real
data with assigned orientations; it does **not** include fractional update.

**Phase A -- make the reconstructor testable on real data.**

1. **Batching and memory.** [DONE, verified on real data] Remove `y_planes`,
   fuse the accumulate into the batch loop (5.1). Blocking defect: nothing runs
   at real dataset size until this is done, and it is also the accumulate/solve
   split everything later needs.
2. **Cheap prerequisites.** [1.1 DONE; 5.4 DONE incl. `ft_map_ctf_kernel` form,
   `test=pcg_recon` green; setup re-run done -- no win, see below] Remove
   `calibrate_kernel`'s particle pass (1.1) and hoist the `absT2` `(h,k)` lookup
   tables plus adopt the transcendental-free CTF form (5.4). Both pay for
   themselves immediately and 1.1 is a hard prerequisite for phase B. The 5.4
   hoist precomputes `spafreqsq`, `ang` (`atan2`) and `shell` (`sqrt`) once over
   the fixed `lims2` disk in `build_hk_luts`, shared by `absT2_plane` and
   `build_transfer`; the LUT expressions are bit-identical but the adopted
   `ft_map_ctf_kernel` form is numerically-equivalent-not-bit-identical (see
   5.4), so it is validated by the operator stages (stage-4 corr `0.97559`,
   stage-5/6 agreements), not the residual trace. Per-phase setup
   instrumentation and a 3-way accumulation sub-split (read / prep / scatter)
   were added to the commander to localize the cost: `scatter |T|^2` is ~15.9 s
   of ~25 s setup, and the post-CTF-form re-run shows no measurable change
   (15.9 vs 15.1 s, within noise) -- the cost is the KB-window lattice scatter,
   not CTF eval, so that is the next target.
   *Not originally scoped but done here:* the colouring race and its fix, the
   merged `test=pcg_recon`, and a `dx/x` early-stop in `solve_core`
   (`PCG_XTOL = 1.5e-2`, exits alongside `rtol`; fired at iter 27 in the stage-8
   run). See "Where this stands".
3. **Symmetry** (2). Coordinate replication **DONE** (see 2.3): implemented in
   `reconstructor_pcg`, wired into `reconstruct3D_pcg` (guard lifted), and gated
   by `test=pcg_recon` stage 9 (in-operator c2 build produces the same `Khat`
   and `b` as a c1 build of the expanded set).
   Still open: the lattice-exact permutation path gated against replication, and
   the asymmetry diagnostic (2.2, 2.4).

At the end of phase A the command reconstructs symmetric particles from a real
project at real dataset sizes, and can be compared against production
`reconstruct3D` on the same particles.

**Phase B -- nail fractional update in isolation.**

4. Standalone trailing emulator (6): module, test commander, and the six
   assertions. Includes box growth (3.5) and the warm-start measurement that
   the cost argument depends on.

Phase B needs nothing from phase A except the accumulate/solve split, so the two
can overlap if convenient. It deliberately does not touch any workflow.

**Phase C -- workflow integration.** Only after A and B report.

5. Even/odd half sets + FSC; `box_crop`/`smpd_crop` and the NU-driven `sqlp`.
6. Particle weights.
7. Multi-state and the single-object memory scheme (4.2); distributed
   accumulator reduction (4.4).
8. Project write-back and postprocess handoff.

Note the reordering against earlier drafts: even/odd was previously step 2 on
the grounds that nothing is measurable without it. That still holds for
*resolution* claims, but not for the phase A goal -- a symmetric reconstruction
at scale can be compared against production directly, map to map and time to
time, without a half-set split. Half sets move to phase C where FSC-based
acceptance actually starts.

## 9. Acceptance gates

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

## 10. Open questions

1. Does the axial-view multiplicity (2.7) matter on real symmetric data?
2. Is the ~8 s residual rim cost really type-bound dispatch on `win_wraps`?
   Unverified; see "Where this stands". If the module-level change does not
   move setup, instrument setup per phase instead of guessing again.
3. Should the wrapping rim be removed entirely by giving the accumulator a
   margin, as production's `expand_cmat` does, rather than serializing it? That
   would restore full parallelism and is the same mechanism that defeated the
   earlier packed-accumulator attempt (4.3). Bigger change; revisit if the rim
   cost stays material.
4. Is the signed-permutation argument (2.2) right in practice -- does the
   permutation path reproduce coordinate replication to round-off on `d2`, and
   does `get_sym_rmat` return exactly-integer matrices for the lattice-exact
   groups or only near-integer ones? If the latter, snapping to the nearest
   signed permutation needs a tolerance and a justification.
5. Does warm-started PCG actually converge in 2-3 iterations under trailing, or
   does the changing operator break the warm start's usefulness? The whole cost
   argument in 3.3 depends on this and it is unmeasured.
6. Is the kernel adequate across the `box_crop` values `refine3D_auto` actually
   uses? Its ~3% error was measured at box 256 only. The smallest boxes are no
   longer in scope, but `autoscale` still varies the box within the target
   regime.
7. How should `lambda` and `RHO_FLOOR_FRAC` scale with `uf` and with `M`? Both
   are calibrated against the data term's absolute scale, which both mechanisms
   change.
8. *(resolved in 3.5: exact zero-pad, no rescaling, because `autoscale`
   preserves `box * smpd`. Retained here only so the resolution is on record.)*
9. *(withdrawn: the small-`box_crop` early stages stay on gridding, so the
   cropped-plane cache question of 5.3 is closed rather than open.)*
10. `sig2` at 258 MB per 500 k particles (5.2) -- worth replacing with a group
   index plus per-group spectra, or is it below the noise floor?
11. Is squared conditioning (5.3.1) ever the binding constraint? If preconditioned
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
