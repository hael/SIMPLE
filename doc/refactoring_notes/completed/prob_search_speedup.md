# The probabilistic search: paying less for the same search

Status: completed 14 October 2026 in an interactive session (not the planned
unattended run); the report is `prob_search_speedup_report.md` beside this
file. Phases 0-3 landed behind one key, `prob_fast` (ruling 7, which
supersedes the three keys of the plan below); Phase 4 was skipped on its
measurements; the subspace-map table of the neighbourhood search was added on
the default path (ruling 8). The last batch is written but not yet built or
run (report, "Open"). Section 11 is the record. Written 2026-10-13; revised
14 October 2026 after review (section 13).

The probabilistic pre-step (`prob_align` and its `prob_tab` parts) was 46% of
a three-state ribosome `solve3D` run
(`doc/implementation_notes/rejected/particle_caching_3D.md`). This plan
times it, then makes it cheaper:

| Phase | Change | Kind | Key |
| --- | --- | --- | --- |
| 0 | Benchmark files, baseline, thread scaling (3.4) | exact | none |
| 1 | Particle-frame coarse grid in the multi-state refinement (3.1(a)) | approximate | `prob_grid_ptcl` |
| 2 | Particle-frame seed in the single-state table (3.1(b)) | approximate | `prob_seed_ptcl` |
| 3 | Quadratic shift step after the coarse grid (3.2) | approximate | `prob_grid_fit` |
| 4 | Blocked table loop, if Phase 0 or 3 calls for it (3.3) | exact | none |
| 5 | Documents and close | | |

Each approximate change has its own default-off developer key and is
validated on its own against the baseline. The keys are temporary: after
review the maintainer either promotes a method, removing its key, or removes
the method with its key.

## 1. Aim

- **Faster pre-step, same search.** An exact change keeps the table's
  contents, up to the order of random draws inside a thread (3.3).
- **A breakdown of the pre-step's time.** The pre-step gets benchmark files
  of its own, so its time is split into named parts instead of being one
  number in the refine3D strategy file.
- **Measurement over argument.** The table-loop change is decided by a
  measurement.

## 2. What the pre-step costs today

Measured (the ribosome chunk, consolidated layout, 29,101 particles, three
states, `nspace` 1,000, 5 parts of 16 threads, box 336): the pre-step is
3,033 of 6,572 s over 78 iterations, 37 s per iteration at `box_crop` 132,
50 s at 160 and 138 s of a 190 s iteration at 182. Those are the run's
`refine=prob` stages. In today's `solve3D`, `refine=prob` runs in the middle
stages (3–5). Single-state runs use `prob_neigh` before and after them;
multi-state runs stop after stage 5 by default (`solve3D_policy.md` sections
3 and 6). Inside the pre-step only the whole is timed (`rt_prob` in
`simple_refine3D_strategy.f90`).

Counted from the code, per sampled particle, at `box_crop` 182 (`nk` 60
shells, `pftsz` 254, `nrots` 508):

- **The multi-state table.** `fill_tab_range` in `simple_eul_prob_tab.f90`
  (`process_particle_with_shift`) uses one thread per particle. The shift
  seed is single-state only (`l_sh_first`), so with more than one state
  every one of the 3,000 references is scored at zero shift by
  `gen_prob_likelihood_euclid_val`. That is three steps:
  - the memoised contraction `ft_ctf2 · ft_ref2 − 2 ft_ptcl_ctf ·
    conj(ft_ref)` over `nk × (pftsz+1)`;
  - one real inverse FFT of length `nrots`;
  - a bounded heap for the best `nsample` rotations
    (`sample_likelihood_dist`, O(n log nsample)).

  Call the contraction plus the inverse FFT the base kernel: about 0.32
  Mflop, reading 248 KB of reference and 248 KB of particle memos. 3,000 base
  kernels per particle.
- **The single-state table.** The seed comes from `estimate_shift_seed`: one
  coarse grid and L-BFGS-B per particle against the particle's previous
  reference. In alternating mode the search starts from the best rotation at
  zero shift, not from the incoming one. The seed is a reference-frame shift, and every reference is scored
  at it, so `gen_euclid_crvec` takes its shifted-reference branch. That
  branch builds the shift matrix, multiplies the reference by it, and runs
  `nk` complex FFTs of length `nrots` before the contraction. It costs about
  5.3 base kernels per reference, `nspace` × 5.3 per particle.
- **The refinement** (`refine_best_refs`). It refines the `projs_ns` best
  references, `projs_ns = athres/180 × nspace`, where `athres` is the smaller
  of `prob_athres` (10°) and the sampled particles' mean angular move: 17 to
  55 references. Each gets `grad_shsrch_minimize`
  (`simple_pftc_shsrch_grad.f90`) in alternating mode with `coarse_init` on.
  The single-state case turns `coarse_init` off by passing `xy_in`. Per
  refined reference:
  - one base kernel at zero shift;
  - the coarse grid: `coarse_num_steps`² = 25 calls of `gen_objfun_vals`,
    each shifting the reference (`gen_shmat4aln`, then `nk` complex FFTs on
    `plan_fwd1_many`) before the contraction. About 1.7 Mflop each (5.3 base
    kernels), the FFTs four fifths of it.
  - up to `max_evals` = 5 alternations of L-BFGS-B (`maxits_sh` = 60), each
    followed by one shifted `gen_objfun_vals`. Every L-BFGS-B evaluation
    (`gen_corr_grad_for_rot_8` → `gen_euclid_grad_for_rot_8`) copies the
    reference to double precision (244 KB) and sums residual and gradient
    directly over `nk × pftsz`: about 0.6 Mflop in double.
  - The search object counts the evaluations per run, E
    (`profile_objective_evals`, `profile_gradient_evals`), and has a getter
    and a reset (`get_profile`, `reset_profile`), but nothing calls them.
    Alternations are not counted.

  Per refined reference that is 1 + 25 × 5.3 + alternations × (2E + 5.3) base
  kernels: about 250 to 450, for two to three alternations and E of 20 to 50.
  Per particle it is 4,000 to 25,000, against 3,000 for the table.

**Consistency check.** 5,000 particles × about 13,000 base kernels × 0.32
Mflop is 2 × 10¹³ flop per iteration. Over 138 s on 80 threads that is 1.9
Gflop/s per thread, a plausible rate for these loops. So the pre-step is
these two loops plus 10 to 15 s per iteration of process launches,
reprojection-model reads and memoisation, the table write and read, and
`ref_assign`.

**Where the grid runs.** It runs once per particle in:
- the matcher (`grad_shsrch_first_obj` in `simple_strategy3D_srch.f90`);
- the 2D search;
- the 2D table (`simple_eul_prob_tab2D.f90`);
- every single-state seed.

It runs once per refined reference in the multi-state 3D table, and in the
multi-state path of the neighbourhood table (`simple_eul_prob_tab_neigh.f90`).
The neighbourhood table is out of scope.

## 3. Target design

### 3.1 Particle-frame shifts in the table (approximate; two keys)

Two changes, both in the 3D table's `fill_tab_range`, each under its own key:
- (a) the coarse grid of the multi-state refinement (`prob_grid_ptcl`);
- (b) the seeded scoring of the single-state table (`prob_seed_ptcl`).

They share the polar calculator's builder of shifted particle memos. They
never act on the same particle: (a) needs more than one state, (b) exactly
one.

**No exact particle-side form.** The kernels shift the reference in its own
polar frame (`shmat(:,k) * pfts_refs(:,k,iref,ieo)` in
`gen_raw_euclid_vals_impl` and `gen_euclid_crvec`), then correlate over every
in-plane rotation by FFT. At rotation r, a reference-frame shift s therefore
acts as the particle-frame shift R_r s. The code converts this way itself:
`grad_shsrch_minimize` rotates the shift by `get_rot(irot)` when `sh_rot` is
set, and `materialize_seed_shift` rotates the seed by the assigned in-plane
angle. A fixed particle shift is a different reference-frame shift at every
rotation, so a fixed bank of particle shifts reproduces the reference-frame
grid at one rotation only. Restricting the change to fixed-angle searches
gains nothing: the table uses the alternating object, and the fixed-angle
`coarse_search` scores single rotations by direct sums with no FFT to save.
Caching the grid's 25 shift matrices would be exact, but it is worth about 5%
of a grid point and is not pursued.

**Exact pointwise.** For one particle shift t and one rotation r, the
particle-side value equals the reference-side value at the reference-frame
shift that t corresponds to at r, to single-precision rounding:
|X·S(t) − CTF·R|² = |X − CTF·R·S(−t)|², with |S| = 1 and a real CTF, so
`ft_ctf2` and the particle's weighted power are unchanged. Both changes rest
on this identity, and a test pins it.

**(a) The multi-state coarse grid (`prob_grid_ptcl`).**

- **Lattice.** The grid's points are recorded shifts u, in the particle
  frame like every recorded shift: today's 5 × 5 lattice (spacing h =
  2·trs/5, at 0, ±0.4·trs and ±0.8·trs), 25 points. Converted to the
  reference frame at a rotation, a corner point (|u| = 1.13·trs) can fall
  outside the search box; the start is then clamped into the box.
- **Bank.** A search object built with the option makes, when its particle
  changes in `set_indices`, one bank entry per lattice point. Entry u holds
  the particle's `X·CTF` shifted by −u and transformed over rotation, as
  `ft_ptcl_ctf` is in `memoize_ptcls`. It scores as the reference shifted by
  u in the particle frame. The sign lives in the builder only. Size: 25 ×
  (pftsz+1) × nk complex, 3.1 MB per thread at `box_crop` 182.
- **Grid.** Each entry is scored by a variant of the unshifted branch (Euclid
  and CC) that takes the entry in place of `ft_ptcl_ctf(:,:,i)`: one base
  kernel per lattice point, all rotations. The grid keeps the lowest
  (u*, r*) and every point's loss vector (for 3.2).
- **Start.** The search starts at rotation r* from the reference-frame shift
  that u* corresponds to there (the inverse of the `sh_rot` conversion),
  clamped into the box. It then continues as today, with today's acceptance
  against the zero-shift cost.
- **Scope.** `fill_tab_range` passes the option at construction only when the
  key is on and the table has more than one state. Every other
  `new_alternating` caller is unchanged by construction. A
  one-grid-per-particle consumer would gain nothing, because its bank costs
  what its grid costs.
- **What changes.** The lattice samples the particle's own shift, the same at
  every rotation, with the same 25 points as today. Today's grid samples a fixed reference shift, which is a
  different particle shift at every rotation. So the start of L-BFGS-B, and
  through it the recorded distances, may change.
- **Cost.** Per refined reference: 25 base kernels instead of 25 × 5.3 ≈ 133.
  Per particle: 25 ramps and FFT sets (about 25 × 4.3 base kernels' worth),
  shared by the 17 to 55 refined references.

**(b) The single-state seeded table (`prob_seed_ptcl`).** `estimate_shift_seed`
finds a reference-frame shift s at the rotation r_s where its search ends
(`minimize` returns it in `irot`, or 0 with a zero shift). Today every reference is scored at s as a
shift of the reference, and every refined reference starts from s (`xy_in`).
An unrefined assigned candidate gets s rotated by its own in-plane angle
(`materialize_seed_shift`). Away from r_s, s is a different particle shift.
With the key on, the seed is used throughout as one particle-frame shift u_s,
which is s converted at r_s (what `sh_rot` gives). Three things change
together, so that every recorded shift matches how its candidate was scored:

- **Scoring.** The particle's `X·CTF` is shifted once by −u_s into a
  one-entry bank (the builder of (a)). Every reference is scored with the
  unshifted-kernel variant taking that entry, in place of the
  shifted-reference branch (Euclid and CC). A seed below `SHERRSQ` uses the
  plain memo, as today.
- **Refinement start.** A refined reference at rotation r starts from u_s
  converted to the reference frame at r.
- **Assignment.** The seed table holds u_s.
  - In `ref_assign`, `assign_current_ref` records u_s for an unrefined
    candidate, whatever its in-plane angle, instead of calling
    `materialize_seed_shift`. The rule is a small module procedure of
    `simple_eul_prob_tab.f90`, so it can be tested without a project.
  - `materialize_seed_shift` stays as it is for the 2D and neighbourhood
    tables.
  - `prob_tab` and `prob_align` read the same command line, so the table file
    format does not change.

The seed's estimation is unchanged. At r_s, every score, start and recorded
shift equals today's; elsewhere they differ. Cost: `nspace` × 5.3 base
kernels per particle become `nspace` + about 4.3, so with a non-zero seed
about a fifth of today's.

### 3.2 A quadratic shift step from the grid (approximate; `prob_grid_fit`)

**The table's refinement can be made cheaper, not dropped.** `ref_assign`
consumes the refined distances, because `refine_best_refs` replaces the
`projs_ns` entries of `loc_tab`:
- with more than one state, the state is the argmin distance over all
  references (`assign_greedy_state_labels`);
- the projection is drawn from the distances within that state
  (`assign_refs_for_state`).

The matcher then takes state, projection and in-plane index as given
(`srch_prob` in `simple_strategy3D_prob.f90`) and only polishes the pose.

**Scope.** The step runs in the multi-state refinement, in the objects
3.1(a) scopes, after the coarse grid. The grid keeps its 25 loss vectors in
either frame: the reference-frame grid with `prob_grid_ptcl` off, the
lattice of 3.1(a) with it on.

**The step.**

1. **Neighbourhood.** Take the grid minimum g*, its rotation r*, and the
   3 × 3 stencil around g*, all at r*. If g* is on the edge of the 5 × 5
   lattice, there is no fit (9 of the 25 points have a full stencil). A
   stencil shifted inside the lattice for an edge g* was tried and rejected
   (section 11).
2. **Fit.** Fit the six-parameter least-squares quadratic over the nine
   values, in units of h: f(a,b) = c + gᵀ(a,b) + ½ (a,b) H (a,b)ᵀ, in closed
   form on the regular stencil.
3. **Curvature.** H must be positive definite, with λmin/λmax ≥ 0.01;
   otherwise there is no fit.
4. **Step.** δ = −H⁻¹g must lie in the centre cell, |δa|, |δb| ≤ ½;
   otherwise there is no fit. A stationary point outside it contradicts g*
   being the grid minimum.
5. **Limits.** The fitted shift, converted to the reference frame at r*,
   must lie in the search box; otherwise there is no fit.
6. **Value.** Evaluate exactly once at the fitted shift over all rotations
   (about 5.3 base kernels); the rotation is its argmin. The fit is kept only
   if its cost is below the grid minimum's; otherwise there is no fit.

"No fit" falls back to today's L-BFGS-B alternations from the grid minimum.
Either way the result passes today's acceptance against the zero-shift cost.
The search object counts fits, and fallbacks per reason (steps 1, 3, 4, 5
and 6).

**Second form.** It is built only if the fit alone fails the case
acceptance of section 10 (Phase 3): the fit, then L-BFGS-B from the fitted
shift with `maxits` 10 and one alternation.

**Cost.** On the fitted fraction, about 5.3 base kernels replace
alternations × (2E + 5.3) ≈ 100 to 300 per refined reference. The fallback
fraction keeps today's cost.

### 3.3 The table loop (exact; Phase 4)

The table loop streams, per particle and thread, the whole memoised
reference bank of the particle's parity (3,000 × 248 KB = 744 MB), at under
one flop per byte.

**Blocking.** A thread takes B particles and, for each reference, evaluates
the B pairs before moving on. A reference's memos are then loaded once per
block. The per-pair work is unchanged.

**Candidate storage comes with it.** Today's buffers are append-only, and
`prob_candidate_buffer%append_or_replace` finds a particle's earlier
candidate by scanning back while the particle index matches. Reference-major
blocking interleaves particles, so refined candidates would be appended as
duplicates. `read_tab_to_glob` would then refuse the table, because every
particle must have exactly `nrefs` candidates. So blocking uses block-local
dense storage:
- per thread, arrays over (`nrefs`, B) for distance and in-plane index (as
  `dists_refs` and `inpl_refs` today), plus shift and `has_sh` for refined
  entries;
- the score loop fills them reference-major;
- then, particle by particle, `refine_best_refs` writes its refined entries
  in place, and the particle's `nrefs` candidates are emitted with `append`
  in reference order.

The stream keeps today's layout and contents, and the 3D table no longer
calls `append_or_replace`. Memory: `nrefs` × B × 24 bytes per thread, 1.2 MB
at 3,000 references and B = 16.

**Scope.** Every table scored with the unshifted kernel: the multi-state
table, `doshift=no`, and the single-state table with `prob_seed_ptcl` on,
which needs one seed entry per block particle (B entries per thread).
The single-state table with that key off shifts the reference per pair, so a
block would share no memos there, and it keeps today's loop.

**Determinism** (ruled: default path, no key). Within a thread, the order of
the `ran3` draws that pick the in-plane rotation changes, so the picks are
the same in distribution, not draw for draw.

**Product form.** Writing the contraction as BLAS matrix products over
blocks of particles is not in this run. A later plan takes it up if blocking
leaves the score loop the larger part.

### 3.4 Benchmark files

- **`PROB_TAB_BENCH_ITER###_PART##.txt`**, written by `exec_prob_tab` in the
  matcher's format (`CONTEXT`, `TIMINGS`, `% accounted for`):
  - **Wall-clock entries**, measured outside parallel regions: project read,
    reprojection-model read, reference memoisation, particle preparation
    (summed over batches), the table fill (the parallel region of
    `fill_tab_range`, summed over batches), the table write. `% accounted
    for` uses only these and the total.
  - **A separate thread-seconds block.**
    - Contents: seed estimation, score loop and refinement, each summed over
      threads; their sum; and the load balance, which is that sum divided by
      the fill region's wall time × `nthr`.
    - Per-thread accumulators read `system_clock` into thread-local
      variables. They do not use `tic`/`toc`, which write the module
      variable `last_time_point`.
    - These are labelled thread-seconds and never enter `% accounted for`.
  - **Context:** references, sampled particles, `projs_ns`, alternations (a
    new counter in `grad_shsrch_minimize`), grid evaluations; the search
    objects' counts from `get_profile`, summed after the loop. E is the
    gradient count, one per L-BFGS-B evaluation; the objective count adds
    `nrots` per all-rotation evaluation. With the keys on, the file also
    records fits and fallbacks per reason; Phases 1 and 3 add those lines.
- **`PROB_ALIGN_BENCH_ITER###.txt`**, written by `exec_prob_align`, all wall
  time: sampling and project write, scheduling (launch to last part
  finished), table reads, `ref_assign`, assignment write.
- The file names go in `simple_refine3D_fnames.f90`, beside
  `refine3D_bench_fname`.

## 4. What does not change

- **With every key off (the default):** the search, the table's contents (up
  to the in-plane draw order once blocking is in), and the table file format.
- **With or without the keys:**
  - the estimation of the single-state seed;
  - the matcher;
  - every one-grid-per-particle consumer;
  - the 2D and neighbourhood tables;
  - the multi-state rule of `refine3D_policy.md`, that no seed from the
    particle's previous state scores other states.
- **The shift convention:** every recorded shift is in the particle frame,
  as after `sh_rot`.

## 5. The keys

**Superseded (ruling 7, 14 October 2026): one key, `prob_fast`, for all three
changes.** The text below is the original plan.

`prob_grid_ptcl`, `prob_seed_ptcl` and `prob_grid_fit`, each
`(yes|no){no}`, one per approximate change. They are registered as
`prob_athres_taper` was (commit `cfbc56bf3`):
- the parameter in `simple_parameters.f90` and its registration in
  `simple_parameters_parse.f90`;
- a `UI_VIS_DEVELOPER` entry on refine3D and on solve3D (in `new_solve3D`,
  not in `add_solve3D_base_inputs`, which the stream emulators share).

solve3D's stages copy its command line into refine3D's
(`simple_solve3D_utils.f90`), and refine3D's reaches `prob_tab` through
`prob_align`. The keys are read in `fill_tab_range`, and `prob_seed_ptcl`
also in `ref_assign`. No key stays: on the maintainer's decision after review
a validated method becomes the only path and its key is removed (release 4
keeps no backward compatibility); a method that is not promoted is removed
with its key.

## 6. Tests

The new tests run on one thread on small synthetic polar transforms, within
the fast gate's budget.

- **Particle-frame grid** (Phase 1), in `simple_pftc_inplane_tester.f90`:
  - *Identity.* On random particle and reference transforms with CTF, for
    lattice points u and rotations r (every r for a few u), the bank entry's
    loss L at r and the reference-side loss at the converted reference-frame
    shift agree within the tester's existing route tolerance,
    5e-4·(1+|L|), under both objectives.
  - *Clamping.* A converted start outside the box is clamped into it; one
    inside is unchanged.
  - *Recovery.* Synthetic particles with known shift and rotation, on and off
    the lattice, noise-free and noisy. Grid plus L-BFGS-B recovers them as
    well as the reference-frame grid, at the tolerances of
    `test_euclid_routes_and_recovery`.
  - *Option off.* Without the option there is no bank and today's result.
- **Particle-frame seed** (Phase 2):
  - *At r_s* (inplane tester). Scoring with the one-entry bank gives today's
    shifted-reference loss at r_s, within 5e-4·(1+|L|), under both
    objectives.
  - *Start* (inplane tester). The start converted from u_s at r_s equals s.
  - *Assignment* (3D table tester, below). With the key on, the recording
    rule gives u_s whatever the in-plane angle; with it off, today's
    `materialize_seed_shift` value.
- **Quadratic step** (Phase 3), in the same tester:
  - *Guards.* An edge minimum, a non-positive-definite or ill-conditioned H,
    a step outside the cell, a shift outside the box, and a worse exact
    value. Each takes the fallback, with a result identical to today's path.
  - *Accuracy.* Noise-free, over a range of true shifts, the fitted shift
    lies within 0.2 pixel of L-BFGS-B from the same grid point.
- **3D table tester.** A new `simple_eul_prob_tab_tester.f90` (the 3D table
  has none). Phase 2 creates it with the assignment case. It is registered as
  the sub-suite "3D probability table" after "2D probability table I/O" in
  `suites_pftc_align2D3D`, and listed in that suite's `suite=` help in the
  test UI.
- **Blocked loop** (Phase 4), in the 3D table tester:
  - *Placement.* The blocked score loop is a module procedure of
    `simple_eul_prob_tab.f90` over the polar calculator, a block of
    particles and the references.
  - *Identity.* With the in-plane sample size forced to 1 (a tiny
    `prob_athres`), `sample_likelihood_dist` picks the argmin. It still draws
    one `ran3` number, which cannot change the pick. The blocked and
    per-particle paths must then emit identical candidate streams, field by
    field, with exactly `nrefs` candidates per particle, and the stream must
    read back with `read_tab_to_glob`.
- **Benchmark files.** `% accounted for` above 90% of the wall time; the load
  balance printed.

## 7. Documents

- `doc/policies/3D/refine3D_policy.md`: the pre-step's benchmark files; each
  key while it exists (what it changes, developer key pending promotion).
- `.github/skills/simple-main-pftc/SKILL.md`: the shifted particle memos as
  part of the memoisation design, if a particle-frame key is promoted. The
  text is proposed in the report and applied on approval.
- On completion this plan moves to `doc/refactoring_notes/completed/` with
  its report.

## 8. Consequences

- **Memory per thread:**
  - `prob_grid_ptcl`: 3.1 MB of bank;
  - `prob_grid_fit`: 25 loss vectors (51 KB);
  - `prob_seed_ptcl`: one 122 KB entry, B with blocking;
  - blocking: 1.2 MB.
- **Keys off:** nothing changes but the benchmark files and, after Phase 4,
  the loop order.
- **Estimates by the counts of section 2:**
  - Phase 1 removes a quarter to two fifths of the multi-state refinement.
  - Phase 2 cuts the single-state table to about a fifth, in the `refine=prob`
    stages only.
  - Phase 3 replaces 100 to 300 base kernels per refined reference with about
    5 on the fitted fraction.
  - Phase 4 pays only if the loop is bandwidth-bound.
  - Phase 0 replaces these estimates with numbers.

## 9. Rulings

All given 14 October 2026.

1. **Dropping the table's continuous refinement: withdrawn.** The refined
   distances decide state and projection (3.2); the build order of 3.2 takes
   its place.
2. **Per-state seeds for more than one state: refused.** The multi-state
   rule of `refine3D_policy.md` stands, and the multi-state refinement keeps
   a grid per refined reference.
3. **The search counters: void.** They have `get_profile` and
   `reset_profile`, and the benchmark files print them.
4. **The `ft_ctf2 · ft_ref2` term shared per micrograph: dropped.** With patch
   CTF (`ctfpatch=yes`, the default), extraction gives every particle its own
   defocus (`pix2polyvals` in `simple_extract_strategy.f90`), so the term
   differs per particle.
5. **The single-state seeded table in the particle frame: accepted, under its
   own key** (`prob_seed_ptcl`, Phase 2). It is validated, kept or removed
   independently of `prob_grid_ptcl`.
6. **Blocking's change of the in-plane draw order: accepted on the default
   path, without a key**, on the identity test of section 6.
7. **One key** (14 October 2026, after the measured runs of section 11):
   `prob_grid_ptcl`, `prob_grid_fit` and `prob_seed_ptcl` are removed;
   `prob_fast` switches on all three, each in its own scope (the grid and the
   fit in the multi-state table, the seed in the single-state table). This
   replaces ruling 5's separate key.
8. **The subspace-map table on the default path, without a key** (14 October
   2026): `sym_dist_table` is a refactoring of the pair loop with the same
   distances.

## 10. Execution

An unattended arun run, `prob_search_speedup`, on the Dell (Oracle Linux
8.10, 24 cores). Nothing is committed; the maintainer reviews the complete
diff. Every phase records itself in section 11. Evidence goes under the
run's `scratch/`, kept projects and benchmark files under `scratch/keep/`.

**Build.** As recorded in `global_sigma2_only.md` (section 13, Phase 0):
- GCC 15.2 (gcc-toolset-15), CMake 3.26.5, `OPENBLAS_CORETYPE=Haswell`, and
  Python 3.12 as `python3` for the gate's checkers;
- `./compile_debug.sh` (Debug, runs the fast gate);
- an out-of-tree Release build with tests (`BUILD_TESTS=ON`,
  `USE_ARCHOPT=ON`, installed in its own tree) for every case;
- every case sets `SIMPLE_PATH` to the Release build and puts its `bin` and
  `scripts` first on `PATH`.

**Data.** Every case starts from a fresh copy, named `bgal.simple`, of
`/data2/bgal_addon/control/bgal_addon.simple` (beta-galactosidase after 2D
classification: 5,513 particles, box 256 at 1.275 Å per pixel, 50 selected
class averages), or from a project kept in `scratch/keep/`. The originals
are never written. `scratch/run_case.sh <case>` holds the exact command
lines and prepares each run directory, as in that run.

**Cases.** Distributed runs use `nparts=4 nthr=6`. Cases run one at a time on
the otherwise idle machine.

| Case | Starts from | Command line | Exercises |
| --- | --- | --- | --- |
| S1 | the control project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6` | the start of P1 |
| S3 | the control project | S1's line with `nstates=3` | the start of P3 (`cavg_ini_ext` with `nstates=3` needs particles carrying three states, so P3 cannot start from S1) |
| P1 | S1's final project | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` | stages 4–8; the single-state table in stages 4–5 (`refine=prob`), `prob_neigh` in 6–8 |
| P3 | S3's final project | P1's line with `nstates=3` | stages 4–5, both `refine=prob`: the multi-state table and refinement |
| T | a copy of P3's final project | `simple_exec prg=refine3D projfile=bgal.simple vol1=… vol2=… vol3=… mskdiam=180 pgrp=d2 nstates=3 refine=prob objfun=euclid update_frac=0.2 maxits=1 box_crop=… lp=… trs=… nspace=… nthr=N mkdir=no`, shared memory, N = 1, 2, 4, 8, 16 | the table fill's thread scaling |
| W | | `simulated_workflow_1jxy`, `simulated_workflow_6vxx`, `solve3D_addon`, each alone: `ctest -R '^<name>$' -V` in the Release build directory | ground-truth errors of single-state `solve3D` |

- **S3.** If S3 leaves a state empty, P3 cannot start (`cavg_ini_ext` needs
  every state populated). Rerun S3 once with `reseed_states=yes` and record
  it; if a state is still empty, stop.
- **T.** T takes `vol1` to `vol3` from P3's final state volumes. It takes
  `box_crop`, `lp` and `trs` from P3's last stage, which its
  `solve3D_manifest.txt` records (crop, low-pass, shift limit), and that
  stage's `nspace` from `simple_solve3D_controller.f90`. It runs in shared
  memory because `prob_align` without `nparts` fills the table in process,
  so `PROB_TAB_BENCH` holds the whole fill. `mkdir=no` runs it in place on
  the copy.
- **W.** W's three entries run single-state `solve3D`, so they reach only
  Phase 2's key, and Phase 4's loop with it. Phase 0 confirms from their logs
  that they reach `refine=prob`. If Phase 0 finds a high-level entry that
  runs the multi-state table, it is recorded and run in Phases 1 and 3 under
  W's rule. Otherwise those phases are judged on P3 alone, and the report
  says so.

**Measures per case.** Record the command line, the wall time, the final FSC
0.143 resolution per state, and the peak memory. From the benchmark files,
per prob-stage iteration, record:
- the fill region's wall time;
- the thread-seconds of seed, score loop and refinement;
- the load balance, E, the alternations and the grid evaluations;
- with a key on, the fits and fallbacks.

**Acceptance of a repeated case.** A case's Phase 0 value is the mean of its
two Phase 0 runs, and its spread is their difference. W runs once, so its
Phase 0 value is that run.
- P1 and P3: the final FSC 0.143 resolution of every state is no worse than
  its Phase 0 value by more than the larger of the spread and 3%.
- W: `simulated_workflow_1jxy` and `simulated_workflow_6vxx` pass, with
  whole-volume correlation against the truth ≥ 0.97 × Phase 0's and truth
  FSC 0.143 resolution ≤ 1.03 × Phase 0's.
  - `solve3D_addon` is recorded, not judged: its pose-error excess varies by
    about 1° between runs of one build (`global_sigma2_only.md`, section 13).
  - A W entry that fails on the base build is recorded and left out, not a
    stop.
- Control runs and key-off runs are recorded and must meet acceptance; they
  are not judged for speed.

**Faster.** A change counts as faster when all of these hold:
- its target measure is below the Phase 0 value by more than the Phase 0
  spread;
- that drop is at least 5%;
- the case's wall time is not longer than the longer Phase 0 run.

A target measure is thread-seconds summed over the parts, averaged per
iteration within each `refine=prob` stage, and summed over those stages
with Phase 0's iteration counts as weights. That way a run that stops
earlier or later does not count as faster or slower.

The targets:
- Phases 1 and 3: P3's refinement thread-seconds.
- Phase 2: P1's score-loop thread-seconds.
- Phase 4: P3's score-loop thread-seconds.

**A rejected change.** A change that fails acceptance or is not faster is
not kept:
1. the phase saves it with `git diff` (read-only) as
   `scratch/keep/phase_<N>_rejected.diff`;
2. it reverts the change and its key with
   `patch -R -p1 < scratch/keep/phase_<N>_rejected.diff`, since git commands
   that write are refused;
3. it records the numbers and the reason in section 11;
4. the run continues with the next phase.

A phase that comes later and needs code from a rejected phase (the shifted
memo builder) writes that code itself.

**Stop and write a stop report** when:
- the base build, the fast gate, or S1, S3, P1, P3 or T fails on unchanged
  code;
- an identity test fails for a reason other than a bug in the new code
  (never loosen its tolerance);
- a change needs a file outside the file table;
- a case cannot run as written.

**Gates** in every phase that changes source:
- `./compile_debug.sh` (the Debug build and the fast gate, every unit suite);
- the Release build;
- `python3 scripts/check_test_registry.py .` whenever tests or their
  registration change.

**Not in this run:**
- the neighbourhood table and the product form;
- the generated indexes under `doc/code_overview/`;
- the skill file (its text goes into the report);
- the ribosome timing on the production server, which is the maintainer's;
- promoting or removing a kept key, which is the maintainer's after review.

### 10.1 Phase 0: benchmark files, baseline, thread scaling

1. Find how SIMPLE builds on the Dell; build the base, Debug with the fast
   gate and Release.
2. Write the benchmark files of 3.4 and the alternation counter, the only
   source changes of this phase. Use SIMPLE's conventions (`L_BENCH_GLOB`;
   `tic`/`toc` outside parallel regions only). Rebuild both.
3. Run S1 and S3 once each and keep them; P1 and P3 twice each (the spread);
   each W entry once; T at the five thread counts.
4. Record in section 11:
   - the measures of every case;
   - the score loop's efficiency at 16 threads in T (its thread-seconds at 1
     thread divided by those at 16);
   - P1's score-loop thread-seconds per reference;
   - which W entries reach `refine=prob`, and with how many states.

Exit:
- the builds and the fast gate pass;
- the benchmark files appear in every prob-stage iteration with `% accounted
  for` above 90%;
- every case is recorded;
- S1, S3 and P3's final project are kept.

### 10.2 Phase 1: particle-frame grid (`prob_grid_ptcl`)

Source:
- the key (section 5);
- the shifted memo storage, builder and kernel variants;
- the grid option and the converted start (3.1(a));
- the tests of section 6.

Runs:
- P3 twice with `prob_grid_ptcl=yes` and every other key off (with W's
  multi-state entry, if Phase 0 found one);
- P3 once with every key off;
- P1 once as a control, unchanged by construction.

Exit: the gates; the identity, clamping, recovery and option-off tests pass;
acceptance holds; the change is faster. Otherwise it is rejected.

### 10.3 Phase 2: particle-frame seed (`prob_seed_ptcl`)

Source:
- the key;
- the one-entry bank on the shifted memo builder (3.1(b));
- the converted refinement start;
- the recording of u_s in `ref_assign`;
- the tests of section 6, including the new 3D table tester and its
  registration.

Runs, with `prob_seed_ptcl=yes` and every other key off: P1 twice and W
once. Then P1 once with every key off. P3 once as a control.

Exit: the gates; the seed tests pass; acceptance holds; the change is faster.
Otherwise it is rejected.

### 10.4 Phase 3: quadratic step (`prob_grid_fit`)

Source:
- the key;
- steps 1 to 6, the fallback and the counters (3.2);
- the tests of section 6.

Runs:
- P3 twice with `prob_grid_fit=yes` and every other key off (with W's
  multi-state entry, if Phase 0 found one). These are the runs that judge the
  step.
- If Phase 1 kept its key: P3 once with both keys on. This checks that the
  two combine and is not judged.
- Record the fallback fractions per reason.

If the fit alone fails the case acceptance, build the second form of 3.2
and run the same cases. If that fails too, the step is rejected.

Exit: the gates; the guard and accuracy tests pass; acceptance holds; the
change is faster. Otherwise it is rejected.

### 10.5 Phase 4: blocked table loop (conditional)

Phase 4 runs only if either holds:
- Phase 0's efficiency at 16 threads is below 0.75;
- after Phase 3, P3's score-loop thread-seconds exceed its refinement
  thread-seconds. Use P3's last run with every kept key on, or Phase 0's runs
  if none was kept.

Otherwise it is skipped, with the reason in section 11.

Source:
- block-local storage and the blocked loop in `fill_tab_range`, for the
  scope of 3.3, with B a module parameter;
- if Phase 2 kept its key, seed storage for B entries;
- the blocked-loop case in the 3D table tester (section 6).

Runs:
- T at 16 threads with B = 8, 16 and 32, every key off; choose the fastest B.
- P3 twice with every key off.
- If Phase 2 kept its key: P1 and W once with `prob_seed_ptcl=yes`.

Exit: the gates; the stream-identity test passes; acceptance holds; the
change is faster. Otherwise it is rejected.

### 10.6 Phase 5: documents and close

1. Write the documents of section 7, except the skill file.
2. Run P1 and P3 once more on the final build: once with every kept key on,
   once with every key off. These runs are recorded and must meet
   acceptance.
3. Write the report, `doc/refactoring_notes/completed/prob_search_speedup_report.md`.
   It reads on its own and covers:
   - what changed and why;
   - the evidence against Phase 0, in tables;
   - every decision and where it is recorded;
   - the rejected changes and their saved diffs;
   - the proposed skill text;
   - what is left open.
4. Move this plan beside the report.

### File table

| File | Change | Phases |
| --- | --- | --- |
| `src/main/commanders/simple/simple_commanders_prob.f90` | Benchmark files of `exec_prob_tab` and `exec_prob_align`; wall timers outside the parallel regions (0); the fit, fallback and block lines (1, 3, 4). | 0 1 3 4 |
| `src/defs/simple_refine3D_fnames.f90` | `prob_tab_bench_fname`, `prob_align_bench_fname`. | 0 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | Thread-seconds and counters (0); the grid option and the step under their keys, multi-state only (1, 3); the seed in the particle frame, including u_s in `ref_assign` (2); the blocked loop with block-local storage (4). | 0 1 2 3 4 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | Alternation counter (0); the grid option, its bank and the converted start (1); the quadratic step, its guards, fallback and counters (3). | 0 1 3 |
| `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`, `simple_polarft_memo.f90`, `simple_polarft_corr.f90`, `simple_polarft_geom.f90` | Per-thread storage and builder of shifted particle memos; the kernel variants taking such a memo (Euclid and CC); seed storage for B entries (4). | 1 2 4 |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90` | `prob_grid_ptcl`, `prob_seed_ptcl`, `prob_grid_fit`. | 1 2 3 |
| `src/main/ui/simple/simple_ui_refine3D.f90`, `simple_ui_solve3D.f90` | Developer entries for the keys. | 1 2 3 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | The grid, seed, guard and accuracy cases. | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` (new) | The assignment case (2); the blocked-loop stream case (4). | 2 4 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | "3D probability table" in `suites_pftc_align2D3D`. | 2 4 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | The new sub-suite in the `suite=` help of `unit_pftc_align2D3D`. | 2 4 |
| `doc/policies/3D/refine3D_policy.md` | Section 7. | 1 2 3 5 |
| `doc/refactoring_notes/completed/prob_search_speedup.md`, `prob_search_speedup_report.md` | The plan moved on completion, and the report. | 5 |

## 11. Progress

### 14 October 2026, written without a build (interactive session, no Dell)

Source written by Claude from a laptop session that could not compile or run
anything (no Fortran compiler, no data). The maintainer built it, ran the
cases below and reported back (see "Measured results" at the end of this
entry). The unit gate was not re-run after the CC reference fix in
`test_particle_bank` (see "Open"). Nothing is committed.

- **Phase 0 (source written, unbuilt).** Built on `PROB_TAB_BENCH`, which commit
  `77a8df9a8` (particle I/O) had already added with the wall-clock entries; the
  note below that no `prob_tab_bench` names exist is stale. Added: the
  thread-seconds block, the context counters and the load balance
  (`prob_fill_profile` in `simple_eul_prob_tab_utils.f90`, filled by
  `fill_tab_range`, read through `get_fill_profile` before `kill`); the
  alternation, grid-point and shift-search counters in `pftc_shsrch_grad`
  (`get_profile` takes them as optional arguments); `PROB_ALIGN_BENCH`
  (`prob_align_bench_fname`, `write_prob_align_bench`). Thread timing uses
  `omp_get_wtime` into thread-indexed accumulators.
- **Phase 3 (source written, unbuilt).** `prob_grid_fit`, 3.2 as specified:
  `try_grid_fit`, the pure `fit_quadratic_step` and `grid_stencil_exists`,
  counters per reason, key in `simple_parameters*.f90`, developer entries on
  refine3D and solve3D, bench lines, tests `test_grid_fit_step` and
  `test_grid_fit_search` in `simple_pftc_inplane_tester.f90`. Choice not fixed
  by the plan: the quadratic is fitted to the loss, `-log(objective)` for
  Euclid (the objective is `exp(-loss)`) and `-objective` for cc. The closed
  form was checked against `numpy.linalg.lstsq`. The second form of 3.2 is not
  built. Measured: see "Measured results".
- **Phases 1 and 2 (source written, unbuilt): the identity of 3.1 is exact only
  with the Nyquist bin.** A numpy emulation of the kernels shows that the
  particle-side bank entry is `X*CTF*conj(shmat4aln(t))` with
  `t = s * rotmat2d(angtab(irot))`, the conversion `sh_rot` already uses, and
  that with the angular Nyquist bin included on both sides it equals the
  reference-frame loss at every rotation to 1e-16. The Euclid kernels sum
  `p = 1..pftsz` into a `pftsz+1` spectrum and so drop that bin, which differs
  between the two formulations: about 8% of `1+|L|` for white-noise transforms,
  about 2e-4 for a smooth angular spectrum. Ruling (maintainer, 14 October
  2026): include the Nyquist bin in the new kernels (exact against the true
  shifted loss, slightly different from today's shifted-reference kernel).
  Written to that ruling:
  - `polarft_calc`: `memoize_shifted_ptcl` (shifted `X*CTF` transformed over
    rotation, written into a caller-owned array in the calling thread's FFT
    buffers), `gen_objfun_vals_bank` and `gen_likelihood_val_bank` (Euclid and
    CC, `p = 1..pftsz+1`). `SHERRSQ` moved to the module so the table can use it.
    Deviation from the file table: the bank is owned by its user (the search
    object, the table's per-thread seed array), not by the calculator, so
    nothing is added to `new`/`kill`.
  - `pftc_shsrch_grad`: `grid_ptcl` (26 entries: 25 lattice points and a scratch
    entry for the fit), built lazily when the particle changes, the start
    converted by `shift_ptcl2ref` and clamped into the box, the quadratic step
    on the lattice (the fitted lattice shift is converted at its rotation,
    scored from the scratch entry, and the box test is made at the rotation it
    is kept at), `get_bank_builds`.
  - `eul_prob_tab`: `prob_grid_ptcl` passed to the search objects in the
    multi-state table; `prob_seed_ptcl`: the seed converted with
    `shift_ref2ptcl` at the rotation `estimate_shift_seed` ends at, one bank
    entry per particle when above `SHERRSQ`, the refined start converted per
    reference and clamped, `record_seed_shift` for the unrefined assigned
    candidate; bench lines for the keys and the memos built.
  - Tests: `test_shift_frames`, `test_particle_bank` (every rotation against the
    direct sum, the shifted-reference kernel up to the Nyquist bin at 5e-2, the
    likelihood route, a seed at its own rotation) and `test_grid_ptcl_search`
    in the in-plane tester; `simple_eul_prob_tab_tester` ("3D probability
    table") for `record_seed_shift`. Not written: a table-level test of the
    seed scoring (needs a project).
  - Measured: see "Measured results".
- **Phase 4 (not written; skipped, 14 October 2026).** Its criterion was not
  met (see "Measured results"), and the one measurement that could have
  justified it says the score loop is compute-bound. Same two-state case with
  both keys, one part of 6 threads (`nparts=1`, 13_refine3D): `score loop
  s/(particle*ref)` 1.93e-5 s, against 2.44e-5 s with 4 parts of 6 threads
  running together (24 threads; the master reports 20). The 26% difference is
  the size of the oversubscription (24 threads on 20), so little or none of it is
  memory contention, and blocking could recover at most about 6 of the score
  loop's 28 thread-s (about 1 s of 12 s per part). Rejected for the cost: block
  storage and a new candidate path in the 3D table for at most that. Revisit
  only if the loop is ever measured memory-bound (efficiency below 0.75 at 16
  threads on a machine without oversubscription).

- **`prob_fast` (source written, unbuilt; ruling 7).** The only key: the three
  separate keys are removed from the parameters, the registry and the refine3D
  and solve3D interfaces, and `fill_tab_range`, `ref_assign` and the bench
  writer read `prob_fast` (grid and fit when `nstates > 1`, seed when the
  single-state table has a shift search). The measured runs above were made
  with the separate keys: `prob_seed_ptcl=yes` is now `prob_fast=yes` in a
  single-state run, `prob_grid_ptcl=yes prob_grid_fit=yes` is `prob_fast=yes`
  in a multi-state run; `prob_grid_ptcl` alone has no equivalent any more.

#### Measured results (14 October 2026, bgal bench data, commit `5e1e6aee` plus the working tree)

Both runs on the maintainer's machine, 6 threads per part, 4 parts. Keys not
named were off.

**Single state, `solve3D`, keys off against `prob_seed_ptcl=yes`.**
- Whole run 631.1 s to 505.9 s. FSC 0.5 / 0.143 at 4.29 / 3.75 A in both.
- Table pre-steps over iterations 1-8: 198.1 s to 82.0 s. The score loop fell
  from 3462 to 907 thread-s (-74%). Seed estimation and refinement are
  unchanged and gradient evaluations per search are about equal. The fraction
  of memos sampled is 0.75-1.0.
- `prob_neigh` time is unchanged (117 s of table time, summed over the slowest
  part of each iteration). Where the 505.8 s goes, from `SOLVE3D_BENCH` and the
  per-iteration `REFINE3D_STRATEGY_BENCH` (base run 631.0 s in brackets):

  | | modified | [base] |
  |---|---|---|
  | pre-step, `prob` iterations 1-8 | 82.2 s | [198.1] |
  | pre-step, `prob_neigh` iterations 9-17 | 159.9 s | [159.7] |
  | refine3D parts and merge (search, per-iteration partial reconstruction) | 142.6 s | [150.3] |
  | strategy setup, includes the reprojection model | 39.0 s | [39.6] |
  | assembly and postprocess inside the iterations | 26.0 s | [26.2] |
  | sigma2 commit | 0.4 s | [0.4] |
  | refine3D stages, sum of the rows above | 450.1 s | [574.3] |
  | final reconstruction, includes the sigma pass of iteration 19 (7.8 s) | 48.5 s | [48.4] |
  | start state and epilogue | 7.0 s | [7.9] |

  The `prob_neigh` pre-step is the largest block left (32% of the run) and no key
  touches it. Of its 160 s, 117 s is table time; the rest is process launch and
  the table reads. In iterations 12-17 (`nspace` 5000) a part spends 3.5 s of
  15.5 s in startup/setup and 11.1 s in the fill. Whether the particle-frame
  grid applies to the `prob_neigh` table has not been looked at.

**Two states, `refine3D`, 2 iterations, `update_frac=0.2`, `lp=8`, `pgrp=d2`,
`vol1`/`vol2` the even and odd half-maps of one single-state volume** (5000
references, 910 sampled particles, 138 refined references per particle).

| | keys off | `prob_grid_ptcl` | `+ prob_grid_fit` |
|---|---|---|---|
| table step, iterations 1 / 2 | 34.4 / 33.5 s | 17.3 / 18.5 s | 16.3 / 15.3 s |
| table fill, wall (part 1, iteration 1) | 28.6 s | 14.7 s | 12.7 s |
| refinement, thread-s | 141 | 53 | 40 |
| score loop, thread-s | 27.0 | 27.1 | 27.9 |
| gradient evaluations per search | 12.2 | 11.9 | 6.1 |
| load balance | .982 | .916 | .904 |
| whole run | 96.9 s | not logged | 60.4 s |

- `prob_grid_fit`: the fit was kept for 15206 of 31602 searches (48%). Fallbacks:
  edge 13114 (41%), step 2073 (6.6%), value 1144 (3.6%), curvature 65 (0.2%),
  search box 0. Alternations per search fell from 1.00 to 0.52.
- Search quality, runs not paired (the 20% sample differs between runs):
  score 0.637 / 0.640 / 0.644 (iteration 1) and 0.548 / 0.553 / 0.556
  (iteration 2); state populations 814/96, 815/95, 815/95 and 1632/188,
  1632/188, 1643/177; state 1 FSC 0.5 at 10.2 A in iteration 2 in all three.
- State 2 has about 95-188 particles and a cFAR of 0.001-0.003. Its resolution
  at FSC 0.143 varies between runs (13.6, 21.8 and 23.3 A in iteration 1; 14.8,
  29.7 and 20.4 A in iteration 2) and says nothing about the keys until the
  spread of the keys-off run itself is known.
- "Continuous in-plane runs improved" rose from 67-71% to 80% and 86% with
  `prob_grid_fit`. Not explained. Hypothesis: the fit leaves the table's shift
  slightly less polished and the post-assignment in-plane route finds more to
  improve. The score is unchanged.
- What is left, with every key on (part 1, iteration 1, wall 12.7 s of the
  fill): the score loop is 28 of 68 thread-s (41%) and refinement 40 (59%). A
  search whose fit is kept costs about no gradient evaluations and a fallback
  about 12, so the 41% edge fallbacks are about half the refinement left. The
  load balance fell from 0.98 to 0.90-0.92 (about 1.2 s idle at the tail of a
  14.6 s region).
- Phase 4 is not triggered by its own criterion (score loop 28 thread-s is
  below refinement 40 thread-s; the 16-thread efficiency was never measured).
  The maintainer asked for the multi-state gain regardless. The deciding
  measurement is cheap: `score loop s/(particle*ref)` at 24 concurrent threads
  is 2.4e-5 s; at 6 threads (`nparts=1`) a lower figure means the loop competes
  for memory and blocking can help, an equal one means it is compute-bound and
  blocking cannot.
- **Tried and rejected: a stencil shifted inside the lattice for an edge grid
  minimum** (14 October 2026; reverted, the diff is saved in
  `prob_search_speedup_edge_stencil_rejected.patch`). Estimate: half of the 41%
  edge fallbacks would convert. Measured (14_refine3D against 12_refine3D, same
  case, both keys): about 14% convert (1816 of 12854 edge cases in part 1,
  iteration 1); the rest fall back on curvature (56 to 4537) or step (2046 to
  6448) instead. Gradient evaluations per search 6.1 to 5.1, fill wall within
  the part-to-part noise (12.15 / 11.58 s against 11.64 / 11.89 s, slowest part
  of iterations 1 / 2), whole run 58.9 s against 58.7 s, score and populations
  unchanged. Unit gate: a kept edge fit was 0.23 px from L-BFGS-B (bound 0.2).
  Reading: an edge grid minimum mostly means the optimum lies beyond the
  lattice, where the loss keeps falling outward and no quadratic fits, so the
  fallbacks are genuine. The edge fallback stays; the fit's remaining cost is
  the fallbacks' L-BFGS-B runs.
- Written, unbuilt: `fill_tab_range` runs the shifted fill with
  `schedule(dynamic)` when `prob_grid_ptcl` or `prob_grid_fit` is on (keys off:
  static, unchanged). Expected: the idle tail falls from about 1.2 s to the
  cost of one particle (0.35 thread-s). The random in-plane picks are drawn in
  a different order per thread, as for any change of assignment of particles to
  threads.
  Measured (same two-state case, both keys, 12_refine3D against 10_refine3D):
  load balance 0.84-0.97 to 0.984-0.989 in all eight parts; slowest part's fill
  wall 12.76 / 12.82 s to 12.24 / 11.67 s (iterations 1 / 2); table step
  16.3 / 15.3 s to 15.2 / 14.4 s; whole run 60.4 s to 58.9 s. Score 0.649 /
  0.551 (against 0.644 / 0.556), populations 808/102 and 1617/203.

**Open.**
- Repeat the keys-off two-state run once to measure its own spread (state 2).
- `unit_pftc_align2D3D`: `test_particle_bank` with the CC reference fix passes
  (differences 2e-7), confirmed on 14 October 2026. The only failure in that
  run was the 0.2 px bound for the edge-stencil experiment above, since
  reverted; the gate has to be run once more on the reverted tree.
- A state left empty by volassemble on iteration 1 stopped the run with
  `fsc_state02.bin does not exist`. Found when identical `vol1`/`vol2` tied in
  `ref_assign` (`minloc` picks state 1) and state 2 got no sampled particles.
  Cause: `determine_dropped_states` / `carry_forward_dropped_state` in
  `simple_commanders_rec_distr.f90` write nothing, so a state with no earlier
  volume or FSC has nothing to carry, while the driver
  (`distr_execute_iteration` in `simple_refine3D_strategy.f90`) takes the
  project's state population as the test, and with `update_frac<1` the
  unsampled particles keep their state-2 labels, so state 2 counts as populated
  and its FSC is read. Existing behaviour, not caused by these keys.
  Fixed (written, unbuilt, not run): in the driver a populated state whose map or
  FSC file does not exist keeps its input volume as the reference (`vol<state>`
  in the command line and job description untouched), logs
  `>>> STATE n HAS NO RECONSTRUCTION AND NO EARLIER MAP; KEEPING ITS INPUT VOLUME
  AS REFERENCE`, and skips the FSC plot, the registration in `out` and the
  postprocess; the shared-memory `inmem_finalize_run` registers a state only when
  both files exist. Test: the original command (`vol1=vol2=recvol_state01.mrc`,
  `nstates=2`, `maxits=2`) must now end normally with the log line in each
  iteration.

Original note, before any phase ran. Verified against the code on 14 October 2026:
- no `prob_tab_bench` or `prob_align_bench` names exist;
- `coarse_search_opt_angle` still shifts the reference through
  `gen_objfun_vals` (the fixed-angle `coarse_search` uses
  `gen_corr_for_rot_8`);
- the search counters have `get_profile` and `reset_profile`, with no
  caller;
- `rt_prob` is the pre-step's only timing;
- solve3D's stage controller runs `refine=prob` in stages 3–5.

### 14 October 2026, the `prob_neigh` pre-step (single-state solve3D, iterations 9-17)

The largest block left in the single-state run (160 s of 506 s; section 11,
"Measured results"). Late solve3D stages run `refine=prob_neigh` with
`prob_neigh_mode=state` (`fill_tab_subspace_range`).

**Found: an O(nspace²) toolbox build in every process.** `build_general_tbox`
with `l_neigh` builds the subspace map from a table of `sym_dists` over every
(subspace, full) direction pair, with `nspace_sub` scaled with `nspace`
(126 per 2500, `apply_refine3D_search_overrides`): 252 x 5000 pairs at
`nspace` 5000, each forming `ori` objects and `nsym` symmetry copies. The
bench files show it as the `startup/setup` of every `prob_tab_neigh` part:
0.85 s at `nspace` 2500, 3.5 s at 5000 (x4.1 for x2 directions). The same build
runs in the `prob_align_neigh` master and in every refine3D part of those
iterations (the refine3D scheduler time grows by 2.7 s from 2500 to 5000), so
an iteration at 5000 pays it about three times on its critical path, roughly
10 s, and 60-70 s of the 506 s run. It grows with the symmetry order (d2 has 4
operations; I has 60).

- **Written, unbuilt (default path, no key; ruling 8):** `sym%sym_dist_table`
  forms the symmetry copies, projection normals and in-plane unit vectors once
  per orientation (in parallel; 60 copies per orientation in I) and the table
  from dot products, with the operands, branches and tie folding of
  `sym_dists_1`; the builder calls it instead of the pair loop. Cost at
  `nspace` 5000 in d2: 20,000 symmetry copies and 5 million dot products,
  against 5 million `ori` symmetry copies before. Test `test_sym_dist_table` (symmetry suite, `unit_ori`): table against
  `sym_dists` pair by pair in all nine point groups of the tester (1e-4 deg)
  and the same nearest neighbours in both directions. Possible difference: a
  last-bit difference in a dot product could move a full direction to an
  equally near subspace; the test checks the maps.
- Expected: `startup/setup` of the neigh parts from 3.5 s to a fraction of a
  second at 5000; the neigh pre-step and the refine3D scheduler of iterations
  12-17 each shorter by about 3 s.

**Instrumented (written, unbuilt): the neighbourhood fill's split.** The
`PROB_TAB_NEIGH_BENCH` files had no thread-seconds. Both neighbourhood fills
(`fill_tab_subspace_range`, `fill_tab_stoch_range`) now time the seed, the
score loop (the subspace representatives and the pooled neighbourhood) and the
refinement per thread, count the references scored (`nscored`, a new field of
`prob_fill_profile`), and pass the search objects' counters
(`record_neigh_fill_profile`); `exec_prob_tab_neigh` hands the profile to the
bench writer, which prints `scored references (all)` and the score cost per
scored reference for a neighbourhood table.

**Not done yet: the particle-frame seed in the neighbourhood table.** In the
single-state `state` mode every score of a particle (252 representatives and
the pooled neighbourhood) uses the same shift seed, the pattern
`prob_seed_ptcl` exploits in the full table. Worth building only if the score
loop is a large part of the neigh fill; the instrumented run answers that.

## 12. Code map

After the global sigma2 commit `2307ea96f`, the kernels' shell weight is
`sigma2_wk(k,ih)` / `sigma2_wk8(k,ih)`, two global spectra formed by
`memoize_sigma2_weights`. `ft_ctf2 · ft_ref2` appears in
`gen_raw_euclid_vals_impl`, `gen_euclid_crvec` and
`gen_raw_euclid_grad_at_angle`.

| Area | Files |
| --- | --- |
| Pre-step orchestration | `src/main/strategies/parallelization/simple_refine3D_strategy.f90` (`rt_prob`), `src/main/commanders/simple/simple_commanders_prob.f90` (`exec_prob_align`, `exec_prob_tab`) |
| Stage modes | `src/main/solve/simple_solve3D_controller.f90` (`set_refine3D_mode_policy`) |
| The table | `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` (`fill_tab_range`, `estimate_shift_seed`, `refine_best_refs`, `ref_assign`, `read_tab_to_glob`), `simple_eul_prob_tab_utils.f90` (`sample_likelihood_dist`, `calc_num2sample`, `prob_candidate_buffer`, `materialize_seed_shift`) |
| The assignment's consumer | `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90` (`srch_prob`) |
| The shift search | `src/main/pftc/simple_pftc_shsrch_grad.f90` (`grad_shsrch_minimize` and its `sh_rot` conversion, `coarse_search_opt_angle`, `get_profile`) |
| The kernels | `src/main/pftc/simple_polarft_corr.f90` (`gen_raw_euclid_vals_impl`, `gen_euclid_crvec`, `gen_corrs`, `gen_prob_likelihood_euclid_val`, `gen_prob_likelihood_cc_val`), `simple_polarft_geom.f90` (`gen_shmat4aln_8`), `simple_polarft_memo.f90` (`memoize_ptcls`, `memoize_refs`) |
| Test registration | `src/main/commanders/test/simple_commanders_test_class.f90` (`suites_pftc_align2D3D`), `src/main/ui/simple_test/simple_test_ui_class.f90`, `scripts/check_test_registry.py` |
| Precedents | `prob_athres_taper` (commit `cfbc56bf3`) for a key; `doc/refactoring_notes/completed/global_sigma2_only.md` section 13 for the Dell build, data and command lines |

## 13. Review of 14 October 2026

The review found six problems in the first version. Each was checked against
the code and accepted.

| Finding | Resolution |
| --- | --- |
| The particle-side grid is not exact across rotations (the reference is shifted in its own frame before the FFT over rotations). | 3.1: a different grid behind a key, with a pointwise identity test and validation. |
| The matcher cannot replace the table's refinement (`ref_assign` uses the refined distances; `srch_prob` treats the assignment as given). | 3.2 makes the refinement cheaper; the option to drop it is withdrawn. |
| Blocking breaks `append_or_replace` and the reader's `nrefs` check. | 3.3: block-local dense storage and particle-major emission. |
| The benchmark mixed thread-seconds with wall time. | 3.4: wall entries outside parallel regions; thread-seconds apart. |
| The quadratic fit had no guards. | 3.2: six steps with a defined fallback, each tested. |
| Approximate changes need default-off keys. | One key per change (section 5). |

Found in the same re-check:
- the counters already have getters;
- the single-state table pays a shifted reference per pair (section 2,
  Phase 2);
- under patch CTF, ruling 4's premise fails;
- the multi-state neighbourhood table also grids per refined reference (out
  of scope).
