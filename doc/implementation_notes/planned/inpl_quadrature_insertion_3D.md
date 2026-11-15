# In-plane quadrature insertion in 3D

Status, 2026-11-15: delivered; the key removed and the insertion made the default (Hans, decision 1). Was, 2026-11-14:
delivered, not validated at scale; the run's report is
[inpl_quadrature_insertion_3D_report.md](inpl_quadrature_insertion_3D_report.md) and the maintainer's runs of
section 7 follow. Both design reviews incorporated, including the per-view window-prior centre correction (section
10). Base:
`06d3a4de9` (first written at `b2a9d544f`; the three commits between them touch no file of the file table except
an unrelated paragraph of the solve2D policy). Supersedes: none.
Depends on: the 2D in-plane quadrature (`inpl_quadrature`, `prob_insert_inpl`, `INPL_POSTERIOR_BUDGET`;
`d6cf4ae25` to `b2a9d544f`; validated: yes, Hans's beta-galactosidase and TRPV1 class averages, 2026-11-13);
the 3D weighted insertion (`topk`, `top_posterior_extras`, `prob_insert_cone`, `prob_insert_poses`; validated:
yes, production, default `topk=3` since 2026-11-01); the continuous in-plane solve of the probabilistic route
(`refine_assignment_continuously`, `refine_extra_pose_continuously`; validated: yes, production).
Companion: [the 2D plan](inpl_quadrature_insertion_2D.md). Its rulings on the row, the temperature, the phase
and the continuous angle are applied here to every pose a particle is inserted at; they are restated below so
that this note reads on its own.

## Terms

- **Row (table row).** A particle's losses against every reference of the probability table (`eul_prob_tab`,
  `eul_prob_tab_neigh`), one in-plane index and one shift per reference.
- **Committed pose.** The pose the assignment draws from the particle's row at `prob_beta` times the physical
  temperature and the search commits to the project after its continuous solve. It is the hypothesis the next
  iteration starts from.
- **Extras (weighted insertion, `topk`).** Up to `topk - 1` other refined entries of the committed state's block,
  within the insertion cone of the committed pose, of highest posterior at the physical temperature
  (`top_posterior_extras`). The particle enters the iteration's reconstruction at the committed pose and at its
  extras, with weights summing to 1. The final and bootstrap reconstructions use the committed pose alone.
- **Insertion cone.** The poses whose rotation from the committed pose (smallest over the point group) times
  the mask radius, plus the shift difference, move the rim by at most the finer of the map's FSC=0.5
  resolution and the matching low-pass limit (`insert_cone_within`).
- **In-plane row.** The Euclidean losses `L_r`, `r = 1..nrots`, of one reference against the particle at one
  shift over every grid rotation (`gen_euclid_dists`, `gen_euclid_dists_bank`).
- **Physical temperature `T`.** `posterior_temperature(wsqsum, pftsz) = pi wsqsum / (2 pftsz)`: the posterior
  over a row is `p_r` proportional to `exp(-T (L_r - L_min))`. The table's draws run at `prob_beta` times `T`
  (`ptcl_posterior_temp`). The extras' weights and the 2D quadrature use `T` itself.
- **Partition sum `Z`.** `Z = sum_r exp(-T (L_r - L_min))` over the valid rotations of an in-plane row. `Z = 1`
  for a row with one plausible angle, about `n` for a row with `n` equally plausible angles.
- **View weight `m_k`, `m_c`.** The share of a particle that inserted pose `k` (its direction, all angles
  together) receives under the key: the in-plane marginal of decision 4, normalised over the particle's poses.
  `m_c` is the committed pose's.
- **Quadrature.** The `K` rotations holding the points `u_j = (j - 1 + v) / K`, `j = 1..K`, of the row's cumulative
  posterior, each at weight `1/K`, coincident points merged (`inpl_quadrature`). `v` is a keyed per-particle phase.
  A row whose mass sits in one rotation collapses to one entry.
- **Solved pose.** A pose after its continuous joint solve of shift and fractional in-plane angle.
- **`n_eff`.** `exp(entropy)` of a distribution: the effective number of entries it spreads over.
- **`m_max`.** `ceiling(pi mskdiam / lp_band)`, the highest angular harmonic the band resolves at the rim; the 2D
  readout prints it beside `K`.

## 1. Why

The 2D quadrature was built against radial artifacts in top-view class averages of small membrane proteins.
With `K = 64`, Hans called the beta-galactosidase averages (5,513 particles, 90 classes) almost perfect and the
TRPV1 averages perfect (2026-11-13), and the next day made it a fixed part of every probabilistic 2D pass
(`b2a9d544f`). The run's readout showed 1.1 to 1.8 entries per particle, 40 to 90 % of rows collapsed and mean
`n_eff` 1.0 to 1.4 (2D report, section 8). The loop it acts on exists in 3D. A particle whose in-plane angle the
data do not decide enters its half-map at an angle its own noise chose, and the half then favours that angle
for that particle. This self-reinforcement was measured on PfCRT with `prob_post_diag=yes` (pose-marginalization
note, section 12, 2026-10-22, at an older base). The own half's advantage at the inserted pose averaged
+0.81 nats in stage 3, +2.30 in stage 4 and +4.71 in stage 5, against a control of 0.00 ± 0.03. It was largest on
the sparsely populated views: in stage 5, 5.5 nats in the lowest population quartile and 3.3 in the highest.
Those are the views a small membrane protein is short of. The 3D weighted insertion does not spread the
angle. Each direction of a row carries one in-plane angle, the extras are other directions, and the cone tests
the whole pose, angle included. This plan carries the 2D quadrature to every pose a 3D particle is inserted at.
Whether it improves 3D maps is Hans's test after delivery (decision 2). The history of the weighted insertion
before the cone (section 2, last bullet) is why the feature arrives behind a key that is off by default.

## 2. What the code check changed

Checked at `b2a9d544f` with `rg`, the review's points re-checked at `06d3a4de9`. These premises are fixed by the
code.

- **A refined entry's in-plane index is its row's best angle, not a draw.** The fill scores every reference at a
  hot-drawn rotation (`score_ref`, `gen_posterior_val` at `ptcl_posterior_temp`). It then refines the shift of the
  `projs_ns` best references with alternating search objects (`new_alternating`, `opt_angle = .true.`), and
  `grad_shsrch_minimize` re-selects the in-plane argmax over every rotation at each shift. A refined entry
  therefore holds the minimum of its in-plane row at its refined shift. Only unrefined entries hold the hot draw.
  The extras are refined entries (`select_extras` only considers them), so an extra's table weight already
  ranks its direction by its in-plane minimum. This corrects the design discussion of 2026-11-14, which assumed
  the extras' angles were random draws.
- **The cone includes the angle.** `insert_cone_within` composes the in-plane and direction matrices of both
  poses and bounds the geodesic angle, smallest over the point group.
- **Extras take the continuous route before the committed pose.** `strategy3D_prob%srch` calls
  `resolve_extra_poses` (each extra through `refine_extra_pose_continuously` from its own cell and shift, then
  `candidate_pose` and `resolve_inplane_e3`) and then `refine_assignment_continuously` for the committed pose.
  The extras keep the table's weights after their solve; the policy records Hans's reason (2026-11-01): the
  solve is local and bounded, so the losses move little and alike.
- **The extras' weights carry a shift prior** (found by the review). Under the polar shift map
  (`shift_sigma > 0`, `l_centring`), a refined entry's loss for the insertion is `insert_pose_loss`: the data term
  plus `prior_w |du|^2 / T`, `prior_w = map_prior_weight(shift_sigma, smpd) = 1 / (2 (shift_sigma/smpd)^2)`, `du`
  the increment from the stored shift (dense table, `process_particle_with_map`; neighbourhood table,
  `select_extras`). Under the window search with `shift_prior=yes` (the default; off under the map), the table's
  shift refinement minimises the data term times `exp(-w |u - u0|^2)`, `w = 1 / (2 sigma2_prior T)` with
  `sigma2_prior` the iteration's `sigma2_prior_model` (`shift_prior_sigma2`) and `u0` the actual start of that
  reference's refinement, converted back to the stored frame by `set_prior_centre`. Both tables rotate the
  particle's seed into the reference's starting-angle frame, then clip each component to `[-trs, trs]` before
  passing it as `xy_in` to `minimize`. Thus the centre is `R(irot_start)^T clip(R(irot_start) seed, -trs, trs)`,
  which can differ by view; it is zero without the single-state seed. The refined loss the extras are weighed
  by holds this penalty. An in-plane row at a fixed shift has the prior constant across its rotations:
  it cancels within the row and in `Z`, but not between poses at different shifts.
- **Shift frames.** The table's refinement returns the shift in the stored frame (`minimize` with `sh_rot`:
  `matmul(x, R)`, the `shift_ref2ptcl` conversion). The strategy hands the joint solve `R(irot)` times that
  increment (`matmul(sh, transpose(rotmat))`), the optimizer's frame at the cell. `minimize_joint` returns
  `matmul(x, R)` again (`sh_rot=.true.`). The fill memoizes the particle at the stored-frame increment
  (`memoize_shifted_ptcl` after `shift_ref2ptcl`). The row of a solved pose is therefore
  `memoize_shifted_ptcl(iptcl, sh)` with `sh` the solve's returned increment. The 2D plan proved the bank row's
  frame (`check_bank_particle_frame`). The 3D strategy's conversion is not covered by that check, so Phase 3
  asserts it.
- **Insertion.** `calc_3Drec`, through `update_state_half_rec`, inserts the committed plane at `wbase * wc`, then
  each extra (`insert_extra_poses`) from the committed plane re-shifted by the shift difference
  (`fplane_shift`, one phase ramp per extra), skipping slots of weight 0. `insert_plane_oversamp` writes every
  point-group copy of each insertion.
- **PCG.** `pcg_extra_entries` lays every extra of positive weight out as one more operator entry after the
  committed ones. `accumulate_pcg_extras` copies each extra's observation into a scratch plane, and
  `accumulate_pcg_state_half` sizes that scratch by the most extras any 500-particle batch holds
  (`MAXIMGBATCHSZ = 500`). At the base that is at most `500 * (topk - 1)` planes. With `K` entries per pose it
  would be unbounded in practice. Capping the scratch is not enough: `pcg_extra_entries` builds one `oris` record
  per entry for the whole half, and `reconstructor_pcg%prep_particles` allocates a rotation matrix, a shift, CTF
  parameters and two weights per entry, all before the first batch accumulates. **Confirmed before development
  (2026-11-14, at `06d3a4de9`): the production path reads none of that metadata after accumulation.** In
  refine3D and solve3D, shared memory and distributed alike, the matcher's `write_partial_recs` runs
  `execute_rec3D_pcg_worker`. Each (state, half) gets its own `reconstructor_pcg` in
  `accumulate_worker_state_half`: `prep_particles`, `begin_accum`, `accumulate_batch` per batch, then
  `write_raw_accum` and `kill`. `write_raw_accum` writes only the B and D accumulators, with the particle count
  and mass passed by the caller. The solve is the master's (`assemble_refine3D_pcg` →
  `execute_rec3D_pcg_distributed_master`, in both modes). It builds its own object from the raw files
  (`begin_reduction`, no `prep_particles`), finalises with `end_accum` and solves with
  `set_op_mode(PCG_OP_KERNEL)`. The only routines that read the per-entry metadata after accumulation are
  `apply_normal_matrixfree` and `apply_adjoint_all`. Both stop on an object without `prep_particles`, and
  production reaches neither: they serve the matrix-free reference operator and the one-shot `solve(y_planes)`,
  which only the PCG testers call (FLEX's `solve` is another type, `flex_pcg_t`). So the metadata can be prepared per chunk of entries. Repeated `prep_particles`
  calls leave the open accumulators untouched. `accumulate_rhs_density` adds a batch entry by entry, and each
  lattice cell is written by one thread at a time in entry order, so accumulating the same entries in chunks
  gives B and D byte for byte equal to one call.
- **Layout.** `prob_insert_poses%setup_range` allocates `nptcls * (kmax - 1)` slots at the first batch. It takes
  each batch's extra counts and weights from the table. The matcher calls it only when `topk > 1`
  (`pipeline_batch`), and the strategy sees the structure only when it has particles.
- **Readout.** `prob_readout` already holds the 2D fields `inpl_kmax`, `inpl_mmax`, `inpl_nptcls`, `inpl_npose`,
  `inpl_nsingle` and `inpl_neff_sum`. They are written, read and summed for every route. Only the 2D branch of
  `report_prob_readout` prints them (`inpl_report_counts`). The 3D tables close their assignment after the
  last batch's search (`close_assignment` in `refine3D_exec`), so the search's counts can reach the record.
- **Gate inputs.** `l_prob_align_mode` is true for `refine=prob|prob_neigh`. The neighbourhood table turns the
  extras off for `prob_neigh_mode=shc|snhc` (`open_assignment` of `eul_prob_tab_neigh`). solve3D copies its
  command line into every refine3D stage (`prep_class_command_lines`), so a solve3D key reaches each stage.
  The multi-state final assignment sets `topk=1`, `maxits=1` and `volrec=no` (`run_multistate_final_assignment`).
  Single-state solve3D stages 1-2 are `prob_neigh` `shc`, stages 3-5 `prob` and stages 6-8 `prob_neigh` `state`
  (`PROB_REFINE_STAGE = 3`, `PROB_NEIGH_REFINE_STAGE = 6`). Multi-state runs `shc` before stage 3.
- **Reusable pieces.** `inpl_quadrature` returns the entries, weights and `n_eff` but not `Z` or `L_min`.
  `RND_FAM_INPLQ` is the quadrature's draw family, and `rnd_draw_key` takes an optional `reference`. 2D allocates
  its per-thread bank in `prep_strategy2D_inpl_bank`; `s3D` has no bank. The search's per-thread work counts
  already flow through `strategy3D_cont%add_work(cont_work(ithr))` in `choose_and_run_strategy`.
- **Exits and memory.** The weighted insertion's margin on the resolution exits is `l_margin =
  l_prob_align_mode .and. topk > 1` (`simple_convergence`). The matcher's memory estimate counts
  `64 * max(1, topk)` bytes per particle for the insertion (`simple_matcher_memory`).
- **History that bounds the design** (refine3D policy, "The weighted insertion"). Before the cone, extras over
  the whole row at the physical temperature made the map the posterior expectation from a stage's first
  iteration: every stage-3 map was bad on PorF and PfCRT, while the late stages gained. The quadrature is also a
  physical-temperature expectation over a whole row, but over one coordinate and conditional on the drawn
  direction. In 2D that was enough. In 3D it is untested, and stage 3 is the first thing Hans's runs read.

## 3. Contract

- **Off state.** With `inpl_posterior=no`, the default, the run is the base's byte for byte at one thread. No
  in-plane row is evaluated, no bank or reservation is allocated, the layout is `setup_range`'s at the base,
  gridding and PCG run the base's insertion code, `l_margin` is as at the base, and the readout's new fields are
  zero. This is the gate no ruling may waive.
- **The committed pose is the hypothesis.** For identical incoming references, the table, its draws and keyed
  streams, the extras' selection, every continuous solve, the project's pose, `corr`, shift, sigma2 residual
  and convergence statistics are unchanged at either key value. Under the key, the references change from the
  second iteration on, so "unchanged" holds for one pass, not for a run.
- **Concentrated limit.** A particle at `topk = 1` whose committed row collapses enters at its solved pose at
  weight 1, as with the key off. With extras, a particle whose rows all collapse enters at the base's poses with
  the marginal weights of decision 4. It differs from the key-off run in its weights only.
- **Mass and determinism.** The weights of a particle sum to 1, so the reconstruction's mass and the trailing
  chain's population accounting are unchanged. The phase is keyed on the particle and the reference, and each
  particle writes its own reserved block. The entries therefore do not depend on the thread count, the batch or
  the part.
- **Modes equivalent.** Shared memory and distributed parts apply the same rule to the same particles. Nothing
  in it depends on the execution mode.
- **Final maps from committed poses.** The final and bootstrap reconstructions, `reconstruct3D` from the project
  and the multi-state final assignment use committed poses, as with `topk > 1` today.

## 4. Decisions and rulings

1. (2026-11-14, Hans) **The switch.** A temporary key `inpl_posterior=yes|no`, default `no`, on `refine3D` and
   `solve3D`. solve3D forwards it to every stage through its copied command line (section 2). After Hans's read
   it becomes the shared constant or is removed, as `topk_inpl` went in 2D (`b2a9d544f`).
   (2026-11-15, Hans) The key is removed: the insertion runs in every pass decision 3 names that reconstructs
   (`volrec=yes`) under `objfun=euclid`, so the multi-state final assignment needs no override.
2. (2026-11-14, Hans) **Testing.** The plan carries the off-state byte identity, unit assertions with planted
   violations, and smoke runs. The scale test is Hans's, after delivery. The report names the runs and the
   readout to read. The plan closes on delivery.
3. (2026-11-14, Hans) **Where it runs.** Every pass with a posterior table except the hill-climb
   neighbourhoods: `refine=prob`, and `refine=prob_neigh` with `prob_neigh_mode=state|geom`. In solve3D that is
   single-state stages 3 to 8 and the multi-state stages from 3. It does not run in `prob_neigh_mode=shc|snhc`
   (stages 1-2), in `refine=shc`, in the multi-state final assignment (which reconstructs nothing; the commander
   sets `inpl_posterior=no` there beside `topk=1`), or in the final and bootstrap reconstructions. It is
   independent of `topk` and of the cone: it runs at `topk=1` and before the first FSC exists. Derived once as
   `l_inpl_posterior` in the parameter phases, as `l_shift_prior` is.
4. (2026-11-14, Hans; the shift-prior term added after the review, section 10) **View weights: the in-plane
   marginal.** Under the key, the weight of each inserted pose `k` (the committed pose and its extras) is
   recomputed from the rows at the solved shifts, keeping the shift prior the insertion applies today (section 2,
   the shift-prior bullet):

   ```text
   log m_k = log Z_k - T L_min,k - P_k        (then minus the largest log m_k, exponentiated, normalised)
   P_k = prior_w * |du_k - c_k|^2
   ```

   Under the polar shift map, `prior_w = map_prior_weight(shift_sigma, smpd)`, `du_k` is the pose's solved
   increment from the stored shift and `c_k = 0`, which is `T` times `insert_pose_loss`'s prior term. Under the
   window search with `shift_prior=yes`, `prior_w = 1 / (2 sigma2_prior_model)` and `c_k` is that pose's actual
   stored-frame prior centre, after rotating and clipping the seed as described in section 2 (zero without the
   single-state seed). The table retains that centre with each refined entry and exposes it, with the width,
   to the strategy for the committed pose and extras. Capture it at the refinement's start: its final in-plane
   index may differ, so neither the unchanged particle seed nor clipping at the final angle reproduces it.
   The subsequent continuous solve leaves this retained centre unchanged for the marginal weight; the search,
   selection and solves themselves are unchanged. Otherwise `P_k = 0`. The conditions are the existing insertion path's, and
   one helper computes `P_k` for both routes. The product with the quadrature, `m_k / K` per point, is then the
   joint posterior over view and angle at the solved shifts under the same prior, restricted to the kept poses.
   With the prior retained, the table's weight differs from `m_k` by the factor `Z_k` and the solve's small moves
   only (section 2, first bullet), so it favours poses whose in-plane row is sharp. The selection of the extras
   stays the table's (`top_posterior_extras`); only the kept set is reweighted. After reweighting, an extra under
   `PROB_INSERT_WMIN = 0.01` is dropped and the rest renormalised; the committed pose is always kept. With the key
   on, this replaces the 2026-11-01 rule that "the weights stay the table's". With the key off, that rule stands.
   What it does to the amount of averaging: it moves weight toward poses whose angle is uncertain, and every
   inserted entry carries its own share of the joint posterior. Usually that lowers the largest single entry,
   but not always. Example: a pose with two equally plausible angles fits best, and a one-angle pose has half its
   likelihood (a loss higher by `ln 2 / T`). Marginal gives 0.4, 0.4, 0.2 (`n_eff` 2.87); table weights give
   1/3 each (`n_eff` 3.0). With equal likelihoods the order reverses: marginal 1/3 each, table 1/4, 1/4, 1/2
   (`n_eff` 2.83). The readout compares the committed view's marginal share `m_c` with the table's committed
   weight (decision 10).
5. (proposed; the 2D plan's rulings on the row, the temperature, the selection, the phase and the continuous
   angle, per pose) **The product rule.** For each inserted pose, after its continuous solve:
   - the whole in-plane row at the pose's solved shift, at the physical temperature
     `temperature_or_argmin(posterior_temperature(wsqsum_ptcl_now(iptcl), pftsz))`;
   - the `K = INPL_POSTERIOR_BUDGET = 64` point quadrature, with the phase
     `v = rnd_keyed(rnd_draw_key(prob_rnd_scope(which_iter, .false.), iptcl, RND_FAM_INPLQ, reference=iref))`,
     so that two poses of one particle take independent phases;
   - the point or points in the solved bin placed at the solved continuous angle (`resolve_inplane_e3`) and merged
     into one entry;
   - a row collapsed to one entry inserts the pose at its solved angle at its weight, whichever bin that is;
   - the committed pose not inserted beside its own quadrature (it is one sample of the same posterior), so the
     committed angle's weight is `m_c` times the solved bin's point count over `K`, which can be 0;
   - every entry of pose `k` at weight `m_k` times its point count over `K`, inserted at the pose's direction,
     shift and the entry's angle.
6. (proposed) **Where the rows are evaluated.** In `strategy3D_prob%srch`, on the particle's thread:
   - each extra's row in `resolve_extra_poses` after `refine_extra_pose_continuously`;
   - the committed row after `refine_assignment_continuously`;
   - the weights (decision 4) and the fill of the particle's reserved block (decision 8) after all rows.
   A row is `memoize_shifted_ptcl(iptcl, sh, s3D%inpl_bank(:,:,ithr))` then
   `gen_euclid_dists_bank(iref, iptcl, bank, dists)`, with `sh` the stored-frame increment the solve returns
   (section 2). It falls back to `gen_euclid_dists(iref, iptcl, sh, dists)` when `|sh|^2 <= SHERRSQ`, as in 2D.
   Without continuous refinement (no shifts, or an invalid solve), it is taken at the table's increment. If any
   pose of a particle has no valid row, the whole particle enters as with the key off (table weights, one angle
   per pose), so marginal and table weights are never mixed within a particle. A particle needs only its per-pose
   results (entries, `L_min`, `ln Z`, `n_eff`), not the rows themselves. `inpl_quadrature` gains optional `lnz`
   and `dmin` outputs; the 2D callers do not pass them. The cone is unchanged.
7. Withdrawn 2026-11-14 after the review (section 10): no point-group merge, because the search already removes
   the symmetric duplicates. The directions are a spiral over the asymmetric unit (`build_refspiral`), so a
   symmetry mate of a pose never appears as another direction. Each direction's in-plane row covers the full
   360 degrees, though. The only poses that row duplicates are those of a direction lying exactly on a symmetry
   axis: the top view along a `Cn` axis, where in-plane angles `360/n` apart are the same pose, and the views
   along the 2-fold axes of `Dn`. Every other direction's row peaks are different poses. Merging at one in-plane
   step (the first version) would also have merged distinct poses near an axis, which the review showed the rim
   resolves (a `C4` view 0.4 degrees off the axis, `dang = 1` degree: 0.57 degrees, about one pixel at a 100-pixel
   rim). Merging only the exact on-axis duplicates would change no map, because every insertion writes all
   point-group copies anyway, and would save only extra insertions for those few directions. A pose has at most
   `K` entries either way, so the worst case is the same. Every entry is inserted as it is, and the REC3D bench's
   poses inserted shows the cost.
8. (proposed; revised after the review) **Layout: a worst-case reservation.** Under the key, `prob_insert_poses`
   reserves, at the first batch and for the whole part, one block per pose of every particle: the committed pose
   and up to `topk - 1` extras, so `max(1, topk)` blocks per particle. A block holds its entry count, the solved
   bin, and up to `K` packed 32-bit entries (the rotation index in the low 16 bits, the point count in the high
   bits; `nrots < 2^16` and the count is at most `K`). That is `max(1, topk) * (8 + 4 K)` bytes per particle,
   792 bytes at `topk = 3` and `K = 64`, 79.2 MB per 100,000 particles of a part, on top of today's extras slots.
   The reservation is never exceeded: a quadrature gives at most `K` entries. It does not grow, copy or compact.
   The pose weights live where they live today: an extra's `w` holds its `m_k`, and the committed pose's `wc`
   holds `m_c`, the committed view's weight (decision 10). An entry's weight is its pose's weight times its count
   over `K`. The solved bin's entry is inserted at the pose's continuous angle (the project's `e3` for the
   committed pose, the extra's `euler(3)`), every other entry at its grid angle `360 - rot(irot)`. A block with
   one entry is inserted at the pose's continuous angle whatever its bin (decision 5). A block with
   count 0 (a particle that fell back, decision 6) leaves the particle to the key-off insertion. The table's
   extra count stays apart from the blocks. With the key off nothing is reserved and `setup_range` is the base's.
9. (proposed; revised after the review) **Insertion paths and bounds.** With the key off, gridding and PCG run
   the base's code (`insert_extra_poses`, `pcg_extra_entries`, `accumulate_pcg_extras`) unchanged. Under the key:
   - gridding: one new loop per particle, over its poses, then over each pose's entries. The committed plane is
     re-shifted once per pose by the pose's shift minus the committed one (none for the committed pose), and each
     entry is one `insert_plane_oversamp` at its angle and weight;
   - PCG: the entries (the committed pose's and the extras') are expanded in chunks of at most `MAXIMGBATCHSZ`
     entries. Each chunk gets its own `oris` records and weights, its `prep_particles`, its observation copies and
     one `accumulate_batch` from index 1. The per-entry metadata, the records and the plane scratch are bounded by
     one chunk whatever the part's size, `topk` or `K`, which is less than the base holds for a half of more than
     `MAXIMGBATCHSZ` particles. The B and D accumulators span the half as today and come out byte for byte as one
     call over the same entries would make them (section 2). The contributor count and mass published with them
     stay the particles'.
   Section 2 records the check made before development: nothing on the production path reads per-entry metadata
   after accumulation. The REC3D bench file gains the poses inserted next to the particle images inserted.
10. (proposed; revised after the review) **Readout.** The 3D route fills the existing `inpl_*` fields
    (`inpl_kmax = K`, `inpl_mmax`, `inpl_nptcls` particles with a row, `inpl_npose` entries inserted,
    `inpl_nsingle` rows collapsed, `inpl_neff_sum` over rows) and two new ones:
    - `inpl_nrows`, the rows evaluated;
    - `inpl_committed_view_weight_sum`: the weight the committed view (its direction, all angles together)
      receives after the in-plane-marginal reweighting of decision 4, `m_c`, summed over the table's particles,
      with 1 for a particle without extras. It is the counterpart of `ins_wc_sum`, the table's committed weight
      that `PROB_TAB_INSERT` prints, so the two lines show how far the reweighting moves the particle's mass
      between its committed direction and its extras. It is a view weight, not the
      committed angle's inserted weight (`m_c` times the solved bin's count over `K`, which at `topk = 1` with a
      flat row is 1/64 or 0 while `m_c` is 1). The angular spreading shows in the collapsed share and `n_eff`, so
      the committed angle's weight is not reported.
    The new fields are written, read and summed with the others. `PROB_TAB_INPL (budget K, m_max M)` follows
    `PROB_TAB_INSERT` in the 3D branch. It prints rows and entries per particle, the collapsed share of rows, the
    mean `n_eff`, and the committed view weight beside `PROB_TAB_INSERT`'s table value. The 2D line is unchanged.
    The counts are set on the table before `close_assignment` (`set_inpl_counts`, on `eul_prob_tab` and copied
    into both tables' records). Revised by the run in Phase 2 (launch ruling R2): each particle's search puts its
    counts into the particle's own slot of the pass (`add_inpl_counts`), and the matcher sums the slots in particle
    order. Sums per thread, the `add_work` pattern first proposed, depend on which particles the dynamic schedule
    gives each thread, so the double-precision sums of a part's record could differ between repeats at more than one
    thread (found by the Phase 2 review). The rows, entries, collapsed rows and `n_eff` are counted over the
    particles that enter by their rows. A particle that falls back to the insertion without the key (decision 6)
    adds its table weight to the committed view weight sum, and a particle without an assignment adds 1.
11. (proposed) **Exits.** `l_margin = l_prob_align_mode .and. (topk > 1 .or. l_inpl_posterior)`. The quadrature
    averages the draw's noise out of each half, as the extras do, so the halves agree a shell early for the same
    reason.
12. (proposed; revised after the review) **Memory estimate.** Under the key, the matcher's estimate adds the
    reservation of decision 8 exactly, `max(1, topk) * (8 + 4 K)` bytes per particle of the part. It is a worst
    case, not an expectation, so the first iteration cannot exceed it. The alignment context's share of the
    budget then grows by that amount, and the pass buffer takes the remainder, as it does today. Over budget,
    the existing warning is printed. The PCG chunk of decision 9 is bounded by `MAXIMGBATCHSZ` entries, the order
    of one observation batch the reconstruction already holds, so it adds no estimate term.
13. (proposed; revised by the run in Phase 1 under launch ruling R2) **Objective.** The key requires
    `objfun=euclid`, as the probability tables already do (`gen_posterior_val`). The parameter phases refuse
    `inpl_posterior=yes` with another objective on a pass the key acts on: a `refine3D` pass in one of decision 3's
    modes. That is also where `l_inpl_posterior` is derived; elsewhere it is false and the key is ignored. Revised
    because the children of solve3D and refine3D copy the parent's command line, key included, and switch it to
    `objfun=cc`: the start-volume reconstruction of `simple_refine3D_strategy.f90` keeps `refine=prob`, solve3D's
    reconstructions (`simple_solve3D_utils.f90`) and the external-reference pose initialisation (`refine=greedy`)
    do the same. Refusing the key whatever the program and mode would stop every key-on run at such a child.
    solve3D itself sets `objfun=euclid` before its parameters are made.
14. (2026-11-14, Hans) **PCG is supported and exercised.** The PCG reconstruction takes the quadrature's
    entries as it takes the extras (decision 9), and the plan runs it. solve3D reconstructs by gridding unless
    `rec_backend=pcg` is given, and then from `PCG_START_STAGE = 6` (`set_refine3D_backend_policy`), so the
    gridding smoke runs never reach `pcg_extra_entries` and `accumulate_pcg_extras`. The smoke therefore includes
    a solve3D run with `rec_backend=pcg`, in shared memory and with `nparts=2`, key on. The same command with the
    key off is compared byte for byte against the base at one thread, which proves that the PCG path is unchanged
    in the off state.
15. (2026-11-14, Hans) **Multi-state is implemented here and tested separately afterwards.** The enabled path
    covers every state count. The rows take the multi-state reference index
    (`(istate - 1) * nspace + iproj`), the extras stay in the committed state's block, and the marginal weights run
    over the committed state's poses only, so the particle never enters another state's map. The multi-state
    stages from 3 run it (decision 3). This plan's smoke runs are single-state. Multi-state is covered here by
    the unit assertions and the off-state identity only; its scale test is a separate step after delivery.

Defaults taken without a question: the key name `inpl_posterior`; `l_inpl_posterior`; the field names
`inpl_nrows`, `inpl_committed_view_weight_sum`; `add_inpl_counts`, `set_inpl_counts`, `prep_strategy3D_inpl_bank`; reuse of
`RND_FAM_INPLQ` (its comment widened to both routes); `K` shared with 2D.

## 5. Not in this plan

- The table, its draws, the extras' selection and the cone. A table row carrying the in-plane marginal would
  change the search, which is an experiment of its own (searches are experiments, measured on restarts).
- Shift marginalisation and shift extras; class or state extras; soft state membership.
- The Cartesian route (`refine3D_local refine=cont`): it has no polar row and no `topk`.
- solve3D_chunks, the solve3D add-on and the stream's 3D stages: the key reaches them only where their command
  lines carry it. Their command builders are not edited, and their behaviour under the key is untested.
- 2D, the class averager, `make_cavgs` and the 2D readout line.
- A `K` tied to `m_max`. When a stage's `m_max` reaches `K`, the readout shows it; the 2D grid residual (the
  points snap to the rotation grid, so cancellation is exact only where `K` divides `nrots`) applies unchanged.
- The exact angular marginalisation by a polar circular convolution per ring (2D plan, section 5).
- A default of `yes`, and the removal or conversion of the key (decision 1, after Hans's read).
- Larger-scale testing (decision 2) and the simulated workflow's command lines: its `POSE RECOVERY` lines are the
  known-pose readout, but carrying the key into the workflow is not in this plan.
- Multi-state smoke and scale testing (decision 15): a separate step after delivery.
- Skill edits (`simple-refine3d`, `simple-main-strategies-probabilistic`): one line each is proposed in the
  report and applied only with approval.

## 6. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 baseline | none | Release and Debug builds, fast gate, both simulated workflows twice at one thread, the PCG smoke twice at one thread, smoke dry runs, thread-count probe |
| 1 plumbing | identical (key off) | `test=identity` at one thread against Phase 0; the PCG smoke at key off byte-identical to the base; fast gate; unit suites; static gates |
| 2 enabled path | identical at key off, result-changing at key on | identity at key off again (workflows and the PCG smoke); smoke at key on in Release, Debug and `nparts=2`, gridding and PCG; thread gate per Phase 0 |
| 3 assertions | tests-only | each planted violation fails its check, then passes; `scripts/check_test_registry.py`; fast gate |
| 4 handoff | none | the report and the runs named for Hans |

## 7. Acceptance

- **Identity (the gate).** `simple_test_exec test=identity` over `test=simulated_workflow` (`suite=1jxy` and
  `suite=6vxx`, each running solve2D and solve3D) at one thread with `SIMPLE_SEED` fixed, base against the build
  with the key off, compared in the directories the runs executed in (copies fail, 2D report finding 6).
  Exclusions are declared from Phase 0's two base runs and never added later (*measured*; the 2D plan's two base
  runs needed none). An item that cannot be made identical is restored and recorded, never rerun until it passes.
  This gate proves the off state only.
- **Thread gate (Phase 2).** One `refine3D refine=prob` iteration on the smoke project at `nthr=1` and `nthr=16`,
  three repeats each, key on. If Phase 0 shows the base's partial reconstructions byte-identical across these
  thread counts, the key-on partials must be too (*measured*). If not, the base's spread is recorded and the
  key-on difference must not exceed it.
- **Assertions of the enabled path (Phase 3)**, each planted with its violation first, in existing suites:
  - `refine3D in-plane state` (`simple_strategy3D_inplane_tester`, `unit_pftc_align2D3D`). On the tester's
    fixture (it runs `ctf=no`) with a non-zero stored shift and a known angle: after
    `refine_assignment_continuously` and `refine_extra_pose_continuously`, the row at the returned increment has
    its minimum within one bin of the solved angle. Its value at the solved bin equals the shifted-reference
    kernel at `shift_ptcl2ref(sh, rot)` within the single-precision spread measured on that pair before the plant
    (*measured*, recorded in Progress). The bank kernel's own frame under an angle-dependent CTF is the 2D check
    `check_bank_particle_frame` (mirror tester), which stays as it is. Plants: the increment passed unrotated;
    the sign reversed.
  - `3D probability table` (`simple_eul_prob_tab_tester`). The marginal weights on synthetic rows (a one-angle,
    a flat and a two-peak row) equal `exp(log Z - T L_min - P) / sum` within `1e-6` relative (*derived*:
    double-precision sums of at most `3 nrots` terms). The committed pose is kept below `PROB_INSERT_WMIN`, an
    extra under it is dropped and the rest renormalised, the entry weights of a particle sum to 1, and `topk=1`
    gives the committed pose weight 1. The worked example of decision 4 is reproduced (0.4, 0.4, 0.2). Plants:
    `Z` omitted; the `PROB_INSERT_WMIN` rule not re-applied; the committed pose dropped.
  - Same suite, the shift prior (decision 4). Two poses with identical rows, one at the prior's centre and one two
    prior sigmas from it, take 0.881 and 0.119 (`1 : exp(-2)`), under the map (centre: the stored shift, width
    `shift_sigma`) and under the window with `shift_prior=yes` (centre: each pose's retained prior centre,
    variance `sigma2_prior_model`). With neither prior active they take 0.5 each. Also exercise a seed
    `(0.9*trs, 0.9*trs)` at starting angles 0 and 45 degrees: the first centre is unchanged; the second is
    `(trs/sqrt(2), trs/sqrt(2))` after clipping. Check the captured centres and marginal weights against this
    closed form, including a refinement whose final in-plane index differs from its starting index. The rows'
    synthetic examples above cannot detect a missing or miscentred prior; these cases can. Plants: the prior
    omitted; the window prior centred at the stored shift; the unchanged seed reused for all views; the centre
    reconstructed using the final angle.
  - Same suite, the layout and its reservation (decision 8). A batch whose every pose has a flat row (a broad C1
    posterior) fills every block with `K` entries, which fits the reservation exactly, and the reserved bytes equal
    the memory estimate's term. A batch with fills from 0 to `K` per block decodes to the same rotation indices and
    counts as written, and the table's extra counts are kept. With the key off the layout equals the base
    `setup_range`'s on the same table, counts and first slots exactly, and nothing is reserved. Plants: a count
    packed into the low bits; the reservation sized by `topk - 1` blocks.
  - Same suite, the committed view weight against the committed angle weight (decision 10). At `topk = 1` with a
    flat row, the readout counts `m_c = 1` while the committed angle's inserted weight is its solved bin's count
    over `K`. Plant: the field fed with the angle's weight.
  - `weighted insertion` (`simple_reconstructor_tester`, `unit_reconstruction`). The key-on loop's volume for a
    particle with entries at two poses equals inserting each entry separately, with `fplane_shift` computed per
    entry from the committed plane, bit for bit in the same order. Entries at one pose with weights summing to 1
    equal one insertion at weight 1 within `1e-6` relative (*guessed*, frozen here). Plant: the first pose's
    shifted plane carried over to the second pose.
  - `rec3D backend` (`simple_rec3D_strategy_tester`, beside `test_pcg_extra_entries`). The key-on PCG path with
    chunks smaller than a batch's entries expands the same entries, in the same order and with the same weights,
    as one flat list, and accumulates each exactly once. The B and D accumulators equal one flat accumulation of
    that list byte for byte, at one and at several threads (each cell takes its contributions in entry order,
    section 2). The published contributor count is the particles'.
    Plant: an off-by-one at the chunk boundary.
  - `simple_convergence_tester`: `topk=1` with `l_inpl_posterior` takes the weighted insertion's margin, and the
    key off keeps the single draw's rules. Plant: the margin ignores the key.
  - Same `3D probability table` suite: the readout record's round trip (write, read, add) carries the new fields.
    Plant: a field missing from the sum.
- **Repository gates, unchanged.** Fast gate (`scripts/run_fast_gate.sh`), `test=unit_pftc_align2D3D`,
  `test=unit_reconstruction`, the static gates (`check_cli_chain.py`, `check_test_registry.py`, `check_descr.py`,
  `check_guide_imports.py`, `check_tester_linkage.py`), the complete UI validation, a clean Debug build.
- **Smoke, not a test of the effect.** solve3D on a copy of the 6vxx workflow's project with `inpl_posterior=yes`:
  - gridding: shared memory in Release, shared memory in Debug (bounds checks on the new arrays), and `nparts=2`;
  - PCG (decision 14): `rec_backend=pcg` in shared memory and with `nparts=2`, which reaches the PCG worker from
    stage 6;
  - each completes;
  - `PROB_TAB_INPL` appears from stage 3 on and not in stages 1-2;
  - the parts' counts are summed in the master's line.
  Phase 0 fixes the exact command lines and their wall times.
- **PCG off-state identity.** The PCG smoke command with the key off, at one thread with `SIMPLE_SEED` fixed,
  gives the same bytes as the base on every image and binary output (*measured*: Phase 0 runs it twice on the
  base and declares its exclusions then, as for the workflows).
- **Reported, not gated.** On the smoke run at key on and off: the REC3D accumulation time and poses inserted;
  the `PROB_TAB_INPL` lines (rows and entries per particle, collapsed share, `n_eff`, the committed view weight
  against `PROB_TAB_INSERT`'s); the stage-3 and final FSC lines; the peak resident memory of a part at key on and
  off beside the estimate's reservation. Paths go in Progress. The 6vxx fixture cannot say whether the maps
  improve.
- **Validation of the effect** is Hans's (decision 2). The report names the runs, all with the solve3D test
  protocol (`nsample=2500`, a common `solve3D_cavgs` start, `cavg_ini_ext=yes`):
  - PorF (`~/Processing/PorF`), ten restarts per arm, key on and off, paired by seed;
  - PfCRT, ten restarts per arm;
  - streptavidin (`D2`: symmetric views keep all their peaks' entries, so read the REC3D poses inserted);
  - one small C1 membrane protein (SLC or Tat), where the top views' in-plane angle is the question.
  The readout to read, in order: stage-3 maps and the stage-3 exit iteration; final FSC 0.143; `PROB_TAB_INPL`
  per iteration; REC3D accumulation. A suggested reading rule for Hans to adopt or replace before the runs: the
  arm is not worse if it matches the control's restart count at the control's final resolution range, and its
  stage-3 exits are no later in more than two of ten paired seeds.

## 8. Common rules and stop conditions

- No parameter but `inpl_posterior` (decision 1). It is registered in `simple_parameters.f90`, parsed, checked
  and derived in the phases, and present in the UI of refine3D and solve3D. `simple_args.f90` and the UI default
  values are generated.
- No floor, tolerance or assertion relaxed unless this plan or a ruling says so. Fix as you go, with a file-table
  row before the edit. Every file changed outside the table is listed in the report. Dead code found on the way
  is deleted.
- Thread policy: per-thread work is sized by the team in use (`omp_get_max_threads()` read before allocation).
  There is no allocation per particle in the search loop. Byte identity is judged at one thread; the thread gate
  of section 7 covers the parallel region.
- Compile-time policy: the new components are allocatable; imports are lean `use ..., only:`; no new CMake targets
  or per-file options.
- Builds and runs happen on the dell under arun (Hans, 2026-11-14, launch ruling R1; the machine notes and
  build trees are the run's `spec/phases.md`). Never into the checkout's `build/`. Long runs go in the
  background. Nothing is committed. Master is not frozen during the run; the run does not sync with it.
- Sync: no other run is known to edit these files on 2026-11-14. If master moves past `06d3a4de9` in a file of
  the table before launch, re-check section 2 against the new commit and record the differences in Progress.
- **Stops.** The off-state identity fails and cannot be restored. The frame assertion fails (the row disagrees
  with the kernel at the solved pose): stop, do not adjust the tolerance. A file outside the table needs a
  change that alters behaviour outside the plan's scope (launch ruling R3 names the scope; a row may be added
  for a change the design implies within it). A maintainer's decision (1 to 4, 14, 15) contradicts the code;
  the proposed decisions may be revised by the run within the contract (launch ruling R2). Decision 4 turns out to need the table's selection changed, or
  the table cannot retain and expose each refined pose's actual prior centre and the prior width to the strategy. Section 2's PCG premise
  (nothing on the production path reads per-entry metadata after accumulation) no longer holds at the launch
  commit.
- **Gates no ruling may waive.** The off-state identity; the fast gate and the unit suites; the Debug smoke at
  key on.

## 9. Phases

- **Phase 0, baseline.** Release and Debug builds of `06d3a4de9` in the launch environment (section 8); the fast
  gate; machine, compiler and load recorded; the base build kept. `test=simulated_workflow`, both suites, twice at
  one thread with a fixed seed, for the identity exclusions. The gridding and PCG smoke commands run by hand on the
  base without the key, to check the UI's required keys and record their wall times; the PCG command twice more at
  one thread with `SIMPLE_SEED` fixed, for its identity exclusions. If a smoke run exceeds one hour, a smaller
  fixture is chosen and named. The thread-count probe of section 7 on the base. The base's `PROB_TAB_INSERT` lines
  and REC3D bench recorded. Exit: Progress holds the build records, the run directories, both exclusion lists
  (possibly empty), the smoke command lines with their times, and the thread probe's result.
- **Phase 1, plumbing (identical with the key off).**
  - The key, its parse, check, derivation and UI.
  - The final multi-state assignment's `inpl_posterior=no`.
  - The reservation of decision 8 (allocated only under the key), retention of each refined pose's actual
    stored-frame prior centre, and the table's getters for that centre and the iteration's prior width.
  - The readout fields, the 3D line and `set_inpl_counts`.
  - `l_margin`, the memory estimate, the REC3D bench's poses inserted, and the bank allocation routine (used only
    under the key).
  Exit: `test=identity` at one thread against Phase 0 gives 0 differing files outside the declared exclusions;
  the PCG smoke at key off gives the base's bytes (section 7); the fast gate, `unit_pftc_align2D3D`,
  `unit_reconstruction` and the static gates pass; the UI validation passes; `inpl_posterior=maybe` and
  `inpl_posterior=yes objfun=cc` are refused with their messages.
- **Phase 2, the enabled path.**
  - In `strategy3D_prob`: the rows after each solve (decision 6), the quadrature and its phase (decision 5),
    the marginal weights with the shift prior (decision 4), the fill of the reserved blocks, and the per-thread
    counts.
  - `inpl_quadrature`'s optional `lnz` and `dmin`.
  - The key-on insertion paths of decision 9: the gridding loop and the chunked PCG expansion.
  - The policy text, in the same phase: a subsection "In-plane posterior insertion" in
    `doc/policies/3D/refine3D_policy.md` after the weighted insertion, its sentence on the 2D route updated; the
    key's reach and the final assignment in `doc/policies/3D/solve3D_policy.md`; one sentence in the 2D policy's
    in-plane section pointing to the 3D counterpart.
  Exit: identity at key off against Phase 0 again (0 differing), for the workflows and the PCG smoke; the smoke
  runs of section 7 complete, gridding and PCG, with `PROB_TAB_INPL` from stage 3 and none in stages 1-2; the
  thread gate passes; the fast gate and both unit entries pass; the reported numbers and paths are in Progress.
- **Phase 3, the assertions.** The checks of section 7, each written with its violation planted first and shown
  to fail, then the plant removed. Exit: every planted violation fails its check; the suites pass;
  `scripts/check_test_registry.py` and the fast gate pass; the touched testers equal their pre-plant content
  apart from the new checks.
- **Phase 4, handoff.** The report `inpl_quadrature_insertion_3D_report.md` beside this note, in the order of
  `doc/how2s/how2write_agent_plans.md` section 8: what was learned, files changed with counts, decisions taken
  under delegation, deviations, the identity and thread evidence, Hans's runs with their command lines and the
  readout to read, the proposed skill lines, and a commit message. Exit: the report exists; the tree is left
  uncommitted with the changed files listed.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | `inpl_posterior`, its parse, check, refusal with `objfun` other than euclid, `l_inpl_posterior` (decision 3) | 1 |
| `src/main/ui/simple/simple_ui_refine3D.f90`, `simple_ui_solve3D.f90` | the key | 1 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | `inpl_posterior=no` in `run_multistate_final_assignment` | 1 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | `prob_insert_poses` reservation, blocks and table count kept apart; refined-entry prior-centre retention; readout fields and the 3D line; `inpl_quadrature` outputs; the marginal-weight helper with the shift-prior term (Phase 3 under launch ruling R3: the seeded window refinement and the particle counts as testable helpers; plants, reverted) | 1, 2, 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90` | `set_inpl_counts`; the counts into each table's record at `close_assignment`; capture each refined pose's actual prior centre and expose it with the prior width (Phase 3 under launch ruling R3: both tables call the seeded window refinement helper) | 1, 3 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | layout gated on `topk > 1 .or. l_inpl_posterior`; per-thread counts; `m_max`; the bank's allocation | 1, 2 |
| `src/main/strategies/search/simple_strategy3D_alloc.f90` | `s3D%inpl_bank`, `prep_strategy3D_inpl_bank`, release in `clean_strategy3D`; the per-thread per-pose work of the rows (Phase 2, row extended under launch ruling R3) | 1, 2 |
| `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90` | the rows, quadratures, marginal weights, block fill, `add_inpl_counts` (Phase 3 under launch ruling R3: the particle counts through the helper) | 2, 3 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | if needed: access to the solved increment for the row (Phase 3 under launch ruling R3: plants only, reverted) | 2, 3 |
| `src/main/strategies/search/simple_matcher_3Drec.f90` | the key-on gridding loop; poses inserted in the REC3D bench (Phase 3 under launch ruling R3: `insert_inpl_entries` public, on the point group; plants, reverted) | 1, 2, 3 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | the key-on chunked PCG expansion; the poses inserted in its REC3D bench (Phase 1, row extended under launch ruling R3) (Phase 3 under launch ruling R3: the chunk loop as a routine over planes in memory; plants, reverted) | 1, 2, 3 |
| `src/main/simple_convergence.f90` | `l_margin` (decision 11) (Phase 3 under launch ruling R3: plant only, reverted) | 1, 3 |
| `src/main/strategies/search/simple_matcher_memory.f90` | the estimate (decision 12) | 1 |
| `src/utils/math/simple_rnd.f90` | the `RND_FAM_INPLQ` comment covers both routes | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | if needed: the master's readout merge of the new fields | 1 |
| `src/main/strategies/search/simple_strategy3D_inplane_tester.f90` | the frame assertion | 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | marginal weights, shift prior, layout and reservation, view weight against angle weight, readout round trip | 3 |
| `src/main/volume/simple_reconstructor_tester.f90` | the key-on gridding loop and weighted entries | 3 |
| `src/main/strategies/parallelization/simple_rec3D_strategy_tester.f90` | the chunked PCG expansion; the key-on gridding loop and weighted entries (moved here from the reconstructor tester: a volume submodule importing the matcher would invert the layering) | 3 |
| `src/main/simple_convergence_tester.f90` | the margin under the key | 3 |
| `src/main/strategies/search/simple_matcher_memory_tester.f90` | if the estimate's checks need the key's term | 3 |
| `doc/policies/3D/refine3D_policy.md`, `doc/policies/3D/solve3D_policy.md`, `doc/policies/2D/solve2D_policy.md` | the living description (Phase 2 bullet) | 2 |
| `.github/skills/simple-refine3d/SKILL.md`, `.github/skills/simple-main-strategies-probabilistic/SKILL.md` | one line each, with approval only | 4 (if approved) |
| `doc/implementation_notes/planned/inpl_quadrature_insertion_3D_report.md` | the report | 4 |

## 10. Review (2026-11-14) and what it changed

Two design reviews checked this proposal against the code at `06d3a4de9`. The first made four findings;
the second confirmed the symmetry, memory and readout revisions and identified a remaining detail in the
window-prior centre. All findings are incorporated below and in the contract, phases and acceptance checks.
This section is the review record; the separate review document has been removed. No build or runtime
validation was performed for either review; the implementation gates remain pending.

1. **The marginal weights dropped the shift prior.** The extras are weighed today with a prior on the shift:
   `insert_pose_loss` under the polar shift map, and the prior folded into the refined loss under the window
   search with `shift_prior=yes` (section 2). An in-plane row at one shift has the prior constant across
   rotations, so it cancels in `Z` but not between poses at different shifts. Without it, equal rows at increments
   0 and two prior sigmas would weigh 0.5 each instead of 0.881 and 0.119. Changed: decision 4's formula carries
   `P_k` under the existing path's conditions, centres and widths, and section 7 asserts the equal-row,
   unequal-shift case under both routes. The second review found that the window prior is centred on the
   rotated, component-clipped seed converted back to the stored frame, not always the original particle seed.
   A seed `(0.9*trs, 0.9*trs)` is unchanged at zero degrees but becomes approximately
   `(0.707*trs, 0.707*trs)` at a 45-degree starting angle after clipping. Changed: decision 4 retains `c_k`
   per refined pose at its refinement's start, and section 7 checks the clipped case and a changed final angle.
   Evidence: `refine_best_refs` in `simple_eul_prob_tab.f90` (lines 887–893 at the base),
   `refine_best_neighbors` in `simple_eul_prob_tab_neigh.f90` (1514–1519), and `set_prior_centre` and its call
   in `minimize` in `simple_pftc_shsrch_grad.f90` (283–299 and 572–588). The map prior and zero-seed case are
   unaffected; search and selection are unchanged.
2. **A one-step tolerance merges poses that are not equivalent.** Below one in-plane step does not mean
   equivalent under the point group: a near-axis view keeps a separation the rim resolves. Changed: decision 7 is
   withdrawn. There is no merge. Hans pointed out that the search covers only the asymmetric unit, so the only
   duplicates left are the in-plane angles of directions exactly on a symmetry axis. Merging those would change
   no map and would not lower the worst case.
3. **The memory budget was a guess where the feature makes the worst case likely.** A broad posterior fills all
   `K` points of every pose, so the estimate of 4 entries per pose could be 16 times low. Doubling and copying
   added peaks, and capping the PCG plane scratch left the per-entry PCG records and metadata unbounded.
   Changed: decision 8 is a worst-case reservation in a compact layout (`max(1, topk) * (8 + 4 K)` bytes per
   particle, with no growth or compaction). Decision 12 counts it exactly. Decision 9 expands the PCG entries in
   chunks of at most `MAXIMGBATCHSZ`, with metadata per chunk. At Hans's request, that the production path reads
   no per-entry metadata after accumulation was confirmed before development rather than left to Phase 1
   (section 2). The same reading showed that chunked accumulation is byte-identical to one call. The key-off paths
   are the base's code, untouched. Section 7 fills the reservation with a broad C1 posterior.
4. **The readout mixed the committed view's weight with the committed angle's.** `m_c` is a view weight; the
   committed angle's inserted weight is `m_c` times its bin's share and can be 0 at `topk = 1` with a flat row.
   Changed: decision 10's field is `inpl_committed_view_weight_sum`, the view weight compared with the table's
   committed weight (first named `inpl_mc_sum`; renamed because the abbreviation did not say what it holds). The
   angle's weight is not reported, and section 7 asserts the `topk = 1` broad-row case.

## Progress

(one entry per phase, written by the run)

### Phase 0, baseline (2026-11-14)

No source file changed; this entry is the only edit. All paths below are under `/data2/agent_runs/inpl3D/scratch/keep/`
(the run's kept scratch, `keep/` below); one line per run is in `keep/summary.tsv`, the session's journal in
`/data2/agent_runs/inpl3D/state/journal_phase_0.md`.

- **Base.** The run's HEAD is `85213734e`. It differs from this plan's base `06d3a4de9` in this plan file alone
  (`git diff --stat 06d3a4de9 85213734e`: one file, 607 lines added), so section 2 needed no re-check. Launch ruling R1
  (builds and runs on the dell under the run driver) is the one applied here; section 8 already carries it.
- **Machine and compiler.** The dell: Oracle Linux 8.10, 24 cores (8 performance cores, CPUs 0-7, and 16 efficiency
  cores), 62 GB, no GPU. GNU Fortran 15.2.1 (gcc-toolset-15), CMake 3.26.5, Python 3.12 for the build's checkers,
  `OPENBLAS_CORETYPE=Haswell`. The maintainer's checkout `~/src/SIMPLE/build` is configured the same way (Release,
  tests on, `USE_ARCHOPT`, OpenMP); it was read, never built or run. Load averages: 0.2 before the builds, up to about
  8 while the run's own jobs ran, otherwise idle; each timing below names its load in `keep/summary.tsv`.
- **Build trees.** Four trees, each configured by cmake directly with `-DBUILD_TESTS=ON -DNICE=OFF -DUSE_ARCHOPT=ON
  -DFORTLS=OFF -DPython3_EXECUTABLE=/usr/bin/python3.12`, built with make and installed into themselves: the base
  (`git archive` of HEAD in `keep/src_base`) in Release (`keep/build_base`, 131 s at 12 jobs) and Debug
  (`keep/build_base_debug`, 64 s), and the run's tree in Release (`repo/build`, 130 s) and Debug (`repo/build_debug`,
  66 s), identical to the base but for this note. The fast gate (`scripts/run_fast_gate.sh`, 15 suites) passed on all
  four: 15 of 15, 10 s in Release and 38 s in Debug (`keep/logs/fastgate_*.log`).
- **Environment pitfall found.** The login environment sets `SIMPLE_PATH` to the maintainer's `~/src/SIMPLE/build`, and
  the job scripts that the local queue writes (even for one part) start `$SIMPLE_PATH/bin/simple_private_exec`. The
  first launch of the workflows therefore ran the maintainer's workers; it was stopped within minutes, before any
  comparison, and its directories deleted. Every run of this plan now sources `keep/env.sh <tree>`, which points
  `SIMPLE_PATH` and `PATH` at the tree under test and removes the maintainer's build from `PATH`.
- **Identity baseline.** `simulated_workflow` (each suite runs a simulated dataset from movies through solve2D and
  solve3D at one thread) ran twice on the base Release build with `SIMPLE_SEED=7`, in `keep/identity/base1/<suite>` and
  `keep/identity/base2/<suite>`: 1jxy 822 and 845 s, 6vxx 663 and 631 s, peak memory 2.4 GB, all four ending normally.
  `test=identity` of `base1` against `base2`, in the directories the runs executed in: 1jxy 223 files the same, 6vxx
  186 the same, none differing (`keep/identity/cmp_base1_base2/<suite>/identity_report.txt`). **The exclusion lists
  of both suites are empty.** `base1` is kept as the reference of every later comparison; `base2`'s outputs were deleted
  after the comparison (its logs are kept).
- **Smoke fixture.** `keep/fixture/fixture_6vxx.simple`: the 6vxx workflow's solve2D project of `base1` (87 particles,
  point group c3, box 192, 1.3 A per pixel), copied with its stacks and class averages (the movies hard-linked) and
  relocated with `simple_exec prg=update_project old_root=<base1 tree> new_root=<fixture tree> projfile_out=...`.
  `print_project_info` opens it. `update_project` leaves six side-file keys of the micrograph records
  (`mc_starfile`, `mcmodel`, `thumb`, `mceps`, `ctfdoc`, `ctfjpg`) pointing at `base1`; a solve3D run of the fixture
  under `strace` opened, probed or stat'ed no file under `base1` (`keep/smoke/fixturecheck_strace_summary.txt`), so the
  fixture is self-contained for these runs. The thread probe uses a second relocated project,
  `keep/fixture/fixture_6vxx_refined.simple`, the 6vxx workflow's solve3D result.
- **Smoke commands** (base Release, no key, each in `keep/smoke/<label>/` on its own copy of the fixture, `SIMPLE_SEED=7`;
  solve3D's required keys are `pgrp`, `mskdiam` and `nthr` with a project, and the other keys are the workflow's):
  `simple_exec prg=solve3D projfile=<copy> mkdir=yes pgrp=c3 mskdiam=180 gold_std_stage=9` plus
  - gridding, shared memory: `nthr=16`, 92 s, peak 2.46 GB;
  - gridding, distributed: `nparts=2 nthr=8`, 129 s, master 1.73 GB, largest part 2.21 GB;
  - PCG (preconditioned conjugate-gradient reconstruction), shared memory: `rec_backend=pcg nthr=16`, 131 s, 2.57 GB;
  - PCG, distributed: `rec_backend=pcg nparts=2 nthr=8`, 176 s, master 2.34 GB, largest part 2.21 GB.
  All ended normally, well under the hour, so the fixture stands. PCG takes over at stage 6 as section 4 says (its
  lines start after "ENTERING REFINE3D STAGE 6"). FSC=0.143 resolution (the half-map agreement criterion) at the end of
  stage 3 is 10.85 A in all four, at the end of stage 8 4.89, 4.80, 4.38 and 4.38 A. Peak memory is the processes' own
  high-water mark (`/usr/bin/time -v` and a sampler of `VmHWM`; `rss_peaks.txt` per run).
- **PCG off-state identity, exclusions.** The PCG shared-memory command at `nthr=1` twice
  (`keep/smoke/base_pcg_t1_a`, `_b`; 626 and 537 s, run side by side). Every image and binary output (`.mrc`, `.mrcs`,
  `.jpg`, `.bin`; 58 files) is byte-identical by `cmp`. What differs, from this comparison alone: the project file in
  `projinfo solve3D_run_id` (made from the process number) and the job record's wall-clock `date` and `time`; the
  `.provenance` files and the manifest's artifact lines in the run directory's path; the manifest's `run_id` and
  `checksum`; the PCG `*_diagnostics.txt` files in their timing lines; `metadata.json` (process number, timestamps) and
  `BENCH/` (timings). `keep/pcg_identity.sh <tree> <runA> <runB>` freezes these exclusions: everything else byte for
  byte, the project segment by segment with the run's path normalised. Base pair: 158 files the same, 0 differing
  (`keep/smoke/pcg_identity_base_a_b.txt`); against the `nthr=16` run it reports 88 differing, so the check can fail.
  `base_pcg_t1_a` is kept as the reference for Phases 1 and 2.
- **Thread probe.** One `refine3D refine=prob objfun=euclid maxits=1 mkdir=no` iteration on the refined fixture
  (`keep/probe.sh`; the run directory holds the fixture's `fsc_state01.bin`, so the insertion cone and the extras are
  active, and the start volume is the workflow's stage-3 map), three times at `nthr=1` pinned with `taskset -c 0-7`
  (22.6 s each) and three times at `nthr=16` (4.6-4.8 s). The partial reconstructions are
  `recvol_state01_part1_{even,odd}.mrc` and `rho_recvol_state01_part1_{even,odd}.mrc` (the density, SIMPLE's own format
  with a 64-byte header). The three single-thread runs are byte-identical, and so are the three 16-thread runs, but
  **one thread and 16 threads differ**: in the even and odd sums 42 and 27 of 7.08 million values differ, by at most
  2.5e-8 and 6.0e-9 of the volume's largest value; in the even and odd densities 14 and 14 of 3.58 million values,
  by at most 2.1e-10 of the largest. Headers are equal, and the committed poses of all 87 particles (the `ptcl3D`
  records) are identical. The noise spectra that the iteration bootstraps (`sigma2_state.bin`, 89 bytes) also differ
  between the thread counts. The gridding insertion adds into each lattice cell in a fixed order whatever the thread
  count (`insert_plane_oversamp`), so the spread most likely enters through the whitening of the planes by those spectra.
  Measured by `keep/mrc_spread.py` (standard library, reads the outputs). Few extras are inserted on this fixture
  (`PROB_TAB_INSERT`: 1.00 poses per particle); under the key the committed pose's own quadrature still exercises the
  new insertion path. **The thread gate of Phase 2, frozen here (*measured*):** at key on, the repeats at each thread
  count byte-identical among themselves; the `ptcl3D` records identical between `nthr=1` and `nthr=16`; and for each
  partial, the largest difference between the thread counts over the volume's largest value no larger than the
  base's for that kind of file, 2.465e-8 for the sums and 2.134e-10 for the densities (the larger of the two halves:
  which half takes the larger value is chance). The count of differing values is reported, not gated, because the key
  inserts more entries per particle and so touches more cells.
- **Base readout.** The `PROB_TAB_INSERT` lines with their stages and iterations, the REC3D bench records (the
  reconstruction's per-pass benchmark, `BENCH/REC3D_BENCH_PART<n>.txt`) and each run's command line, wall time and
  memory are in `keep/base_readout/`. Over the whole solve3D run, shared-memory gridding inserted 7,482 particle images
  in 0.97 s of accumulation (86 passes); PCG in shared memory 7,221 images in 4.29 s (68 gridding and 15 PCG passes).
  Stage 3's first table iteration inserts 1.69 poses per particle at a committed weight of 0.72.
- **Exit.** The build records, the run directories, both exclusion lists (empty for the workflows; the PCG smoke's as
  above), the smoke command lines with their times and the thread probe's result are recorded here. Files changed: this
  plan only.

### Phase 1, plumbing (2026-11-14)

The key exists and every structure the enabled path needs is in place, while the key off reproduces the base byte
for byte. Nothing inserts differently yet: the rows, the quadrature and the new insertion loops are Phase 2. Paths are
under `/data2/agent_runs/inpl3D/scratch/keep/` (`keep/` below); the session's journal is
`/data2/agent_runs/inpl3D/state/journal_phase_1.md`.

- **Files changed (16 source files, lines added/removed):** `simple_parameters.f90` (2/0), `_parse.f90` (1/0),
  `_phases.f90` (14/0); `simple_ui_refine3D.f90` (8/0), `simple_ui_solve3D.f90` (9/0); `simple_commanders_solve3D.f90`
  (2/0); `simple_eul_prob_tab_utils.f90` (154/7), `simple_eul_prob_tab.f90` (86/3), `simple_eul_prob_tab_neigh.f90`
  (21/4); `simple_strategy3D_matcher.f90` (37/4), `simple_strategy3D_alloc.f90` (15/1), `simple_matcher_3Drec.f90`
  (24/3), `simple_matcher_memory.f90` (6/2); `simple_rec3D_pcg_strategy.f90` (6/3); `simple_convergence.f90` (3/3);
  `simple_rnd.f90` (1/1); and this plan. `simple_refine3D_strategy.f90` needed no change: the master merges the
  parts' readout records through the generic record routines, which carry the new fields.
- **What it does.**
  - The key `inpl_posterior=yes|no` (default `no`) is declared, bound as a yes/no key and offered on refine3D and on
    solve3D. On solve3D it sits in its own descriptor, not in the inputs it shares with solve3D_chunks.
    `l_inpl_posterior` is derived once in the parameter phases: the key is `yes`, the program is `refine3D`, the
    particles are 3D, and the mode is `refine=prob` or `prob_neigh` with `prob_neigh_mode=state|geom`.
  - The multi-state final assignment sets `inpl_posterior=no` beside `topk=1` when its parent's command line carries
    the key. Without the key the command line stays the base's.
  - The tables keep each refined entry's window-prior centre: the start of its shift refinement, converted to a
    stored shift at the rotation the refinement starts from. It is taken by the shift search's own rule, before the
    search can move to another rotation. That rule's centre is private to the search object, and `src/main/pftc/`
    is out of scope, so the tables apply the same double-precision formula (`prior_window_centre`). Extras carry the
    centre through the insertion table. The committed pose's centre goes into a per-particle array on the table,
    which is allocated only under the key. `get_prior_centre` and `get_shift_prior_sigma2` expose the centre and
    the iteration's prior width.
  - Under the key, the matcher lays the resolved poses out even at `topk=1`. At the first batch it reserves
    `max(1,topk)` blocks per particle of the part (`reserve_inpl`). Each block holds an entry count, a solved bin
    and 64 packed 32-bit entries: the rotation index in the low 16 bits and the point count above. That is
    `max(1,topk)(8 + 4 x 64)` bytes per particle, the term the memory estimate now adds under the key.
  - The matcher allocates the per-thread shifted-particle banks (`prep_strategy3D_inpl_bank`) and computes `m_max`.
    It sums the per-thread counts in thread order and sets them on the table (`set_inpl_counts`) before
    `close_assignment`.
  - The readout record gains `inpl_nrows` and `inpl_committed_view_weight_sum`, written, read and summed with the
    others. A `PROB_TAB_INPL` line follows `PROB_TAB_INSERT` in the 3D report when the budget is set. Without the key
    the record keeps its defaults and no line is printed.
  - `l_margin` takes the key (decision 11).
  - The REC3D bench records "poses inserted" for gridding and PCG alike: each particle image plus each extra of
    positive weight, the count both insertion paths use.
  - The draw family's comment covers both routes.
- **Decisions taken under delegation.**
  - Decision 13 revised (launch ruling R2; text in section 4): the objective is refused only on a refine3D pass the
    key acts on. Children of solve3D and refine3D copy the key and run under `objfun=cc`.
  - Ruling R3 applied once: the PCG strategy's table row gained Phase 1, for its bench line, before the edit.
  - A design review of the diff by a separate reader found no risk to the off state. Its findings were applied:
    - the PCG bench line;
    - the neighbourhood seed fallback clears its centre;
    - the committed centre's getter stops if the centres are missing under the key;
    - the counts are sized by the search's team where the other per-thread counters are;
    - shorter comments.
  - It also suggested refusing `solve3D inpl_posterior=yes objfun=cc` early. That is moot: solve3D sets
    `objfun=euclid` before its parameters are made.
- **Exit evidence** (all on the final Phase 1 code, Release `repo/build` and Debug `repo/build_debug`, both built
  with no compiler warning).
  - **Identity at one thread** (`test=identity` against `base1`, the declared exclusions being none):
    - 6vxx: 186 files the same, 0 differing;
    - 1jxy: 223 files the same, 0 differing. Its second solve3D runs two parts, so the extended readout record
      was written, read and merged.
    - Reports: `keep/identity/cmp_base1_p1/<suite>/identity_report.txt`; runs `keep/identity/p1/` (670 and 804 s).
  - **PCG smoke, key off, `nthr=1`, `SIMPLE_SEED=7`** (`keep/smoke/p1_pcg_t1`, 609 s) against the base's
    `base_pcg_t1_a` with `keep/pcg_identity.sh`: 158 files the same, 0 differing (`keep/smoke/pcg_identity_base_p1.txt`).
  - An earlier pass of the same checks on the code before the review fixes also gave 0 differing
    (`keep/identity/p1_first/`).
  - **Fast gate:** 15 of 15 on both trees (`keep/logs/p1/fastgate_final_*.log`).
  - **Unit suites:**
    - Release: `unit_pftc_align2D3D` 11 s and `unit_reconstruction` 128 s;
    - Debug: 25 s and 133 s;
    - all pass (`keep/p1_units/final_*.log`).
  - **Static gates:** `check_cli_chain`, `check_test_registry`, `check_descr`, `check_guide_imports` and
    `check_tester_linkage` pass (`keep/logs/p1/static_final_*.log`). The first gate run failed `check_cli_chain`:
    it reads a key's choices from the first parenthesised `a|b` group of the declaration's comment, and the comment
    named `state|geom` first. The comment was reworded.
  - **UI validation:** the build's complete UI JSON validation and the fast gate's `unit_ui` pass.
  - **Refusals** (`keep/p1_refusals/final_*.log`):
    - `refine3D ... inpl_posterior=maybe` stops with
      `ERROR! inpl_posterior=maybe: the key takes yes or no` (solve3D the same);
    - `refine3D refine=prob objfun=cc inpl_posterior=yes` stops with
      `ERROR! inpl_posterior=yes needs objfun=euclid: the in-plane posterior is the Euclidean one`;
    - `refine3D refine=greedy objfun=cc inpl_posterior=yes`, a mode the key does not act on, runs.
  - **Key-on sanity, not a gate** (Debug, bounds checks on): one `refine3D refine=prob` iteration with
    `inpl_posterior=yes` at `topk=3` and at `topk=1` completes; no `PROB_TAB_INPL` line yet (no rows);
    REC3D poses inserted 87 (`keep/p1_keyon_sanity`).
- **For Phase 2.** These have no reader yet and must be used or deleted:
  - the entry packing helpers;
  - `prob_insert_poses%block`;
  - the centre getters;
  - `s3D%inpl_bank`;
  - the per-thread counts the search must add into (`add_inpl_counts`).

  The `PROB_TAB_INPL` committed view weight is divided by the table's particles; a particle without a row or extras
  counts 1.

  With the key on, Phase 1 already widens `l_margin` and reserves memory, so key-on runs before Phase 2 are not
  meaningful.

### Phase 2, the enabled path (2026-11-14)

With `inpl_posterior=yes`, every pose a particle is inserted at now enters the reconstruction at a quadrature of its
in-plane posterior: the committed pose and the weighted insertion's extras, in gridding and in PCG (the
preconditioned conjugate-gradient reconstruction). With the key off, the run is still the base's byte for byte.
Paths are under `/data2/agent_runs/inpl3D/scratch/keep/` (`keep/` below); the session's journal is
`/data2/agent_runs/inpl3D/state/journal_phase_2.md`.

- **Files changed in this phase (lines added/removed against the phase's start):**

  | File | Added/removed | What changed |
  | --- | --- | --- |
  | `simple_strategy3D_prob.f90` | 112/5 | the rows, quadratures, view weights, block fill, counts |
  | `simple_eul_prob_tab_utils.f90` | 98/19 | `inpl_quadrature`'s `lnz`/`dmin`, `inpl_view_weights`, `inpl_prior_penalty`, the grid-angle table and the entry decoders of the reservation |
  | `simple_strategy3D_alloc.f90` | 18/4 | per-thread, per-pose work beside the bank |
  | `simple_strategy3D_matcher.f90` | 19/6 | grid angles for the reservation, counts per particle |
  | `simple_matcher_3Drec.f90` | 61/4 | the key-on gridding loop `insert_inpl_entries`, poses counted under the key |
  | `simple_rec3D_pcg_strategy.f90` | 161/3 | `accumulate_pcg_inpl_entries`, `pcg_inpl_particle_entries` |
  | `doc/policies/3D/refine3D_policy.md` | 47/1 | new subsection "In-plane posterior insertion"; the sentence on the 2D route |
  | `doc/policies/3D/solve3D_policy.md` | 6/0 | the key's reach and the final assignment |
  | `doc/policies/2D/solve2D_policy.md` | 3/1 | pointer to the 3D counterpart |

- **How it works.**
  - After each pose's continuous solve, the strategy evaluates that pose's in-plane row: the Euclidean losses over
    every grid rotation, at the solved increment in the stored frame. It uses the thread's shifted-particle bank, or
    the unshifted row below `SHERRSQ`.
  - The row's 64-point quadrature uses a phase keyed on the particle and the full reference index. Its entries are
    packed into the pose's block of the reservation. The solved bin and a collapsed row are inserted at the pose's
    continuous angle, every other entry at its grid angle. The grid angles are copied from the polar calculator
    before it is released.
  - The view weights are `log m_k = ln Z_k - T L_min,k - P_k`. Under the polar shift map, `P_k` is the map's prior
    on the increment. Under the window search with `shift_prior=yes`, it is the prior on the increment from the
    refinement's retained centre.
  - An extra under 0.01 is dropped, and the committed pose is always kept.
  - If any pose of a particle has no valid row, the whole particle falls back to the insertion without the key.
  - Gridding inserts each entry at `wbase x m_k x count/64`, from the committed plane re-shifted once per extra
    pose. PCG expands the same entries in the same order.
  - PCG lays the entries out per observation batch in chunks of at most `MAXIMGBATCHSZ`. Each chunk has its own
    records, `prep_particles`, plane copies and `accumulate_batch` from index 1.
- **Decisions taken under delegation.**
  - Decision 10 revised (launch ruling R2; text in section 4). The counts are kept per particle and summed in
    particle order, because sums per thread depended on the dynamic schedule (found by the review below).
  - Ruling R3 applied: the row of `simple_strategy3D_alloc.f90` gained Phase 2 before the edit. The per-pose
    results of a particle's rows need per-thread work arrays sized by the team, and the strategy object, made per
    particle, cannot hold them without an allocation per particle.
  - A design review of the Phase 2 diff by a separate reader found no risk to the off state. Its findings were
    applied:
    - counts per particle;
    - a particle without blocks always enters at its committed pose in PCG, and a half that would accumulate no
      entry stops;
    - a zero-sum guard in the view weights;
    - corrected comments.
  - **The thread gate: the probe frozen in Phase 0 missed, and a stricter probe passes. The maintainer should read
    this item.**
    - The frozen probe at key on (`keep/probe/p2on_*`, `keep/probe/thread_gate_p2on.txt`): the three repeats at each
      thread count are byte-identical, and so are the committed poses between the counts.
    - Between one and 16 threads, the largest difference over the volume's largest value is above the base's in
      two files:
      - even sums: 2.773e-8 against 2.465e-8;
      - even densities: 2.402e-10 against 2.134e-10;
      - the odd sums (6.040e-9) and odd densities (2.134e-10) are equal to the base's.
      The key-on run inserts 98 planes for 87 particles where the base inserts 87.
    - The base's spread comes from the noise spectra (sigma2) that each probe run bootstraps from the particles.
      That bootstrap differs between thread counts in 89 bytes (Phase 0) and runs before the search.
    - A probe starting from a finished single-thread iteration reuses that iteration's committed sigma2, so both
      thread counts start from identical inputs (`keep/probe/sig_template`, `keep/probe_fixed_sigma.sh`,
      `keep/probe/thread_gate_fixed_sigma.txt`). In it, the four partial reconstructions are byte-identical across
      all six runs of each arm, `nthr=1` pinned and `nthr=16`, three repeats each:
      - the base build with the key off;
      - the run's build with the key on, exercising the path at 1.60 rows and 10.51 entries per particle and 914
        poses inserted.
    - By the plan's own rule (section 7: when the base's partials are byte-identical across thread counts, the
      key-on partials must be too), the key-on path passes. The gate is judged on this probe. That is a stricter
      criterion, not a tolerance relaxed.
    - The same fixed-input probe also checks the contract that the search does not depend on the key: from identical
      references, the base with the key off and the run with the key on commit the same 87 poses. Only `cfar`, a
      map statistic written after assembly, differs; the committed sigma2 is identical.
    - The Phase 0 probe design is revised accordingly and listed among the report's deviations.
- **Exit evidence** (final Phase 2 code, Release `repo/build` and Debug `repo/build_debug`, both built without a
  compiler warning).
  - **Identity at key off** (`test=identity` against `base1`; reports in `keep/identity/cmp_base1_p2/<suite>/`):
    - 6vxx: 186 files the same, 0 differing;
    - 1jxy: 223 the same, 0 differing.
  - **PCG smoke at key off, `nthr=1`:** 158 the same, 0 differing against the base (`keep/smoke/pcg_identity_base_p2.txt`).
  - The same checks on the code before the review fixes also gave 0 differing (`keep/identity/p2_first/`).
  - **Fast gate:** 15 of 15 on both trees (`keep/logs/p2/fastgate_final_*.log`).
  - **Unit suites:** `unit_pftc_align2D3D` and `unit_reconstruction` pass on both trees (`keep/p2_units/`).
  - **Smoke at key on** (solve3D on the fixture, the command lines of Phase 0 plus `inpl_posterior=yes`; logs and
    readouts in `keep/smoke/p2on_*` and `keep/p2_readout/`). Every run completes:
    - gridding, shared memory, Release (85 s) and Debug with bounds checks (267 s, no runtime error);
    - gridding, `nparts=2` (120 s);
    - PCG, shared memory (133 s) and `nparts=2` (173 s).

    In all of them `PROB_TAB_INPL` appears from stage 3 on (12, 15, 4, 3, 3 and 8 lines in stages 3 to 8) and never
    in stages 1-2. The two-part masters print the parts' counts summed: 81 particles in stage 3, the same as shared
    memory.
  - **The thread gate,** as above.
- **Reported numbers, key on against key off** (same Release build, sequential, load 4-10 from these runs alone;
  `keep/p2_readout/summary.txt`):

  | Route | Wall time | REC3D accumulation | Images / poses inserted | FSC=0.143 final |
  | --- | --- | --- | --- | --- |
  | Gridding, shared memory | 85 s on, 90 s off | 1.34 s on, 0.97 s off | 7,656 / 12,248 on, 7,482 / 7,680 off | 4.71 A on, 4.89 A off |
  | Gridding, two parts | 120 s on, 127 s off | 1.97 s on, 1.42 s off | 12,701 poses on, 7,777 off | 4.89 A on, 4.80 A off |
  | PCG, shared memory | 133 s on, 130 s off | 4.72 s on, 4.30 s off | | 4.46 A on, 4.38 A off |
  | PCG, two parts | 173 s on, 175 s off | 9.88 s on, 10.17 s off | | 4.38 A both |

  - Stage-3 FSC=0.143 is 10.85 A in every run.
  - First stage-3 iteration, key on:
    `PROB_TAB_INSERT` 1.69 poses per particle, committed weight 0.72;
    `PROB_TAB_INPL` 81 particles, 1.74 rows and 6.68 entries per particle, 1.4 % of rows collapsed, mean
    `n_eff` 2.9, committed view weight 0.72 against the table's 0.72.
  - The last stage's rows mostly collapse: 92 % collapsed, 1.08 entries per particle.
  - Peak resident memory: 2.33 GB on against 2.28 GB off in shared memory; the largest part 2.11 GB on against
    2.10 GB off. The estimate's reservation is `max(1,topk)(8 + 4 x 64)` = 792 bytes per particle: 69 kB for the
    fixture's 87 particles and about 35 kB per part.
  - The 87-particle fixture cannot say whether maps improve (decision 2).

### Phase 3, the assertions (2026-11-14)

Every check of section 7 now exists. Each was made to fail by a violation planted in the code it guards, then
passed again with the plant removed. Paths are under `/data2/agent_runs/inpl3D/scratch/keep/` (`keep/` below);
the session's journal is `/data2/agent_runs/inpl3D/state/journal_phase_3.md`.

- **Files changed in this phase (lines added/removed against the phase's start):**

  | File | Added/removed | What changed |
  | --- | --- | --- |
  | `simple_eul_prob_tab_tester.f90` | 269/1 | five tests |
  | `simple_strategy3D_inplane_tester.f90` | 196/2 | two tests |
  | `simple_rec3D_strategy_tester.f90` | 307/2 | two tests |
  | `simple_convergence_tester.f90` | 61/0 | one test |
  | `simple_matcher_memory_tester.f90` | 13/0 | the estimate's term |
  | `simple_eul_prob_tab.f90` | 23/11 | refactor |
  | `simple_eul_prob_tab_neigh.f90` | 3/14 | refactor |
  | `simple_eul_prob_tab_utils.f90` | 27/1 | refactor |
  | `simple_strategy3D_prob.f90` | 3/7 | refactor |
  | `simple_matcher_3Drec.f90` | 7/6 | refactor |
  | `simple_rec3D_pcg_strategy.f90` | 52/33 | refactor |
  | this plan | 8/8 | file-table rows |

  `simple_strategy3D_srch.f90` and `simple_convergence.f90` carried plants only; both are back to their state at
  the phase's start.
- **Deviation: behaviour-neutral production refactors (the plan's change table calls Phase 3 "tests-only").** A check
  can catch a violation planted in production code only if it runs that code, and the code each check needed was
  private or reached only through a file read. Under launch ruling R3, rows were added to the file table before the
  edits. The refactors:
  - the seeded window refinement of both tables, the start, its prior centre and the shift search together, became
    one routine (`seeded_window_refine` in `simple_eul_prob_tab.f90`), which the neighbourhood table calls too;
  - the particle's readout counts became `inpl_particle_counts` in the utilities, which the strategy calls;
  - the gridding loop `insert_inpl_entries` became public and takes the point group instead of the builder;
  - the PCG chunk loop became `accumulate_pcg_inpl_batch`, public, over planes in memory, with the chunk capacity as
    an argument (production passes `MAXIMGBATCHSZ`).

  They change no result:
  - the key-off identity against `base1` holds after them: 6vxx 186 files the same and 1jxy 223, 0 differing; the
    PCG smoke 158 the same, 0 differing (`keep/identity/cmp_base1_p3/`, `keep/smoke/pcg_identity_base_p3.txt`);
  - the fixed-input key-on probe of Phase 2 gives the same four partial reconstructions byte for byte, and the same
    `PROB_TAB_INPL` line, as before them (`keep/probe/p3_fs_on_t1` against `keep/probe/fs_on_t1_1`).
- **Deviation: placement.** The key-on gridding check is in `simple_rec3D_strategy_tester.f90`, in the same
  `unit_reconstruction` entry, not in `simple_reconstructor_tester.f90`. The reconstructor tester is a submodule of a
  volume module, and importing the matcher's reconstruction module from it would make a lower layer depend on a
  higher one. The prior-centre capture check is in the in-plane tester (`unit_pftc_align2D3D`, suite
  `refine3d_in_plane_state`), which holds the polar-calculator fixture a real shift search needs.
- **The checks and their plants** (each plant is a patch in `keep/plants/`, indexed in `keep/plants/README.txt` with
  the suite to run; its failures are in `keep/plants/<plant>.run.log`):
  - **The frame of the solved increment** (`test_solved_row_frame`). On a white reference, rotated and shifted
    about the particle's axes, with a non-zero stored shift, the row at the increment that both continuous solves
    return is checked three ways:
    - its minimum lies within one bin of the solved angle;
    - its value at that bin equals the shifted-reference kernel at `shift_ptcl2ref(sh, rot)`. The measured spread on
      eight solves is 7.6e-8 in Release and 1.5e-7 in Debug, relative; the tolerance is frozen at 1e-6;
    - it beats the rows at the negated increment and at the increment rotated either way by the solved angle.

    Plants, in `strategy3D_srch`:
    - the joint solve's increment returned unrotated: 11 failures, among them "the row's minimum lies within one bin
      of the solved angle" and "the returned increment's row beats the negated and the rotated frames";
    - its sign reversed: 14 failures.
  - **The view weights** (`test_inpl_view_weights`). Rows with one angle, flat, and two peaks, against
    `exp(ln Z - T L_min) / sum` computed in double precision, within 1e-6 relative. Also checked:
    - the committed pose is kept below the 0.01 floor (`PROB_INSERT_WMIN`), and an extra under it is dropped with the
      rest renormalised;
    - the entries sum to 1, and a lone pose takes weight 1;
    - the worked example: 0.4, 0.4 and 0.2, and a third each when the likelihoods are equal.

    Plants:
    - `ln Z` omitted: 11 failures, "the view weight equals exp(ln Z - T L_min) / sum" and the worked example;
    - the floor not re-applied: 4, "an extra under WMIN is dropped";
    - the committed pose dropped: 7, "the committed pose is kept below WMIN".
  - **The shift prior** (`test_inpl_prior_weights`). Identical rows at the prior's centre and two prior sigmas from
    it take 0.881 and 0.119, under the polar shift map and under the window with each pose's own centre; without a
    prior they take 0.5 each. Plants:
    - the prior omitted: 4 failures;
    - the window prior centred at the stored shift: 2 failures, "the window's prior: the pose at its centre takes
      0.881".
  - **The prior centre** (`test_seeded_window_centre`, through `seeded_window_refine`, a real shift search):
    - the seed (0.9 trs, 0.9 trs) at 0 degrees is its own centre;
    - at 45 degrees the start is clipped and the centre is (trs, trs)/sqrt(2);
    - a refinement whose final rotation differs from its start keeps the start's centre (the test asserts the change
      of rotation).

    Plants:
    - the unchanged seed reused as the centre: 2 failures;
    - the centre taken at the final angle: 3, "the centre is the start's, not the final rotation's".
  - **The layout and the reservation** (`test_inpl_reservation_layout`, `simple_matcher_memory_tester`):
    - without the reservation nothing is allocated, and the layout's counts, first slots and weights are the table's
      either way;
    - there are `max(1,topk)` blocks per particle, and their bytes equal the memory estimate's term, which the memory
      tester pins as the estimate's exact increase under the key (792 bytes per particle at `topk=3`);
    - a batch of flat rows fills every block with 64 entries;
    - fills of 0 to 64 entries decode to the rotations and counts written.

    Plants:
    - the count packed into the low bits: 3 failures, "fills of 0 to K entries decode to the rotations and counts
      written";
    - the reservation sized by `topk - 1` blocks. Its first run corrupted the heap, because the test went on filling
      past the undersized blocks, so the check was fixed to stop once the block count is wrong. It then failed
      cleanly: "max(1,topk) blocks per particle expected=15 actual=10" and "the reserved bytes are the memory
      estimate's term".
  - **The committed view weight against the committed angle's weight** (`test_inpl_view_against_angle_weight`). At
    `topk=1` with a flat row, the view weight counted is 1, while the committed angle's inserted weight is 1/64. Plant:
    the field fed with the angle's weight, 1 failure, "the committed view weight is 1 expected=1.00 actual=0.0156".
  - **The readout round trip** (`test_inpl_readout_round_trip`). All in-plane fields survive a part's record, and the
    master's merge sums them. Plant: `inpl_nrows` missing from the sum, 1 failure.
  - **The key-on gridding loop** (`test_inpl_gridding_entries`). On a committed pose and two extras at their own
    shifts, the loop equals every entry inserted on its own from the committed plane re-shifted for that entry, bit
    for bit. Entries at one pose whose weights sum to 1 equal one insertion at weight 1, within 1e-6. Plant: the first
    extra's shifted plane carried over to the next pose, 1 failure, "the key-on loop equals each entry inserted on its
    own ... bit for bit".
  - **The chunked PCG accumulation** (`test_pcg_inpl_chunks`). Chunks of three entries, in one observation batch and
    in two, give the raw B and D accumulators of one flat accumulation of the same entries byte for byte, at one
    thread and at four. The flat list holds one entry per inserted pose entry. Plant: off by one at the chunk boundary
    (a chunk's last entry dropped), 4 failures, at both thread counts.
  - **The margin** (`test_posterior_exits_inpl_posterior`). At `topk=1`, the key takes the weighted insertion's
    margin: exit at the fourth iteration, and none at 7.2 A against an 8 A limit. Without the key the single draw's
    rules hold: exit at the third. Plant: the margin ignoring the key, 2 failures.
- **A false start, for the record.** The first batch of plants ran the previous binary: the driver built but did not
  install, and the tests read `build/bin`. It reported nothing caught. The driver was fixed to install and every plant
  was rerun. None of those first results is used.
- **Exit evidence:**
  - every plant failed its check and was reverted. The planted files equal their clean state, and all 16 patches
    replay on the clean tree (`patch --dry-run`);
  - the suites pass on both trees: `unit_pftc_align2D3D` (6,304 assertions), `unit_reconstruction` (975) and
    `unit_cart_align3D` (1,488), Release and Debug (`keep/p3_units/final_*.log`). The five touched suites pass again
    after the last rebuild (`keep/p3_units/end_*.log`);
  - `scripts/check_test_registry.py` passes;
  - the fast gate gives 15 of 15 on both trees (`keep/logs/p3/fastgate_end_*.log`);
  - the testers' earlier content is unchanged apart from the new checks and their entries in each `run_all_*` list.

### Phase 4, handoff (2026-11-14)

The report is `inpl_quadrature_insertion_3D_report.md`, beside this plan, in the order of
`doc/how2s/how2write_agent_plans.md` section 8:
- what was learned;
- the changes phase by phase, pointing to the driver's phase diffs;
- the decisions taken under delegation;
- the deviations, the thread gate's revised probe first;
- the evidence against the base;
- what is left open, with the maintainer's runs, their command lines and the readout to read;
- the proposed skill lines;
- the delivered file list;
- a commit message.

This plan stays in `planned/` because the maintainer's scale runs follow delivery. Its status line now reads
"delivered, not validated at scale", and no README row changed.

Files changed in this phase: this plan (the status line and this entry) and the new report. No source file changed,
and both build trees were up to date: nothing recompiled (`keep/logs/p4/make_check*.log`). The fast gate passed 15 of
15 on the Release and the Debug tree (`keep/logs/p4/fastgate_*.log`), and the tree is left uncommitted with its
changed and new files listed in the report.
