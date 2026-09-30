# flex_pca-style class splitting (`cls_split pca_mode=flex`): plan

Status: plan, 2026-09-30 (revised the same day: no even/odd machinery; fixed `ncls` subclasses and
fixed `neigs` rank per parent class; EM fits the basis, the subclasses are placed once afterwards).
Target branch: master. Nothing implemented yet.

## 0. What cls_split does today, and the answer to "does it weight?"

`cls_split` (`src/main/strategies/parallelization/simple_cls_split_strategy.f90`, commander in
`simple_commanders_denoise.f90`) runs per parent 2D class:

1. `transform_ptcls` (`simple_classaverager_restore.f90:780`): reads the class members, optionally
   phase-flips, shifts and rotates each particle into the class frame by KB gridding, returns
   real-space images and their plain average.
2. automask from the parent average, noise normalisation, soft mask.
3. `make_pcavecs` -> diffusion-map or kPCA embedding (`make_split_embedding`), ICM rank selection.
4. Euclidean distance matrix in the embedding -> k-medoids, k chosen by silhouette over
   `nsubcls_min..nsubcls_max` (or fixed `ncls`). Every class is split into at least 2.
5. Hard labels written to `ptcl2D/class` (global subclass) and `ptcl2D/cluster` (parent);
   `os_cls2D` rebuilt with `cluster`/`pop`.
6. A trailing distributed `make_cavgs` regenerates `cls_split_cavgs.mrc` from the hard labels.

**It does not build weighted class averages.** `cavger_update_sums` calls `accumulate_fplane`
without the optional `weight` argument (`simple_classaverager_restore.f90:342,345`), so every member
enters its subclass at unit weight. The embedding also has no noise or CTF model: phase-flipped,
masked real-space pixels with a Euclidean/kNN affinity. Sub-class membership is hard and
silhouette-driven, and a homogeneous class cannot answer "do not split".

## 1. Goal

A per-class **projection-aware low-rank covariance fit** (the flex_pca model, in 2D with the
in-plane pose fixed) that yields, per parent class:

- a latent embedding `z_i` with posterior precision, from a CTF- and sigma2-weighted likelihood;
- a fixed rank `neigs` (existing cls_split key); no rank gate;
- exactly `ncls` subclasses per parent (existing cls_split key): k-means centres on the final
  embedding, soft weights `w_is` from a tied-covariance GMM with those `ncls` components
  (`gmm_state_weights`, `nstates=ncls`); no auto placement, pruning or merging;
- **weighted sub-class averages**: `A_s = sum_i w_is C_i^* y_i / (sum_i w_is |C_i|^2 + eps)`; the 2D
  analogue of `reconstruct_flex_weighted_states`. No even/odd halves: the classes are small and the
  halves only served the rank gate and a per-average `res`;
- hard labels (argmax) kept for the existing project contract, and the weight table delivered
  through the `flex_weights` store so downstream tools can use the soft assignment.

No new user-facing keys: `ncls` and `neigs` are the two inputs, both already on cls_split. Poses are not refined
(the low-rank family is a poor aligner; flex2D branch, proposal section 15.1).

## 2. Model (per parent class, Cartesian half-plane, pose fixed)

For member `i` in the class frame (rotated by `-e3`, shifted by `-shift`, as `transform_ptcls`
does, but **without** phase flipping and kept in Fourier space):

```
y_i(h,k) = C_i(h,k) [ mu(h,k) + sum_q U_q(h,k) z_iq ] + n_i(h,k),   n_i ~ N(0, sigma2_i(|k|))
```

`C_i` is the particle's CTF evaluated on the rotated lattice: same defocus, astigmatism angle
rotated by `e3` (sign to be pinned by a one-particle test against explicit rotation). Weights
`W_i(h,k) = |C_i|^2 / sigma2_i(k)`.

E-step (flex_pca's MAP embedding, `_em_embed`): `A_i = I + U^H W_i U`, `b_i = U^H W_i (y_i/C_i - mu)`
(formally: `b_i = U^H C_i^* (y_i - C_i mu)/sigma2`), `m_i = A_i^-1 b_i`, `S_i = A_i^-1`.

M-step (flex2D branch, proposal section 6, transplanted from polar to Cartesian): for every
coefficient `(h,k)` a `(K+1)x(K+1)` weighted least-squares system

```
G(h,k) = sum_i W_i [1 m_i^T; m_i  S_i + m_i m_i^T],   rhs(h,k) = sum_i W_i (y_i/C_i) [1; m_i]
[mu; U] = G^-1 rhs        (Tikhonov ridge on the basis block, relative to G(0,0))
```

over all members, `K = neigs` fixed.

The fit is expectation maximization (PPCA EM, flex_pca's probe loop): `mu` starts from the plain
weighted average and `U` from random projections; prior-free (ALS) probe iterations first, then the
PPCA E-step/M-step alternation with rescaling (flex2D findings: random-projection init is 1000x too small; ALS is
scale-invariant so `U` is rescaled by the latent second moment). Iteration count: flex_pca's
`n_probe_iters` upper bound with the convergence stop.

Fit band: coarse. State differences are low-resolution features and the fit degrades as the band
marches (flex2D 15.11: labels held at 20-30 A beat 8 A). Default `lp_fit` derived, not typed:
the coarser of flex_pca's derived band rule and the parent's own `res` from `os_cls2D`; the
restoration uses the full band. Fit at `box_crop` matched to `lp_fit`, restore at `box`.

Sigma2: the project's canonical sigma2 state (what `make_cavgs ml_reg=yes` reads via
`cavger_read_euclid_sigma2`); if absent, run the same `ensure_canonical_sigma_state` path flex_pca
uses (move it out of `simple_commanders_flex_pca.f90` into a shared helper) rather than a silent
white fallback (memory: the white-noise fallback invalidated real-data runs).

## 3. States, weights, averages

- Placement, once, after the basis EM has converged: `kmeans_latent_targets` (`simple_flex_pca_targets.f90`)
  with `ncls` centres on the MAP embedding, then `gmm_state_weights` (`simple_flex_pca_gmm.f90`)
  with `nstates=ncls` and those centres for the responsibilities. The subclasses are not part of the
  EM (no mixture prior on `z`); `nsubcls_min`/`nsubcls_max` are ignored on this path and `ncls>=2`
  is required.
- Weighted averages: per state, numerator `sum_i w_is C_i^* y_i`, denominator `sum_i w_is |C_i|^2`
  (+ Wiener eps from sigma2 when ml_reg), on the full-box lattice, in the class frame (no
  re-gridding: the members were already rotated once). `res` inherits the parent's `os_cls2D` value.
  Members below the responsibility floor (`RESP_FLOOR`, as flex_pca) contribute nothing.
- Delivery: `cls_split_cavgs.mrc` in global subclass order, registered with `add_cavgs2os_out` like
  today; `os_cls2D` rows with `cluster`, `pop` (hard), `neff` (soft), `res` (parent's), `state=1`; `ptcl2D/class` = argmax subclass, `ptcl2D/cluster` = parent (unchanged contract).
  Weights: one `flex_weights_state_NNN.bin` per global subclass via `flex_weights_deliver`
  with `os_ptcl2D` (the store is keyed on box/smpd/layout digest, so a 2D field works; check that
  `flex_weights_consumable` does not assume a `vol_flex` sibling; if it does, add a `which_imgkind`
  argument rather than a second store).

## 4. Architecture and ownership (master)

The flow is the repo's: `ui -> exec -> commander -> strategy -> domain`. Nothing new is added
above the strategy; the new code is one domain module plus one extension of the class subsystem.

| layer | file | owns | must not |
|---|---|---|---|
| ui | `src/main/ui/simple/simple_ui_cluster2D.f90` (`new_cls_split`) | `pca_mode` gains `flex`; help text says `ncls` and `neigs` are required on that path | anything else |
| exec | `src/main/exec/simple_exec_denoise.f90` | unchanged: `cls_split` -> `commander_cls_split` | |
| commander | `src/main/commanders/simple/simple_commanders_denoise.f90` (`exec_cls_split`) | cmdline normalisation and defaults before `params%new`; the canonical sigma2 state on the master only (workers consume it); the four lifecycle calls; skipping the trailing `make_cavgs` when `pca_mode=flex` | reading part files, numerics |
| strategy | `src/main/strategies/parallelization/simple_cls_split_strategy.f90` | roles (shmem/master/worker), class partitions, the part-file protocol, the merge, all project writes (`ptcl2D`, `os_cls2D`, out-segment registrations of the cavg stack and the weight files). `split_one_parent_class` becomes a dispatcher: the existing body moves to `split_class_diffmap`, the new `split_class_flex` prepares the inputs (planes, CTFs, sigma2 rows via the class and sigma2 subsystems) and calls the domain | the model, the EM, the averaging |
| domain | `src/main/flex/simple_flex_cls_split.f90` (new) | `flex_cls_model` and its EM (`fit`), `embed`, `place_states` (k-means + fixed-K GMM through `simple_flex_pca_targets` / `simple_flex_pca_gmm`), `restore_states` (weighted averages). Inputs are arrays: complex planes, CTF planes, sigma2 rows, `ncls`, `neigs`. Outputs: `z`, `weights(ncls,N)`, `labels(N)`, `cavgs(ncls)` as Fourier images | importing `parameters`, `builder`, `sp_project`, `qsys_env`, or any part I/O (the 3D driver `simple_flex_pca_model` does; the 2D module stays pure so the unit test needs no project) |
| class | `src/main/class/simple_classaverager_restore.f90` (`transform_ptcls`) | the one gridding kernel that puts a member into the class frame; gains an optional Fourier-plane output (unflipped) and the rotated CTF parameters per member. The real-space path is unchanged | knowing about flex |
| sigma2 | `src/main/sigma2/simple_sigma2_bootstrap.f90` | `ensure_canonical_sigma_state` lifted here from the flex_pca commander. Today three private copies exist (flex_pca commander, rec3D strategy, cluster2D strategy `prepare_canonical_sigma_update`); cls_split calls the lifted one, the other copies can converge later, out of scope | |
| weights store | `src/main/flex/simple_flex_weights_state.f90` + `sp_project%add_flex_weights2os_out` | reused as is with `os_ptcl2D`; if `flex_weights_consumable` turns out to assume a `vol_flex` sibling, the fix is a `which_imgkind` argument in the flex module, not a project change | a second store |
| test | `production/tests/simple_test_flex_cls_split.f90` (new) | exercises the domain module alone on synthetic planes (auto-globbed, no CMake edit) | project fixtures |
| docs | `src/main/flex/README.md`, this note | | |

Reuse from the 3D flex code, exactly this and no more: `kmeans_latent_targets`
(`simple_flex_pca_targets.f90`), `gmm_state_weights` (`simple_flex_pca_gmm.f90`), the small
Hermitian solver `solve_real_spd_complex` (`simple_flex_reconstructor_latent_ops.f90`), and the
`flex_weights` store. Not reused: the EM driver, E-step, M-step and their submodules (all bound to
volume projection and the polar bank), `simple_flex_pca_rec3D`, `simple_flex_pca_rounds` and the
flex_pca strategy (its stage/round scheme and producer-owned part I/O are what this design avoids).
The per-coefficient EM in 2D is ~300 lines of new code; the flex2D branch's
`simple_flex2D_model.f90` is the reference for it and is not imported.

Reused unchanged elsewhere: `discrete_read_imgbatch`/`prepimgbatch`, the KB gridding loop,
`ctf%eval`, `load_sigma2_groups`, `add_cavgs2os_out`, the cls_split master/worker scheduling and
merge.

## 5. Parts: the standard SIMPLE distribution, not flex_pca's rounds

The existing cls_split strategy already follows the usual pattern and is kept as is:

1. The master (`nparts>1`, no `part=`) partitions the parent classes into `nparts` lists
   balanced by population (`prepare_class_partitions`), writes `cls_split_classes_partNN.txt`,
   and submits `simple_exec prg=cls_split ... part=NN class_assignment=<list>` through `qsys_env`.
2. Every worker runs the same shared-memory code (`run_local_split`) on its classes, in the
   master's directory, and writes only part files: today `cls_split_class_map_partNN.txt` and
   `cls_split_assignments_partNN.txt`; the flex path adds `cls_split_weights_partNN.bin`
   (pind, parent, local sub, `ncls` weights) and `cls_split_cavgs_partNN.mrc` (local sub order).
   A worker never touches the project.
3. The master waits on the job-finished flags, merges the part files once
   (`merge_worker_outputs`: global subclass numbering, then labels, weights and the concatenated
   cavg stack in that order), writes the project, registers the stack and the weight files, and
   cleans up. One submission, one merge, no rounds, no stages, no master-side numerics.
4. `nparts=1` runs the same `run_local_split` in process and skips the part files.

The strategy owns every part file, read and write, as it does today. Cost per class per EM
iteration is `N * ncoeff * (K+1)^2` (N = 1000, box_crop 128, K = 8: ~7e8 flops), negligible
against the particle read, so a class is never split across workers.

## 6. Validation (batch2, never the Mac)

1. Unit test (item 7), no data.
2. One-particle convention test: rotated-CTF sign, plane conventions against `transform_ptcls`'
   real-space output (inverse FFT of the returned plane must equal the returned image).
3. 10028 8000-particle reference run (`/mnt/beegfs/elmlund/afan/10028/flex2D_tests/ref.simple`):
   `ncls=2 neigs=4` on every class; do the two sub-averages differ visibly, and how do their
   weights split (a 50/50 split on a homogeneous class is the null result to expect).
4. Synthetic 50/50 mix (EMD-8440 vs 8445, `t4_mix5050`): abinitio2D's 60 classes are all ~50/50;
   score per-subclass purity after `cls_split pca_mode=flex` against the existing
   `pca_mode=diffusion_maps` on the same run. This is the decision test.
5. EMPIAR-10076 with the published per-particle labels (`labels/real_recall.py`): purity/NMI per
   subclass vs the diffusion-map path.

## 7. Risks

- **View vs state.** Within a parent class the residual out-of-plane spread can be the first
  dominant component, and a fixed `ncls` split will then partition on view. The current cls_split
  has the same confound. Mitigation: report the component's dominant resolution band and leave the
  split, note the confound in the manifest.
- **Small classes.** With `ncls` and `neigs` fixed, a 40-member class is fitted and split whatever
  its content; classes below `3 x (neigs + 1)` members are skipped (as cls_split skips pop <= 2).
- **Rotated-CTF sign and the gridding conventions**: covered by validation step 2 before anything
  else is run.
- **Weight store consumers** were written for `vol_flex` siblings; verify item 3's delivery does not
  break `flex_weights_consumable` for 3D runs.

## 8. Open decisions

- `cls_split pca_mode=flex` (recommended: reuses scheduling, merge, and the project contract) versus a
  separate `flex_cls_split` program (cleaner ownership, duplicated master/worker plumbing).
- Whether the weighted averages should also be re-centred per subclass (cluster2D centres its
  averages every iteration; cls_split does not).
- Whether to expose the soft weights to `map2ptcls`/selection downstream, or leave the hard labels
  as the only consumed output for now.
