# refine3D_states_cavgs: design note

**Date:** October 8, 2026
**Status:** Design agreed on October 9, 2026; implementation starts with the
controller extraction. Revised on October 8 after an independent
review (Codex), and on October 9 after Hans's rulings on flex initialization,
the PCG backend, the flex settings (flex keeps its own cross-validation and its
code is not changed) and the handling of the even and odd class averages,
which follows `solve3D_cavgs` exactly. Written from a read of GitHub
master at commit `89bbc49` and of the `refine3D_states_cavgs` worktree.
Nothing was compiled or run.

## Summary

`refine3D_states_cavgs` is a proposed program that sorts the sub-class averages
written by `cls_expansion` into conformational states. It is a careful analysis
of the dynamics that `cls_expansion` resolves, aimed at high-quality
class-average state maps. It is not an initial-model step.

The program works on a temporary class-average work project and treats the
even and odd averages exactly as `solve3D_cavgs` does: each sub-class enters as
two independent rows, one per half, and each row is aligned and assigned a
state on its own. At the end the better-scoring row of each sub-class gives
its pose and state.

First `flex_pca` maps the rows to discrete states. Then the `refine3D_states`
frequency march refines those states, exactly as particle-level
`refine3D_states` refines its flex initialization. Reconstruction uses the PCG
backend throughout.

The program uses no nonuniform (NU) machinery and no automasking. Its product
is a state label and a pose per sub-class average in `cls3D`, handed back to
the particles, plus class-average state maps. Particle-level `refine3D_states`
then continues from the handed-back labels; it skips its own `flex_pca` step
and owns the final particle maps.

## Rulings of 2026-10-08 and 2026-10-09

Hans ruled on the following. Each ruling is carried into the sections below.

| Ruling | Consequence in this note |
| --- | --- |
| Handle the even and odd class averages exactly as `solve3D_cavgs` does | Two independent rows per sub-class; the better-scoring row wins at the end |
| Do not change `solve3D_cavgs`; it works very well | `solve3D_cavgs` stays as it is; the work-project builder here is new code |
| Start with flex: flex maps the class averages to discrete states, which the march refines as `refine3D_states` does | Flex initialization is in scope; random seeding is gone |
| Use PCG everywhere; flex already does by default. PCG did not work when tried in `solve3D_cavgs`, for a reason not remembered, but it should be useful on sparse class-average data | `rec_backend=pcg` for flex's basis solve, its state maps, the pose-seeding checkpoint map, the march and the final maps |
| Flex's covariance band is the march's `lpstart` | `lp=lpstart` on the flex step |
| Keep flex's own cross-validation and leave the flex code unchanged (this replaces the earlier ruling to turn the cross-validation off) | Flex runs with its defaults, set only through command-line keys |
| `nstates >= 3` | As for particle `refine3D_states` with flex initialization |
| Flex's mean map is a consensus map reconstructed from the class averages | The checkpoint map of the pose-seeding pass |
| No automasking in this workflow for now | PCG solves on the spherical `mskdiam` support |
| Full box, with whatever box ceiling flex uses by default | Flex has no ceiling of its own below the native box, so the native box is used |
| Enlarge the phantom | The FLEX phantom fixture gets size and state-count arguments |

How `solve3D_cavgs` handles the halves, which this program copies
(`simple_commanders_solve3D.f90`, lines 205-224 and 270-287;
`solve3D_cavgs_policy.md`, sections 3 and 6):

- The even and odd averages of each class are two rows of the work project,
  each in its own half.
- Every row is aligned and, with `nstates > 1`, assigned a state on its own.
  The even half map of a state is built from the even averages in it, the odd
  half map from the odd ones, so a class whose halves sit in different states
  feeds half its data to each.
- For each class, the pose, state and score of the better-scoring row go into
  `cls3D`. A disagreement between the halves is not resolved; it is measured
  and logged as the even/odd state overlap and the even/odd angular distance.

No record says why PCG failed in `solve3D_cavgs`. The code shows only that the
`solve3D` controller keeps every class-average stage on gridding
(`set_refine3D_backend_policy` in `simple_solve3D_controller.f90` sets
`gridding` whenever `l_cavgs` is true). The validation below compares the two
backends on this program, and `rec_backend=gridding` stays accepted as an
override.

## Scope and non-goals

The base `refine3D` and reconstruction code is not touched.

In scope:

- A new program `refine3D_states_cavgs`, wired through ui, exec and commander.
- A states controller extracted from the `refine3D_states` commander. It owns
  the march and the flex initialization, each under a typed policy.
- A class-average work-project builder.
- A strict-full option on the pose initializer.
- A project routine for the hand-back of states to the particles.
- Input: a project after `cls_expansion` that also holds a consensus map
  (`vol`, state 1) and particle poses from a prior 3D refinement.
- Output: state, pose and score per sub-class average in `cls3D`; a state label
  per particle in `ptcl3D`; class-average state maps registered as `vol_cavg`.

Not in scope:

- Any change to `solve3D_cavgs`.
- Running the base `refine3D` on the `cls3D` segment. The march runs on
  `ptcl3D` rows of a temporary work project.
- Final particle maps. Those stay with particle-level `refine3D_states`; this
  program ships class-average maps, as `solve3D_cavgs` does.
- The weighted mode (`m_estimator`) beyond what flex's own state maps use, NU
  filtering, automasking and solvent masks or priors.
- Any change to the flex production code under `src/main/flex`; `flex_pca` is
  driven through command-line keys only.
- Any change to `cls_expansion`, or any behaviour change in particle
  `refine3D_states` (whose flex initialization keeps gridding and the derived
  working box).
- Moving `classify3D_refs` onto the shared controller. It can follow later.

## Work project

The work project is built exactly as in `solve3D_cavgs` (section 3 of
`solve3D_cavgs_policy.md`). `cls_expansion` already writes the class averages
in that form: `cls_expansion_cavgs.mrc`, which it registers as the `cavg`
stack, with `cls_expansion_cavgs_even.mrc` and `cls_expansion_cavgs_odd.mrc`
beside it.

- It preserves project info, compute environment and job-process metadata.
- It registers the even and odd class-average stacks as two stacks.
- It expands each class average into an even and an odd `ptcl3D` row without
  CTF, carrying class and state. For `N` sub-averages, rows `1..N` are the even
  averages (`eo=0`) and rows `N+1..2N` the odd ones (`eo=1`).
- It owns a transient canonical sigma2 state file and never inherits the input
  project's.
- It copies no out segment, so every map a child needs is passed to it
  explicitly.
- It is deleted, with that file, when the run completes.

Everything downstream then sees an ordinary `ptcl3D` lineage, so Euclidean
matching, `reconstruct3D`, `flex_pca` and `calc_final_rec` run unchanged.

The builder is new code that follows the construction in `exec_solve3D_cavgs`;
`solve3D_cavgs` is not touched and could adopt the builder later. The commit is
not shared either: `solve3D_cavgs` recreates `cls3D`, which here would erase
the `cluster` and `accept` fields that `cls_expansion` wrote.

## Flex initialization on the class averages

Flex runs once, after pose seeding, on the work project's rows. It turns them
into discrete states: hard labels in `ptcl3D` and one map per state. The march
then refines them. This is the particle path's initialization with different
settings, and the controller owns both.

| Setting | Particle `refine3D_states` | `refine3D_states_cavgs` |
| --- | --- | --- |
| Basis solve (`rec_backend`) | `gridding` | `pcg`, flex's own default |
| State maps (`rec_states_backend`) | `gridding` | `pcg`, with flex's solver budget (`maxits_pcg=4`, `rtol=0`) |
| Working box (`box_crop`) | derived from `lp` (flex default) | the native box, set explicitly; `lp` alone would derive a smaller box |
| Covariance band (`lp`) | flex default, 16 A | the march's `lpstart`, 10 A by default |
| Cross-validation | flex default (paired engine) | flex default (paired engine) |
| Mean map (`vol1`) | project consensus map | the pose-seeding checkpoint map |
| PCG solve support | none (gridding) | spherical `mskdiam` |
| State count (`npreimages`) | `nstates`, at least 3 | `nstates`, at least 3 |
| Population floor (`min_state_frac`) | 0.1 | 0.1 |
| Minimum effective state size (`min_neff`) | flex default, 2000 | 20, flex's lower bound |
| Distribution | as configured | shared memory |

Flex labels each row on its own, so the even and odd rows of a sub-class can
start in different states, as they can after the random start of
`solve3D_cavgs`. The march refines every row and the better row wins at the
end.

The mean map is the consensus map the pose-seeding pass reconstructs from the
class averages themselves (step 4 of the workflow). It has the same
observation model as the rows: no CTF and the same restoration. It is passed
explicitly as `vol1`; the reconstruction made it at the native box. The pose
initializer writes it under the state-1 map name
(`refine3D_state_vol_fname(1)`), and flex's delivered maps reuse that name, so
the commander first copies it to a name of its own.

Flex runs in shared memory. Its distributed master hands each worker a
contiguous slice of rows, and a worker refuses a slice that holds only one half.
With the even averages first and the odd ones second, `nparts=2` would give one
all-even part and one all-odd part. The image count is small, and every worker
would hold the full-box M-step volumes anyway, so distribution saves no memory.

The population floor decides the state count. With `min_state_frac > 0`, flex
places `nstates` states and aims to keep each at or above the floor; it warns,
but goes on, when it cannot. It runs neither its two-gate state merge nor its
bandwidth cross-validation. Rows it cannot place (outside every kernel, or in
placement clusters below the floor) get a uniformly random state, logged as
`randomized`, and the march then moves them. The weights it publishes are
one-hot. It then drops every state whose population is below `min_neff`; the
rows of a dropped state go to state 0.

`min_neff` must be set here. Its default is 2000 particles, and flex uses
`max(20, min(nrows, min_neff))`, so with fewer than 2000 rows the minimum is
the row count itself. No state can reach it, and flex then keeps only the two
best-supported states whatever `nstates` asks for. The class-average policy
sets `min_neff=20`, flex's own lower bound, so a state with fewer than 20 rows
is dropped (unless fewer than two states would remain, when flex keeps the two
best-supported). For the 10% floor not to fall below those 20 rows, the run
needs at least 200 rows, which is 100 selected sub-averages at the default
`min_state_frac`. Flex's own minimum (100 active rows, 20 in each half) is
lower.

The maps flex delivers (`reconstruct3D m_estimator=flex`, PCG, native box) are
the march's starting references, as in the particle path. The flex state
weight set is then discarded; the march works on hard labels only.

When `cls3D` already carries multi-state labels, flex is skipped. The labels are
copied to both rows of each sub-class and the startup maps are reconstructed
from them; the work project holds no state maps to reuse.

Flex's M-step memory grows with the cube of the box and with
`ncomp(ncomp+1)/2` coupled component pairs (up to 55 at the default
`neigs=10`), however few images there are: at box 256 it is 64 times what it is
at box 64. Flex logs the GB of the coupled normal matrix and of the PCG kernels
and right-hand sides at the start of every M-step. (Its `(box/64)^3` footprint
note prints only on the derived-box path, so not here.) The per-image work on
class averages is small; the per-iteration volume work is not.

## Flex's own cross-validation on this layout

Flex runs with its own cross-validation, and the flex code is not changed:
every setting in the table above is a command-line key. Under this program's
settings its validation devices behave as follows.

| Device | Behaviour here |
| --- | --- |
| Paired engine: two disjoint fits split by row index modulo 4, matched, merged at the accumulators, each component weighted by `2c/(1+c)` from the cross-fit match `c` (below 0.143 its eigenvalue is scaled by 1e-6), a cross-fit FSC ridge, one polish iteration | on; each fit needs at least 100 rows |
| Even/odd half-sets inside each fit: the half-set FSC sets the Wiener filter of every M-step | on; the halves are the even and odd rows, as in the march |
| Number of mixture components `K` in the latent deconvolution, chosen on two disjoint halves of the rows | `K <= nrows/2000`, so `K = 1` below 4000 rows |
| Bandwidth cross-validation (`nbins > 1`) | off by default, and off under the population floor |
| Two-gate state merge on trial half maps | runs only with `preimage_auto=yes`, never under the population floor |

The paired engine stops when either fit gets fewer than 100 rows
(`run_flex_pca_paired`). That needs about 200 active rows, which the run's
minimum of 100 sub-averages already gives, but only about: the split is by row
index, so the preflight counts each fit's rows with flex's own rule
(`flex_pca_half_of`, public in `simple_flex_pca_stages.f90`) and stops before
any work if a fit would be short.

The split also decides where the two rows of a sub-class go. With the even
averages in rows `1..N` and the odd ones in `N+1..2N`, rows `i` and `N+i` land
in the same fit when `N` is a multiple of 4, in opposite fits when `N` leaves
remainder 2, and half and half otherwise. In the same fit, the two fits see
disjoint sets of sub-classes. In opposite fits, they see the same sub-classes
through their even and odd averages, so their independence rests on the
even/odd split, the assumption this program and `solve3D_cavgs` already make
for the FSC. Either way nothing fails; the cross-fit agreement can read
optimistic, as the FSC does.

The split-half reliability prior in flex's embedding is hard-wired off
(`l_relprior = .false.` in `simple_flex_pca_embed.f90`), so its branches are
dead. Like the dead lines in `exec_solve3D_cavgs`, they are left for a separate
cleanup, because this work does not touch the flex code.

## Evidence contract: FSC as in solve3D_cavgs, no NU

The FSC compares the half map built from the even rows with the one built from
the odd rows, as in `solve3D_cavgs`. It enters the program in these places:

| What the FSC may drive | Where | Setting |
| --- | --- | --- |
| ML regularization of the state maps | PCG: the regularized pair is the closed-form `P_tau` replay from the FSC of the base pair (`shrink_by_ml_prior` in `simple_reconstructor_pcg.f90`; `reconstruct3D_pcg_policy.md`, "ML two-map contract") | `ml_reg=yes`, from the first block |
| Final maps through the ML route, with canonical sigmas reused or bootstrapped | `calc_final_rec` in `simple_final_rec.f90`, PCG branch | follows `ml_reg` |
| Low-pass of the final snapshot at FSC 0.143 | `solve3D_state_fsc_lowpass`, called by `write_final_rec_outputs` | never finer than `lpstop` |
| FSC and resolution records per iteration | volume assembly | always, diagnostic |
| Wiener filter of every flex M-step, within each of flex's two fits | flex probe fit | flex default |

| What stays off | Where it would enter | Closed by |
| --- | --- | --- |
| NU filter bank, NU masks and low-pass promotion | `filt_mode=nonuniform*`: `nonuniform_filter_state`, `refresh_matching_lp_from_project` | `filt_mode=none`, fixed |
| Matching low-pass taken from the FSC | `filt_mode=fsc` or `uniform` sets `l_lpauto`, `simple_parameters_phases.f90` | `filt_mode=none`; the march sets `lp` per block |
| Automasking, envelope-masked FSC and the density-envelope PCG support | `automsk`, `envfsc` | both `no`, fixed; PCG solves on the spherical `mskdiam` support |
| PCG solvent mask and solvent prior | `pcg_mskfile`, `pcg_solvent` | unset, refused |
| Postprocessing of the final maps | `postprocess_states` in `simple_rec3D_service.f90` | `l_postprocess=.false.` |

The second table is what `solve3D_cavgs` already fixes for its class-average
route (`filt_mode=none`, `automsk=no`, `envfsc=no`, no postprocessing), plus
the PCG support settings. The backend is where the two programs differ:
`solve3D_cavgs` stays on gridding.

With the low-pass set per block, the matcher puts the merged state map into
both reference slots (`read_mask_filter_refvols`,
`simple_matcher_refvol_utils.f90`), so both halves are matched against the
merged reference, as in `solve3D_cavgs`. With `ml_reg=yes` that reference is
regularized at assembly, and the block's `lp` limits the matching band.
`simple_refine3D_stage_plan` needs only box, sampling, `lpstart` and `lpstop`.

The FSC is not gold-standard, as in `solve3D_cavgs`: both halves are matched
against the merged reference. With class expansion as it stands today, the
halves also share the fitted basis, the labels, the in-plane alignment and the
Wiener prior, so the curve will read optimistic.

Avoiding a resolution claim does not remove that bias. The FSC is converted to
an SSNR that enters the regularization, so an optimistic curve
under-regularizes the state maps, and inside flex the basis.

Decided on 2026-10-08: ML regularization is on from the first block all the
same, because the march starts from good orientations. The comparison against
`ml_reg=no` stays in the validation list as a check on that choice.

Sibling sub-averages do not share particles, because the delivered stack is
restored from hard labels.

## Workflow

The program builds a class-average work project, seeds poses, lets flex map
the rows to states, runs the shared states controller there, and writes
labels, poses and class-average maps back in one final step.

```text
INPUT PROJECT                                  TEMPORARY WORK PROJECT

cls_expansion output            -- build -->   2N rows, even and odd averages
  sub-average stack, even/odd                    no CTF, own transient sigma2 state
  cls2D selection, ptcl2D class                          |
  consensus vol, ptcl3D poses                  seed poses
                                                 one cc pass against the consensus map
  (not modified while                            checkpoint map from the class averages
   the run proceeds)                                     |
                                               flex_pca, PCG, native box, lp = lpstart,
                                                 its own cross-validation
                                                 hard labels and state maps
                                                         |
                                               refine3D_states frequency march, PCG
                                                 lp set per block, even/odd FSC, no NU
                                                         |
                                               better row per sub-class
                                                         |
written back                    <-- commit --  state maps at native sampling, PCG
  cls3D: state, pose, score                      no postprocessing, low-pass snapshot
  ptcl3D: state only, poses kept
  out: vol_cavg state maps
```

The input project is read once at build and written once at the end; every
FSC the run computes stays in the work project.

1. **Preflight.** The run stops before any work unless all of these hold (the
   flex-fit count runs right after the work project is built, before any
   refinement):
   - the sub-average stack and its `_even` and `_odd` companions exist, share
     box and sampling, and hold the same number of images as `cls2D` and
     `cls3D` have rows;
   - the consensus `vol` of state 1 has the same box and sampling;
   - `ptcl3D` is congruent with `ptcl2D`, and every particle active in
     `ptcl3D` is active in `ptcl2D` (particle poses are not read; they are
     left as they are);
   - `mskdiam` is given (the PCG support);
   - at least 100 sub-averages are selected;
   - with two or more populated `cls3D` labels on the selected sub-averages,
     a given `nstates` equals their number;
   - otherwise `nstates >= 3`; each of flex's two fits gets at least 100
     rows, counted with `flex_pca_half_of`; and flex's population floor fits,
     by flex's own check with its rounding:
     `nstates * nint(min_state_frac * nrows) <= nrows`, with `nrows = 2N`.
2. **Snapshot the particle selection.** The `ptcl2D` states of the input
   project are recorded before anything else, for the hand-back in step 9.
   The `ptcl2D` selection, not `ptcl3D`, so that a rerun on the output starts
   (ruled on 2026-10-09). The rerun takes the labels path, so a sub-average the
   earlier run left at state 0 stays out, with its particles, unless its
   `cls3D` state is reset first.
3. **Build the work project** with the builder: `2N` rows, selection from the
   `cls2D` state flags, its own transient sigma2 state file. All children run
   from a separate work command line with the work project file and
   `mkdir=no`; the public command line is never mutated.
4. **Seed poses.** One fixed-reference `cc` pass aligns every row against the
   input consensus map at `lpstart`, through
   `initialize_poses_against_external_references` in its new strict-full mode.
   The same service bootstraps sigma2 from the images and reconstructs a
   checkpoint map from the aligned rows (PCG). That checkpoint map is the
   consensus map from the class averages that flex uses as its mean. The
   service clears `sampled` and `updatecnt` before its pass; in the
   strict-full mode it clears and writes them again after its validation, so
   the coverage check in step 6 counts the march and not the seed pass.
5. **Map to states with flex.** Flex runs on the work project with the
   settings of the flex section: hard labels in `ptcl3D` (rows of a dropped
   state at 0), one PCG map per state at the native box. Its weight set is
   discarded. The number of states that survive becomes the run's `nstates`.
   A project whose `cls3D` already carries two or more populated labels on the
   selected sub-averages skips flex and copies them, numbered consecutively
   (an empty label is dropped), to both rows of each sub-class; the pose
   seeding then keeps no mean map.
6. **Run the shared states controller** with the class-average policy. It
   starts from the flex maps (or reconstructs startup maps from copied labels),
   runs the frequency march under `pose_policy` in strict full-batch mode on
   the PCG backend, and requires every active row to be updated by the march.
7. **Pick the better row.** For each sub-class, the active row with the better
   score gives the pose, state and score, as in `solve3D_cavgs`; a sub-class
   with no active row gets state 0. The even/odd state overlap and angular
   distance are logged, and so are the state populations. A state that ended
   empty is reported, not treated as a failure (see Decisions).
8. **Reconstruct the final maps** from all active rows as assigned, with the
   ending `solve3D_cavgs` uses on the PCG backend: `calc_final_rec` at native
   sampling without postprocessing (the controller defines `ml_reg`
   explicitly, which its ML route requires), then the raw and low-pass
   snapshot maps, named by the consecutive state numbers of step 9, which are
   assigned right after the pick.
9. **Commit, hand back and clean up, as one final write of the input
   project.**
   - The populated states are numbered consecutively, and the old-to-new
     mapping is logged.
   - `cls3D` in place: each committed sub-average's state, pose and score are
     written into its existing row; the `cluster` and `accept` fields from
     `cls_expansion` are kept.
   - `ptcl3D`: particles not selected in `ptcl2D` stay at state 0; selected
     particles of a committed sub-class take its state, the others go to
     state 0; poses are left as they are.
   - `out`: the `vol_cavg` entries are replaced as a set.

   Then the work project and its sigma2 state are deleted. Flex's diagnostics
   (eigenvolumes, latent readouts, its manifest) stay in the run directory.
   Nothing of flex is registered in the input project: its rows are the work
   project's class averages, not the input particles.

Particle-level `refine3D_states` can then be run on the same project. It finds
populated multi-state labels, reconstructs state maps from them and proceeds
without `flex_pca`.

| Setting | Value | Status |
| --- | --- | --- |
| `filt_mode` | `none` | fixed, anything else is refused |
| `ml_reg` | `yes` | from the first block; `ml_reg=no` turns it off |
| `envfsc`, `automsk`, `combine_eo` | `no` | fixed, other values are refused |
| `objfun`, `sigma_est` | `euclid`, `global` | fixed, as in `solve3D_cavgs` |
| `rec_backend` | `pcg` | seed checkpoint, flex, march and final maps; `gridding` accepted as an override |
| `mskdiam` | none | required, the PCG support |
| `pcg_mskfile`, `pcg_solvent` | unset | refused |
| `nsample` | all active rows | fixed, user values are refused |
| `update_frac`, `cohort_sampling` | unset | fixed, user values are refused |
| `trail_rec`, `trail_seed` | `no`, unset | fixed, other values are refused |
| `pose_policy` | `global` | overridable, `local` or `global` |
| `nstates` | none | required, at least 3; flex may deliver fewer, at least 2 |
| `min_state_frac` | 0.1 | overridable; flex's population floor |
| `lpstart`, `lpstop` | 10 A, 6 A, as `refine3D_states` | overridable; `lpstart` is also flex's band |
| flex working box, `min_neff`, distribution | native box, 20, shared memory | fixed |
| GUI publication | from the input project only | fixed |

All of these live in the typed controller policy, not in command-line state.

## Refactoring steps

Six code changes carry the program. The controller extraction must leave
particle `refine3D_states` unchanged. The flex production code is not changed.

1. **Shared states controller.** The march and the flex initialization now
   live in contained procedures of `exec_refine3D_states`. They move into a
   module, working name `simple_refine3D_states_controller`, with typed
   policies:
   - the controller owns the flex initialization, sampling and frequency
     planning, state-count resolution from labels, startup state maps, stage
     execution, the coverage read, the missing-update pass and the
     state-overlap read;
   - the march policy carries every fixed setting, the backend, the update
     scheme (sampled or strict full batch) and the GUI switch;
   - the flex policy carries the command-line settings of flex: the backends,
     the working box, the band, `min_neff`, the mean map and the distribution,
     with the particle path's values as one instance and this program's as the
     other;
   - each commander keeps what is specific to its lineage and calls the flex
     step and the march in its own order. The particle commander keeps the
     postprocessed ending. The class-average commander keeps build, pose
     seeding, the better-row pick, finalization, state renumbering and commit.

   Large objects such as `parameters` and `cmdline` are allocatable components,
   as the compile-time policy requires.
2. **Class-average work-project builder.** A new module, working name
   `simple_cavgs_work_project`, builds the work project as `exec_solve3D_cavgs`
   does: stack lookup, even/odd companions, project info copy, transient sigma2
   state path, row expansion, and the `nparts` cap at the row count. It also
   provides the even/odd agreement diagnostics (the equivalents of `conv_eo`
   and `conv_eo_states`). `exec_solve3D_cavgs` is not touched. Beyond
   `solve3D_cavgs`, the builder mirrors the selection into the `ptcl2D`
   states, checks that each companion stack holds one average per class, and
   counts the state overlap only over sub-classes with an active row. Its unit
   test is `simple_test_exec test=unit_project suite=cavgs_work_project`.
   (Done on 2026-10-09.)
3. **Strict-full option on the pose initializer.** An optional argument of
   `initialize_poses_against_external_references` sets the sample to all
   active rows, lifts the 100,000-row cap, requires every row to be updated,
   and clears `sampled` and `updatecnt` again after its validation. Without it
   the three existing callers behave as before. The checkpoint map it already
   returns is flex's mean. (Done on 2026-10-09: the optional argument
   `l_all_rows`; the counters are cleared after the checkpoint reconstruction,
   so the checkpoint itself runs exactly as before.)
4. **New commander `commander_refine3D_states_cavgs`.** It owns workflow steps
   1 to 4 and 7 to 9, and calls the controller for steps 5 and 6. (Done on
   2026-10-09. The controller policy gained the flex settings of the flex
   section, the strict full batch, the GUI switch and the iteration offset;
   every new field defaults to the particle path's behaviour.)
5. **Project routines** in `simple_sp_project_ptcl.f90`. (Done on 2026-10-09
   as one routine, `map_cls3D_states2ptcls`: the state-only hand-back that
   respects a supplied active mask and leaves poses alone, tested in
   `simple_test_exec test=unit_project suite=project_records`. The in-place
   `cls3D` update is five setter calls per row and stays in the commander.)
6. **Registration** in the UI and the exec router. (Done on 2026-10-09: a
   developer-visible program in `simple_ui_heterogeneity.f90`, counted in
   `simple_ui_visibility_tester.f90`.)

| File | Change |
| --- | --- |
| `simple_refine3D_states_controller.f90` (new, location to confirm) | controller, march and flex policies |
| `simple_commanders_refine3D.f90` | `exec_refine3D_states` calls the controller; new commander |
| `simple_cavgs_work_project.f90` (new, location to confirm) | class-average work-project builder, even/odd diagnostics |
| `simple_external_reference_pose_initialization.f90` | optional strict-full argument |
| `simple_sp_project.f90`, `simple_sp_project_ptcl.f90` | state-only hand-back `map_cls3D_states2ptcls` |
| `simple_ui_heterogeneity.f90` | program `refine3D_states_cavgs` |
| `simple_exec_refine3D.f90` | router case `refine3D_states_cavgs` |
| `simple_flex_pca_application_tester.f90` | phantom size and state-count arguments (test code; defaults unchanged) |

Dead code noticed on the way. The OpenMP-offload reconstructor
(`simple_reconstructor_openmpoffload.f90`, whose `calc_3Drec_gpu` has no caller
anywhere) is deleted in this work, as a diff of its own, together with its
anchor in `separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`.
In `exec_solve3D_cavgs`, which this work does not touch: both state reads in
the `imgkind` branch, which the `cls2D` states read just after overwrite; the
`noise_norm=no` default, a key only the `normalize` program reads; the
`filt_mode=none` default, which repeats an unconditional setting a few lines
above; and the three branches after the first in the `set_lplims_from_frcs`
dispatch, unreachable because `lpstart` and `lpstop` are always defined by
then. Those are left for a separate cleanup of that file.

Alternatives rejected:

- Entering `exec_refine3D_states` through an in-process command-line key, as
  `solve3D_addon` enters `exec_solve3D`. The class-average path differs in pose
  seeding, flex settings, update scheme, counter handling, the better-row pick,
  ending and commit, which is too much to hide behind a mode branch and a
  key outside the argument vocabulary.
- Aligning and assigning one merged average per sub-class and inserting its
  even and odd averages into the half maps (considered on 2026-10-09). It would
  need a new mode in base reconstruction, both backends, and its halves would
  share every pose. `solve3D_cavgs`'s independent halves are proven on class
  averages, including multi-state runs.

The cost is real: about 600 lines of contained procedures over host variables
move, and the particle path must stay identical. `exec_classify3D_refs`
carries second copies of the coverage and missing-update helpers and can adopt
the controller later.

## Decisions

All design decisions are made. Settled by the rulings: the even/odd handling
of `solve3D_cavgs` (independent rows, the better row wins, disagreement
logged), `solve3D_cavgs` untouched, flex initialization with the settings of
the flex section, flex's own cross-validation with no change to the flex code,
PCG everywhere with `gridding` as an override, `nstates >= 3`, the
class-average mean map, no automasking, the native box, the shared controller
owning flex and the march, strict full-batch operation, the counter reset after
pose seeding, the in-place commit, inactive particles staying inactive, and ML
regularization from the first block. Label balance is no longer a question,
because flex replaces the random labels.

Decided on 2026-10-09:

| Decision | Ruling |
| --- | --- |
| Snapshot low-pass | At the FSC 0.143 point, never finer than `lpstop`, because the FSC reads optimistic. |
| Empty states | Never a failure, as everywhere else in SIMPLE: `refine3D` carries an emptied state forward or drops its out-segment entries (`simple_refine3D_strategy.f90`), and `calc_final_rec` skips a state with no particles. The run reports empty states and completes; at the hand-back the populated states are numbered consecutively, because particle `refine3D_states` requires every state label to be populated (the state-count check of `simple_refine3D_states_controller.f90`). If only one state stays populated, the run completes and says so. |
| Active particles without a committed sub-class | State 0, with the count logged and the original mask kept so a later pass can re-admit them. The kept mask is the `ptcl2D` selection: a rerun hands back over it, and its preflight requires only that every particle active in `ptcl3D` is active in `ptcl2D`. A particle switched off in `ptcl3D` alone takes its sub-average's state again in a rerun; a sub-average left at state 0 stays out on the rerun's labels path unless its `cls3D` state is reset. |
| Frequency band | As `refine3D_states`: 10 to 6 A by default, both overridable. The alignment from the earlier class-average analysis is already good. `lpstart` is also flex's band. |
| Reference low-pass shape | Cosine edge at `lp`, as `refine3D_states`. This applies only with `ml_reg=no`. |

Active particles without a committed sub-class come from three sources: parents
too small for `cls_expansion` to split, deselected sub-classes, and sub-classes
whose rows both fell in a state flex dropped. Class 0 in `ptcl2D` also holds
particles that were already inactive, which is why the selection is
snapshotted first.

Per-particle labels are a seed, not a result. The `cls_expansion` note states
that the latent orders members but does not label individual particles. The
particle-level march must stay free to move them, which it is.

Current behaviour of `cls_expansion` that bears on this program (since
2026-10-12). `ncls` is the total number of sub-classes, budgeted over
class-average groups and classes; a class too small to split is kept whole, so
some rows are whole-class averages (a parent with a single row in `cls2D`). The
final restoration pass reuses the split's soft weights on the refined alignment:
the averages are soft, `neff` is the sum of a sub-class's weights, `sep` is the
split's, and the particle labels stay the split's hard, most likely sub-class.
A sub-class can end with no hard members (`pop` 0, `neff` > 0); it keeps its row
and average. `cls2D` is rebuilt with state 1 on every `cls3D` row, so fresh
output always goes through flex here.

## Validation and documentation

Compilation and runtime tests are performed by the user; nothing here has been
built.

The program's automated test builds on the generated FLEX phantom
(`create_flex_pca_phantom_fixture` and `build_flex_pca_phantom_project` in
`simple_flex_pca_application_tester.f90`), enlarged. Today it holds 2000
particles of two planted states at box 64 with CTF. The fixture gets a size
argument and a state-count argument whose defaults keep today's instance, so
the flex suites and `simple_rec3D_service_tester.f90`, which also uses it, are
unchanged. This test asks for three planted states, because the program needs
`nstates >= 3`, and for enough particles that `cls_expansion` delivers at least
150 selected sub-averages, half again the minimum. The box stays 64, so the
full-box cost is tested only on the manual benchmark. The test runs
`cls_expansion` and this program end to end; it is a long test and runs
nightly, not in the build. The EMD-8440/8445 mix used for the `cls_expansion`
benchmarks stays as a secondary, manual benchmark.

Done on 2026-10-09 (uncompiled). The fixture routines take optional `nptcls`,
`nstates` (2 or 3, with `truth_c`) and `nviews` arguments; without them the
two-state instance is the same, file for file. `nviews > 0` repeats `nviews`
spiral directions, each as one pose, so a direction is an exact 2D class (no
in-plane rotation, no shift) and the project builder writes it as the
`ptcl2D` class; `cls_expansion` then needs no 2D classification before it. The
third state moves the lobe to the third corner of a near-equilateral triangle
with the other two. The test is a stage of the high-level `flex_pca_blobs`
gate, added as checks to that suite rather than as a new CTest entry (the
process budget is an owner decision): 3 600 particles in 60 views (60 per
view, 20 per state), `cls_expansion ncls=3 lp=12` (the states differ in one
lobe of 4 A sigma, which the default 30 A fit band would not see), then this
program with `nstates=3` against the truth mean as consensus map. It gates the
contract (`cls2D` unchanged, `cluster` and `accept` kept, poses kept,
consecutive populated states, the hand-back through the 2D classes, one
`vol_cavg` per state, the consensus entry kept, the work project deleted) and
three detection floors: at least half the sub-averages committed, sub-average
state accuracy of at least 0.60, and a median pair pose error of at most 15
degrees. Particle accuracy, member purity and the map-truth correlations are
reported for calibration. `flex_pca_blobs` is a `highlevel` entry, run by an
explicit CTest command, not nightly.

Regression:

- [ ] `solve3D_cavgs` is untouched: its source and outputs are as before.
- [ ] Particle `refine3D_states` gives the same outputs before and after the
  controller extraction.
- [ ] The three existing callers of the pose initializer are unchanged without
  the strict-full argument.
- [ ] The flex production code is untouched; the flex suites pass on the
  fixture's default instance.

Assertions on the phantom:

- [ ] Flex on the class averages: the paired engine with both fits above 100
  rows, the PCG backend for the basis and the maps, the native box, `lp` equal
  to `lpstart`, the checkpoint map as mean, and `nstates=3` delivering three
  states (no occupancy-floor collapse to two).
- [ ] Full-batch coverage: every active row is updated in every iteration, and
  no fractional or trailing flag is set.
- [ ] Counter reset: `sampled` and `updatecnt` are empty between pose seeding
  and the first march iteration.
- [ ] Better row: each committed sub-average carries the pose, state and score
  of its better-scoring active row, and the even/odd overlap is logged.
- [ ] Metadata preservation: `cluster` and `accept` in `cls3D`, and all of
  `cls2D`, are unchanged by the commit.
- [ ] Empty state: a forced empty state does not fail the run; it gets no
  map, it is reported, and the hand-back numbers the populated states
  consecutively, so particle `refine3D_states` starts on the output.
- [ ] Inactive rows: particles not selected in `ptcl2D` are at state 0 after
  the run, and a rerun on the output starts.
- [ ] Shared-memory and distributed runs give equivalent labels, including
  `nparts` larger than the row count; the flex step runs in shared memory in
  both.
- [ ] State recovery against the planted truth: purity of the committed
  sub-averages and of the handed-back particle labels, and orientation errors
  of the committed sub-averages, with and without ML regularization.
- [ ] Backends: the same recovery measures with `rec_backend=pcg` and with
  `rec_backend=gridding`, as a check on the PCG ruling given the earlier
  failure in `solve3D_cavgs`.

Further checks:

- [ ] The policy refuses `filt_mode=nonuniform_lpset`, `automsk=yes`,
  `pcg_mskfile`, `pcg_solvent=yes`, `update_frac`, `nsample`, `trail_rec=yes`,
  `nstates < 3`, a missing `mskdiam`, missing `_even`/`_odd` stacks, fewer
  than 100 selected sub-averages, a flex fit with fewer than 100 rows and a
  population floor that does not fit, before any write.
- [ ] No NU product appears: the run directory holds no `_nu_filt` volumes and
  no NU envelope masks.
- [ ] End to end: particle `refine3D_states` on the output project starts from
  the handed-back labels without `flex_pca`.
- [ ] The input project's out segment keeps its `vol` and `fsc` entries, and
  its `vol_cavg` entries are replaced as a set or not at all.
- [ ] A failure at any step leaves the input project unchanged; the work
  project and its sigma2 state stay in the run directory.

Documentation to add or update once the design is agreed:

- New `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md`, laid out
  like `solve3D_cavgs_policy.md`.
- `refine3D_states_policy.md`: name `refine3D_states_cavgs` as a second source
  of populated labels, and the controller as the owner of the march and the
  flex initialization.
- New page in `doc/algorithms/heterogeneity_analysis/` and a line in its
  README.
- `doc/how2s/how2_process_heterogeneous_datasets.md`: the route
  `cls_expansion`, then `refine3D_states_cavgs`, then `refine3D_states`.
- The `simple-refine3d` skill code map, only with approval.
