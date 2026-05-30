# abinitio2D stream snapshot collapse report

Date: 2026-05-30

## Summary

We have not proven the root cause yet. The changes made so far fall into two
families:

1. Remove or ignore stale 2D/3D project state when a new abinitio2D run starts.
2. Prevent stale first-iteration search state from forcing all particles toward
   a previous class index.

Those changes address real failure modes in stream-derived projects, but they do
not fully explain the newest observation that the starting references themselves
were all zero. If the references on disk are already zero before the search, the
search/objective layer will only expose the problem: every class comparison is
nearly identical, and the old tie handling can then collapse everything into the
last class encountered.

The most important remaining distinction is:

- If `start2Drefs.mrc` is already all zero immediately after reference
  initialization, the bug is upstream of `shc_inpl` and PFTC reference
  polarization.
- If `start2Drefs.mrc` is nonzero on disk but PFTC sees zero references, the bug
  is in reference preparation, most likely stale `cls2D`/`ptcl2D` state causing
  some refs to be skipped.

## Observed symptom

On stream-derived projects with previous selections/classification state, a
fresh 2D run initializes but the first iteration assigns all particles to the
last class, for example class index 100. One diagnostic score from the first
abinitio2D iteration was:

```text
>>> SCORE [0,1] AVG/SDEV/MIN/MAX: 0.368 0.000 0.368 0.368
```

That score is a strong clue. For the Euclidean objective, comparing particles to
zero references can produce a nearly constant value around `exp(-1) = 0.367879`.
A zero standard deviation across classes means the class-dependent term has
effectively disappeared. That is consistent with zero references, or references
that were never polarized/prepared and therefore behave as zero in the objective.

The previous greedy update condition used `>=`. With identical class scores,
that condition walks through ties and keeps replacing the best class with the
latest equal-scoring class. That makes the collapse to the last class a
downstream effect of the constant-score condition, not by itself the source of
the zero scores.

## Snapshot inspection

The failing project snapshot inspected was `~/snapshot_1.simple`, copied to
`/private/tmp/snapshot_1.simple` for local SIMPLE project readers.

Facts from the snapshot:

- `nptcl2D = 320406`
- `ncls2D = 150`
- `nstk = 2009`
- Active/selected particles: `62075`
- Inactive particles: `258331`
- Active particles with stale 2D search state: `62075`
- Particles with positive class assignment: `62075`
- Particles assigned to the last class in the snapshot: `375`
- Particles with class index greater than `ncls2D`: `0`
- Particles with nonzero `corr`: `62075`
- Particles with nonzero `e3`: `0`
- `stkind` values are in range.
- `indstk` is not present as a real particle field; printed zeros are just the
  default value. `get_stkname_and_ind` can still fall back through `stkind`,
  `fromp`, and `top`.
- All stack rows had `ctf=yes`.
- `cls2D` contained 30 zero-population classes:
  `4 6 11 13 24 27 37 55 59 61 63 66 80 83 84 85 90 94 102 105 110 112 113 116 117 120 121 128 140 148`

Interpretation:

- The snapshot definitely contains stale 2D state in `ptcl2D` and stale
  aggregate `cls2D` state from a previous stream/classification step.
- The provided snapshot does not support the earlier suspicion that bad
  `stkind` reindexing is the immediate cause. The `stkind` values look valid.
- The provided snapshot does not support the CTF-flag suspicion. All stack rows
  are CTF-enabled.
- The project metadata alone is not enough to verify whether the source stack
  images or generated starting-reference MRC files are nonzero.

## Last three commits

### `a5cebb3b` - removing stale 2D artifacts before fresh abinitio2D run

This was the first broad cleanup attempt in
`src/main/commanders/simple/simple_commanders_abinitio2D.f90`.

What it did:

- Added `cleanup_stale_run_artifacts`.
- Deleted `FRCS_FILE` and `ABINITIO2D_FINISHED` before the new run.
- Changed the fresh-run particle cleanup to:
  `delete_2Dclustering(keepshifts=.false., keepcls=.false.)`.
- Explicitly killed `os_cls2D` and `os_cls3D`.
- Replaced `write_segment_inside(params%oritype, params%projfile)` with a full
  `spproj%write(params%projfile)` so killed class segments would be omitted from
  the binary project file.

Why we tried it:

- The first hypothesis was that stale class segments and stale particle
  alignment fields were being left in the project and read by the fresh run.

Why it was not enough:

- The collapse still happened.
- It was too broad: it rewrote the whole project and changed abinitio2D behavior
  outside the narrow first-iteration search problem.
- It was later reverted in `77c272d9`.

### `2d967ac6` - ignoring output fields in CTF-heteroegenous merge projects

This changed merge behavior in `src/fileio/simple_projfile_utils.f90` and the
merge test in `src/main/project/simple_project_merge_tester.f90`.

What it did:

- Treated only `mic`, `stk`, `ptcl2D`, `ptcl3D`, and `optics` as mergeable data
  segments for `merge_selected_project_files`.
- Ignored analysis-product segments `cls2D`, `cls3D`, and `out`.
- Stopped remapping/copying class-output rows from those ignored segments.
- Reset copied `ptcl2D` clustering with
  `merged_proj%os_ptcl2D%delete_2Dclustering(iptcl2D_glob)`.
- Updated tests to expect empty merged `cls2D`, `cls3D`, and `out`, and reset
  `ptcl2D` class assignments.

Why we tried it:

- Stream/selection workflows can merge projects that already carry old 2D class
  products. Carrying `cls2D`, `cls3D`, and `out` forward into a new analysis is
  dangerous because they look like valid analysis state but no longer match the
  intended fresh run.

Why it was not enough:

- The failing snapshot inspected here still has `cls2D`, `cls3D`, and `out`.
  It is a stream snapshot, not necessarily a project produced by the modified
  merge path.
- The snapshot's `stkind` mapping looks valid, so the collapse is unlikely to be
  explained by merge `stkind` reindexing in this specific case.

### `77c272d9` - trying to get 2D going from stream projects

This reverted the broad abinitio2D cleanup and moved the fix into 2D search
preparation.

What it did:

- Reverted the full-project rewrite and class-segment killing in
  `simple_commanders_abinitio2D.f90`.
- Returned the abinitio2D fresh-run cleanup to the narrower:
  `spproj_field%delete_2Dclustering`
  followed by `write_segment_inside(params%oritype, params%projfile)`.
- In `src/main/strategies/search/simple_strategy2D_alloc.f90`, introduced a
  fresh-start condition for `startit <= 1`, first iteration, not `continue=yes`,
  and not fill-in mode.
- On fresh start, `prep_strategy2D_glob` now allocates all classes as eligible
  instead of trusting stale `cls2D` populations or stale `ptcl2D` class fields.
- On fresh start, `prep_strategy2D_batch` now ignores the old particle class
  when generating the stochastic search order, so `put_last` does not move a
  stale previous class to the end.
- In `src/main/strategies/search/simple_strategy2D_srch.f90`, fresh start now
  ignores stale `os%get(class)` when seeding `prev_class`.

Why we tried it:

- The user observation was that search-order generation puts the current class
  index last to prevent cycling, and every particle appeared to have the last
  class index assigned. Ignoring stale previous class on the first fresh
  iteration directly addresses that mechanism.

Why it was not enough:

- It protects the search order and initial `prev_class`, but does not guarantee
  that references on disk are nonzero.
- It also did not initially protect PFTC reference preparation from stale
  project state.

## Current uncommitted changes

The current working tree adds a small helper and applies the same fresh-start
definition consistently in search and PFTC prep.

### Shared fresh-start helper

File: `src/main/strategies/search/simple_strategy2D_alloc.f90`

Added public helper:

```fortran
logical function is_fresh_2D_start( params, which_iter )
    class(parameters), intent(in) :: params
    integer,           intent(in) :: which_iter
    is_fresh_2D_start = params%startit <= 1 .and. which_iter <= params%startit &
        &.and. trim(params%continue) /= 'yes' .and. .not. params%l_fillin
end function is_fresh_2D_start
```

This removes duplicated fresh-start logic from the search modules.

### Search allocation and previous-class setup

Files:

- `src/main/strategies/search/simple_strategy2D_alloc.f90`
- `src/main/strategies/search/simple_strategy2D_srch.f90`

Current behavior:

- Fresh start gives every class a nonzero synthetic population for search
  eligibility.
- Fresh start does not use stale `ptcl2D%class` when constructing per-particle
  search order.
- Fresh start does not seed `prev_class` from stale `ptcl2D%class`.

This addresses the user-observed mechanism where stale current class could be
placed last in the search order for every particle.

### PFTC reference preparation

File: `src/main/strategies/search/simple_matcher_pftc_prep.f90`

Current behavior:

- Imports `is_fresh_2D_start`.
- Computes:

```fortran
l_fresh_start     = is_fresh_2D_start(params, which_iter)
has_been_searched = (.not. l_fresh_start) .and. (.not.build%spproj%is_virgin_field(params%oritype))
```

Why this matters:

- Previously, stale `ptcl2D` search fields could make the project look
  non-virgin on the first fresh iteration.
- If PFTC reference prep believes the run has already been searched, it can use
  existing class populations to decide which references to prepare.
- In the inspected snapshot, `cls2D` has zero-population classes. Those could be
  skipped even though the fresh run has just generated a full reference stack.
- A skipped reference can look like a zero reference to the objective function.

This is downstream of starting-reference generation. It explains zero
references inside the search/PFTC layer, not an already-zero `start2Drefs.mrc`
file on disk.

### Greedy tie handling

File: `src/main/strategies/search/simple_strategy2D_greedy.f90`

Changed:

```fortran
if( inpl_corr >= corr )then
```

to:

```fortran
if( inpl_corr > corr )then
```

Why this matters:

- With all class scores equal, `>=` makes the last equal score win.
- With strict `>`, ties no longer deterministically walk to the last class.

This is a guard against collapse under degenerate scores. It does not fix the
underlying zero-reference or constant-score cause.

## What we reverted or rejected

### Broad abinitio2D project rewrite

The earlier change that killed `cls2D`/`cls3D` and rewrote the entire project was
reverted because it did not fix the bug and changed too much state at the
commander level.

### CTF/sigma suspicion

We considered whether CTF or sigma handling could be causing the constant score.
The inspected snapshot argues against this as the primary cause:

- All stack rows had `ctf=yes`.
- The score pattern is more directly explained by zero references or skipped
  reference preparation.

No CTF code change remains in the current diff.

### Shift carryover

We considered stale shifts early, but shifts alone should not make every class
score exactly identical with zero standard deviation. The `0.368` constant-score
observation points more strongly at zero references or a degenerate objective
input.

## What explains starting references being all zero?

There are two different "zero reference" cases, and they imply different bugs.

### Case 1: the MRC starting-reference file is already zero

If `start2Drefs.mrc` is all zero immediately after `init_cluster2D_refs` or
abinitio2D `inirefs`, the current search/PFTC fixes cannot explain it. The bug
is upstream in reference initialization.

Relevant code paths:

- `abinitio2D` defaults `cls_init=rand` if no value is supplied.
- For a normal fresh run with no input `refs`, `cluster2D` should call
  `init_cluster2D_refs`.
- `init_standard_refs` with `cls_init=rand` calls:
  `noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)`.
- `noise_imgfile` calls `img%ran` for every image before writing, so the output
  should not be all zero.

Therefore, if `cls_init=rand` and no input `refs` are defined, all-zero
`start2Drefs.mrc` means one of these is true:

- `noise_imgfile` was not called.
- The reference file was overwritten after `noise_imgfile`.
- The run is not actually using the default no-refs path.
- The file being inspected is stale from an earlier failed run, not the one just
  generated.

There is also one suspicious direct path:

- `abinitio2D%inirefs` normalizes user-supplied references into
  `start2Drefs.mrc`.
- If the incoming `refs` argument already equals `start2Drefs.mrc`, then
  `copy_imgfile(refs, params%refs, ...)` becomes a copy from a file to itself.
- `copy_imgfile` opens the source for read and then opens the destination for
  write. It does not guard against source and destination being the same path.
- A self-copy can plausibly truncate or corrupt the file before reads complete,
  producing zero or invalid references.

That self-copy path is a plausible explanation if the stream/NICE command line
or a resumed failed job is carrying `refs=start2Drefs.mrc` into abinitio2D.
It is not yet proven because we have not captured the stage-1 command line.

### Case 2: the MRC file is nonzero, but PFTC/search sees zero references

This is what the current uncommitted `prep_pftc4align2D` change targets.

In the inspected snapshot, stale `ptcl2D` state makes the project non-virgin and
stale `cls2D` has zero-population classes. If PFTC reference preparation trusts
that stale state during a fresh start, it can skip classes that are valid for
the new run. Those skipped classes then behave as zero references in the
objective.

This case explains:

- valid-looking starting reference files,
- zero or unprepared PFTC references,
- flat `0.368` Euclidean scores,
- and collapse to the last class through old tie behavior.

It does not explain an already-zero `start2Drefs.mrc` file.

## Why clean runs behave differently

A clean project that has never had 2D classification typically lacks stale
`ptcl2D` class/corr fields and stale `cls2D` populations. That makes these
checks naturally fall into the "first iteration, all classes are valid" path:

- `is_virgin_field('ptcl2D')` is true.
- `ptcl2D%class` is absent or zero.
- `cls2D` is absent or empty.
- Search order has no meaningful previous class to place last.
- PFTC prep has no stale populations to use for skipping refs.

The stream snapshot is different:

- Active particles already have previous class/corr state.
- `cls2D` exists and includes zero-population classes.
- `out` and `cls3D` also exist.
- The project can look like a continuation even when the user intends a fresh
  abinitio2D/cluster2D run.

That explains why stale-field guards are necessary. It still does not fully
explain zero starting-reference files.

## Verification so far

Completed:

- Inspected the provided snapshot metadata with a temporary Fortran reader.
- Confirmed `stkind` ranges look valid.
- Confirmed all stack rows in the snapshot have `ctf=yes`.
- Confirmed stale active-particle 2D search/classification state exists.
- Confirmed stale `cls2D` zero-population classes exist.
- Built the current code with:
  `cmake --build build --target SIMPLE3.0.0`
- Ran `git diff --check`.

Not completed:

- We have not reproduced the full failing run from the snapshot because the
  project file alone does not include all image-stack data needed to prove the
  reference images are nonzero.
- We have not yet instrumented the exact reference files at the moment they are
  created.
- We have not captured the actual stage-1 `cluster2D` command line to verify
  whether `refs` is defined and whether it points to `start2Drefs.mrc`.

## Recommended next diagnostic

Add temporary logging or a debug helper at two boundaries.

First boundary: immediately after abinitio2D/cluster2D reference initialization.
Log:

- `params%cls_init`
- whether the incoming command line defines `refs`
- `params%refs`, `params%refs_even`, `params%refs_odd`
- for each reference stack: number of images, dimensions, min, max, mean,
  standard deviation, and sum of squares per image

Second boundary: immediately after PFTC reference polarization in
`prep_pftc4align2D`. Log:

- whether `is_fresh_2D_start` is true
- whether `has_been_searched` is true
- class populations used for reference prep
- count/list of references with zero polar sum of squares

This should separate the two failure modes cleanly:

- Disk refs are zero: fix reference initialization or stale `refs` command-line
  propagation.
- Disk refs are nonzero but PFTC refs are zero: fix PFTC prep/skipping logic.

The most suspicious new lead from the zero-starting-reference observation is the
self-copy possibility in `copy_imgfile` when an inherited `refs` argument already
points at `start2Drefs.mrc`.
