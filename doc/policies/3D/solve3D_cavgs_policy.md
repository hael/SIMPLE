# Solve3D Cavgs Policy

This document records the current policy for `solve3D_cavgs`, the
class-average route for de novo map determination (ab initio 3D
reconstruction coupled with initial 3D refinement). It is separate from
the particle-based [solve3D_policy.md](solve3D_policy.md) policy and
from the base [refine3D_policy.md](refine3D_policy.md) policy.

## 1. Scope

`solve3D_cavgs` determines a 3D map de novo from selected 2D or 3D class
averages. It creates a temporary project containing even and odd class-average
entries, runs staged `refine3D` over that temporary project, maps the selected
class orientations back to the original project, and optionally writes
validation artifacts.

Only refinement and reconstruction children are distributed. The
`solve3D_cavgs` master itself rejects direct worker execution with `part`.

## 2. Defaults

The route sets:

- `sigma_est=global`
- initial `oritype=out`, then `ptcl3D` for the temporary project
- `bfac=0`
- `filt_mode=none`
- `automsk=no`

Canonical sigma persistence is used for single-state and multi-state runs.

When unset, it supplies:

- `mkdir=yes`
- `objfun=euclid`
- `overlap=0.95`
- `prob_athres=90`
- `cenlp` from the solve3D controller default
- `imgkind=cavg`
- `noise_norm=no`
- `lpstart=20`
- `lpstop=8`
- `gauref=yes`

NU filtering and automasking are intentionally disabled for this route.

## 3. Temporary Project

The workflow writes a temporary project named
`solve3D_cavgs_tmpproj.simple`.

The temporary project:

- preserves project info, compute environment, and job-process metadata
- registers the even and odd class-average stacks as two stacks
- expands each class average into an even entry and an odd entry
- copies class/state/orientation metadata into both entries
- sets even/odd flags and stack indices for the temporary `ptcl3D` segment

The temporary project never inherits the input project's canonical sigma path.
It registers a workflow-local transient canonical state
file for the temporary class-average particle lineage. That file is rebuilt at
startup and removed with the temporary project after successful completion.

The staged matcher and subsequent standalone reconstruction children resolve
that registered state through the shared sigma-group loader. Canonical loads
validate the committed state against the temporary project's native grid,
ordered particle layout, and grouping policy before reconstruction begins.

The temporary project is deleted at the end of the workflow.

## 4. Inputs and Low-Pass Schedule

`imgkind=cavg` reads state labels from `cls2D`. `imgkind=cavg3D` is accepted
as an input mode, but the current implementation still retrieves the selected
class-average stack and then uses class state information for the temporary
entries.

The class-average stack plus `_even` and `_odd` companion stacks must exist.
If all class-average states are zero, the workflow stops.

The number of particles in the temporary project is twice the number of class
averages. If distributed execution requested more partitions than temporary
entries, `nparts` is reduced to the number of even/odd class-average entries.

Stage low-pass limits are derived from class FRCs by default. Explicit
`lpstart_ini3D` and `lpstop_ini3D` override that schedule together; supplying
only one is an error.

## 5. Staged Refinement

The number of ini3D stages is capped by `solve3D_nstages_ini3D_max()`. A user
`nstages` value can shorten the route up to that cap.

`nstates` alone selects the mode, as in particle `solve3D`: `nstates=1` is
the single-state run and `nstates > 1` refines independent states from the
start. `solve3D_cavgs` has no `multivol_mode` input and refuses it on the
command line; the docked mode (one model up to a split stage, then a random
split of the class averages into states) was removed on 2026-10-06.

Before staged refinement, `rndstart` randomizes orientations, zeros shifts,
randomizes states with balanced uniform labels for multi-state runs, and
reconstructs starting volumes. Starting volumes and half maps are renamed to
the standard `refine3D` start-volume names, including `_unfil` copies of the
half maps.

Each stage is configured through the shared solve3D stage controller with
`l_cavgs=.true.`. In this mode:

- `envfsc=no`
- `filt_mode=none`
- `automsk=no`
- `update_frac` is deleted from emitted refine3D child commands
- `snr_noise_reg` is set from the stage controls
- early Gaussian reference filtering remains available through stage policy

For multi-state runs, class-average stages 1 and 2 retain
`refine=prob_neigh` with `prob_neigh_mode=shc`. They do not adopt the direct
`refine=shc` startup used by particle-based multi-state `solve3D`.

Because `l_cavgs=.true.` removes `update_frac`, `refine3D` derives full-update
mode and disables trailing reconstruction effectively, even when the shared
stage controller emits `trail_rec=yes`.

At the symmetry-search stage, the workflow runs the shared symmetry handling
used by the solve3D workflows.

### 5.1 State reseeding (`reseed_states`, default `no`)

In a multi-state run, a state can empty out and never come back. Each
iteration gives every class-average half the state it fits best, so a state
with more members reconstructs better and draws more. The probabilistic search
also drops a state of 5 halves or fewer (`eul_prob_tab`), and an empty state
takes no members back. `reseed_states=yes` (added 9 October 2026) lets a state
that has emptied start again, at each stage boundary before the last stage
(`reseed_state_labels`, `simple_solve3D_utils`):

- **A state is weak** when it holds no more than 5 halves or less than 2% of
  the selected ones. Both constants are provisional (`RESEED_MIN_POP`,
  `RESEED_MIN_FRAC`).
- **Each weak state, the emptiest first, takes members of the most populated
  state.** These are its worst-fitting classes, ranked by the mean score of
  their halves in that state. Both halves of a class move together. It takes
  them until it holds an equal share of the selected halves, or until half the
  donor's halves have moved, whichever comes first. Moved halves keep their
  poses.
- **Nothing is relabelled** when a reseeded state would hold 5 halves or fewer,
  or a state would stay empty: there are too few classes.
- **The next stage starts from rebuilt volumes.** They are reconstructed from
  the new labels by the stage-boundary reconstruction (`calc_rec`), so that
  stage's `_stageNN` volumes are the reseeded ones. The relabelling and the
  populations are logged.
- **The collapse check comes after the reseeding.** A reseeded run is not
  exited early (`exit_collapse`), and the restart driver (`nrestarts_collapse`)
  sees the states it ends with. A state that empties in the last stage is not
  reseeded.

The worst-fitting classes are also the noisiest, so a reseeded state can
start from weak signal. A small but genuine state that falls under 2% after a
stage is reseeded too. Validate the option before making it a default.

## 6. Mapping Back

After staged refinement, the temporary `ptcl3D` and `out` segments are read
back. Temporary stack-index and even/odd fields are removed.

For each original class, the even or odd temporary entry with the better
correlation is selected. Its correlation, projection index, Euler angles,
2D shift, and state are copied into the original project's `cls3D` segment.

The resulting `cls3D` orientations are mapped back to particles with
`map2ptcls`.

## 7. Validation and Outputs

When the route runs the maximum ini3D stage count, it produces validation
artifacts:

- final original-sampling reconstructions without automatic postprocess
- final raw and low-pass diagnostic volumes
- `vol_cavg` entries in `os_out` for final class-average volumes
- `final_oris.txt`
- reprojections of the final low-pass volumes
- an alternating `cavgs_reprojs.mrc` stack
- a shifted class-average stack registered as `cavg_shifted`
- optional class-average ranking by cavg-vs-reprojection correlation

Even/odd convergence diagnostics report average angular distance. Multi-state
runs also report even/odd state overlap.
