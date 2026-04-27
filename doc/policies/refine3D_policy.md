# Refine3D Policy

## Scope

This document defines the current architectural policy for `refine3D` in SIMPLE.

It covers:

- command and strategy ownership
- probabilistic pre-alignment and particle-update flow
- builder lifetime and stage-dependent derived state
- the boundary between particle-domain work and volume-domain work
- the role of explicit assembly pathways in shared-memory and distributed execution

This document is intentionally broader than a single implementation note. It is the top-level policy for how the `refine3D` workflow is divided conceptually and in code.

## Core design rule

`refine3D` is a layered fixed-point workflow with a strict boundary:

- particle-domain work stays in probabilistic alignment, search strategies, matcher preparation, orientation update, and partial reconstruction
- assembled-reference work stays in the explicit assembly pathways and their postprocessing policy/helpers

In other words:

- `refine3D` decides and updates particle poses
- assembly pathways build and postprocess state references from those decisions

That split should be preserved.

## Public policy

The current public workflow for `refine3D` is:

1. initialize execution mode and iteration state
2. optionally run probabilistic pre-alignment
3. run the particle-update matcher/search pass
4. write partition-local Cartesian partial reconstructions or polar partial sums
5. assemble Cartesian volumes or polar references through the explicit assembly pathway
6. persist updated orientation and volume artifacts for the next iteration or downstream programs

The same policy applies to:

- shared-memory `refine3D`
- distributed-master `refine3D`

Only the process-launch mechanism differs between those routes. The scientific workflow and ownership boundaries should remain aligned.

## Execution ownership

### `simple_commanders_refine3D`

`simple_commanders_refine3D.f90` owns:

- the entrypoint for `refine3D`
- top-level defaults
- selection of execution strategy through `create_refine3D_strategy`

It should remain thin. It is not the place for low-level search logic or assembled-volume policy.

### `simple_refine3D_strategy`

`simple_refine3D_strategy.f90` owns:

- iteration control
- shared-memory versus distributed-master execution policy
- scheduler interaction
- iteration counters and run-finalization bookkeeping
- orchestration of pre-alignment, matcher execution, partial-reconstruction writing, and assembly-command dispatch

This layer may thread command-line state and workflow state across steps, but it should not absorb numerical postprocessing logic.

### Builder lifetime and derived state

The `builder` owns derived execution state, not durable workflow state. Its
contents are valid only for the command-line and parameter policy used to build
that instance. This matters because `refine3D` stage policy can change derived
state such as reference-grid size, cropped-box geometry, PFTC frequency range,
symmetry-expanded projection grids, masks, and strategy work arrays.

Shared-memory `refine3D` must therefore rebuild the strategy toolbox at each
iteration after the iteration command line has the settled stage policy for that
iteration. It should not keep one builder instance alive across iterations and
then repair it with ad hoc signature checks. Any particle-state changes made
during initialization, such as random initial orientations or cleanup of
sampling counters, must be written to the project before the first per-iteration
rebuild so the freshly built toolbox sees the intended state.

Distributed `refine3D` naturally follows the same policy because workers and
assembly commands are launched as fresh command executions. Shared-memory mode
should mirror that lifecycle explicitly: persistent handoff state lives in the
project, command-line parameters, and documented artifacts; per-iteration
builder state is disposable.

When rebuilding, `which_iter` and other iteration-local counters may be threaded
through the build command line, but the stage-level `startit` must remain
available in `params`. Planning predicates such as final-stage-iteration checks
depend on the original stage interval, not on a single-iteration child-command
view.

### `simple_commanders_prob` and probability-table modules

`simple_commanders_prob.f90` together with:

- `simple_eul_prob_tab.f90`
- `simple_eul_prob_tab_neigh.f90`

own the probabilistic pre-alignment phase.

That phase is responsible for:

- sampling particles for update
- generating partition-local probability tables
- aggregating those tables across partitions
- writing a single assignment artifact consumed by the matcher

This phase is particle-domain work, not volume-domain work.

### `simple_strategy3D_matcher`

`simple_strategy3D_matcher.f90` owns the core particle-update pass through `refine3D_exec`.

That includes:

- reference preparation
- search-strategy dispatch
- candidate evaluation
- orientation/state/in-plane/shift update
- Euclidean sigma update during search when applicable
- writing partial reconstructions or polar partial sums for downstream assembly when instructed by the strategy

This is the execution center of particle-domain refinement.

### Assembly Commanders

`simple_commanders_rec_distr.f90` currently exposes two explicit assembly
commanders:

- the Cartesian assembly commander, whose `execute` procedure is
  `exec_cartesian_assembly` (`commander_cartesian_volassemble`)
- the polar assembly commander, whose `execute` procedure is
  `exec_polar_assembly` (`commander_polar_volassemble`)

`exec_cartesian_assembly` owns Cartesian volume assembly. It reduces
partition-local Cartesian reconstruction updates, handles even/odd volumes,
does gridding correction, writes merged state volumes, updates per-particle
FSC-derived resolution metadata, and runs shared-memory postprocessing such as
automasking and nonuniform filtering. Cartesian matching still uses projected
polar central sections, so this path always projects the prepared Cartesian
volumes into `POLAR_REFS_even.bin` and `POLAR_REFS_odd.bin` for the next
matcher/probability-table pass. That projection is accounted for in the
assembly benchmark as `polar ref projection`.

`exec_polar_assembly` owns polar reference assembly for `polar=yes`,
`polar=direct`, and `polar=obsfield`. The matcher writes partition-local polar
partial sums. Polar assembly then calculates polar populations, reduces the
partial sums, dispatches to the common-line restore path for `polar=yes` or the
direct restore path for `polar=direct|obsfield`, and writes the updated
`POLAR_REFS.bin`, `POLAR_REFS_even.bin`, and `POLAR_REFS_odd.bin` triplet.
Its benchmark reports this work as `polar reference assembly`.

This policy applies to both shared-memory and distributed `refine3D`. In both
execution modes, `simple_refine3D_strategy.f90` calls the polar assembly
commander for polar modes and the Cartesian assembly commander for non-polar
Cartesian volume assembly. The strategy sets the legacy force-assembly key
when `volrec=yes` or when any polar mode is active, then deletes that key after
assembly.

The matcher-side reference contract is file based. `prob_tab`,
`prob_tab_neigh`, and `refine3D_exec` require polar central sections to be
available through `POLAR_REFS*` before they prepare references. They do not
reproject Cartesian volumes themselves. For Cartesian reconstruction this file
set is generated by `exec_cartesian_assembly`; for polar restoration modes it
is generated by `exec_polar_assembly`. `polar_ref_sections_available` accepts
either a merged `POLAR_REFS.bin` file or a complete even/odd pair. The reader
accepts all three forms currently produced by the assembly code: merged-only,
even/odd-only, or merged plus even/odd. The availability check is a header
contract, not an existence-only check: the stored number of references must
match `nspace * nstates`, the stored `pftsz` must match the current run, and
the stored high-frequency interpolation limit must be readable by the current
polar Fourier calculator. When `FRCS_FILE` is absent during the bootstrap pass,
matcher preparation creates neutral in-memory FRCs rather than using a silent
all-zero FRC model.

Stage scheduling may emit `nspace_next` when the next stage will use a larger
reference grid than the current matcher iteration. `nspace_next` is a
forward-looking assembly hint; it is not the grid used by the current particle
matching pass. On the final planned iteration of a stage, Cartesian assembly
may reproject the assembled volumes with `nspace_next` so the emitted
`POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` files match the next stage.
`polar=obsfield` may do the same because the observations are collected into a
3D structure before reprojection. Legacy `polar=yes` and `polar=direct` keep
their emitted references on the current matching grid; when trailing
reconstruction averages across a grid increase, previous state-local
projections are remapped to the nearest current projection within the same
state. Reconstruction-only child command lines must delete `nspace_next`
because the value belongs to the refine3D-to-assembly handoff, not to plain
volume reconstruction.

There is one explicit bootstrap exception. If polar reference files are missing
but all `vol1..volN` inputs exist, `simple_refine3D_strategy.f90` calls
`write_initial_polar_ref_sections` before probabilistic alignment or matching.
That helper uses the live `params`, `builder`, and `cmdline`, projects the
starting volumes with `read_mask_filter_reproject_refvols`, and writes the
initial `POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` pair. It does not create a
second temporary builder. Because this bootstrap path can run before matcher
setup calls `set_bp_range3D`, reference-volume reprojection must enforce the
same PFTC range contract before constructing `polarft_calc`: the requested
search high-frequency index cannot exceed the cropped interpolation limit.
If the reference files are missing and the full
starting-volume set is not available, the run must first create those references
through reconstruction/assembly; shared-memory polar initialization errors in
that case.

Symmetry search inside `abinitio3D_cavgs` is a content-changing handoff for
polar reference sections. After the symmetrized map becomes the next `vol1`,
existing `POLAR_REFS*` files are stale even when their headers remain
dimension-compatible, so they must be invalidated and rebuilt from the
symmetrized reference model.

Multi-state polar references are stored state-major: state 1 occupies local
projection slots `1:nspace`, state 2 occupies `nspace+1:2*nspace`, and so on.
All polar insertion, mirroring, common-line restoration, direct restoration,
FSC/FRC bookkeeping, and ML regularization preserve that state-local reference
block structure. Common-line restoration is intra-state only; cross-state
common lines are not physical.

## Particle-domain versus volume-domain boundary

The most important architectural boundary in `refine3D` is the distinction between particle-domain and volume-domain work.

### Particle-domain work

Particle-domain work includes:

- particle sampling for update
- probabilistic candidate generation
- state/projection/in-plane/shift assignment
- reference matching; 3D matching consumes `POLAR_REFS*` central-section files
- sigma updates during search
- writing partition-local Cartesian partial reconstructions or polar partial sums

The key code owners are:

- `simple_commanders_prob.f90`
- `simple_eul_prob_tab*.f90`
- `simple_strategy3D_matcher.f90`
- `simple_strategy3D_*.f90`
- `simple_matcher_3Drec.f90`

### Volume-domain work

Volume-domain work includes:

- summing partial Cartesian reconstructions into state volumes
- reducing partition-local polar partial sums into state-major polar references
- even/odd averaging and merged-volume generation
- gridding correction
- automask generation and reuse
- nonuniform filtering
- artifact writing that downstream FSC and reference preparation consume

The key code owners are:

- `simple_commanders_rec_distr.f90`
- `simple_volume_postprocess_policy.f90`
- `simple_reconstructor_eo.f90`

Mixing these responsibilities makes the workflow harder to reason about and risks divergence between execution modes.

## Probabilistic refinement policy

`refine3D` supports probabilistic search-related modes, but the workflow is not a monolithic soft-assignment EM implementation.

The current policy is:

- probabilistic pre-alignment may build candidate tables and select assignments before the main matcher pass
- the main matcher consumes that information and applies hard orientation/state updates to the working orientation model
- downstream reconstruction and volume assembly operate on those assigned particle updates

This distinction should be preserved in documentation and code review:

- "probabilistic pre-alignment" is accurate
- "hard-assignment particle update" is accurate
- "fully soft volume-integrated EM" is not an accurate description of the current `refine3D` implementation

## Artifact policy

The workflow depends on conventional artifacts that act as handoff points between layers.

Examples include:

- assignment maps written by probabilistic alignment
- partition-local Cartesian partial reconstructions
- partition-local polar partial sums for `polar=yes|direct|obsfield`
- `POLAR_REFS.bin` / `POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` central-section files consumed by 3D matcher/search preparation
- state volumes and even/odd volumes
- state-specific automasks such as `automask3D_stateNN.mrc`

These artifacts are not incidental. They are part of the current contract between workflow phases and should be treated as stable unless a coordinated policy change is made.

## Why Explicit Assembly Pathways Own Shared Reference Work

Keeping assembled-reference work in explicit assembly pathways is good policy because it:

- removes duplicated post-assembly logic from separate `refine3D` routes
- keeps the expensive shared-memory volume work in one place
- improves behavioral parity between shared-memory and distributed execution
- keeps `refine3D_strategy` focused on orchestration rather than low-level volume handling
- makes benchmarking more meaningful because the same step owns the same costs

This is the right direction and should not be rolled back.

## Current implementation status

The current implementation has already moved in the right direction.

In particular:

- `simple_refine3D_strategy.f90` now orchestrates probabilistic pre-alignment, matcher execution, bootstrap polar-reference projection, and assembly-command dispatch
- `simple_strategy3D_matcher.f90` is the core particle-update path
- `simple_commanders_rec_distr.f90` has distinct `commander_cartesian_volassemble` and `commander_polar_volassemble` types
- `simple_volume_postprocess_policy.f90` factors automask and nonuniform-filter decisions out of `exec_cartesian_assembly`

That policy helper is important. It is now the owner of the postprocessing decision table:

- `automask_action` says whether the state automask is regenerated, reused, or ignored
- `automask_tight` preserves the exact `automsk=tight` command-line policy value
- `nu_mask_source` says whether nonuniform filtering should use the freshly generated automask, an existing compatible automask, or the spherical fallback
- compatibility checks for state automasks remain centralized in the same module

## Recommended boundary going forward

The recommended split is:

- `simple_commanders_refine3D.f90`
  Entry point and defaults only.
- `simple_refine3D_strategy.f90`
  Iteration control, execution-mode orchestration, scheduler interaction, and workflow threading.
- probability-table and matcher/search modules
  Particle-domain candidate generation, assignment, and partial reconstruction.
- `simple_volume_postprocess_policy.f90`
  Decide automask cadence, compatibility, nonuniform-mask precedence, and related artifact policy.
- assembly commanders
  Execute the plan for one iteration and state set.
- numerical modules
  Implement automasking, FSC masking, and nonuniform filtering details.

## Near-term improvements

- add end-to-end regression coverage that compares shared-memory and distributed `refine3D` artifact sets

## Rules to preserve during refactors

- Do not move assembled-volume postprocessing back into `refine3D_strategy`.
- Do not treat the assembly commanders as distributed-only helpers.
- Do not merge probabilistic particle-update logic with volume postprocessing logic.
- Do not describe the current `refine3D` implementation as if it were a single undifferentiated Bayesian engine.
- Do not reuse shared-memory builder-derived state across iterations when stage policy can change; rebuild the toolbox instead.
- Do preserve parity between shared-memory and distributed workflows.

## Workflow summary

1. `commander_refine3D` selects the execution strategy.
2. `refine3D_strategy` manages iteration state, execution mode, and per-iteration shared-memory builder rebuilds.
3. probabilistic pre-alignment may sample particles and write an assignment map.
4. `refine3D_exec` performs particle-domain search, update, and writes Cartesian partial reconstructions or polar partial sums.
5. Assembly commanders assemble state volumes for `polar=no` or state-major polar references for `polar=yes|direct|obsfield`.
6. output artifacts are persisted for the next iteration and for downstream FSC/reference consumers.

That is the current `refine3D` policy and the architectural model future changes should follow.
