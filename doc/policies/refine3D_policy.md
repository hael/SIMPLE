# Refine3D Policy

## 1. Purpose and Scope

This document defines the current architectural policy for `refine3D` in SIMPLE.

It establishes:

- the architectural boundary between particle-domain work and assembled-reference work
- ownership of command, strategy, matcher, assembly, and postprocessing layers
- iteration and builder-lifecycle rules
- reference-generation, handoff, and assembly policy
- stable artifact contracts between workflow phases
- refactor rules that preserve shared-memory and distributed parity

This document is intentionally broader than an implementation note, but it is not a line-by-line implementation specification. Its purpose is to define the architectural model, the durable workflow contracts, and the review standards that future changes must satisfy.

Where this document describes a current mechanism, file contract, or exception path, that material should be read as an implementation contract in service of the policy, not as a license to dissolve the higher-level architectural boundaries.

---

## 2. Architectural Policy

### 2.1 Core design rule

`refine3D` is a layered fixed-point workflow with a strict architectural boundary:

- particle-domain work stays in probabilistic alignment, search strategy, matcher preparation, orientation update, and partial reconstruction
- assembled-reference work stays in the explicit assembly pathways and their postprocessing policy/helpers

In other words:

- `refine3D` decides and updates particle poses
- assembly pathways build and postprocess state references from those decisions

That split is the primary policy boundary and must be preserved.

### 2.2 Public workflow contract

The public workflow for `refine3D` is:

1. initialize execution mode and iteration state
2. optionally run probabilistic pre-alignment
3. run the particle-update matcher/search pass
4. write partition-local Cartesian partial reconstructions or polar partial sums
5. assemble Cartesian volumes or polar references through the explicit assembly pathway
6. persist updated orientation and volume artifacts for the next iteration or downstream programs

The same scientific workflow applies to both:

- shared-memory `refine3D`
- distributed-master `refine3D`

Only the process-launch mechanism differs between those routes. The workflow, contracts, and ownership boundaries must remain aligned.

### 2.3 Review standard

Changes to `refine3D` should be reviewed against the following questions:

- Does the change preserve the particle-domain versus assembled-reference boundary?
- Does it keep orchestration separate from numerical postprocessing and assembled-reference construction?
- Does it preserve a single clear source of truth for current matching references?
- Does it maintain parity between shared-memory and distributed execution?
- Does it preserve explicit artifact and handoff contracts between workflow phases?

If the answer to any of those questions is no, the change is presumptively a policy change rather than a routine implementation update and should be evaluated as such.

---

## 3. Execution and Ownership Policy

### 3.1 `simple_commanders_refine3D`

`simple_commanders_refine3D.f90` owns:

- the `refine3D` entry point
- top-level defaults
- execution-strategy selection through `create_refine3D_strategy`

This layer should remain thin. It is not the place for low-level search logic or assembled-volume policy.

### 3.2 `simple_refine3D_strategy`

`simple_refine3D_strategy.f90` owns:

- iteration control
- shared-memory versus distributed-master execution policy
- scheduler interaction
- iteration counters and run-finalization bookkeeping
- orchestration of probabilistic pre-alignment, matcher execution, partial-reconstruction writing, and assembly-command dispatch

This layer may thread command-line state and workflow state across steps, but it must not absorb numerical postprocessing logic or assembled-reference implementation detail beyond what is needed to dispatch the correct assembly pathway.

### 3.3 Probabilistic pre-alignment ownership

`simple_commanders_prob.f90`, together with:

- `simple_eul_prob_tab.f90`
- `simple_eul_prob_tab_neigh.f90`

owns the probabilistic pre-alignment phase.

That phase is responsible for:

- sampling particles for update
- generating partition-local probability tables
- aggregating those tables across partitions
- writing a single assignment artifact consumed by the matcher

This phase is particle-domain work, not volume-domain work.

### 3.4 Matcher ownership

`simple_strategy3D_matcher.f90` owns the core particle-update pass through `refine3D_exec`.

That includes:

- reference preparation
- search-strategy dispatch
- candidate evaluation
- orientation, state, in-plane, and shift update
- Euclidean sigma update during search when applicable
- writing Cartesian partial reconstructions or polar partial sums for downstream assembly when instructed by the strategy

This is the execution center of particle-domain refinement.

### 3.5 Assembly ownership

`simple_commanders_rec_distr.f90` exposes two explicit assembly commanders:

- the Cartesian assembly commander, whose `execute` procedure is `exec_cartesian_assembly` (`commander_cartesian_volassemble`)
- the polar assembly commander, whose `execute` procedure is `exec_polar_assembly` (`commander_polar_volassemble`)

`exec_cartesian_assembly` owns Cartesian volume assembly. It reduces partition-local Cartesian reconstruction updates, handles even/odd volumes, performs gridding correction, writes merged state volumes, updates per-particle FSC-derived resolution metadata, and runs shared-memory postprocessing such as automasking and nonuniform filtering.

Cartesian matching still uses projected polar central sections. In `refine3D`,
Cartesian assembly owns the next-iteration reprojection handoff just as polar
assembly does: after the assembled half volumes have been restored, trailed,
low-resolution inserted, and FSC-filtered according to the Cartesian policy, the
assembly commander refreshes `POLAR_REFS_even.bin` and
`POLAR_REFS_odd.bin` for the next matcher or probability-table pass. The next
pass should consume that compatible handoff instead of reprojecting the same
volumes again. Volume reprojection is the bootstrap fallback for a fresh/random
starting volume or for a missing/incompatible handoff. Standalone or terminal
`reconstruct3D` calls assemble only volumes and do not own matcher handoff state.

`exec_polar_assembly` owns polar reference assembly for `polar=yes` and for
the current `polar=obsfield` direct-polar-handoff benchmark. The legacy polar
worker writes partition-local polar partial sums for `polar=yes`. For
`polar=obsfield`, reconstruction uses the normal `reconstruct3D` worker and
`calc_3Drec`, writing the same partition-local Cartesian partial
reconstructions used by `polar=no`; polar assembly consumes those partials
directly. The active path does not use an `fgrid_obsfield` observation-field
accumulator. Polar assembly then:

- calculates polar populations
- reduces polar partial sums or Cartesian reconstruction partials
- dispatches to common-line polar-reference normalization for `polar=yes` only
- dispatches to Cartesian-style volume restoration and direct polar projection for `polar=obsfield`
- writes the updated polar-reference handoff files. `polar=yes` writes the
  merged/even/odd triplet, while the direct-polar-handoff benchmark writes the
  projected even/odd pair; readers synthesize the merged reference when only
  that pair exists.

`polar=obsfield` is not a new reconstruction method or a new scientific model.
It is a performance investigation that keeps the `polar=no` Cartesian
accumulation and restoration semantics while emitting polar matcher references
for the next iteration. Its correctness baseline is therefore `polar=no`; any
divergence must be treated as an implementation defect or as evidence that the
performance experiment is not viable, not as a new method with separate quality
criteria.

`polar=obsfield` is currently a Cartesian-reconstructor accumulation path
followed by Cartesian-style restoration and direct polar projection. It uses the
normal `calc_3Drec` path, including `prep_imgs4rec`, `update_rec`, and the
normal partial reconstruction files, then reduces those partials in
`exec_polar_assembly`. `exec_polar_rec3D_worker` is reserved for legacy
`polar=yes` polar partial sums and must not be used for `polar=obsfield`.
Assembly performs the same core
sampling-density correction, FSC prior handling, and optional trailing blend as
the `polar=no` volume assembly path, but does not run the post-assembly automask
or nonuniform filtering stage. It also does not compute or consume common lines;
no common-line redistribution, weighting, or normalization is part of this
direct-polar-handoff benchmark. The final polar matcher references are generated
with the current assembly-time `nspace` and `pftsz`; when the strategy promotes
`nspace_next` or `pftsz_next`, that promoted grid defines the new reprojection
model for the next iteration.

The Cartesian restoration sequence is implemented once in
`restore_cartesian_state_from_partials` and shared by `polar=no` and the
`polar=obsfield` direct-polar-handoff benchmark. New changes to ML ordering,
trailing reconstruction, FSC-prior use, low-resolution even/odd insertion, or
`lpset` handling must be made in that shared restoration helper rather than
copied into separate assembly branches.

The direct-polar-handoff restoration sequence mirrors `polar=no`: reduce the current
Cartesian half-map reconstruction partials, use previous half volumes as the
trailing/FSC prior when `trail_rec=yes`, apply the Cartesian sampling-density
correction, and apply any fractional trailing blend in volume space. The
reference handoff then follows the same `lp`-set versus gold-standard split as
Cartesian refinement. In gold-standard mode, assembly inserts the merged
low-resolution region into the even/odd half volumes to keep independently
matched references docked, applies trailing to the half volumes, and projects the
separate even/odd references. In `lp`-set mode, assembly does not perform
low-resolution even/odd insertion; if trailing is active, it first trails the
even and odd maps independently, merges the trailed maps into the state volume,
and projects that merged trailed volume once for both even and odd polar
reference slots. That final projection is a performance handoff for the next
matcher iteration, not a separate reconstruction model.

Its benchmark labels must name the operations actually being timed: Cartesian
partial reconstruction accumulation, Cartesian partial reduction,
Cartesian-style restoration, direct polar projection, resolution metadata
update, and `POLAR_REFS*` handoff writing. Labels should not describe these
steps as obsfield insertion, obsfield normalization, or obsfield restoration
unless an actual observation-field accumulator is reintroduced. These timings
exist to answer the performance question; they do not establish a separate
reconstruction policy.

This policy applies to both shared-memory and distributed `refine3D`. In both execution modes, `simple_refine3D_strategy.f90` calls the polar assembly commander for polar modes and the Cartesian assembly commander for non-polar Cartesian volume assembly.

Shared-memory refinement still sets the legacy `force_volassemble` key when `volrec=yes` or when any polar mode is active, then deletes that key after assembly. Distributed refinement does not rely on that key; its workers receive the explicit write-partial-reconstruction decision from the distributed strategy path.

### 3.6 Postprocessing policy ownership

`simple_vol_pproc_policy.f90` owns the postprocessing decision table for Cartesian assembly output.

That includes:

- `automask_action`, which decides whether the state automask is regenerated, reused, or ignored
- `automask_tight`, which preserves the exact `automsk=tight` command-line policy value
- `nu_mask_source`, which decides whether nonuniform filtering uses a freshly generated automask, an existing compatible automask, or the spherical fallback
- centralized compatibility checks for state automasks

Postprocessing-policy decisions belong here, not in `refine3D_strategy`.

---

## 4. Iteration and Builder-Lifecycle Policy

### 4.1 Builder lifetime

The `builder` owns derived execution state, not durable workflow state.

Its contents are valid only for the command-line and parameter policy used to build that instance. This matters because stage policy can change derived state such as:

- reference-grid size
- cropped-box geometry
- PFTC frequency range
- symmetry-expanded projection grids
- masks
- strategy work arrays

### 4.2 Shared-memory rebuild rule

Shared-memory `refine3D` must rebuild the strategy toolbox at each iteration after the iteration command line has the settled stage policy for that iteration.

It must not keep one builder instance alive across iterations and then repair it with ad hoc signature checks.

Any particle-state changes made during initialization, such as random initial orientations or cleanup of sampling counters, must be written to the project before the first per-iteration rebuild so the freshly built toolbox sees the intended state.

Likewise, each shared-memory iteration must persist updated orientations to the active project segment before the next iteration rebuilds the toolbox. `algndoc` output is an iteration artifact, but the rebuild path reads the project segment.

### 4.3 Distributed parity

Distributed `refine3D` naturally follows the same policy because workers and assembly commands are launched as fresh command executions.

Shared-memory mode must mirror that lifecycle explicitly:

- persistent handoff state lives in the project, command-line parameters, and documented artifacts
- per-iteration builder state is disposable

### 4.4 Iteration counters and stage semantics

When rebuilding, `which_iter` and other iteration-local counters may be threaded through the build command line, but the stage-level `startit` must remain available in `params`.

Planning predicates such as final-stage-iteration checks depend on the original stage interval, not on a single-iteration child-command view.

The same rule applies to distributed workers. A worker may execute only one iteration, but:

- `which_iter` is the current iteration
- `startit` remains the stage start

Worker command lines must not collapse `startit` to `which_iter`, because assembly-output planning such as `nspace_next` and `pftsz_next` promotion depends on the stage interval.

### 4.5 Per-iteration stale-input cleanup

Partition-local assembly inputs are per-iteration artifacts.

Before a matcher iteration writes Cartesian `vol_stateNN_partPP_*` partial reconstructions or polar `cavgs_*_part*.bin` / `ctfsqsums_*_part*.bin` partial sums, the strategy must remove stale files from previous iterations.

This is especially important when the assembly output reference space changes through `nspace_next` or `pftsz_next`. Stale part files with the old reference count, or stale Cartesian partials from a previous workflow phase, must never be accepted as valid assembly input.

---

## 5. Reference and Assembly Policy

This section contains both policy rules and current implementation contracts. The policy-level rules state what must remain true architecturally; the implementation-contract subsections record the currently valid operational mechanisms that realize those rules.

### 5.1 Policy: source of truth for matching references

The matcher-side reference contract is handoff based after assembly.

For normal iterative `refine3D`, assembly produces the reference representation
that the next matcher or probability-table pass should consume. For `polar=no`
and for the current `polar=obsfield` benchmark, that representation is the
`POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` handoff emitted by assembly on the
appropriate reference grid.

A complete `vol1..volN` set is a bootstrap source only when no compatible
assembly handoff exists, or when a content-changing fresh start has deliberately
invalidated the existing handoff. In that case the volume source is reprojected
once to seed `POLAR_REFS*`.

This distinction must remain explicit. Assembly handoff files are authoritative
for the next iterative pass; fresh starting volumes are authoritative only for
the bootstrap pass that creates the first compatible handoff.

The matcher must not silently choose between a stale handoff and a fresh
content-changing volume. The producer of the new content must invalidate stale
`POLAR_REFS*` before consumers ask for references.

### 5.2 Policy: current-source rules by mode

For Cartesian matching without a probabilistic pre-step, the matcher consumes a
compatible assembly-emitted `POLAR_REFS*` handoff when it exists. If the handoff
is missing or incompatible, it reprojects the complete current volumes and
materializes the result so later consumers do not repeat the same fallback work.

For `polar=obsfield`, reconstruction workers write Cartesian partial
reconstructions, and polar assembly converts the restored dense model into the
next `POLAR_REFS*` handoff. The matcher consumes that handoff. It does not
re-enter the volume domain except for the bootstrap fallback from an explicit
starting volume.

For legacy `polar=yes`, polar assembly remains the producer of state-major
common-line normalized polar reference files.

`POLAR_REFS.bin`, `POLAR_REFS_even.bin`, and `POLAR_REFS_odd.bin` are the file
contract for consumers that need central-section references across process
boundaries, including matchers and probabilistic table workers.

### 5.3 Policy: centralized reference production

Reference-section production is centralized in `simple_matcher_refvol_utils`.

Consumers that require a file handoff may request materialization through the
centralized utility path, but they do not define independent reprojection policy.

`prob_tab`, `prob_tab_neigh`, and `prep_pftc4align3D` are consumers. They do not
own live reprojection policy. When they need references, they should first use a
compatible assembly handoff and materialize from volumes only as the bootstrap
fallback.

### 5.4 Policy: gold-standard versus lp-set assembly handoff

Gold-standard refinement and `lp`-set refinement have different reference
handoff contracts.

In gold-standard mode, even and odd references are independent matching sources.
Assembly must keep them separate, and when low-resolution even/odd docking is
needed it must insert the merged low-resolution region into the half volumes
before the next matcher handoff. If trailing reconstruction is active, the even
map is blended with the previous even map and the odd map is blended with the
previous odd map. The next matching handoff consumes those separate trailed
halfmaps or their separately projected polar references.

In `lp`-set mode, the matcher consumes a merged reference model. Assembly must
not do low-resolution even/odd docking insertion for this path. If trailing
reconstruction is active, the current even and odd maps are first trailed against
their corresponding previous halfmaps, then the trailed halfmaps are merged into
the state volume. The next matching handoff is derived from that merged trailed
volume. For `polar=obsfield`, this means projecting the merged trailed volume
once and copying the result into both even and odd polar-reference slots.

This split applies to both `polar=no` Cartesian assembly and the
`polar=obsfield` direct-polar-handoff benchmark.

### 5.5 Implementation contract: materialization and distributed cache behavior

If a volume-derived reprojection model is generated as a bootstrap fallback, it
is also materialized to `POLAR_REFS*` as the initial handoff artifact.

In distributed matching, all workers can generate the same bootstrap
volume-derived reference model, so only the first worker materializes the shared
cache to avoid concurrent writes to the same handoff files.

Generated polar reference sections are filled through the interpolation limit of
the consumer they are meant to serve. For iterative handoff, assembly owns that
range and may promote the output grid through `nspace_next` / `pftsz_next`.

Probabilistic modes consume the same handoff. `prob_tab` and `prob_tab_neigh`
are separate worker programs and do not own a live volume-reprojection path. If a
probabilistic pre-step is launched before any compatible handoff exists,
`prob_align` materializes `POLAR_REFS*` from the available starting volumes
before launching the table workers. After assembly has run, the probability
pre-step should reuse the assembly-emitted handoff.

### 5.6 Implementation contract: file-handoff reuse and validity

Callers that need a file handoff use `ensure_polar_refs_on_disk`, which applies one policy:

- compatible existing `POLAR_REFS*` are reused first
- missing or incompatible files fall back to a complete parsed volume set only when all state volumes exist
- fresh content-changing volume starts must remove stale `POLAR_REFS*` before calling the helper

For file-based polar handoffs, `polar_ref_sections_available` accepts either:

- a merged `POLAR_REFS.bin` file
- a complete even/odd pair

The reader accepts all three forms currently produced by assembly:

- merged-only
- even/odd-only
- merged plus even/odd

Trailing reconstruction trails the previous even and odd references independently. If the merged `POLAR_REFS.bin` companion is absent, the merged companion is derived from the already-trailed even/odd references.

Availability is a header contract, not an existence-only check. The stored metadata must remain compatible with the current run:

- the stored number of references must match `nspace * nstates`
- the stored `pftsz` must match the current run
- the stored high-frequency interpolation limit must be readable by the current polar Fourier calculator

When `FRCS_FILE` is absent during the bootstrap pass, matcher preparation creates neutral in-memory FRCs rather than using a silent all-zero FRC model.

### 5.7 Implementation contract: bootstrap exception

There is one explicit bootstrap exception.

If polar reference files are missing or incompatible and the full starting-volume
set is available, `ensure_polar_refs_on_disk` may materialize the initial
`POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` pair from the parsed state volumes.

This bootstrap path uses the live `params`, `builder`, and `cmdline`. It does not create a second temporary builder.

Matcher preparation also calls this helper as a last local fallback before
reading `POLAR_REFS*`, so a freshly generated random/start volume can seed the
first matching pass without relaxing the stale-reference checks.

`set_bp_range3D` is the owner of the 3D matching frequency-range policy. Callers must settle the PFTC range with `set_bp_range3D` before invoking `ensure_polar_refs_on_disk` or any path that materializes polar references from volumes. The requested search high-frequency index must not exceed the cropped interpolation limit before `polarft_calc` is constructed.

If the reference files are missing, no complete `vol1..volN` set is available, and the full starting-volume set is not available, the run must first create those references through reconstruction and assembly. Shared-memory polar initialization errors in that case.

### 5.8 Policy: stage-transition reference-space rule

Stage scheduling may emit `nspace_next` or `pftsz_next` when the next stage will use a different polar reference grid than the current matcher iteration.

`nspace_next` and `pftsz_next` are forward-looking assembly hints. They are not the grid used by the current particle-matching pass.

Reconstruction-only child command lines must delete `nspace_next` and `pftsz_next` because those values belong to the `refine3D`-to-assembly handoff, not to plain volume reconstruction.

The policy question is the assembly output reference grid: which `nspace` and `pftsz` should the reference model emitted for the next iteration use? That policy is owned by assembly and its handoff contract, not by the current matcher pass.

### 5.9 Implementation contract: current stage-transition behavior

On the final planned iteration of a stage:

- Cartesian assembly may reproject the assembled volumes with `nspace_next` and `pftsz_next` so the emitted `POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` files match the next stage
- `polar=obsfield` may do the same for performance testing because assembly directly projects the polar handoff model from state-local Cartesian partial reconstructions

Legacy `polar=yes` keeps its emitted references on the current matching grid. When trailing reconstruction averages across a grid increase, previous state-local projections are remapped to the nearest current projection within the same state.

`polar=obsfield` must not choose a separate matcher-side polar output grid. Workers write state-local Cartesian partial reconstructions, and assembly extracts the emitted polar reference model on the assembly output grid.

### 5.10 Policy: content-changing handoffs must invalidate stale references

Content-changing handoffs must invalidate stale reference artifacts even when those artifacts remain dimension-compatible.

Symmetry search inside `abinitio3D` and `abinitio3D_cavgs` is one explicit example. After the symmetrized map becomes the next volume set, the polar route must inject that map as an explicit `vol1` source for the next `refine3D` stage. Stale companion even/odd files must not be allowed to override the injected symmetric model.

### 5.11 Implementation contract: multi-state layout

Multi-state polar references are stored state-major:

- state 1 occupies local projection slots `1:nspace`
- state 2 occupies `nspace+1:2*nspace`
- and so on

All polar insertion, mirroring, FSC/FRC bookkeeping, and ML regularization must
preserve that state-local reference block structure. For `polar=yes`, this
includes common-line normalization.
For the direct-polar-handoff benchmark, this includes projecting the restored
half volumes with the assembly-time reference grid.

Common-line normalization is intra-state only and applies to `polar=yes`.
The direct-polar-handoff path does not use it. Cross-state intersections are not physical.

### 5.12 Policy: abinitio3D stage low-pass snapshots

`abinitio3D` and `abinitio3D_cavgs` write stage volume snapshots and matching `_lp` diagnostic snapshots through `simple_abinitio_utils.f90`.

The stage planner's `lpinfo(istage)%lp` controls the staged search/reference schedule. It is not the policy cutoff for the saved `_stageNN_lp.mrc` diagnostic volume.

Saved stage `_lp` volumes must be low-pass filtered to the current state FSC resolution when an FSC file exists. The planned stage LP is only a fallback when no valid FSC-derived resolution is available. This keeps the diagnostic stage volume tied to measured even/odd agreement rather than to the intended search schedule.

---

## 6. Particle-Domain and Volume-Domain Boundary

The boundary introduced in Section 2 governs both code ownership and future refactors.

### 6.1 Particle-domain work

Particle-domain work includes:

- particle sampling for update
- probabilistic candidate generation
- state, projection, in-plane, and shift assignment
- reference matching; 3D matching either derives central sections from a current volume source or consumes `POLAR_REFS*` handoff files
- sigma updates during search
- writing partition-local Cartesian partial reconstructions or polar partial sums

The key code owners are:

- `simple_commanders_prob.f90`
- `simple_eul_prob_tab*.f90`
- `simple_strategy3D_matcher.f90`
- `simple_strategy3D_*.f90`
- `simple_matcher_3Drec.f90`

### 6.2 Volume-domain work

Volume-domain work includes:

- summing partial Cartesian reconstructions into state volumes
- reducing partition-local polar partial sums or Cartesian partial reconstructions for direct polar handoff
- even/odd averaging and merged-volume generation
- gridding correction
- automask generation and reuse
- nonuniform filtering
- artifact writing that downstream FSC and reference preparation consume

The key code owners are:

- `simple_commanders_rec_distr.f90`
- `simple_vol_pproc_policy.f90`
- `simple_reconstructor_eo.f90`

Mixing these responsibilities makes the workflow harder to reason about and risks divergence between execution modes.

---

## 7. Artifact and Handoff Policy

The workflow depends on conventional artifacts that act as explicit handoff points between layers.

Stable examples include:

- assignment maps written by probabilistic alignment
- partition-local Cartesian partial reconstructions
- partition-local polar partial sums for `polar=yes`
- `POLAR_REFS.bin` / `POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` central-section files consumed by 3D matcher and search preparation
- state volumes and even/odd volumes
- state-specific automasks such as `automask3D_stateNN.mrc`

These artifacts are not incidental implementation leftovers. They are part of the current contract between workflow phases and should be treated as stable unless a coordinated policy change is made.

---

## 8. Probabilistic Refinement Policy

`refine3D` supports probabilistic search-related modes, but the current workflow is not a monolithic soft-assignment EM implementation.

The current implementation policy is:

- probabilistic pre-alignment may build candidate tables and select assignments before the main matcher pass
- the main matcher consumes that information and applies hard orientation and state updates to the working orientation model
- downstream reconstruction and volume assembly operate on those assigned particle updates

This distinction must be preserved in documentation and code review:

- “probabilistic pre-alignment” is accurate
- “hard-assignment particle update” is accurate
- “fully soft volume-integrated EM” is not an accurate description of the current `refine3D` implementation

Future changes may alter the scientific model, but they should not be described as already implemented before the ownership, workflow, and artifact contracts actually change.

---

## 9. Why Explicit Assembly Pathways Own Shared Reference Work

Keeping assembled-reference work in explicit assembly pathways is good policy because it:

- removes duplicated post-assembly logic from separate `refine3D` routes
- keeps expensive shared-memory volume work in one place
- improves behavioral parity between shared-memory and distributed execution
- keeps `refine3D_strategy` focused on orchestration rather than low-level volume handling
- makes benchmarking more meaningful because the same step owns the same costs

This is the right direction and should not be rolled back without an explicit policy decision.

---

## 10. Current Implementation Status

The current implementation has already moved in the right direction.

In particular:

- `simple_refine3D_strategy.f90` now orchestrates probabilistic pre-alignment, matcher execution, bootstrap polar-reference projection, and assembly-command dispatch
- `simple_strategy3D_matcher.f90` is the core particle-update path
- `simple_commanders_rec_distr.f90` has distinct `commander_cartesian_volassemble` and `commander_polar_volassemble` types
- `simple_vol_pproc_policy.f90` factors automask and nonuniform-filter decisions out of `exec_cartesian_assembly`

That postprocessing policy helper is important because it centralizes the decision table for automask cadence, compatibility, and nonuniform-mask precedence.

This section is descriptive rather than normative. It records the present implementation status against the policy above.

---

## 11. Recommended Boundary Going Forward

The recommended long-term split is:

- `simple_commanders_refine3D.f90`  
  Entry point and defaults only.

- `simple_refine3D_strategy.f90`  
  Iteration control, execution-mode orchestration, scheduler interaction, and workflow threading.

- probability-table and matcher/search modules  
  Particle-domain candidate generation, assignment, and partial reconstruction.

- `simple_vol_pproc_policy.f90`  
  Decide automask cadence, compatibility, nonuniform-mask precedence, and related artifact policy.

- assembly commanders  
  Execute the plan for one iteration and one state set.

- numerical modules  
  Implement automasking, FSC masking, and nonuniform filtering details.

This section is guidance for future code movement within the policy boundaries already established above.

---

## 12. Near-Term Improvements

- add end-to-end regression coverage that compares shared-memory and distributed `refine3D` artifact sets
- add a targeted review checklist or test matrix keyed to the policy questions in Section 2.3
- keep implementation-contract details in this document grouped under explicit labels so future policy changes and mechanism changes are easier to separate

---

## 13. Rules to Preserve During Refactors

- Do not move assembled-volume postprocessing back into `refine3D_strategy`.
- Do not treat the assembly commanders as distributed-only helpers.
- Do not merge probabilistic particle-update logic with volume postprocessing logic.
- Do not introduce a second ambiguous source of truth for current matching references.
- Do not describe the current `refine3D` implementation as if it were a single undifferentiated Bayesian engine.
- Do not reuse shared-memory builder-derived state across iterations when stage policy can change; rebuild the toolbox instead.
- Do preserve parity between shared-memory and distributed workflows.

---

## 14. Workflow Summary

1. `commander_refine3D` selects the execution strategy.
2. `refine3D_strategy` manages iteration state, execution mode, and per-iteration shared-memory builder rebuilds.
3. probabilistic pre-alignment may sample particles and write an assignment map.
4. `refine3D_exec` performs particle-domain search and update, then writes Cartesian partial reconstructions or polar partial sums.
5. assembly commanders assemble state volumes for `polar=no` or state-major polar references for `polar=yes|obsfield`.
6. output artifacts are persisted for the next iteration and for downstream FSC and reference consumers.

That is the current `refine3D` policy and the architectural model future changes should follow.
