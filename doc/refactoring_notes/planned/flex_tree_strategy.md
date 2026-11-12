# FLEX tree strategy: contract and incremental migration

## Current increment

The tree strategy now executes the sequential **node FLEX -> multistate refine3D
-> split** loop. These stages call the existing numerical programs through
caller-supplied command objects. The strategy owns the tree and orchestration;
`simple_flex_tree_project_gateway` owns stage-project preparation, model artifact
references, and state-derived child payloads. No numerical implementation was copied.

There is no `mode=tree` CLI dispatch yet. Final project/report/checkpoint publication
is still explicitly unimplemented: `finalize_run()` raises an error. `get_tree()`
returns a copy of the current hierarchy for the caller to inspect. Existing FLEX
execution and all legacy tree code remain unchanged. The legacy reference is
`simple_hetero_tree.f90` at hetero-tree-wip commit
`bda4b66c03f6bbb3107b04c29ae91036237fe9f4`.

## Lifecycle and ownership

The future commander normalizes defaults, then initializes parameters once through
`build%init_params_and_build_spproj`. The tree strategy borrows `params`, `build`
and `cline`; it only owns its tree and work-directory reference. Input validation
requires matching, nonempty 2D/3D particle fields, positive box/sampling metadata,
`oritype=ptcl3D`, and an active 3D selection. The root stores an absolute input
project reference and optional consensus `vol1`; construction does not write the
project or allocate the numerical toolbox.

The caller supplies concrete existing commanders through the abstract
`commander_base` interface. This avoids importing a concrete commander back into
a strategy that the future FLEX commander will itself call:

```fortran
! Caller declarations (plus its existing params, build and cline):
type(commander_flex_pca) :: flex_runner
type(commander_refine3D) :: refine_runner
type(flex_tree_strategy) :: strategy
type(flex_tree) :: tree

call strategy%initialize(params, build, cline)
call strategy%execute(params, cline, flex_runner, refine_runner)
tree = strategy%get_tree()
! Publication is a later increment; do not call finalize_run() yet.
call strategy%cleanup()
call build%spproj%kill
nullify(build%spproj_field)
```

The project-only builder must be released explicitly; `kill_general_tbox` is
guarded by `general_tbox_exists` and does not release this lifecycle on its own.
Each child command owns and releases its own parameters, builder and numerical
strategy. After each normal command return, the tree strategy restores the
working directory and OpenMP/thread-count context. Fatal stage errors retain the
stage files and terminate through the existing SIMPLE error handling.

## Stage and membership contract

Execution requires a fresh `flex_tree/` directory, with
`node_000001/flex/flex.simple` and `node_000001/refine/refine.simple`, and similarly
for later stable node IDs. Existing run/stage directories are refused; resume is
not implemented. Stages run sequentially; the existing commands can use distributed
workers when `nparts > 1`. There is no independent-node scheduler yet.

Each stage project retains every row of the original input in the original order.
Only its selected population remains active in both 2D/3D fields; `pind=i` records
the row in the tree's source project. No separate tree membership file or stack
compaction is introduced. This costs a full project copy per stage but makes the
initial identity contract explicit. Original source files are not overwritten.
Project names and selections change, so inherited output registrations and the
canonical sigma2 registration are removed; the existing commands rebuild their
own stage outputs and noise state.

1. **Node FLEX:** select the node's state(s), reset its active 3D labels to one
   population, and call the normal FLEX command with its mean map. The commander
   supplies its normal fit defaults; explicitly supplied FLEX controls are
   forwarded, and `rec_states=yes` requests candidate maps. The frozen fitting
   project, embedding (latents/eigenvalues/precisions), mean, registered weight
   manifest, state-target table and run manifest are retained in `node%models`.
2. **Refine:** count populated FLEX states and copy their maps/poses to a separate
   refinement project. Relabel occupied states contiguously, preserving map order,
   then pass exactly that count as `nstates` and its maps as `vol1...volK` to
   `refine3D`. The initial policy is `refine=prob_neigh`,
   `prob_neigh_mode=state`, `multivol_mode=input_oris_refine`, no centering or
   trailing reconstruction, and ordinary hard-label reconstruction. It runs up to
   10 iterations unless the supplied `maxits` overrides that cap. Refinement does
   not overwrite the FLEX fitting project or its artifacts.
3. **Split:** count populated positive states in the completed refinement project.
   Zero gives `TREE_JUNK`; one gives `TREE_LEAF`; K >= 2 gives `TREE_SPLIT` and
   exactly K children. Each child stores that frozen project, one state selector,
   its registered map/FSC, and the existing paired half maps when present.
   Zero/unassigned rows remain in the frozen project, outside the children.

The indexed-tree mutation happens only after preparing all child payloads. Parent
decisions are assigned after adding children because `add_child` invalidates an
old parent decision. Children enter the same breadth-first loop. A node with fewer
than 100 selected rows in either mod-4 fit stops as a leaf, retaining its inherited
map and membership; the existing paired engine requires that minimum. An occupied
state is a workflow branch here, not a new claim of independently validated
conformational novelty. Scientific split/merge assessment remains a later layer.

## Migration checklist

Keep this table current as each operation is implemented and verified. Nothing
is silently retired; record a deliberate replacement or removal here. Existing
source is the reference, not an automatic requirement to preserve every defect.

| Responsibility / legacy entry points | Destination | Status |
| --- | --- | --- |
| Parameter/project initialization and root construction | Tree strategy lifecycle | Implemented; static checks passed, runtime pending |
| `node_project`, `compact_proj`, `expand_proj`, row/state mapping | Project gateway; canonical `pind` and frozen fit-project provenance | Full-row masked stage copies implemented; compaction not introduced |
| `start`, `advance`, mean reconstruction and initial/focused FLEX | Strategy stage orchestration using existing single-run FLEX/reconstruction strategies | Existing FLEX and multistate refine3D commands invoked; focused/site variants pending |
| `oneflex_inputs`, parent-model reuse and lineage settings | Typed stage inputs and model provenance | Pending |
| `cheap_children`, `falloff_mask`, `half_variance`, site/island discovery | FLEX state services and existing image/volume kernels | Pending |
| `decide`, `state_split`, `state_mergeup`, children/leaves/junk | Scientific assessment services; controller applies accepted results | Populated-state splitting and leaf/junk handling implemented; scientific assessment and merges pending |
| Site batches, `fs_sync`, ranking, refinement, guards, factorial outputs | Site orchestration and numerical site services | Pending |
| `submit*`, `step_status`, `watch`, `core_room`, `busy_update` | Existing qsys plus strategy job/dependency/resource accounting | Pending |
| `read_table`, `write_table`, site tables and restart/recovery | Versioned scientific checkpoint and separate execution ledger | Pending |
| `submit_subtree`, `merge_subtree` | Explicit delegation policy under the same resource/identity contract | Pending |
| `finish`, coverage/duplicate checks, `final_project`, final state hierarchy | Delivery/gateway and state-similarity calculation using the unified tree type | Pending |
| `final_maps`, tree Wiener filtering, landscape/reports | Existing volume operations and delivery services | Pending |
| `prune_step`, `prune_node`, `prune_maps` | Retention based on referenced artifacts and unfinished readers | Pending |
| UI/typed mode parameter and commander branch | Existing FLEX entry chain; no selector module | Pending |

## Validation criteria and review

Use language-server diagnostics, descriptor/source-unit and FLEX DAG checks,
whitespace/line-length checks, and source review of stage interfaces and ownership.
Compilation and runtime tests remain for the user unless explicitly requested.

Source review covers populated/gapped state labels, 0/1/K outcomes, row-preserving
selection, maps matching compacted labels, isolated fitting/refinement projects,
ordinary versus weighted refinement, one-part shared-memory execution, and parent
decision ordering. All these checks are static; they do not establish runtime or
numerical correctness.

Before CLI integration, run native single-node and multi-node cases, small-child
termination, exhausted states, single/distributed parity, project identity and
coverage checks, and failure propagation. Later increments still need scientific
split/merge gates, complete fitted-basis/noise metadata registration, final output
publication, transactional restart, resource accounting and independent-node
scheduling. Existing stage directories must be retained until those delivery and
retention contracts are implemented.

The current 391-line strategy and 236-line gateway pass fortls 3.2.2 diagnostics
(no errors or warnings), FLEX dependency checks, descriptor checks and source-style
checks. Navigation was regenerated in scratch. At the user's request, both modules
and their 150-module dependency closure compiled successfully with GNU Fortran
14.3.0, with no diagnostics from either new module.

A fresh copy of master at fa9b82b5723e66883f36e8a891837526090b5516, with only these
four review files overlaid, also completed the CMake Release build and production
executable links. The build used GNU Fortran 14.3.0, Apple Clang 21 and the system
macOS linker; the bundled linker could not parse the installed SDK. That selection
was confined to the temporary build configuration. Existing build/link warnings
remain; no repository build configuration or source fix was needed. The final
configuration had BUILD_TESTS=OFF and FORTLS=OFF. The build's complete UI JSON
validation passed. No workflow runtime tests or scientific jobs were run.
