# Abinitio3D Cavgs Reject Policy

This document records the current policy for `abinitio3D_cavgs_reject`, the
class-average rejection workflow based on restarted multi-state
`abinitio3D_cavgs` runs. It should be read alongside
[abinitio3D_cavgs_policy.md](abinitio3D_cavgs_policy.md) and
[../microchunk_and_rejection/model_cavgs_rejection.md](../microchunk_and_rejection/model_cavgs_rejection.md).

The implementation lives in `exec_abinitio3D_cavgs_reject` in
`src/main/commanders/simple/simple_commanders_abinitio.f90`. The executable
entry point is registered in `src/main/exec/simple_exec_abinitio3D.f90`, and
the UI definition is in `src/main/ui/simple/simple_ui_abinitio3D.f90`.

## 1. Scope

`abinitio3D_cavgs_reject` selects one good consensus class-average state and
rejects all other class averages. It does this by:

1. evaluating class-average quality with the shared cavg-quality model backend;
2. launching `nrestarts` independent short `abinitio3D_cavgs` runs;
3. reading the restart state labels after stage 2;
4. mapping randomized state labels into one common label space;
5. voting a consensus state label for each class average;
6. choosing the quality-best populated consensus state as good;
7. selecting the corresponding best-state volume from every restart;
8. docking those restart volumes to the best-scoring selected restart volume;
9. averaging the docked volumes into a consensus volume;
10. mapping the resulting binary selection back into the project.

The restart runs may use two or three ab initio states, but the final project
selection is always binary:

- `state=1` means selected or good;
- `state=0` means rejected or bad.

Internally, `final_state=2` is used in the consensus report for active
consensus states that are not selected.

The workflow also writes a docked consensus 3D volume:

```text
abinitio3D_cavgs_reject_consensus_vol.mrc
```

## 2. Public Inputs and Defaults

The route is master-only. Supplying `part` is an error.

The route accepts `nstates=2` or `nstates=3`; values outside that range are an
error. When unset, `nstates=2`.

The route always runs the restart children through the first two
`abinitio3D_cavgs` stages. If the user supplies `nstages` with any value other
than `2`, the command stops. The command then sets `nstages=2` before parsing
parameters.

Point-group symmetry is not a public input to this route. The command sets
`pgrp=c1` and deletes `pgrp_start` before parameter parsing because only C1 is
used for the first two `abinitio3D_cavgs` stages.

When unset, the route supplies:

- `nrestarts=3`
- `mkdir=yes`
- `quality_model=chunk_default_v2`
- `prune=no`

`nrestarts` must be at least 1. In the UI, `nrestarts` is grouped with
`nthr` under compute controls because it controls the number of independent
restart jobs, analogous to the number of parts used by other SIMPLE workflows.

## 3. Quality Evaluation

The command reads the input project and uses the number of `cls2D` rows as the
class-average count. A project with no `cls2D` entries is an error.

The original `cls2D` `state` array is saved before any restart labels are
processed. Original states less than or equal to zero are treated as inactive
class averages during restart-label sanitization and consensus voting.

Class-average images are read through the shared class-average stack reader.
The stack size must match the `cls2D` count. The command then calls
`evaluate_cavg_quality` with the selected `quality_model` and `mskdiam`.

The quality model is initialized from the built-in preset named by
`quality_model`. If `infile` is supplied, the model file is read after preset
initialization and overrides the built-in model specification.

The quality backend produces:

- scalar `quality_scores`, where higher is better;
- automatic model states in `quality_auto_states`;
- quality-cluster labels in `quality%labels`;
- hard-reject and feature diagnostics used later by the feature table.

The automatic model states are diagnostic only in this workflow. The final
accept/reject decision is made from consensus state labels plus mean quality
score per consensus state.

## 4. Restart Execution

Each restart runs in its own folder:

```text
abinitio3D_cavgs_reject_restart_001
abinitio3D_cavgs_reject_restart_002
...
```

The completion marker is:

```text
ABINITIO3D_CAVGS_REJECT_FINISHED
```

The command copies the original project file into each restart folder and
removes any existing completion marker before submission. Absolute paths to the
restart project files and marker files are constructed from the original
working directory, restart folder, and project basename.

Each child command line is a sparse copy of the wrapper command line with these
restart-specific settings:

- `prg=abinitio3D_cavgs`
- `projfile=<project basename>`
- `mkdir=no`
- `nstates=<2 or 3>`
- `nstages=2`
- `verbose_exit=yes`
- `verbose_exit_fname=ABINITIO3D_CAVGS_REJECT_FINISHED`

The wrapper-only and output-routing arguments are removed from the child
command line:

- `nrestarts`
- `nparts`
- `numlen`
- `dir_exec`
- `outdir`
- `quality_mode`
- `quality_model`
- `filetab`
- `fname`
- `infile`

The command submits each restart asynchronously through the queue environment.
The child keeps shared-memory parallelization through the parsed `nthr` value,
while the wrapper controls parallelism across restarts by launching multiple
asynchronous jobs.

After all submissions, the command watches the absolute marker-file paths. If
any marker is missing after the watcher returns, the workflow stops.

## 5. Restart Label Collection

For each finished restart, the command reads only the `cls3D` segment from the
restart project. The restart `cls3D` count must match the original `cls2D`
count.

Restart labels are read from `cls3D%state` and sanitized per class:

- if the original `cls2D` state was less than or equal to zero, the restart
  label is set to `0`;
- if the restart label is outside `1:nstates`, it is set to `0`;
- otherwise, the restart label is kept.

The sanitized raw labels are retained in `restart_labels` and written to the
consensus report.

## 6. State-Label Correspondence

Restart 1 defines the reference label space. Its sanitized labels are copied
directly into `mapped_labels(1,:)`.

For every later restart, the command enumerates all state-label permutations
and chooses the mapping with the highest agreement to restart 1. For restart
`r` and candidate permutation `P`, the score is:

```text
score(P, r) =
    count over class averages i where
        restart_labels(1,i) > 0
        restart_labels(r,i) > 0
        restart_labels(1,i) == P(restart_labels(r,i))
```

For `nstates=2`, the two permutations are tested. For `nstates=3`, all six
permutations are tested. The selected permutation maps raw restart labels into
the restart-1 consensus label space. Label `0` remains `0`.

If two permutations have the same score, the first permutation encountered by
the ascending enumeration wins.

## 7. Consensus Voting

For each active original class average, the command counts mapped restart
labels:

```text
votes(label, class) =
    number of restarts whose mapped label for class is label
```

Labels outside `1:nstates` do not vote. Original classes with state less than
or equal to zero are skipped and keep consensus state `0`.

The consensus state is the label with the largest vote count. Ties are broken
in favor of the mapped label from restart 1 when that label has the same vote
count as the current best count. If an active class has no valid votes, the
current implementation falls through to the first consensus label because all
vote counts are zero.

The full vote vector is retained. The project annotation stores only the
maximum vote count as `cavgs_reject_votes`; the per-state vote counts are
written in the consensus report.

## 8. Good/Bad State Assignment

After consensus voting, the command chooses one consensus state as good using
the cavg-quality scores.

For each consensus label, it computes:

```text
mean_quality(label) =
    average quality_score over class averages where
        original_state > 0
        consensus_state == label
```

Only populated consensus labels are eligible. The populated label with the
highest mean quality is selected as the good consensus state. If two populated
labels have identical mean quality, the lower-numbered label wins because the
implementation updates the winner only on a strict improvement.

Final internal states are assigned as:

- `0` for inactive classes with consensus state `0`;
- `1` for classes in the good consensus state;
- `2` for classes in every other active consensus state.

The project-facing binary selection is:

```text
selection_state = 1 if final_state == 1
selection_state = 0 otherwise
```

For `nstates=3`, this means one consensus state is selected and both remaining
active consensus states are rejected.

## 9. Project Mapping

The command calls `spproj%map_cavgs_selection(selection_states)`.

That project method:

- requires the selection array length to match the `cls2D` row count;
- sets each `cls2D` `state` to the corresponding selection value;
- creates or resizes `cls3D` to match `cls2D` when needed;
- sets each `cls3D` `state` to the corresponding selection value;
- when both `ptcl2D` and `ptcl3D` are present, maps each class state to all
  particles whose `ptcl2D%class` equals that class index.

After mapping, `abinitio3D_cavgs_reject` annotates `cls2D` with:

- `quality`
- `accept`
- `quality_cluster`
- `cavgs_reject_consensus`
- `cavgs_reject_votes`

If `cls3D` has the same row count as `cls2D`, the same annotations are written
to `cls3D`.

If `prune=yes`, `spproj%prune_particles` is called after selection mapping and
annotation. The project is then written back to `projfile`.

The consensus volume is registered in `os_out` as a `vol_cavg` entry with
`state=1` before the project is written.

## 10. Consensus Volume

After the good consensus state is known, the command identifies one restart
state volume per restart. The raw state selected for restart `r` is the raw
state whose mapped label equals the good consensus state.

For each restart, the command computes the mean cavg-quality score over active
class averages whose mapped restart label equals the good consensus state. The
restart with the highest such mean is used as the reference volume. Restarts
with no active class averages in the good mapped state cannot become the
reference.

The selected restart volume path is:

```text
abinitio3D_cavgs_reject_restart_NNN/vol_stateXX.mrc
```

where `XX` is the raw restart state corresponding to the good consensus state.

Every non-reference selected restart volume is docked to the reference with an
asynchronous `dock_volpair` child job. Docking jobs run in folders named:

```text
abinitio3D_cavgs_reject_dock_NNN
```

and signal completion with:

```text
ABINITIO3D_CAVGS_REJECT_DOCK_FINISHED
```

The docking child receives the selected reference volume as `vol1`, the
selected target restart volume as `vol2`, the stage-volume sampling read from
the reference volume header as `smpd`, the default docking band-pass range
`hp=100` and `lp=15`, the parsed `mskdiam`, the parsed `nthr`, and `mkdir=no`.

The reference volume is used as-is. Docked target volumes are written as:

```text
abinitio3D_cavgs_reject_dock_NNN/consensus_docked_restart_NNN.mrc
```

After all docking markers appear, the command reads the reference and docked
target volumes, requires matching dimensions and sampling, averages them with
equal weight, and writes:

```text
abinitio3D_cavgs_reject_consensus_vol.mrc
```

The volume report `abinitio3D_cavgs_reject_consensus_volume.txt` records the
reference restart, selected raw state per restart, per-restart state quality
mean and population, selected volume path, docked volume path, docking LP/HP,
and per-dock report path.

## 11. Output Files

The command writes the same selected/rejected stack style used by
`model_cavgs_rejection` apply mode:

- `quality_selected_cavgs.mrc`
- `quality_rejected_cavgs.mrc`

It also writes:

- `cavgs_quality_features.txt`
- `abinitio3D_cavgs_reject_consensus.txt`
- `abinitio3D_cavgs_reject_consensus_volume.txt`
- `abinitio3D_cavgs_reject_consensus_vol.mrc`

The feature table is written through `write_cavg_quality_feature_table` with
`manual_states=selection_states`, so the accepted/rejected states in the table
are the final consensus-derived binary selection.

The consensus report records:

- class index;
- original state;
- consensus state;
- internal final state;
- binary selection state;
- one vote count column per restart state;
- quality score;
- quality-model automatic state;
- raw restart label for every restart;
- mapped restart label for every restart.

## 12. Failure Conditions

The workflow stops when:

- `part` is supplied;
- `nstates` is not `2` or `3`;
- `nstages` is supplied with a value other than `2`;
- `nrestarts < 1`;
- the input project has no `cls2D` rows;
- the saved original-state array does not match the `cls2D` row count;
- the class-average stack size does not match the `cls2D` row count;
- any restart completion marker is missing after the watcher returns;
- any restart `cls3D` row count differs from the original `cls2D` count;
- no populated consensus class can be chosen as good;
- no raw restart state can be mapped to the good consensus state;
- a selected restart state volume is missing;
- any asynchronous docking completion marker is missing;
- any docked consensus volume is missing;
- docked volumes have inconsistent dimensions or sampling.

## 13. Current Limits

State-label correspondence is based on label agreement with restart 1. It does
not compare maps, class-average projections, or inter-state volumes. Volume
docking is applied only after the good consensus state has been chosen.

The cavg-quality model selects which consensus state is good by mean class
quality. It does not directly override individual class votes except through
that state-level good/bad assignment. The same class-average quality scores are
used to choose the reference restart volume for docking.

The workflow supports at most three restart states. The final project state is
binary regardless of the restart state count.

The route is fixed to C1 because it exits after stage 2 of
`abinitio3D_cavgs`.
