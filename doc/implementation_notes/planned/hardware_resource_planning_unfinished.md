# Hardware resource planning: status and remaining work

The external Python planner calibrates motion correction, CTF estimation,
picking, solve2D and solve3D using installed SIMPLE binaries. It writes JSON
recommendations; SIMPLE does not consume them automatically (`auto_apply=false`).
Extraction is excluded. See the [usage guide](../../../scripts/resource_planner/README.md).

## Current behavior

- Quick calibration provides initial synthetic recommendations. Full calibration
  increases workload, repetitions and thread coverage; scientific gates remain required.
- Slurm comparisons use native partition workers and master merging. Worker
  placement, allocation identity, completion and scientific outputs are checked.
  Picking remains single-node. Node trials hold partition count fixed and use one
  worker per node; they do not exhaustively optimize workers and partitions.
- Optional storage probes report synchronized writes and verified first/repeated
  reads. Slurm storage results cover only the coordinator node.
- Failed trials retain diagnostics and return a nonzero exit status.

## Verified scope

All five quick workflows passed on the desktop and CSB. A bounded CSB comparison
passed for motion correction, CTF estimation, solve2D and solve3D on one and two
nodes, including assembled-output checks. The Python regression suite passed
70 tests. These checks establish synthetic behavior, not production optima.

## Remaining work

- Complete full-workload and four-node validation; qualify representative real data.
  Quick solve3D uses a supplied volume and does not qualify de novo reconstruction.
- Measure concurrent per-node memory; current distributed estimates conservatively
  sum process peaks. Extend partition, NUMA, GPU and storage-contention coverage.
- Add particle CTF correction and streaming calibration if requested.
- Design opt-in native JSON consumption separately. Explicit job parameters must
  take precedence; incompatible profiles must leave defaults unchanged, and restart
  partition layouts must be preserved.

Historical design details are in the
[completed implementation record](../completed/hardware_resource_planning_completed.md).
