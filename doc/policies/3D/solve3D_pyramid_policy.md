# solve3D_pyramid Policy

This document records the contract of `solve3D_pyramid`: multistate 3D by class-average
assignment, in one pass, on a 2D-classified project. It does offline what stream p07 does
session-long, with p07's rules and none of its stream machinery. The rationale of the method
(why sub-class averages, why the check, why the rejection by fit) is in
[stream_3D_ingestion_policy.md](../stream/stream_3D_ingestion_policy.md), sections 9 and 10. The
design record is
[solve3D_pyramid_plan_2026-10-15.md](../../refactoring_notes/planned/solve3D_pyramid_plan_2026-10-15.md).

Status (2026-10-15): implemented, not yet compiled or run.

## 1. Scope

`simple_exec prg=solve3D_pyramid projfile=<project> nstates=<n> nthr=<n>` makes the states'
reference volumes once, from class averages of the first `nptcls_refs` selected particles. It
then assigns every selected particle to a state by matching sub-class averages against those
references. Last, one `refine3D refine=cont` iteration refines the particles' poses with the
states fixed.

It has no watching, publications, retries, fallbacks or GUI. Every failed step is a hard exit,
after the step's log and files are on disk and the report records the failure. Nothing is
retried, skipped or replaced. The one lenience kept is that a pool model that keeps no class of
the first set means every populated class is used, with a warning. That outcome is a verdict on
the classes, not a failure of a step.

## 2. Command Line

| Key | Status |
| --- | --- |
| `projfile` | required: a 2D-classified project (section 3) |
| `nstates` | required, at least 2 (`NSTATES_MIN`) |
| `nthr` | required; every step's threads |
| `nptcls_refs` | optional, default 100,000 (`NPTCLS_REFS`): the first set's selected particles at most |
| `use_model` | optional, `yes` (default) or `no`: the pool model selects the first set's classes |
| `mskdiam` | optional; default the project's class averages' (`get_mskdiam('cavg')`). With neither, the program exits |
| `nparts` | optional; makes `solve2D`, both `cls_expansion` runs and `solve3D_cavgs` distributed when above 1. `refine3D` always runs distributed, with `nparts` parts (1 by default). `assign_cavgs_states` runs in shared memory |
| `mkdir` | `yes` (forced); `mkdir=no` is refused, since the steps live in the execution directory |

`nspace` (the local refinement's grid, default 2500) must be even. Every other setting is a
constant of `simple_state_assignment`, shared with p07 (section 6).

## 3. Input

- `ptcl2D` holds a 2D classification (refine2D, solve2D or the like), and `ptcl3D` has as many
  rows.
- The class averages are registered in `os_out` (one `cavg` entry), with one class average per
  `cls2D` entry, their even and odd halves and the `frc2D` FRCs on disk (`usable_job_2D`).
- A particle selected in one particle segment only is deselected in both, and this is logged.

## 4. Steps

Each step runs in-process in a folder of the execution directory. Its commander executes on a
child command line built by the shared `set_*_cline` builders, then the program returns to the
execution directory. The two exceptions are the first set's cut and the check, which run in the
program itself.

| # | Folder | What | Exit when |
| --- | --- | --- | --- |
| 1 | `refs_subset/` | the first set: whole stacks in arrival order up to `nptcls_refs` selected particles (`cut_whole_stacks`), the other particles deselected in a copy of the project | fewer than `MIN_PTCLS_PER_STATE` (5) × `nstates` |
| 2 | `refs_subset/solve2D/` | `solve2D`, `NCLS2D_FIRST` (150) classes, `cache=yes` | no project, or no usable 2D (`usable_job_2D`) |
| 3 | `refs_subset/quality_selection/` | with `use_model=yes`, the pool model's selection of those classes, mapped to their particles; the selected and rejected stacks with JPEGs (`select_route_classes`) | never |
| 4 | `refs_subset/cls_expansion/` | `cls_expansion` of the kept classes into `expansion_budget` (500) sub-classes in all | no project, or no usable 2D |
| 5 | `refs_subset/solve3D_cavgs/` | `solve3D_cavgs`, `nstates`, c1, `prune=no`, 20 → 12 Å forced, `prob_athres_taper=yes` | no project |
| 6 | `state_refs/` | the check (`check_state_refs`): the states that hold class averages renumbered from 1, their volumes, JPEGs, the class averages per state and the record `state_refs.txt` | no state, a state without a volume, or the check fails |
| 7 | `assign/quality_selection/` | the pool model's verdict on the project's classes, logged beside the fit, not applied | never |
| 8 | `assign/cls_expansion/` | `cls_expansion` of every selected particle's classes (the project's own 2D) into 500 sub-classes | no project, or no usable 2D |
| 9 | `assign/` | `assign_cavgs_states` against `state_refs/` (15 Å, 1000 directions per state) | no project |
| 10 | `assign/` | the rejection by fit, the ambiguity flags and the mapping onto the particles (`map_assignment`); the per-state, rejected and ambiguous stacks with JPEGs (`write_matched_cavgs`) | no candidate sub-class, or other particles |
| 11 | `refine3D/` | the seeds (`set_refine_seeds`: the project's half-sets, the nearest of `nspace` directions), the particles' states renumbered to the populated ones (`write_refine_project`), `refine3D refine=cont maxits=1 objfun=cc pgrp=c1` | no particle with a state, no project, or no volume registered |
| 12 | (the execution directory) | the result written over the execution directory's project; the report | — |

The check passes when the same-state even/odd median angular distance is at most 50°, the state
agreement is at least 70 %, no two states correlate above 0.97 after docking in both hands, and
at least two states each hold 5 % of the class averages (`refs_usable`). A failed check leaves
`state_refs/` with the record, the volumes and the per-state class averages for inspection.

When `nptcls_refs` is at least the selected particles, the first set is the whole project. The
route then classifies and expands every particle, and step 8 expands them again from the
project's own 2D. The result is correct, but the second expansion repeats work.

## 5. Outputs

- **The execution directory's project** (the input's name, copied in by `mkdir=yes`; this is
  where `simple_exec` records the job). Every particle has its state (0 for the fit's
  rejections) and its refined pose. Each refined state's volume and FSC are registered (`vol`,
  `fsc`), in the states of the references (the renumbering of step 11 is mapped back).
- **`state_refs/`**: the references, their record, the reprojection and orientation JPEGs, and
  the class averages each was made from.
- **`solve3D_pyramid_report.txt`**: a text table, rewritten after each step and on a failure.
  It holds the status, the source and command line, the selected particles and the first set,
  and the route's classes and sub-classes. It holds the check: medians, agreement, each pair's
  docked correlation, the class averages per state and the verdict. It holds the assignment:
  the candidates, the fit's rejections and cap, the ambiguous sub-classes, the median, MAD and
  threshold, and the particles assigned and deselected, with the fit's and the model's
  rejections side by side. Per state it lists the sub-class averages, the particles and the
  FSC 0.143 and 0.5 resolutions. Last come the wall times of the steps.

## 6. Shared code

The method is not copied from p07. The stateless parts live in `simple_state_assignment`
(`src/main/solve/`), which both programs use:

- the constants: `MIN_PTCLS_PER_STATE`, `NCLS2D_FIRST`, `NCLS_EXPANDED`, `LPSTART_CAVGS3D`,
  `LPSTOP_CAVGS3D`, the check's (`EO_ANGDIST_MAX`, `EO_OVERLAP_MIN`, `STATE_FRAC_MIN`,
  `NSTATES_MIN`, `REFS_PAIR_CC_MAX`, `DOCK_HP_REFS`), the matching's (`LP_ASSIGN`,
  `NSPACE_ASSIGN`) and the fit's (`FIT_NMAD`, `FIT_REJECT_MAX`, `MARGIN_AMBIGUOUS`);
- the rules: `refs_usable`, `reject_by_fit`, `flag_ambiguous`, `compact_states`,
  `expansion_budget`, `count_active_classes`, `fit_mskdiam`, `median_of`, `usable_job_2D`,
  `cut_whole_stacks`;
- the steps on projects: `score_cavgs_with_model`, `select_route_classes`,
  `model_rejected_rows`, `check_state_refs`, `docked_pair_cc`, `read_refs_record`,
  `map_assignment`, `write_matched_cavgs`, `set_refine_seeds`, `write_refine_project`,
  `map_refine_result`, `write_state_cavgs`;
- the jobs' command lines: `set_solve2D_route_cline`, `set_clsexp_cline`, `set_cavgs3D_cline`,
  `set_assign_cline`, `set_refine_cont_cline`. The callers add the threads and the parts.

A change of any of these changes p07 and `solve3D_pyramid` together. Their rules are tested on
fixed numbers in `unit_reconstruction` "state assignment rules".

## 7. Tests

- `unit_reconstruction` "state assignment rules" (`simple_state_assignment_tester`): the check,
  the renumbering, the sub-class budget, the mask diameter, the fit rules, the first set's cut at
  whole stacks, the class averages by state, and the jobs' command lines.
- `unit_stream` "solve 3D": the shared steps through p07's stage (the attempts, the
  assignment, the local refinement's seeds and result); `lib_stream` "solve 3D reference
  docking": the docking of the check.
- High level, nightly: `solve3D_pyramid_blobs` (`simple_commanders_test_highlevel_pyramid`).
  It uses the three-state phantom of 3,600 particles in 60 exact views that
  `refine3D_states_cavgs` uses, with its class averages made from the exact classes
  (`make_cavgs`), and runs end to end with `nstates=3`. It checks the report, the references
  and the result. It checks the assignment's purity against the planted states, each reference
  docked onto its closest truth, and each refined map docked onto its particles' majority truth
  (the hand free). Its floors are provisional, to be set from the first run.
  - The phantom is one stack ordered by direction, so a first set of a third would hold a cap
    of 20 directions. The test therefore takes the whole phantom as the first set (the
    default `nptcls_refs`).
