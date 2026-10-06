# Reference Generation Policy

How the stream makes the picking references that reference-based picking uses, who may publish
them, and what a change must preserve.

## 1. Scope

- The initial analysis (p03): `src/main/stream/stages/simple_stream_stage_initial_analysis.f90`,
  driven by `src/main/commanders/stream/simple_commanders_stream_p03_initial_analysis.f90`.
- The master's wiring of p03 and p04: `src/main/commanders/stream/simple_commanders_stream_p00_master.f90`.
- Reference picking's use of the references (p04): `src/main/stream/stages/simple_stream_stage_refpick.f90`.

The particle sieve that p03 drives is governed by
`doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`; the class-average quality selection by
`doc/policies/sieving_and_rejection/model_cavgs_rejection.md`.

## 2. Routes to the references

There are three, and exactly one supplies the references of a run.

1. **Given references.** When the master is started with `pickrefs=<stack>`, p03 is not run
   (skipped) and p04 picks with that stack.
2. **A selection in the GUI (manual route).** The user selects class averages of p03's cycle 1 or
   cycle 2 in the GUI; the master forwards the selection (`pickrefs_selection`, `pickrefs_cycle`)
   to p03, which publishes the selected class averages.
3. **The 3D route (automatic).** p03 runs its plan to the end (section 3) and publishes
   reprojections of an ab initio volume.

## 3. The automatic route

p03 runs a fixed plan of two cycles over the first preprocessed micrographs:

1. **Cycle 1** on the first `NMICS_PLAN(1)` (100) accepted micrographs, or on fewer once
   preprocessing will hand on no more (below): picking with
   `segdiam_bin_picker`, which decides the diameter bins and the box; extraction; `solve2D`;
   class-average selection; the mask diameter estimated from the selected class averages
   (section 3.1). The selected class averages go to the GUI. A pick that accepts no diameter
   bin decides no box: nothing is extracted (the "all" set does not start either), a warning is
   logged, and cycle 1 picks again once `NMICS_PLAN(1)` more micrographs are in, or once
   preprocessing hands on no more and new micrographs came since (follow-up plan, decision 4).
2. **The "all" set**, as soon as the bins are known: every imported project up to
   `NMICS_PLAN(2)` (500) micrographs is picked and extracted with the same bins and box, and its
   particles are fed to a coarse-only particle sieve. The cap is checked before each project, so
   one pass never goes past it.
3. **Cycle 2**, once the sieve has every particle: the sieve's chunks are combined; `solve2D`;
   class-average selection
   (the class averages go to the GUI); class balancing; `solve3D_cavgs`; the choice of a
   state; its reprojections, rescaled to the particle sampling, are published. The 2D status of
   cycle 2 counts cycle 2's particles.

The sieve's final ingestion is set once every extraction of the "all" set is collected and
either `NMICS_PLAN(2)` micrographs are picked, or preprocessing will hand on no more and every
imported micrograph is picked. "No more" is preprocessing's `STREAM_IDLE` (no new movie for 15
minutes, its sets done) or `STREAM_FINISHED` marker in its folder, counted once a watch made a
settle time after it was first seen has found nothing new; a master that skips preprocessing
(`dir_preprocess`) marks the existing folder finished. Right after setting it, p03 cycles the
sieve once more, so the leftover chunk exists before cycle 2 asks whether the sieve is done. A
session with fewer micrographs than the plan's therefore still reaches references; p03 never
takes the trigger back, and the wait is logged once.

The job settings are fixed in the stage, except the 3D route's and the jobs' resources, which
are inputs (decision 20): `nstates_pickrefs`, `nstages_pickrefs`, `lpstop_pickrefs`,
`nspace_pickrefs`, `nthr3D_pickrefs`, and `nrestarts_collapse`, `lpstart_ini3D`, `lpstop_ini3D`
under their `solve3D_cavgs` names. The master offers them at developer visibility and forwards each
only when given (`nthr3D_pickrefs` comes from its resources table otherwise); p03's commander gives
the defaults below and checks `nstates_pickrefs` ≥ 2, `lpstart_ini3D` > `lpstop_ini3D` ≥
`lpstop_pickrefs` > 0 and positive counts. The resources (`nthr2D`, `nparts`, `nchunks`,
`nthr3D_pickrefs`) come from the master's resources table (`doc/policies/stream/README.md`).

| Job | Settings |
|---|---|
| `solve2D` (both cycles) | `ncls` = particles / `nptcls_per_cls`, clamped to 10..100; `nsample` = max(2000, particles/5 rounded up to 1000); `lpstop=8`; `mskdiam`: `mskdiam_box` in cycle 1, the estimate in cycle 2 (section 3.1); `center=yes`; `autoscale=yes`; `sigma_est=global`; `nthr2D` (16) threads, `nparts` (1) |
| sieve of the "all" set | coarse only (`single_pass=yes`); `nchunks` (4) chunks of `nthr2D` (16) threads, on `SIMPLE_STREAM_REFGEN_PARTITION`; `mskdiam_box`, also as the scoring mask; no starting low-pass (the chunks' `solve2D` derives it) |
| class balancing | the selected class averages replicated in proportion to population up to `TARGET_NCLS` (501) rows, written to a project of its own (`balance_classes/all/all_balanced.simple`) on which `solve3D_cavgs` runs; the cycle 2 project keeps the class averages the GUI shows, so a selection of cycle 2 is read against them |
| `solve3D_cavgs` | `nstates_pickrefs` (3), `nstages_pickrefs` (3), `nrestarts_collapse` (3), `lpstart_ini3D` (100), `lpstop_ini3D` (20), `lpstop_pickrefs` (8), `pgrp=c1`, `prune=no`, the estimated mask diameter, `nthr3D_pickrefs` (16) threads, `nparts` only when above 1: `nparts` on the command line makes refine3D distributed even at 1, and a distributed iteration that empties a state stops the next one, while in shared memory `solve3D_cavgs` handles the collapse after the stage. Its restart driver names the run whose result stands in `SOLVE3D_CAVGS_FINAL_DIR`, which p03 reads |
| state choice | among the populated states with a volume (`choose_state`, `state_vetoes`; follow-up plan, decisions 11-14). Each binarised volume (low-passed to 20 Å, Otsu twice with both thresholds from the voxels inside the mask, so the solvent outside does not drive them; `vol_shape_descr`) gives its largest component's share of the foreground inside the mask. A state passes the vetoes when it is one object (that share at least `STATE_DOMINANT_FRAC`, 0.8, so a speck or a small split does not fail it), holds at least `STATE_POP_FLOOR` (10%) of the candidates' population, and its FSC 0.143 resolution (its FSC file in the 3D result) is within `STATE_RES_FACTOR` (1.5) times the best candidate's; a state without an FSC file is not vetoed on resolution. Among those, the most distinct projection directions of its classes (`os_cls3D` `proj`), then the largest population, then the lowest state. When none passes, the same order over every candidate, with a warning to the log and the GUI (the 3 October order, fewest components first, preferred compact junk and is gone). When the directions cannot be counted, the population breaks the ties. Every key of every state is logged. The three constants are provisional until a validation run on datasets with a known answer sets them. Each state's binarised volume and components are written as `vol_binarized_stateNN.mrc` and `vol_cc_stateNN.mrc` |
| reprojection | `nspace_pickrefs` (50), `pgrp=c1`, the estimated mask diameter, `nthr3D_pickrefs` threads; a job on the local queue |

Class-average selection in both cycles is the chunk quality model (`score_project_cavgs`)
followed by the class compatibility filter, trained and applied on the same selection. Its mask
diameter is the one of the cycle's `solve2D`.

### 3.1 The mask diameter

1. **Cycle 1 masks with the box's default**, `mskdiam_box = (box − COSMSKHALFWIDTH) · smpd` for
   the picker's box and the micrographs' pixel size: the value `parameters` derives when no mask
   is given. It is passed explicitly, since `solve2D` requires `mskdiam`.
2. **Cycle 1's selected class averages give the estimate.** These are the classes that the
   quality model and the compatibility filter keep. The estimate is generous, and it uses the
   measure and the rule `make_pickrefs` applies to its references:
   - each selected class average is automasked (`automask2D` with `gen_pickrefs`' `ngrow`,
     `winsz`, `amsklp` and `edge`, which its commander defaults to `make_pickrefs`' 3, 5, 20 Å
     and 6), which gives the diameter of its largest connected component;
   - the largest of these diameters is taken, so that every view fits;
   - it is widened by two soft edges, rounded to an even box, capped at the class averages' box
     and multiplied by `MSK_EXP_FAC` (1.2) (`automask2D_mskdiam`);
   - the result is capped at `mskdiam_box`.

   Example: a largest automask of 150 Å at 1.3 Å/px gives a 128 px particle box (166 Å) and a
   mask of 200 Å. The box default is usually wider, because the picker makes the box 1.0 to 1.5
   times the largest diameter of its accepted bins. For large particles, where the factor nears
   1.0, the cap applies. The estimate is logged and sent to the GUI as the mask diameter.
3. **When cycle 1 selects no class, the estimate is `mskdiam_box`**, with a warning. A selected
   class whose automask finds no object counts as a disc nearly the size of the box
   (`automask2D`'s fallback). That puts the estimate at or near the box default.
4. **Cycle 2 and 3D use the estimate.** That covers cycle 2's `solve2D` and class-average
   selection, `solve3D_cavgs`, the volume shape descriptors and the reprojection.
5. **The sieve of the "all" set masks with `mskdiam_box`.** It starts as soon as the bins are
   known, before cycle 1 has selected classes.
6. **A restart estimates again.** A restarted p03 without published references clears its working
   folders and runs its plan from cycle 1 (`doc/policies/stream/restart_policy.md`).
7. **Particle sieving (p05) does not use this estimate.** It reads the mask diameter that
   `make_pickrefs` (p04) writes to `moldiam.txt`, measured on the published references by the
   same rule.

## 4. Output contract

1. The references are one stack, `INITIAL_ANALYSIS_PICKREFS` (`selected_references.mrcs`,
   `simple_defs_stream`), in p03's folder (`INITIAL_ANALYSIS_JOB_NAME`). The master gives p04
   `pickrefs=../initial_analysis/selected_references.mrcs` built from the same two constants; no other
   name may be used on either side.
2. They are published once per run by `publish_pickrefs`: a complete stack is renamed into place,
   so p04, which polls for the file, never reads a partial stack. Nothing writes the file
   in place.
3. Once published they are final. `publish_pickrefs` refuses to publish over them, and no other
   code writes the file.
4. Publishing sends the stack and its sprite sheet (`selected_references.jpg`) to the GUI as the
   picking references.
5. p04 reads the file once: `make_pickrefs` (`ncls=10`, `nrots=12`, `mirr=yes`), a job on the
   local machine in p04's `make_pickrefs` folder, makes its picking templates at the pixel size
   of the first upstream project with an accepted micrograph, low-passed to 0.15 of the largest
   class-average diameter within [15, 30] Å (`template_lowpass`), and writes `moldiam.txt`. A
   `box_extract` given on the command line is used as given, with a warning below 128 px. p04
   renames the templates and then `moldiam.txt` into its folder, so `moldiam.txt` exists only
   once the templates do. Its `box_for_extract` is the extraction box and its mask diameter is the
   one the particle-sieving stage reads. The upstream projects wait in p04's watcher until the
   templates are in place.
6. p03's reprojection of the chosen state is a job on the local machine too, in the 3D result's
   `reproject` folder; its reprojections are rescaled and published from there.

## 5. Precedence

1. **The user's selection wins until references are published.** Each p03 pass reads the GUI's
   updates before the cycle steps, so a selection pre-empts a 3D result collected in the same
   pass.
2. **A selection that publishes references ends p03** at once, wherever the plan is. The rest of
   the plan is skipped; jobs already submitted (the "all" extractions, `solve2D`,
   `solve3D_cavgs`) are cancelled when the stage stops, and no result of theirs is published.
3. **A selection that publishes nothing is ignored** with a warning, and p03 keeps running: the
   cycle has no class averages yet, or no index names one of its classes.
4. **A selection larger than the update holds** (`MAX_PICKREFS_SELECTION`, 500) is dropped whole
   by the master (`doc/policies/stream/ipc_policy.md`).
5. **Published references are final, including across restarts**: a restarted p03 that finds
   them sends them to the GUI again and is finished at once (`restore_pickrefs`).
6. **After the 3D route has published, no selection can replace the references**: p03 has stopped
   (its last status is `terminating`), and the master sends updates only to running stages.

## 6. Change rules

- Keep one publisher (`publish_pickrefs`) and one name (`INITIAL_ANALYSIS_PICKREFS`).
- Keep the GUI updates ahead of the cycle steps in `iterate`.
- A change that lets references change after publication must also make p04 re-make its templates
  for the sets submitted afterwards, and must say how particles picked with the old references
  are treated downstream.
- A change of the plan or the job settings updates section 3.
- Cycle 1 and the sieve mask with `mskdiam_box`. The mask of cycle 2 and 3D comes only from
  cycle 1's selected class averages and never exceeds `mskdiam_box`.
- p03 and `make_pickrefs` measure with `automask2D` (p03 with `gen_pickrefs`' automasking
  inputs, defaulted to `make_pickrefs`' values), and widen with one routine,
  `automask2D_mskdiam`. A change of the rule changes both and updates section 3.1.
- Tests: `test_gui_selection_ends_stage`, `test_published_pickrefs_are_final`,
  `test_estimate_mskdiam` and `test_choose_state` in
  `src/main/stream/stages/simple_stream_stage_initial_analysis_tester.f90`
  (`unit_stream`, "initial analysis"); `test_automask2D_mskdiam` in
  `src/main/image/simple_image_msk_tester.f90` (sub-suite "masks").

## 7. Known gaps

- **Class averages are duplicated** to weight `solve3D_cavgs` (review M1). Each copy carries
  its class's full population, the cycle-2 project is rewritten with the 501-row table, and copies
  are aligned and assigned to states independently. The review proposes per-row weights in the
  reconstruction instead (R2).
- **The references are reprojections of a multistate model** (M2). The settings are inputs and
  the result path is fixed (R4, done on 5 October 2026), but their defaults are those of
  3 October and no validation has set them.
- **The mask estimate rests on cycle 1's micrographs** (about 100). A view that is rare there,
  or a larger particle that appears only later, can be cut by the mask of cycle 2 and 3D. The
  generous rule is the only margin, and nothing estimates again.
- **The "all" sieve masks with the box default**, which is wider than the particle, because it
  starts before the estimate exists.
- **The state choice's fractions are guesses.** `STATE_CC_MIN_FRAC` and `STATE_POP_FLOOR` (10%
  each) wait for a validation run; a speck larger than 10% of the particle still makes a second
  object.
- **The quality selection is fitted and applied on the same classes** (M5).
- **Jobs left running after a user selection** run until the stage stops, when they are
  cancelled; a job still queued then, or on another host without a scheduler id, is not
  (`simple_qsys_job_record`).
- **No validation** compares the routes' references on a known dataset.
