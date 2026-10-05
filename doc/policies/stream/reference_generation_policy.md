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

1. **Cycle 1** on the first `NMICS_PLAN(1)` (100) accepted micrographs: picking with
   `segdiam_bin_picker`, which decides the diameter bins and the box; extraction; `solve2D`;
   class-average selection; the mask diameter estimated from the selected class averages
   (section 3.1). The selected class averages go to the GUI.
2. **The "all" set**, as soon as the bins are known: every imported project up to
   `NMICS_PLAN(2)` (500) micrographs is picked and extracted with the same bins and box, and its
   particles are fed to a coarse-only particle sieve.
3. **Cycle 2**, once the sieve has every particle (final ingestion is set when the "all" set is
   picked and extracted): the sieve's chunks are combined; `solve2D`; class-average selection
   (the class averages go to the GUI); class balancing; `solve3D_cavgs`; the choice of a
   state; its reprojections, rescaled to the particle sampling, are published.

The job settings are fixed in the stage:

| Job | Settings |
|---|---|
| `solve2D` (both cycles) | `ncls` = particles / `nptcls_per_cls`, clamped to 10..100; `nsample` = max(2000, particles/5 rounded up to 1000); `lpstop=8`; `mskdiam`: `mskdiam_box` in cycle 1, the estimate in cycle 2 (section 3.1); `center=yes`; `autoscale=yes`; `sigma_est=global`; `nthr=16`, `nparts=1` |
| sieve of the "all" set | coarse only (`single_pass=yes`); `nchunks=4`, `nmics=100`, `nthr=16`; `mskdiam_box`; no starting low-pass (the chunks' `solve2D` derives it) |
| class balancing | the selected class averages replicated in proportion to population up to `TARGET_NCLS` (501) rows |
| `solve3D_cavgs` | `nstates=3`, `nstages=3`, `nrestarts_collapse=3`, `lpstart_ini3D=100`, `lpstop_ini3D=20`, `lpstop=8`, `pgrp=c1`, `prune=no`, the estimated mask diameter, `nthr=16` |
| state choice | among the populated states with a volume (`choose_state`): the fewest connected components of the binarised volume (low-passed to 20 Å, Otsu twice; `vol_shape_descr`), so a single object wins; then the most distinct projection directions of its classes (`os_cls3D` `proj`); then the largest population; then the lowest state. A volume with no component ranks last; when the directions cannot be counted, the population breaks the ties |
| reprojection | `nspace=50`, `pgrp=c1`, the estimated mask diameter |

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
   - each selected class average is automasked (`automask2D` with its defaults: `ngrow=3`,
     `winsz=5`, `amsklp=20`, `edge=6`), which gives the diameter of its largest connected
     component;
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
6. **A restart estimates again.** A restarted p03 runs its plan from cycle 1 and keeps nothing
   (`doc/policies/stream/restart_policy.md`).
7. **Particle sieving (p05) does not use this estimate.** It reads the mask diameter that
   `make_pickrefs` (p04) writes to `moldiam.txt`, measured on the published references by the
   same rule.

## 4. Output contract

1. The references are one stack, `OPENING2D_PICKREFS` (`selected_references.mrcs`,
   `simple_defs_stream`), in p03's folder (`OPENING2D_JOB_NAME`). The master gives p04
   `pickrefs=../opening_2D/selected_references.mrcs` built from the same two constants; no other
   name may be used on either side.
2. They are published once per run by `publish_pickrefs`: a complete stack is renamed into place,
   so p04, which polls for the file, never reads a partial stack. Nothing writes the file
   in place.
3. Once published they are final. `publish_pickrefs` refuses to publish over them, and no other
   code writes the file.
4. Publishing sends the stack and its sprite sheet (`selected_references.jpg`) to the GUI as the
   picking references.
5. p04 reads the file once: `make_pickrefs` (`ncls=10`, `nrots=12`, `mirr=yes`) makes its picking
   templates at the micrographs' pixel size and writes `moldiam.txt`, whose `box_for_extract` is
   the extraction box and whose mask diameter the particle-sieving stage reads.

## 5. Precedence

1. **The user's selection wins until references are published.** Each p03 pass reads the GUI's
   updates before the cycle steps, so a selection pre-empts a 3D result collected in the same
   pass.
2. **A selection that publishes references ends p03** at once, wherever the plan is. The rest of
   the plan is skipped; jobs already submitted (the "all" extractions, the sieve's chunks,
   `solve2D`, `solve3D_cavgs`) run on unattended, and no result of theirs is published.
3. **A selection that publishes nothing is ignored** with a warning, and p03 keeps running: the
   cycle has no class averages yet, or no index names one of its classes.
4. **A selection larger than the update holds** (`MAX_PICKREFS_SELECTION`, 500) is dropped whole
   by the master (`doc/policies/stream/ipc_policy.md`).
5. **Published references are final, including across restarts**: a restarted p03 that finds
   them sends them to the GUI again and is finished at once (`restore_pickrefs`).
6. **After the 3D route has published, no selection can replace the references**: p03 has stopped
   (its last status is `terminating`), and the master sends updates only to running stages.

## 6. Change rules

- Keep one publisher (`publish_pickrefs`) and one name (`OPENING2D_PICKREFS`).
- Keep the GUI updates ahead of the cycle steps in `iterate`.
- A change that lets references change after publication must also make p04 re-make its templates
  for the sets submitted afterwards, and must say how particles picked with the old references
  are treated downstream.
- A change of the plan or the job settings updates section 3.
- Cycle 1 and the sieve mask with `mskdiam_box`. The mask of cycle 2 and 3D comes only from
  cycle 1's selected class averages and never exceeds `mskdiam_box`.
- p03 and `make_pickrefs` measure with `automask2D` and its defaults (`simple_default_clines`),
  and widen with one routine, `automask2D_mskdiam`. A change of the rule changes both and
  updates section 3.1.
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
- **The references are reprojections of a three-state model** at `nspace=50` (M2). The settings
  are literals, not inputs. p03 finds the 3D result by taking the highest-numbered
  `<n>_solve3D_cavgs` directory. The review proposes inputs and a fixed result path (R4).
- **The mask estimate rests on cycle 1's micrographs** (about 100). A view that is rare there,
  or a larger particle that appears only later, can be cut by the mask of cycle 2 and 3D. The
  generous rule is the only margin, and nothing estimates again.
- **The "all" sieve masks with the box default**, which is wider than the particle, because it
  starts before the estimate exists.
- **Every connected component counts in the state choice.** A speck of noise above the threshold
  makes a volume two components, and it then loses to a single-object volume with far fewer
  views. The count is of the whole box, not of the mask. Each state's binarised volume and
  components are written under one name (`vol_binarized.mrc`, `vol_cc.mrc`), so only the last
  state's remain to check a choice against.
- **The quality selection is fitted and applied on the same classes** (M5).
- **Jobs left running after a user selection** cannot be cancelled (`qsys_async_job` has no
  cancel).
- **No validation** compares the routes' references on a known dataset.
