# Processing heterogeneous data sets in SIMPLE

This guide describes a practical SIMPLE workflow for removing junk, finding
discrete compositional or conformational states, and optionally examining
continuous motion. It is written for a user who has not used the SIMPLE
command line before.

The recommended path is:

```text
picked particles
      |
      +-- not processed by SIMPLE stream --> particle_sieving
      |
      +-- already processed by SIMPLE stream --> skip offline sieving
      |
      v
one or more abinitio2D runs
      |
      +-- optional: abinitio2D_chunks + learned class-rejection model
      |
      v
initial abinitio3D cleanup, normally with three/four states
      |
      v
choose the state(s) to retain
      |
      +-- one state --> continue that state with abinitio3D state=N
      |
      +-- merged states --> select/merge, then single-state abinitio3D
      |
      v
refine3D_states
      |
      v
select each useful state --> refine3D_auto
      |
      +-- optional: flex_pca for remaining continuous heterogeneity
```

The first half of this workflow is **junk removal**. Do not interpret every
class or volume from those runs as a biological state. The second half starts
only after a clean particle set has been established.

## 1. Before running anything

### Check that SIMPLE is available

```bash
which simple_exec
simple_exec prg=list
```

If `which` prints nothing, load the SIMPLE environment used at your site before
continuing.

### Understand a SIMPLE command

SIMPLE commands consist of a program name followed by `key=value` arguments:

```bash
simple_exec prg=abinitio2D projfile=my_project.simple ncls=100 \
  mskdiam=180 nparts=4 nthr=16
```

There must be no spaces around `=`. In this example:

- `projfile` is the input SIMPLE project.
- `ncls` is the requested number of 2D classes.
- `mskdiam` is the particle mask diameter in Angstroms.
- `nparts` is the number of distributed partitions or worker processes.
- `nthr` is the number of shared-memory threads used by each process.

The values of `nparts` and `nthr` depend on the computer or scheduler. Ask the
local SIMPLE administrator for suitable values. Do not copy the example CPU
values without checking the available resources.

Most workflows create a new numbered execution directory and write a new
`.simple` project there. At every step:

1. read the final lines of the log;
2. record the execution directory and output project;
3. inspect the result; and
4. use that output project as `projfile` for the next step.

Never assume that the original project was modified. Keeping a small text file
with every command and its accepted output project prevents most lineage
mistakes.

The commands below are templates. Replace every value in angle brackets,
including the angle brackets themselves.

## 2. Remove picking junk

### 2.1 Decide whether offline sieving is needed

- If the particles have already been processed by **SIMPLE stream**, do not run
  offline sieving again. Start at section 2.2.
- If the particles were picked or imported without SIMPLE stream processing,
  run `particle_sieving` first.

The input project must contain its micrographs, particle stacks, and particle
metadata:

```bash
simple_exec prg=particle_sieving \
  projfile=<EXTRACTED_PROJECT.simple> \
  nparts=<CONCURRENT_CHUNKS> nthr=<THREADS_PER_CHUNK>
```

Here `nparts` means the number of sieving chunks run at the same time, not MPI
partitions within one chunk. Each concurrent chunk uses `nthr` threads. Start
conservatively if memory is limited. The standard coarse and fine sieving
settings are used when their advanced parameters are omitted.

**Checkpoint:** inspect the retained particle count and the class averages.
The output should be substantially cleaner without losing convincing,
well-resolved particle views.

### 2.2 Run ab-initio 2D classification

Run `abinitio2D` on the sieved or stream-cleaned project:

```bash
simple_exec prg=abinitio2D \
  projfile=<SIEVED_OR_STREAM_PROJECT.simple> \
  ncls=<NUMBER_OF_CLASSES> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

Inspect the final ranked class-average stack and reject classes that are
obviously ice, carbon, aggregates, empty boxes, or unrecognizable noise. Run
ab-initio 2D a second time on the retained particles. More than one run is
useful because genuine views should recur, whereas unstable junk classes tend
not to.

The easiest manual-selection route is to use `e2display.py` from EMAN2:

1. Open the final `cavgs_iterNNN_ranked.mrcs` stack from the `abinitio2D`
   execution directory.
2. In the image-stack window, use the middle-mouse menu and `Del` to mark the
   bad classes.
3. Use `Save` from the same menu to write the remaining good classes to a new
   stack, for example `selected_cavgs.mrcs`.
4. Map that subset back to its particles:

```bash
simple_exec prg=map_cavgs_selection \
  projfile=<ABINITIO2D_PROJECT.simple> \
  stk2=<SELECTED_CAVGS.mrcs> prune=yes
```

`map_cavgs_selection` identifies the selected classes by correlation, so the
saved subset does not have to preserve the original class order. Use the
project from the generated numbered `_selection/` directory for the next 2D
run.

If EMAN2 is unavailable, a text selection remains possible. Create a file
containing one `1` (keep) or `0` (reject) per 2D class, in project class order,
and apply it with:

```bash
simple_exec prg=selection \
  projfile=<ABINITIO2D_PROJECT.simple> oritype=cls2D \
  infile=<CLASS_KEEP_FLAGS.txt> prune=yes
```

Use the resulting project from the generated numbered `_selection/` directory
for the next 2D run.

**Checkpoint:** continue only when the retained classes have recognizable
particle features, cover the expected range of views, and no longer contain a
large obvious junk population. Prefer reproducibility over choosing classes
merely because they look sharp.

### 2.3 Optional: classify independent chunks

For a large data set, independent chunk classifications can expose rare junk
or minority populations that disappear in a single global classification:

```bash
simple_exec prg=abinitio2D_chunks \
  projfile=<CLEANING_PROJECT.simple> \
  mskdiam=<MASK_DIAMETER_A> nptcls_per_cls=500 nchunks=0 \
  nchunks_in_parallel=<PARALLEL_CHUNKS> \
  nparts=<PARTITIONS_PER_CHUNK> nthr=<THREADS>
```

`nchunks=0` lets SIMPLE choose enough balanced chunks to target roughly 100
classes per chunk. Peak process use is approximately
`nchunks_in_parallel * nparts`, so set both values deliberately.

A data-set-specific learned rejection model can make repeated chunk cleanup
consistent. Training is an advanced, supervised step: it requires projects in
which the final 2D classes have first been reviewed and labelled by a person.
The complete procedure is in the appendix. Never train a model on its own
unreviewed predictions.

## 3. Perform the initial 3D cleanup

For most heterogeneous data sets, start with three states:

```bash
simple_exec prg=abinitio3D \
  projfile=<SELECTED_CLEAN_2D_PROJECT.simple> \
  nstates=3 pgrp=<POINT_GROUP> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

Use `pgrp=c1` when no point-group symmetry is justified. Do not impose symmetry
only to make a map appear cleaner.

Inspect all three volumes, their particle populations, their directional
coverage and generated re-projections. At this stage the three-state run is
primarily a cleanup tool. A state can represent junk, damaged particles, a
preferred-view failure, or an unstable reconstruction rather than a biological
conformation.

For a large, high-contrast complex with an already clean particle set, it can
be reasonable to use one state and postpone heterogeneity analysis:

```bash
simple_exec prg=abinitio3D \
  projfile=<SELECTED_CLEAN_2D_PROJECT.simple> \
  nstates=1 pgrp=<POINT_GROUP> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

Use this shortcut only when the single-state result is stable and there is no
substantial junk population left to isolate.

### Choose the acceptable initial states

Record which state or states contain the particles to retain. If one state is
sufficient, no separate selection command is needed: section 4 uses SIMPLE's
direct state-continuation mode.

To merge two or more acceptable states into one retained particle group:

```bash
simple_exec prg=selection \
  projfile=<THREE_STATE_ABINITIO3D_PROJECT.simple> \
  oritype=ptcl3D states=<COMMA_SEPARATED_STATES> prune=yes
```

For example, use `states=1,3` to merge states 1 and 3. Run every independent
selection from the original multi-state project, not from a project that has
already been pruned. Keep the original project so that the decision can be
revisited.

## 4. Establish a clean consensus model

If one state was retained from the multi-state cleanup, a complete new
reference-free `abinitio3D` run is unnecessary. Continue directly from that
state:

```bash
simple_exec prg=abinitio3D \
  projfile=<MULTISTATE_ABINITIO3D_PROJECT.simple> \
  state=<STATE_NUMBER> pgrp=<POINT_GROUP> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

This mode internally selects and prunes the requested state, preserves its
particle poses, reconstructs a same-lineage starting map, and resumes the
single-state search at stage 5 with nonuniform filtering. It is therefore a
gentle continuation of the selected solution rather than another full
ab-initio search. Do not supply `vol1` or `nstates`: the state-continuation
path owns preparation of the starting reference.

If two or more states were merged in section 3, there is no single state map
that represents the merged particle group. In that case, establish a new
single-state consensus with:

```bash
simple_exec prg=abinitio3D \
  projfile=<MERGED_STATE_SELECTION_PROJECT.simple> \
  nstates=1 pgrp=<POINT_GROUP> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

If the initial cleanup itself used `nstates=1`, that project already provides
the clean consensus and this section can be skipped. In all three cases, the
resulting single-state project supplies the common map and pose scaffold for
the actual heterogeneity analysis. Inspect its map, half-map agreement,
angular coverage, and particle count before continuing.

## 5. Separate conformational or binding states

Run `refine3D_states` from the clean single-state project:

```bash
simple_exec prg=refine3D_states \
  projfile=<SINGLE_STATE_ABINITIO3D_PROJECT.simple> \
  nstates=<NUMBER_OF_STATES> pose_policy=global \
  pgrp=<POINT_GROUP> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

Do not pass `vol1` through `volN` to this program. It derives same-lineage
states from the consensus project. For a state-0/1 input, the default
`flex=yes` uses `flex_pca` to initialize the states before their joint
refinement.

The pose policy controls how much orientation can change while states compete:

- `global` is the default and performs a full orientation search.
- `local` permits a limited change around each consensus direction.
- `fixed` keeps the projection direction fixed while optimizing state,
  in-plane angle, and shift.

Start with `global` unless the consensus orientations are already trusted and
the scientific question specifically calls for more constrained
classification.

**Checkpoint:** useful states should contain enough particles, retain broad
view coverage, show interpretable structural differences, and remain similar
when the analysis is repeated or the requested state count is changed. Treat
very small, poorly oriented, or irreproducible states as suspect.

## 6. Extract and refine each useful state

For every state worth retaining, run a separate selection from the original
`refine3D_states` project:

```bash
simple_exec prg=selection \
  projfile=<REFINE3D_STATES_PROJECT.simple> \
  oritype=ptcl3D state=<STATE_NUMBER> prune=yes
```

Then refine that particle group as a single structure:

```bash
simple_exec prg=refine3D_auto \
  projfile=<SELECTED_STATE_PROJECT.simple> \
  vol1=<REFINE3D_STATES_DIRECTORY/recvol_stateNN.mrc> \
  pgrp=<POINT_GROUP> mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

Replace `NN` with the zero-padded original state number: state 1 is `01`,
state 2 is `02`, and so on. Use the map produced by the same
`refine3D_states` run as the selected particles. Passing `vol1` explicitly
is important and guarantees that `refine3D_auto` starts from the matching
state map. Because this reference volume and these poses have the same SIMPLE
lineage, leave the advanced `ref_pose_init` parameter at its default `none`.

This refinement includes a methodology similar to what is usually referred to
as auto- or envelope masking. SIMPLE derives a density envelope internally and
suppresses low-density solvent outside that envelope when preparing matching
references and filtering the reconstruction; no input mask is required. This
procedure can be switched off by setting the optional parameter `automsk=no`.

Do not refine several states together. Each selected state gets its own
selection project, `refine3D_auto` run, and validation record.

## 7. Optional continuous-heterogeneity analysis

Discrete states do not always describe a flexible complex well. After a good
consensus or state-specific refinement, `flex_pca` can reveal reproducible
low-dimensional variability and reconstruct representative volumes:

```bash
simple_exec prg=flex_pca \
  projfile=<REFINED_PROJECT.simple> \
  neigs=10 npreimages=16 mskdiam=<MASK_DIAMETER_A> \
  nparts=<PARTITIONS> nthr=<THREADS>
```

The default `neigs=10` requests ten covariance components and
`npreimages=16` sets an upper limit of sixteen representative state volumes;
the recovered number can be smaller. The project consensus map is used when
`vol1` is omitted. Begin with the defaults rather than requesting many
components: excess components tend to describe fitting noise and make the
state analysis less stable.

Interpret a motion only when it is supported by reproducible components,
adequate particles and views across the trajectory, and coherent changes in
the reconstructed volumes. A smooth-looking sequence by itself is not proof
of continuous biology.

## 8. A minimal command worksheet

Fill this in before starting and update it after every accepted run:

```text
Input picked/stream project:
Mask diameter (A):
Point group (use c1 if none):
nparts:
nthr:

Accepted sieving project (or STREAM):
Accepted abinitio2D project, run 1:
Accepted abinitio2D project, run 2:
Initial abinitio3D cleanup project:
State retained from cleanup, or states merged:
Merged-state selection project, if applicable:
State-continuation or single-state consensus project:
refine3D_states project:

State | selected project | matching recvol_stateNN.mrc | refine3D_auto project
------|------------------|----------------------------|----------------------
      |                  |                            |
```

## 9. Common mistakes

- **Using the wrong project:** always copy the exact output project path from
  the accepted previous run.
- **Running offline sieving after SIMPLE stream:** this repeats a cleanup stage
  unnecessarily and can remove useful particles.
- **Treating cleanup states as final biology:** the initial three-state
  `abinitio3D` is primarily for separating useful particles from junk.
- **Repeating ab initio unnecessarily:** when retaining one cleanup state, use
  `abinitio3D state=N`; do not start over or also supply `vol1`.
- **Skipping the clean consensus:** whether continued from one state or rebuilt
  after merging states, `refine3D_states` needs a common map and pose scaffold.
- **Selecting from an already selected project:** make every independent state
  extraction from the original multi-state project.
- **Using the wrong map for a selected state:** pass the matching
  `recvol_stateNN.mrc` explicitly to `refine3D_auto`.
- **Over-interpreting tiny states or PCA components:** require particle support,
  view coverage, map quality, and repeatability.
- **Oversubscribing the machine:** total resource use can multiply across
  chunks, partitions, and threads.

## Appendix: training a data-set-specific 2D rejection model

This optional route is useful when many comparable `abinitio2D_chunks`
outputs must be screened consistently.

### A. Manually label representative projects

Review the final classes in several representative chunk projects. Mark good
classes as state 1 and rejected classes as state 0. Both good and bad examples
must be present, and the labels—not model predictions—are the training truth.

### B. Export training features

For each reviewed project, run:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=analyze quality_context=chunk \
  projfile=<MANUALLY_LABELLED_PROJECT.simple> \
  mskdiam=<MASK_DIAMETER_A>
```

This writes `cavgs_quality_training.txt` without changing the selection. Build
a plain-text SIMPLE file table containing one absolute training-file path per
line, for example:

```text
/data/run_01/cavgs_quality_training.txt
/data/run_02/cavgs_quality_training.txt
/data/run_03/cavgs_quality_training.txt
```

All files in one training set must use the same quality context.

### C. Learn the model

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=learn \
  filetab=<TRAINING_FILETAB.txt> \
  fname=<DATASET_QUALITY_MODEL.txt>
```

Keep the model together with a note identifying the data set, training
projects, manual reviewer, and date. A model trained for one specimen or
imaging regime should not silently become a universal model.

### D. Apply the learned model

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply infile=<DATASET_QUALITY_MODEL.txt> \
  projfile=<ABINITIO2D_CHUNK_PROJECT.simple> \
  mskdiam=<MASK_DIAMETER_A> prune=yes
```

Inspect the accepted and rejected class-average reports before using the
pruned output. If the model makes systematic mistakes, correct the manual
labels, expand the representative training set, retrain, and record the new
model separately.
