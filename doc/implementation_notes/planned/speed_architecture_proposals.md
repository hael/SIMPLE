# Speed and architecture: proposed next changes

Proposal note, 2026-10-18. **Nothing here is decided or implemented.** It
collects the changes proposed after the run `speed_mac`, so that they can be
ranked apart from what was done. What was done, and how it measures against
the SIMPLE standards, is in `identical_speed_run_review.md` (beside this
note) and in the run's own report, `identical_speed_run_report.md`.

Each proposal says what it is, why (the measurement or finding behind it),
the expected gain, whether results stay byte-identical, its size, and what it
needs first. Gains are estimates unless a measurement is named.

## A. Finishing the current change

Small changes to the code of the run `speed_mac`, found by the review. None
changes a result; all can be checked by the unit tests and, on the test
machine, by the one-thread identity cases. Proposed as one short follow-up
(about half an hour of work and the fast gate) on top of the run's diff, so
that test jobs started on the run's diff stay valid.

1. **Stop on a short or unreadable dense part file.** Today the read checks
   end of file with `fileiochk`, which only prints, and the master's table is
   no longer filled first, so a truncated file would leave undefined
   distances in the assignment. Stop on any read error and compare the file
   size with the size its header implies.
2. **Check the active states of a dense part file.** A row's index now means
   "reference among the active states" on the worker that wrote it; the old
   format let the master check each reference. Write the active state list
   after the header and compare it on read.
3. **Raise part-read errors after the parallel loop**, with a status per part,
   instead of stopping from inside the OpenMP region.
4. **Free volume assembly's handed-on project early.** Copy it after the
   reconstructors are freed and release it after its last use, or hand back
   only the fields assembly changes. Saves about one full project (about 1 GB
   at 750,000 particles) at the master's peak.
5. **Delete the dead code the review found**: `gen_corr_grad_for_rot_8` (its
   only caller is now a test), the per-thread buffers `kcorrs` and `pft_r1_8`
   (never used), and the release of the kept builder in `distr_initialize`
   (cannot fire); fix the stale comments and the tester header that names the
   run.
6. **Guard the "load the reference before evaluating" contract** of the new
   shift-search evaluators: keep the loaded reference's index in the thread's
   buffers and stop on a mismatch.
7. **Tidy the table types**: the dense format word private to the dense
   table; `write_rows` overridden with an error in the neighbourhood type,
   where it would overwrite the file header; 64-bit offsets for the refined
   lists.
8. **Documents**: the skills `simple-main-nu-filt` and
   `simple-nonuniform-regularization` still describe the scratch-file filter
   (an example in the first would not compile); the planned note
   `next_release/calpha_nu_evidence_mask.md` calls a deleted routine. Skill
   edits need the maintainer's approval.

## B. The commanders used as library routines

**What is there.** Three commanders gained typed entry points beside
`execute(cline)`: `commander_prob_align%execute_from_project` (and the
neighbourhood one), `commander_volassemble%execute_hand_on` and
`commander_commit_sigma2%execute_from_project`. The distributed `refine3D`
strategy calls them with in-memory projects, so that they copy the project
instead of reading the file. Each entry point and the old `execute` call one
private body in the commander module.

**Why it is odd.** In SIMPLE a commander is a command object: it is reached
through `execute(cline)`, builds its own `parameters` and builder, and owns a
named program's workflow. Called from a strategy with typed arguments, the
three commanders become library functions, and about 520 lines of volume
assembly workflow stay in a commander module while a strategy drives them.
SIMPLE already has the right shape for this need: the PCG assembly's body is
`execute_rec3D_pcg_distributed_master(params, build, ...)` in
`strategies/parallelization`, and the strategy calls it with the master's own
builder.

**Options.**

| option | what | identical? | cost |
| --- | --- | --- | --- |
| keep, then restructure (recommended) | keep the entry points for now (they are thin wrappers over private bodies); move the bodies to strategy routines in the next step of section C, "one project owner" | yes | none now; the move is part of that step |
| restructure now | move the three bodies into strategy routines that take the master's `(params, build)` and the few settings that differ today through the command line (threads, `combine_eo`, `update_frac`, state) as explicit arguments; `execute(cline)` becomes a thin wrapper that builds `params` and the builder and calls the routine | yes | 500-900 lines, mostly moved; delays the merge; test jobs started now would need repeating |
| revert the hand-offs | the three commanders read the project file again as before; keep every other change of the run | yes (the base did this) | small; gives up about 4 s per iteration at 750,000 particles (one full read costs 1.3 s, a copy 0.28 s, measured) |

**Recommendation.** Keep the entry points through the current testing and
merge, and make the restructuring the first architectural step (section C,
step 1). Restructuring removes the copies altogether: the strategy routines
work on the master's own builder, so there is nothing to copy and nothing to
hand back. If the pattern should not reach master even briefly, revert the
hand-offs instead; that is the quick, clean option, at a small cost in speed.

## C. Architectural steps for performance

Where the time goes, from the measurements so far: on PfCRT (single state,
16,832 particles, 10 partitions of 8 threads) the probability tables were 44%
of the run, the search 19%, volume assembly 13% (the NU filter 330 s), the
pre-step outside the tables 12%, and reading and preparing particles twice
per iteration about 9.5%. On msp1 (three states, 141,838 particles) the
tables were 68%. On the ribosome run the table pre-step was 46% and the
particle pipeline 3-4%. These baselines predate the runs `prob_speed`,
`ms_speed` and `speed_mac`, which have since cut much of the table time.

Ranked by expected gain on large runs.

1. **One project owner on the master.** The master builds its toolbox once
   and materialises the reprojection model on its own builder; the bodies of
   volume assembly, `prob_align` and the sigma commit become strategy or
   domain routines on `(params, build)` (this settles section B); the master
   reads only what workers write and writes the project once per iteration.
   - Gain: about 3-6 s per iteration at 750,000 particles, growing with the
     particle count; one to two fewer full projects at the master's peak;
     the kept builder's hand-made steps and the copies go.
   - Identical: yes. Size: 500-900 lines, mostly moved.
   - In the way: volume assembly's internal procedures share about 40 local
     variables; assembly's settings differ from the master's only through the
     command line, and `parameters` may not be copied; the shared-memory
     path's rule that volume assembly rereads the project file
     (`publish_before_file_assembly`, also in the strategies skill).
2. **The worker fold.** One worker process per partition and iteration fills
   its table part, signals the master, waits for the assignment, then
   searches and reconstructs with the references and particles it already
   holds.
   - Gain: the second process start, project read, reference build and
     particle preparation per iteration (about 9.5% on PfCRT), and a queue
     round per iteration on clusters.
   - Identical: achievable if the search starts from the state a fresh process
     would have. Size: medium to large.
   - In the way: the job scheduler blocks until jobs end, so it needs separate
     schedule, poll and wait steps; the table partitions and the search
     partitions must coincide; memory per worker rises (table, references and
     particles live together). Easier after step 1.
3. **Faster label optimisation in the NU filter.** The new timers show the
   label optimisation as the largest part of the filter (8.7 of 14.9 s on the
   three-state identity case). Reuse the setup's smoothing at each level
   (drops 7 of 35 smoothings at 8 labels), skip all-zero lines and block the
   strided passes of the tent smoothing, and give the Potts sweeps packed
   voxel lists per colour.
   - Gain: likely 2-3 times on the optimisation; measure its share on PfCRT
     first with the new split timers.
   - Identical: yes. Size: about 150 lines.
4. **A transposed master table and contiguous reads.** Store the master's
   dense table particle-major, so the assignment streams instead of striding
   by the number of references; read runs of consecutive particles in one
   statement; let workers write per-particle refined counts so the master
   places refined entries directly and in parallel.
   - Gain: the assignment and the reads on large multi-state tables; measure
     the assignment time on msp1 first.
   - Identical: yes. Size: about 150 lines.
5. **Batched, GPU-ready kernels.** Pure kernels over explicit arrays
   (references, particle times CTF, shifts, rotations) that score many
   references per call; on the CPU first by batching references per particle.
   - Gain: the score loop is the largest single item in the tables, and this
     interface is what the GPU release needs.
   - Identical: probably not (FFT plans and precision change). Size: large.
6. **The compact record for the 2D table.** The same record, parallel read
   and coverage check as the 3D dense table.
   - Gain: `solve2D` table reads and memory. Identical: yes. Size: 300-400
     lines.
7. **Shift-independent work once per search.** Compute the correlation's
   shift-independent term once per minimisation; build the shift matrix
   inside the kernel loops instead of writing it and reading it back.
   - Identical: yes for the correlation objective, rounding-level for the
     Euclidean one. Size: small to medium.
8. **The last two heap trampolines**, in the bounded and power samplers (2D
   sampling, `greedy_sub`), with the array sampler of the run.
   - Identical: yes. Size: about 120 lines.
9. **The NU filter as an object** instead of module state, so that states
   could be filtered concurrently.
   - Identical: no (sums whose rounding depends on the thread count; FFT
     plans). Size: 1-2 weeks. Needs a ruling on the multi-state memory
     policy, which forbids holding several states' volumes.
10. **Persistent workers** that live for the run and exchange only model,
    sample and assignment.
    - Gain: perhaps 5-15% per iteration on queued clusters, less on one
      machine. Identical: possible, risky. Size: large; after steps 1 and 2.
11. **A streaming assignment** that keeps only the best candidates per
    reference from each part.
    - Gain: the master's table memory by orders of magnitude on large
      multi-state runs. Identical: no; needs a ruling. Size: large.
12. **Single-precision weights in the score kernels.**
    - Gain: twice the vector width in the score loop. Identical: no (rounding
      level); judged by the truth suite. Size: small.

## D. Measurements that would settle the ranking

- **PfCRT and msp1 on the merged code**, with the new split timers of volume
  assembly: how much of the run is now tables, search, NU filter (and its
  label optimisation), the master's project handling and the assignment.
- **The master's peak memory** in a distributed `refine3D` at a large
  particle count, to size section A item 4 and step 1.
- **Identity on the x86 build** of the test machine (the run's checks were
  on the Mac's ARM build).

## E. Suggested order

1. Section A, now, while the first test jobs run on the run's diff.
2. Section B decided: keep and restructure in step 1 (recommended), or
   revert the hand-offs before the merge.
3. After the measurements of section D: steps 1 and 3 (both identical; step 1
   also settles section B), then step 2.
4. Steps 4, 6, 7 and 8 as small identical items whenever convenient.
5. Steps 5, 9, 10, 11 and 12 as the GPU-release and results-changing work,
   each with its own ruling and truth-suite validation.
