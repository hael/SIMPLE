# Accelerating refine=cont: fewer evaluations, a cheaper kernel, and a benchmark

Implementation note, 2026-10-19, completed 2026-10-21. **All seven stages are
done: built, checked, and measured against the baseline. The measured gains
are in section 2.** The maintainer's rulings are in section 8. The findings come from a static review of `master` at
`4da5cb9c0`. The savings below are work counts read from the code. The rule of
the local-workflow plan applies: do not infer a speedup percentage from removed
calls (`doc/refactoring_notes/planned/refine3d_local_workflow.md:264-267`).
Every stage is therefore measured with the benchmark of section 6.

This is also a learning note for Maz. Each stage says why it comes where it
does, what to run and what to look for, so the work can be followed and
repeated without the conversation it came from.

How to read this note:

- Section 2 is the plan, as seven stages. Stage 2 ends with the baseline
  record, a table to fill in.
- Sections 3 and 4 hold the findings the stages implement.
- Sections 5 to 7 hold the measurement and test detail.

Companion notes:

- `cont_kb_gather_performance.md`: the fused gather and its list of
  future work.
- `doc/refactoring_notes/planned/refine3d_local_workflow.md`: the
  workflow-level startup savings. This note covers the per-particle cost inside
  `strategy3D_cont`.
- `identical_speed_run*.md` and `single_state_speed_week*.md`: the polar speed
  work this note compares against.

The engine is shared, so every item also applies to the polish
(`pose_cont=yes`) and to the proposed `refine3D_local`.

## 1. The route and where it stands

One particle runs one transaction:

```text
strategy3D_cont%srch_cont
  -> cartft_pose_opt%refine -> refine_pose -> refine_joint            (LM loop)
       -> cartft_calc%pose_normal_terms     objective, 5-gradient, 5x5 Gauss-Newton matrix
            -> gather_packed_bspline3_window_grad 27-sample B-spline gather and 3 derivatives, per pixel
  -> euclid_sigma2%calc_sigma2 -> cartft_calc%sigma_contribution     (euclid only)
```

The search time of a particle is the number of full-disk passes it makes,
times the cost of one pass. Section 3 reduces the first factor and section 4
the second.

Done for this route so far:

- **Fused 3x3x3 KB gather.** Median worker orientation-search time -36.8 to
  -37.5%, wall time -20.5 to -21.9%, on beta-gal
  (`cont_kb_gather_performance.md:142-148`).
- **Joint stage alone.** 37% fewer LM iterations and 38% less transaction time
  than `shift_then_joint`
  (`doc/refactoring_notes/completed/pose_cont_refactoring.md:1615-1618`).
- **Dynamic scheduling of the matcher's particle loop.** Commit `aea5b1846`,
  `simple_strategy3D_matcher.f90:235`.
- **A continuous gather (2026-10-20, after stage 5).** The quadratic B-spline
  replaced the KB stencil in the gather (record after stage 5).

None of the speed-week kernel changes reached `src/main/cftc`. The speed-run
documents name `refine=cont` only as the identity case CR1. These polar items
have a cont counterpart:

| Speed-week item (polar) | Cont counterpart | Section |
| --- | --- | --- |
| Reuse of the last accepted objective (`simple_pftc_shsrch_grad.f90:59-64`) | Redundant seed, final and sigma passes | 3.2 |
| Dynamic OpenMP schedule of the per-particle loops | Done at `4da5cb9c0` | - |
| Half the in-plane rotations stored (`nrots = 2*pftsz`, `simple_polarft_core.f90:60-61`) | Full-disk pixel loops | 4.1 |
| Evaluation counters and thread-seconds in the prob-table BENCH records | None: the stage counts never leave the optimizer | 5 |
| Replaced routine kept as an oracle (`reference_assign`, the pftc in-plane tester) | No oracle for the 5x5 matrix | 7 |

**Cost evidence.** The only per-particle record is from beta-gal: 5,513
particles, box 256, `nthr=12`, an NCC objective, and the build before the fused
gather. It shows 136,968 LM iterations, about 25 per particle, with a mean
transaction time of 0.217 s and a maximum of 0.562 s
(`pose_cont_refactoring.md:1597-1604`). Every LM iteration that evaluates makes
one full-disk pass of the gather.

## 2. Development stages

The agent builds each stage, runs the checks it names and reports the results.
The maintainer reviews and decides what is committed and pushed, and when.
Stages are not gated on commits; several may go into one larger commit.

**What is measured (ruled 2026-10-20).** The working tree, with every
uncommitted change in it, is the code under development, so the builds and
measurements use the full tree. That includes changes outside this work, such
as PCG reconstruction work in progress. For an A/B pair the tree is
snapshotted once (`git diff HEAD`, kept as
`~/src/SIMPLE_cont_bench/patches/tree_snapshot_<time>.patch`) and both
variants are built from that one snapshot. A pair then differs only in the
change being measured. The stage 2 baseline (builds A and B) predates this
rule and stays the reference of the stage 2 and stage 3 records. The order is: measure first, then protect the kernel with tests,
then make the changes that keep results identical, and only then the changes
that alter them. Every gain can then be attributed to one stage.

| Stage | What | Results | Check |
| --- | --- | --- | --- |
| 1 | Count the work | Identical | CR1 identity check (section 6.2) |
| 2 | Benchmark workflow and baseline | No change in production code | Baseline recorded in this note |
| 3 | Remove work without changing results | Identical | CR1 identity check; benchmark |
| 4 | Converge in fewer iterations | Path and endpoints change | Robustness profile on truth; benchmark |
| 5 | Kernel oracle and test gaps | Tests only | New tests pass on the current kernel |
| 6 | Cheaper evaluations | Rounding changes | Stage 5 oracle; truth profile; benchmark |
| 7 | Close-out | Documents only | - |

**Measured gains** (orientation search at 16 threads, wall seconds; each row
from the build its record names, so rows from different batches hold to about
5%):

| After | Beta-gal search | Passes per solve | ns per pixel-evaluation | Polish euclid search | Truth |
| --- | --- | --- | --- | --- | --- |
| Stage 2 baseline (B) | 19.3 | 13.75 | 215 | 24.1 | the before-profile |
| Stage 3 (D): duplicate passes removed, `padsc` folded | 16.1 | 11.75 | 210 | 20.6 | identical |
| Stage 4 (S4c): curvature scale, predicted-reduction stop | 9.4 | 6.68 | 215 | 11.7 | inside the before-profile's spread |
| Continuous gather (G) | 5.4 | 5.70 | 143 | 6.1 | rotation error down in every polish seed |
| Stage 6a (H): half plane and normal-term arithmetic | 2.40 | 5.70 | 64 | 2.68 | unchanged |
| Stage 6c (K): sigma2 from the accepted evaluation | 2.01 | 4.70 | 65 | 2.28 | unchanged |
| Stage 6b (L): separable gather | 1.90 | 4.70 | 62 | 2.08 | unchanged; gate 20 of 20 |

Over the whole plan, the beta-gal search fell from 19.3 to 1.9 s (-90%) and
the simulated polish from 24 to 2.1 s (-91%), with the same map and truth
metrics inside the before-profile's spread. Stages 1, 2, 5 and 7 change no
production result: counters, benchmark, tests and documents.

### Stage 1: count the work

**Why.** Today the bench record reports one number for the whole Cartesian
search, `match3D orientation search`
(`simple_strategy3D_matcher.f90:373`). A later change can speed the search up
in two ways: fewer passes per particle (section 3), or cheaper passes
(section 4). A single time cannot tell the two apart, and the policy forbids
inferring a speedup from removed calls. The first change therefore makes the
code report how much work it does.

- **What:** the counters of section 5, added to the REFINE3D bench record.
  These include passes per particle, LM iterations, ns per pixel-evaluation and
  the actual/predicted ratio that diagnoses section 3.1.
- **Results:** identical. Only BENCH content changes, and identity excludes
  BENCH files.
- **Check:** the CR1 identity check (section 6.2), run in step 2.2.
- **Done when:** one run reports these figures for each partition.

**Status (2026-10-19): done.** It was committed in `0c9542e50` and checked in
step 2.2: the CR1 identity holds, and the counters cost nothing measurable.
The changes:

- `simple_cartft_calc.f90`: each slot counts the pixels of its shells
  (`get_ptcl_npix`), and `get_ref_bytes` gives the bytes of the reference
  transforms.
- `simple_cartft_pose_opt.f90`:
  - the solve counts its calculator evaluations (`npasses`);
  - each stage counts its evaluated trials (`ntrials`, beside `nattempted`,
    which includes the joint stage's bound rejections) and sums actual over
    predicted at accepted steps;
  - `cartft_pose_work` sums solves and writes the bench section.
- `simple_strategy3D_cont.f90`: `add_work` adds the particle's solve, with its
  slot's pixel count and the sigma2 residual pass.
- `simple_strategy3D_matcher.f90`:
  - one work record and the busy time per thread;
  - a `THREAD-SECONDS` block for every 3D search, following `match2D`;
  - the `CONTINUOUS POSE WORK` section for a Cartesian pass;
  - the reference bytes in the context;
  - the Cartesian reference read moved from startup to reference preparation.
- `refine3D_policy.md` section 11.1 describes the new sections.

### Stage 2: benchmark workflow and baseline

**Why the baseline comes first.** An optimization is only a claim until it has
been compared with what came before: the same work, on the same machine, with
the same build settings. The baseline is that "before". It is measured on the
code as it is now, with the stage 1 counters in place and no optimization, and
every later stage is judged against it. It has to be taken before stage 3,
because once a change lands the tree no longer builds the unoptimized code.

A baseline is more than one row of timings. It has four parts, each answering
its own question:

| Part | The question it answers | Steps |
| --- | --- | --- |
| The baseline executable | What do we rerun next to each new stage? | 2.1 |
| The identity and overhead checks | Did stage 1 change a result, or slow the search? | 2.2, 2.3 |
| The benchmark figures | What does a particle cost today, in work and in time? | 2.4, 2.5 |
| The before-profile on truth | How good are the poses today, and how much do they vary from seed to seed? | 2.6 |

**What was implemented (2026-10-19; built and checked, commits `b2d8f34e5` and `920b0e410`):**

- `simple_test_exec test=cont_refine3D_bench`, the benchmark program of
  section 6.3:
  - the code is in `simple_commanders_test_highlevel_cont.f90`;
  - it is registered in `simple_test_ui_highlevel.f90` and
    `simple_test_exec_highlevel.f90`;
  - it is listed in section 1.3 of `test_environment_policy.md`.

  `scripts/check_test_registry.py` passes with it.
- `scripts/plot_refine3d_bench.py` summarizes the stage 1 sections.
  `bench_report.txt` gains a "Thread-seconds (3D search)" block and a
  "Continuous pose work" block, and `bench_records.csv` gains their values. The
  totals and ratios are recomputed from the sums over the parts. This was
  checked on a synthetic record.

#### 2.1 Keep two executables

- **Build A:** `86bf97213`, the commit before stage 1.
- **Build B:** the stage 2 commit.
  - Stage 1 is in `0c9542e50`.
  - Stage 2 adds the benchmark program and the plot-script summary and changes
    no production code, so B is still the unoptimized code with the counters.
  - The benchmark program exists only in B, which is why B, not `0c9542e50`,
    is the baseline.

`0c9542e50` also carries solve3D stage-progression changes. They act only
through `lpexit`, which only the solve3D controller sets, so they reach neither
`refine=cont` nor the CR1 cases.

Build both in Release with the usual compile script, each into an install
directory of its own, and write each commit hash beside it. Leave B unchanged
until the work is finished.

- **Why B:** every later stage is compared with B in pairs run on the same
  machine, not with numbers recorded weeks earlier. Hardware state, other jobs
  and file caches change over weeks; two runs made back to back share them.
- **Why A:** stage 1 is meant to change no result and cost no time, and only A
  can show that.

#### 2.2 Check that stage 1 changed no result

1. With B, run the unit suite of the area:

   ```bash
   simple_test_exec test=unit_cart_align3D
   ```

2. Run the CR1 identity check (section 6.2). Run the gate at one thread with A
   and with B, each in an empty directory of its own. This takes about 35
   minutes and 9 GB of memory per build.

   ```bash
   simple_test_exec test=cont_refine3D_1jxy nthr=1
   ```

   Then compare the two `metrics.tsv` files the runs leave in their
   directories:

   ```bash
   diff run_A/metrics.tsv run_B/metrics.tsv
   ```

   Every row must be equal.

**Why one thread.** With several threads, floating-point sums are added in an
order that changes from run to run, so two runs of the same build can differ in
the last digits. At one thread a run is reproducible, so any difference comes
from the code.

#### 2.3 Check what the counters cost

Stage 1 adds a timer call per particle and a few integer additions per LM step.
That should be invisible, but it has to be shown, or the baseline already
contains a slowdown.

**What to run.** One order-balanced pair on the beta-gal checkpoint (section
6.3): A, B, B, A, each run in a fresh copy of the checkpoint.

**What to compare.** The median worker orientation-search time,
`match3D orientation search` in `BENCH/REFINE3D_BENCH_PARTppp.txt`. The
difference between A and B must lie within the spread of the two A runs.

**Why order-balanced.** Each build runs once early and once late. A machine
that warms up or slows down during the session then favours neither build.

#### 2.4 Cover the realistic SNR range

**The range.** Cryo-EM particles have an SNR of 0.02 or better. The default,
`snr=0.05`, is a typical value. Step 2.5 also runs the polish regime at 0.02,
the noisiest realistic case. The convergence problem of section 3.1 grows as
the particle correlation falls, so this is where it is largest.

**A trap to avoid.** Do not tune `snr` to reproduce a correlation measured with
a different objective. The beta-gal record's 1 - NCC of 0.83 (correlation
0.17) was measured on sigma-whitened data, where the noisy high shells weigh
as much as the low ones. Today's `objfun=cc` weights every pixel equally, so on
the same data its correlation is much higher: 0.79 at `snr=0.05` in the smoke
runs below. Chasing 0.17 with the cc objective pushed `snr` to 0.0005, 40
times noisier than any real data set. That mistake was made, and corrected,
while this note was written.

The program refuses to overwrite `cont_refine3D_bench/`, so move that
directory away between runs.

#### 2.5 Run the benchmark with B

**What to run.**

```bash
simple_test_exec test=cont_refine3D_bench nthr=8
```

Replace `nthr=8` with the machine's core count. Then, after moving
`cont_refine3D_bench/` away, run the polish regime at the noisiest realistic
SNR:

```bash
simple_test_exec test=cont_refine3D_bench suite=polish snr=0.02 nthr=8
```

**The cases.** It runs four cases: `polish_euclid`, `polish_cc`,
`capture_euclid` and `capture_cc`. Each is one `refine3D refine=cont` pass in a
directory of its own.

**Summarize each case.** For example:

```bash
python3 scripts/plot_refine3d_bench.py cont_refine3D_bench/polish_euclid
```

**What to record.** Take these from `bench_summary/bench_report.txt` and
`metrics.tsv` for each case:

- search thread-seconds, search wall time, peak RSS and Cartesian reference
  bytes;
- passes per solve, LM iterations per solve, the solves at the iteration cap,
  the rejected fraction of trial evaluations and the actual/predicted
  histogram;
- ns per pixel-evaluation;
- median rotation and shift error after the pass, the fraction in the basin and
  the mean `corr_cart`.

**Beta-gal.** Run the checkpoint with B three times, from fresh copies, and
record the same figures. Real data has no truth, so there are no pose errors.

**How to read the first numbers:**

- **Redundant passes (section 3.2).** Let N be the trial evaluations per solve.
  Passes per solve near N+4 (euclid) or N+3 (cc) confirm them.
- **The CC diagnosis (section 3.1).** cc's accepted steps centred near
  `2 - cc` on the actual/predicted histogram confirm it. At realistic SNR
  that is 1.2-1.35, a mild effect.
- **The open euclid question (section 8).** The euclid ratio and rejected
  fraction decide whether the euclid matrix needs scaling too.

#### 2.6 Take the before-profile on truth

**Why a profile, not a run.** Stages 4 and 6 change results. A single run on
simulated data is one draw, and its error moves from seed to seed. These stages
are therefore judged by the distribution of the pose errors over many seeds,
before and after. The "before" distribution is taken now, with B.

**The gate.**

```bash
simple_test_exec test=profile target=cont_refine3D_1jxy nseeds=20
```

It runs the gate at 20 seed offsets. For every metric it writes the median, the
spread, the worst value and the failure rate.

**The benchmark simulation.** Run the polish regime at five more seed offsets,
k = 1 to 5, each in its own directory. Offset 0 is the run of step 2.5.

```bash
simple_test_exec test=cont_refine3D_bench suite=polish seed_offset=1 nthr=8
```

Keep all the profile directories together with build B.

**Done when:** the baseline record below is filled in.

#### Baseline record

Taken 2026-10-19 by the agent with `~/src/SIMPLE_cont_bench/run_baseline.sh`.
Builds A and B are detached worktrees there, each with its own `build/`, and
every run directory is under `runs/`. The beta-gal checkpoint is a read-only
copy of `~/Processing/bgal/bench/5_solve3D`: `bench.simple`, whose stack path
is absolute, and `recvol_state01{,_even,_odd}.mrc`. Each beta-gal run starts
from a fresh copy of it.

| Item | Value |
| --- | --- |
| Build A commit | `86bf97213` |
| Build B commit (stage 2) | `920b0e410`: stage 2 with the metric fix. The timed beta-gal runs and the CR1 identity check ran on `b2d8f34e5`, whose solver is the same (`920b0e410` changes only the benchmark's metric and a `solve3D` `nu_stage` check). Work counts on the simulated cases are identical on both. |
| Machine, threads | FWM-S116118, Apple M2 Ultra (16 performance + 8 efficiency cores, 192 GB); timed runs `nthr=16`, alone on the machine |
| SNR | 0.05; polish also at 0.02 |
| `unit_cart_align3D` (B) | PASS, 761 of 761 assertions; fast gate 15/15 for A and for B |
| CR1 identity, A against B | Identical: the 60 check values of `metrics.tsv` are byte-identical, and both gates pass (`runs/cr1_A`, `runs/cr1_B`, about 40 minutes each at one thread) |
| Counter cost (beta-gal, A B B A) | Orientation search A 19.34, 19.41 s; B 19.33, 19.46 s (a third B run 19.12 s). The difference is inside the spread: no measurable cost. |
| Before-profile directories | `runs/profile_cr1_B/profile`: 20 of 20 seeds pass (table below). `runs/bench_default` and `runs/bench_polish_seed{1..5}`: six seed offsets of the polish regime (table below). The superseded runs with the coarse metric are in `runs/superseded_metric`. |

Simulated cases: 10 000 particles, box 256 at 0.73 A, reference reconstructed at
the truth poses. The polish band is shell 62 (lp 3 A) and the capture band
shell 23 (lp 8 A). Beta-gal: 4 199 active particles of 5 513, box 256 at
1.275 A, band shell 88 (3.7 A), D2, euclid.

| Case | Search thread-s | Search wall (s) | Passes / solve | LM its / solve | Solves at the cap | Rejected fraction | Actual/predicted mode bin (share) | ns / pixel-eval | Median rot (deg) | Median shift (A) | Mean corr_cart |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| polish_euclid, snr 0.05 | 380.3 | 24.1 | 14.88 | 10.88 | 0 | 0.41 | 1.25-1.75 (0.32) | 208 | 0.0845 | 0.029 | 0.67 (`exp(-L)`) |
| polish_cc, snr 0.05 | 384.0 | 24.3 | 14.43 | 11.43 | 0 | 0.38 | 1.25-1.75 (0.38) | 217 | 0.0864 | 0.030 | 0.78 |
| capture_euclid, snr 0.05 | 67.5 | 4.3 | 17.69 | 19.70 | 33 | 0.16 | 0.75-1.25 (0.47) | 220 | 0.250 | 0.095 | 0.83 (`exp(-L)`) |
| capture_cc, snr 0.05 | 61.5 | 3.9 | 15.71 | 19.91 | 15 | 0.16 | 0.75-1.25 (0.64) | 226 | 0.321 | 0.133 | 0.89 |
| polish_euclid, snr 0.02 | 376.0 | 23.8 | 14.78 | 10.78 | 0 | 0.42 | 1.25-1.75 (0.31) | 207 | 0.127 | 0.045 | 0.56 (`exp(-L)`) |
| polish_cc, snr 0.02 | 429.3 | 27.1 | 16.40 | 13.40 | 0 | 0.36 | 1.25-1.75 (0.51) | 213 | 0.130 | 0.046 | 0.65 |
| beta-gal, B (median of 3) | 305.2 | 19.3 | 13.75 | 9.75 | 0 | 0.48 | 0.25-0.75 (0.64) | 215 | - | - | - |

Seeds were 2 degrees and 1 A off for polish, and 15 degrees and 2.6 A off for
capture. Every solve ended "accepted". Memory: the Cartesian references take
1.08 GB (two halves at 0.54 GB); peak RSS is 6.9-7.8 GiB for the simulated
cases and 5.1 GiB for beta-gal. The full histograms are in each case's
`BENCH/bench_summary/bench_report.txt`.

**Before-profile of the gate** (`test=profile target=cont_refine3D_1jxy nseeds=20`,
8 threads each, build B). Each value is the median over 20 seeds, with its
spread, the worst seed and the floor. No seed failed.

| Metric | Median | Spread | Worst | Floor |
| --- | --- | --- | --- | --- |
| a (euclid) median rotation error (deg) | 0.943 | 0.035 | 1.003 | 1.5 |
| a (euclid) median shift error (px) | 0.071 | 0.002 | 0.075 | 0.2 |
| b (cc) median rotation error (deg) | 0.229 | 0.003 | 0.237 | 1.5 |
| b (cc) median shift error (px) | 0.059 | 0.001 | 0.061 | 0.2 |
| c (distributed) difference to a (deg) | 0 | 0 | 0 | 0.01 |
| e (auto, `pose_cont=only`) median rotation error (deg) | 0.980 | 0.306 | 1.709 | 7.5 |
| f (polish) symmetry-aware error (deg) | 0.423 | 0.019 | 0.474 | 0.93 |
| g (auto with polish) symmetry-aware error (deg) | 0.176 | 0.022 | 0.564 | 1.55 |
| h (continuation) symmetry-aware error (deg) | 0.174 | 0.030 | 0.258 | 1.55 |

**The pose-error metric of the benchmark was too coarse.** It took the `acos`
of a single-precision trace, the gate's formula. Near zero that resolves only
about 0.002 degrees, so at the polish errors (about 0.085 degrees) the median
moved in steps of a few percent: five seed offsets gave 0.0839 degrees to four
digits. The benchmark now uses the chord form in double precision,
`2 asin(||R1 - R2||_F / (2 sqrt 2))`, which resolves the seeds (0.0858 against
0.0873 for euclid). The simulated cases and the seed runs were therefore
retaken with build B at `920b0e410`, the commit of this fix; the case table
above holds the retaken values.

- **Work counts:** identical to the first run, as expected with the solver
  unchanged.
- **Search thread-seconds:** within 1.4% of the first run for every case.
- **One outlier:** the polish euclid run of seed offset 3 took 435.7
  thread-seconds against about 380 for the other offsets, with the same work
  (passes per solve 14.82-14.89 for every offset). That is a reminder that
  single timings are not comparable, and that later stages are timed in
  order-balanced pairs.
- **Not affected:** the beta-gal runs and the CR1 identity check do not use
  this metric.

**Polish regime over six seed offsets** (offset 0 is `runs/bench_default`;
build B, `snr` 0.05, 10 000 particles). This is the before-profile that stages
4 and 6 compare against.

| Metric | Median | Min | Max |
| --- | --- | --- | --- |
| euclid median rotation error (deg) | 0.0844 | 0.0834 | 0.0854 |
| euclid median shift error (A) | 0.0292 | 0.0288 | 0.0293 |
| cc median rotation error (deg) | 0.0861 | 0.0850 | 0.0873 |
| cc median shift error (A) | 0.0299 | 0.0297 | 0.0301 |
| euclid mean `corr_cart` (`exp(-L)`) | 0.6722 | 0.6717 | 0.6726 |
| cc mean `corr_cart` | 0.7789 | 0.7786 | 0.7792 |

What the baseline says:

- **The fixed passes are 21-28% of all passes.** They are 3 or 4 of 13.75-17.7
  per solve (section 3.2, stage 3a).
- **36-48% of trial evaluations are rejected** in the polish cases and on
  beta-gal. That is the terminal phase of the stopping rule (section 3.3), the
  largest single waste, and stage 4's main lever.
- **Real data answers the open euclid question.** On beta-gal, 64% of euclid's
  accepted steps have actual/predicted below 0.75: its Gauss-Newton matrix is
  too small, the over-step that section 3.1 predicts for a filtered,
  attenuated reference. The simulation, whose reference is an unfiltered
  reconstruction, does not show it, so stage 4 is judged on beta-gal as well.
- **cc's steps are about 1.3-1.5 times too short at realistic SNR** (mode bin
  1.25-1.75). The scaling of section 3.1 is a modest gain.
- **The capture regime iterates longer** (about 20 LM iterations, some of them
  without an evaluation because the bound guard rejected the proposal), and a
  few solves reach the cap. Its time is small next to the polish band's, since
  shell 23 has a seventh of the pixels of shell 62.
- **About 207-225 ns per pixel-evaluation on every case,** the figure stage 6
  must lower.
- **Work counts are identical across repeated beta-gal runs** (13.75 passes
  per solve in all three), so a change in them between builds is a real
  change, not noise.

#### First observations (smoke runs, 2026-10-19)

These are not the baseline. The runs used the polish regime, 500 particles and
8 threads, with the reference reconstructed at the truth poses, on build
`0c9542e50` plus the uncommitted stage 2. They checked that the program and the
counters work.

| `snr` | Case | Mean `corr_cart` | LM iterations / solve | Passes / solve | Rejected trials | At the cap | Accepted steps by actual/predicted bin (0.25-0.75 / 0.75-1.25 / 1.25-1.75 / 1.75-2.5 / >= 2.5) | ns / pixel-eval |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.05 | euclid | 0.68 (`exp(-L)`) | 10.1 | 14.1 | 0.45 | 0 | 0.12 / 0.37 / 0.18 / 0.15 / 0.17 | 193 |
| 0.05 | cc | 0.78 | 9.7 | 12.7 | 0.42 | 0 | 0.05 / 0.34 / 0.28 / 0.17 / 0.16 | 195 |
| 0.02 | euclid | 0.56 (`exp(-L)`) | 12.6 | 16.6 | 0.44 | 0 | 0.42 / 0.25 / 0.13 / 0.09 / 0.10 | 192 |
| 0.02 | cc | 0.65 | 9.7 | 12.7 | 0.40 | 0 | 0.05 / 0.33 / 0.29 / 0.18 / 0.16 | 194 |

What they show:

- **The pass count matches section 3.2 exactly.** Passes per solve are LM
  iterations + 3 for cc and + 4 for euclid (the sigma2 pass). Every trial is
  evaluated: no proposal was cut by a bound.
- **The fixed passes are a quarter of the work.** They are 24% of cc's passes
  and 28% of euclid's at these iteration counts, so stage 3a is worth more
  than its size suggests.
- **About 40-45% of trial evaluations are rejected, in every case.** These
  are mostly the last steps of a solve, where the 1e-6 tolerance meets
  rounding noise and up to eight consecutive rejections end it (section 3.3).
  At realistic SNR this is the largest single waste.
- **The CC under-stepping of section 3.1 is mild at realistic SNR.** With a cc
  of 0.65-0.78, the theory predicts actual/predicted near `2 - cc`, that is
  1.2-1.35, and the histograms centre there: steps about 1.3 times too short,
  not 1/cc.
- **Euclid over-steps at `snr` 0.02.** Its histogram peaks below 0.75, so its
  matrix is too small, as section 3.1 predicts for a noisy, attenuated
  reference. A reconstruction from 500 particles is far noisier than one from
  10 000, so the baseline decides the open euclid question of section 8.
- **The cost per pixel-evaluation is about 192-195 ns whatever the SNR or the
  objective.** Sections 3 and 4 act on independent factors.

Two diagnostic runs outside the realistic range are kept only for what they
taught.

- **Against the noise-free truth map, at `snr` 0.05 to 0.0003** (200
  particles, cc), iterations grew from 13.3 to 39.4 as cc fell from 0.79 to
  0.12. Most solves at the lowest SNR reached the 40-iteration cap and still
  ended "accepted". This is the mechanism of section 3.1. It led to two
  counter changes:
  - the cap count;
  - the actual/predicted histogram. Steps with a near-zero predicted reduction
    gave ratios in the hundreds and made a mean meaningless (18.8 with a
    spread of 241).
- **Euclid against the truth map put 77% of its accepted steps in the 1.75-2.5
  bin.** This was an artefact of the map's arbitrary amplitude scale, which is
  why the reference is now a reconstruction.

### Stage 3: remove work without changing results

**Status (2026-10-20): done, shipped as variant D of the record below.**

- **Duplicate passes removed** (section 3.2, rows 1 and 2). The transaction
  evaluates the seed once with `pose_normal_terms` and hands the terms to
  `refine_joint`, which takes them instead of evaluating again
  (`seed_objective`, `seed_gradient`, `seed_hessian`). `refine_joint` returns
  the objective at its endpoint (`objective_out`) instead of having it
  re-evaluated.
- **`padsc` folded into the reference** (section 4.4). `set_ref` stores
  `padf**3` times the transform, and the gather no longer multiplies.
- **Results:** identical. A reused evaluation is the same computation on the
  same inputs, and scaling by 8 is exact.
- **Check:**
  - per-particle poses identical to B;
  - the CR1 identity check (section 6.2);
  - `unit_cart_align3D` and the fast gate;
  - the benchmark: two fewer passes per solve.

The planned 3b, a full centred band-limited reference with a branch-free
gather, was built and measured too (variant C below). On its own it is slower
at production bands, so it moved into stage 6a, to be measured together with
the Friedel half-plane. Its patch, against `920b0e410`, is kept at
`~/src/SIMPLE_cont_bench/patches/stage3b_centred_cube.patch`. It includes the
bit-identity oracle of the centred gather and the in-band property test.

#### Stage 3 record (2026-10-19)

Two variants were built on build B's commit `920b0e410` and measured against B,
in order-balanced pairs (B, X, X, B) at 16 threads, alone on the machine:

- **C:** 3a and 3b as planned. The duplicate seed and final passes are gone,
  `padsc` is folded in, and the reference is the full, centred, band-limited
  cube with a branch-free gather.
- **D:** 3a and the `padsc` fold only, on today's packed layout.

**Identity.** Both are bit-identical to B.

- The per-particle pose fields (Euler angles, shifts, `corr_cart`, status) of
  the project are identical for all 5 513 beta-gal records and all 10 000
  records of every simulated case.
- Two B runs are identical to each other too, which shows the runs are
  deterministic.
- The comparison flags the 4 199 refined particles when a refined project is
  compared with its input, so it is sensitive.
- The CR1 check values are byte-identical (C).
- The area suite passes, 773 of 773, including two new checks:
  - the centred gather equals the packed gather, scale aside, bit for bit
    (`simple_cartesian_fourier_tester`);
  - a band-limited reference predicts the same bits inside its band as a
    Nyquist one (`simple_cartft_calc_tester`).

**Work.** Both remove exactly two passes per solve: 13.75 to 11.75 on
beta-gal, 14.88 to 12.88 (euclid) and 14.43 to 12.43 (cc) on polish.

**Time** (orientation search, s; ns per pixel-evaluation):

| Case | B | C (cube) | D (packed) |
| --- | --- | --- | --- |
| beta-gal, band 88 | 19.0-19.4; 212-216 | 21.6; 281 (+12%) | 16.1; 210 (-15.4%) |
| polish euclid, band 62 | 23.9; 207 | 24.5; 246 (no gain) | 20.6; 206 (-13.8%) |
| polish cc, band 62 | 23.8-24.4; 212-218 | 24.3; 252 (no gain) | 20.5; 213 (-13.6%) |
| capture, band 23 | 3.9-4.3; 220-226 | 3.4-3.8; 219-227 (-12%) | not run |

**Why the cube is slower.** The cube costs 15-31% more per pixel-evaluation,
and the penalty grows with the cube: none at band 23, where it fits in cache.
A central section is centrosymmetric. In the packed half-complex layout, a
pixel and its Friedel mate therefore fetch the same memory. The full cube
stores both halves separately, which doubles the working set of a full-disk
pass. The branch-free fetch saves less than that costs.

**What follows for the order of work.** With the Friedel half-plane of stage
6a, a pass visits only half the disk. Both layouts then touch the same
footprint, and the cube's branch-free fetch should win. The cube also needs
less memory: 0.73 GB against 1.08 GB on beta-gal, 0.26 GB against 1.08 GB on
polish.

**Decided 2026-10-20:** stage 3 ships as D, and the cube moves into stage 6a,
measured together with the half-plane (section 8).

### Stage 4: converge in fewer iterations

**Status (2026-10-20): done on the working tree, shipped as variant S4c of the
record below; the after-profile is within the before-profile's spread.** These changes alter results and were ruled wanted for release 4, for
`refine=cont` and the polish alike.

- **What:**
  - scale the CC Gauss-Newton matrix by the correlation (section 3.1);
  - the stopping rule (section 3.3).

  The euclid matrix changes only if the baseline shows its over-step
  signature (section 3.1). At realistic SNR the stopping rule is the larger
  lever: in the smoke runs, 40-45% of trial evaluations were rejected.
- **What shipped** (variant S4c of the record):
  - **The curvature scale.** Both objectives' Gauss-Newton matrices are
    multiplied by the regression slope of the data on the model, weighted by
    `|f|^2` as the matrices are, and clamped to [0.05, 20]
    (`curvature_scale` in `simple_cartft_calc.f90`). For cc it is the slope
    of the normalized data on the normalized model, which equals the
    correlation `c` when every shell has the same slope. For euclid it is
    `Re<X,M>/|M|^2` with the same weights. This corrects euclid's over-step
    on beta-gal as well, so section 3.1's euclid question is settled by the
    same change.
  - **The predicted-reduction stop.** A proposal whose predicted reduction is
    below `1e-6` of the objective ends the stage without being evaluated
    (`LM_PREDICTED_RELATIVE_TOL`). The accepted-reduction stop stays as it
    was, and both tolerances stay at `1e-6`.
  - **The benchmark's capture regime** bounds the rotation at 20 degrees
    instead of 15 (record, finding 3).
- **Tests added or changed:**
  - `test_curvature_scale` (calculator tester): data filtered by shell, at the
    evaluated pose, scale both kernels' blocks by the weighted slope. The
    fixture asserts that it separates the weighted slope from the unweighted
    one.
  - E12 (`test_tiny_accepted_reduction_stop`, pose optimizer tester): a first
    proposal predicted below the threshold ends the stage unevaluated, with
    the pose unchanged and the stage bound-rejected, since its step was
    capped.
  - The cumulative-guard test (E12, joint) accepts a stage that ends before
    the eighth rejection: growing damping shrinks the rejected proposals until
    their predicted reduction stops the stage.
- **Check:**
  - a robustness profile on truth, before and after: `test=profile` on
    `cont_refine3D_1jxy`, and the benchmark simulation over its seeds;
  - the benchmark's iterations per particle and its time;
  - beta-gal FSC within the run-to-run spread.
- **Done when:** the tolerance values are set from the profile, and the truth
  metrics are no worse than the before-profile's spread.

#### Stage 4 record (2026-10-20)

Every build is a snapshot of the working tree on `e560a49fe`, in a worktree of
`~/src/SIMPLE_cont_bench` (patches in `patches/`). Runs use 16 threads, alone
on the machine.

| Build | What it adds |
| --- | --- |
| S3 | Stage 3 only: the reference of this record |
| S4a | The curvature scale, with an unweighted slope |
| S4b | S4a and the predicted-reduction stop |
| S4c | S4b with the `\|f\|^2`-weighted slope: **what shipped** |
| S4e | S4c with a Gauss-Newton decrement stop in place of both stops (dropped, finding 4) |

**Beta-gal** (orientation search in s, the mean of two order-balanced runs;
the work counts of the two runs are identical):

| Build | Search | Passes per solve | Trials per solve | Rejected | Accepted at 0.25-0.75 | Accepted at 0.75-1.25 |
| --- | --- | --- | --- | --- | --- | --- |
| S3 | 16.46 | 11.75 | 9.75 | 0.48 | 0.64 | 0.17 |
| S4a | 16.39 | 11.46 | 9.46 | 0.48 | 0.61 | 0.21 |
| S4b | 12.51 | 8.69 | 6.69 | 0.41 | 0.72 | 0.22 |
| S4c | 9.68 | 6.68 | 4.68 | 0.31 | 0.45 | 0.42 |
| S4e | 12.61 | 8.66 | 6.66 | 0.39 | 0.40 | 0.38 |

Every build gives the same map: FSC 3.752 A at 0.143 and 4.352 A at 0.5, cFAR
0.558. Against stage 2's build B (19.3 s, 13.75 passes per solve), stages 3
and 4 together halve the orientation search.

**Simulated polish** (seed offset 0, one run each; errors are medians over
10 000 particles, rotation in degrees and shift in A):

| Build | cc search | cc passes | cc errors | euclid search | euclid passes | euclid errors |
| --- | --- | --- | --- | --- | --- | --- |
| S3 | 20.6 | 12.43 | 0.0864 / 0.0297 | 20.6 | 12.88 | 0.0845 / 0.0288 |
| S4a | 18.3 | 10.81 | 0.0863 / 0.0296 | 19.7 | 11.95 | 0.0847 / 0.0289 |
| S4b | 11.1 | 6.43 | 0.0862 / 0.0296 | 12.5 | 7.57 | 0.0848 / 0.0289 |
| S4c | 10.3 | 6.04 | 0.0864 / 0.0296 | 11.7 | 7.11 | 0.0849 / 0.0288 |
| S4e | 15.1 | 8.70 | 0.0864 / 0.0296 | 16.4 | 9.78 | 0.0848 / 0.0288 |

Every error lies inside the before-profile's spread over six seed offsets
(euclid 0.0834-0.0854 degrees, cc 0.0850-0.0873 degrees).

**Capture with the rotation bound at 20 degrees** (seed offset 0, S3 against
S4c):

| Case | Search | Passes per solve | Rotation error | Shift error |
| --- | --- | --- | --- | --- |
| euclid | 3.56 to 2.67 (-25%) | 14.88 to 10.96 | 0.2732 / 0.2729 | 0.0856 / 0.0855 |
| cc | 3.16 to 2.74 (-13%) | 12.85 to 10.91 | 0.3471 / 0.3482 | 0.0995 / 0.0993 |

Per particle the two builds agree: the median distance between their
endpoints is 0.005 degrees and 0.001 A for euclid, and each build is closer to
the truth for half of the particles. No solve reaches the iteration cap.

What the record says:

1. **An unweighted slope does not correct euclid on real data.** On beta-gal
   it stays near 1 (S4a hardly changes the histogram). The overall amplitude
   of the reference matches the data's, but the matrix takes its curvature
   from the high frequencies, where a filtered reference is weaker than the
   data, and there the slope is above 1. Weighting the slope by `|f|^2`, as
   the matrix is, halves the over-steps: accepted steps below 0.75 fall from
   64% to 45%, and those at 0.75-1.25 rise from 17% to 42%. A scalar slope
   cannot follow every axis and shell, so 45% remain; a finer correction is
   left open.
2. **The predicted-reduction stop is the larger lever.** From S4a to S4b it
   removes 2.8 passes per solve on beta-gal and 4.4 on polish, with no change
   in the map or the truth metrics. It also removes the high-ratio tail
   (accepted steps at or above 2.5 fall from 9.5% to 0.4%): those were the
   tiny final steps.
3. **The capture regime measured the bound, not the solver.** Its seeds lay
   exactly 15 degrees from the truth, and `athres_cont` was 15 degrees, so the
   truth lay on the cumulative bound, and about half of the particles had
   their noisy optimum just outside it. With that bound, every scaled build
   gave worse capture endpoints than S3 (euclid shift error 0.109-0.113 A
   against 0.095 A, its 90th percentile 0.26-0.28 A against 0.185 A), and the
   particles that moved ended at a worse objective. S4e, which has none of the
   damped stops, gave the same endpoints, so the cause was not an early stop.
   With the bound at 20 degrees the builds agree particle by particle (table
   above), so the difference was where each path met the bound. The gate test
   E25 (`cont_refine3D_1jxy`) has the same 15-against-15 design; it is a stage
   5 item.
4. **Why the decrement stop was dropped.** The Gauss-Newton decrement,
   `g^T H^-1 g / 2`, is the reduction the undamped step predicts, and damping
   cannot shrink it. S4e stopped on it alone, to rule out stops forced by
   damping. A trace of 36 beta-gal solves (a throwaway build, T) shows why
   that is worse:
   - most solves converge in 2 to 5 accepted steps with actual/predicted
     near 1, and end on the decrement;
   - near the optimum the objective is not quadratic at about `1e-6` of
     itself. A step of 0.015 degrees, predicted to gain `1.2e-6`, lost
     `2.5e-6`. Steps of `1e-5` to `1e-8` degrees changed the objective by
     `1e-7` to `1e-10` of itself, with random sign;
   - with only the decrement, such a solve keeps trying: one took 26
     iterations for 11 accepted steps of no value. Solves pressed against the
     cumulative bound creep to the iteration cap (3 192 to 3 464 of 10 000 in
     the old capture regime).

   The damped step never predicts more than the undamped one, so whenever
   the decrement is below the tolerance, so is the damped prediction. The
   predicted-reduction stop therefore also ends every converged stage, and the
   decrement test adds nothing. Finding 3 shows that its stops forced by
   damping cost nothing measurable.
5. **The tolerance stays at `1e-6`.** In the traces, the Gauss-Newton step
   that remains at a decrement of `1e-5` is 0.05-0.08 degrees, the size of the
   polish precision (0.085 degrees); at `1e-6` it is 0.02-0.03 degrees. A
   step tolerance in physical units (section 3.3) is not needed.

**After-profile of S4c** (taken with the shipped build; the before-profile is
section 2.6's, with build B, whose results stage 3 left unchanged):

- **Beta-gal, S3 against S4c** in one order-balanced batch (S3, S4c, S4c,
  S3): 16.18 and 16.20 s against 9.34 and 9.38 s, -42%. The map is unchanged.
- **Polish over six seed offsets** (0 to 5), against the same offsets of the
  before-profile:

  | Metric | Before: median (min-max) | After: median (min-max) | Paired difference: median (min..max) |
  | --- | --- | --- | --- |
  | euclid median rotation error (deg) | 0.0844 (0.0834-0.0854) | 0.0844 (0.0834-0.0853) | -0.0001 (-0.0002..+0.0004) |
  | euclid median shift error (A) | 0.0292 (0.0288-0.0293) | 0.0291 (0.0288-0.0294) | +0.0000 (-0.0001..+0.0001) |
  | cc median rotation error (deg) | 0.0861 (0.0850-0.0873) | 0.0861 (0.0851-0.0871) | +0.0001 (-0.0001..+0.0001) |
  | cc median shift error (A) | 0.0299 (0.0297-0.0301) | 0.0298 (0.0296-0.0300) | -0.0001 (-0.0001..+0.0001) |
  | euclid mean `corr_cart` | 0.6722 | 0.6722 | 0 to five digits |
  | cc mean `corr_cart` | 0.7789 | 0.7789 | 0 to five digits |

- **The gate over 20 seeds** (`test=profile target=cont_refine3D_1jxy
  nseeds=20`, 8 threads): 20 of 20 pass, and every metric of the
  before-profile's table stays inside its spread:

  | Metric | Before: median (spread) | After: median (spread) |
  | --- | --- | --- |
  | a (euclid) median rotation error (deg) | 0.943 (0.035) | 0.937 (0.034) |
  | a (euclid) median shift error (px) | 0.071 (0.002) | 0.071 (0.001) |
  | b (cc) median rotation error (deg) | 0.229 (0.003) | 0.229 (0.003) |
  | b (cc) median shift error (px) | 0.059 (0.001) | 0.058 (0.001) |
  | c (distributed) difference to a (deg) | 0 | 0 |
  | e (auto, `pose_cont=only`) median rotation error (deg) | 0.980 (0.306) | 0.976 (0.326) |
  | f (polish) symmetry-aware error (deg) | 0.423 (0.019) | 0.424 (0.018) |
  | g (auto with polish) symmetry-aware error (deg) | 0.176 (0.022) | 0.173 (0.022) |
  | h (continuation) symmetry-aware error (deg) | 0.174 (0.030) | 0.181 (0.033) |

- **Checks of the working tree:** the incremental build, the fast gate (15 of
  15), and `unit_cart_align3D` (772 checks) pass.

The S3 and S4c worktrees carry the 20-degree capture bound on top of their
snapshots; the S4e and T worktrees are kept for reference only.

### Stage 5: kernel oracle and test gaps

**Status (2026-10-20): done on the working tree, tests only (record below).**

- **What:**
  - the independent checks of section 7 that stage 6 needs: an independent
    `J^H J`, the current `pose_normal_terms` kept as the reference, half-plane
    against full-disk sums, and the 1-against-`num_procs` thread rule;
  - the test gaps of section 7: large tester objects as `class(T), allocatable`,
    a negative control in `cont_refine3D_1jxy`, and a gate stage on the
    production `joint` route.
- **Results:** tests only.
- **Done when:** the new tests pass on the current kernel. They define the
  target stage 6 must meet.

#### Stage 5 record (2026-10-20)

**What was built:**

- **The Gauss-Newton oracle** (`test_gauss_newton_oracle`, calculator tester).
  At a pose where the data equal the model, the residual term of the Hessian
  vanishes and the curvature scale of stage 4 is exactly one. Both kernels'
  matrices are then the normal matrix of the model's Jacobian:
  - the 5x5 matrix of `pose_normal_terms` is compared with one built from
    per-sample finite differences of predicted planes (`predict`), finalized
    as the kernel finalizes cc and euclid;
  - the 2x2 matrix of `shift_normal_terms` is compared with second differences
    of objective values alone;
  - for cc and euclid, without CTF, with CTF and phase-flipped.

  Agreement: 8.3e-5 (5x5) and 3.4e-6 (2x2) of the diagonal scale, against a
  tolerance of 2e-3 derived from the single-precision samples and the step
  sizes.
- **The thread rule.** N26 (the calculator: objective, gradient, residual,
  both normal-term blocks and the pose transaction of 24 particles) and N27
  (the strategy's batch through `srch`) run once at one thread and then 20
  times at `omp_get_num_procs()` threads (24 on the M2 Ultra), compared bit for
  bit. Both pass. They replace a fixed 3-thread team run once.
- **Large tester objects.** The calculator, strategy and 1JYX testers hold
  their `parameters`, `builder`, `cmdline` and `sp_project` objects as
  `class(T), allocatable` locals (15 procedures) and take `class(T)` dummies
  (6 helpers), as `compile_time_policy.md` asks.
- **The gate** (`cont_refine3D_1jxy`):
  - **a negative control for every E25 stage** (a, b, c, i, j): the same
    rotation error and basin fraction against the truth shuffled among the
    particles (each particle against its successor's pose) must fail the
    floors. They come out at 130 degrees and 0. The simulation has no shift,
    which shuffling keeps, so the shift error has no such control;
  - **stages i and j:** stages a and b on the production route, the joint
    stage alone, with the same floors. They run last, so the earlier stages'
    random draws, and with them their results, stay as they were.

**Finding: the gather jumps where a sample changes stencil.** The first form
of the oracle took second differences of the objective in rotation as well.
At a step of 1e-3 radians they came out up to 26% above the kernel's matrix
(2% with CTF); at 5e-4 they agree with it to 1e-5. Between the two steps a few
samples cross a half-integer offset, where the 3x3x3 stencil moves to the next
grid cell and the interpolated value jumps (`test_stencil_switch` measures the
jump; the KB window is not zero at its edge). Each crossing adds its jump
squared over the step squared. So:

- the kernel's matrix is right, and the oracle now uses per-sample
  differences: where the forward and backward differences of a sample
  disagree, the side that crossed carries the jump over the step, and the
  other is the derivative of the sample's piece at the pose. At the
  oracle's step 14 samples of each plane take a one-sided difference, and the
  test asserts that some do;
- the objective is piecewise smooth in rotation, with jumps large enough to
  change a curvature estimate over 0.06 degrees by a quarter. This is a likely
  cause of what stage 4 found near the optimum: the objective is not
  quadratic at about `1e-6` of itself, and steps of 0.015 degrees lost what
  they were predicted to gain (stage 4 record, finding 4);
- a gather without the jump would make the objective smooth and the last LM
  steps reliable: a window that vanishes at its edge, or a stencil that does
  not switch. That changes results, so it went to the maintainer, who ruled
  it in (the continuous gather, below).

**Moved to stage 6:** the kept reference (today's `pose_normal_terms`, kept to
compare the rewrite against) and the half-plane against full-disk check. Both
compare against code that exists only once the rewrite begins; a copy of the
code under test, made now, would check nothing.

**Not changed:** E25's 15-against-15 bound (the gate holds floors, and the
benchmark makes the comparisons, with its 20-degree capture regime), the
euclid reference scale and the small-angle metric (section 7).

**Checks:** the incremental build; the fast gate (15 of 15);
`unit_cart_align3D` (780 checks); `lib_cart_align3D` (55 checks, 27 s); the gate
once (371 s at 8 threads, every row passing); and the gate's profile over 20
seed offsets at 8 threads:

- 20 of 20 pass;
- 48 of the 60 metrics it shares with the stage 4 profile are identical (stages
  a-d and f). The other 12 come from stages e, g and h, the `refine3D_auto`
  and polar runs, which are not reproducible at 8 threads; each stays inside
  its spread (stage e's median rotation error 0.976 to 0.981 degrees, spread
  0.33);
- the floors of stages i and j are admissible by the floor rule, far from their
  bounds: stage i rotation median 0.932 degrees, spread 0.029, worst 1.009;
  shift worst 0.072 px; stage j rotation worst 0.237 degrees; basin 1.000 in
  every run of both;
- the shuffled controls come out at 128.8-135.2 degrees, with a basin fraction
  of at most 0.001.

### The continuous gather (ruled 2026-10-20)

**Status (2026-10-20): done on the working tree (record below).** Stage 5
found that the gather jumps where a sample's nearest grid point changes; the
maintainer ruled to make it continuous.

**The cause.** The gather takes three taps per axis around the nearest grid
point, with KB weights normalized per axis. At a half-integer offset the
stencil moves to the next cell: one tap at distance 1.5 leaves and another
enters. The KB window is 1/3 at its edge, about 1% of its centre, so the
leaving tap carries about 0.7% of the weight, and the value jumps by that much
of the difference between the two outer coefficients. The derivative jumps
too.

**The choice.** A window whose value and slope vanish at 1.5 makes the value
and the derivative carry across every switch, on the same three taps. Two such
windows were compared in a one-dimensional model of the gather: the padf=2
transform of an object in the central half of the padded box, sampled off grid
and compared with the exact transform of the tapered object, the taper being
the stencil's own on-grid response, as `get_ptcl_taper` makes it.

| Window | Value error | Derivative error | Jump at a switch |
| --- | --- | --- | --- |
| KB (before) | 4.1e-3 | 6.4e-2 | 1.1e-2 |
| KB less its edge value and slope | 1.1e-2 | 1.4e-1 | 0 |
| Quadratic B-spline | 3.2e-3 | 3.2e-2 | 0 |

The quadratic B-spline fits the three taps exactly (support [-1.5, 1.5]), is
continuous with its first derivative, and its weights sum to one at every
offset, so the per-axis normalization goes. It is also more accurate than the
KB stencil in this model, and its derivative twice as accurate. Its weights
are three quadratics per axis instead of nine degree-14 polynomials.

**Keeping the polar link.** The particle carries the stencil's taper (O5 of the
pose_cont refactoring), and the polar and Cartesian representations share one
sigma2 store, which the taper's noise-power factor scales. A B-spline taper on
the particle would shift that factor by about 8% at a typical mask, which the
sigma agreement tests (N14, N33) bound at 4-9%. So the particle keeps the KB
taper, and the reference is weighted instead: `set_ref` multiplies the padded
volume, per axis, by the KB stencil's response over the B-spline stencil's
(`ref_weight`, 0.92 at the native box edge). On grid, the B-spline gather of
the weighted transform equals the KB gather of the unweighted one exactly
(circular convolution), so the model is the polar branch's there. The model
puts this at the same accuracy against the KB-tapered target as the plain
B-spline (3.1e-3, derivative 3.2e-2).

**Changes:** `gather_packed_bspline3_window_grad` and
`bspline3_stencil_envelope_1d` (`simple_cartesian_fourier.f90`) replace
`gather_packed_kb3_window_grad`; `cartft_calc` gathers with it and weights its
references by `ref_weight`. The KB machinery of PCG and the polar branch is
unchanged. The KB interpolation policy has a section on it (section 6a).

**Tests:**

- `test_fused_bspline3_gather_oracle` (Cartesian Fourier tester): the fused
  gather against the generic gather given the B-spline stencil built from the
  spline's piecewise form (3.2e-8 and 3.9e-8 of the largest coefficient), and
  across half-grid switches, where the value moves by its derivative times the
  step (ratio 1.06) and the derivative by at most 1% of its smoothness bound;
- E10: the gather equals the projector's KB stencil on grid (1.7e-7 and 2.7e-7
  relative), and off grid within the two kernels' interpolation errors (1.3e-2
  against a bound of 2e-2);
- the Gauss-Newton oracle takes second differences of the objective in
  rotation as well; they agree with the kernel's matrix to 4.1e-5, where the
  KB gather gave 26%;
- the sigma agreement tests N14, N31 and N33 pass unchanged (N33's worst shell
  at 0.98 of its bound of 3).

#### Record (2026-10-20)

G is the stage 5 working tree with the continuous gather, against S4c, whose
cont code stage 5 left unchanged. Timed alone, 16 threads.

| Case | Search (s) | Passes per solve | Rejected trials | ns per pixel-evaluation | Truth |
| --- | --- | --- | --- | --- | --- |
| beta-gal | 9.42 to 5.36 (-43%) | 6.68 to 5.70 | 0.31 to 0.05 | 215 to 143 | FSC 3.752 A, cFAR 0.558, unchanged |
| polish cc | 10.29 to 5.45 (-47%) | 6.04 to 4.68 | 0.27 to 0 | 219 to 150 | 0.0864 to 0.0860 deg, 0.0296 to 0.0295 A |
| polish euclid | 11.60 to 6.13 (-47%) | 7.11 to 5.72 | 0.27 to 0 | 210 to 138 | 0.0849 to 0.0849 deg, 0.0288 to 0.0288 A |
| capture cc | 2.74 to 1.24 (-55%) | 10.91 to 7.30 | 0.34 to 0 | 228 to 155 | 0.348 to 0.344 deg |
| capture euclid | 2.68 to 1.29 (-52%) | 10.96 to 7.96 | 0.31 to 0 | 221 to 148 | 0.273 to 0.271 deg |

Beta-gal times are the means of two order-balanced runs, whose work counts are
identical; the simulated cases are one run each at seed offset 0.

What it says:

- **The jumps caused every rejected trial on simulated data.** With a
  continuous gather no trial is rejected, and 73-89% of accepted steps land in
  the 0.75-1.25 bin. On beta-gal 5% remain, with 52% of accepted steps short of
  their prediction: real data still over-steps a little.
- **Each evaluation is a third cheaper,** from the weights alone.
- **Against stage 2's baseline,** the beta-gal orientation search went from
  19.3 s to 5.4 s (-72%), and the polish cases from about 24 s to 5.5-6.1 s.

**After-profile** (G, against the stage 5 profile and the polish runs of S4c):

- **The gate over 20 seed offsets** at 8 threads: 20 of 20 pass, every
  metric of stages a-j stays inside the stage 5 profile's spread, and every
  floor stays admissible (stage a rotation worst 1.016 degrees, stage i 0.991,
  against 1.5). Stage a's median moved from 0.937 to 0.950 degrees, a third of
  its spread; the controls are unchanged.
- **Polish over six seed offsets** (0 to 5), paired against the same offsets
  with the KB gather:

  | Metric | KB: median (min-max) | B-spline: median (min-max) | Paired difference: median (min..max) |
  | --- | --- | --- | --- |
  | euclid median rotation error (deg) | 0.0844 (0.0834-0.0853) | 0.0843 (0.0831-0.0849) | -0.0002 (-0.0005..-0.0000) |
  | euclid median shift error (A) | 0.0291 (0.0288-0.0294) | 0.0291 (0.0288-0.0294) | -0.0000 (-0.0001..+0.0000) |
  | cc median rotation error (deg) | 0.0861 (0.0851-0.0871) | 0.0859 (0.0846-0.0867) | -0.0003 (-0.0005..-0.0001) |
  | cc median shift error (A) | 0.0298 (0.0296-0.0300) | 0.0298 (0.0295-0.0300) | -0.0000 (-0.0001..+0.0001) |
  | euclid mean `corr_cart` | 0.6722 | 0.6732 | +0.0010 in every seed |
  | cc mean `corr_cart` | 0.7789 | 0.7791 | +0.0003 in every seed |

  The rotation error falls in every seed, by 0.2-0.3%, and the model fits the
  data a little better (higher `corr_cart` in every seed).
- **Checks of the working tree:** the build, the fast gate (15 of 15),
  `unit_cart_align3D` (788 checks), `lib_cart_align3D` (55 checks).

### Stage 6: cheaper evaluations

**Status (2026-10-21): done on the working tree (record below).** 6a, 6b and
6c are in. The centred cube was measured again, on the half plane, and
dropped.

- **6a. One rewrite of `pose_normal_terms` and its siblings.** It covers the
  Friedel half plane, loop bounds set by `kto`, the separable shift phase and
  the normal-term arithmetic (sections 4.1, 4.2, 4.3 and 4.5).
  `shift_normal_terms`, `sigma_contribution` and the particle power of
  `fill_ptcl` follow the half plane. The particle power is summed in double
  precision in the kernel's order, so the CC normalization matches the
  kernel's sums. The slots store `transfer` as real.
- **The centred cube** (sections 4.4 and 4.5) was measured on top of the half
  plane. It is slower on beta-gal, so it was dropped again (record).
- **6b. SIMD across taps.** The gather contracts its 27 taps one axis at a
  time (x, then y, then z). It forms the packed/Friedel indices once per axis
  instead of once per tap. The sums are reordered.
- **6c. sigma2 from the accepted evaluation** (section 3.2, row 3).
  - `pose_normal_terms` accumulates the whitened `|r|^2` per shell when asked.
  - `cartft_pose_opt` keeps it for the evaluation that chose the committed
    pose: the seed, the shift endpoint on the seam route, or the last
    accepted joint step.
  - `calc_sigma2` unwhitens it per shell (`sigma_from_shell_residual`:
    `sigma2 res / (2 npix)`, with the shell's full-disk pixel count).
  - The separate `sigma_contribution` pass goes.
- **Results:** rounding changes.
- **Kept reference.** The full-disk kernels stay as
  `pose_normal_terms_full_disk` and `shift_normal_terms_full_disk`, and the
  tests compare the half-plane kernels with them.

#### Stage 6 record (2026-10-21)

Each build is a snapshot of the working tree (`patches/`):

- **G:** the stage 5 tree with the continuous gather.
- **H:** G with 6a.
- **I:** H with the centred cube.
- **J:** I with 6c.
- **K:** H with 6c, the kept tree before 6b.
- **L:** K with 6b, the working tree.

All runs were timed alone at 16 threads. Beta-gal times are the means of
order-balanced pairs. The simulated cases are one run each at seed offset 0.

**6a (G to H).**

| Case | Search (s) | ns per pixel-evaluation | Passes per solve | Truth |
| --- | --- | --- | --- | --- |
| beta-gal | 5.50 to 2.40 (-56%) | 146.7 to 64.0 | 5.70, unchanged | FSC 3.752 A, cFAR 0.558, unchanged |
| polish cc | 5.39 to 2.23 (-59%) | 147.9 to 61.2 | 4.68, unchanged | 0.0860 deg, 0.0295 A, unchanged |
| polish euclid | 6.12 to 2.68 (-56%) | 138.2 to 60.2 | 5.72, unchanged | 0.0849 deg, 0.0288 A, unchanged |
| capture cc | 1.24 to 0.40 (-68%) | 154.8 to 49.6 | 7.30, unchanged | 0.344 deg, unchanged |
| capture euclid | 1.30 to 0.44 (-66%) | 148.5 to 49.3 | 7.96, unchanged | 0.271 deg, unchanged |

Work counts and truth metrics agree to the digits the bench prints. The half
plane halves the pixels. The rest takes a further 13-36% off each pixel:

- the phase tables;
- the 15 distinct Hessian entries;
- the cross-product Jacobians;
- the `kto` bounds.

**The centred cube (H to I).**

| Case | Search (s) | ns per pixel-evaluation | Reference memory |
| --- | --- | --- | --- |
| beta-gal | 2.35 to 2.77 (+18%) | 62.7 to 74.0 | 1.08 to 0.73 GB (peak RSS 5.48 to 4.47 GiB) |
| polish cc | 2.27 to 2.27 | 62.2 to 62.2 | - |
| polish euclid | 2.61 to 2.59 | 58.9 to 58.3 | - |
| capture cc | 0.41 to 0.39 | 50.3 to 48.2 | - |
| capture euclid | 0.45 to 0.42 | 50.7 to 46.7 | - |

At beta-gal's band (shell 88 of 128) the cube saves a third of the reference
memory but costs 18% of the search time. On the small simulated box it breaks
even. As in stage 3, the packed layout wins at a production band, so it stays.
I's snapshot (`patches/tree_snapshot_20261021_1042_I.patch`) holds the cube.

**6c (H to K).** One pass fewer per euclid solve; cc is unchanged.

| Case | Search (s) | Passes per solve | Truth |
| --- | --- | --- | --- |
| beta-gal | 2.38 to 2.01 (-16%) | 5.70 to 4.70 | FSC 3.752 A, cFAR 0.558, unchanged |
| polish euclid | 2.61 to 2.28 (-13%) | 5.72 to 4.72 | unchanged |
| capture euclid | 0.45 to 0.41 | 7.96 to 6.96 | unchanged |
| polish cc | 2.27 to 2.28 | 4.68 | unchanged |
| capture cc | 0.41 to 0.41 | 7.30 | unchanged |

- **The cube tree gave the same** (I to J): beta-gal 2.77 to 2.37 s.
- **The sigma2 store** that beta-gal writes agrees with the separate pass's to
  rounding. Of its 4199 x 128 entries, 44% differ, by at most 8.6e-6 relative
  (99th percentile 9.7e-7).
- **The ns per pixel-evaluation rises by 2%** (63.4 to 64.9). The per-shell
  sum costs a little, and the per-solve costs are spread over fewer passes.

**6b (K to L).**

| Case | Search (s) | ns per pixel-evaluation | Truth |
| --- | --- | --- | --- |
| beta-gal | 2.10 to 1.90 (-9%) | 67.8 to 61.5 | FSC 3.752 A, cFAR 0.558, unchanged |
| polish cc | 2.33 to 2.14 (-8%) | 64.0 to 58.8 | unchanged |
| polish euclid | 2.34 to 2.08 (-11%) | 63.9 to 56.8 | unchanged |
| capture cc | 0.42 to 0.38 (-10%) | 52.8 to 48.0 | unchanged |
| capture euclid | 0.41 to 0.37 (-10%) | 52.4 to 48.2 | unchanged |

- **The pairs.** Beta-gal is the mean of two batches of order-balanced pairs,
  eight runs in all. The simulated cases come from one batch. In these
  batches the machine ran K 5% slower than in the morning's (2.10 against
  2.00 s), with a backup and another session's tests running. The pairs share
  that load.
- **A standalone model of the gather** (padded box 512, band 88, random
  planes, 16 threads) chose the form:
  - the contraction took the gather from 145 to 121 thread-ns per sample with
    the reference out of cache, and from 76 to 62 in cache;
  - collecting the 27 taps first and contracting after was slower (127 and
    69).

  In the pass, consecutive LM evaluations of a particle sample nearly the same
  plane, so the reference is mostly in cache, and the pass gains 9-11%.

**Against stage 2's baseline,** the beta-gal orientation search went from
19.3 s to 1.9 s (-90%), and the polish cases from about 24 s to about 2.1 s.
Against G, the tree stage 6 started from, beta-gal went from 5.50 s to 1.90 s
(-65%). These comparisons span batches, so they hold to about 5%.

**Checks of the working tree:**

- the build and the fast gate (15 of 15);
- `unit_cart_align3D` (801 checks) and `lib_cart_align3D` (55 checks);
- `check_descr.py` and `check_test_registry.py`;
- `test_half_plane_kernels`: the half-plane kernels against the full-disk
  ones, on Hermitian noise. It covers cc and euclid, with and without CTF, and
  bands [2,6] and [0,6]. The objective, the gradient and the Gauss-Newton
  matrix agree within 6.3e-8, 1.1e-7 and 1.3e-7 relative (bound 1e-5).
- `test_sigma_from_shell_residual`: the conversion against
  `sigma_contribution` at the same pose, bands from shell 2 and from DC:
  identical.
- `test_committed_shell_residual`: for an accepted solve, a rejected one and
  the seam route, the kept residual is bit-identical to a fresh evaluation at
  the returned pose.
- The stage 5 oracles:
  - the gather oracle within 9.5e-8 and 1.3e-7 of the largest coefficient
    (bound 3.8e-6), with switch ratios 1.005 and 0.009;
  - the Gauss-Newton oracle within 4.5e-5 (rotation) and 3.1e-6 (shift);
  - E10 within 1.6e-7 and 2.2e-7 on grid, and 1.3e-2 off grid.
- **The gate over 20 seed offsets on L,** at 8 threads:
  - 20 of 20 pass;
  - every one of its 94 metrics has its median inside the min-max of G's
    after-profile;
  - every floor stays admissible: stage a rotation worst 1.016 degrees and
    stage i 0.991, against 1.5; stage e 1.98 against 7.5; stage h continued
    0.295 against 1.53;
  - stages a-d and i-j agree with G's to four digits, and the
    multi-iteration stages e-h move within their spread.

### Stage 7: close-out

**Status (2026-10-21): done.**

- The Cartesian section of `refine3D_policy.md` (continuous Cartesian pose
  policy) covers the reference layout, the continuous gather, the half-plane
  evaluation, the curvature scale, the stopping rule, the passes per solve,
  sigma2 from the accepted evaluation and the counters.
- The `production/CMakeLists.txt` comment of `lib_cart_align3D` says 1 000
  particles and that the fixture files are removed.
- The measured gains of every stage are in section 2, and the note moved to
  `completed/`.

## 3. Fewer passes per particle (the LM driver)

### 3.1 The CC Gauss-Newton matrix lacks the correlation factor

**What the code computes.** `finalize_cc_normal_terms`
(`simple_cartft_calc.f90:958-981`) turns the sums into three things:

- the objective `f = 1 - <x^, n(theta)>`, where `x^ = x/|x|` and `n = m/|m|` is
  the unit model;
- its gradient `-<x^, J_n>`;
- the matrix `H_GN = J_n^T J_n`, computed as
  `hessian/prediction_power - mg mg^T/prediction_power**2`.

**Why it is too large.**

1. The exact Hessian is `-<x^, d2n>`.
2. Because `|n| = 1`, we have `<n, d2n> = -|J_n|^2`.
3. Write `x^ = c n + sqrt(1-c^2) u`, where `c` is the correlation and `u` is the
   noise direction. The noise term averages to zero.
4. So `H ~ c J_n^T J_n`, and `H_GN` is too large by a factor `1/c`.

**Consequence.** Each LM step is `c` times the Newton step. The pose error
shrinks only by `(1 - c)` per iteration, and the objective gap by `(1 - c)^2`.

**The beta-gal record fits this.**

- At `c ~ 0.17`, with the 1e-6 relative stop of section 3.3, the contraction
  predicts about 18 accepted steps.
- The run averaged about 25 LM iterations per particle, over a spread of
  correlations.
- The pre-refactor refiner had the same finalize and the same LM constants:
  lines 684-706 and 43-47 of
  `git show 214fad5d6^:src/main/volume/simple_cartesian_pose_refiner.f90`.

**But that `c` is a whitened correlation.** The old objective was a normalized
correlation on sigma-whitened data, where the noisy high shells weigh as much
as the low ones, so its correlation was low. Today's `objfun=cc` weights every
pixel equally, and on the same data its correlation is far higher. In the
smoke runs (stage 2) it is 0.65-0.78 at SNR 0.02-0.05, with about 10 LM
iterations per solve, and the steps are only about 1.3 times too short. The fix
is still right, but at realistic SNR its gain is modest. The large factor
returns for any objective whose correlation is low, such as one on whitened
data.

**Signature to look for with the counters of section 5.** An undamped step
predicts a reduction of `(1/2) g^T H_GN^-1 g`, and the true quadratic gives
`(1 - c/2) g^T H_GN^-1 g`. Accepted steps should therefore show
`actual/predicted ~ 2 - c`: almost no rejections, and most accepted steps in
the 1.75-2.5 bin of the histogram for a correlation below 0.25.

**Change (stage 4).** In `finalize_cc_normal_terms`, multiply the matrix by
`c = cross_real/normalizer`, with a small positive floor. This is a one-line
change. The optimum is unchanged, but the path changes, and with it the
endpoints at tolerance level.

**As built (stage 4 record).** The factor is not `c` itself but the
regression slope of the normalized data on the normalized model, weighted by
`|f|^2`; it equals `c` when every shell has the same slope. The weighting
matters because the matrix takes its curvature from the high frequencies, and
the correlation falls with frequency.

**`objfun=euclid` is different.** It is the production objective of
`refine3D_auto` and of the proposed `refine3D_local`.

- `L = sum |m - x|^2 / P`, with `H_GN = 2 J^T J / P`. This is right on average
  when the reference amplitude matches the data.
- If the reference is `beta` times the data signal, then `H_GN/H ~ beta`.
- A filtered reference (`beta < 1` at high k) therefore over-steps and pays in
  rejections and damping.

Stage 2 measures the ratio and the rejection rate. A euclid change follows
only if they show over-stepping.

**Stage 4 result.** They do on beta-gal, and an unweighted slope does not
help: it stays near 1, because the reference's overall amplitude matches the
data's. The slope weighted by `|f|^2`, the same correction as cc's, halves the
over-steps (stage 4 record, finding 1).

### 3.2 Redundant full-disk passes

| Pass | Where | Duplicate of | Change | Results |
| --- | --- | --- | --- | --- |
| Seed `objective_gradient` | `simple_cartft_pose_opt.f90:208` | The seed `pose_normal_terms` of `refine_joint` (`:457`) on the joint route | Pass the seed evaluation into `refine_joint`, or return it | Identical (stage 3a) |
| Final `objective_gradient` | `:246` | The last accepted evaluation inside `refine_joint`, or the seed | Return the accepted objective from `refine_joint` | Identical (stage 3a) |
| `sigma_contribution` (euclid) | `simple_strategy3D_cont.f90:151-155` -> `simple_cartft_calc.f90:740` | The committed pose, already evaluated | Accumulate per-shell whitened `\|r\|^2` in `pose_normal_terms`, keep it for the accepted pose, and unwhiten by sigma2 per shell. Pixel counts per shell are fixed by the box. | Rounding (stage 6c) |

**Totals.** With N trial evaluations:

- A euclid particle makes N+4 full-disk passes; a cc particle makes N+3.
- The test route `shift_then_joint` adds one more pass (`:219`).
- `sigma_contribution` also uses the derivative gather and discards the
  derivatives.

### 3.3 The stopping rule is far below the noise

**Current rule.** `LM_ACCEPTED_RELATIVE_TOL = 1e-6`
(`simple_cartft_pose_opt.f90:38`; the exit is at `:563-564`). A step ends the
solve only when it reduces an objective near 0.83 by less than about 8e-7.

**For comparison:**

- The polar continuous in-plane solve stops at 1e-4 (`JOINT_IMPROVE_REL_TOL`,
  `simple_pftc_shsrch_grad.f90:505`).
- `check_cont_conv` counts a particle as stable at 0.5 degrees and 1 A.

**Changes (stage 4):**

- Add a step tolerance in physical units: a small fraction of `cont_conv_rot`
  and `cont_conv_shift`.
- Skip the evaluation of a trial whose predicted relative reduction is already
  below the tolerance. Such a trial would be accepted and would end the loop,
  so skipping it saves the last pass.

Both change results. The tolerance values are set from the stage 4 profile,
together with section 3.1, since a correct matrix makes fewer, larger steps.

**As built (stage 4 record).** Only the second change, with the tolerance kept
at `1e-6`. The proposal's prediction is the damped one, which also ends a
converged stage, because the undamped prediction is larger. A test on the
undamped prediction alone (the Gauss-Newton decrement) was tried and dropped
(finding 4), and the step tolerance in physical units was not needed
(finding 5).

The smoke runs of stage 2 make this the larger lever at realistic SNR: 40-45%
of all trial evaluations are rejected, in both objectives.

### 3.4 Objective-only trials (not planned)

This is listed as future work in the gather note.

- A value-only gather is cheaper, but an accepted trial then needs a second
  pass for the derivatives.
- It is worth doing only if the counters show a real rejection rate.
- After stage 4, nearly every step is expected to be accepted.

## 4. Cheaper passes (the kernel)

### 4.1 Friedel symmetry

`pose_normal_terms` loops over the full disk (`simple_cartft_calc.f90:677-680`).
Every term is Hermitian between `(h,k)` and `(-h,-k)`:

- The observation is `expand_ft` of a real image
  (`simple_matcher_2Dprep.f90:147-148`), so `X(-h,-k) = conj X(h,k)`.
- The transfer `|C|/sqrt(sigma2)` is even.
- The shift phase conjugates.
- The gather is Hermitian: `V(-loc) = conj V(loc)` and
  `dV(-loc) = -conj dV(loc)`. The rotation lever arm also changes sign, so the
  rotation Jacobians conjugate.

So every summand takes the same value at the two Friedel mates: `|r|^2`,
`Re(J^H r)`, `Re(J_a^H J_b)` and the CC sums.

**Change (stage 6a).** Loop over `h > 0`, plus the half line `h = 0, k > 0`,
with weight 2. Visit the DC pixel once if the band includes shell 0. This
halves the pixel work of every pass, including `sigma_contribution` (its
per-shell sums double the same way) and `shift_normal_terms`.

**Details:**

- Rounding changes, because the gather of the mirrored stencil sums in reverse
  order.
- The Nyquist pixels `(+-box/2, 0)` and `(0, +-box/2)` matter only when
  `kto = box/2`. Handle them explicitly.
- The particle slots can then store the half plane, which halves their memory.
  `transfer` is real but stored as complex (`:37`); storing it real halves it
  again.

### 4.2 The shift phase

**Current cost.** Two double-precision transcendental calls per pixel
(`:688-689`). The same code appears in `shift_normal_terms`,
`sigma_contribution` and `predict`.

**Change (stage 6a).** The phase separates as
`exp(2 pi i (h s1 + k s2)/box) = e_h(s1) e_k(s2)`. Build two 1-D tables per
evaluation: `2(2 kto + 1)` sin/cos pairs instead of about `pi kto^2`. Rounding
changes.

### 4.3 Normal-term arithmetic

In `:684-720`, `m = wp val` is the model and `wp = phase transfer`. All of the
following belong to stage 6a.

**Hessian entries.** All 25 entries are accumulated (`:703-708`, `:713-719`),
but only 15 are unique.

**Shift Jacobians.** They are `i f m` (`:699`). It follows that:

- the shift components of `model_gradient` are identically zero;
- `H(4:5,4:5) = sum f_a f_b |m|^2`;
- `H(rotation a, shift) = -sum f Im(conj(J_a) m)`;
- the shift gradient is `-sum f Im(conj(m) x)`.

`f_h` depends on `h`, and `f_k` is fixed along a row. These terms are therefore
row sums that are multiplied by k at the row end.

**Rotation Jacobians.** The three rotation Jacobians (`:685-697`) are the
components of `wp (dval x loc)`. One complex cross product replaces building
`dloc` and three `sum(cmplx(dval,dp)*dloc(:,axis))` calls.

**Objective branch.** The `select case(objfun)` sits inside the pixel loop
(`:700`). Either hoist it, or accumulate one set of sufficient statistics and
finalize both objectives from it. The set is: `sum |m|^2`, `sum Re conj(x) m`,
`sum Re J^H m`, `sum Re J^H x` and `J^H J`.

**Smaller items:**

- The CC branch recomputes `sum |x|^2` on every evaluation (`:710`). For a CC
  slot, `p%power` already holds it (same shells, weight 1).
- The per-pixel `matmul(real([h,k,0],dp), rotmat)` (`:681`) can be computed as
  `h r1 + k r2`, where `r1` and `r2` are the padded rotation rows formed once per
  evaluation. The gather note made keeping the `matmul` a non-goal of that
  change only.

### 4.4 The gather

The gather is `simple_cartesian_fourier.f90:141-212`.

**KB weights (stage 6b; superseded on 2026-10-20 by the continuous gather,
whose B-spline weights are three quadratics per axis).**

- Per pixel there are three elemental calls of `apod_fast_value_deriv` on
  3-vectors (`:166`). That is nine scalar degree-14 Horner evaluations, each with
  an `error stop` parameter check (`simple_kbinterpol.f90:169-196`).
- The Release flags have no LTO (`cmake/CompilerConfig.cmake:107`). Calls that
  cross module boundaries are therefore not inlined; this covers both the KB
  calls and the per-pixel gather call.

There are two options. Both must stay inside the 32-epsilon oracle contract of
the gather note:

- one 9-lane evaluation (12-lane padded) in the gather's own module;
- a table of the three weights and derivatives per fractional offset, which
  needs three lookups per pixel.

**Scaling (stage 3a).** `padsc` (= 8) multiplies the value and the three
derivatives at every pixel (`simple_cartft_calc.f90:915-916`). Fold it into the
reference when it is built. Because 8 is a power of two, the result is
bit-identical, barring overflow.

**Indexing (stage 6a; ruled 2026-10-19: a branch-free gather; moved from stage 3 on 2026-10-20).**

- Today each sample costs three wrap-table lookups, the packed/Friedel branch
  and two conditional index adds (`simple_cartesian_fourier.f90:183-211`).
- On the reference of section 4.5, each 3-tap run along x is three contiguous
  complex values at a fixed offset. There is no wrap, no Friedel branch and no
  conditional index.
- Keeping today's accumulation order keeps the result bit-identical.
- Vectorizing across the three taps reorders the sums, so it belongs to
  stage 6b.
- **Outcome (stage 6 record):** the cube lost again on the half plane and
  was dropped; 6b contracts the taps one axis at a time on the packed layout.

### 4.5 Loop bounds and the reference layout

**Loop bounds (stage 6a).** The pixel loops visit the full `(box+1)^2` square and
test `shells` (`:677-680`). Bound them by `kto`, and by `h >= 0` under section
4.1.

**Reference layout (stage 6a; ruled 2026-10-19; moved from stage 3 on
2026-10-20).** The reference becomes a full (unpacked) complex cube. Measured
alone in stage 3 it is slower at production bands. It doubles the working set
of a full-disk pass, because the packed layout lets a pixel and its Friedel
mate share memory (stage 3 record). So it is adopted only together with the
half-plane, if the measurement there says so. On the half plane it was 18%
slower on beta-gal (stage 6 record), so it was not adopted.

- **Extent.** It is centred on the origin and covers the band of the reference
  file (`kto`, the header's band). Its half-side is `padf kto + 2` voxels per
  axis. A pixel of shell at most `kto` has `|loc| < padf (kto + 0.5)`, and the
  3-tap stencil reaches one voxel beyond `nint(loc)`, so the cube holds every
  sample a pass can need.
- **Construction.** It is filled once from the padded transform, as copies for
  `h >= 0` and exact conjugates for `h < 0`. It therefore holds exactly the
  values the packed gather fetches today.
- **Memory per state and half.** `(2(padf kto + 2) + 1)^3` complex values,
  against today's 0.54 GB at `box_crop` 256 (`pose_cont_refactoring.md:325`):
  - smaller below about 0.78 of Nyquist; at beta-gal's shell 87 of 128 it is
    0.35 GB;
  - up to about twice today's size at Nyquist, 1.1 GB.

  The benchmark reports reference memory against this. The figures matter more
  under the `autoscale=no` default proposed for `refine3D_local`
  (`refine3d_local_workflow.md:107`).
- **Consequences:**
  - `predict`, a test instrument that today covers shells up to `box/2`, stops
    at the reference's band.
  - `check_eval_args` refuses a slot band beyond the reference band. With that
    guarantee, the per-pixel `inside` test can go.

**Precision.** Keep the single-precision gather and the double accumulation;
the gather note rejects lower precision. Single-precision per-pixel products
with double accumulation are possible, but they change LM decisions at the
margin. They are not planned.

## 5. Counters (stage 1)

When `l_cart_refine` is set, the REFINE3D bench record gains a
`CONTINUOUS POSE WORK` section, summed over the partition's solves. It reports:

- solves by `CARTFT_*` status;
- LM iterations, trial evaluations (evaluated trials, not `nattempted`),
  accepted steps and bound hits;
- calculator passes (the solve's evaluations plus the sigma2 residual pass) and
  pixel-evaluations (passes times the pixels of the slot's shells in the full
  disk);
- passes and LM iterations per solve;
- the solves whose joint stage reached the iteration cap (they still end
  "accepted");
- a histogram of actual/predicted at accepted steps, in the bins 0.25-0.75,
  0.75-1.25, 1.25-1.75, 1.75-2.5 and >= 2.5. A histogram rather than a mean,
  because steps whose predicted reduction is near zero give ratios that ruin a
  mean. The bin around `2 - cc` is the signature described in section 3.1;
- the search thread-seconds per pixel-evaluation, in ns.

Elsewhere in the record:

- **Thread-seconds.** The search's thread-seconds sit in a `THREAD-SECONDS`
  block written for every 3D search, as `match2D` does. They are the busy time
  per particle, so the time per pixel-evaluation covers the whole particle
  search.
- **Reference memory.** The context gives the bytes of all Cartesian reference
  transforms.
- **Reference read.** The read, padding and transform of `cart_refvols_*.bin`
  move from "startup/setup" to "reference preparation".

**Implementation.** Each thread keeps its own work record and busy time, and
they are merged after the parallel loop. Nothing is written per particle. This
mirrors the THREAD-SECONDS block and the evaluation counters of the prob-table
records.

## 6. Benchmark workflow

### 6.1 Correctness

Use the E25 fixture:

- 1JYX, box 144, 1.3 A;
- 1000 particles from `simulate_particles`;
- CTF with defocus 1.5 +- 0.5 and astigmatism 0.1;
- SNR 1;
- every pose perturbed by 15 degrees and 2 px.

Run it through two entry points:

- `simple_test_exec test=lib_cart_align3D`, at strategy level;
- `simple_test_exec test=cont_refine3D_1jxy`, through the production programs.
  It runs euclid, cc, `nparts=2`, auto, the polish and the polar continuation.

### 6.2 Identity

The CR1 identity check is how the speed runs judged `refine=cont` (case CR1 of
`identical_speed_run.md`). Run `simple_test_exec test=cont_refine3D_1jxy nthr=1`
with each build, each in an empty directory of its own. Then compare the
`metrics.tsv` each leaves there: every check value must be equal.

`test=identity` cannot judge this entry, for two reasons:

- its registry row declares no outputs;
- the gate works in a directory named after its process id, and deletes it when
  it passes.

BENCH files are never part of an identity comparison, so the stage 1 counters
keep identity.

### 6.3 Performance (stage 2)

#### The simulated benchmark, `cont_refine3D_bench`

**Form (ruled 2026-10-19).** A `simple_test_exec` program outside the unit
suite.

- It is reachable as `simple_test_exec test=cont_refine3D_bench`, but it has no
  CTest entry, like `rec3D_backends` (`test_environment_policy.md`
  section 1.3).
- It measures. It asserts nothing about time and holds no metric to a floor.

**What it does:**

1. **The truth.** It builds the embedded 1JYX model at box 256 and 0.73 A.
   This is E25's 187 A extent on a finer grid.
2. **The particles.** It simulates `nptcls` particles (default 10 000) with
   E25's CTF spread, at `snr` (default 0.05, a typical particle SNR; real
   particles are 0.02 or better, step 2.4).
3. **The project and the reference.** It writes a project with the truth CTF,
   the truth poses and the fields a discrete workflow leaves. The reference is
   `reconstruct3D` of that project, in `reference/`, with its own half-maps.
   This is the reference a converged refinement leaves, on the particles'
   amplitude scale.

   Why not the truth map: euclid's Gauss-Newton matrix is right only when the
   reference has the data's scale (section 3.1). The pdb2mrc map has an
   arbitrary scale next to noise-normalized particles. A first smoke run
   against it put 77% of euclid's accepted steps in the 1.75-2.5 bin, which
   was an artefact of the scale, not production behaviour.
4. **The cases.** For each regime it perturbs every truth pose by exactly the
   regime's error. It then runs one `refine3D refine=cont` pass per objective
   against that reference. The pass runs on the native grid (no
   `box_crop`) and the production `joint` route, with `trs` 6.5 A (E25's 5
   pixels):

   | Regime | Pose error | Band | `athres_cont` | Stands for |
   | --- | --- | --- | --- | --- |
   | polish | 2 degrees, 1 A | lp 3 A | the default, 10 degrees | the error a discrete search leaves, in a finishing band |
   | capture | 15 degrees, 2.6 A | lp 8 A | 15 degrees | E25 |

5. **The report.** For every case, `metrics.tsv` gets:
   - the seed's and the result's median rotation error (degrees) and median
     shift error (A);
   - the fraction inside the basin width at the band;
   - the mean `corr_cart` and the improved fraction;
   - the settings.

   Recording the truth metrics with the timings means a faster arm that
   converges less cannot pass unnoticed.

**Keys:**

| Key | Effect |
| --- | --- |
| `suite=polish\|capture\|list` | one regime, or print the regimes (default: both) |
| `objfun=euclid\|cc` | one objective (default: both) |
| `nptcls`, `snr` | the simulation |
| `nparts` | distributes each pass |
| `nthr` | threads |
| `seed_offset` | moves every seed (handled by the front door of `simple_test_exec`) |

**Outputs.** Everything goes into `./cont_refine3D_bench/`, and the program
refuses to start when that directory exists. It holds:

- `metrics.tsv`;
- the fixture;
- one directory per case (`polish_euclid`, `polish_cc`, `capture_euclid`,
  `capture_cc`), each with its project and its `BENCH/`.

#### Real data: the frozen beta-gal checkpoint

The workload is 5,513 particles at box 256 with D2 symmetry, the run of
`pose_cont_refactoring.md` (lines 1569-1576). Each run starts in a fresh copy of
the checkpoint directory:

```bash
simple_exec prg=refine3D mkdir=no vol1=recvol_state01.mrc mskdiam=180 pgrp=d2 objfun=euclid refine=cont maxits=1 nthr=12
```

The comparison standard, the one the fused gather was held to
(`cont_kb_gather_performance.md:101-106`):

- three order-balanced pairs from fresh copies;
- judged on the median worker orientation-search time;
- memory within 5%, and matching work counts for the identical stages.

### 6.4 Why the SNR must be realistic

Iteration counts depend on the correlation (section 3.1).

- E25 runs at SNR 1 (`simple_strategy3D_cont_1jyx_tester.f90:50`), far above
  the realistic 0.02-0.1.
- The pose-optimizer and pose-strategy unit fixtures are noise-free
  observations from the calculator's own `predict`
  (`simple_cartft_pose_opt_tester.f90:127`).

At a correlation near 1 the Gauss-Newton matrix is nearly exact. No existing
cont test can therefore show the problem of section 3.1, and timings taken from
them would understate the iteration count.

## 7. Gaps against the test and benchmark standard (stage 5)

| Area | Library standard | Cont route |
| --- | --- | --- |
| Work counters | Prob tables write thread-seconds and evaluation counts | Orientation-search time only; the stage counts never leave the optimizer (stage 1) |
| Kept harness | Developer benchmark, `BENCH/`, plot script | The gate deletes its fixture, `BENCH/` included, on pass (`simple_commanders_test_highlevel_cont.f90:26`); the fused-gather evidence stayed local to one checkout; the joint-route timing is "not measured yet" (stage 2) |
| Route under test | Production defaults | Stages a-c and the library suite force `shift_then_joint` (`simple_commanders_test_highlevel_cont.f90:495`, `simple_strategy3D_cont_1jyx_tester.f90:464`); production uses `joint`. **Stage 5:** stages i and j run a and b on the production route |
| Fast-path oracle | Replaced routine kept and compared | Only the gradient is checked, against finite differences; the 5x5 matrix and `shift_normal_terms` have no independent check. **Stage 5:** `test_gauss_newton_oracle` checks both; the kept routine comes with the rewrite (stage 6) |
| Threads | 1 against `omp_get_num_procs()` threads, 20 repeats, bit for bit (speed2D rule) | A fixed 3-thread team, run once. **Stage 5:** N26 and N27 follow the rule |
| Negative control | Required for every gate | None in `cont_refine3D_1jxy`. **Stage 5:** every E25 stage against the shuffled truth |
| Reference scale | The objective's assumptions met | The gate's euclid stages (a, c) run against the pdb2mrc truth map, whose amplitude scale does not match the noise-normalized particles. Their median rotation error is 0.94 degrees, against 0.23 degrees for cc (stage b) on the same data (before-profile, stage 2). The benchmark reconstructs its reference for this reason. |
| Small-angle metric | Resolves the tested differences | The gate's `acos` of a single-precision trace resolves about 0.002 degrees near zero; enough for its floors near 1 degree, not for comparing stages |
| Bounds clear of the truth | A bound under test does not sit on the answer | E25 seeds the poses exactly 15 degrees from the truth with `athres_cont` 15, so the truth lies on the cumulative bound and about half of the particles end against it. The benchmark's capture regime had the same design and measured the bound rather than the solver (stage 4 record, finding 3); it now uses 20 degrees |
| Large tester objects | `class(T), allocatable` | Plain `type(parameters)`/`type(builder)` locals (`simple_strategy3D_cont_tester.f90:75-76`). **Stage 5:** allocatable in the three cont testers |
| Comments | Current | `production/CMakeLists.txt:243-245` says 5,000 particles and kept TSVs; the code uses 1,000 and removes its files (stage 7) |

The oracle checks stage 6 needs:

- an independent `J^H J`, built from per-pixel centred finite differences of the
  predicted model, compared with the 5x5 matrix of `pose_normal_terms`, for cc
  and euclid under all three CTF flags;
- the current `pose_normal_terms` kept, either test-local or as the slow path,
  and compared with the rewrite;
- half-plane sums compared with full-disk sums at the same pose;
- 1 thread against `omp_get_num_procs()` threads, 20 repeats, bit for bit.

## 8. Rulings (2026-10-19 and 2026-10-20)

- **The convergence changes (sections 3.1 and 3.3) are wanted for release 4,**
  for `refine=cont` and the polish alike (stage 4).
- **The benchmark is a `simple_test_exec` workflow outside the unit suite**
  (section 6.3, stage 2).
- **The gather is branch-free,** on a full, centred, band-limited reference
  (sections 4.4 and 4.5).
- **2026-10-20: stage 3 ships as variant D,** the pass removal and the
  `padsc` fold on the packed layout. The cube moves into stage 6a with the
  Friedel half-plane. Measured alone, it was 12% slower on beta-gal (stage 3
  record).

Settled by the stage 4 measurements (2026-10-20), for the maintainer to
confirm:

- **the stopping tolerances stay at `1e-6`,** with the predicted-reduction
  stop added (stage 4 record, finding 5);
- **euclid is corrected like cc,** by the `|f|^2`-weighted slope of the data
  on the model (finding 1);
- **the benchmark's capture regime bounds the rotation at 20 degrees,** so the
  bound no longer sits on the truth (finding 3). The gate's E25 keeps its
  15-against-15 design: it holds floors, and the comparisons are the
  benchmark's (stage 5).

- **2026-10-20: the gather is continuous.** The quadratic B-spline replaces
  the KB stencil on the same three taps, with the reference weighted so that
  on grid it is the polar branch's model (record after stage 5).
