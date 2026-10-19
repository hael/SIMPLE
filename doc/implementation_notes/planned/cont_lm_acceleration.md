# Accelerating refine=cont: fewer evaluations, a cheaper kernel, and a benchmark

Implementation note, 2026-10-19. **Stages 1 and 2 are implemented but not yet
compiled; stages 3 to 7 are proposals.** The maintainer's rulings of 2026-10-19
are in section 8. The findings come from a static review of `master` at
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

- `completed/cont_kb_gather_performance.md`: the fused gather and its list of
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
            -> gather_packed_kb3_window_grad 27-sample KB gather and 3 derivatives, per pixel
  -> euclid_sigma2%calc_sigma2 -> cartft_calc%sigma_contribution     (euclid only)
```

The search time of a particle is the number of full-disk passes it makes,
times the cost of one pass. Section 3 reduces the first factor and section 4
the second.

Done for this route so far:

- **Fused 3x3x3 KB gather.** Median worker orientation-search time -36.8 to
  -37.5%, wall time -20.5 to -21.9%, on beta-gal
  (`completed/cont_kb_gather_performance.md:142-148`).
- **Joint stage alone.** 37% fewer LM iterations and 38% less transaction time
  than `shift_then_joint`
  (`doc/refactoring_notes/completed/pose_cont_refactoring.md:1615-1618`).
- **Dynamic scheduling of the matcher's particle loop.** Commit `aea5b1846`,
  `simple_strategy3D_matcher.f90:235`.

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

Each stage is one commit, or a short series of commits. The agent builds the
stage and runs the checks it names, and reports the results. The maintainer
reviews and decides what is committed and pushed. The order is: measure first, then protect the kernel with tests,
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

**Status (2026-10-19): implemented and committed in `0c9542e50`, not yet
compiled.** Step 2.2 checks it. The changes:

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

**What was implemented (2026-10-19, not compiled):**

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
| Build B commit (stage 2) | `b2d8f34e5` |
| Machine, threads | FWM-S116118, Apple M2 Ultra (16 performance + 8 efficiency cores, 192 GB); timed runs `nthr=16`, alone on the machine |
| SNR | 0.05; polish also at 0.02 |
| `unit_cart_align3D` (B) | PASS, 761 of 761 assertions; fast gate 15/15 for A and for B |
| CR1 identity, A against B | Identical: the 60 check values of `metrics.tsv` are byte-identical, and both gates pass (`runs/cr1_A`, `runs/cr1_B`, about 40 minutes each at one thread) |
| Counter cost (beta-gal, A B B A) | Orientation search A 19.34, 19.41 s; B 19.33, 19.46 s (a third B run 19.12 s). The difference is inside the spread: no measurable cost. |
| Before-profile directories | `runs/profile_cr1_B/profile`: 20 of 20 seeds pass (table below). `runs/bench_polish_seed{1..5}`: to be retaken, because the pose-error metric was quantized (below). |

Simulated cases: 10 000 particles, box 256 at 0.73 A, reference reconstructed at
the truth poses. The polish band is shell 62 (lp 3 A) and the capture band
shell 23 (lp 8 A). Beta-gal: 4 199 active particles of 5 513, box 256 at
1.275 A, band shell 88 (3.7 A), D2, euclid.

| Case | Search thread-s | Search wall (s) | Passes / solve | LM its / solve | Solves at the cap | Rejected fraction | Actual/predicted mode bin (share) | ns / pixel-eval | Median rot (deg) | Median shift (A) | Mean corr_cart |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| polish_euclid, snr 0.05 | 381.6 | 24.1 | 14.88 | 10.88 | 0 | 0.41 | 1.25-1.75 (0.32) | 209 | 0.084 | 0.029 | 0.67 (`exp(-L)`) |
| polish_cc, snr 0.05 | 378.7 | 24.0 | 14.43 | 11.43 | 0 | 0.38 | 1.25-1.75 (0.38) | 214 | 0.088 | 0.030 | 0.78 |
| capture_euclid, snr 0.05 | 66.7 | 4.2 | 17.69 | 19.70 | 33 | 0.16 | 0.75-1.25 (0.47) | 218 | 0.250 | 0.095 | 0.83 (`exp(-L)`) |
| capture_cc, snr 0.05 | 61.1 | 3.9 | 15.71 | 19.91 | 15 | 0.16 | 0.75-1.25 (0.64) | 225 | 0.320 | 0.133 | 0.89 |
| polish_euclid, snr 0.02 | 375.6 | 23.8 | 14.78 | 10.78 | 0 | 0.42 | 1.25-1.75 (0.31) | 207 | 0.128 | 0.045 | 0.56 (`exp(-L)`) |
| polish_cc, snr 0.02 | 429.1 | 27.1 | 16.40 | 13.40 | 0 | 0.36 | 1.25-1.75 (0.51) | 213 | 0.131 | 0.046 | 0.65 |
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
0.0873 for euclid). The simulated truth metrics in the case table above, and the
benchmark seed runs, are therefore retaken with build B rebuilt at the commit
of this fix. The solver is unchanged, so the timings and work counts stand; the
beta-gal runs and the CR1 identity check do not use this metric.

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

- **3a. Duplicate passes and the `padsc` fold.**
  - Drop the duplicate seed and final passes (section 3.2, rows 1 and 2).
  - Fold `padsc` into the reference (section 4.4).
- **3b. Branch-free reference and gather.** Use the full, centred,
  band-limited reference with a branch-free gather (sections 4.4 and 4.5, per
  ruling). Keep today's accumulation order.
- **Results:** identical. Removed passes produce bitwise-equal values. The
  centred cube holds exactly the values the packed gather fetches today, as
  copies and exact conjugates. Scaling by 8 is exact.
- **Check:**
  - the CR1 identity check (section 6.2);
  - `unit_cart_align3D`;
  - the gather oracle of `simple_cartesian_fourier_tester` extended to the new
    layout;
  - the benchmark: two fewer passes per particle (3a), a lower cost per
    pixel-evaluation (3b), and reference memory as in section 4.5.

### Stage 4: converge in fewer iterations

These changes alter results and were ruled wanted for release 4, for
`refine=cont` and the polish alike.

- **What:**
  - scale the CC Gauss-Newton matrix by the correlation (section 3.1);
  - the stopping rule (section 3.3).

  The euclid matrix changes only if the baseline shows its over-step
  signature (section 3.1). At realistic SNR the stopping rule is the larger
  lever: in the smoke runs, 40-45% of trial evaluations were rejected.
- **Check:**
  - a robustness profile on truth, before and after: `test=profile` on
    `cont_refine3D_1jxy`, and the benchmark simulation over its seeds;
  - the benchmark's iterations per particle and its time;
  - beta-gal FSC within the run-to-run spread.
- **Done when:** the tolerance values are set from the profile, and the truth
  metrics are no worse than the before-profile's spread.

### Stage 5: kernel oracle and test gaps

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

### Stage 6: cheaper evaluations

- **6a. One rewrite of `pose_normal_terms` and its siblings.** It covers the
  Friedel half plane, loop bounds set by `kto`, the separable shift phase and
  the normal-term arithmetic (sections 4.1, 4.2, 4.3 and 4.5).
- **6b. The KB weights** (section 4.4), plus any SIMD across taps that
  reorders sums.
- **6c. sigma2 from the accepted evaluation** (section 3.2, row 3).
- **Results:** rounding changes.
- **Check:**
  - the stage 5 oracles within tolerance;
  - `unit_cart_align3D`;
  - the truth profile unchanged within its spread;
  - the benchmark's ns per pixel-evaluation and passes per particle.

### Stage 7: close-out

- Update the Cartesian section of `refine3D_policy.md` to cover the objective
  scaling, the stopping rule, the counters and the reference layout.
- Fix the stale `production/CMakeLists.txt:243-245` comment.
- Record the measured gains of every stage in this note, and move the note to
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

**`objfun=euclid` is different.** It is the production objective of
`refine3D_auto` and of the proposed `refine3D_local`.

- `L = sum |m - x|^2 / P`, with `H_GN = 2 J^T J / P`. This is right on average
  when the reference amplitude matches the data.
- If the reference is `beta` times the data signal, then `H_GN/H ~ beta`.
- A filtered reference (`beta < 1` at high k) therefore over-steps and pays in
  rejections and damping.

Stage 2 measures the ratio and the rejection rate. A euclid change follows
only if they show over-stepping.

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

**KB weights (stage 6b).**

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

**Indexing (stage 3b; ruled 2026-10-19: a branch-free gather).**

- Today each sample costs three wrap-table lookups, the packed/Friedel branch
  and two conditional index adds (`simple_cartesian_fourier.f90:183-211`).
- On the reference of section 4.5, each 3-tap run along x is three contiguous
  complex values at a fixed offset. There is no wrap, no Friedel branch and no
  conditional index.
- Keeping today's accumulation order keeps the result bit-identical.
- Vectorizing across the three taps reorders the sums, so it belongs to
  stage 6b.

### 4.5 Loop bounds and the reference layout

**Loop bounds (stage 6a).** The pixel loops visit the full `(box+1)^2` square and
test `shells` (`:677-680`). Bound them by `kto`, and by `h >= 0` under section
4.1.

**Reference layout (stage 3b; ruled 2026-10-19).** The reference becomes a full
(unpacked) complex cube.

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
(`completed/cont_kb_gather_performance.md:101-106`):

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
| Route under test | Production defaults | Stages a-c and the library suite force `shift_then_joint` (`simple_commanders_test_highlevel_cont.f90:495`, `simple_strategy3D_cont_1jyx_tester.f90:464`); production uses `joint` |
| Fast-path oracle | Replaced routine kept and compared | Only the gradient is checked, against finite differences; the 5x5 matrix and `shift_normal_terms` have no independent check |
| Threads | 1 against `omp_get_num_procs()` threads, 20 repeats, bit for bit (speed2D rule) | A fixed 3-thread team, run once |
| Negative control | Required for every gate | None in `cont_refine3D_1jxy` |
| Reference scale | The objective's assumptions met | The gate's euclid stages (a, c) run against the pdb2mrc truth map, whose amplitude scale does not match the noise-normalized particles. Their median rotation error is 0.94 degrees, against 0.23 degrees for cc (stage b) on the same data (before-profile, stage 2). The benchmark reconstructs its reference for this reason. |
| Small-angle metric | Resolves the tested differences | The gate's `acos` of a single-precision trace resolves about 0.002 degrees near zero; enough for its floors near 1 degree, not for comparing stages |
| Large tester objects | `class(T), allocatable` | Plain `type(parameters)`/`type(builder)` locals (`simple_strategy3D_cont_tester.f90:75-76`) |
| Comments | Current | `production/CMakeLists.txt:243-245` says 5,000 particles and kept TSVs; the code uses 1,000 and removes its files (stage 7) |

The oracle checks stage 6 needs:

- an independent `J^H J`, built from per-pixel centred finite differences of the
  predicted model, compared with the 5x5 matrix of `pose_normal_terms`, for cc
  and euclid under all three CTF flags;
- the current `pose_normal_terms` kept, either test-local or as the slow path,
  and compared with the rewrite;
- half-plane sums compared with full-disk sums at the same pose;
- 1 thread against `omp_get_num_procs()` threads, 20 repeats, bit for bit.

## 8. Rulings (2026-10-19)

- **The convergence changes (sections 3.1 and 3.3) are wanted for release 4,**
  for `refine=cont` and the polish alike (stage 4).
- **The benchmark is a `simple_test_exec` workflow outside the unit suite**
  (section 6.3, stage 2).
- **The gather is branch-free,** on a full, centred, band-limited reference
  (sections 4.4 and 4.5, stage 3b).

Still open, to be decided from measurements:

- the stopping tolerances (the stage 4 profile);
- how the euclid matrix is corrected. The stage 2 baseline shows that it over-steps on beta-gal
  (64% of accepted steps below 0.75), so stage 4 covers euclid too; the form of the correction is open.
