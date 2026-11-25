# Reproducible phase-randomized FSC

Status, 2026-11-25: complete and ready for the maintainer's commit review. The maintainer reported Hans's approval on 2026-11-25; decisions 2–4 are accepted as implemented and Phase 5 is closed. Phase 0 baseline, Phase 1 transfer, Phase 2 GCC16.2 Debug validation, Phase 3 macOS GCC15.2 Debug/Release qualification and both tester follow-ups are complete. Phase 4 is complete by the explicitly authorized reuse of verified research experiments and source comparison. Base: `89b1b4817`; research source: `07afabbbf` plus the preserved eight-file patch. The final scope is nine Fortran files, one test-policy document and this completed note; the original eight-file transfer is historical scope. Applied recovery stashes were verified and dropped with user approval; all validation logs and the unfixed Phase 0 source/builds remain retained. Approval of the implementation does not authorize staging, committing or pushing; commit and push each require separate user approval.

This is the single living implementation record, including the final report below. There are no open implementation phases or acceptance placeholders. Reported scientific limits and work outside this fix remain explicit. The format guide is `doc/how2s/how2write_agent_plans.md`; the report is kept here under the repository's single-note rule.

## 1. Why

Phase randomization used persistent OpenMP worker RNG streams. Reseeding the calling thread could therefore change the corrected FSC curve on identical half-maps, shifting its resolution cutoff and the reference supplied to nonuniform filtering. The original GCC15.2 baseline shows it: three repeats at 1/2/4/16/48 threads give three passing serial cases, and all twelve multithreaded cases reproduce the same-seed randomized/corrected-FSC failure. The fix generates one serial per-call seed, then keys Fourier phases by that seed, state, half and coefficient coordinates. Replaying the same seed and call sequence no longer depends on worker RNG history.

This change concerns FSC repeatability. It does not change the FSC correction equation, mask, onset or threshold policy, average multiple realizations, or modify an optimizer. Different seeds still give different phase realizations. Hans's suspected sporadic LM failure remains a separate investigation.

## 2. Contract

1. Reproduce same-seed randomized/corrected curves independently of worker RNG history and thread count, within the tested numerical tolerance (`TOL` = 2e-6, the existing tolerance). Preserve raw/masked curves and input maps/mask.
2. Retain the existing seed/reseed lifecycle. Generate a fresh realization on consecutive calls without reseeding; calls without an explicit run seed remain stochastic. This preserves reseeding behavior, not pre-fix RNG stream positions: an active FSC randomization now consumes one calling-thread intrinsic draw instead of per-coefficient calling-thread/worker draws. Subsequent unkeyed RNG consumers can receive different values, even under the same run seed. Same-seed reproducibility requires the same call sequence; identity with the old code or across arbitrary call-sequence changes is not promised.
3. Reserve a separate FSC scope in the existing keyed RNG (`RND_PROG_FSC` = 3, scope `rnd_scope(RND_PROG_FSC,RND_PASS_DIAG)` = 52). Even and odd halves use different keys; phases remain uniform, amplitudes remain unchanged, and coefficients below the onset remain untouched.
4. Keep the new `image%phase_rand` key optional. Other callers retain their existing path. No new CLI control, CMake target or CTest entry.
5. Off-state invariant: with `phase_key` absent, `phase_rand` runs only the original `ran3() * TWOPI` branch. The callers in `simple_commanders_imgops.f90` and `simple_procimgstk.f90` do not pass a key and are unchanged. Negative control: test "legacy no-key serial phase randomization remains repeatable".
6. Preserve current master changes and all existing acceptance tolerances.
7. Seed determinism requires `random_number` to be called from a serial context after `seed_rnd()`. `phase_rand_fsc` is not to be called from inside a parallel region. Callers: `simple_halfmap_diagnostics.f90` (no OpenMP in the file, checked in review), `simple_commanders_resolest.f90` (checked in the Phase 0 caller audit), and the serial-context test calls in `simple_fsc_tester.f90`.

Review points for the maintainer: the serial 53-bit seed conversion, scope separation, half identities, logical Fourier coordinates (h ≥ 0 in the stored half; negative k and l enter as two's-complement low words in `rnd_key_lanes`) and the OpenMP private key (`coefficient_key`).

## 3. Decisions and rulings

1. (2026-11-25, user) Transfer the verified eight-file patch into the master working tree for manual review only. Only the eight listed source/test files and this note belong to the transfer; the pre-existing `.gitignore` edit (`graphify-out/`) is unrelated. No optimizer changes, no new HPC campaign, no commit or push.
2. (Hans approved as implemented; approval reported by the maintainer on 2026-11-25) One serial seed per `phase_rand_fsc` call from `random_number`, converted as `int(u * 2^53, int64)`, so consecutive calls without reseeding give fresh realizations.
3. (Hans approved as implemented; approval reported by the maintainer on 2026-11-25) Keys: `scope` = FSC scope, `reference` = state, `family` = 1 for even and 2 for odd, `iteration`/`particle`/`ordinal` = logical h/k/l. The accepted implementation keeps numeric half families 1 and 2 in the separate FSC scope; their values coincide with `RND_FAM_SAMPLE` and `RND_FAM_INPL` in alignment scopes. Adding named FSC family constants is not required for this approved change.
4. (Hans approved as implemented; approval reported by the maintainer on 2026-11-25) Same-seed FSC curves now differ from pre-fix runs, so this is result-changing (section 5). Different seeds keep legitimate sampling variation; no threshold or onset rule changes.
5. (2026-11-25, plan) Failures are preserved; no tolerance is widened to obtain a pass. Rejected evidence (earlier failed-launcher and wrong-input-stage attempts) stays preserved and excluded from scientific acceptance.
6. (2026-11-25, plan) Prior research qualification does not qualify a newer integration target. If affected source changes again, retrace and rerun the relevant checks.
7. (2026-11-25, user) Reuse the completed research experiments instead of repeating the HPC campaign for this FSC fix. This amends decision 6 only for Phase 4's real-data evidence: preserve the experiments' original source/build identities, connect them to the unchanged FSC implementation by source comparison, and use Phases 0 and 3 for current-target local validation. Do not relabel historical jobs as current-master or current-snapshot Linux runs. If the relevant numerical path changes, reassess reuse.
8. (2026-11-25, scope clarification implementing decision 7) The frozen five-seed study is a measured diagnostic of seed variation, not a statistical acceptance band. Its original purpose explicitly declares no inference floor or quality verdict. Replace the unexecuted before-profile/sample-size placeholder with that descriptive evidence; do not invent a retrospective pre-registered pass rule. Same-seed identity, raw/masked invariance and frozen-input checks remain required and were met. Duplicated LM pairs and the sporadic-LM investigation are outside this plan.
9. (2026-11-25, user) Correct the FSC tester's seed-offset handling, remove fast-tier timing and assert real OpenMP coverage; document the RNG stream-consumption impact and rerun GCC15.2 Debug/Release validation. Production behavior and files remain unchanged. Use the repository-convention fixed three-thread team, independent of machine/core settings, and preserve the caller's OpenMP settings. At that stage this tester follow-up did not approve production decisions 2–4; their subsequent approval is recorded in decision 11.

10. (2026-11-25, user) Apply the tester, documentation, comment and style review fixes; retain `real(TWOPI,dp)`, field-by-field keys, fixture noise and the absent optional parallel-region guard. Rerun GCC15.2 Debug/Release and report only the delta since manual review. This authorizes no production numerical change, commit or push.

11. (2026-11-25, maintainer report of Hans’s approval) Accept production decisions 2–4 and the implementation as reviewed. The maintainer moved this note to `completed/` and authorized documentation closure after the final review. Phase 5 is complete with the final report and proposed commit message below. This approval does not authorize production changes, staging, commits, pushes or baseline cleanup.

Defaults taken without a question: test suite name `FSC phase repeatability` with selector `fsc_phase_repeatability` in the existing `unit_reconstruction` area; tester tolerances 1e-5 (random-phase difference between seeds), 1e-4 (distinct key realizations) and 2e-6 (repeatability); the field names `phase_key` and `coefficient_key`.

## 4. Not in this plan

- Optimizer changes, the FSC correction equation, mask, onset or threshold policy, averaging of multiple realizations.
- The sporadic LM failure investigation, including the five duplicated LM pairs. The completed ten-run cross-optimizer pilot is supplementary context and does not resolve that separate bug.
- New CLI control, CMake target or CTest entry; research probes, HPC launchers, datasets, compiled files and private development history.
- Commit, push and staging: the maintainer decides.
- Final map quality and the threshold rule. The five-seed study uses one newly frozen bootstrap, not the unavailable historical bootstrap, and does not establish either.
- Global seeds also affect other draws than FSC.
- Queue check closed at final local review, 2026-11-25: among the other registered worktrees, `codex-lm-diagnostics` carries the same FSC production implementation in `simple_fsc.f90`, `simple_image_filt.f90` and `simple_rnd.f90`, differing only in spacing and the FSC field-map comment. The other inspected worktrees have no dirty changes to those three files. No competing numerical change or scope retirement was found. This is a local source/worktree check, not a claim about unobserved remote activity; later integration must preserve the shared fix without applying it twice.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 Baseline on current base | tests-only (measurement) | complete locally: GCC15.2 Debug/Release builds, installations and 15/15 gates; six serial reproducer passes, twelve expected multithreaded failures; production spread and caller audit recorded below |
| 1 Recover and transfer | result-changing (code) | complete: eight files equal the research patch apart from whitespace |
| 2 Validate GCC16 Debug | result-changing | complete: fast gate, reverse reconstruction, three threaded FSC repeats |
| 3 GCC15.2 and Release qualification | result-changing | complete locally on macOS: GCC15.2 Debug/Release builds, 15/15 gates, install, latest cleanup: 2,009/2,009 reverse reconstruction, 996/996 reverse image and five 55/55 threaded FSC repeats per build |
| 4 Reuse verified real-data evidence | result-changing evidence; no new code | complete by decisions 7–8: four same-seed pairs have identical saved FSC curves/final-map hashes; corrected five-seed probe meets repeat/input invariants; current FSC source comparison recorded |
| 5 Maintainer review | reported, not gated | complete: Hans’s approval reported by the maintainer on 2026-11-25; final review and report complete; commit/push remain unauthorized |

## 6. Acceptance

Tags: *measured*, *existing* (a tolerance already in the repository), *guessed*.

| Quantity | Value | Tag |
| --- | --- | --- |
| Same-seed repeat difference (unmasked, masked, randomized, corrected) | within `TOL` = 2e-6, printed differences zero | existing |
| Raw and masked curves across seeds | within `TOL` | existing |
| Frozen APOF probe vs native bootstrap, raw/masked curves | 5.97e-8 observed against 2e-6 | measured |
| Different-seed randomized-FSC difference in tester | greater than 1e-5 | guessed |
| Distinct half/key realizations in tester | greater than 1e-4 of amplitude | guessed |
| Amplitude preservation in tester | 2e-6 relative | existing |
| Fast-gate time | advisory 60 s budget; 138.89 s observed on GCC16 Debug | existing budget, measured |
| Build resources | at most five compile jobs, tests serial with at most four OpenMP threads and BLAS1; build-process RSS 0.33 GiB sampled against a 12 GiB laptop ceiling | existing |

Identity cases: same seed and call sequence, one, two, four, sixteen and 48 threads (research) and actual three-thread teams (latest current-target tester), with unchanged inputs byte for byte. Exclusions: none declared; different seeds are expected to differ.

Measured seed variation, reported descriptively under decisions 7–8: the frozen APOF bootstrap study has five seeds plus one repeated seed (N = 5). Maximum per-shell ranges are 0.042057 for randomized FSC and 0.028740 for corrected FSC. FSC0.5 is 4.378947 Å for all five seeds. FSC0.143 and the NU-bank handoff are 3.670588 Å for seeds 20261123/26/27 and 3.466666 Å for 20261124/25. Shell 69 straddles 0.143, explaining the four-shell cutoff change. The repeated seed is byte-identical; raw/masked curves are seed-invariant; frozen inputs are unchanged. The observed variation is distinct from the same-seed defect. N = 5 does not establish tail probabilities, a general uncertainty band or final-map quality. No retrospective statistical gate or pre-fix-curve identity is claimed.

Reported, not gated (maintainer reads): APOF LM cFAR changed from 0.9675 to 0.9596 between historical and fixed runs. Historical final resolution endpoints were unchanged, but source snapshots and random realizations also differed; this is an unresolved quality comparison, not evidence of a regression-free scientific result. The controlled synthetic baseline/fix comparison independently establishes the RNG defect and its correction. Artifacts: local `codex_logs/fsc_phase_repro_20261124/` (including `report.md` and `realdata_review/evidence/`) and `codex_logs/fsc_master_transfer_20261125/`; `codex_logs/` is gitignored, so these are not in the repository. The final report below preserves the essential findings in this repository note and explicitly identifies the underlying artifacts as local evidence.

## 7. Common rules and stop conditions

Rules: no new parameter; no floor, tolerance or assertion relaxed; fix as you go with a file-table row first; compilation at most five jobs; tests serial with at most four OpenMP threads and BLAS1; no full particle refinement on the laptop; research probes and launchers are excluded.

Stops: unexpected build/test failure, changed test inputs, nonfinite curves, raw/masked seed dependence, a differing same-seed control, or a file outside the table. Preserve failures; do not widen tolerances to obtain a pass.

Gates no ruling may waive: the same-seed control, raw/masked seed invariance, input isolation, and the off-state (no-key) test. Waiving GCC15.2 qualification on this exact snapshot is Hans's explicit decision, not an agent ruling.

## 8. Phases

- **Phase 0, baseline on the current base (complete locally, 2026-11-25).** On `89b1b4817` plus the dirty tree minus the eight files (or a base build kept for comparison): build Release and Debug and record the machine (compiler, cores, load); run the fast gate; run the new FSC tester against the unfixed `phase_rand_fsc` to record the same-seed failure on this snapshot (the twelve research failures are on older GCC15.2 builds); run the production entries that touch FSC twice with a fixed seed (`simple_commanders_resolest.f90` and `simple_halfmap_diagnostics.f90` paths); check by `rg` that `phase_rand_fsc` is not called inside a parallel region at the pinned commit, including `simple_commanders_resolest.f90`; confirm the command lines against the UI's required keys. Exit: build logs, a gate log, a failing-baseline log for the FSC tester, two-run spread of the fixed-seed entries, and the caller check recorded.
- **Phase 1, recover and transfer (complete).** Recover the verified eight-file patch, inspect current master, preserve unrelated edits and adapt formatting-only conflicts (three files). Leave the code and this note uncommitted for manual review. Exit: eight files changed as in the table; the exact eight-file change matches the research patch apart from whitespace; original patch SHA256 `f8b1a45462604f692961a59ac5f978237b20eeb8499a6bac2e9a64fe84b62a65`; nothing staged, committed or pushed (met).
- **Phase 2, validate the transferred target on GCC16.2 Debug (complete).** Incremental build/install, all fifteen fast-gate suites, reverse reconstruction, three fresh-process FSC repeats. Exit: all pass (met; numbers in Progress).
- **Phase 3, GCC15.2 and Release qualification (complete locally on macOS, 2026-11-25).** Build and install on GCC15.2 and as Release; run the fifteen-suite gate, `unit_reconstruction suite=fsc_phase_repeatability`, and `unit_reconstruction` reverse reconstruction; three fresh-process threaded FSC repeats; record the gate time. Exit: all suites pass with 0 failed assertions on both builds, the gate time is recorded (and, if over the 60 s budget, whether the overrun predates this change), and the logs are named in Progress. If the maintainer waives this phase, that ruling is recorded here by number.
- **Phase 4, reuse of completed real-data experiments (complete, 2026-11-25; decisions 7–8).** Reuse the preserved beta-gal/APOF C/R and LM L/L2 campaign (eight parents, four same-solver repeat pairs) and the corrected frozen APOF pre-ML probe (five seeds plus one repeated seed). Verify saved-shell comparisons for iterations 1–3 and final reconstruction, paired final-map hashes, probe repeat identity, raw/masked seed invariance and frozen-input records. Record the original source/build identities and the comparison to the current FSC implementation. Treat seed spread descriptively as section 6 specifies. Exit: preserved evidence and source comparison support those invariants, current-target validation remains recorded in Phases 0/3, and the reuse decision and limits are explicit (met; Progress below). No new HPC campaign or duplicated-LM study is part of this phase.
- **Phase 5, maintainer review (complete, 2026-11-25).** The maintainer reported Hans’s approval of the reviewed implementation, closing decisions 2–4, and moved this note to `completed/`. Final source/policy review found no remaining code defect, architecture violation or new scientific decision requiring approval. Documentation closure and the final report, including a proposed commit message, complete this phase. Commit and push still require separate explicit authorization.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/utils/math/simple_rnd.f90` | Reserve the FSC RNG scope (`RND_PROG_FSC`). | 1 |
| `src/main/image/simple_image.f90` | Declare the optional phase key. | 1 |
| `src/main/image/simple_image_filt.f90` | Generate keyed phases when a key is supplied. | 1 |
| `src/utils/filter/simple_fsc.f90` | Generate the per-call seed and separate half keys. | 1 |
| `src/utils/math/simple_rnd_tester.f90` | Check scope separation alongside existing RNG vectors. | 1 |
| `src/utils/filter/simple_fsc_tester.f90` (new) | Test correction math, seed/thread repeatability and input isolation; seed-offset, timer removal, team coverage, canonical seed cleanup and named-test follow-ups; direct image contracts now live in the image tester. | 1, follow-up |
| `src/main/commanders/test/simple_commanders_test_class.f90` | Register the tests in the existing reconstruction area. | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | Add the suite selector. | 1 |
| `src/main/image/simple_image_tester.f90` | Own the direct keyed/no-key image phase contracts, with a deterministic impulse fixture. | review follow-up |
| `doc/policies/test_environment_policy.md` | Add FSC phase repeatability to the reconstruction suite inventory. | review follow-up |
| `doc/implementation_notes/completed/fsc_phase_randomization_transfer_plan.md` | This note. | all |

Formatting follow-up applied on 2026-11-25: align the test-class import, procedure imports and key assignments in `simple_fsc.f90`, and phase-key declarations in the image interface/filter. Production executable tokens are unchanged from the manually reviewed snapshot.

## Progress

**Research, 2026-11-23/24 (earlier hosts, GCC15.2 and GCC16 research builds).** Source `07afabbbf` plus the eight-file patch.

| Test | Observed result |
| --- | --- |
| Original GCC15.2 baseline: three repeats at 1/2/4/16/48 threads | Three serial cases pass; all twelve multithreaded cases reproduce the same-seed randomized/corrected-FSC failure. |
| Fixed GCC15.2 build `55697434` | All fifteen gate suites and all fifteen FSC cases pass. |
| Fixed local GCC16 Debug research build | All fifteen gate suites pass; reverse reconstruction passes 2,009/2,009 assertions. Gate time 143 seconds exceeds the advisory 60-second budget. |
| Same-seed beta-gal/APOF C/R and LM L/L2 repeats | All eight parents and 240 workers succeed. Every saved FSC shell in iterations 1-3 and final reconstruction matches within each repeat pair; final maps are byte-identical. |
| Frozen APOF bootstrap, five seeds plus one repeated seed | Repeated seed is byte-identical; raw/masked curves are seed-invariant; frozen inputs are unchanged. |

The frozen APOF probe used production's unfiltered pre-ML half-map pair and density mask. Its raw/masked curves agree with the native bootstrap within 5.97e-8 (existing tolerance 2e-6). Bootstrap `55697680` and corrected probe `55697779` succeeded. Earlier failed-launcher and wrong-input-stage attempts are preserved and excluded from scientific acceptance. Seed-variation numbers are in section 6.

**Phase 1, 2026-11-25, master working tree.** Base `89b1b4817`. Eight files changed as in the file table; master formatting changes required adapting three files; the transferred algorithm and tests are unchanged and the change matches the research patch apart from whitespace, preserving unrelated current-master content. Nothing is staged, committed or pushed. Authorized by the user on 25 November.

**Phase 2, 2026-11-25, local GCC16.2 Debug.** Incremental build/install and all fifteen fast-gate suites pass (138.89 seconds; advisory 60-second warning). Reverse reconstruction passes 2,014/2,014 assertions in 23 groups (the research build had 2,009, a different source target). Three fresh-process FSC repeats each pass 60/60 checks with an actual four-thread team; all printed same-seed differences are zero and saved curves agree across processes. Compilation used at most five jobs, tests ran serially with at most four OpenMP threads and BLAS1; sampled build-process RSS was 0.33 GiB, below the 12 GiB laptop ceiling. No full particle refinement ran locally. Logs: `codex_logs/fsc_master_transfer_20261125/` (`build.log`, `install.log`, `fast-gate.log`, `fsc-repeat-{1,2,3}.log`, `reverse-reconstruction.log`, `registry.log`, `validation.json`). Meaning: the fix holds on this snapshot under GCC16 Debug; GCC15.2 qualification was outstanding at that stage and was subsequently completed in Phase 3 and the follow-ups below. Historical GCC15.2 evidence remains in the research table above.

**Phase 0, 2026-11-25, local macOS baseline (complete).** User authorized stashing all current master changes, preserving logs and completing the local baseline with GCC15.2. Master was clean throughout Phase 0 at `89b1b4817c2adfc1132210b742973516adb7e25d`. Fix/plan/untracked tester/previous `.gitignore` edit: stash `cdb4e26811f868e8795dd1b236bd659a7c086410`. The graph outputs exposed when `.gitignore` was stashed are separately preserved in stash `338dad6dd1b1a79571a2d1704676a2866cc75e68`. Both original validation directories remain in place; an archive of all 5,186 files was verified by SHA256 and retained with this phase's evidence.

Baseline worktree: `/private/tmp/simple-fsc-phase0-20261125/source` (detached at the exact base). Machine: macOS 26.7.3 arm64, 11 logical CPUs, 18 GiB RAM. C, C++ and Fortran all report Homebrew GNU 15.2.0 from `/opt/homebrew/opt/gcc@15.2.0/bin`. Independent `build_debug` and `build_release` trees are siblings of the baseline source; the pre-existing master GCC16 Debug tree was preserved. No new compiler or package was installed. Existing dependency archives were used to build compiler-specific OpenBLAS/gtk-fortran caches. Compilation used at most five jobs; tests ran one process at a time, OpenMP capped at four, BLAS at one. Exact compiler paths, flags, source hashes and commands are retained.

| Baseline check | Observed result |
| --- | --- |
| Untouched Debug build, fast gate, install | PASS; 15/15 suites, 141.51 s; advisory 60 s warning |
| Untouched Release build, fast gate, install | PASS; 15/15 suites, 51.94 s |
| Baseline-compatible FSC reproducer: Debug and Release, 1/2/4 actual threads, three fresh processes each | 6/6 serial cases pass 49/49 assertions; all 12 multithreaded cases fail exactly the two expected randomized/corrected repeatability assertions (47/49 pass, nonzero exits preserved) |
| Serial and explicit worker-stream-reset controls | All printed curve differences zero in all 18 cases |
| Raw/masked invariance, correction equation, finite curves, input isolation and no-onset fallback | All corresponding reproducer assertions pass |
| Production `fsc automsk=density` and `nu_filt3D automsk=order` (half-map diagnostics path) | Two fresh processes per entry at 1/4 threads in each build: all 16 normal stops; saved curves finite, original input SHA256 unchanged |
| Production paired-curve spread | All one-thread pairs match. Debug four-thread `fsc`: max randomized delta 0.05799585, corrected delta 0.04691760; Debug four-thread NU randomized delta 0.05758105; Release four-thread NU randomized delta 0.05855625. Release four-thread `fsc` pair happens to match; retained without rerunning to seek a failure. All raw/masked paired deltas are zero. NU CLI saves raw/masked/randomized curves; it does not save the corrected curve. |
| Production caller audit | `exec_fsc` and the gridding `evaluate_halfmap_pair` paths run outside inspected OpenMP regions. Volume assembly uses a sequential state loop; nearby parallel restoration loops end before diagnostics. PCG diagnostics do not call phase randomization. Pinned caller-chain references/hashes retained. |

Method clarification: first build/gate/install the completely untouched base. Then add the original baseline-compatible research tester (with optional fixture export only) and its registration/UI selector in this isolated worktree. This avoids depending on the stashed fix's new RNG constant or keyed image interface. Only the three already-listed test/registration files are instrumented; all four production/interface files remain byte-identical to the base. The production CLI checks use the pristine installed executables and the same exported 48-cubed synthetic half-map pair (`smpd=1`, `mskdiam=34`, seed 20261124). Required UI keys and both numerical paths were checked against source and exercised successfully. This is a small-map local baseline, not a real-data refinement or scientific-quality comparison.

One environment deviation: the first full Debug build reached 100%, then `/usr/bin/time -l` failed because the sandbox blocked `sysctl kern.clockrate`. The entire log was preserved; only the resource-reporting wrapper was replaced by elapsed-time measurement. An incremental build then succeeded before any gate or instrumentation. No production source or tolerance was changed to obtain a pass.

Evidence: `codex_logs/fsc_phase0_local_20261125/` in the canonical checkout (local, gitignored), including `phase0-results.json`, `measurements.json`, `acceptance-before-measurements.json`, `compiler-versions.txt`, `build-configuration.json`, `caller-audit.json`, `harness.patch`, complete `logs/`, all 18 synthetic logs and all 16 production run directories. At Phase 0 completion the source worktree retained only the measurement harness and an execution-note copy; the transferred fix was still stashed. The later restoration is recorded below. Meaning: the unfixed exact base builds and passes its existing gates on macOS GCC15.2, while the controlled reproducer confirms the RNG defect in both Debug and Release. This does not qualify the stashed fix on GCC15.2 or replace Linux/HPC or real-data gates. No main-repository commit, push, server operation or full particle refinement occurred.

**Phase 1 restoration, 2026-11-25 (complete).** Following the user request after Phase 0, applied fix stash `cdb4e26811f868e8795dd1b236bd659a7c086410` to master and restored graph-output stash `338dad6dd1b1a79571a2d1704676a2866cc75e68`. Both stashes were retained at restoration; their later verified cleanup is recorded below. Verified all eight restored source/test files byte-for-byte against their stash contents, including the untracked FSC tester; the unrelated `.gitignore` edit is also restored. Updated this master plan with the observed Phase 0 results above. No numerical source changes were made during restoration; no commit or push occurred. Phase 1's recovery/transfer exit remains satisfied, and the earlier Phase 2 GCC16 Debug evidence applies to the identical restored source.

**Preparation for Phase 3, 2026-11-25.** Phase 0 is already logged above with its build/gate logs, 18 reproducer cases, 16 production invocations and unchanged-input checks. Before stash cleanup, compared all eight tracked changed files (including `.gitignore`), the untracked FSC tester and all 986 graph-output files byte-for-byte against their stash objects. The only intentional difference is this updated plan; preserved its original stashed version with the verification receipt. With explicit user authorization, dropped only applied FSC stash `cdb4e26811f868e8795dd1b236bd659a7c086410` and graph stash `338dad6dd1b1a79571a2d1704676a2866cc75e68`; unrelated Cartesian stash `83e1907dbc46a636456c35ada472af3f4556f68b` remains intact. Receipts: `codex_logs/fsc_phase3_local_20261125/stash-verification-before-drop.json`, `stash-cleanup.txt`, `plan-original-stash.md`.

User authorized local Phase 3. Frozen SHA256 manifest covers 1,219 candidate source/build files, including the untracked FSC tester. Exact GCC15.2 C/C++/Fortran paths are explicit. New trees: `build_gcc15_debug/` and `build_release/`; preserve existing GCC16 `build/` and all Phase 0 source/build/evidence trees. At most five compile jobs; one test process at a time, OpenMP capped at four and BLAS at one. Required exits are both fifteen-suite gates, installation after each passing gate, reverse `unit_reconstruction`, and three fresh-process FSC repeats per build with 60/60 assertions, actual four-thread teams and identical saved full-shell records. No source or tolerance changes, HPC operations, commits or pushes are part of this execution. Status at execution start: running; superseded by the completed Phase 3 results below. Logs and live stage: `codex_logs/fsc_phase3_local_20261125/` (`stage.txt`, `results.json`, `driver.log`).

**Phase 3, 2026-11-25, local macOS GCC15.2 fixed-code qualification (complete).** Supersedes the running status above. Built the exact restored candidate with matching GCC/G++/GNU Fortran 15.2.0 frontends, in independent `build_gcc15_debug/` (Debug) and `build_release/` (Release) trees. The existing compiler-specific OpenBLAS and gtk-fortran caches were reused. The source SHA256 manifest of 1,219 files matches before and after all validation. No source, tolerance, compiler policy or CTest registration was changed during this phase.

| Fixed-code check | GCC15.2 Debug | GCC15.2 Release |
| --- | --- | --- |
| Full build | PASS | PASS |
| Fast gate | 15/15 suites, 141.86 s | 15/15 suites, 53.30 s |
| Installation after passing gate | PASS | PASS |
| Reverse-order `unit_reconstruction` | 2,014/2,014 assertions, 0 failures | 2,014/2,014 assertions, 0 failures |
| Three fresh-process `suite=fsc_phase_repeatability nthr=4` runs | 60/60 assertions each, actual four-thread teams | 60/60 assertions each, actual four-thread teams |
| Serial, ordinary-reseed and worker-stream-reset controls | All printed same-seed deltas zero | All printed same-seed deltas zero |
| Saved full-shell FSC records across the three processes | Equal (96 records per run) | Equal (96 records per run) |

The Debug gate's advisory 60-second overrun predates the fix: unfixed Phase 0 was 141.51 s, fixed Phase 3 was 141.86 s. Release was 51.94 s before and 53.30 s after; neither exceeds the budget. These times document the gate budget, not a performance claim. No failed stage or assertion was observed in Phase 3. Builds used at most five jobs; tests ran serially with OpenMP capped at four and BLAS at one.

Baseline preservation: the four numerical/interface source files and four installed executable hashes in the Phase 0 tree remain unchanged. The existing GCC16 `build/` was not reconfigured or rebuilt. All 986 restored graph-output files were verified again after qualification. The two applied FSC/graph stashes are removed after verified restoration; the unrelated Cartesian stash remains intact. The master index is empty, HEAD remains `89b1b4817`, and no commit, push, HPC operation or full particle refinement occurred.

Evidence: local, gitignored `codex_logs/fsc_phase3_local_20261125/`. `results.json` contains every command, exit code, elapsed time, assertion count, thread measurement and retained shell records. Logs: `configure-{debug,release}.log`, `build-{debug,release}.log`, `fast-gate-{debug,release}.log`, `install-{debug,release}.log`, `reverse-reconstruction-{debug,release}.log`, and `fsc-{debug,release}-repeat{1,2,3}.log`. Exact compiler/cache identities and flags are in the three `*-version.txt` files and `CMakeCache-{debug,release}.txt`; source identity is in `candidate-source-sha256.json` and `candidate-before.patch`, with the new tester preserved separately. Final checks: `verification-final.json`, `baseline-preservation-check.json`; stash receipts are listed in Preparation above.

Meaning: the fixed code satisfies Phase 3's local macOS GCC15.2 Debug and Release gates, without relaxing any acceptance criterion. This qualifies the fixed candidate rather than the unfixed baseline. No current-snapshot Linux result or real-data scientific-quality result is claimed.

**Phase 4, 2026-11-25, reuse of completed real-data evidence (complete).** User agreed to reuse the completed experiments rather than repeat them, and excluded duplicated LM pairs as part of a separate bug investigation. Decisions 7–8 supersede the earlier requirement to rerun on this base and judge a new seed-spread before-profile. The numerical invariants below are retained; no failed assertion or tolerance was changed.

| Reused evidence | Original source/build | Verified result and scope |
| --- | --- | --- |
| Same-seed beta-gal/APOF C/R and LM L/L2 campaign | `a4788cb44` plus keyed-FSC patch; L-BFGS-B qualification `55692551`, LM qualification `55693462`; parent arrays `55693468` and `55693470` | Eight parents and 240 workers completed. All four repeat pairs have zero differences at every saved shell in iterations 1–3 and final reconstruction (128 shells beta-gal, 96 APOF); paired final-map SHA256 values match. Original project hashes unchanged. |
| Corrected frozen APOF pre-ML probe | `07afabbbf` plus FSC patch; qualification `55697434`, bootstrap `55697680`, corrected probe `55697779` | Six 96-shell tables: repeated seed 20261123 is byte-identical; raw/masked curves match across five seeds. Frozen-input audit records all five MRC files unchanged and matching bootstrap sources. Randomized/corrected maximum per-shell ranges are 0.04205675609/0.02873989944. Production bootstrap raw/masked agreement is within 5.97e-8 against existing 2e-6 tolerance. |
| Ten-run APOF solver comparison | `07afabbbf` plus FSC patch; array `55697950`, collector `55697951` | Supplementary downstream context across five seeds and two solvers. All ten parents/300 workers completed. This is not a same-solver duplicated-LM study and is not used as a gate for the separate LM bug. |

Source comparison: reconstructed the four numerical/interface files from both preserved patches (`hpc/keyed-source.patch` on `a4788cb44`, and `five_seed_review/probe03/transfer-production-07af.patch` on `07afabbbf`) in memory. Their executable Fortran token sequences match current `simple_fsc.f90`, `simple_rnd.f90`, `simple_image.f90` and `simple_image_filt.f90`, with quoted string contents preserved and formatting/comments excluded. Current `simple_halfmap_diagnostics.f90` and `simple_volassemble_strategy.f90` also match the later research base by this comparison. This establishes reuse of the FSC implementation, not identity of the entire historical refinement executable. Other source changes exist between the research bases and master; current-target synthetic and build/integration validation is supplied by Phases 0 and 3.

Local evidence was checked without HPC access: the saved comparison JSON has zero full-shell deltas for four pairs at four stages, and all four paired final-map hash records agree. Independently read all six 96-shell probe tables, confirmed byte-identical repeated-seed files and exact raw/masked columns, and reproduced the two reported per-shell ranges. Large maps remain remote; their recorded hashes were inspected locally, not freshly recomputed remotely. Earlier wrong-stage probe `55697681` and failed-launcher attempts remain excluded and preserved.

Evidence locations (local, gitignored; historical source identities retained):

- Same-seed campaign review: `.codex/diagnostics/polar-optimizer-realdata-ab/20261124-fsc-repeats-01/result.md`; underlying comparisons, map-hash records and logs: `codex_logs/fsc_phase_repro_20261124/realdata_review/` (`summary.json`, `final_map_hashes.json`, `scientific_parent_comparison.json`, `evidence/`).
- Corrected probe review: `.codex/diagnostics/polar-optimizer-realdata-ab/20261124-five-seed-apof-preml03/result.md`; underlying tables, probe source and frozen-input audit: `codex_logs/fsc_phase_repro_20261124/five_seed_review/probe03/evidence/` (`seeds/`, `fresh_input_audit.json`, `probe_source.f90`).
- Supplementary ten-run review: `codex_logs/fsc_phase_repro_20261124/apof_seed_campaign/review/report.md`; individual decisions and chronology: `.codex/branch-context/codex-fsc-phase-repro.md`.

Meaning: Phase 4's FSC-repeatability evidence is complete by explicit reuse, not by new execution on master `89b1b4817`. The five-seed spread remains descriptive. Final quality, APOF cFAR comparison, the unavailable historical bootstrap and the unrelated sporadic-LM defect remain outside this acceptance claim. No files in the historical evidence were changed, and no new experiment, HPC operation, commit or push occurred.

**Phase 4 handoff closed.** The then-outstanding maintainer review and final report were completed in Phase 5 (decision 11). No further HPC campaign is required by the revised Phase 4 scope.

**Tester policy follow-up, 2026-11-25 (complete; decision 9).** The user authorized fixing the three tester-policy gaps without changing production behavior. Replaced direct intrinsic reseeds with `set_fixed_seed` from `simple_test_utils` and direct key bases with `offset_seed`; each base receives the offset once. The diagnostic worker-stream reset still uses `random_seed` deliberately to set individual intrinsic worker states without racing shared library seed state, but its base now comes from `offset_seed` in serial code and its words use bounded 64-bit arithmetic. Removed `system_clock`, stored elapsed times and `FSC_TIMING` output from the fast tester. Main FSC and direct keyed-image comparisons request a fixed three-thread team and each asserts the actual team size when OpenMP is enabled. Thread count and dynamic settings are restored on exit. Existing numerical tolerances and fixture data remain unchanged. Two new team assertions make the OpenMP suite 62 checks; earlier 60-check results above describe the previous tester.

The production RNG impact is now explicit in contract 2: one intrinsic seed draw replaces per-coefficient worker draws, so old RNG stream positions and subsequent unkeyed draws may differ. This is documentation of the existing result-changing fix, not a further production change.

Renewed validation of the exact revised tester passed on the existing GCC15.2 caches; builds were incremental, with at most five compiler jobs and serial test processes. Installation followed each successful fast gate. OpenMP was capped at four; the tester's requested and asserted team is three; BLAS was bounded to one.

| Follow-up check | GCC15.2 Debug | GCC15.2 Release |
| --- | --- | --- |
| Incremental build, fast gate, install | PASS; 15/15 suites, 140.17 s | PASS; 15/15 suites, 51.78 s |
| Reverse `unit_reconstruction` | 2,016/2,016 assertions | 2,016/2,016 assertions |
| Three fresh FSC repeats, `seed_offset=0` | 62/62 checks each; actual three-thread team | 62/62 checks each; actual three-thread team |
| Two fresh FSC repeats, `seed_offset=1` | 62/62 checks each; actual three-thread team | 62/62 checks each; actual three-thread team |
| Offset arithmetic and coverage | Effective seeds 20261124 / 2674696885; offset changes randomized curves while raw/masked curves remain unchanged | Same |
| Same-seed repeatability | All printed deltas zero; 96 full-shell records agree across repetitions at each offset | Same |
| Negative control, `OMP_THREAD_LIMIT=1` | Expected exit 1: exactly the two team-size assertions fail (60/62 pass) | Same |

The Debug advisory overrun already existed in Phase 0 (141.51 s); no budget or numerical tolerance was relaxed. The one-thread-cap failures are intentional coverage controls, not failed qualification gates, and their logs/fixtures are retained. There is no `system_clock`, timing field, `FSC_TIMING` output or direct `seed_rnd_fixed` call in the fast tester. External validation tooling alone records elapsed time. An initial compiler-path regular-expression error in that external driver was corrected before any build/test ran; its original log/result are preserved as `driver-initial-preflight.log` and `initial-preflight-result.json`.

All 1,218 existing source/build files, including every production file, remain byte-identical to the pre-follow-up snapshot. Only the already-untracked FSC tester and this living plan were edited. The four numerical/interface source files and four installed executables of the unfixed Phase 0 baseline also still match their retained hashes; preserve that baseline, its builds and evidence until the user authorizes cleanup after approval and submission to origin. No changes were made to historical experiment records, no HPC operation was performed, and nothing was staged, committed or pushed.

Evidence: local, gitignored `codex_logs/fsc_test_policy_20261125/` (`results.json`, `qualification-contract.json`, `candidate-source-sha256.json`, `unchanged-source-before.json`, `tester-before.f90`, `tester-after.f90`, `tester-policy.patch`, `plan-before.md`, `plan-after.md`). Logs: `build-{debug,release}.log`, `fast-gate-{debug,release}.log`, `install-{debug,release}.log`, `reverse-reconstruction-{debug,release}.log`, and the offset/repeat/negative-control `fsc-*.log` files. Exact cached compiler/build flags are retained in `CMakeCache-{debug,release}.txt`. This supersedes the earlier tester's 60-check qualification for the current test source; production behavior and the explicit Phase 4 evidence-reuse decision remain unchanged. Phase 5 was outstanding at this historical checkpoint and is now closed by decision 11.

**Manual-review cleanup, 2026-11-25 (complete locally).** User authorized tester, documentation, comment and style fixes after manual review. Reset the canonical test seed with `fix_test_seed()` before leaving the FSC fixture; organize its assertions into named `test_<what>` routines. Move the seven direct image-phase assertions to `simple_image_tester.f90` in the existing image sub-suite, with a deterministic centred impulse that has Fourier amplitude at every shell. The FSC map/noise fixture remains the same: only kind/PI spellings use repository constants. The moved image test also resets the canonical seed and restores OpenMP team/dynamic settings. Name the independent FSC scope expectation (52 = 3×16+4), add the FSC RNG field-map comment, and update the policy suite inventory. No new suite, target or CLI control. Preserve all assertion guarantees and tolerances.

Pre-existing numerical behavior for Hans: both the former `ran3` path and the keyed path randomize conjugate coefficients independently on the h = 0 plane (also the stored Nyquist plane for even boxes); they do not enforce Hermitian symmetry there. This repeatability change does not alter that behavior.

Explicit exclusions for this cleanup: retain `real(TWOPI,dp)`, field-by-field key construction, the fixture noise hash and the documented serial-call contract without an `omp_in_parallel()` guard. Snapshot/delta and local validation logs: `codex_logs/fsc_review_cleanup_20261125/`. No commits or pushes.

Latest validation on the cleanup snapshot, macOS GCC/G++/GNU Fortran 15.2.0, incremental builds/install with five compile jobs and serial test processes (OpenMP cap four, BLAS one):

| Check | Debug | Release |
| --- | --- | --- |
| Fast gate | 15/15; 140.60 s, advisory budget warning | 15/15; 51.50 s |
| Reverse reconstruction | 2,009/2,009 | 2,009/2,009 |
| Reverse image area | 996/996 | 996/996 |
| Focused image suite at offsets 0 and 1 | 392/392 each | 392/392 each |
| Focused FSC: three fresh processes at offset 0, two at offset 1 | 55/55 each; actual three-thread teams | 55/55 each; actual three-thread teams |
| One-thread-cap negative controls | FSC 54/55, image 391/392; only the intended team assertions fail, exit 1 | Same |

All same-seed FSC deltas are zero; all 96 printed shell records agree across repeats at each offset and with the corresponding pre-cleanup record at offset 0. Offset 1 moves random-phase draws while leaving raw/masked curves invariant. The seven image assertions moved areas; they were not removed, explaining reconstruction's change from 2,016 to 2,009 assertions. Seed cleanup follows the existing `simple_eul_prob_tab_tester` pattern; reverse order checks pass in both affected areas. The Debug warning predates this cleanup (previous follow-up: 140.17 s; Phase 0: 141.51 s). No numerical tolerance or gate budget changed.

Evidence is local and gitignored: `codex_logs/fsc_review_cleanup_20261125/results.json`, `new-since-manual-review.patch`, `before/`, `before-sha256.json`, `candidate-source-sha256.json`, `production-token-verification.json`, compiler version files, cached configuration/flag files and all command logs. This record supersedes the previous tester qualification for the current test source. The four production files retain identical executable tokens, including OpenMP directives, after excluding ordinary comments/spacing; Phase 4's numerical-source evidence reuse remains applicable. The unfixed Phase 0 source, baseline builds and evidence remain retained. No HPC/Linux validation is claimed for this cleanup. Phase 5 and decisions 2–4 are now closed by the approval recorded in decision 11.

**Phase 5 closure, 2026-11-25.** The maintainer reported Hans’s approval and requested completion of this document. Final review checked source ownership, the optional image API, keyed scope/coordinate mapping, caller serial context, test registration and lifecycle/seed cleanup against the plan and repository policies. Registry, FLEX DAG, source headers, submodule linkage, guide imports, CLI chain and whitespace checks passed. Compared with the latest GCC15.2-qualified source manifest, the only source drift was two import-spacing changes in `simple_fsc.f90`; executable tokens are unchanged. No runtime qualification was repeated for this documentation-only closure. No production or tester code was changed, no evidence was deleted, and no commit or push was made.

## Final report

### What was learned

The defect is reproducible on the unfixed current base: across Debug and Release, six serial baseline cases pass while all twelve ordinary multithreaded same-seed cases fail in randomized/corrected FSC. Raw/masked curves and inputs remain invariant. One serial per-call seed plus coefficient keys removes dependence on OpenMP worker RNG history while retaining fresh realizations without reseeding. Same-seed replay requires the same call sequence; it does not preserve old RNG stream positions or old corrected curves.

The latest local GCC15.2 Debug/Release checks show zero same-seed deltas and repeated shell-record identity, including explicit worker-stream controls and actual three-thread teams. Historical real-data evidence supplies four identical repeat pairs and a repeated frozen-probe seed; reuse is justified by the production token comparison, not by claiming those jobs ran on current master.

### What changed

| Stage | Source/document scope and result |
| --- | --- |
| Phase 0 | Three measurement-only tester/registration/UI files in the isolated baseline worktree; four production/interface files unchanged. Preserve that source, both baseline builds and evidence. |
| Phase 1 | Eight Fortran files transferred: four production/interface files and four tester/registration/UI files. No optimizer or FSC correction-policy change. |
| Phases 2–3 | Build/install and qualification; no additional production implementation. |
| Phase 4 | Evidence reuse and source comparison; no new experiment or code. |
| Tester-policy follow-up | One tester plus this note: offset-aware seeds, no real timers and asserted fixed teams. |
| Manual-review cleanup | Eight already-scoped/newly-listed Fortran files, the test-policy inventory and this note: test ownership/organization, canonical seed cleanup, comments and formatting. Production executable tokens unchanged. |
| Phase 5 | Final review and closure of this note, moved by the maintainer to `completed/`; no code change. |

Final commit scope is nine Fortran files (four production/interface, five tester/registration/UI), one test-policy document and this note: eleven files total. The unrelated pre-existing `.gitignore` edit is excluded. The file table above is the authoritative list.

### Decisions and deviations

Hans’s approval, as reported by the maintainer, accepts decisions 2–4: the serial 53-bit seed, separate FSC scope/half/coordinate keys and the result-changing effect. Decisions 7–8 authorize historical real-data reuse and descriptive reporting of the five-seed spread; no retrospective statistical gate was invented. Decisions 9–10 authorize the tester corrections and explicitly retain `real(TWOPI,dp)`, field-by-field keys, the fixture noise hash and the absence of a parallel-region guard. Decision 11 closes review without authorizing Git actions.

The original eight-file scope expanded through approved tester-policy cleanup to add the existing image tester and policy inventory; their rows were added before editing. Direct image assertions moved to their owning area with an impulse fixture; the FSC noise fixture and all tolerances were preserved. No source file outside the final table was changed for this fix. The final report stays in this single note instead of creating the companion report proposed in the original template. Debug’s advisory gate overrun was present on the baseline and was not hidden or waived by changing the budget. Failed launchers, the wrong-stage research probe and intentional negative controls remain preserved and distinguished from passing qualification.

### Evidence against the baseline

Latest macOS GCC/G++/GNU Fortran 15.2.0 Debug and Release each pass all 15 fast-gate entries, 2,009 reverse reconstruction assertions, 996 reverse image assertions, two 392-assertion focused image runs at offsets 0/1 and five 55-assertion focused FSC runs. Requested and observed teams are three; all same-seed differences are zero. Across offsets, random-phase draws change while raw/masked curves remain invariant. One-thread-cap controls exit 1 solely at the intended team assertions. Gate times are 140.60 s Debug and 51.50 s Release, versus baseline 141.51/51.94 s; these are budget observations, not a performance claim.

Reused research comparisons establish identity of every saved FSC shell at iterations 1–3 and final reconstruction, plus matching recorded final-map hashes, in four same-solver pairs. The frozen APOF five-seed probe has a byte-identical repeated seed, invariant raw/masked curves and unchanged frozen inputs. Maximum per-shell randomized/corrected ranges are 0.04205675609/0.02873989944. This is N = 5 descriptive variation, not a general uncertainty bound.

These are **local evidence** records under gitignored `codex_logs/` and `.codex/`, with original source/build identities retained in Progress. Latest qualification: `codex_logs/fsc_review_cleanup_20261125/results.json` and its source manifests/logs. Baseline: `codex_logs/fsc_phase0_local_20261125/`. Historical real-data comparisons and probe tables: the Phase 4 evidence locations above. This repository note preserves the review summary; the raw artifacts are not included in the commit. Large remote-map hashes were inspected from retained records, not freshly recomputed remotely. No current-master Linux/HPC run is claimed.

### Remaining scope and authorization

No implementation decision or acceptance gate remains open for this FSC fix. Reported, not gated: the historical APOF LM cFAR comparison (0.9675 versus 0.9596), descriptive seed-dependent cutoff variation and the pre-existing lack of enforced Hermitian symmetry on stored boundary planes. These are not claims of regression-free scientific quality. Data-dependent limits: the historical bootstrap is unavailable and final-map quality is not established by the replacement frozen probe. Separate work: sporadic LM failures and duplicated-LM campaigns; no follow-up run is required to close this plan.

Retain the unfixed Phase 0 source, builds and all validation/experiment evidence until approval and submission to origin have occurred and the maintainer separately authorizes cleanup. No staging, commit, push or cleanup is authorized by completion of this plan.
