# Adopting X's Execution Layer in SIMPLE

**Status:** Outline for developer review, 2026-09-02

**Depends on:** nothing. Independent of, and prior to, the staged DDD
refactor (`staged_domain_driven_design_refactor.md`).

**Delivery:** small commits on `master`, each behavior-preserving unless it
adds an opt-in key (default off). The user compiles and runs the gates.

## 1. What X is and what we take from it

X (`/Users/elmlundho/src/X`) is our cross-linking MS engine. It borrowed
SIMPLE's commander, exec-router, qsys, and OpenMP patterns, then modernized
the command line and the distribution model. Those modernizations come back
here. Read in X: `docs/execution_layer.md`, `docs/architecture_principles.md`,
`AGENTS.md` (OpenMP rules), and the sources `xlms_cmd`, `xlms_cmd_spec`,
`xlms_exec_strategy`, `xlms_qsys`, `xlms_qsys_ctrl`, `xlms_commander_search`.

Ground rule: every mechanism lands inside the SIMPLE layer that already owns
the concern (`cmdline`, `ui_program`, the executables, `qsys`, the
parallelization strategies). No parallel registry, configuration path, or
executable.

## 2. Command object and CLI

| X | SIMPLE today | Decision |
|---|---|---|
| Typed command entries (int/real/str/flag tag), growable, strict typed getters | `cmdline`: 100 fixed slots, `real(dp)` + `string` per entry, type sniffed, `get_iarg = nint(rarg)` | Adopt, in place |
| `to_string`: exact round trip; the command object is the unit of invocation and of distribution | `gen_job_descr` to `chash` to script text; a second serialization | Adopt; one serialization, `prg=key=val` grammar kept |
| Parse returns status + message; the driver maps it to exit 2, no backtrace | `cmdline%parse` has eleven bare `stop`s and `THROW_HARD`s; not usable as a library | Adopt |
| Overlength argument refused with a named error | Silent truncation at `XLONGSTRLEN` | Adopt |
| One spec per subcommand = whitelist + validator + defaulter + help; options declare unit, inclusive/exclusive range, closed value set; NaN/Infinity rejected; help prints exactly what parse enforces | `ui_program`/`ui_param` (about 14k lines) hold type, choices, units, required-ness, display defaults; parse whitelist is a separate generated `simple_args.f90`; no ranges anywhere; non-finite reals pass; range checks scattered | Adapt: put ranges/strictness/finiteness on `ui_param`, make the registered program the parse authority, render the domain in help and JSON, retire `simple_args.f90` |
| Defaults declared once (typed constant + string form) and written by the spec into the command object; presence-semantic keys deliberately undefaulted | Four default sources; `default_audit.py` generates 3k lines to recover display defaults; CLI/GUI divergence accepted by policy | Adapt: per program, literal `if(.not. defined) set` overrides become declared defaults; `derived` marker for data-dependent ones; explicit presence-semantic list |
| No results read back out of the command object | `cline%get_iarg('endit')` after refine3D in five places | Adopt for new code; migrate old side channels when the workflow is touched |
| `--key value` grammar, one binary, specs declared inside the driver | `prg=key=val`, four executables, registry in `ui` | Reject all three (external contracts; second registry) |

## 3. Distributed execution

| X | SIMPLE today | Decision |
|---|---|---|
| One `exec_role(cline)`: worker if `part`, master if a concurrency bound, else shared-memory; same commander for all roles | The role decision is written locally in each of the fifteen strategy factories (three variants that disagree on edge cases), in `set_shmem_flag`, and in the `parameters` worker flag | Adopt one resolver used by all of them (section 3.1); the strategy modules and their scheduling, merge, and shared-memory code stay; keep `simple_private_exec` as the worker entry |
| Exit status captured before any pipe; exit-code file always; success sentinel only on success; controller sees running/done/failed; bounded retries; watchdog; run dir kept on failure | Batch scripts pipe through `tee`; exit-code files only in streaming; sentinel touched by the worker Fortran itself, so a crash polls forever; no retries | Adopt; retries and watchdog as opt-in keys, default off |
| Slots + part-size floor; part count derived; queue topped up (self-balancing) | `nparts`, `ncunits` user knobs; refill loop already exists; part count is persistence-visible (algndocs, partial recs, chunks, restart) | Adapt as opt-in only; derived count recorded for restart; pilot on preprocessing |
| Chunk commit: tmp file, atomic rename, manifest (hash + row count); merge refuses uncommitted or mismatched parts | Parts written directly under final names; presence is trusted | Adopt per workflow; hash mandatory for metadata docs, opt-in for large binaries |
| Merged distributed output byte-identical to serial, gated on two split geometries; every parallel path has a bit-identity test | Equivalence required in principle; many reductions are order-dependent | Adopt with a declared class per path: `bit_identical` or `governed_tolerance`; new kernels reduce lane-locally and merge in fixed order |
| Typed `qsys_settings` value object; unset field omits the directive; CLI-pure, no environment fallbacks | `chash` queue descriptor; `SIMPLE_QSYS*` env fallbacks; `compenv` persisted in the project | Adopt the typed object with precedence cline > compenv > env; reject CLI purity (compenv is a NICE contract) |
| Golden-script test per backend | None | Adopt first, before touching script generation |
| Rejected in X: persistent worker pool, coarray, SLURM arrays, memory estimator | Serve streaming and scale here | Keep; must not regress |

### 3.1 The role resolver, and what stays

The fifteen modules under `strategies/parallelization` (about 12.5k lines)
are not the target. Their `distr`/`master` types own script generation
through `qsys_env`, scheduling, and merging; their `inmem`/`shmem`/`worker`
types own the shared-memory loops. That is where X says parallelization
code belongs. Only the role decision at the top of each factory is replaced.

Today that decision has three variants:

| Variant | Where |
|---|---|
| `nparts` defined and `part` not defined | cluster2D, refine3D, calc_pspec, rec3D (gridding) |
| the same, and `fromp`/`top` not defined | ctf_estimate, motion_correct, pick, extract, reextract, preprocess, gen_pspecs_and_thumbs |
| `nparts > 1` and not worker; worker is `part` defined (denoise_project, cls_split) or `part`, or `fromp`+`top` with `nparts > 1` (make_cavgs); rec3D PCG adds `nparts > 1` | denoise_project, cls_split, make_cavgs, rec3D (PCG) |

Plus `set_shmem_flag` (`nparts` absent or 1 means shared memory, and it
deletes `nparts`) and `parameters%l_distr_worker` (`part` defined),
published as a global read by the error handler, jiffys, the sentinel
writer, the polar memo, and both matchers.

The variants disagree: with `nparts=1`, cluster2D and refine3D build a
distributed master with one part while make_cavgs and abinitio2D go shared
memory; with `fromp`/`top` set and no `part`, the first group says master
and the second says worker. E0 tabulates which shapes occur in production
and records the current answer for each.

E4 then adds one function in the parallelization layer:

```fortran
role = exec_role(cline)   ! ROLE_WORKER  if part defined
                          ! ROLE_MASTER  if nparts defined, > 1, not worker
                          ! ROLE_SHMEM   otherwise
```

and each factory becomes a `select case` on it, mapping roles onto the
concrete types it already has (two-role modules map worker and shared
memory onto `inmem`, as now; rec3D keeps its backend branch). No strategy
body, deferred interface, or hook changes. `set_shmem_flag` becomes
`exec_role(cline) == ROLE_SHMEM` with the `nparts` deletion kept as an
explicit normalization; `parameters` derives `l_distr_worker` from the same
function. The global flag stays; removing globals is a separate program.

E4's failure semantics land in `qsys_ctrl`/`qsys_env` beneath the strategies,
so a `distr` strategy sees a controller that reports a failed part without
per-strategy edits. E6's chunk commit touches the merge helpers the `distr`
strategies call, one workflow at a time.

Not proposed: collapsing the fifteen strategies into one role-dispatching
commander (X could, having one workflow; SIMPLE's workflows differ in
partition unit, merge rule, and restart state), or unifying the two-role and
three-role module shapes.

## 4. OpenMP hot-path rules (new and migrated code)

- Nothing allocates, constructs objects, dispatches on `class(...)`, does
  formatted I/O (internal writes included), or passes internal procedures
  inside a parallel region. Per-thread workspaces sized before the region,
  selected by `omp_get_thread_num()+1`.
- Immutable shared state separate from per-thread state; `intent(in)` self on
  the hot call.
- Hot kernels are pure procedures on plain contiguous arrays behind a thin
  validating wrapper.
- Engine-lifetime invariants travel as one validated context object, not a
  long argument list and never a back-reference from a workspace.
- Derived types with default initializers are resident on allocation; size
  per-worker buffers by what is retained, not the population.
- Serial optimization first; counters first, then the knife: a structured
  per-thread counter type, summed, drives every optimization decision.
- Every parallel path ships its identity or equivalence test.

## 5. Tests and build

- ctest labels (`unit`, `regression`, `platform`), a fast inner loop with a
  time budget, `WILL_FAIL` for abort-path tests, data-dependent tests
  registered only when a cache variable names the data.
- Benchmark writes a dated, git-stamped report only after its determinism
  gates pass; peak RSS per phase in the report.
- Zero-warning tree in Debug and Release, warnings-as-errors option for CI,
  FP traps in Debug, `fix_warnings.py` (X's port of our perl script).
- Finalizers on resource-owning types only; keyword-only optionals at call
  sites; no `==` on reals.

## 6. Steps

Each step has a characterization gate captured on unchanged `master` first.

| Step | Content | Gate |
|---|---|---|
| E0 | Characterize: `simple_test_cmdline` accept/reject matrix; golden scripts for local/SLURM/PBS/SGE/LSF/persistent-worker/coarray; job-description round trip; role table (command-line shape x factory, including `nparts=1` and `fromp`/`top` without `part`, section 3.1); inventory of missing-key overrides and command-object read-backs | Passes on unchanged master; each role disagreement has a recorded decision |
| E1 | Typed `cmdline` (tag, growable, strict getters); `to_string` as the only job description; status-returning parse; exit codes in the executables; overlength refusal | E0 matrix unchanged except documented exit-code changes; golden scripts byte-identical |
| E2 | Ranges/strictness/finiteness on `ui_param`; registered program is the parse authority; help and JSON render the domain; retire `simple_args.f90` after a superset check | Enforced domain equals advertised domain per program; JSON gains fields only; no default changes |
| E3 | Per program: literal missing-key overrides become declared defaults; `derived` marker; presence-semantic key list; `default_audit.py` reports nothing for migrated programs | CLI behavior unchanged; JSON defaults unchanged except where the audit previously failed |
| E4 | `exec_role` replacing the factory heuristics, `set_shmem_flag`, and the `parameters` worker test (factories only; strategy bodies untouched); exit-code capture, sentinel on success, tri-state parts; opt-in retries and watchdog; run-dir retention; typed `qsys_settings` | Every factory selects the same concrete type as before for every shape in the E0 role table; golden-script diffs reviewed; local failure-path `WILL_FAIL` test; overlays pass existing gates |
| E5 | Opt-in part-size floor deriving `nparts`; derived count recorded for restart; `ncunits=0` fits the local machine | Preprocessing pilot equals fixed-`nparts` runs; restart re-derives the same split; default path unchanged |
| E6 | Per-part commit manifest and refusing merge; manifest-first resume for workflow commanders | Two split geometries merge byte-identical; truncated part refused; mismatched manifest refuses resume |

E0 to E3 touch `cmdline`, `ui`, the executables, and `params` parse only.
E4 to E6 touch `qsys` and `strategies/parallelization` and need the platform
gates (local end to end; schedulers at script level until a cluster run).

## 7. Questions for developers

1. Exact real formatting in `to_string`: do we want bit-exact values in job
   descriptions (a behavior change for workers) or the current rounding?
2. Which programs go first in E3? Proposal: preprocessing programs, whose
   overrides are few and literal.
3. Retry and watchdog defaults: off (behavior-preserving) or on with
   conservative values?
4. Should the derived part count (E5) ever become the default for
   preprocessing once piloted?
5. Hash policy for large binary partials in E6: size-only, sampled, or full.
6. Are there consumers of `simple_args.f90` or the `chash` job description
   outside the tree (NICE, scripts) that E1/E2 would break?

## 8. Not in scope

The `prg=key=val` grammar; merging executables; removing the persistent-worker,
coarray, array, or subproject paths; a derived part count by default; any
change to the `.simple` format or to `compenv` semantics; retrofitting
existing kernels to section 4.
