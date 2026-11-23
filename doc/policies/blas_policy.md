# BLAS Policy

Scope: every BLAS or LAPACK call in SIMPLE-owned code, the bundled OpenBLAS that
provides both, and the OpenMP context of those calls. Written against
`0a8b0903e`; library facts from the OpenBLAS 0.3.34 source and probes on the M2
Ultra (2026-11-23). Section 7 lists where the code does not follow it yet.

## 1. The library

`cmake/Dependencies.cmake` builds OpenBLAS 0.3.34 from a SHA-256-pinned archive:
static, LP64, with LAPACK, `NUM_THREADS=128`, `NUM_PARALLEL=1`. It is cached in
`.deps/numerics/<key>/` (system, CPU family, compilers, flavour) and rebuilt when
its recipe stamp changes. The flavour is OpenMP (on SIMPLE's libgomp pool) with
GNU C and Fortran: always on Linux, on macOS when AutoDetectGCC picks Homebrew
GCC; with AppleClang as C compiler (also a macOS configure that sets only the
Fortran compiler) it is sequential. CMake's dependency summary prints it.

Rules: no system BLAS, LAPACK, Accelerate, vendor library or second BLAS; no
ILP64 and no `integer(kind=8)` BLAS argument; a recipe change follows section 6.

## 2. What a call does

| Called from | OpenMP flavour | sequential flavour |
| --- | --- | --- |
| an active team (two or more threads, any nesting level) | calling thread | calling thread |
| serial code, a one-thread region, a region whose `if` is false | threaded by size, up to min(current `omp_get_max_threads()`, its value at process start, 128) | calling thread |

- **Start-time cap.** OpenBLAS reads `omp_get_max_threads()` at program load. A
  process started with `OMP_NUM_THREADS=1` never threads a product, whatever
  `nthr` is later set to (probed: 1.00 cpu/wall; 7.7 with the variable unset).
  After a `fork()` the next threaded product re-reads it; stream stages are
  forked children of the master, so they inherit its launch environment.
- **One-thread bounds** (from serial code): real `gemm` m·n·k ≤ 262 144 generic,
  `sgemm` < 125 000 on Neoverse V2, both covered by `SGEMM_ONE_THREAD_MNK =
  124 999` (`simple_linalg`); complex `gemm` m·n·k ≤ 32 768. On arm64 a real
  `gemm` with one row or column goes to `gemv` (bound ≈ 25 000 on Neoverse,
  460 800 generic).
- **Concurrency.** Each thread inside a call holds one of 256 buffer slots; past
  them OpenBLAS warns "precompiled NUM_THREADS exceeded" and may crash. Teams up
  to 128 threads are covered.
- An element's last bits depend on the kernel (the build host's CPU), the call's
  shape, its column (`cgemm`: not its row within a full call, probed on x86) and
  whether the call was threaded. Identity is per build and machine only.

Never call `openblas_set_num_threads(_local)` (in-team calls would then start
nested teams) or rely on `OPENBLAS_*` variables; never call BLAS from a non-OpenMP
thread; size teams that call BLAS by the team that runs, at most 128.

## 3. Rules at a call site

**Placement.** Parallelism belongs to SIMPLE's loops, not to the BLAS.

- *In-team, always.* A team loops over independent items; each item's products
  run on its thread with that thread's work arrays (`bwork(ithr)`, the shift-map
  work), sized by the running team. A wide product is cut into chunks of a fixed
  size and the chunks go to the team (the picker's batches of 256 windows, flex's
  direction assignment in chunks of 128 particles, HHCNN's templates one class a
  pass): its speed then does not depend on the launch environment, nor its
  result on the team. Never leave a wide product to OpenBLAS's own threads, and
  never serialize in-team calls (a critical region): they are independent.
- The one exception is `flex_traj_denoise`, which measures once per run whether
  the library threads (`probe_blas`) and, if so, calls its products from serial
  code.
- A call in a one-thread or `if(.false.)` region is serial code; `matmul` stays the
  intrinsic (no `-fexternal-blas`); a BLAS replacement of loops is adopted on an
  interleaved timing at production sizes and 24 threads.

**Identity.** A product whose result must not depend on the team, batch or
caller has a shape fixed by the problem's geometry, never by the team or how
full a batch is: pad partial batches (`score_block_posterior` always multiplies
`PFTC_BLOCK_NB` × `PFTC_BLOCK_NR` = 16 × 32, absent rows and columns zero), keep
an item in a fixed column, and, if a serial caller can reach it, cut it into
blocks below the one-thread bound (`shell_sums`) or state that it must run in a
team. Size products in `int64`. The BLAS arithmetic is not frozen against the
loops it replaced: a change is judged by tolerance, team identity and the
workflows.

**Interfaces and arguments.** `simple_linalg` owns every explicit BLAS/LAPACK
interface (reference argument order, `integer(kind=4)`, explicit kinds,
assumed-size dummies, intents); callers `use simple_linalg, only:`. No interface
in a consumer module and no implicit call. LAPACK goes through the wrappers
(`eigh`, `jacobi`, `matinv`, `svdcmp`, `svdfit`, `qr_solve`); extend one before
adding a direct call. Pass sub-blocks by first element and the parent's leading
dimension or a contiguous pointer view, never a non-contiguous section. In a team,
errors set a status raised after the region.

**Comments and constants.** The comment at a call states, in itself, its context
(in-team or measured), what fixes its shape and what it keeps identical
(`shell_sums` is the model). A bound or shape relied on is a named parameter with
its reason (`SGEMM_ONE_THREAD_MNK`, `PFTC_BLOCK_NB/NR`).

## 4. Why: performance

In-team calls of the hot shapes on the M2 Ultra (16 P + 8 E cores; each thread
calls on its own arrays; medians of three interleaved passes, 2026-11-23),
aggregate GFLOP/s:

| Shape | Library | 1 thread | 16 threads | 24 threads |
| --- | --- | --- | --- | --- |
| block kernel, `cgemm` 16 × 32 × 80 | OpenBLAS | 73 | 1 026 (0.86 per thread) | 1 029 |
| | Accelerate | 60 | 131 | 122 |
| shift map, `sgemm` 512 × 8 × 24 | OpenBLAS | 46 | 557 (0.74) | 573 |
| | Accelerate | 422 | 1 853 | 2 033 |

In-team OpenBLAS scales; threading inside BLAS at these sizes (0.1–0.3 Mflop a
call) would only add fork-join cost. Accelerate runs real products on the
matrix coprocessor, one per core cluster, so its team throughput saturates at
the unit count; in a team its complex product is about 8× slower than
OpenBLAS's, its real product about 3.4× faster. On M4/M5 OpenBLAS 0.3.34 uses
NEON kernels for Fortran calls (its SME code serves only row-major CBLAS calls).
A critical region around HHCNN's products made its training step at 24 threads
as slow as at one (4.07 s against 0.27 s without, nk 25, width 32).

## 5. Testing

The fast gate runs at `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=4`: team-against-one-
thread checks work, but by the start-time cap no serial product is ever threaded,
so serial-against-team checks (`test_blocked_shell_sums`) only mean something by
hand: `env -u OMP_NUM_THREADS build/bin/simple_test_exec test=unit_pftc_align2D3D
suite=polar_shift_map_kernel`. A new or changed call site has: an independent
oracle check (loops or `matmul`) including partial blocks and strided operands;
for in-team use a forced team against one thread, bit for bit, varying batch fill
and item place where it claims independence of them; a full-machine stress run
by hand. Benchmarks: `polar_matching_bench` (`suite=fills`, `shiftmap`,
`block_prod`), `run_hhcnn_bench`.

## 6. Changing the recipe

A new version, flag, flavour rule or `NUM_THREADS` changes SIMPLE's arithmetic and
needs the owner's agreement, then on x86 Linux and arm64 macOS: the in-team,
start-time-cap and concurrency probes; the bounds re-read (`interface/gemm.c`,
`gemv.c`); the fills' row invariance and the shift maps' team identity; the fast
gate, `simulated_workflow suite=6vxx`, identity files re-baselined; the CMake and
`simple_linalg` comments, `installation.md` and this policy updated together.

## 7. Open items

1. HHCNN's thread-identity suite was removed as flaky on the remote CI
   (`0f01b5c66`); without the critical region of `2712df0bb` it passes 10 of 10
   at the full machine on the M2 Ultra, and a production-size training step is
   identical at 1 and 24 threads. It is not registered again.
2. `simple_hhcnn`, `simple_cavg_cnn`, `simple_opt_lbfgsb` declare their own
   interfaces; `flex_traj_denoise` and `flex_probe_fit_estep` call `sgemm`,
   `sgemv`, `dgemv`, `ssyrk` with none.
3. `SIMPLE_BLAS_TEAM_SAFE` is always defined; the non-team-safe branches are dead.
4. Under the start-time cap `probe_blas` takes the serial route with fewer
   threads than the team when a job starts with 2 to `nthr` - 1 threads;
   `test_environment_policy.md` explains the gate's BLAS limit wrongly.
5. `score_block_posterior` passes the complex bound when nk > 32 and relies,
   unstated, on team or one-thread callers; `shell_sums`' one-magnitude path can
   pass `sgemv`'s Neoverse bound from a serial caller. The `cgemm` row invariance
   is not probed on arm64.
6. The cache key has no CPU target (a `.deps` shared by different x86 generations
   reuses the first node's kernels); `lapack_stop` exits with status 0.
