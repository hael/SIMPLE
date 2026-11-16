# Porting SIMPLE's two hot spots to the Apple Silicon GPU

Proposal, 2026-11-16. Not implemented. Written for review by engineers on Apple's MLX team and by
SIMPLE developers. It belongs to the next release: release 4 is CPU-only by ruling
(`speed_architecture_proposals.md`, "batched, GPU-ready kernels").

[One-page briefing for Samy Bengio](apple_silicon_gpu_brief.md).

## Summary

SIMPLE is an open-source cryo-EM image processing program written in modern Fortran, with OpenMP
threading. It builds and runs on macOS. Two computations dominate 3D refinement: the probabilistic
search, which scores every particle image against every reference projection, and the reconstruction
of the 3D map from the aligned images.

The plan moves the dense part of the search to the GPU through Apple's C API for MLX (mlx-c), and keeps
the irregular, double-precision parts on the CPU cores. Because CPU and GPU share memory on Apple
Silicon, the GPU can score the next batch while the CPU threads refine the current one, using separate
particle buffers and shared read-only references. This avoids transfers between separate device
memories when buffers are allocated or wrapped appropriately. Layout conversion and a fallback from
external-buffer wrapping can still copy data within shared memory; both must be measured.

The search maps onto MLX well. Its core can be written as one batched complex matrix product, followed
by a batched inverse real FFT and one small custom kernel. Reconstruction maps less cleanly. Its central
step scatters values into a 3D grid, which on a GPU means atomic additions, and those collide with
SIMPLE's rule that results must not depend on thread count. Our own timings also put most of the
reconstruction cost outside the scatter step. Reconstruction therefore starts with a profile on the
target machine before any kernel is chosen.

## Terms

- **Particle**: one noisy 2D projection image of the molecule. A data set has 10^5 to 10^6 of them.
- **Reference**: a 2D projection of the current 3D map in one viewing direction. A typical run has
  2,500 directions per state and one to four states (conformations), so up to 10,000 references.
- **Polar Fourier transform**: SIMPLE stores the central Fourier region of each image on polar
  coordinates: `nk` rings (shells) of `pftsz` angular samples each. A typical size is `nk` = 61 and
  `pftsz` = 182. Only half of each ring is stored; the other half is its complex conjugate.
- **In-plane rotation**: rotating an image about the viewing axis. On the polar grid this is a cyclic
  shift along the angular axis, so all `nrots = 2*pftsz` rotations can be scored at once with an FFT
  along that axis.
- **CTF**: the contrast transfer function of the microscope, a real, per-particle filter.
- **sigma2**: the noise power per shell. The score weights each shell by `k/sigma2(k)`.
- **Even/odd halves**: particles are split into two halves, each matched against references from its
  own half-map.
- **Batch**: SIMPLE currently processes particles in batches of 50 per CPU thread, so about 1,200
  particles at 24 threads. GPU tile dimensions will be fixed independently of the CPU thread count,
  with masked padding for short tiles, to keep the GPU arithmetic stable when CPU scheduling changes.

## Where the time goes today

Measured on a ribosome benchmark (78 refinement iterations, 6,572 s; distributed over 5 processes of
16 threads, 3 states of 1,000 directions), recorded in
`doc/implementation_notes/rejected/particle_caching_3D.md`:

| part | share of the run |
| --- | ---: |
| probabilistic pre-step: the search table, plus its process launches and file reads and writes | 46% (73% of a late iteration) |
| assembly of the 3D maps (reading and summing partial volumes, statistics, filtering) | 15% |
| reprojection model: volumes read, masked, filtered, projected to references, written | 11% |
| particle preparation | 1.2% |
| partial reconstruction: particles read, prepared, inserted into the grid, written | 1.7% |

How the pre-step splits between the table itself and its overhead was not measured in that run.

On a beta-galactosidase run (506 s; 4 processes of 6 threads), the final full-size reconstruction took
48.5 s and assembly inside the iterations 26 s. The two search table steps took 82 s (the global
search table, `prob`) and 160 s (the local table, `prob_neigh`, of which 117 s was table time and the
rest process launch and file reads) (`prob_search_speedup.md`). In small test runs, inserting
7,400 to 12,700 poses into the grid took 1 to 2 s in total with the default gridding route and 4 to
10 s with the iterative route (`inpl_quadrature_insertion_3D_report.md`).

So the search is the clear first target. For reconstruction, the time sits in assembly, the full-size
final map and the reference projection, more than in the insertion loop itself. Both measured runs
were distributed; part of their cost is file traffic between processes, which a shared-memory run on
one Mac does not have.

## Hot spot 1: the probabilistic search

### What it computes

For a particle `X` and a reference `R`, the score at in-plane rotation `r` is a weighted Euclidean
distance between the particle and the CTF-filtered, shifted, rotated reference:

```
L_r = sum_k w_k |X_k - CTF_k * rot_r(S_s * R_k)|^2 / sum_k w_k |X_k|^2
```

SIMPLE evaluates it for all rotations at once in the Fourier domain of the angular axis. The particle
side and the reference side are each transformed once and kept ("memoized"):

```
A(p,k,i) = FFT_c2c(X * CTF)              B(p,k,j) = FFT_c2c(R)
C(p,k,i) = FFT_r2c(CTF^2)                D(p,k,j) = conj(FFT_r2c(|R|^2))

acc(p; i,j) = sum_k w_k(i) [ C(p,k,i) * D(p,k,j) - 2 * A(p,k,i) * conj(B(p,k,j)) ],  p = 1..pftsz+1
dist(r; i,j) = 1 + c2r_p(acc)(r) / (2 * nrots * wsqsum(i))
```

`p` is the angular frequency bin, `i` the particle, `j` the reference. The transforms run along the
angular axis over the full ring of `nrots` samples, rebuilt from the stored half, and keep bins
`1..pftsz+1`. `A` and `B` are complex-to-complex transforms; `C` and `D` are real-to-complex. The
inverse `c2r` is unnormalized, and distances above a fixed cutoff are set to a huge value. All memoized
arrays are single precision complex. The weights `w_k = k/sigma2(k)` depend on the particle's half; the
reference half follows the particle's half.

The search does not take the best rotation. For each (particle, reference) pair it draws one rotation
from the tempered distribution `exp(-T * (dist_r - min_r dist))`. The random numbers come from a
counter-based hash generator. Its key holds the run seed, the program and pass, the iteration, the
particle, the reference and a draw ordinal. The draw is therefore independent of thread count and
evaluation order. The hash uses only 32-bit integer operations and ports bit for bit to the GPU. The
uniform variate it yields (53 bits as a double) and the cumulative sum of the draw are double
precision on the CPU; the precision section below covers what that means on the GPU.

Each pair's result is two numbers: the distance and the drawn rotation. They fill one row per particle
over all references. A cheap per-particle step then draws one reference from a tempered softmax over
the row.

### The full table step

1. **Dense score.** Every particle of the batch against every reference, with the draw. In the default
   configuration the score is at zero shift. With the shift map switched off, single-state runs score
   at a per-particle seed shift instead; the shift is applied to the particle side, so only `A`
   changes. Cost: about 0.3 million floating-point operations per pair, so `nptcls * nrefs` pairs
   dominate. This is the part this plan moves to the GPU.
2. **Shift screen.** For a pool of the best references per particle, a shift map built from
   Fourier-Bessel sums and one batched 2D inverse FFT. A later GPU candidate.
3. **Shift refinement.** For the best few references, a bounded quasi-Newton search (L-BFGS-B) over
   the 2D shift, in double precision, with branching and early exits. Stays on the CPU.

### The dense score as a matrix product

For a fixed frequency bin `p` and a fixed half, the contraction over shells is a complex matrix
product. Concatenate the two terms along the shell axis:

```
P_p (particles x 2nk) = [ w_k * C(p,k,i) , -2 * w_k * A(p,k,i) ]
Q_p (2nk x refs)      = [ D(p,k,j) ; conj(B(p,k,j)) ]
acc_p = P_p @ Q_p       (particles x refs), batched over the pftsz+1 bins
```

With `nk` = 61 the inner dimension is 122. The weights depend on the particle's half and each half
reads its own references, so there is one product per half per bin, each over about half the batch.

A CPU version of this layout, using BLAS `cgemm` per bin, ran 3.3 to 3.8 times faster than the then
production per-pair routine on an M3 Pro, and about 2.4 to 3.1 times faster than its separately compiled
replica ([October layout benchmark](../../completed/polar_matching_layout_bench_report.md)). These are
historical kernel timings, not a speedup prediction for current production or a whole run. The benchmark
found arithmetic to be the limiting factor on that machine, which motivates testing the GPU.

Before serving as the GPU reference, the benchmark needs two corrections. Its replica and matrix
layouts omit the last angular bin (`p = pftsz+1`), which production includes. Their `pair_finish` also
uses `ran3()` while the production layout uses the keyed generator. Every layout must include all bins
and use identical variates keyed by the production scope, particle and reference. Then remeasure the
CPU baselines against current production.

On MLX the step becomes:

1. `mlx_matmul` on complex64 arrays of shape `(bins, particles, 2nk)` and `(bins, 2nk, refs)`, giving
   `(bins, particles, refs)`, once per half.
2. Transpose to `(particles, refs, bins)` and `mlx_fft_irfft` along the last axis, length `nrots`,
   with the normalization mode set so the transform is unscaled (`MLX_FFT_NORM_FORWARD`), as the CPU
   route expects. The view alone leaves bins strided; MLX's GPU FFT makes a contiguous input copy.
   Include that copy in the timing and peak-memory budget.
3. A custom kernel (`mlx_fast_metal_kernel_new` / `mlx_fast_metal_kernel_apply`) per (particle,
   reference) row: convert to distances with the cutoff, find the minimum, form the tempered weights,
   run the keyed draw, write distance and rotation.

Steps 2 and 3 are candidates for fusion into one custom kernel. Each row is a 364-point inverse real
transform plus a reduction; register use, threadgroup memory and the strided reads must be profiled.
Fusing avoids materializing the real rotation array, about 17.5 GB for a 1,200-particle batch against
10,000 references. The complex contraction is itself 17.6 GB, or 8.8 GB per half with 600 particles.
Both the fused and unfused routes therefore tile over references and particles. Fusion is adopted only
if it improves the complete tile pipeline over the library FFT plus a separate draw kernel.

### Layout

The memoized arrays are stored with the frequency bin fastest: `(bins, shells, reference, half)`. The
matrix product needs the reference operand as `(half, bins, 2nk, refs)`, with `B` conjugated and
stacked under `D`. This is a reordered copy, built once per iteration when the references change; the
October benchmark took under a second to build its version on the M3 Pro. The original arrays stay,
because the shift refinement reads them.
Copy-free wrapping therefore applies to the reordered operand, which SIMPLE allocates for this purpose,
not to the existing arrays. The particle operand is reordered the same way per batch.

### Sizes

At `pftsz` = 182 and `nk` = 61, one reference costs 89 KB per memoized array (`B` and `D`) per half.
10,000 references with both halves take 3.6 GB, and the reordered copy as much again. The particle side
(`A` and `C`) of a 1,200-particle batch takes about 214 MB, before packing and double buffering. The
output per batch is 1,200 x 10,000 pairs x 8 bytes, about 96 MB. Sizes here use decimal GB and MB.

For a tile of 600 particles from one half and 512 references, the complex contraction occupies about
450 MB, the FFT's contiguous input copy another 450 MB, and the real rotation output about 447 MB.
Budget for their overlapping lifetimes, kernel work arrays, allocator caching, both particle buffers
and SIMPLE's existing volumes and memos. Measure peak memory for the complete pipeline; being below
1 GB per intermediate does not bound their combined footprint.

### What stays on the CPU

The shift refinement and the per-particle assignment stay on the CPU. The custom kernel's outputs are
arrays that MLX allocates. After evaluation the CPU reads them in place through `mlx_array_data_*`,
with no copy, because memory is shared.

The batch loop becomes a two-stage pipeline: the GPU scores batch `n+1` while the OpenMP threads refine
batch `n`. This needs two changes on the CPU side. Today the table step handles each particle whole,
dense score, then screen, then refinement, inside one threaded loop; it splits into a batch-wide dense
phase and a per-particle phase. And the particle memo holds one batch at a time; it becomes two
buffers, so that one batch can be refined while the next is scored. Submit the next GPU graph with
`mlx_async_eval` before starting CPU refinement; wait for its outputs before consuming or reusing that
buffer. A lazy graph that has only been constructed has not launched GPU work. Read-only references
can be shared throughout; particle-buffer reuse waits for both CPU and GPU readers.

### Precision

The contraction is already single precision on the CPU. Double precision is used for the normalizer,
conversion to distances, uniform variate, exponential argument and evaluation, total weight, draw
target and cumulative sum. `draw_posterior_row` rounds the exponential weights to single precision
before accumulating them in double precision. The `arg < 80` cutoff, invalid-row fallback and CDF
comparison convention are also part of the contract.

Apple GPUs do not support native double precision. Candidate implementations use compensated
single-precision summation or two-float ("double-single") arithmetic for the totals and CDF; the
cost and accuracy must be measured over the 364 values. Compensating the sum alone does not remove
errors in the distances, exponentials or representation of the 53-bit variate.

On GPUs with neural accelerators (MLX's `is_nax_available`: macOS 26.2 or later and a recent GPU
generation), MLX routes float32 matrix products to the accelerators only when TF32 is enabled, which
is the default, and never routes complex types there (`backend/metal/matmul.cpp`). Complex64 products
therefore run at full float32 precision without the accelerator. Splitting into real products can use
the accelerator at TF32 precision; with `MLX_ENABLE_TF32=0` they keep full precision and skip it. On
GPUs without accelerators the setting does not affect ordinary matrix products. Benchmark complex64,
real products at TF32 and real products at full precision, and gate each against the CPU
([MLX numerical precision](https://ml-explore.github.io/mlx/build/html/usage/precision.html)).

The October CPU matrix-product layout changed distances by at most 4.1e-7 relative and the best rotation
in 3 of 320,000 pairs. That result does not establish the accuracy of GPU posterior draws. Compare
normalized CDFs and the variate representation separately, then explain each changed draw using the
measured displacement of a CDF boundary. Bit-identical CPU/GPU draws are not required at those
boundaries; repeatability on one fixed GPU configuration is a separate gate.

## Hot spot 2: reconstruction and reference projection

These are two mirror-image operations on the same 3D Fourier grid. Insertion scatters each particle's
2D Fourier plane into the 3D grid. Reference projection gathers 2D central sections, on polar
coordinates, from the 3D grid. Both use a Kaiser-Bessel interpolation window of 3 x 3 x 3 voxels.

### Reference projection (gather)

For every reference direction, SIMPLE interpolates the polar samples from a 2x padded 3D transform.
Each output sample reads 27 voxels; no two outputs write the same place. This is a natural GPU kernel
without write races. Repeatability and CPU/GPU numerical agreement still need separate checks. The
reprojection model that contains it is 11% of the ribosome run,
including file reading and writing and the volume filtering. It also feeds the search directly: run on
the GPU, it writes the polar references, then `mlx_fft_fft` and `mlx_fft_rfft` along the angular axis
form `B` and `D` in the reordered layout the matrix product reads. The references then never leave
GPU-visible memory. The projector already has an OpenMP target-offload branch for this loop, which
shows the kernel's shape.

### Insertion (scatter)

For each particle and each symmetry copy, every Fourier component of the plane adds a weighted complex
value (CTF times data) and a real weight (CTF squared) to the 27 voxels around its rotated position.
SIMPLE's CPU version avoids races without atomics. Threads take lines of one plane that are far enough
apart never to touch the same voxel, so results are identical at any thread count. A GPU must insert
many planes at once, and planes in different orientations overlap everywhere. The options are:

- **Floating-point atomic additions.** Fast, but the order of the additions changes from run to run,
  so the last bits change. This breaks SIMPLE's thread-identity gate unless the gate is relaxed to a
  tolerance for the GPU route.
- **Fixed-point integer atomic additions.** Integer addition does not depend on order, so results are
  reproducible with a fixed quantization rule and proven bounds against overflow. The required dynamic
  range may demand 64-bit addition; support and a safe scaling scheme must be established first.
- **Tiled accumulation.** Each threadgroup owns a block of the grid in threadgroup memory, gathers the
  contributions that fall inside it, then writes the block once. Determinism also requires a fixed
  contribution order and reduction tree within each block; exclusive ownership alone is insufficient.

The measured insertion cost is small today: 1.7% of the ribosome run. It grows with the number of
poses inserted per particle (several weighted poses per particle with `topk`, up to a 64-point in-plane
quadrature with the in-plane posterior). It is worth porting only if the profile below says so.

### What the earlier GPU insertion teaches

SIMPLE has a batched GPU insertion prototype written by Cyril Reboul in May and
June 2026, for NVIDIA GPUs through OpenMP target offload (`volume/simple_reconstructor_openmpoffload.f90`,
`calc_3Drec_gpu`). It ran outside the production call chain; its timings are not recorded in the tree.
Its batching and overlap provide a starting point for Apple Silicon:

- **Coarse work items over a batch of planes.** Up to 500 prepared particle planes are packed into
  contiguous arrays (complex plane, CTF-squared plane, two rows of the rotation matrix, half flag). One
  GPU thread takes one (particle, symmetry copy) pair and inserts the whole plane. There is no
  per-voxel scheduling.
- **The grids stay resident.** The even and odd grids go to the device once per pass and come back once
  at the end. Per batch only the planes and rotations travel.
- **CPU and GPU overlap.** The kernel is launched without waiting. While it runs, the CPU threads read
  and prepare the next batch of particles; then the host waits and swaps. This is the same two-stage
  pipeline this plan proposes for the search.
- **A specialised kernel.** Fixed 3-tap window, a polynomial form of the Kaiser-Bessel weight, the
  rotated position advanced incrementally along each line, Friedel symmetry so that only half of each
  plane is stored, and zero components skipped.

With suitable shared allocations on Apple Silicon, the grids and packed planes can be accessed by both
CPU and GPU without the explicit device transfers. Packing, synchronization and buffer ownership remain.

One part must not carry over. Threads working on different planes write the same voxels at the same
time, and the kernel has no atomic operations or other protection. Updates are lost, mostly at low
frequency, where every plane passes through the same voxels, and the result changes from run to run.
The loss is partly hidden because lost contributions drop out of the data sum and the weight sum alike,
which the division by the weights then evens out; this is an inference from the code, not a
measurement. The colouring loop inside the kernel (every third line) does not prevent this: within one
thread it runs in sequence, and it cannot separate different planes. The Apple version keeps the
structure and replaces the unprotected update with one of the three options above: atomics, which are
fast but order-dependent; fixed-point integer atomics; or tiled accumulation.

### What the earlier GPU reference projection teaches

The same offload work produced the reference projector that production uses today,
`fproject_polar_batch_opt` (`image/simple_projector_pft_batch.f90`, called through
`vol_pad2ref_pfts_opt` once per state and half). It is one loop with two directive sets: built with
offload it runs on an NVIDIA GPU, otherwise the same loop runs on the CPU threads. Its features are the
ones a Metal kernel needs:

- **One work item per output sample.** The loop runs over (direction, shell, angular sample), collapsed
  into one index space. Each item computes one polar Fourier sample of one reference and writes it to
  its own place. There are no write conflicts or reductions between work items. Each item still sums
  27 contributions; compiler contraction and arithmetic choices can change its last bits across devices.
- **Everything precomputed into flat tables.** The polar sample coordinates, scaled to the padded grid,
  are one table `(2, pftsz, nk)`. Each direction contributes the first two rows of its rotation matrix,
  `(2, 3, ndirections)`. The kernel does six multiplications to place a sample, then a 27-voxel
  weighted sum.
- **A device-callable weight function.** The Kaiser-Bessel weight is a pure function marked for the
  device (`apod_fast_device`, a polynomial), with per-axis normalisation of the three weights. The
  Metal kernel needs the same function, with the same coefficients, to match the CPU.
- **Half the work by symmetry.** With mirror projection on (`mirr_proj=yes`), only one direction of
  each mirror pair is projected. The other is derived from it by a mirror operation on the polar
  arrays (`mirror_ref_bank`).

What changes on Apple Silicon is the data movement. The offload version copies the whole expanded,
padded volume to the device and the references back on every call: twice per state, every iteration.
With suitable shared buffers both explicit device copies can be avoided. The output can also feed the
next step on the GPU directly: the angular FFTs that form `B` and `D`, and the reordering into the
layout of the matrix product. The references are then made, transformed and scored without leaving
GPU-visible memory. The CPU keeps the
steps before the projection (reading, masking and filtering the volume, padding and transforming it),
unless the profile shows they matter.

This kernel is a low-risk first port. It is a pure gather without write races, has a CPU twin in
production to compare against, and feeds the search. The steps below therefore move it ahead of
reconstruction.

### Assembly, the final map and the iterative route

Assembly sums the partial grids, computes Fourier shell statistics, applies the Wiener filter, divides
by the sampling density, runs one 3D inverse FFT and corrects for the interpolation window. Except for
the FFT these are simple loops over voxels. The final reconstruction repeats this at full box size (up
to about 512^3). The iterative route (`rec_backend=pcg`) applies a Toeplitz operator built on 3D FFTs at
each iteration. All of these are candidates for `mlx_fft_rfftn` and element-wise MLX operations.
Part of the assembly time on the ribosome run is file reading and writing in the distributed mode,
which a shared-memory run on one Mac does not have.

## Integration with SIMPLE

- **Build.** A CMake option, on macOS arm64 only, links mlx-c (version 0.8.0, built against MLX
  0.32.4). A Fortran module declares the mlx-c functions with `iso_c_binding`. mlx-c already exposes a C
  ABI with opaque handles and integer error codes, so no C++ reaches Fortran. Callbacks for
  `mlx_compile` are plain C function pointers (`mlx_closure_new_func`), so a Fortran `bind(c)` procedure
  passed with `c_funloc` serves. The Metal kernel sources are text, kept as separate `.metal` files and
  read at start-up. The x86 Linux build is unaffected and stays the CPU performance target.
- **Memory.** `mlx_array_new_data_managed` reaches the core constructor, which attempts to wrap the
  pointer as a Metal buffer. If wrapping returns no buffer, the constructor copies the data and calls
  the supplied destructor immediately; otherwise the destructor runs when the wrapped memory is
  released. The lifetime contract must handle both paths without premature freeing or double freeing.
  Metal requires a page-aligned address and extent within one VM region
  ([buffer requirements](https://developer.apple.com/documentation/metal/mtldevice/makebuffer%28bytesnocopy%3Alength%3Aoptions%3Adeallocator%3A%29)).
  MLX passes the array's logical byte count to Metal, so aligning only the allocation address is not
  sufficient; the wrapping experiment must establish a supported padded-storage/view arrangement.
  Check pointer identity and destructor timing. `mlx_get_active_memory` counts successfully wrapped
  buffers too, so memory counters and timing alone cannot prove that wrapping was copy-free.
  Fortran writes only after all readers have completed and the ownership protocol permits reuse.
- **Execution.** One GPU stream runs the scoring pipeline. Invoke the pipeline for each new batch,
  submit both result arrays with `mlx_async_eval`, then run CPU work on the preceding batch.
  `mlx_array_eval` on the outputs or `mlx_synchronize` marks completion before the CPU consumes them.
  Reuse a compiled callable per fixed tile shape if `mlx_compile` improves dispatch overhead; do not
  expect an already evaluated output graph to recompute after a host write to its input memory.
  Release obsolete graphs and handles before recycling buffers, and keep outputs alive while the CPU
  holds pointers into them. Measure overlap under CPU load as well as isolated GPU throughput.
- **Prior art in the tree.** SIMPLE keeps GPU code for discrete NVIDIA GPUs, by the maintainer's
  ruling. `gen_many2many_euclids_cufft` (OpenMP target offload with
  cuFFT) is a benchmark of the unshifted dense contraction: it uses the even references only and has no
  draw. `simple_reconstructor_openmpoffload.f90` is the batched offload insertion described under
  hot spot 2; neither has a production caller. The production reference projector has an offload
  branch. The Apple route is written beside this code and does not replace it.

## Correctness gates

Every step passes gates before it is wired in. Numerical agreement between backends and repeatability
within a backend are separate requirements. The tolerances below are proposed initial acceptance
targets, not measured GPU results.

- **Score agreement.** Add the GPU layout to `simple_test_exec test=polar_matching_bench` after all
  CPU layouts include the last angular bin and identical keyed variates. Compare complete distance
  rows with both current production and the CPU matrix-product layout. Initial target for valid
  distances: `abs(gpu-cpu) <= 1e-6 + 1e-6*abs(cpu)`. Check cutoff decisions, invalid rows, both halves,
  seed shifts and the angular Nyquist bin explicitly. On the fixed benchmark, changes of the best
  rotation must be at most 1 in 10,000 pairs; report the competing minima's gaps as well.
- **Draw agreement.** Hash words must match bit for bit. Check flat rows, tied minima, temperature
  zero, sharp posteriors, invalid entries, the exponential cutoff and variates near CDF boundaries.
  Initial target for the normalized CDF is `max(abs(F_gpu-F_cpu)) <= 1e-5` per row. Measure variate
  error separately. Every changed draw must have its CPU variate within the measured CDF error plus
  variate error of a CPU CDF boundary; count all changes and report the maximum errors. This measures
  the effect of the entire scoring and draw pipeline, including exponentials and accumulation.
- **Projection and reconstruction agreement.** Compare polar samples, data grids and weight grids
  separately with CPU results, using mixed absolute/relative tolerances to handle zeros. For projection,
  start with `abs(gpu-cpu) <= 1e-6*S + 1e-6*abs(cpu)`, where `S` is the maximum CPU sample magnitude
  for the state and half. Fix reconstruction tolerances from CPU reference and independent synthetic
  fixtures before qualification. Compare map scale and voxel errors as well as shell correlation;
  correlation alone would miss a uniform scaling error.
- **Workflow, small simulated sets** (up to 2,000 particles, known orientations): orientation errors
  against the ground truth match the CPU route within its spread across prescribed seeds. Resolution
  estimates are not used at this size; they are too noisy.
- **Workflow, large real sets** (the ribosome and beta-galactosidase benchmarks): final resolution and
  maps match the CPU route within its spread across prescribed seeds. Agree on the seed set and
  acceptance ranges before comparing GPU runs.
- **Repeatability.** For fixed inputs, seed, device, software versions and GPU tile shapes, repeated
  runs must be bitwise identical. Test separately across CPU thread counts (including 11 and 24),
  host batch boundaries, short tiles and distributed partitions, preserving global draw keys and
  fixed GPU arithmetic shapes. Do not infer these properties from repeatability at one shape.
  CPU/GPU comparisons use the numerical and scientific gates above. Any relaxation for floating-point
  atomic reconstruction requires an explicit change to this contract before that route is adopted.

## Steps

**Profile on the target machine.** Run the ribosome and beta-galactosidase benchmarks in
shared-memory mode on a Mac Studio with SIMPLE's bench timers on. Split the search into dense score,
shift screen and shift refinement, and reconstruction into insertion, assembly, final map and
projection. Done when the numbers fix the order of the later steps.

**CPU reference and review kit for Apple.** Correct the benchmark's missing angular bin and its
`ran3()` draws, validate against current production, then export a fixed case with binary input arrays,
shapes, layouts, normalizers, temperatures, draw keys, variates, full distance rows and CPU results.
Supply a short C program on mlx-c that runs the tiled dense score alone, with separate timings for
packing, wrapping/copying, matrix products, FFT input conversion, FFT and drawing. Record hardware,
library revisions, precision settings and peak memory. The kit lives outside the SIMPLE repository
and is prepared before production integration so MLX engineers can reproduce and tune the workload
without building SIMPLE.

**Binding and memory.** The Fortran module for mlx-c. Wrap one reference array, run one matrix product,
read the result back in place. Verify aliasing through pointer identity and destructor timing, exercise
the copying fallback, and check buffer reuse over asynchronous submissions. Use memory counters and
timings to quantify the cost, not to establish aliasing. Done when the ownership protocol is verified
and copy-free wrapping is confirmed, or the copy's cost is measured and accepted.

**Dense score in the benchmark.** Add the matrix product, inverse FFT and distance-and-draw kernel as
a layout of the corrected `polar_matching_bench`. First use the library FFT and a separate draw kernel;
test fusion if the profile justifies it. Tile both routes, including their intermediate lifetimes.
Done when score, draw and repeatability gates pass. Report complete pipeline time per pair, peak memory
and packing cost against the refreshed CPU matrix-product layout on the Mac and current production on
the x86 test machine. Keep isolated-kernel and end-to-end timings separate.

**Dense score in production.** Split the table step into its batch-wide and per-particle phases,
double the particle memo, and wire the GPU route in with the two-stage CPU/GPU pipeline. Done when
the workflow gates pass.

**Reference projection on the GPU.** A Metal port of `fproject_polar_batch_opt`: the same flat tables,
one work item per polar sample, the same weight polynomial. It writes directly into the reordered
arrays the dense score reads after angular FFTs. Done when the projection agreement and repeatability
gates pass. It needs no change to the search pipeline, so it can also be done earlier, right after the
binding step, as the first kernel to exercise the build, the wrapping and the custom-kernel route.

**Reconstruction.** Start from the structure of the earlier offload insertion: batches of packed
planes, one thread per plane and symmetry copy, resident grids, CPU preparation of the next batch
overlapped with the kernel. Replace its unprotected scatter with a safe one. Then the 3D FFTs and the
iterative route, in the order the profile gives.

**Shift screen.** The Fourier-Bessel sums and the batched 2D FFT, if the profile shows it matters once
the dense score is on the GPU.

## Questions for the MLX team

The first three questions define the initial experiment; the others support tuning and later ports.

1. **Dense-scoring layout.** MLX already implements complex64 GPU matrix products, using three real
   matrix products within its complex kernel. For 183 products of `(600 x 122) @ (122 x 512)` per half,
   what tile shapes and layouts would you start with? On GPUs with neural accelerators, complex types
   are excluded from the accelerator route and float32 uses it only with TF32. Would you keep complex64
   at full precision, or split into real products and accept TF32 for the contraction?
2. **External-buffer ownership.** The core attempts copy-free wrapping and has a copying fallback.
   What allocation, padded-storage/view and lifetime pattern do you recommend for Fortran-owned
   buffers? What host-write and graph-reuse protocol is supported between asynchronous evaluations?
   Would MLX-owned storage exposed to the host be preferable?
3. **Numerical agreement and repeatability.** What reproducibility guarantees apply to GPU matrix
   products and FFTs on a fixed device and software version? Can we preserve them across host batch
   partitions using fixed padded tiles? What precision controls and error checks would you recommend
   for the posterior draw, including its exponentials, 53-bit variate and cumulative sum?
4. **FFT layout and fusion.** The GPU planner supports the factors 2, 7 and 13 in our 364-point transform.
   How does that length perform at this batch size, including the copy to make bins contiguous? Would
   you fuse the inverse transform and draw, or retain the library FFT? What are the threadgroup-memory,
   register and input-layout tradeoffs? Other angular grid sizes are possible.
5. **Custom kernels from C.** Is `mlx_fast_metal_kernel_new` / `mlx_fast_metal_kernel_apply` the intended
   route for production kernels called from a host program in another language? How are compiled
   kernels cached? What exactly does `atomic_outputs` enable? Is 32-bit float atomic addition fast on
   current GPUs? We understand that
   Metal offers 64-bit atomics only as minimum and maximum, not addition; is that right, and what do you
   recommend for deterministic accumulation into a 3D grid?
6. **Asynchronous execution.** Does submitting both outputs with `mlx_async_eval`, running the CPU
   phase and waiting at buffer handover provide the intended overlap? Which synchronization and
   handle-lifetime practices avoid hidden serialization or premature buffer reuse?
7. **Overhead and compilation.** What is the per-operation dispatch cost when a graph is built one call
   at a time through mlx-c? Is `mlx_compile` with C callbacks the recommended way to build a fixed
   callable per fixed tile shape, invoked with each new batch?
8. **CPU and GPU together.** With the CPU cores busy on their own work, how much does the GPU kernel
   slow down from memory bandwidth contention? Should multi-gigabyte, long-lived arrays be kept resident
   with `mlx_set_wired_limit`?
9. **Large 3D FFTs.** How do `mlx_fft_rfftn` and `mlx_fft_irfftn` perform on 3D grids up to 512^3, and how
   much extra memory do they need?
10. **Accumulation precision.** For reconstruction and posterior sums that exceed ordinary float32
    accuracy, would you recommend compensated summation, two-float arithmetic or CPU execution?
    Are there planned capabilities that should influence this boundary?
11. **Profiling.** Which tools show why a custom kernel underperforms on M-series GPUs? Is
    `mlx_metal_start_capture` with Xcode's GPU tools the usual route?

## MLX implementation references

The proposed build pairs mlx-c 0.8.0 with MLX v0.32.4, the tag selected by mlx-c's CMake configuration.
The local core checkout inspected on 2026-11-16 was `fc344c986c4169d99f76b5d95ae6b4a58a359a63`,
82 commits beyond that tag; mlx-c was `a26f2d584a14ec619ae8db9a09be03becd47598a`. The wrapping path,
complex matrix kernel and supported FFT radices described here are also present at v0.32.4. Pin and
report the actual built revisions for the review kit and all measurements.

- [C managed-array wrapper](https://github.com/ml-explore/mlx-c/blob/a26f2d584a14ec619ae8db9a09be03becd47598a/mlx/c/array.cpp)
  and [core constructor with copy fallback](https://github.com/ml-explore/mlx/blob/v0.32.4/mlx/array.cpp).
- [Metal external-buffer allocation and accounting](https://github.com/ml-explore/mlx/blob/v0.32.4/mlx/backend/metal/allocator.cpp).
- [Complex matrix kernel](https://github.com/ml-explore/mlx/blob/v0.32.4/mlx/backend/metal/kernels/steel/gemm/mma.h).
- [GPU FFT planner and input-layout copy](https://github.com/ml-explore/mlx/blob/v0.32.4/mlx/backend/metal/fft.cpp).
- [C asynchronous evaluation](https://github.com/ml-explore/mlx-c/blob/a26f2d584a14ec619ae8db9a09be03becd47598a/mlx/c/transforms.h).

## Appendix: code map

Paths are under `src/main/`.

| what | where |
| --- | --- |
| per-pair score, all rotations | `pftc/simple_polarft_corr.f90`, `gen_euclid_crvec` and `euclid_shell_sum` |
| memoized particle and reference spectra | `pftc/simple_polarft_memo.f90` |
| polar array layout | `pftc/simple_polarft_calc.f90`, `pftc/simple_polarft_core.f90` |
| table step, dense score, shift screen, refinement | `strategies/search/probabilistic/simple_eul_prob_tab.f90`, `fill_tab_range` |
| tempered draw, keyed random numbers | `strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`, `draw_posterior_row`; `src/utils/math/simple_rnd.f90` (outside `src/main/`) |
| shift map | `pftc/simple_polarft_shiftmap.f90` |
| shift refinement (L-BFGS-B, double precision) | `pftc/simple_pftc_shsrch_grad.f90` |
| batch loop | `strategies/search/simple_strategy3D_matcher.f90`, `pipeline_batch` |
| plane preparation for insertion | `image/simple_image_ctf.f90`, `gen_fplane4rec` |
| insertion | `volume/simple_reconstructor.f90`, `insert_plane_oversamp` |
| assembly | `strategies/parallelization/simple_volassemble_strategy.f90` |
| iterative route | `volume/simple_reconstructor_pcg.f90`, `opt/simple_pcg_solver.f90` |
| reference projection | `image/simple_projector_pft_batch.f90` |
| kernel benchmark | `commanders/test/simple_polar_matching_bench.f90` |
| GPU offload code for NVIDIA (kept, no callers) | `pftc/simple_polarft_corr.f90`, `gen_many2many_euclids_cufft`; `volume/simple_reconstructor_openmpoffload.f90` |
