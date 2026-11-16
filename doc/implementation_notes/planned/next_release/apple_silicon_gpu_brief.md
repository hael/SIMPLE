# Accelerating cryo-EM image matching with MLX

For discussion with Samy Bengio on 20 November 2026

**Why this matters.** Cryo-electron microscopy (cryo-EM) determines 3D structures of proteins from
hundreds of thousands of noisy 2D images. Today this needs GPU clusters. We maintain SIMPLE, an
open-source cryo-EM program, and we would like a high-memory Mac Studio to do the same work on a
desk. The step that dominates our run time looks like a fit for MLX.

**What we are asking for.** If one engineer could look at the three questions below, and later at a
small standalone program that reproduces the workload, that would be ideal. Written answers are fine.

**The computation.** We compare every image in a batch with every one of up to 10,000 reference images,
at 364 rotation angles. Per batch this is three steps:

1. A batched complex matrix product: 183 independent products of `(600 x 122) @ (122 x 512)` in
   complex64, one per frequency, for each tile of 600 images and 512 references.
2. An inverse real FFT of length 364 along the frequency axis, for every (image, reference) pair. This
   needs the product's output transposed, so that each pair's 183 values are contiguous.
3. A small per-pair kernel: find the minimum of the 364 values, form softmax weights with a given
   temperature, and sample one angle using a uniform random number that we supply. Output: one value
   and one index per pair.

The CPU cores then do irregular follow-up work on the best pairs of one batch while the GPU processes
the next. Our code is Fortran; we would call MLX through mlx-c.

**What we know so far.** On the CPU, writing step 1 as matrix products made this step 3.3–3.8 times
faster on an M3 Pro, and the CPU version is limited by arithmetic, not memory bandwidth. In one
large benchmark, the stage containing this matching took 46% of the total run time, including some
file and process overhead; we will measure the matching alone on a Mac Studio. MLX supports complex64
matrix products on the GPU, and its GPU FFT handles 364 = 4·7·13 natively.

**What concerns us.**

- **Memory.** One tile's product output, the transposed copy for the FFT and the FFT output take about
  1.35 GB together, so we must tile. The references take 3.6 GB and stay resident.
- **Sharing buffers.** We want the GPU to read arrays that Fortran allocated, without copying. MLX can
  wrap external memory but falls back to a copy when the buffer does not meet Metal's requirements.
- **Numbers.** We need results that agree closely with the CPU, and that repeat exactly on one machine.
  The sampling step is sensitive: a tiny difference in the weights can pick the neighbouring angle.

**Three questions**

1. Which layout and tile sizes would you choose for steps 1 and 2? Should the FFT and the per-pair
   kernel be fused into one custom Metal kernel? On GPUs with neural accelerators, MLX uses them for
   float32 only with TF32 and never for complex types. Would you keep complex64 at full precision, or
   split into real products and accept TF32?
2. What is the supported way to share Fortran-allocated buffers with mlx-c, and to overlap GPU work on
   one batch with CPU work on another (submission, waiting, buffer reuse)?
3. Are GPU matrix products and FFTs bit-for-bit repeatable on one machine and software version? What
   would you do for the softmax and the cumulative sum of the sampling step, given no double precision
   on the GPU?

**Next step on our side.** We will export one real case and write a short C program on mlx-c that runs
the three steps alone, with our CPU results to compare against, timings and peak memory. Your team
could run it without building SIMPLE.

The detailed technical proposal is attached (`apple_silicon_gpu_port.pdf`).
