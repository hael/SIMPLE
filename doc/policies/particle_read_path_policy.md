# Particle Read Path Policy

How SIMPLE reads particle images from stack files: many small per-micrograph
stacks in a stream project, or a few large ones. The layout on disk is never
changed for reading's sake. No particle is rewritten into another file, and no
cache of particles exists. The design and its measurements are in
`doc/implementation_notes/completed/particle_io_few_particle_stacks.md`.

## Bit-exact

A consumer receives exactly the bytes and values it always received, in every
pixel mode. Only how the bytes are fetched changes. Any change to the read path
is tested image for image against a single read of each image
(`simple_stack_io_tester.f90`).

## The read session

Every pass that reads particles reads them through one read session
(`stack_read_session`, `src/fileio/simple_stack_read_session.f90`). A pass is a
search, a reconstruction, a probability table, power spectra, class sums or a
FLEX batch loop. In the matchers the session is wrapped by `ptcl_read_session`
(`simple_matcher_ptcl_io.f90`). Nothing reads a particle by opening its stack
on its own.

- **One open per stack per pass.** The session decodes a stack's header once
  per open. The stack read last in a batch stays open for the next batch.
- **Plain read-only descriptors.** Each read stream opens, reads and closes its
  own stacks inside the parallel region. The opens run concurrently; there are
  no Fortran units and no existence checks before an open.
- **Coalesced spans.** The images needed from one stack are read in spans, one
  positioned request each. A span may bridge at most `SPAN_GAP_BUDGET` (2 MiB)
  of unneeded bytes between needed images. The raw scratch of a session, over
  all streams, is capped at `SPAN_SCRATCH_CAP` (64 MiB).
- **The prefetch window.** Told the next batch, the session opens up to
  `PREFETCH_MAX_STACKS` (64) of its stacks. A producer thread then reads up to
  `PREFETCH_WINDOW_BYTES` (128 MiB) ahead while the caller computes.
  Prefetching changes no byte read.
- **The stream count** is a property of the storage, not of the CPU count
  (`io_read_nstreams`): it is probed per device, and the existing
  `SIMPLE_IO_NSTREAMS` overrides it.
- **Pixel modes.** The descriptor path decodes bytes, 16-bit integers, float32
  and float16. Any other mode is read through the image file reader
  (`dstack_io`), one image at a time.
- **Failures stop the run.** A missing or truncated stack, a short read, a box
  that differs from the pass's, an index outside the stack, or running out of
  descriptors each stop the run, naming the stack. Interrupted reads are
  retried.

Stack rows that share one file are grouped by file, so the file is opened once
for all of them (`project_stack_indexing_policy.md`, reader policy).

## The particle buffer

A 3D search fills a pass-sized buffer of decoded images before its batch loop.
The reconstruction of the same pass, which runs after the assignment barrier,
copies its images from that buffer instead of reading them again. The buffer
lives for one pass in one process
(`separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`).
In shared-memory 2D, the probability table (`prob_tab2D`, run in the master
process) fills the buffer with the iteration's sampled particles, the
matcher of the same iteration copies its batches from it, and the matcher
releases it. Batches are copied out of the buffer in parallel.
Its budget is half the task memory a queue system grants; otherwise it is the
node's physical memory (`MemTotal` on Linux, `hw.memsize` on macOS) divided by
twice the number of processes the node runs.
A share above the budget, or an allocation that fails, leaves the buffer empty,
and the batches are read from the stacks as before. The outputs are the same
either way.

## Benchmarks

With `L_BENCH_GLOB`, every file of a reading pass reports the following for
that pass: particle reading time, apart from preparation; stacks opened; read
requests; bytes requested, read and prefetched; prefetched bytes used; and
images copied from the buffer. The refine3D search file also reports the
buffer's budget, its size and its outcome.

A gain is claimed only from the same case measured before and after.
