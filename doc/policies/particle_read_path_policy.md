# Particle Read Path Policy

How SIMPLE reads particle images from stack files: many small per-micrograph
stacks in a stream project, or a few large ones. The layout on disk is never
changed for reading's sake: reading rewrites no particle into another file. The
2D particle cache (`2D/particle_cache_policy.md`) and `consolidate_particles`
(`project_stack_indexing_policy.md`) are separate products, made on request.
The design and its measurements are in
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
FLEX batch loop. Each wraps the session in `ptcl_read_session`
(`simple_matcher_ptcl_io.f90`). None of these passes opens a particle stack on
its own. Outside the session are the reads of the 2D particle cache
(`ptcl_cache_read_batch`) and the single-image utilities: `ctf_correct`
(`read_imgbatch`), `write_classes`, the random and time-series reference
selection of `cls_init=ptcl` (`simple_procimgstk.f90`), the project thumbnails,
and the byte copy of `consolidate_particles`. The session takes pairs of stack
name and image index; `ptcl_read_session` resolves a particle to them with
`sp_project%get_stkname_and_ind` (`project_stack_indexing_policy.md`, reader
policy).

- **One open per stack per pass.** The session decodes a stack's header once
  per open. The stack read last in a batch stays open for the next batch.
- **Plain read-only descriptors.** Each read stream opens, reads and closes its
  own stacks inside the parallel region. The opens run concurrently; there are
  no Fortran units and no existence checks before an open.
- **Coalesced spans.** The images needed from one stack are read in spans, one
  positioned request each. A span may bridge at most `SPAN_GAP_BUDGET` (2 MiB)
  of unneeded bytes between two needed images, holds no more unneeded images
  than needed ones, and stays within its stream's share of the scratch. The raw
  scratch of a session, over all streams, is capped at `SPAN_SCRATCH_CAP`
  (64 MiB).
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
(`3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`).

When the pass reconstructs at a cropped box (`box_crop < box`, a polar
search), the buffer holds records instead of decoded images. The search
reads its batches from the stacks. Each particle's alignment preparation
normalizes and transforms the full box. The unshifted Fourier crop of that
transform, back in real space, is the observation the reconstruction would
make from the same image (`prep_rec_observation`, before its taper). The
batch keeps it at `box_crop` in the particle's slot. The reconstruction (both
backends) takes its observations from these records. It does not read the
stacks or transform the full box again. The outputs are bit for bit those of
the raw buffer, and the buffer is `(box/box_crop)^2` smaller. One exception:
the PCG backend's own observation image transforms with several threads when
`box_crop` is 200 or more, and a record is transformed on one thread, so there
the two may differ by rounding. A particle without a record is read from the
stacks.

In shared-memory 2D without the particle cache, the probability table
(`prob_tab2D`, run in the master process) fills the buffer with the
iteration's sampled particles, the matcher of the same iteration copies its
batches from it, and the matcher releases it. Batches are copied out of the buffer in parallel. When the 2D
particles are cropped (`box_crop < box`, uncached, outside the stream's 2D
pool), the table opens the buffer for records instead: each particle's
in-memory cache entry (noise-normalized at the full box, Fourier-cropped,
unshifted, real space at `box_crop`), put before the table's alignment
preparation transforms its copy. The matcher copies a batch's records, which
are its cropped particles for the class sums and, as cached particles, for
alignment, and neither reads the stacks nor prepares the full box. The
buffer is then `(box/box_crop)^2` smaller.
Its budget is half the task memory a queue system grants; otherwise it is the
node's physical memory (`MemTotal` on Linux, `hw.memsize` on macOS) divided by
twice the number of processes the node runs.
A share above the budget, or an allocation that fails, leaves the buffer empty,
and the batches are read from the stacks as before. The outputs are the same
either way.

## Benchmarks

With `L_BENCH_GLOB`, every file of a reading pass reports the following for
that pass: particle reading time, apart from preparation; stacks opened; read
requests; bytes requested, read, prefetched and advised to the kernel;
prefetched bytes used; and images copied from the buffer. The refine3D search
file also reports the buffer's budget, its size and its outcome.

A gain is claimed only from the same case measured before and after.
