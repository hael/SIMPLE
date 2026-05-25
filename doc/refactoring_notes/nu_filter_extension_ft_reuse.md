# NU filter high-resolution extension FT reuse

## Summary

The high-resolution nonuniform filter extension currently regenerates each
challenged shell as an independent filtered even/odd pair. This is safe and
simple, but it repeats the same invariant preparation for every shell:

- copy the input even/odd halfmaps
- taper the real-space volume edges
- transform both halfmaps to Fourier space
- apply one shell-specific Butterworth filter
- inverse transform and write the scratch cache pair

A future CPU optimization could retain the prepared Fourier-space even/odd
volumes for the duration of a shell-walk session. Each challenged shell would
then copy from those prepared FT volumes, apply only the new shell filter,
inverse transform, and proceed with the existing objective calculation.

## Proposed Shape

1. Add an extension-session helper local to the NU filter implementation.
   It owns the tapered, Fourier-transformed even/odd input volumes and is
   created by the multi-shell extension loop.

2. Add a generator path that takes the prepared FT pair plus a `cutoff_find`.
   It should produce the same filtered scratch pair as the current
   `generate_single_filtered_pair` path.

3. Keep the current single-shell generator as a fallback.
   `extend_nu_filter_highres_shell_next` should remain usable without an
   explicit session object.

4. Preserve cache safety.
   Any challenged shell should still delete stale `nu_filter_cache_*_k_*.mrc`
   files before writing fresh scratch products. The optimization is about
   avoiding repeated FT preparation, not about trusting old cache files.

5. Gate the optimization by memory policy.
   Keeping two prepared FT volumes alive may improve CPU saturation during long
   shell walks, but it increases peak memory. For fine-sampled, distributed
   refinements where memory limits MPI part counts, the current lower-memory
   path may be preferable.

## Tradeoff

This is primarily a parallel throughput optimization, not a RAM optimization.
It should reduce repeated copy/taper/FFT work and help long high-resolution
extension walks keep the CPUs busier. It should not be enabled blindly for
memory-bound jobs such as very fine-sampled ApoF-like cases.
