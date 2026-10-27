# BLAS threading: an OpenMP BLAS on Linux

Implementation note, 2026-10-26. **Implemented in `cmake/Dependencies.cmake`; the configure logic was simulated
with stand-in libraries, not yet configured or built against a real toolchain.**

## Why

SIMPLE calls BLAS in two ways:

- from serial code, expecting the library to use the cores: the picker's correlation products
  (`simple_pickref_corr_batch`, `gemm_tn`), LAPACK eigen/SVD/inverse in `simple_linalg`, the wide products of
  `trajectory_denoise pca_mode=flex` when the BLAS threads;
- inside its own OpenMP teams, expecting the library to stay on the calling thread: the flex probe E-step
  (`sgemv`, `dgemv`, `ssyrk`), the per-thread products of `trajectory_denoise pca_mode=flex`.

A serial BLAS (RHEL's plain `libopenblas`, reference BLAS) runs the first kind on one core; a pthreads OpenBLAS
(Ubuntu's default `libopenblas0-pthread`, the bundled build until now) oversubscribes the cores in the second kind
unless `OPENBLAS_NUM_THREADS=1` is set, which then makes the first kind serial. An OpenMP build of OpenBLAS does
both right: outside a team it uses the OpenMP thread count (SIMPLE's `nthr`), inside one the calling thread.

## What configure does (Linux, GNU C and Fortran, OpenMP on)

1. A system OpenMP OpenBLAS first: `libopenblaso` (Fedora/RHEL/Rocky `openblas-devel`) or
   `/usr/lib/<arch>/openblas-openmp/libopenblas` (Debian/Ubuntu `libopenblas-openmp-dev`), used for BLAS and
   LAPACK when a probe program links both from it and `openblas_get_parallel()` returns 2.
2. Otherwise `find_package(BLAS/LAPACK)` as before, then the probe on what was found. Kept: OpenBLAS OpenMP
   build, MKL, FlexiBLAS with an OpenMP backend. Anything else (serial or pthreads OpenBLAS, a BLAS the probe
   cannot classify): the bundled OpenBLAS, now the OpenMP build, and ARPACK built against it (a system ARPACK
   would load the system BLAS too). Configure prints why and which package would avoid the build.
3. The bundled OpenBLAS (also when no BLAS is found) is the OpenMP build whenever C and Fortran are GNU (one
   libgomp); its own build folder (`_deps/openblas-omp`), so existing build trees rebuild it, linked by full path
   so the binaries relink.
4. `SIMPLE_BLAS_TEAM_SAFE` is defined when the linked BLAS is known to stay on the calling thread inside a team
   (system or bundled OpenMP OpenBLAS, MKL, FlexiBLAS OpenMP); the dependency summary prints it.

Left alone: macOS (Accelerate), other compilers, cross-compiling, an explicit `BLA_VENDOR` or `BLAS_LIBRARIES`,
`-DSIMPLE_BLAS_OPENMP=OFF`. CI installs `libopenblas-openmp-dev` instead of `libopenblas-dev`, so CI does not
build OpenBLAS.

## Verified

Configure simulated in a C-only project (`gcc` 13, CMake 3.28, Ubuntu's `--as-needed` linker default) with
stand-in shared libraries: the probe returns `openblas:0/1/2`, `flexiblas:<backend>`, `other`, and empty when
LAPACK is missing; no BLAS -> bundled OpenMP OpenBLAS and ARPACK; serial OpenBLAS -> bundled OpenBLAS and ARPACK
with the explanation; `libopenblaso` in the prefix -> used directly, nothing bundled; `SIMPLE_BLAS_OPENMP=OFF` ->
the serial system BLAS kept; re-configure gives the same result.

## Still to do

- A real configure and build on Linux (the maintainer's RHEL-type box: expect `libopenblaso` to be picked if
  `openblas-devel` is installed) and on macOS (unchanged path).
- A run of the picker and of `trajectory_denoise pca_mode=flex` with the new BLAS, watching core use.
- Whether to also use `SIMPLE_BLAS_TEAM_SAFE` elsewhere (the flex E-step already relies on it implicitly).
