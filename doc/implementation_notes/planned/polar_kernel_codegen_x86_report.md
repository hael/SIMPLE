# Polar kernels on x86: the report

Report of Phases 3 and 4 of the plan `polar_kernel_codegen.md` (same directory), the x86 follow-up of the Mac run
committed as ee5b68be1. Unattended run `kernel_x86` on the Dell (Intel Core Ultra 9 285: 8 performance cores with
their own L2, 16 efficiency cores in clusters of four sharing an L2; gfortran 15.2.1, `-O3 -march=native`, system
OpenBLAS 0.3.15, FFTW 3.3.5), 2026-10-22/23. The only production file changed is
`src/main/pftc/simple_polarft_corr.f90`. Nothing is committed.

## The answer

- **Why x86 was slower.** With the sums in helpers with plain array arguments, gcc can prove that two shells may be
  processed in one pass and applies unroll-and-jam to the shell loop. On the Mac that gave two shells per pass with
  de-interleaving loads and stores, a gain. On x86 with AVX2 gcc vectorizes the jammed loop with lane permutes
  (36 `vpermd` per pass), reloads of permute constants and register spills: 17 instructions per bin and shell
  against 15 for the old loop. Permutes are slow on the efficiency cores most of all, hence HEAD's euclid kernel at
  0.81x of the old one there and 0.92x on a performance core. The old loops could not be jammed, because their
  accumulator is a pointer component; they were already vectorized at full AVX2 width.
- **The form kept (c14).** The euclid sums call a helper once per shell (`euclid_shell_sum`), the cc sums call one
  once per two shells (`cc_shells_sum`). The helpers still take plain array arguments, so the compiler knows unit
  stride and no overlap (what the Mac needed), but there is no shell loop around a bin loop in one procedure left
  to jam. Each element's expression and precision are unchanged; the distance conversion of ee5b68be1 stays.
- **Its speed against the old form (f3ccb5fbd) here:** euclid 1.15x on a performance core, 1.13x on an efficiency
  core, 1.00-1.28x in teams of 12 and 24 threads; cc 1.39x and 1.31x alone, 1.22-1.49x in teams. Against HEAD the
  euclid kernel is 1.24x faster on a performance core and 1.39x on an efficiency core.
- **Identity.** Distances differ from the old form's by at most one ulp, two in the shifted-reference branch and
  in the bank's exp(-loss) values; the best rotation of a pair never changes; 6 of 320,000 lists of the 20 lowest
  rotations change order (HEAD: 1 for euclid, 4 for cc). cc distances (1 - cc): at most one ulp, 99.9 % of the
  vectors byte-identical.
- **The Mac.** Speed is judged on the Dell, which is closer to the machines SIMPLE normally runs on; the Mac only
  has to build and pass its checks (maintainer, 2026-10-23). Whether c14 keeps the Mac's euclid gain of ee5b68be1
  (its euclid now runs one shell per pass there too) is not measured.

## What changed, phase by phase

- **Phase 3** (`review/phase_3.diff` under `~/agent_runs/kernel_x86`): the plan gained Phases 3-4, the judging rule,
  a file table and the Phase 3 progress (diagnosis, every candidate's timings, identity). The benchmark driver
  `simple_commanders_test_highlevel_matchbench.f90` times its middle team at half the threads (24 -> 24, 12, 1 on
  the Dell; 11 -> 11, 5, 1 on the Mac as before). No production change.
- **Phase 4** (`review/phase_4.diff`): `simple_polarft_corr.f90` gets the c14 form: `euclid_sums` calls the new
  `euclid_shell_sum` once per shell; `cc_sums` calls the new `cc_shells_sum` once per two shells (once more for an
  odd last shell). Callers, argument passing and the distance conversion are HEAD's. This report and the plan's
  Progress complete the phase.

## Why HEAD was slower here

Read from gfortran's vectorization and loop reports (`-fopt-info-vec-all`, `-fopt-info-loop-optimized`) and the
assembly of the hot loops, with the build's flags:

- The old euclid loops vectorize 4 complex values per iteration (32-byte vectors), behind a run-time overlap check.
  The accumulator `crvec1(ithr)%c` is a pointer component of unknown stride, so it is loaded and stored 8 bytes at
  a time: 61 instructions per 4 bins.
- HEAD's `euclid_sums` contains a clean loop of 50 instructions per 4 bins, but it runs only for the odd last
  shell. For more than two shells the code enters the unroll-and-jammed version (gcc reports "applying unroll and
  jam with factor 2" for the shell loops of `euclid_sums` and both of `cc_sums`): 276 instructions per 8 bins x
  2 shells, of them 36 `vpermd` (cross-lane permutes), 16 `vblendps`, 18 reloads of permute constants and 24 stack
  accesses for spilled registers. The cc loop has 20 `vpermd` per iteration; the old cc loop none.
- No packing copies (`_gfortran_internal_pack`) in either form's hot path; the FFTW buffer is 32-byte aligned, so
  HEAD's vector stores never split a cache line; the distance conversion vectorizes without masked stores and is not
  involved (byte-identical; alone it makes the old sums 2.7 % faster on a performance core).
- What does not help: `!$omp simd` on the bin loop and `!GCC$ unroll 1` on the shell loop leave the jam in place;
  a local `pointer, contiguous` alias of the accumulator does not change gfortran's code at all (still strided).
- Writing two shells per pass by hand removes the permutes and is the fastest form on one thread (c6, c12, c13:
  euclid 1.15x), but it reads eight memory streams per pass (four arrays, two shells) and falls to 0.88-0.96x at
  12 threads. The likely reason, not proven: cores that share an L2 run out of hardware prefetch streams. Forms
  that read at most four streams per pass (one euclid shell, or two cc shells of two arrays) keep their gain in a
  team. This is why c14 treats the two objectives differently.

## The candidates in every setting

How measured: each form is HEAD's file with only the sums and their callers replaced, built in one scratch tree;
the benchmark's production layouts `per_pair` (euclid) and `per_pair_cc` on the cases 1000x1 (1,000 references)
and 2500x4 (10,000). P1 and E1: one thread pinned to cpu3 (performance) and cpu20 (efficiency), best of three after
a warm-up, median of two or three pinned rounds (repeatable within 1 %). T12 and T24: unpinned teams as production
runs them; they move 10-30 % from launch to launch, so six rounds in rotating order (c14 four, HEAD two), median.
Nothing else ran during the timings. Speed = old time / form time, 1000x1 / 2500x4; the old row gives µs per pair
and thread; team cells give the number of rounds. Forms without team timings were already outside the rule on one
thread (c4, c5), share their euclid loop with a timed form (c6 with c12 and c13), equal the reference (c1), or are
diagnostics (d1); c7 is the old-loop version of the two-shell form.

euclid:

| form | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| old (f3ccb5fbd) | 14.39 / 14.42 | 23.64 / 23.65 | 24.10 / 30.34 (6) | 31.02 / 38.00 (6) |
| c1 old loops, new conversion | 1.028 / 1.029 | 1.008 / 1.008 | - / - | - / - |
| c2 HEAD helpers | 0.924 / 0.927 | 0.811 / 0.809 | 0.833 / 0.834 (2) | 1.019 / 0.829 (2) |
| c3 old loops, contiguous pointer alias | 1.023 / 1.025 | 1.002 / 1.001 | 0.964 / 0.969 (6) | 1.022 / 1.007 (6) |
| c4 helpers, contiguous assumed shape | 0.927 / 0.927 | 0.802 / 0.806 | - / - | - / - |
| c5 helpers, local work array | 0.929 / 0.929 | 0.814 / 0.801 | - / - | - / - |
| c6 helpers, two shells per pass | 1.155 / 1.150 | 1.092 / 1.096 | - / - | - / - |
| c7 old loops, two shells per pass | 1.090 / 1.077 | 1.063 / 1.055 | - / - | - / - |
| c9 helper call per shell | 1.123 / 1.127 | 1.119 / 1.120 | 0.997 / 1.042 (6) | 1.165 / 1.079 (6) |
| c12 euclid two shells, cc per shell | 1.154 / 1.151 | 1.103 / 1.101 | 0.892 / 0.959 (6) | 0.896 / 1.039 (6) |
| c13 euclid two shells, cc two shells per call | 1.156 / 1.153 | 1.039 / 1.033 | 0.897 / 0.924 (6) | 0.882 / 0.989 (6) |
| c14 euclid per shell, cc two shells per call | 1.148 / 1.150 | 1.127 / 1.126 | 0.999 / 1.062 (4) | 1.284 / 1.094 (4) |
| d1 HEAD sums, old conversion (diagnostic) | 0.911 / 0.910 | 0.803 / 0.803 | - / - | - / - |

cc:

| form | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| old (f3ccb5fbd) | 17.61 / 17.68 | 20.39 / 20.50 | 25.51 / 34.25 (6) | 28.37 / 39.10 (6) |
| c1 old loops, new conversion | 0.998 / 1.001 | 1.002 / 1.001 | - / - | - / - |
| c2 HEAD helpers | 1.176 / 1.176 | 0.919 / 0.922 | 0.989 / 1.119 (2) | 0.984 / 1.110 (2) |
| c3 old loops, contiguous pointer alias | 1.182 / 1.182 | 1.282 / 1.281 | 0.947 / 1.135 (6) | 1.181 / 0.996 (6) |
| c4 helpers, contiguous assumed shape | 1.149 / 1.172 | 0.916 / 0.918 | - / - | - / - |
| c5 helpers, local work array | 1.183 / 1.183 | 0.908 / 0.907 | - / - | - / - |
| c6 helpers, two shells per pass | 1.410 / 1.427 | 0.997 / 0.996 | - / - | - / - |
| c7 old loops, two shells per pass | 1.055 / 1.059 | 1.013 / 1.010 | - / - | - / - |
| c9 helper call per shell | 1.350 / 1.341 | 1.383 / 1.386 | 1.110 / 1.172 (6) | 1.141 / 1.202 (6) |
| c12 euclid two shells, cc per shell | 1.347 / 1.339 | 1.370 / 1.380 | 1.263 / 1.146 (6) | 1.076 / 1.209 (6) |
| c13 euclid two shells, cc two shells per call | 1.386 / 1.363 | 1.297 / 1.284 | 1.472 / 1.254 (6) | 1.258 / 1.304 (6) |
| c14 euclid per shell, cc two shells per call | 1.387 / 1.398 | 1.305 / 1.307 | 1.491 / 1.218 (4) | 1.345 / 1.327 (4) |
| d1 HEAD sums, old conversion (diagnostic) | 1.186 / 1.187 | 0.919 / 0.922 | - / - | - / - |

Compiled but not timed (their code is HEAD's or c6's, checked in the assembly): c8 (c6 with cc as explicit bin
loops: cc still jammed), c10 (`!$omp simd`), c11 (`!GCC$ unroll 1`).

The tightest cell for the kept form, euclid 1000x1 at 12 threads, was timed eight more times interleaved with old
and c9: c14 1.040x (24 threads 1.052x); pooled with the six rounds 1.043x (1.057x). The cc gain of all forms that
pass the memoized FT(REF) directly (c2, c3 and every helper form) partly comes from no longer copying it into
`cmat2_many` when there is no shift.

**The rule.** Acceptable (within 2 % of the old form in all four settings, both objectives, and within two ulp):
c9 and c14. c3 misses (euclid 0.964 and cc 0.947 at 12 threads), c12 and c13 miss (two-shell euclid in teams),
HEAD misses. Chosen per objective: euclid, one helper call per shell (c9 and c14 share it); cc, two shells per
helper call (c14's, faster than c9's in every setting). Both tell the compiler that the accumulator and the memo
slices are contiguous and do not overlap. Together that is c14.

## Identity

The benchmark's `suite=identity`: 32 particles, the production kernels' per-rotation vectors at one thread,
compared bit for bit with the old form's. For cc the judged distances are those the table scores, 1 - cc (the Mac
run's reading); the raw correlations cancel to about 1e-2 and are reported in ulp of each vector's largest |cc|.
2500x4 (320,000 pairs); 1000x1 gave the same bounds or tighter and no reordered list.

| form | euclid max ulp (dists / exp values) | euclid vectors identical | 20-lowest lists reordered | cc 1-cc max ulp | raw cc, ulp of vector max | cc vectors identical | 20-lowest lists reordered |
|---|---|---|---|---|---|---|---|
| c1 | 0 / 0 | 100 % | 0 | 0 | 0 | 100 % | 0 |
| c2 HEAD, c4, c5, d1 | 2 / 2 | 74.4 % | 1 | 1 | 4 | 77.3 % | 4 |
| c3 | 0 / 0 | 100 % | 0 | 1 (shifted only) | 6.5 | 100 % unshifted | 0 |
| c6 | 2 / 2 | 77.5 % | 4 | 1 | 4 | 77.5 % | 4 |
| c7 | 1 / 2 | 67.8 % | 6 | 1 | 2 | 100 % | 0 |
| c9, c14 | 2 / 2 | 67.8 % | 6 | 1 | 2 | 99.9 % | 0 |
| c12, c13 | 2 / 2 | 77.5 % | 4 | 1 | 2 | 99.9 % | 0 |

Two ulp occur only in the shifted-reference branch and in the bank's exp(-loss) values. The best value of a pair
never changes its rotation in any form. c14's euclid bits equal c9's and its cc bits c13's.

## Checks of the final tree (Phase 4)

On a Release build of the repository with the c14 form (`compile_clean.sh`, `-O3 -march=native`):

| check | result |
|---|---|
| fast gate | 15/15 passed |
| `simple_test_exec test=unit_pftc_align2D3D` | 9/9 groups, 718/718 assertions |
| `simple_test_exec test=unit_cart_align3D` | 5/5 groups, 801/801 assertions |
| CTest `simulated_workflow_1jxy` | passed (785 s), 17/17 metrics pass |
| benchmark `suite=check`, 24 threads | every layout's 24-thread repeats equal their own one-thread run; production at 24 threads equals one thread |
| benchmark `suite=identity` against the old form | exactly c14's bits (table above); five 24-thread repeats of all four kernels equal one thread, both objectives and cases |
| re-timing against the old form | below |

Re-timing of the repository build against the old form (one pinned round, two team rounds in opposite order,
nothing else running), speed = old time / new time, 1000x1 / 2500x4:

| objective | P1 | E1 | T12 | T24 |
|---|---|---|---|---|
| euclid | 1.150 / 1.148 | 1.126 / 1.125 | 1.129 / 1.050 | 1.245 / 1.074 |
| cc | 1.376 / 1.395 | 1.323 / 1.308 | 1.113 / 1.117 | 1.170 / 1.288 |

Evidence, all under `~/agent_runs/kernel_x86/scratch/`: `logs/build_repo_p4.log`, `p4/*.log`, `p4/check_table.md`,
`p4/metrics_1jxy.tsv`, `p4/ident_final_vs_old.txt`, `logs/retime_table.txt`; Phase 3's timings `time/<round>/`,
tables `logs/mdtable.txt`, identity `logs/ident_table_*.txt`, compiler reports and assembly `diag/`.

## Decisions taken under delegation

All recorded in the journals (`~/agent_runs/kernel_x86/state/journal_phase_3.md`, `journal_phase_4.md`) and, where
they shape the result, in the plan's Progress.

- **Waiting for an idle machine.** Another run's jobs held all 24 cores until 02:35; no timing ran before 02:47.
  Builds, compiler reports and the identity runs (bits, not time) ran before, niced.
- **The benchmark's middle team** is half the threads (test side), so that `nthr=24` times the 12-thread team the
  run's instructions ask for; the Mac's `nthr=11` still gives 11, 5, 1.
- **cc identity** judged on the distances the table scores, 1 - cc, as the Mac run did; raw correlations reported
  in ulp of each vector's largest |cc|.
- **Team timings** judged on the median of six unpinned rounds (they move 10-30 % between launches), plus eight
  interleaved launches for the closest cell. Production runs teams unpinned, so pinning them was not an option.
- **Candidates beyond the run's list.** The diagnosis pointed at the jam, so forms that avoid it were added (c6-c14);
  c8, c10 and c11 were dropped after their compiler reports showed HEAD's or c6's code.
- **Which forms got team timings:** those not already ruled out on one thread or equal to a timed form (see the
  table's note).
- **The choice** made per objective as the rule says; the two winners combine into the complete form c14.
- **A Release build** of the repository (`compile_clean.sh`, fast gate included) rather than the debug script, so
  that the re-timing compares like with like; the unit suites and the simulated workflow ran on it.
- **Disk:** 48 GB of identity vectors deleted after their comparisons were written to text; the old form's
  vectors kept under `scratch/keep/ident_old`.

## Deviations from the plan

- The plan (the Mac run's) had no Phases 3-4 and no file table; both were added from the run's instructions, with
  the judging rule fixed before any measurement.
- The benchmark driver change (test side) is the only change outside the plan's production file and documents.
- The re-timing of the chosen form is one pinned round and two team rounds against the old form, not one run, for
  the same noise reason as above.
- `simulated_workflow_1jxy` was checked for passing (17/17 metrics), not compared byte for byte with a run of the
  old build; the identity suite covers the kernels themselves.

## Open for the maintainer

- **Real runs.** All numbers are from the synthetic benchmark. What share of a real dense table's time the kernels
  take on the Dell was measured by `ms_prob5check` for HEAD (score loop 52-68 % of the fill time); a real run with
  c14 was not made.
- **The memory-stream explanation** of the two-shell euclid's team loss is a hypothesis: `perf` is not installed on
  the Dell, so there are no counter measurements.
- **Not touched:** `gen_raw_euclid_vals_impl` still accumulates through the pointer component (as the Mac report
  already noted); it is not on the dense table's path.

## Commit message

```text
pftc: per-shell helper calls for the euclid and cc sums

With the shell loop and the bin loop in one helper, gcc unroll-and-jams two
shells; on x86 (AVX2) the jammed loop needs lane permutes and register spills,
and the euclid kernel ran at 0.81-0.92x of the old inline loops. The euclid
sums now call a helper once per shell and the cc sums once per two shells, still
with plain array arguments (unit stride, no overlap). On the Dell euclid is
1.13-1.15x and cc 1.31-1.39x the old form on one core, and at least as fast in
teams; distances change by at most 2 ulp. The benchmark's middle team is half
the threads.
```

## Appendix: the source of every form

Every form is HEAD's `src/main/pftc/simple_polarft_corr.f90` (ee5b68be1) with the routines listed under it in place of
HEAD's routines of the same name; a routine shown as "deleted" is removed. The other routines of the file, the
distance conversion included, are HEAD's. The diagnostic d1 is the opposite: f3ccb5fbd's file with HEAD's six
routines of c2 (and nothing else of HEAD). A routine whose source equals one printed earlier is referred to.
The repository now holds c14 with two lines added to the doc comments of `euclid_sums` and `cc_sums` (why the
shell loop is split into helper calls); the code is the same.

### c2: HEAD (ee5b68be1): explicit-shape helpers, called with the memo sections

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(X.CTF) x FT(S.REF)* ), the memoized FT(REF) when there is no shift
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call cc_sums(self%pftsz+1, self%nrots, self%nk, self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, .true.,&
                &self%crvec1(ithr)%c)
        else
            call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), .true.,&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call euclid_sums(self%pftsz, self%pftsz+1, self%nrots, self%nk, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        else
            call euclid_sums(self%pftsz, self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        integer :: i, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        call euclid_sums(self%pftsz+1, self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i),&
            &self%ft_ref2(:,:,iref,ieo), bank, self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih), self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, bank, self%ft_ref(:,:,iref,ieo), .true., self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk
                acc = acc + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            do kk = 1,nk
                acc = acc + a(:,kk) * b(1:n1,kk)
            enddo
        endif
    end subroutine cc_sums
```

### c1: the old inline loops of f3ccb5fbd (with HEAD's distance conversion); helpers deleted

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        logical  :: needs_shift
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! ========================================================================
        ! Prepare FFT(S.REF) if needed
        ! ========================================================================
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        needs_shift  = shift_mag_sq > SHERRSQ
        if (needs_shift) then
            ! Generate shift matrix
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Shift reference directly into cmat2_many
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            ! Calculate FFT(S.REF)
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        else
            ! no shift, memoized reference is used
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = &
                self%ft_ref(:,self%kfromto(1):self%kfromto(2),iref,ieo)
        endif
        ! ========================================================================
        ! Single IFFT #1: IFFT( sum_k FT(CTF2) x FT(REF2) )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! ========================================================================
        ! Single IFFT #2: IFFT( sum_k FT(X.CTF) x FT(S.REF)* )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + &
                self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,kk))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        complex(sp) :: c
        real(dp)    :: wk
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
        else
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            do k = self%kfromto(1), self%kfromto(2)
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

- `euclid_sums`: deleted

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        complex(sp) :: c
        real(dp)    :: wk
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            wk = self%sigma2_wk8(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                    &kind=kind(self%crvec1(ithr)%c))
            enddo
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, k, kk, k0, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + bank(:,kk) * conjg(self%ft_ref(:,k,iref,ieo))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

- `cc_sums`: deleted

### c3: old inline loops with a local `pointer, contiguous` alias of the FFT buffers; FT(REF) read directly without a shift

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer             :: shmat(:,:)
        complex(sp), pointer, contiguous :: acc(:), sref(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! contiguous aliases of the per-thread FFT buffers: unit stride without a call
        acc     => self%crvec1(ithr)%c
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        acc = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            acc = acc + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(X.CTF) x FT(S.REF)* ), the memoized FT(REF) when there is no shift
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        acc = cmplx(0.,0.,kind=sp)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            sref => self%cmat2_many(ithr)%c
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                acc = acc + self%ft_ptcl_ctf(:,k,i) * conjg(sref(1:self%pftsz+1,kk))
            enddo
        else
            do k = self%kfromto(1), self%kfromto(2)
                acc = acc + self%ft_ptcl_ctf(:,k,i) * conjg(self%ft_ref(:,k,iref,ieo))
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer             :: shmat(:,:)
        complex(sp), pointer, contiguous :: acc(:), sref(:,:)
        complex(sp) :: c
        real(dp)    :: wk
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        ! contiguous aliases of the per-thread FFT buffers: unit stride without a call
        acc          => self%crvec1(ithr)%c
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            sref => self%cmat2_many(ithr)%c
            acc  = cmplx(0.,0.,kind=sp)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(sref(p,kk))
                    acc(p) = cmplx(acc(p) + wk * c, kind=sp)
                enddo
            enddo
        else
            acc = cmplx(0.,0.,kind=sp)
            do k = self%kfromto(1), self%kfromto(2)
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    acc(p) = cmplx(acc(p) + wk * c, kind=sp)
                enddo
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

- `euclid_sums`: deleted

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        complex(sp), pointer, contiguous :: acc(:)
        complex(sp) :: c
        real(dp)    :: wk
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        acc  => self%crvec1(ithr)%c
        acc  = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            wk = self%sigma2_wk8(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                acc(p) = cmplx(acc(p) + wk * c, kind=sp)
            enddo
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        complex(sp), pointer, contiguous :: acc(:)
        integer :: i, ithr, k, kk, k0, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        acc  => self%crvec1(ithr)%c
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        acc = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            acc = acc + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        acc = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            acc = acc + bank(:,kk) * conjg(self%ft_ref(:,k,iref,ieo))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

- `cc_sums`: deleted

### c4: HEAD's helpers with assumed-shape `contiguous` dummies

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(X.CTF) x FT(S.REF)* ), the memoized FT(REF) when there is no shift
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call cc_sums(self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, .true.,&
                &self%crvec1(ithr)%c)
        else
            call cc_sums(self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), .true.,&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call euclid_sums(self%pftsz, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        else
            call euclid_sums(self%pftsz, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Contiguous dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, ctf2, ref2, xctf, ref, wk, acc )
        integer,                  intent(in)    :: np
        complex(sp), contiguous,  intent(in)    :: ctf2(:,:), ref2(:,:), xctf(:,:), ref(:,:)
        real(dp),    contiguous,  intent(in)    :: wk(:)
        complex(sp), contiguous,  intent(inout) :: acc(:)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,size(ctf2,2)
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        integer :: i, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        call euclid_sums(self%pftsz+1, self%ft_ctf2(:,:,i),&
            &self%ft_ref2(:,:,iref,ieo), bank, self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih), self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        call cc_sums(bank, self%ft_ref(:,:,iref,ieo), .true., self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the angular bins of acc, contiguous dummies as in euclid_sums
    subroutine cc_sums( a, b, l_conj, acc )
        complex(sp), contiguous, intent(in)    :: a(:,:), b(:,:)
        logical,                 intent(in)    :: l_conj
        complex(sp), contiguous, intent(inout) :: acc(:)
        integer :: kk, n1
        n1  = size(acc)
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,size(a,2)
                acc = acc + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            do kk = 1,size(a,2)
                acc = acc + a(:,kk) * b(1:n1,kk)
            enddo
        endif
    end subroutine cc_sums
```

### c5: HEAD's helpers accumulating into a local work array, copied once

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Accumulated in a local array, copied to acc once.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c, w(n1)
        integer     :: kk, p
        w = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                w(p) = cmplx(w(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
        acc = w
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, accumulated in a local array
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: w(n1)
        integer :: kk
        w = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk
                w = w + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            do kk = 1,nk
                w = w + a(:,kk) * b(1:n1,kk)
            enddo
        endif
        acc = w
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c6: HEAD's helpers, two shells per pass written out

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    !! Two shells per pass over acc, each still rounded to single precision in shell order.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c1, c2
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk-1,2
            do p = 1,np
                c1 = ctf2(p,kk) * ref2(p,kk)
                c1 = c1 - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                c2 = ctf2(p,kk+1) * ref2(p,kk+1)
                c2 = c2 - 2.0 * xctf(p,kk+1) * conjg(ref(p,kk+1))
                acc(p) = cmplx(cmplx(acc(p) + wk(kk) * c1, kind=sp) + wk(kk+1) * c2, kind=sp)
            enddo
        enddo
        if( mod(nk,2) == 1 )then
            do p = 1,np
                c1 = ctf2(p,nk) * ref2(p,nk)
                c1 = c1 - 2.0 * xctf(p,nk) * conjg(ref(p,nk))
                acc(p) = cmplx(acc(p) + wk(nk) * c1, kind=sp)
            enddo
        endif
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums,
    !! two shells per pass, summed in shell order
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk-1,2
                acc = (acc + a(:,kk) * conjg(b(1:n1,kk))) + a(:,kk+1) * conjg(b(1:n1,kk+1))
            enddo
            if( mod(nk,2) == 1 ) acc = acc + a(:,nk) * conjg(b(1:n1,nk))
        else
            do kk = 1,nk-1,2
                acc = (acc + a(:,kk) * b(1:n1,kk)) + a(:,kk+1) * b(1:n1,kk+1)
            enddo
            if( mod(nk,2) == 1 ) acc = acc + a(:,nk) * b(1:n1,nk)
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c7: old inline loops, two shells per pass written out

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        logical  :: needs_shift
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! ========================================================================
        ! Prepare FFT(S.REF) if needed
        ! ========================================================================
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        needs_shift  = shift_mag_sq > SHERRSQ
        if (needs_shift) then
            ! Generate shift matrix
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Shift reference directly into cmat2_many
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            ! Calculate FFT(S.REF)
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        else
            ! no shift, memoized reference is used
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = &
                self%ft_ref(:,self%kfromto(1):self%kfromto(2),iref,ieo)
        endif
        ! ========================================================================
        ! Single IFFT #1: IFFT( sum_k FT(CTF2) x FT(REF2) )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)) +&
                &self%ft_ctf2(:,k+1,i) * self%ft_ref2(:,k+1,iref,ieo)
        enddo
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c +&
            &self%ft_ctf2(:,self%kfromto(2),i) * self%ft_ref2(:,self%kfromto(2),iref,ieo)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! ========================================================================
        ! Single IFFT #2: IFFT( sum_k FT(X.CTF) x FT(S.REF)* )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            kk = k - k0 + 1
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + &
                self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,kk))) +&
                &self%ft_ptcl_ctf(:,k+1,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,kk+1))
        end do
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c + &
            &self%ft_ptcl_ctf(:,self%kfromto(2),i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,self%nk))
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        complex(sp) :: c, c2
        real(dp)    :: wk, wk2
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            ! two shells per pass, each rounded to single precision in shell order
            do k = self%kfromto(1), self%kfromto(2)-1, 2
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                wk2 = self%sigma2_wk8(k+1,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                    c2 = self%ft_ctf2(p,k+1,i) * self%ft_ref2(p,k+1,iref,ieo)
                    c2 = c2 - 2.0 * self%ft_ptcl_ctf(p,k+1,i) * conjg(self%cmat2_many(ithr)%c(p,kk+1))
                    self%crvec1(ithr)%c(p) = cmplx(cmplx(self%crvec1(ithr)%c(p) + wk * c, kind=sp) + wk2 * c2, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
            if( mod(self%nk,2) == 1 )then
                k  = self%kfromto(2)
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            endif
        else
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            ! two shells per pass, each rounded to single precision in shell order
            do k = self%kfromto(1), self%kfromto(2)-1, 2
                wk  = self%sigma2_wk8(k,ih)
                wk2 = self%sigma2_wk8(k+1,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    c2 = self%ft_ctf2(p,k+1,i) * self%ft_ref2(p,k+1,iref,ieo)
                    c2 = c2 - 2.0 * self%ft_ptcl_ctf(p,k+1,i) * conjg(self%ft_ref(p,k+1,iref,ieo))
                    self%crvec1(ithr)%c(p) = cmplx(cmplx(self%crvec1(ithr)%c(p) + wk * c, kind=sp) + wk2 * c2, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
            if( mod(self%nk,2) == 1 )then
                k  = self%kfromto(2)
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            endif
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

- `euclid_sums`: deleted

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        complex(sp) :: c, c2
        real(dp)    :: wk, wk2
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        ! two shells per pass, each rounded to single precision in shell order
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            kk  = k - k0 + 1
            wk  = self%sigma2_wk8(k,ih)
            wk2 = self%sigma2_wk8(k+1,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                c2 = self%ft_ctf2(p,k+1,i) * self%ft_ref2(p,k+1,iref,ieo)
                c2 = c2 - 2.0 * bank(p,kk+1) * conjg(self%ft_ref(p,k+1,iref,ieo))
                self%crvec1(ithr)%c(p) = cmplx(cmplx(self%crvec1(ithr)%c(p) + wk * c, kind=sp) + wk2 * c2, &
                    &kind=kind(self%crvec1(ithr)%c))
            enddo
        enddo
        if( mod(self%nk,2) == 1 )then
            k  = self%kfromto(2)
            kk = k - k0 + 1
            wk = self%sigma2_wk8(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                    &kind=kind(self%crvec1(ithr)%c))
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, k, kk, k0, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)) +&
                &self%ft_ctf2(:,k+1,i) * self%ft_ref2(:,k+1,iref,ieo)
        enddo
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c +&
            &self%ft_ctf2(:,self%kfromto(2),i) * self%ft_ref2(:,self%kfromto(2),iref,ieo)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            kk = k - k0 + 1
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + bank(:,kk) * conjg(self%ft_ref(:,k,iref,ieo))) +&
                &bank(:,kk+1) * conjg(self%ft_ref(:,k+1,iref,ieo))
        end do
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c +&
            &bank(:,self%nk) * conjg(self%ft_ref(:,self%kfromto(2),iref,ieo))
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

- `cc_sums`: deleted

### c8: c6 with the cc sums as explicit bin loops (compiled only: cc still jammed)

- `euclid_sums`: as c6

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums,
    !! two shells per pass, summed in shell order
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk-1,2
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,kk) * conjg(b(p,kk))) + a(p,kk+1) * conjg(b(p,kk+1))
                enddo
            enddo
            if( mod(nk,2) == 1 )then
                do p = 1,n1
                    acc(p) = acc(p) + a(p,nk) * conjg(b(p,nk))
                enddo
            endif
        else
            do kk = 1,nk-1,2
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,kk) * b(p,kk)) + a(p,kk+1) * b(p,kk+1)
                enddo
            enddo
            if( mod(nk,2) == 1 )then
                do p = 1,n1
                    acc(p) = acc(p) + a(p,nk) * b(p,nk)
                enddo
            endif
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c9: one helper call per shell, euclid and cc

`euclid_sums` (with the helpers it calls):

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc, one shell at a time (euclid_shell_sum)
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            call euclid_shell_sum(np, ctf2(:,kk), ref2(:,kk), xctf(:,kk), ref(:,kk), wk(kk), acc)
        enddo
    end subroutine euclid_sums

    !> One shell of euclid_sums: plain array dummies, unit stride and no overlap
    subroutine euclid_shell_sum( np, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np
        complex(sp), intent(in)    :: ctf2(np), ref2(np), xctf(np), ref(np)
        real(dp),    intent(in)    :: wk
        complex(sp), intent(inout) :: acc(np)
        complex(sp) :: c
        integer     :: p
        do p = 1,np
            c = ctf2(p) * ref2(p)
            c = c - 2.0 * xctf(p) * conjg(ref(p))
            acc(p) = cmplx(acc(p) + wk * c, kind=sp)
        enddo
    end subroutine euclid_shell_sum
```

`cc_sums` (with the helpers it calls):

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, one shell at a time (cc_shell_sum)
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            call cc_shell_sum(n1, a(:,kk), b(:,kk), l_conj, acc)
        enddo
    end subroutine cc_sums

    !> One shell of cc_sums: plain array dummies, unit stride and no overlap
    subroutine cc_shell_sum( n1, a, b, l_conj, acc )
        integer,     intent(in)    :: n1
        complex(sp), intent(in)    :: a(n1), b(n1)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        if( l_conj )then
            acc = acc + a * conjg(b)
        else
            acc = acc + a * b
        endif
    end subroutine cc_shell_sum
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c10: HEAD with `!$omp simd` on the bin loops (compiled only: jam not prevented)

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            !$omp simd private(c)
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk
                !$omp simd
                do p = 1,n1
                    acc(p) = acc(p) + a(p,kk) * conjg(b(p,kk))
                enddo
            enddo
        else
            do kk = 1,nk
                !$omp simd
                do p = 1,n1
                    acc(p) = acc(p) + a(p,kk) * b(p,kk)
                enddo
            enddo
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c11: HEAD with `!GCC$ unroll 1` on the shell loops (compiled only: jam not prevented)

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        !GCC$ unroll 1
        do kk = 1,nk
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            !GCC$ unroll 1
            do kk = 1,nk
                acc = acc + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            !GCC$ unroll 1
            do kk = 1,nk
                acc = acc + a(:,kk) * b(1:n1,kk)
            enddo
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c12: euclid as c6, cc as c9

- `euclid_sums`: as c6

- `cc_sums`: as c9

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c13: euclid as c6, cc two shells per helper call

- `euclid_sums`: as c6

`cc_sums` (with the helpers it calls):

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, two shells per call of cc_shells_sum
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk-1,2
            call cc_shells_sum(n1, ldb, 2, a(:,kk:kk+1), b(:,kk:kk+1), l_conj, acc)
        enddo
        if( mod(nk,2) == 1 ) call cc_shells_sum(n1, ldb, 1, a(:,nk:nk), b(:,nk:nk), l_conj, acc)
    end subroutine cc_sums

    !> One or two shells of cc_sums in one pass over acc, summed in shell order; plain array dummies
    subroutine cc_shells_sum( n1, ldb, ns, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, ns
        complex(sp), intent(in)    :: a(n1,ns), b(ldb,ns)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: p
        if( ns == 2 )then
            if( l_conj )then
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,1) * conjg(b(p,1))) + a(p,2) * conjg(b(p,2))
                enddo
            else
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,1) * b(p,1)) + a(p,2) * b(p,2)
                enddo
            endif
        else
            if( l_conj )then
                acc = acc + a(:,1) * conjg(b(1:n1,1))
            else
                acc = acc + a(:,1) * b(1:n1,1)
            endif
        endif
    end subroutine cc_shells_sum
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c14: euclid as c9, cc as c13: the form kept

- `euclid_sums`: as c9

- `cc_sums`: as c13

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)
