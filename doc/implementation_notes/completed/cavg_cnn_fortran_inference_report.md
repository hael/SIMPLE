# The class-average CNN in Fortran: report

Run of [the plan](cavg_cnn_fortran_inference.md), 2026-11-20, interactive, on the development machine (Ubuntu,
gfortran 15.2, 16 cores) with the Debug build on a second machine (GCC 15.1 toolset, 72 cores). Base `a26f57437`
with the uncommitted `scripts/sscnn/`. Everything is uncommitted for review.

## 1. What was learned

- **The forward pass reproduces PyTorch.** On the fixture the logits differ by at most 5.1e-7 (T1 3.4e-6), with
  `sgemm` and with the explicit loops alike. On three real stacks (205 classes) the probabilities differ by at most
  2.0e-7 (T2 7.5e-6); no group differs and the empty classes agree.
- **There is no transposition.** Written in PyTorch's memory order and read column-major, a weights tensor or an
  image keeps its layout: a Fortran (w, h, c) array is PyTorch's (c, h, w), and SIMPLE's `rmat(x, y)` is (w, h).
  The plan expected to rely on the 8-symmetry average to cancel a transposition; the plain logits match PyTorch's,
  so the parity tests check the orientation directly.
- **The Python resampling had two flaws** (decision 12), found while checking the premises and fixed in
  `scripts/sscnn/cnn_data.py`: the kept frequencies were one bin off for odd resampled sizes, and the transforms'
  origin was the corner pixel, so the centre drifted when the size ratio is not exact. The corrected function is
  exact to 1e-15 against band-limited images at every size. Only stacks resampled to odd sizes or inexact ratios
  change: the transporter's two classes moved in the third decimal (still secondary structure at 0.99).
- **The plan's T2 rule was incomplete.** It counted only the input difference (5.6e-9) and would have failed on the
  network's own rounding. Decision 15 adds half of T1 measured on the real model.
- **Speed is not an issue.** 100 class averages of 320 px take 3.1 s at 16 threads and 31 s on one (the bound was
  60 s at 16). The input preparation's exact operator costs little beside the network.
- **The repository ignores `*.bin` and `*.dat`.** The fixture would have been left out of a commit; it is `*.raw`.

## 2. What changed

| Phase | Files | What |
| --- | --- | --- |
| 0 | 2 | `scripts/sscnn/export_fortran.py` (new); `scripts/sscnn/cnn_data.py` (resampling) |
| 0, 1 | 3 | `testdata/cavg_cnn/` (new fixture: `weights.raw`, `inputs.raw`, `logits.txt`, 0.9 MB) |
| 1 | 4 | `src/main/learn/simple_cavg_cnn.f90` and its tester (new); the suite table; the test UI's `suite=` list |
| 2 | 2 | `src/main/cavg_quality/simple_cavg_cnn_input.f90` and its tester (new) |
| 3 | 5 | `src/defs/simple_defs_cavg_cnn.f90` and `src/main/cavg_quality/simple_cavg_cnn_predict.f90` (new); the commander in `simple_commanders_cavgs.f90`; the route in `simple_exec_cavgproc.f90`; the descriptor in `simple_ui_cavgproc.f90` |
| 4 | 6 | `doc/how2s/how2train_cavg_cnn.md`; `scripts/sscnn/README.md`, `CNN.md`; `library_owners.md`; the class-average quality guide; the notes index |

## 3. Decisions taken under delegation

Decisions 12 to 17 of the plan, each with its reason there:

- 12, the resampling: an exact operator rather than the image class, whose `fft` takes even boxes only and whose
  `pad` weights the spectrum; the Python fixed to match.
- 13, the halves' names: `add2fbody`, as `parameters` names `refs_even`.
- 14, the fixture's form: small inputs, `.raw` files.
- 15, T2: the network's rounding added.
- 16, the program's body: `cavg_cnn_predict_exec(params)` in `cavg_quality/`.
- 17, the orientation: no transposition.

## 4. Deviations

- **Debug build on another machine.** The development machine's disk (1.1 GB free) cannot hold a second build. A
  copy of the working tree was built with `compile_debug.sh` on the second machine (GCC 15.1 toolset), its fast
  gate passed, and after the fixture's rename both suites and the program on two stacks ran there under bounds
  checks.
- **The third negative control** became the average cut to the four rotations (decision 17).
- **Input parity was measured by a temporary harness** in the input tester (removed after its run): the unit
  tests cannot read user data.
- **The second reference stack has no project.** One was made with `import_cavgs` in a scratch folder.
- **Files outside the file table:** `scripts/sscnn/cnn_data.py` (decision 12), `src/defs/simple_defs_cavg_cnn.f90`
  (workflow constants belong in a defs module), `src/main/cavg_quality/simple_cavg_cnn_input_tester.f90`,
  `src/main/ui/simple_test/simple_test_ui_class.f90` (`suite=` help), `scripts/sscnn/CNN.md`, the class-average
  quality guide and the notes index.
- **Not needed:** `learn.inf` and the generated UI files were not edited.

## 5. Evidence

| Gate | Result |
| --- | --- |
| T1 (measured, fixture) | 3.361e-6 |
| fixture parity, sgemm / loops | plain 4.8e-7 / 5.1e-7; averaged 1.2e-7 / 8.2e-8 |
| negative controls | row-major kernel 0.27, stride offset 0.19, four rotations 0.04: all fail the suite |
| unit suites (Release; Debug) | `cavg_cnn` 41 checks, `cavg_cnn_input` 12; pass in both |
| input parity, three stacks | correlation 1.0000000, largest difference 6e-8, empty classes equal |
| T2 (measured, frozen before Phase 3) | 7.466e-6 |
| end to end, three stacks | largest dP 7.0e-8, 2.0e-7, 2.1e-9; 0 group differences in 187 decided classes |
| thread identity | the table byte-identical at 1 and 16 threads, two runs each |
| read-only | the projects' checksums unchanged |
| speed, 100 classes of 320 px | 3.1 s at 16 threads, 31 s on one |
| fast gate, `unit_ui` | 15 of 15 (Release and Debug); 1,042 checks |

## 6. Left open for the maintainer

- **Reported, not gated.** The three prediction tables, in the development machine's
  `~/Data/ss_sim/fortran_ref/tables/` (`cavg_cnn_predictions_{motab,porf,slc}.txt`), with the reference
  probabilities and input dumps beside them. They are the Python predictions, unchanged; the read is whether the
  program is ready to use.
- **Using the predictions** (in `model_cavgs_rejection`, the sieve or the stream) is a later plan, after a
  real-data read of the secondary-structure output, which is still open (`scripts/sscnn/CNN.md`).
- **A trained model file** stays outside the repository (decision 11); users export their own.

## 7. Commit message

```
cavg_cnn_predict: the class-average CNN of scripts/sscnn in Fortran

The network (simple_cavg_cnn: weights file, im2col convolutions with the linked
sgemm or loops, GroupNorm, the 8-symmetry average), its input (simple_cavg_cnn_input:
numpy's Fourier resampling as one exact operator, any size, the centre kept) and the
read-only program cavg_cnn_predict, which writes each class's probabilities and group
and changes nothing in the project. Same probabilities as PyTorch to 2e-7 on three
real stacks; byte-identical at any thread count; 3 s per 100 classes at 16 threads.
scripts/sscnn: export_fortran.py; the resampling's centring fixed for odd sizes.
```
