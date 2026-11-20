# The class-average CNN in Fortran: inference with no new libraries

Status, 2026-11-20: delivered the same day, uncommitted for review (report: `cavg_cnn_fortran_inference_report.md`);
decisions 3, 5, 9 and 11 ruled by the user, 12 to 17 taken under delegation. Base: `a26f57437` plus the uncommitted `scripts/sscnn/` and
`doc/how2s/how2train_cavg_cnn.md` of the same date. Supersedes: none. Depends on: the trained class-average CNN
of `scripts/sscnn/` (validated: on simulated data, yes, AUC 0.995 protein against junk and 0.961 secondary
structure against likely protein on held-out proteins; on real data, protein against junk by eye on two sets,
secondary structure no; `scripts/sscnn/CNN.md`, Results); the image class's Fourier resampling (validated: yes).

## 1. Why

The class-average CNN sorts class averages into junk, secondary structure present and likely protein. Today it
runs only through Python and PyTorch (`scripts/sscnn/cnn_ss.py predict`), which SIMPLE's programs cannot call and
which the stream and the sieve cannot depend on. The network is small: 2.77 M weights, one 3x3 convolution at full
resolution, five residual blocks of 3x3 convolutions with stride 2, GroupNorm, ReLU, a global mean and one linear
layer. Its forward pass is about 1 G multiply-adds per class average (computed from the layer sizes, not measured).
This plan gives SIMPLE that forward pass in Fortran with the libraries the build already links, proves it equal
to PyTorch's, and exposes it as a read-only program. It changes no rejection decision.

## 2. Contract

- **Inference only.** No training, no gradients, no optimiser in Fortran. Models are trained by `scripts/sscnn/`.
- **No new libraries** (ruling, 2026-11-20: existing ones are fine). The new modules use Fortran intrinsics,
  OpenMP, SIMPLE's own modules and what the build already links: the BLAS (`sgemm`, as `simple_hhcnn` uses it) and
  FFTW through the image class. No new dependency, no C, no Python at run time.
- **Read-only.** The new program reads a project and a weights file and writes a table. It never writes the
  project, a state, a `cls2D` field or a stack.
- **Off state.** Nothing existing calls the new modules. `model_cavgs_rejection`, `src/main/cavg_quality/`,
  `src/main/quality/`, the sieve, the stream, solve2D and `simple_hhcnn` are not edited, so their outputs are
  untouched by construction.
- **The average over the eight symmetries of the square is part of the model.** PyTorch reads an MRC image as
  rows by columns and SIMPLE holds it as x by y, which is a transposition. Averaging the logits over the four
  rotations and their mirrors makes the result independent of that. The Fortran forward pass always averages the
  eight; there is no switch. (Decision 17: there turned out to be no transposition; the average stays.)
- **Thread-count independence.** One class average is one single-threaded forward pass; classes are distributed
  over threads, and `sgemm` runs single-threaded inside the team. The table is byte-identical at any thread count.
- **One reading of the weights.** The Fortran reader stops on a file whose magic, version, architecture record or
  tensor sizes differ from what it expects. No silent default for a missing tensor.

## 3. Decisions and rulings

1. (2026-11-20, proposed) **Weights file.** A stream-access unformatted file, little-endian, written by a new
   `scripts/sscnn/export_fortran.py` from `model.pt`: an `int64` magic and version, then `int32` width, input
   channels (1 or 3), box (160), GroupNorm groups (8), number of blocks (5), then `real32` pixel size (2.0), then
   the tensors in forward order, each preceded by its `int32` rank and extents, as `real32` in Fortran (column-major)
   order. The exporter does the reordering; the reader does none.
2. (2026-11-20, proposed) **Where the code lives.** The network is `src/main/learn/simple_cavg_cnn.f90` (type
   `cavg_cnn` with `new`, `read_weights`, `predict`, `kill`; submodule or sibling `simple_cavg_cnn_tester.f90`).
   Input preparation (resample, normalise, the three channels) is `src/main/cavg_quality/simple_cavg_cnn_input.f90`,
   beside the other class-average feature code. Orchestration is a strategy-style routine `run_cavg_cnn_predict(cline)`
   with a thin commander.
3. (2026-11-20, ruled by the user) **Program.** `simple_exec prg=cavg_cnn_predict`, routed by `simple_exec_cavgproc` beside
   `model_cavgs_rejection`. Keys, all existing: `projfile`, `infile` (the weights file), `nthr`, `mkdir`. No new
   command-line parameter.
4. (2026-11-20, proposed) **Output.** `cavg_cnn_predictions.txt` in the execution folder, one line per class in
   SIMPLE's `key=value` form: `class`, `pop`, `res`, `p_junk`, `p_ss`, `p_likely`, `group` (`junk`,
   `secondary_structure`, `likely_protein`, `empty`). An all-zero class average is `empty` and is not passed
   through the network.
5. (2026-11-20, ruled by the user, after the ruling on libraries) **Convolution.** im2col (zero
   padding 1, stride 1 or 2) and one `sgemm` per layer, `real32` throughout, under the same condition
   `simple_hhcnn` uses to call `sgemm` inside a thread team, with its explicit-loop product otherwise. Phase 0
   reads how `simple_hhcnn` decides that and whether its product wrappers can be shared without editing it; if
   they cannot, the new module declares the `sgemm` interface itself, as `simple_hhcnn` does. The 1x1 stride-2
   convolution of a skip is the same product without im2col.
6. (2026-11-20, proposed) **GroupNorm.** 8 groups, epsilon 1e-5 (PyTorch's default; the exporter asserts the
   model's value equals it), statistics accumulated in `real64`, biased variance.
7. (2026-11-20, proposed) **Resampling.** The image class's Fourier pad and clip to
   `nint(box * smpd / 2.0)` pixels, then a centred crop or zero pad to 160 (replaced by decision 12). Normalisation as
   `cnn_data.normalise`: all three channels have the average's mean inside the inscribed circle (radius < box/2)
   subtracted and are divided by the average's standard deviation there.
8. (2026-11-20, proposed) **Class-average files.** The stack and its halves come from the project's `out` segment
   and the helpers of `simple_refine2D_fnames`; no name is built by hand. (Replaced by decision 13.)
9. (2026-11-20, ruled by the user) **Test fixture.** A width-8 model with seeded random weights (about 173,000 values,
   0.7 MB), four 160 px input tensors and PyTorch's logits for them, exported once by
   `scripts/sscnn/export_fortran.py --fixture` and committed under `testdata/cavg_cnn/`. Unit tests never run
   Python.
10. (2026-11-20, proposed) **Speed is reported, with one bound.** A stack of 100 class averages at 16 threads
    must finish in 60 s on the development machine (guessed; Phase 0 cannot measure it, Phase 3 does). If it
    does not, the run stops and reports the profile; it adds no library.

11. (2026-11-20, ruled by the user) **The model comes from a file.** The exported weights file is passed with
    `infile`; no trained model is committed or installed, and the program stops without `infile`.
12. (2026-11-20, under delegation) **Resampling.** Not the image class (its `fft` takes even boxes only and `pad`
    weights the spectrum against aliasing) but numpy's Fourier crop or pad as one exact operator per axis,
    out = Re(R img R^T), in double precision, with the transforms' origin at the centre pixel n/2 (0-based). Its
    Python counterpart had two flaws, fixed in `scripts/sscnn/cnn_data.py`: the kept frequencies were one bin off
    when the resampled size is odd, and the transforms' origin was the corner pixel, so the centre drifted by up to
    half a pixel when the size ratio is not exact. Only stacks resampled to such sizes change (the transporter
    stack: third decimal).
13. (2026-11-20, under delegation) **Halves.** `add2fbody(stk, ext, '_even')` and `'_odd'`, as `parameters` names
    `refs_even`: the `simple_refine2D_fnames` helpers need an iteration number a project's stack does not carry.
14. (2026-11-20, under delegation) **Fixture form.** Inputs of 64, 64 and 96 pixels (the network is fully
    convolutional; 160 px is covered end to end), 0.9 MB in all; `weights.raw` and `inputs.raw`, since the
    repository ignores `*.bin` and `*.dat`.
15. (2026-11-20, under delegation; replaces section 6's T2) **T2** = 2 x (the largest |dP| from feeding the Fortran
    inputs to PyTorch) + 0.5 x T1_real, T1_real being ten times PyTorch's float32-float64 logit difference on the
    real model and the reference inputs (T1's rule on the real network; 0.5 bounds a probability's change when each
    logit moves by T1_real). As written, T2 counted the input difference alone and would have failed on the
    network's own rounding.
16. (2026-11-20, under delegation) **The program's body** is `cavg_cnn_predict_exec(params)` in
    `src/main/cavg_quality/simple_cavg_cnn_predict.f90` (the file table's second home); the commander sets
    `mkdir=yes`, calls `params%new` once and the body. No `run_<x>(cline)` entry: that form is the master steps', and
    this is a whole program (layering policy, section 3).
17. (2026-11-20, under delegation) **Orientation.** The weights are written in PyTorch's memory order and read
    column-major, so a Fortran (w, h, c) array is PyTorch's (c, h, w) and SIMPLE's `rmat(x, y)` is (w, h): there is
    no transposition. The plain logits match PyTorch's; the 8-symmetry average stays part of the model, as in
    `cnn_ss.py`. The third negative control became the average cut to the four rotations.

Defaults taken without a question: the module, type, program and file names above; the fixture's seed (1); the
table's column order.

## 4. Not in this plan

- No change to any rejection decision, preset, default, state or project field. Using the predictions in
  `model_cavgs_rejection`, the sieve or the stream is a later plan, written after the real-data read of Phase 4.
- No training in Fortran, no GPU path, no half-precision or quantised weights.
- No trained model file in the repository (11 MB); only the 0.7 MB fixture.
- No edit of `simple_hhcnn`, the students, `src/main/cavg_quality/` files other than the one new module, or
  `scripts/sscnn/` files other than the new exporter.
- A model trained with `--channels avg` is read and run (1 input channel) but its real-data read is not part of
  this plan.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | measurement | builds, fast gate, PyTorch reference tables |
| 1 | new code, tests-only reach | layer tests and fixture parity in the unit suite |
| 2 | new code, tests-only reach | input-preparation parity |
| 3 | new program, read-only | end-to-end parity; byte identity across thread counts; existing outputs untouched |
| 4 | documentation and report | real-data tables collected; reported, not gated |

## 6. Acceptance

- **Network parity (Phase 1).** On the fixture's four inputs, the Fortran logits against PyTorch's:
  largest absolute difference at most T1. T1 is *measured* in Phase 0 as ten times the largest difference between
  PyTorch's own `float32` and `float64` logits on the fixture, and frozen in the Phase 0 entry before Phase 1
  starts. (Expected order 1e-4; that expectation is *guessed* and is not the criterion.)
- **Product paths (Phase 1).** The `sgemm` path and the explicit-loop path both pass fixture parity at T1.
- **Layer tests (Phase 1).** Each with a hand-computable expectation, exact to 1e-6 (*guessed*, on values of
  order 1): a 3x3 convolution of a delta input returns the flipped kernel at stride 1 and its even samples at
  stride 2; GroupNorm of a constant plane returns the bias; GroupNorm output has mean 0 and variance 1 per group
  with unit weight and zero bias; the eight symmetries applied to an image with one marked corner give eight
  distinct images whose set is closed under transposition.
- **Negative controls (Phase 1).** Each must fail and is then removed: one tensor read in row-major order;
  the stride-2 sampling offset by one pixel; the symmetry average reduced to four rotations with the input
  transposed.
- **Input parity (Phase 2).** For the real model's three reference stacks, the Fortran input tensors against
  `cnn_data.py`'s (written by the exporter with `--dump-inputs`), after undoing the transposition: correlation
  per class at least 0.999 (*guessed*; if SIMPLE's Fourier resampling treats the Nyquist line differently from
  numpy's and this fails, stop: decision 7 is then reopened, not the tolerance).
- **End-to-end parity (Phase 3).** `cavg_cnn_predict` against `cnn_ss.py predict` with the same model on the
  reference stacks: the same group for every class whose PyTorch top-two probabilities differ by more than T2,
  and |ΔP| at most T2 for every class and group. T2 is *measured* in Phase 2 as twice the largest |ΔP| obtained
  by feeding the Fortran input tensors to PyTorch, and frozen before Phase 3's run (replaced by decision 15).
- **Identity (Phase 3).** `cavg_cnn_predictions.txt` byte-identical at `nthr=1` and `nthr=16`, each run twice.
- **Existing behaviour (Phase 3).** The fast gate passes and `git diff --stat` shows no edit outside the file table.
- **Reported, not gated (Phase 4).** The predictions tables for the reference stacks and the run times at 1 and
  16 threads, with paths, for the maintainer's read.

## 7. Common rules and stop conditions

- No new command-line parameter. No tolerance moved after it is frozen. Fix as you go, with a file-table row
  before the edit.
- Compile-time policy: the type holds its weight arrays as allocatable components; module imports with `only:`;
  no new CMake target, no per-file flags.
- Stops:
  - the reference stacks or the trained `model.pt` are absent in Phase 0;
  - `testdata/` cannot take the 0.7 MB fixture (then decision 9 needs a ruling);
  - input parity fails (decision 7);
  - a file outside the table needs an edit;
  - decision 10's bound is missed.
- Gates no ruling may waive: fixture parity at T1; the three negative controls; byte identity across thread counts;
  the read-only contract (the project file's checksum before and after the program).

## 8. Phases

- **Phase 0, baseline.** Build Release and Debug at the base, run the fast gate, record the machine. Confirm which
  unit area runs `run_all_hhcnn_tests` and that `testdata/` is installed for unit tests. Write
  `scripts/sscnn/export_fortran.py`; export the trained model, the fixture and the input dumps; run
  `cnn_ss.py predict` on the three reference stacks (two `solve2D` runs of 100 classes at 1.3 Å and one of
  2 classes at 1.3 Å; the maintainer supplies the paths and the model) and keep the tables. Measure T1.
  Exit: both builds and the gate pass; the weights file, fixture and three reference tables exist; T1 is written
  in the Progress entry.
- **Phase 1, the network.** `simple_cavg_cnn` and its tester: reader, convolution, GroupNorm, residual block,
  global mean, linear layer, the eight-symmetry average, `predict` for a batch of prepared inputs.
  Exit: the layer tests and fixture parity pass in the unit suite in Release and Debug; each negative control was
  shown to fail and removed; the reader stops on a truncated file and on a wrong version (two tests).
- **Phase 2, the input.** `simple_cavg_cnn_input`: read class averages and halves, resample, normalise, flag empty
  classes. Measure T2.
  Exit: input parity holds on the three stacks; T2 is written in the Progress entry.
- **Phase 3, the program.** UI entry, exec route, commander, `run_cavg_cnn_predict`, the table writer.
  Exit: end-to-end parity, identity across thread counts, the project checksum unchanged, the fast gate, the
  speed bound of decision 10, `scripts/check_cli_chain.py` clean.
- **Phase 4, documents and report.** A section in `doc/how2s/how2train_cavg_cnn.md` (export, then
  `cavg_cnn_predict`), a row in `library_owners.md` for the network, the `simple-main-ui`/`cavg` guides if the
  program list there names programs, and `cavg_cnn_fortran_inference_report.md`.
  Exit: the report exists with the parity numbers, the run times and the reported-not-gated tables' paths.

### File table

| File | What | Phases |
| --- | --- | --- |
| `scripts/sscnn/export_fortran.py` (new) | weights, fixture and input dumps from `model.pt` | 0 |
| `testdata/cavg_cnn/` (new) | the fixture: weights, four inputs, PyTorch logits | 0 |
| `src/main/learn/simple_cavg_cnn.f90` (new) | the network: reader and forward pass | 1 |
| `src/main/learn/simple_cavg_cnn_tester.f90` (new) | layer tests, fixture parity, reader failures | 1 |
| `src/main/learn/learn.inf` | description entry, if the folder's listing needs it | 1 |
| unit-suite table that lists the `learn` testers | register the new tester | 1 |
| `src/main/cavg_quality/simple_cavg_cnn_input.f90` (new) | class averages to input tensors | 2 |
| `src/main/strategies/` or `src/main/cavg_quality/` routine `run_cavg_cnn_predict` (new file) | orchestration and table | 3 |
| `src/main/commanders/simple/simple_commanders_cavgs.f90` | thin commander | 3 |
| `src/main/exec/simple_exec_cavgproc.f90` | route | 3 |
| `src/main/ui/simple/simple_ui_cavgproc.f90` | program entry | 3 |
| generated UI defaults and args, if the build regenerates them | regenerated, not hand-edited | 3 |
| `doc/how2s/how2train_cavg_cnn.md`, `scripts/sscnn/README.md` | usage | 4 |
| `.github/skills/simple-modern-fortran/references/library_owners.md` | owner row | 4 |
| `doc/implementation_notes/completed/cavg_cnn_fortran_inference_report.md` (new; plan and report moved to `completed/`) | report | 4 |

## Progress

**Phase 0, baseline (2026-11-20, the development machine: Ubuntu, gfortran 15.2, 16 cores; base `a26f57437` with
the uncommitted `scripts/sscnn/`).**
- Release build, fast gate (15 of 15) and install pass on the base. A Debug build does not fit the local disk
  (1.1 GB free, a build takes 1.3 GB): it was run on a second machine (Phase 1). No production entry changes (the
  modules are new and nothing calls them), so no entry was run twice.
- `scripts/sscnn/export_fortran.py` written. The fixture (`testdata/cavg_cnn/`) and **T1 = 3.361e-6** (*measured*:
  ten times 3.361e-7, PyTorch's largest float32-float64 logit difference on it), frozen here.
- The trained model v3 (46 entries) exported (11 MB, outside the repository). Reference stacks, all at 1.3 A: a
  membrane complex's solve2D run (100 classes, 320 px), a small protein's classes (100, 240 px, no project: one
  made with `import_cavgs`), a transporter's (2, 196 px). Reference tables and input dumps written.
- Found while checking the premises: the Python resampling's two flaws (decision 12), fixed before the references
  were made. The other two stacks' predictions are unchanged by the fix.

**Phase 1, the network.** `simple_cavg_cnn` and its tester (a submodule). 41 checks pass in 0.2 s: the 3x3 and 1x1
convolutions against a brute force of PyTorch's formula (both products), the delta's flipped kernel, GroupNorm's
two closed forms, the symmetries, fixture parity (largest logit difference: sgemm 4.8e-7 plain and 1.2e-7 averaged,
loops 5.1e-7 and 8.2e-8, against T1 3.4e-6), the average's invariance to transposition, predict's probabilities
and its identity in a team of three, the reader's four refusals. Negative controls, each failing the suite and then
removed with the module restored byte for byte: the first block's kernel read row-major (logits off by 0.27), the
stride-2 sampling offset by one pixel (0.19, and the brute-force convolution test), the average cut to the four
rotations (0.04, and the invariance). Debug: a copy of the working tree built with `compile_debug.sh` on the second
machine (GCC 15.1 toolset, runtime checks), its fast gate 15 of 15, and after the fixture's rename both suites again.

**Phase 2, the input.** `simple_cavg_cnn_input` and its tester: 12 checks (exact resampling of band-limited images
at even and odd sizes, down and up, to 1e-9; the centre pixel; the scaling; empty classes; a team of three). A
temporary harness, removed after its run, compared the three stacks with Python's input dumps: correlation 1.0000000
for every class, largest difference 6e-8, the same empty classes (4, 11, 0), orientation as decision 17 says.
**T2 = 7.466e-6** (decision 15: input-caused |dP| 5.6e-9, T1_real 1.49e-5), frozen before Phase 3.

**Phase 3, the program.** Descriptor, route, commander, body. End to end against PyTorch's full-precision
probabilities: largest |dP| 7.0e-8, 2.0e-7 and 2.1e-9 on the three stacks (T2 7.5e-6); no group differs among the
187 decided classes; the empty classes agree. The table is byte-identical at `nthr=1` and `nthr=16`, two runs each.
The projects' checksums are unchanged. 100 classes of 320 px: 3.1 s at 16 threads, 31 s on one (bound 60 s). Fast
gate 15 of 15, `unit_ui` 1,042 checks, the macro scans clean. In the Debug build the program ran both stacks
under bounds checks; its MotAB probabilities differ from PyTorch's by 1.5e-7 and from the Release table by 7.9e-8.

**Phase 4, documents.** The guide (`doc/how2s/how2train_cavg_cnn.md`, a section on predicting inside SIMPLE), the
toolkit's README and CNN.md, two rows in `library_owners.md`, a line in the class-average quality guide; the
report beside this plan. The three tables are kept for the maintainer's read (report, section 6).

What it means, in short: SIMPLE runs the class-average CNN itself, with the same probabilities as PyTorch to 2e-7,
in 3 s per 100 classes, and changes nothing until a later plan decides how to use them.
