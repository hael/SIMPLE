# How to sort 2D class averages with the class-average CNN

`scripts/sscnn/` holds a convolutional network that sorts 2D class averages into three groups:

- **junk**: ice, carbon edges, ethane, aggregates, empty classes;
- **secondary structure present**: protein with visible helices or sheets;
- **likely protein**: protein without visible secondary structure.

It is trained only on simulated data. EMDB maps are run through SIMPLE's own particle simulation and `solve2D`, so the truth of every class is known. One script, `scripts/sscnn/sscnn.py`, does everything: the Python environment, the training set, training with cross-validation, prediction on your class averages, and a summary report.

The model and its validation are described in [scripts/sscnn/CNN.md](../../scripts/sscnn/CNN.md), and the simulation in [scripts/sscnn/README.md](../../scripts/sscnn/README.md).

**Status (2026-11-20):** on held-out proteins, protein against junk is near perfect (AUC 0.995) and secondary structure against likely protein reaches AUC 0.96. Protein-against-junk selections looked good on two real datasets; whether the secondary-structure output carries over to real class averages is still open.

## Requirements

- **SIMPLE**, built. `sscnn.py` uses `--simple`, else `$SIMPLE_PATH`, else the SIMPLE install it sits in, else the repository's `build/`. Only training needs SIMPLE.
- **Python 3.9 or newer.** `sscnn.py` uses the standard library only, and on first use makes its own environment with numpy and PyTorch:
  - location: `~/.cache/simple/sscnn_venv` by default, or `--venv`, or `$SSCNN_VENV`;
  - PyTorch: the CUDA build when `nvidia-smi` finds a GPU, the CPU build otherwise (`--torch`).
- **Internet access:** to PyPI for the environment, and to EMDB and RCSB for training.
- **For training:** a large machine. The full set of 48 entries is 384 SIMPLE runs: about 13 hours on 72 cores (5 entries at once, 14 threads each, about 29 runs an hour), plus about 3 hours of training on one GPU, and about 25 GB of disk.

`make install` copies the scripts to `$SIMPLE_PATH/scripts/sscnn/`. None of them is run by SIMPLE's build or fast gate.

## Predict on your class averages

With a trained model (a `model.pt`, or a model folder from training):

```bash
python3 scripts/sscnn/sscnn.py predict --model-dir MODEL_DIR \
  --stk proj/2_solve2D/cavgs_iter030.mrc --cls proj/2_solve2D/classdoc_ranked.txt
```

- **Input:** the class averages and their `_even` and `_odd` halves, which `solve2D` writes beside them. Any pixel size and box works; the images are resampled to the model's 2 Å and 160 px.
- **`--cls`** takes SIMPLE's class document. It adds each class's population and resolution, and the number of particles in the classes kept.
- **Several stacks:** repeat `--stk` (and `--cls`).
- **A model file elsewhere:** `--model path/to/model.pt`.
- **Output:** written to `MODEL_DIR/predictions/NAME.*`:

  | File | Content |
  | --- | --- |
  | `NAME.tsv` | per class: the three probabilities, the predicted group, population and resolution |
  | `NAME.html` | the class averages grouped by prediction, each captioned with its probability, particles and resolution, with brightness and contrast sliders |
  | `NAME.json` | the record the summary report reads |

The current best model's Fortran weights are in `scripts/sscnn/weights/cavg_cnn_v3.bin` (installed to
`$SIMPLE_PATH/scripts/sscnn/weights/`): model v3, trained 2026-11-20 on 46 simulated proteins with full variant
coverage. The file is kept out of git for now (`.gitignore`), so a fresh clone lacks it: export it from a `model.pt`
as below, or copy it from a colleague.

## Predict inside SIMPLE, without Python

SIMPLE runs a trained model itself, in Fortran with the libraries its build already links. Export the model once:

```bash
python3 scripts/sscnn/sscnn.py setup
~/.cache/simple/sscnn_venv/bin/python -I scripts/sscnn/export_fortran.py --model MODEL_DIR/model.pt --out cavg_cnn_weights.bin
```

then, on any project with class averages:

```bash
simple_exec prg=cavg_cnn_predict projfile=proj.simple infile=$SIMPLE_PATH/scripts/sscnn/weights/cavg_cnn_v3.bin nthr=16
```

- **Output:** `cavg_cnn_predictions.txt` in the execution folder, one line per class: `class`, `pop`, `res`, `p_junk`, `p_ss`, `p_likely` and `group` (`junk`, `secondary_structure`, `likely_protein`, or `empty` for a class without particles).
- **Read-only:** the program reads the project's class averages, their `_even` and `_odd` halves and the class populations, and changes nothing in the project. It makes no selection.
- **The same numbers as Python:** it prepares the images exactly as `cnn_ss.py predict` does and runs the same network. On three real stacks (205 classes) its probabilities differ from PyTorch's by at most 2 × 10⁻⁷, every group and every empty class agrees, and the table is byte-identical at any thread count.
- **Speed:** 100 class averages of 320 px take about 3 s at 16 threads and 31 s on one.
- **The code:** the network is `src/main/learn/simple_cavg_cnn.f90`, the input preparation `src/main/cavg_quality/simple_cavg_cnn_input.f90`, the program's body `src/main/cavg_quality/simple_cavg_cnn_predict.f90`. Tests: `simple_test_exec test=unit_numerics suite=cavg_cnn` and `suite=cavg_cnn_input`, against a small committed model (`testdata/cavg_cnn/`) whose PyTorch outputs they reproduce.

## Label real class averages

To collect labelled real class averages for training, make a labelling page where the data are:

```bash
python3 scripts/sscnn/label_cavgs.py proj/2_solve2D/proj.simple run/3_cavg_cnn_predict --dataset NAME --out labels.html
```

- **Inputs, in any mix:** a project with class averages (a `.simple` file, such as `solve2D`'s); a
  `cavg_cnn_predictions.txt`; or a folder, searched for those tables (a `cavg_cnn_predict` folder, or a
  `scripts/mic_training_generator.py` output with all its runs).
- **The labels:** secondary structure (helices or sheets visible), likely protein, protein with secondary structure
  unsure, junk, and can't tell.
- **`--blind`** leaves the CNN's calls out, so they cannot sway the labels. Otherwise the classes of a prediction
  table are grouped by the CNN's call.
- **`--dataset`** names the data in the labels. Give runs of the same sample the same name: cross-validation will hold
  out whole datasets.
- **Labelling:** the page is one self-contained HTML file. Copy it anywhere and open it in a browser, without a
  server or Python. Labels are kept in the browser as you go. **Export TSV** writes them, one row per class
  (stack, class, population, resolution, label, the CNN's call, curator, time), and **Import TSV** loads an export
  to continue.
- **Label every class of a run**, not only those the CNN gets wrong, and label some runs blind.
- **Not yet used by training:** the training does not read the labels yet.

## Train a model

The whole chain in one command:

```bash
python3 scripts/sscnn/sscnn.py all --set SET_DIR --model-dir MODEL_DIR \
  --stk proj/2_solve2D/cavgs_iter030.mrc --cls proj/2_solve2D/classdoc_ranked.txt
```

or step by step:

```bash
python3 scripts/sscnn/sscnn.py setup
python3 scripts/sscnn/sscnn.py build  --set SET_DIR --jobs 4 --nthr 14
python3 scripts/sscnn/sscnn.py train  --set SET_DIR --model-dir MODEL_DIR
python3 scripts/sscnn/sscnn.py report --set SET_DIR --model-dir MODEL_DIR
```

- **`build`** reads `scripts/sscnn/entries.tsv` (or `--entries`). Per entry it:
  - checks EMDB's reported resolution (it must be better than 4.5 Å);
  - downloads both half maps and the fitted model;
  - makes low-passed and phase-randomised versions of the map at 15, 30 and 50 Å, the negatives;
  - simulates particles in eight runs with junk added, runs `solve2D`, and records the ground truth of every class.

  It is restartable: finished runs are kept, and new entries, or variants missing from an older set, are added on the next run.
- **`train`** collects the classes into `MODEL_DIR/data` (`--recollect` after the set has grown). It then runs 5-fold cross-validation by protein, so a protein is never in both training and validation, and trains the final model, `MODEL_DIR/model.pt`.
- **Options of the underlying scripts:** after `--` for one step, or as quoted strings with `all` (`--build-args`, `--train-args`, `--predict-args`). For example:

  ```bash
  python3 scripts/sscnn/sscnn.py train --set SET_DIR --model-dir MODEL_DIR -- --epochs 30
  ```

- **On a remote machine,** run it detached:

  ```bash
  nohup python3 scripts/sscnn/sscnn.py all ... > sscnn.log 2>&1 &
  ```

  Each step also logs to `SET_DIR/build.log` or `MODEL_DIR/cnn.log`.
- **Adding entries:** append them to `entries.tsv` and never reorder it. Find candidates with `scripts/sscnn/emdb_search.py`, which filters by resolution, half maps, a fitted model and molecular weight. EMDB's weight field is often wrong, so check a candidate's deposited model as well.

## Read the reports

| Page | What to look at |
| --- | --- |
| `MODEL_DIR/summary.html` | the set (each entry's status, size, map-model correlation and runs), the training data, the model's settings, the cross-validation and every prediction, with links |
| `MODEL_DIR/cv_report.html` | the confusion matrix, results per fold and per protein, a breakdown by kind of class (map, low-pass and phase-randomised cutoffs by particle SNR, junk by kind), and galleries of the most confident errors |
| `MODEL_DIR/predictions/NAME.html` | your class averages by predicted group |

- **The per-protein table** is the honest measure of how the model does on a protein it has never seen. Its errors cluster by protein.
- **The SNR-band line** in the cross-validation guards against a shortcut: if the secondary-structure AUC within each particle-SNR band falls well below the overall AUC, the network is reading the noise level instead of the structure.

## Limits

- **The secondary-structure label is the simulation's map, not what shows.** A class from the full map counts as secondary structure even when noise or the view hides its helices. The output therefore means "consistent with a full-resolution protein map" more than "helices visible".
- **Simulated junk is five simple kinds.** Lobed multi-domain complexes (an antigen with several nanobodies, for example) can look like the simulated aggregates.
- **Fixed sampling:** 2 Å per pixel in a 320 Å box. Wider particles are cropped.
- **The halves are needed.** Class averages from other software, which has no halves, need a model trained with `--channels avg`.
