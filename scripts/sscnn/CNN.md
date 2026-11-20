# The class-average CNN

A convolutional network that sorts 2D class averages into three groups:

- **junk**;
- **secondary structure present**;
- **likely protein**: protein without visible secondary structure.

It is trained only on simulated class averages from `sscnn.py build`: EMDB maps run through SIMPLE's own particle simulation and solve2D, so every particle of every class is known.

**Inside SIMPLE:** `simple_exec prg=cavg_cnn_predict` runs a trained model in Fortran (export it with `export_fortran.py`), with the same probabilities as `cnn_ss.py predict`; see [doc/how2s/how2train_cavg_cnn.md](../../doc/how2s/how2train_cavg_cnn.md).

**Status (2026-11-20):** trained on 46 EMDB entries (15,360 simulated classes) with 5-fold cross-validation by protein; see [Results](#results). Protein-against-junk selection has been judged good by eye on two real datasets; whether the secondary-structure output carries over to real data is open.

```text
 simulated set (sscnn.py build)                   real data
 runs/<id>_r1, <id>_r2,                           cavgs_iterNNN.mrc,
 runs/<id>_lp<L>_r1, <id>_pr<L>_r1                its _even and _odd halves
            │                                               │
            ▼                                               ▼
 cnn_data.py                                      cnn_ss.py predict
   every class: average, even, odd                  resampled to 2 Å, 160 px
   label from its particles and the run's map       normalised the same way
   → x.npy, meta.tsv, info.json                             │
            │                                               │
            ▼                                               │
 cnn_ss.py train                                            │
   5-fold cross-validation by protein → cv_report.html      │
   final model on every class         → model.pt ───────────┤
                                                            ▼
                                              preds.tsv, preds.html
                                              (and agreement with labels.json)
```

## Groups and labels

`cnn_data.py` labels every class of every run:

| Group | Rule |
| --- | --- |
| junk | more than half the class's particles are junk |
| secondary structure | otherwise, the run used the full map (`<id>_r<k>`) |
| likely protein | otherwise, the run used a low-passed (`<id>_lp<L>_r1`) or phase-randomised (`<id>_pr<L>_r1`) map |

**Protein fraction:** each class's share of protein particles. It's exact: the run records which particles are junk.

**The secondary-structure label is an assumption.** It comes from the run's map, not from what shows in the class:
- a full-map class is labelled secondary structure even when its helices don't show (low SNR, or helices along the beam);
- a low-passed class is never labelled secondary structure.

**Phase-randomised runs count as low-passed.** Their map is low-passed at L before the phases finer than L are randomised. The fine texture they carry is not secondary structure, which is why they are in the set: they teach the network that fine detail alone isn't enough.

**Visibility numbers** (the helix and strand visibility at 12–6 Å from `sim_truth.py`) are kept beside every class. They are used to find classes whose label is probably wrong, and they can become a measured label later (see [Limits](#limits-and-next-steps)).

**Mixed classes** take the label of their majority. The protein output still learns their protein fraction.

## Input

**Three channels:** the class average and its even and odd halves.
- The halves carry the noise: detail present in both is signal, detail in one only is noise.
- `--channels avg` trains on the average alone, for class averages from software that doesn't write halves.

**Sampling:** 2 Å per pixel, 160 px (a 320 Å box), the simulated set's sampling. Nyquist is 4 Å. Helices are about 10 Å apart, so they are well sampled. Strands are 4.8 Å apart, close to Nyquist, so sheets show mainly as slabs.

**Normalisation:** each channel has the average's mean (inside the inscribed circle) subtracted and is divided by the average's standard deviation. All three channels use the average's numbers, so the halves' extra noise stays visible relative to the average.

**Real data:** `predict` resamples to 2 Å in Fourier space, then crops or pads to 160 px at the centre. Particles wider than 320 Å are cut off.

## Architecture

```text
 3 × 160 × 160   average, even, odd
       │
 stem: 3×3 conv, GroupNorm, ReLU                          32 × 160²   2 Å/px
       │
 stage 1 ┐                                                64 × 80²    4 Å/px
 stage 2 │  residual block:                              128 × 40²    8 Å/px
 stage 3 │   3×3 conv stride 2 ─ GN ─ ReLU ─ 3×3 conv ─ GN 128 × 20²   16 Å/px
 stage 4 │   + 1×1 conv stride 2 ─ GN on the skip         256 × 10²   32 Å/px
 stage 5 ┘   ─ ReLU                                       256 × 5²    64 Å/px
       │
 global average pool                                      256
       │
 dropout 0.3, linear                                      2 logits: protein; secondary structure given protein
```

| Layer | Output | Pixel | Receptive field |
| --- | --- | --- | --- |
| input | 3 × 160 × 160 | 2 Å | — |
| stem | 32 × 160 × 160 | 2 Å | 6 Å |
| stage 1 | 64 × 80 × 80 | 4 Å | 18 Å |
| stage 2 | 128 × 40 × 40 | 8 Å | 42 Å |
| stage 3 | 128 × 20 × 20 | 16 Å | 90 Å |
| stage 4 | 256 × 10 × 10 | 32 Å | 186 Å |
| stage 5 | 256 × 5 × 5 | 64 Å | 378 Å |
| pool, linear | 2 | — | whole box |

**Size:** 2.77 M weights at the default width (`--width 32`: the stem's channels; the stages scale with it).

**Design reasons:**
- **Full-resolution first layer:** 4–10 Å detail is seen before anything is downsampled.
- **Receptive field:** by stage 2 a unit sees 42 Å, a few helices side by side; by stage 5 it sees the whole box.
- **Global pooling:** the answer doesn't depend on where in the box things sit.
- **GroupNorm rather than BatchNorm:** no statistics are carried over from simulated to real data, and it behaves the same at every batch size.

## Outputs

| Output | Meaning | Target | Trained on |
| --- | --- | --- | --- |
| P(protein) | the class is protein | the class's protein fraction (soft) | every class |
| P(ss \| protein) | its secondary structure shows | 1 for secondary structure, 0 for likely protein | protein classes only |

**The three groups:**
- junk = 1 − P(protein);
- secondary structure = P(protein) × P(ss | protein);
- likely protein = P(protein) × (1 − P(ss | protein)).

The predicted group is the largest of the three.

**Why two outputs rather than one three-way output:**
1. Mixed classes give soft targets instead of being forced to one side.
2. The secondary-structure question is only asked of protein, so junk classes don't distort it.
3. The two cut-offs can be set separately. For sieving, for instance, a class might be discarded only when P(junk) > 0.9, since discarding real particles is the costly mistake.

## Training

- **Loss:** binary cross-entropy on P(protein), plus binary cross-entropy on P(ss | protein) for protein classes. Each class is weighted so that the three groups weigh the same; the weights come from the training classes' counts.
- **Optimiser:** AdamW, learning rate 1e-3, weight decay 0.05, batches of 32.
- **Schedule:** one epoch of linear warm-up, then cosine decay to zero over `--epochs` (40).
- **No choice on validation data:** the number of epochs is fixed, with no early stopping and no epoch picked by validation loss. So the cross-validation numbers are not tuned to the validation classes.
- **Augmentation:** for every class in every epoch:
  - one of the 8 rotations by 90° and mirrors (exact; interpolated rotations would blur the 4–8 Å detail);
  - a circular shift of up to 4 px (8 Å);
  - a swap of the even and odd halves;
  - **fewer particles** (half the classes, `--aug-noise 0.5`): independent noise added to each half as if the class had only a fraction f of its particles, f drawn log-uniformly between 0.05 and 1 (`--aug-noise-min`), so down to a twentieth of the population. The noise has the class's own spectrum and radial profile, from (even − odd)/2, and the average gets the mean of the two additions;
  - **unmasked background** (half the classes, `--aug-bg 0.5`): for a masked class (where its noise falls below a quarter of its level), noise at the class's own level added outside the mask, out to a random radius up to the whole box, as in real class averages that aren't masked tightly.

  After these two, the three channels are scaled again as in the data step. Both need the halves, so they're made on all three channels before `--channels avg` keeps only the average. Checked on the local classes: the fewer-particles noise grows the halves' difference by exactly 1/f, and just outside the mask the background noise rises from 0.23 to 0.74 of the level inside.
- **Evaluation:** logits averaged over the 8 rotations and mirrors.
- **Devices:** CUDA (bfloat16 where supported), Apple MPS, or CPU.

## Validation

- **Cross-validation by protein:** the entries are split into `--folds` (5) folds at random, with a fixed seed. Each fold's model trains on the other folds' proteins and predicts its own, so every class is predicted by a model that never saw its protein. With 15 entries that's 3 proteins per fold. (`--split run` makes folds of runs instead; only for testing the code on one protein, because it leaks the protein's shape.)
- **Metrics on those predictions:**
  - the confusion matrix;
  - recall and precision per group;
  - balanced accuracy (the mean recall), accuracy, and log loss of the three group probabilities;
  - AUC of P(protein), protein against junk;
  - AUC of P(ss | protein), secondary structure against likely protein, over protein classes.
- **Spread:** the same per fold and per protein. A protein the model fails on shows up here.
- **By kind of class:** the share predicted in each group for:
  - full-map, low-passed and phase-randomised classes, by cutoff and by particle-SNR band (< 0.003, 0.003–0.01, ≥ 0.01). The band uses `snr_px`, the run's noise per pixel expressed as a 150 Å particle's SNR, so runs from before and after the SNR scaling compare;
  - junk classes by their most common junk kind.
- **Noise-level check:** the secondary-structure AUC within each SNR band. If it is well below the overall AUC, the network is reading the noise level rather than the secondary structure. On the local runs (one protein, seven runs) it did exactly that: AUC 0.32 overall, about 0 within each band, because the one noisy full-map run made "noisy" mean "secondary structure". In the full set every run draws its SNR independently, which should break this; giving each entry's variant runs its full-map runs' SNRs would remove it by design.
- **Error galleries:** for each wrong (label, prediction) pair, the most confident errors, captioned with the helix visibility (protein classes) or the junk kind. Full-map classes with low visibility among the errors are probably label noise, not model errors.
- **Fit:** the log gives each fold's balanced accuracy on its own training classes. A low value means the model is underfitting.
- **Curves:** training and validation loss per epoch and fold, in `cv_curves.tsv`.
- **Final model:** trained on every class with the same settings. The cross-validation numbers are its expected performance.
- **Real data:** `predict --labels` compares with labels made in `ss_trainer.py` (yes: secondary structure visible; no: not visible, whether protein or junk). It reports the AUC of P(secondary structure) for yes against no, and the groups predicted for each label.

## Running it

```
python3 sscnn.py train   --set SET_DIR --model-dir MODEL_DIR [-- options of cnn_ss.py train]
python3 sscnn.py predict --model-dir MODEL_DIR --stk cavgs_iterNNN.mrc [--cls classdoc_ranked.txt] [--name NAME]
python3 sscnn.py report  --model-dir MODEL_DIR --set SET_DIR
```

`sscnn.py` makes the environment first (numpy and PyTorch in `~/.cache/simple/sscnn_venv`, or `--venv`; the CUDA build when there is a GPU). The same steps without it, with that environment's `python`: `python -I cnn_data.py --set SET_DIR --out MODEL_DIR/data`, then `cnn_ss.py train --data MODEL_DIR/data --out MODEL_DIR`, `cnn_ss.py predict --model MODEL_DIR/model.pt --stk ... --out preds.tsv --html preds.html [--cls ...]` and `cnn_ss.py summary --model-dir MODEL_DIR [--set SET_DIR]`.

**Options of `train`** (after `--` with `sscnn.py train`, or `--train-args` with `all`):

| Option | Default | Meaning |
| --- | --- | --- |
| `--folds` | 5 | 0: the final model only |
| `--no-final` | — | cross-validation only |
| `--epochs` | 40 | |
| `--batch`, `--lr`, `--wd` | 32, 1e-3, 0.05 | |
| `--width` | 32 | the stem's channels |
| `--dropout` | 0.3 | |
| `--channels` | `avg+halves` | or `avg` |
| `--aug-noise`, `--aug-noise-min` | 0.5, 0.05 | share of classes given fewer particles' noise; the least particle fraction |
| `--aug-bg` | 0.5 | share of classes given noise outside their mask |
| `--device` | `auto` | `cpu`, `cuda`, `cuda:1`, `mps` |
| `--nthr` | PyTorch's default | CPU threads |

**In `OUT_DIR`:**

| File | Content |
| --- | --- |
| `cv_report.html`, `cv_report.txt` | everything under Validation, with the galleries; the HTML pages (and `predict`'s) have brightness and contrast sliders for the class images |
| `cv_predictions.tsv` | every class: its record, fold, the three probabilities, the predicted group |
| `cv_metrics.json`, `cv_curves.tsv` | the metrics; the loss curves |
| `model.pt` | the final model, its settings and the cross-validation summary (loads with `torch.load(..., weights_only=True)`) |
| `train.log`, `cnn.log` | the logs |
| `data/` | the collected classes: `x.npy` (16-bit), `meta.tsv`, `info.json`. Delete it to collect again |

**Time:** about 5 s per epoch per 140 classes on 16 CPU cores. The 15-entry set has about 3,600 classes, which is about 100 s per epoch on CPU, so 5 folds plus the final model at 40 epochs take about 7 hours. A GPU should take well under an hour.

## Results

Three generations, each cross-validated by protein (every class predicted by a model that never saw its protein), 40 epochs:

| | v1 | v2 | v3 |
| --- | --- | --- | --- |
| proteins / classes | 13 / 2,251 | 32 / 5,310 | 46 / 15,360 |
| variants per entry | one lp, one pr (rotating) | the same | lp and pr at 15, 30 and 50 Å |
| SNR scaled by size, augmentations | no | yes | yes |
| balanced accuracy | 0.871 | 0.876 | 0.911 |
| AUC protein vs junk | 0.995 | 0.976 | 0.995 |
| AUC secondary structure vs likely protein | 0.922 | 0.936 | 0.961 |
| protein classes called junk | 2.8% | 7.9% | 4.7% |
| pr15 classes called secondary structure (high SNR) | — | 58% | 22% |

**v3 in detail** (2026-11-20; MspA and a SARS-CoV-2 RBD–Fab entry were skipped by the map checks):

| Label \ predicted | Junk | Secondary structure | Likely protein | Recall |
| --- | --- | --- | --- | --- |
| junk | 8,727 | 12 | 45 | 0.99 |
| secondary structure | 62 | 1,451 | 230 | 0.83 |
| likely protein | 245 | 209 | 4,379 | 0.91 |

- **Junk:** 99–100% of every junk kind (empty, edge, blob, aggregate, gradient) is recognised.
- **Noise level:** the secondary-structure AUC is 0.89, 0.96 and 0.99 in the particle-SNR bands < 0.003, 0.003–0.01 and ≥ 0.01. It falls only where helices genuinely fade into the noise, so the network doesn't read the noise level.
- **Hard proteins:** the errors cluster by protein. Lowest balanced accuracy:
  - HMG-CoA synthase (EMD-73971), 0.65;
  - serum albumin with four nanobodies (EMD-72279), 0.63: its lobed classes look like aggregate junk;
  - the 20S proteasome (EMD-70078), 0.67;
  - the BAM complex (EMD-71970), 0.67, though its AUCs are 1.00 and 0.98: its classes are ranked right but fall on the wrong side of the cut-off;
  - 14-3-3 (EMD-73503), 0.73.
- **Epochs:** validation loss is lowest near epoch 29 and slightly higher by 40, while each fold fits its own training proteins at 0.96–0.97. More epochs don't help; more proteins do.
- **Real data** (three SIMPLE solve2D runs, no ground truth):
  - **protein against junk:** judged good by eye on a small protein of about 60–80 Å (89 classes) and on a membrane complex (100 classes, most with few particles);
  - **secondary structure:** both classes of a transporter at 8.8 and 10.2 Å are called secondary structure (P ≈ 1.0). No class of the other two is, including their best-resolved ones (7–13 Å by SIMPLE's estimate), which show no rod-like texture by eye.

## Limits and next steps

- **Label noise in secondary structure:** full-map classes whose helices don't show are labelled secondary structure. The error galleries show these. A measured label would replace the run's kind with a visibility threshold, set from your own yes/no labels on a few hundred simulated classes spanning the visibility range.
- **Blobs against low-resolution protein:** at 30–50 Å a protein class is a smooth outline, much like an ethane blob or an aggregate. One class average may not settle it. Whether the outline matches the dataset's other classes would, which needs a model that sees several classes at once.
- **Simulated junk is five simple kinds.** Real ice crystals, cracks and protein aggregates need real labelled examples.
- **Simulation against real data:**
  - Gaussian noise, no beam-induced motion, no conformational heterogeneity, CTF known exactly;
  - nothing here is trustworthy before it holds on labelled real classes.
- **46 proteins are still few shapes.** The per-protein results in the cross-validation show how much a new shape costs.
- **Fixed sampling:** 2 Å and 320 Å; wider particles are cut off.
- **Halves only from SIMPLE.** Class averages from other software need a `--channels avg` model.
- **Rotation:** the network is not rotation-invariant by construction. It learns invariance from the 8 square symmetries and the data, and evaluation averages over those 8.
