# scripts/sscnn: a CNN for 2D class averages

Tools that sort 2D class averages into junk, secondary structure present and likely protein, and the simulated training data behind them:

- **a simulation pipeline** that makes class averages with SIMPLE itself, from EMDB maps (or atomic models), with exact helix and strand ground truth for every class, negatives that keep the particle but lose its secondary structure, and junk classes;
- **the class-average CNN** (`cnn_data.py`, `cnn_ss.py`, described in [CNN.md](CNN.md)), trained on those class averages;
- **a labelling and selection tool** that learns from your labels and exports a selection to SIMPLE.

The guide is [doc/how2s/how2train_cavg_cnn.md](../../doc/how2s/how2train_cavg_cnn.md). `make install` copies this folder to `$SIMPLE_PATH/scripts/sscnn/`. The scripts need numpy and PyTorch, which `sscnn.py` installs into its own environment; none of them is run by SIMPLE's build or fast gate.

## Running everything (`sscnn.py`)

One script does it all: the environment, the simulated set from EMDB, training with cross-validation, prediction on your class averages, and a summary report. From the repository (or `$SIMPLE_PATH/scripts/sscnn/` of an install):

```
python3 sscnn.py all --simple /path/to/SIMPLE/build --set SET_DIR --model-dir MODEL_DIR \
    --stk my_project/2_solve2D/cavgs_iter030.mrc --cls my_project/2_solve2D/classdoc_ranked.txt
```

or step by step:

```
python3 sscnn.py setup                                            # environment only
python3 sscnn.py build   --simple /path/to/SIMPLE/build --set SET_DIR --jobs 4 --nthr 14
python3 sscnn.py train   --set SET_DIR --model-dir MODEL_DIR       # --recollect after the set grew
python3 sscnn.py predict --model-dir MODEL_DIR --stk cavgs_iterNNN.mrc [--cls classdoc_ranked.txt] [--name NAME]
python3 sscnn.py report  --model-dir MODEL_DIR --set SET_DIR       # MODEL_DIR/summary.html
```

| Step | What it does |
| --- | --- |
| `setup` | a virtual environment with numpy and PyTorch (`--venv`, or `$SSCNN_VENV`, default `~/.cache/simple/sscnn_venv`): the CUDA build when `nvidia-smi` finds a GPU, the CPU build otherwise (`--torch`). Skipped when it already works; every other step runs it first |
| `build` | the simulated set (`build_emdb_set.py`), below |
| `train` | the classes collected into `MODEL_DIR/data` (`cnn_data.py`), then 5-fold cross-validation by protein and the final model (`cnn_ss.py train`); see [CNN.md](CNN.md) |
| `predict` | the model on class averages (each `--stk` with its `_even` and `_odd` halves; `--cls` with SIMPLE's class document adds populations, resolutions and the particles kept) into `MODEL_DIR/predictions/NAME.{tsv,html,json}` |
| `report` | `MODEL_DIR/summary.html`: the set (each entry's status and runs), the training data, the model, the cross-validation and every prediction, linked to `cv_report.html` and the prediction pages |
| `all` | `setup`, `build`, `train`, `predict` (when `--stk` is given) and `report` |

- **Options of the underlying scripts:** after `--` for a single step (`python3 sscnn.py train ... -- --epochs 30`), or `--build-args`, `--train-args` and `--predict-args` (quoted strings) with `all`.
- **Requires:**
  - SIMPLE built: `--simple`, else `$SIMPLE_PATH`, else the install this script sits in, else the repository's `build/` (the folder holding `bin/simple_exec`). Build a set with one SIMPLE commit throughout, since class averages depend on solve2D;
  - Python 3.9 or newer (the driver needs nothing else; it makes the environment);
  - internet access to EMDB, RCSB and PyPI.
- **Long runs on a remote machine:** `nohup python3 sscnn.py all ... > sscnn.log 2>&1 &`. Each step's output also goes to `SET_DIR/build.log` or `MODEL_DIR/cnn.log`.

### The set (`build`)

- **Entries:** it builds the 48 entries of `entries.tsv` (`--entries FILE` for another list): 15 larger complexes, 18 proteins of 44–133 kDa, and 15 more added 2026-11-19 (small Fab and nanobody complexes, membrane proteins including a GPCR–Gs complex and the BAM β-barrel complex, larger soluble proteins). Keep the order of `entries.tsv`: with `--coverage rotate` the cutoffs depend on an entry's position. An entry is skipped unless its reported resolution is better than `--max-res` (4.5 Å).
- **Logs:** the whole run's output goes to `OUT_DIR/build.log`, each entry's to `OUT_DIR/logs/<id>.log`, and a line per finished entry to `OUT_DIR/entries.tsv`.
- **Preflight:** before downloading anything, it runs one small run (300 particles, 4 classes, 64 px) on a made-up map through every SIMPLE step the build uses, about 15 s. If SIMPLE does not run on this machine (another version, a missing library, a wrong environment), it stops there with the end of SIMPLE's log and keeps the files in `OUT_DIR/preflight`. `--no-preflight` skips it.
- **Failures:** a failed entry's line in `entries.tsv` and `build.log` ends with the end of its error; `OUT_DIR/logs/<id>.log` has all of it, including the end of the failing SIMPLE program's log. A failed run's folder `OUT_DIR/runs/<name>_work` keeps every log of that run.
- **Restarting:** run the same command again; finished runs are kept, failed entries are tried again.
- **Parallelism:** `--jobs N` runs N entries at once, each with `--nthr` threads. On a large machine several entries at 16 threads use the cores better than one entry at all of them.
- **Disk:** an entry starts only with `--min-free-gb` free (default 3). Each job's peak is about 1 GB (the particle stack of a run) plus up to 1 GB of downloads while it prepares its map. Each run keeps about 50 MB, or about 40 MB with `--compact`.
- **Time:** about 9–10 minutes per run at the defaults on 16 threads here, so about 80 minutes per entry (eight runs). On a 72-core machine with `--jobs 4 --nthr 14`, about 27 runs per hour.

`SIMPLE_QSYS` defaults to `local` and `SIMPLE_EMAIL` to a placeholder; set them if your SIMPLE setup needs others.

## What it makes

Per entry (an EMDB map with both half maps and a fitted model):

| Run | Input map | Particle SNR | Purpose |
| --- | --- | --- | --- |
| `<id>_r1`, `<id>_r2` | the map | drawn log-uniform from 0.001–0.03 | positives; visibility varies with resolution and view |
| `<id>_lp<L>_r1` | low-passed at L Å | drawn | negatives: the particle's shape, no secondary structure |
| `<id>_pr<L>_r1` | phases randomised beyond L Å, amplitudes kept | drawn | negatives: high-resolution texture, no secondary-structure geometry |

L is each of 15, 30 and 50 Å (`--cutoffs`), so an entry has 8 runs: 2 of the map, 3 low-passed, 3 phase-randomised (`--coverage full`, since 2026-11-19; `--coverage rotate` gives each entry one low-passed and one phase-randomised run, their cutoffs cycling between entries, as the first sets had). Entries built under rotate get their missing variants when the set is built again: the maps are downloaded and prepared again. At 30–50 Å only the rough outline survives, as in many real class averages.

Every run simulates, with SIMPLE:

- **protein particles:** 5,000, in random views, with CTF (defocus 1.0–2.5 µm), B-factor 150 Å² and picking error ±6 Å; 30% of them sit on an ice gradient;
- **noise:** the drawn SNR is that of a 150 Å particle. SIMPLE scales its noise with the image's variance, which grows as the fourth power of the particle's diameter, so it gets the drawn SNR × (diameter / 150 Å)⁴. The noise per pixel is then the same for every protein, as the dose makes it in real data, and a small protein stands out less (`--snr-ref-diam`; runs built before 2026-11-19 used the drawn SNR as it was);
- **junk particles:** as many again (`--junk-frac 1.0`: half of all particles), in equal shares of five kinds, with the same CTF and noise level:
  - empty ice;
  - edges (cracks, carbon);
  - dense round blobs (ethane);
  - amorphous aggregates;
  - ice gradients;
- **solve2D** into 60 classes;
- **ground truth** for every class from noise-free particles in the same views.

```text
 entries.tsv ─► build_emdb_set.py (one entry at a time, or --jobs at once)
    download half maps + model ─► emdb_prep.py ─► maps/<id>/_maps_<id>/   full, helix, strand maps, model.json
                                   map_variants.py ─► variants/<id>/_maps_<id>_lp<L>/, _pr<L>/
    build_sim_set.py, per map and run:
      simulate_particles        protein particles
      junk.py check / add       CTF check against SIMPLE (≥ 0.99), junk appended, gradients under 30% of the protein
      import_particles, solve2D class averages, halves, every particle's class and alignment
      simulate_particles x 3    noise-free full, helix and strand particles in the protein particles' views
      sim_truth.py              per-class truth with solve2D's alignment, transfer, targets, visibility numbers
```

## The parts

| File | Role |
| --- | --- |
| `sscnn.py` | the one entry point: environment, set, training, prediction, report |
| `entries.tsv` | the entries: EMDB id, PDB id, group, label |
| `build_emdb_set.py` | per entry: download, `emdb_prep.py`, `map_variants.py`, `build_sim_set.py`, clean-up; `--jobs` in parallel |
| `emdb_prep.py` | an EMDB map and its model as input maps |
| `map_variants.py` | low-passed (`lp`) and phase-randomised (`pr`) variants |
| `build_sim_set.py` | runs: particles, junk, solve2D, ground truth (also takes atomic models with `--models`) |
| `junk.py` | junk particles, gradients under protein particles, and the CTF check |
| `sim_truth.py` | per-class truth, transfer, targets, visibility numbers; `load_truth` reads either storage form |
| `models_to_pdb.py` | an atomic model as PDB files for `pdb2mrc` (atomic-model route) |
| `report.py`, `run_page.py` | HTML pages of runs with their truth |
| `realism.py` | simulated against real class averages on the trainer's measurements |
| `emdb_search.py` | finds EMDB entries with half maps and a fitted model, by resolution, box and molecular weight (`--min-kda`, `--max-kda`, `--pages`) |
| `cnn_data.py`, `cnn_ss.py`, `CNN.md` | the class-average CNN: data, model, training, cross-validation, prediction, summary; its description |
| `ss_trainer.py`, `transfer_check.py` | the labelling and selection tool |
| `label_cavgs.py` | a labelling page for real class averages in the CNN's groups (solve2D projects or `cavg_cnn_predict` tables; blind or with the CNN's calls), exported as TSV; uses `../mic_training_generator.py` |
| `simulate_cavgs.py` | shared code (mmCIF reading, secondary-structure assignment, MRC writing); also a fast numpy-only simulator |
| `export_fortran.py` | a trained model as the weights file of SIMPLE's `cavg_cnn_predict` (Fortran inference, no Python at run time); the unit tests' fixture (`testdata/cavg_cnn/`) |

### Input maps (`emdb_prep.py`)

- **Density:** the mean of the two half maps (unsharpened, unmasked).
- **Clean-up:** background subtracted; a soft sphere inside the source box; low-passed to the target Nyquist.
- **Resampling:** trilinear, onto `--box` × `--smpd` (160 px at 2 Å), centred on the model's electron-weighted centre.
- **Helix and strand maps:** the density within 3 Å of the model's helix and strand atoms. Secondary structure comes from C-alpha geometry: helix where d(i, i+3) is 4.6–5.6 Å and d(i, i+4) 5.7–6.8 Å (runs of 6 or more); strand where an extended residue has an extended partner 4.2–5.5 Å away.
- **Frame check:** the correlation at 8 Å between the map and a density computed from the model; entries below `--min-cc` (0.3) are skipped. So are particles wider than `--max-diam` (82% of the box).
- **Recorded:** the share of the density outside the model (micelle, nanodisc).

### Negative variants (`map_variants.py`)

| Variant | Full map | Helix and strand truth |
| --- | --- | --- |
| `lp<L>` | everything coarser than 1.3 L kept, nothing finer than L (one-sided cosine edge) | low-passed the same way |
| `pr<L>` | the same low-pass, plus phase-randomised detail beyond L, confined to the particle's envelope and rescaled shell by shell to the map's own power | low-passed at L |

The training targets keep only detail finer than 12–15 Å, so both variants' targets are zero.

### Junk (`junk.py`)

| Kind | Signal | CTF |
| --- | --- | --- |
| empty | none | — |
| edge | a tanh step at a random angle, offset ±40 Å, 3–10 Å wide, continuing past the box | yes |
| blob | the projection of a dense sphere, radius 20–60 Å | yes |
| aggregate | 6–20 lumps (8–20 Å) along a random walk, with 10 Å texture | yes |
| gradient | an intensity ramp or a smooth step (30–200 Å wide) | no |

**Amplitudes:** edges at 0.5–3 times the protein signal level, blobs and aggregates twice that, gradients at 0.1–0.6 times the noise sd. Gradients under protein particles are 0.05–0.3 noise sd. **CTF:** SIMPLE's own formula, checked in every run against 20 of SIMPLE's noise-free particles (correlation must reach 0.99).

### Ground truth (`sim_truth.py`)

- **Per-class truth:** each class's average of its particles' noise-free full, helix and strand images, moved by solve2D's alignment. Junk counts as zero.
- **Alignment convention:** found from the data (8 combinations tried on the 10 largest classes) and stored in `conv.json`.
- **Transfer:** H, shell by shell, from the class average against its full truth (units, restoration filter, alignment blur). A class with fewer than 10 protein particles, more than half junk, or a correlation below 0.5 with its noise-free average takes the median H of the run's other classes instead. The truth passed through H is what survives in the class average (`*_on_cavg`).
- **Targets:** `helix_band`, `strand_band`, the surviving truth high-passed from 15 Å (zero) to 12 Å (full), up to Nyquist.
- **Visibility numbers** per class: √(Σ_band |H·T|² / (Σ_band |N|² + 10⁻⁴ Σ |C|²)) for helix and strand, in 12–6 Å and 6–4 Å bands. N is the class average's noise, from (even − odd)/2. The floor (a guess) needs detail to reach about 1% of the class average C's overall contrast, so a band that SIMPLE's restoration emptied doesn't inflate the ratio. No threshold is applied.

### Output per run (`runs/<name>/`)

| File | Content |
| --- | --- |
| `cavgs.mrc`, `cavgs_even.mrc`, `cavgs_odd.mrc` | solve2D's class averages and halves (32-bit MRC) |
| truth | `truth_{full,helix,strand}.mrc` (raw), `truth_{full,helix,strand}_on_cavg.mrc`, `truth_{helix,strand}_band.mrc`; with `--compact` all eight as 16-bit arrays in `truth.npz`. Read either with `sim_truth.load_truth(run_dir, key)` |
| `classes.tsv` | per class: `pop`, `njunk`, `res_A`, `cc_full`, `H_own` (1: its own transfer, 0: the run's median), `H_12_6`, `H_6_4`, `cnr_helix_12_6`, `cnr_helix_6_4`, `cnr_strand_12_6`, `cnr_strand_6_4` |
| `cls2D.txt` | SIMPLE's class records |
| `conv.json`, `junk.json` | the alignment convention; the protein count, each junk particle's kind, the protein particles with gradients |
| `oris.txt.gz`, `ptcl2D.txt.gz` | every particle's true view, shift and CTF; its solve2D class and alignment |
| `run.json` | settings, time, CTF check, the truth step's log |

**To look at runs:**

```
python -I report.py --title T --out r.html --run LABEL=runs/<name>[=maps/<id>/_maps_<id>] [--run ...]
python -I run_page.py runs/<name> [...] --out page.html
```

## The class-average CNN (`cnn_data.py`, `cnn_ss.py`)

Sorts class averages into **junk**, **secondary structure present** and **likely protein**; trained on the simulated set. [CNN.md](CNN.md) describes the labels, input, architecture, training, validation and limits.

```
python3 sscnn.py train   --set SET_DIR --model-dir MODEL_DIR    # then open MODEL_DIR/cv_report.html
python3 sscnn.py predict --model-dir MODEL_DIR --stk cavgs_iterNNN.mrc [--cls classdoc_ranked.txt]
python3 sscnn.py report  --model-dir MODEL_DIR --set SET_DIR    # MODEL_DIR/summary.html
```

## The labelling and selection tool (`ss_trainer.py`)

```
python ss_trainer.py --stk cavgs_iterNNN.mrc --mskdiam A --workdir DIR [--cls cls2D.txt] [--port 8765]
```

Open http://127.0.0.1:8765. It needs `<stk>_even.mrc` and `<stk>_odd.mrc` beside the stack.

- **Label:** `y`, `n`, `u`, `0` (clear); `Enter` shows the average and its halves at 384 px.
- **Train:** a logistic model on six measurements (half FRC at 12–6 Å, rod-ridge strength, helix-band power ratio, its directionality, population, resolution). It needs 3 yes and 3 no labels.
- **Export:** `selected.mrc` for `simple_exec prg=map_cavgs_selection projfile=... stk=... stk2=DIR/selected.mrc`.

The half-based measurements need solve2D classes. The ridge measure does not single out helices: it scores the negative variants higher than real helices.

## Limits

- No conformational heterogeneity, beam-induced motion or preferred orientation.
- Junk is drawn from five simple kinds, with no protein-containing aggregates.
- Strand spacing (4.8 Å) is near the 2 Å pixel's Nyquist limit; sheets show mainly as slabs.
- The CNN's secondary-structure label is the run's kind (full map or not), not what shows in each class.
