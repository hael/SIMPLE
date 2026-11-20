#!/usr/bin/env python3
"""A CNN that sorts 2D class averages into junk, secondary structure present, and likely protein (no visible
secondary structure).

Model: a residual CNN on the class average and its even and odd halves (3 channels, 160 px at 2 A; or the average
alone with --channels avg), full resolution in its first layer, five stride-2 stages, global average pooling and two
outputs:
  P(protein)                         trained on every class against its protein fraction (soft target)
  P(secondary structure | protein)   trained on protein classes only (label 1 against label 2)
from which junk = 1 - P(protein), secondary structure = P(protein) P(ss | protein), likely protein = the rest; the
predicted group is the largest of the three. Each group is weighted equally in the loss.
Training: AdamW, linear warm-up then cosine decay over a fixed number of epochs (no epoch is chosen on validation
data), augmentation by the 8 rotations and mirrors of the square, shifts of up to 4 px and swapping the halves, and
two made from each class's own noise (even - odd)/2, its spectrum and radial profile:
  fewer particles   (--aug-noise 0.5 of the classes) independent noise added to each half as if the class had a
                    fraction f of its particles, f log-uniform in [--aug-noise-min 0.05, 1]
  background        (--aug-bg 0.5) noise at the class's own level added outside its mask, out to a random radius up
                    to the whole box, as in unmasked class averages
then all three channels scaled again as cnn_data.normalise does.
Evaluation averages the logits over the 8 rotations and mirrors.

Subcommands:
  train    cross-validation over --folds folds of whole entries (proteins), so a protein is never in both training
           and validation; then a final model on all classes with the same settings. Writes to --out: model.pt (the
           final model), cv_predictions.tsv (every class predicted by the fold that did not train on it), cv_report.txt
           and cv_report.html (metrics, per-fold and per-entry spread, breakdown by run kind and junk kind, galleries of
           the errors), cv_curves.tsv (training and validation loss per epoch and fold), train.log.
           --folds 0 trains the final model only; --no-final runs the cross-validation only; --split run makes folds of
           runs instead of entries (only for testing the mechanics on one entry: it leaks the protein's shape).
  predict  a model on a class-average stack (resampled to the model's pixel size and box), optionally against labels
           from ss_trainer.py (labels.json: yes = secondary structure visible, no = not) and with SIMPLE's class document
           (--cls: each class's population and resolution, and the particles in the classes kept); writes a TSV, an
           HTML page and a JSON record of the prediction (beside the TSV).
  summary  one page over a whole run: the simulated set (entries, their status, runs per kind), the training data,
           the cross-validation and every prediction in --pred-dir, with links to the detailed pages.

Usage:
  python -I cnn_ss.py train --data DATA_DIR --out RUN_DIR [--folds 5] [--epochs 40] [--batch 32] [--lr 1e-3]
                            [--width 32] [--channels avg+halves] [--device auto] [--nthr 0]
  python -I cnn_ss.py predict --model RUN_DIR/model.pt --stk cavgs.mrc [--labels labels.json] [--cls classdoc.txt]
                              --out preds.tsv [--html page.html]
  python -I cnn_ss.py summary --model-dir RUN_DIR [--set SET_DIR] [--pred-dir DIR] [--out summary.html]
DATA_DIR comes from cnn_data.py. Needs numpy and PyTorch (python3 sscnn.py setup makes them).
"""
import argparse
import base64
import html
import json
import math
import os
import re
import sys
import time

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ss_trainer as T   # noqa: E402
from cnn_data import GROUPS, NAME, normalise, resample   # noqa: E402

SHORT = ('junk', 'ss', 'likely')
SHIFT = 4
SURFACE, TEXT1, TEXT2, HAIR = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'


# ------------------------------------------------------------------ model
def gn(c):
    return nn.GroupNorm(8, c)


class Block(nn.Module):
    def __init__(self, cin, cout, stride):
        super().__init__()
        self.c1, self.n1 = nn.Conv2d(cin, cout, 3, stride, 1, bias=False), gn(cout)
        self.c2, self.n2 = nn.Conv2d(cout, cout, 3, 1, 1, bias=False), gn(cout)
        self.skip = None if cin == cout and stride == 1 else nn.Sequential(nn.Conv2d(cin, cout, 1, stride, bias=False), gn(cout))

    def forward(self, x):
        y = self.n2(self.c2(F.relu(self.n1(self.c1(x)))))
        return F.relu(y + (x if self.skip is None else self.skip(x)))


class Net(nn.Module):
    def __init__(self, cin=3, width=32, dropout=0.3):
        super().__init__()
        w = width
        self.stem = nn.Sequential(nn.Conv2d(cin, w, 3, 1, 1, bias=False), gn(w), nn.ReLU())   # full resolution
        ch = [w, 2 * w, 4 * w, 4 * w, 8 * w, 8 * w]
        self.body = nn.Sequential(*[Block(ch[i], ch[i + 1], 2) for i in range(5)])          # 160 -> 5 px
        self.drop = nn.Dropout(dropout)
        self.head = nn.Linear(ch[-1], 2)       # logits: protein; secondary structure given protein

    def forward(self, x):
        return self.head(self.drop(self.body(self.stem(x)).mean((2, 3))))


def group_probs(lg):
    pp, ps = 1 / (1 + np.exp(-lg[:, 0])), 1 / (1 + np.exp(-lg[:, 1]))
    return np.stack([1 - pp, pp * ps, pp * (1 - ps)], 1)



# ------------------------------------------------------------------ data and training
def load_data(d):
    x = np.load(os.path.join(d, 'x.npy'), mmap_mode='r')
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(d, 'meta.tsv'))]
    meta = [dict(zip(rows[0], r)) for r in rows[1:]]
    return x, meta, json.load(open(os.path.join(d, 'info.json')))


def pick_device(want):
    if want != 'auto':
        return torch.device(want)
    if torch.cuda.is_available():
        return torch.device('cuda')
    if getattr(torch.backends, 'mps', None) and torch.backends.mps.is_available():
        return torch.device('mps')
    return torch.device('cpu')


_RINGS = {}


def rings(n):
    """Integer radius of every pixel of an n x n image, and the pixel count per radius."""
    if n not in _RINGS:
        g = torch.arange(n) - n // 2
        ri = torch.sqrt(g[:, None].float() ** 2 + g[None, :].float() ** 2).round().long()
        _RINGS[n] = (ri, torch.bincount(ri.flatten()).float(), (ri < n / 2))
    return _RINGS[n]


def unit_noise(D, gen):
    """A noise image with D's power spectrum (random phases), of unit standard deviation."""
    mag = torch.fft.rfft2(D).abs()
    ph = torch.rand(mag.shape, generator=gen) * (2 * math.pi)
    N = torch.fft.irfft2(torch.polar(mag, ph), s=D.shape)
    return N / (N.std() + 1e-12)


def noise_aug(v, gen, a):
    """The fewer-particles and background augmentations of one class (3 channels: average, even, odd)."""
    n = v.shape[-1]
    ri, cnt, inner = rings(n)
    avg, ev, od = v[0].clone(), v[1].clone(), v[2].clone()
    D = 0.5 * (ev - od)                                     # the noise of the average
    rms = torch.sqrt(torch.bincount(ri.flatten(), weights=(D ** 2).flatten(), minlength=len(cnt)) / cnt.clamp(min=1))
    if rms.max() <= 0:
        return v
    if torch.rand(1, generator=gen) < a.aug_noise:
        lo = math.log(a.aug_noise_min)
        f = math.exp(lo + (0 - lo) * float(torch.rand(1, generator=gen)))
        amp = math.sqrt(2 * (1 / f - 1))                    # per half: the average's noise variance grows by 1/f
        env = rms[ri]
        add_e, add_o = amp * unit_noise(D, gen) * env, amp * unit_noise(D, gen) * env
        ev, od, avg = ev + add_e, od + add_o, avg + 0.5 * (add_e + add_o)
    if torch.rand(1, generator=gen) < a.aug_bg:
        level = rms[2:n // 4].median()
        r = torch.arange(len(rms))
        out = torch.nonzero((rms < 0.25 * level) & (r > n // 8))
        if len(out):                                        # a masked class: where its noise stops
            r_mask = float(out[0])
            r_out = r_mask + 8 + (0.72 * n - r_mask - 8) * float(torch.rand(1, generator=gen))
            rf = ri.float()
            w = ((rf - r_mask) / 4).clamp(0, 1) * ((r_out - rf) / 4).clamp(0, 1)
            add_e, add_o = level * unit_noise(D, gen) * w, level * unit_noise(D, gen) * w
            ev, od, avg = ev + add_e, od + add_o, avg + 0.5 * (add_e + add_o)
    m, s = avg[inner].mean(), avg[inner].std()
    s = s if s > 0 else 1.0
    return torch.stack([(avg - m) / s, (ev - m) / s, (od - m) / s])


def augment(xb, gen, a):
    out = torch.empty_like(xb)
    for i in range(len(xb)):
        v = xb[i]
        if v.shape[0] == 3 and (a.aug_noise > 0 or a.aug_bg > 0):
            v = noise_aug(v, gen, a)
        v = torch.rot90(v, int(torch.randint(4, (1,), generator=gen)), (1, 2))
        if torch.rand(1, generator=gen) < 0.5:
            v = torch.flip(v, (2,))
        if v.shape[0] == 3 and torch.rand(1, generator=gen) < 0.5:
            v = v[[0, 2, 1]]                                     # the halves are interchangeable
        dy, dx = torch.randint(-SHIFT, SHIFT + 1, (2,), generator=gen).tolist()
        out[i] = torch.roll(v, (dy, dx), (1, 2))
    return out


def d4(x):
    for k in range(4):
        r = torch.rot90(x, k, (2, 3))
        yield r
        yield torch.flip(r, (3,))


def autocast(device):
    if device.type == 'cuda' and torch.cuda.is_bf16_supported():
        return torch.autocast('cuda', dtype=torch.bfloat16)
    return torch.autocast(device.type, enabled=False)


@torch.no_grad()
def logits_of(model, X, idx, device, bs, tta, nch=None):
    model.eval()
    out = []
    for i in range(0, len(idx), bs):
        xb = X[idx[i:i + bs]].to(device).float()[:, :nch]
        with autocast(device):
            lg = torch.stack([model(v) for v in d4(xb)]).mean(0) if tta else model(xb)
        out.append(lg.float().cpu())
    return torch.cat(out).numpy() if out else np.zeros((0, 2))


def loss_fn(lg, pf, lab, w):
    lp = F.binary_cross_entropy_with_logits(lg[:, 0], pf, reduction='none')
    ls = F.binary_cross_entropy_with_logits(lg[:, 1], (lab == 1).float(), reduction='none') * (lab > 0).float()
    return ((lp + ls) * w).sum() / w.sum()


def eval_loss(lg, pf, lab, w):
    t = lambda v: torch.as_tensor(np.asarray(v), dtype=torch.float32)   # noqa: E731
    return float(loss_fn(t(lg), t(pf), t(lab), t(w))) if len(lg) else float('nan')


def group_weights(lab, idx):
    n = np.array([(lab[idx] == k).sum() for k in range(3)], float)
    wk = np.where(n > 0, len(idx) / (3 * np.maximum(n, 1)), 0.0)
    return wk[lab]


def train_model(X, lab, pf, tr, va, a, device, log, tag):
    """One model trained on tr for a.epochs; validation loss per epoch on va (no augmentation, no averaging)."""
    torch.manual_seed(a.seed)
    gen = torch.Generator().manual_seed(a.seed)
    nch = 1 if a.channels == 'avg' else 3
    model = Net(nch, a.width, a.dropout).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=a.wd)
    spe = math.ceil(len(tr) / a.batch)
    total, warm = a.epochs * spe, max(1, spe)
    sched = torch.optim.lr_scheduler.LambdaLR(
        opt, lambda s: (s + 1) / warm if s < warm else 0.5 * (1 + math.cos(math.pi * (s - warm) / max(1, total - warm))))
    w_all = group_weights(lab, tr)
    labt, pft, wt = torch.as_tensor(lab), torch.as_tensor(pf, dtype=torch.float32), torch.as_tensor(w_all, dtype=torch.float32)
    wva = group_weights(lab, va) if len(va) else None
    curve = []
    for ep in range(a.epochs):
        model.train()
        t0 = time.time()
        perm = tr[torch.randperm(len(tr), generator=gen).numpy()]
        tot, nb = 0.0, 0
        for i in range(0, len(perm), a.batch):
            b = perm[i:i + a.batch]
            if len(b) < 2:
                continue
            xb = augment(X[b].float(), gen, a)[:, :nch].to(device)    # augmented on all 3 channels, then cut
            with autocast(device):
                lg = model(xb)
            loss = loss_fn(lg.float(), pft[b].to(device), labt[b].to(device), wt[b].to(device))
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            sched.step()
            tot, nb = tot + loss.item(), nb + 1
        vl = float('nan')
        if len(va) and ((ep + 1) % a.eval_every == 0 or ep + 1 == a.epochs):
            vl = eval_loss(logits_of(model, X, va, device, a.batch * 2, False, nch), pf[va], lab[va], wva[va])
        curve.append((tag, ep + 1, tot / max(nb, 1), vl))
        log(f'{tag} epoch {ep + 1:3d}/{a.epochs}: train loss {tot / max(nb, 1):.4f}, validation loss {vl:.4f}, '
            f'{time.time() - t0:.1f} s')
    return model, curve, nch


# ------------------------------------------------------------------ metrics
def auc(pos, neg):
    return T.auc(pos, neg)


def metrics(P, lab):
    pred = P.argmax(1)
    conf = np.array([[int(((lab == i) & (pred == j)).sum()) for j in range(3)] for i in range(3)])
    rec = [conf[i, i] / conf[i].sum() if conf[i].sum() else float('nan') for i in range(3)]
    prec = [conf[i, i] / conf[:, i].sum() if conf[:, i].sum() else float('nan') for i in range(3)]
    pp = P[:, 1] + P[:, 2]
    ps = np.divide(P[:, 1], pp, out=np.zeros_like(pp), where=pp > 0)
    return {'n': int(len(lab)), 'confusion': conf.tolist(), 'recall': rec, 'precision': prec,
            'balanced_accuracy': float(np.nanmean(rec)), 'accuracy': float((pred == lab).mean()) if len(lab) else float('nan'),
            'auc_protein_vs_junk': auc(pp[lab > 0], pp[lab == 0]),
            'auc_ss_vs_likely': auc(ps[lab == 1], ps[lab == 2]),
            'log_loss': float(-np.mean(np.log(np.clip(P[np.arange(len(lab)), lab], 1e-7, 1)))) if len(lab) else float('nan')}


SNR_BANDS = ('< 0.003', '0.003-0.01', '>= 0.01')


def snr_band(m):
    s = float(m.get('snr_px') or m['snr'])         # the noise per pixel as a 150 A particle's SNR (older data: snr)
    return SNR_BANDS[0] if s < 0.003 else (SNR_BANDS[1] if s < 0.01 else SNR_BANDS[2])


def stratum(m):
    if int(m['label']) == 0:
        return f'junk class, mostly {m["junk_kind"] or "?"}' if m['junk_kind'] else 'junk class'
    kind = 'full map' if m['kind'] == 'map' else f'{m["kind"]}{m["cutoff"]}'
    return f'{kind}, SNR {snr_band(m)}'


def ss_auc_by_snr(P, lab, meta):
    """AUC of P(ss | protein), secondary structure against likely protein, within each particle-SNR band: well below
    the overall AUC means the network reads the noise level rather than the secondary structure."""
    pp = P[:, 1] + P[:, 2]
    ps = np.divide(P[:, 1], pp, out=np.zeros_like(pp), where=pp > 0)
    band = np.array([snr_band(m) for m in meta])
    out = []
    for b in SNR_BANDS:
        s, l = ps[(band == b) & (lab == 1)], ps[(band == b) & (lab == 2)]
        out.append((b, len(s), len(l), auc(s, l)))
    return out


def fmt(v, nd=3):
    return 'n/a' if v is None or (isinstance(v, float) and not np.isfinite(v)) else f'{v:.{nd}f}'


def report_lines(mt, title):
    c = np.array(mt['confusion'])
    out = [title, f'  classes {mt["n"]}; balanced accuracy {fmt(mt["balanced_accuracy"])}, accuracy {fmt(mt["accuracy"])}, '
                  f'log loss {fmt(mt["log_loss"])}',
           f'  AUC protein vs junk {fmt(mt["auc_protein_vs_junk"])}; AUC secondary structure vs likely protein (protein '
           f'classes) {fmt(mt["auc_ss_vs_likely"])}',
           '  confusion (rows: label, columns: predicted)   ' + '  '.join(f'{s:>7s}' for s in SHORT) + '   recall',
           ] + [f'  {GROUPS[i]:43s} ' + '  '.join(f'{c[i, j]:7d}' for j in range(3)) + f'   {fmt(mt["recall"][i], 2)}'
                for i in range(3)] + [
           '  precision' + ' ' * 35 + '  '.join(f'{fmt(p, 2):>7s}' for p in mt['precision'])]
    return out


# ------------------------------------------------------------------ reports
def thumb(img, size=96):
    return 'data:image/png;base64,' + base64.b64encode(T.png(np.asarray(img, np.float64), size)).decode()


VIEW_BAR = ('<div class="bar"><label>Brightness <input type="range" id="br" min="0.2" max="3" step="0.05" value="1">'
            '<span id="brv">1.00</span></label><label>Contrast <input type="range" id="ct" min="0.2" max="4" step="0.05" '
            'value="1"><span id="ctv">1.00</span></label><button id="rs">Reset</button></div>'
            '<script>(function(){var r=document.documentElement.style;'
            'function u(k,v){r.setProperty("--"+k,v);document.getElementById(k+"v").textContent=(+v).toFixed(2);}'
            '["br","ct"].forEach(function(k){var e=document.getElementById(k);e.oninput=function(){u(k,e.value);};});'
            'document.getElementById("rs").onclick=function(){["br","ct"].forEach(function(k){'
            'document.getElementById(k).value=1;u(k,1);});};})();</script>')


def page(title, body):
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(title)}</title><style>'
            f'body{{background:{SURFACE};color:{TEXT1};font:14px/1.45 system-ui,sans-serif;margin:24px;max-width:1200px}}'
            f'h1{{font-size:20px}}h2{{font-size:16px;margin-top:28px}}p,.cap{{color:{TEXT2}}}'
            f'table{{border-collapse:collapse;margin:8px 0}}td,th{{border-bottom:1px solid {HAIR};padding:3px 10px;'
            f'text-align:right;font-variant-numeric:tabular-nums}}th{{color:{TEXT2};font-weight:600}}'
            f'td:first-child,th:first-child{{text-align:left}}.g{{display:flex;flex-wrap:wrap;gap:10px}}'
            f'.t{{width:96px}}.cap{{font-size:11px;line-height:1.25;margin-top:2px}}pre{{font-size:12px}}'
            f'img.cv{{filter:brightness(var(--br,1)) contrast(var(--ct,1))}}'
            f'.bar{{position:sticky;top:0;z-index:2;background:{SURFACE};border-bottom:1px solid {HAIR};padding:8px 0;'
            f'display:flex;gap:24px;align-items:center;color:{TEXT2}}}.bar label{{display:flex;gap:8px;align-items:center}}'
            f'.bar span{{min-width:2.6em;font-variant-numeric:tabular-nums;color:{TEXT1}}}'
            f'</style></head><body><h1>{html.escape(title)}</h1>{VIEW_BAR}{body}</body></html>')


def table(header, rows):
    return ('<table><tr>' + ''.join(f'<th>{html.escape(h)}</th>' for h in header) + '</tr>'
            + ''.join('<tr>' + ''.join(f'<td>{html.escape(str(v))}</td>' for v in r) + '</tr>' for r in rows) + '</table>')


def tile(img, cap, tip):
    return (f'<div class="t"><img class="cv" src="{thumb(img)}" width="96" height="96" title="{html.escape(tip)}" '
            f'style="display:block;border-radius:3px"><div class="cap">{cap}</div></div>')


def cv_html(path, title, intro, mt, folds_rows, entry_rows, strata_rows, X, meta, P, lab, ngal, snr_line=''):
    pred = P.argmax(1)
    c = np.array(mt['confusion'])
    body = f'<p>{html.escape(intro)}</p>'
    body += table(['label \\ predicted'] + list(GROUPS) + ['recall'],
                  [[GROUPS[i]] + [int(c[i, j]) for j in range(3)] + [fmt(mt['recall'][i], 2)] for i in range(3)]
                  + [['precision'] + [fmt(p, 2) for p in mt['precision']] + ['']])
    body += (f'<p>Balanced accuracy {fmt(mt["balanced_accuracy"])}; AUC protein vs junk {fmt(mt["auc_protein_vs_junk"])}; '
             f'AUC secondary structure vs likely protein {fmt(mt["auc_ss_vs_likely"])}.</p>'
             + (f'<p>{html.escape(snr_line)}</p>' if snr_line else ''))
    body += '<h2>Per fold</h2>' + table(['fold', 'entries', 'classes', 'balanced accuracy', 'AUC protein', 'AUC ss'], folds_rows)
    body += '<h2>Per entry</h2>' + table(['entry', 'classes', 'balanced accuracy', 'AUC protein', 'AUC ss'], entry_rows)
    body += '<h2>By kind of class</h2>' + table(['class kind', 'classes', '% junk', '% ss', '% likely'], strata_rows)
    body += ('<h2>Errors</h2><p>Per label and wrong prediction, the classes predicted most confidently. Caption: run, class; '
             'the wrong group\'s probability; helix visibility (12-6 A) or the junk kind; particles (junk).</p>')
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            sel = np.where((lab == i) & (pred == j))[0]
            if not len(sel):
                continue
            sel = sel[np.argsort(-P[sel, j])][:ngal]
            tiles = []
            for k in sel:
                m = meta[k]
                extra = f'vis {float(m["cnr_helix_12_6"]):.2f}' if i > 0 else (m['junk_kind'] or '')
                cap = f'{html.escape(m["run"])} c{m["cls"]}<br>p {P[k, j]:.2f} · {html.escape(extra)}<br>{m["pop"]} ({m["njunk"]})'
                tip = (f'{m["run"]} class {m["cls"]}: label {GROUPS[i]}, predicted {GROUPS[j]}; P = '
                       + ', '.join(f'{GROUPS[g]} {P[k, g]:.2f}' for g in range(3))
                       + f'; snr {m["snr"]}; res_A {m["res_A"]}; helix vis {m["cnr_helix_12_6"]}, strand vis {m["cnr_strand_12_6"]}')
                tiles.append(tile(np.asarray(X[k, 0], np.float64), cap, tip))
            body += f'<h2>{GROUPS[i]} predicted {GROUPS[j]} ({int(((lab == i) & (pred == j)).sum())})</h2><div class="g">{"".join(tiles)}</div>'
    open(path, 'w').write(page(title, body))


def read_class_doc(path, n):
    """Population and resolution per class from a SIMPLE class document (class=, pop=, res= per line, e.g.
    classdoc_ranked.txt or a printed cls2D field); a line without class= is the next class."""
    pop, res = np.full(n, np.nan), np.full(n, np.nan)
    for i, line in enumerate(open(path)):
        v = dict(re.findall(r'(\w+)=(\S+)', line))
        if 'pop' not in v and 'res' not in v:
            continue                                         # blank or other lines
        c = int(float(v['class'])) if 'class' in v else i + 1
        if 1 <= c <= n:
            pop[c - 1] = float(v.get('pop', 'nan'))
            res[c - 1] = float(v.get('res', 'nan'))
    return pop, res


# ------------------------------------------------------------------ subcommands
def plain(v):
    return json.loads(json.dumps(v, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))


def cmd_train(a):
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, 'train.log'), 'a')

    def log(s):
        print(s)
        logf.write(s + '\n')
        logf.flush()
        sys.stdout.flush()
    if a.nthr > 0:
        torch.set_num_threads(a.nthr)
    device = pick_device(a.device)
    xm, meta, info = load_data(a.data)
    nch = 1 if a.channels == 'avg' else 3
    X = torch.from_numpy(np.array(xm))                     # float16, all 3 channels in memory (augmentation needs the halves)
    lab = np.array([int(m['label']) for m in meta])
    pf = np.array([float(m['protein_frac']) for m in meta], np.float32)
    log(f'{time.strftime("%Y-%m-%d %H:%M")} train: {len(meta)} classes ({", ".join(f"{(lab == k).sum()} {g}" for k, g in enumerate(GROUPS))}), '
        f'{len(info["entries"])} entries; device {device}, torch {torch.__version__}; {vars(a)}')
    cfg = {'width': a.width, 'dropout': a.dropout, 'channels': a.channels, 'aug_noise': a.aug_noise,
           'aug_noise_min': a.aug_noise_min, 'aug_bg': a.aug_bg, 'box': info['box'], 'smpd': info['smpd'],
           'groups': list(GROUPS), 'epochs': a.epochs, 'batch': a.batch, 'lr': a.lr, 'wd': a.wd, 'seed': a.seed,
           'data': os.path.abspath(a.data), 'data_counts': info['counts'], 'entries': info['entries'],
           'label_rule': info['label_rule']}
    summary = {}
    if a.folds:
        key = 'entry' if a.split == 'entry' else 'run'
        groups = sorted({m[key] for m in meta})
        if len(groups) < a.folds:
            raise SystemExit(f'{len(groups)} {key}s for {a.folds} folds: use fewer folds' +
                             (' or --split run to test the mechanics' if key == 'entry' else ''))
        order = np.random.default_rng(a.seed).permutation(len(groups))
        fold_of = {groups[g]: i % a.folds for i, g in enumerate(order)}
        fold = np.array([fold_of[m[key]] for m in meta])
        if key == 'run':
            log('WARNING: folds of runs, not entries; the same protein is in training and validation, so these numbers '
                'only test the mechanics')
        P = np.zeros((len(meta), 3))
        curves, folds_rows, fold_metrics = [], [], []
        for f in range(a.folds):
            tr, va = np.where(fold != f)[0], np.where(fold == f)[0]
            model, curve, _ = train_model(X, lab, pf, tr, va, a, device, log, f'fold {f + 1}/{a.folds}')
            curves += curve
            P[va] = group_probs(logits_of(model, X, va, device, a.batch * 2, True, nch))
            mt_tr = metrics(group_probs(logits_of(model, X, tr, device, a.batch * 2, False, nch)), lab[tr])
            mt = metrics(P[va], lab[va])
            fold_metrics.append({'fold': f + 1, 'groups': sorted({meta[i][key] for i in va}), 'validation': mt,
                                 'training_balanced_accuracy': mt_tr['balanced_accuracy']})
            folds_rows.append([f + 1, ' '.join(sorted({meta[i][key] for i in va})), len(va), fmt(mt['balanced_accuracy']),
                               fmt(mt['auc_protein_vs_junk']), fmt(mt['auc_ss_vs_likely'])])
            log(f'fold {f + 1}: validation balanced accuracy {fmt(mt["balanced_accuracy"])}, training (no augmentation) '
                f'{fmt(mt_tr["balanced_accuracy"])}')
        mt = metrics(P, lab)
        entry_rows = []
        for e in sorted({m['entry'] for m in meta}):
            s = np.array([m['entry'] == e for m in meta])
            me = metrics(P[s], lab[s])
            entry_rows.append([e, int(s.sum()), fmt(me['balanced_accuracy']), fmt(me['auc_protein_vs_junk']), fmt(me['auc_ss_vs_likely'])])
        strata = {}
        for k, m in enumerate(meta):
            strata.setdefault(stratum(m), []).append(k)
        pred = P.argmax(1)
        strata_rows = [[s, len(ix)] + [f'{100 * np.mean(pred[ix] == g):.0f}' for g in range(3)] for s, ix in sorted(strata.items())]
        with open(os.path.join(a.out, 'cv_predictions.tsv'), 'w') as f:
            cols = list(meta[0].keys())
            f.write('\t'.join(cols + ['fold', 'p_junk', 'p_ss', 'p_likely', 'predicted']) + '\n')
            for k, m in enumerate(meta):
                f.write('\t'.join([m[c] for c in cols] + [str(fold[k] + 1)] + [f'{v:.4f}' for v in P[k]] + [GROUPS[pred[k]]]) + '\n')
        with open(os.path.join(a.out, 'cv_curves.tsv'), 'w') as f:
            f.write('fold\tepoch\ttrain_loss\tvalidation_loss\n')
            for tag, ep, tl, vl in curves:
                f.write(f'{tag.split()[1].split("/")[0]}\t{ep}\t{tl:.5f}\t{vl:.5f}\n')
        by_snr = ss_auc_by_snr(P, lab, meta)
        snr_line = ('AUC secondary structure vs likely protein within particle-SNR bands (ss / likely classes): '
                    + '; '.join(f'{b}: {fmt(v)} ({ns}/{nl})' for b, ns, nl, v in by_snr)
                    + f'; overall {fmt(mt["auc_ss_vs_likely"])}. If the bands are well below the overall AUC, the network reads '
                    + 'the noise level rather than the secondary structure.')
        intro = (f'{a.folds}-fold cross-validation by {key}; every class predicted by the model that did not train on it. '
                 f'{len(meta)} classes from {len(info["entries"])} entries; {a.epochs} epochs, width {a.width}, channels '
                 f'{a.channels}. Labels: {info["label_rule"]}.')
        lines = [intro, ''] + report_lines(mt, 'All classes (out of fold):') + ['  ' + snr_line, '', 'Per fold:'] + \
            ['  ' + '  '.join(str(v) for v in r) for r in folds_rows] + ['', 'Per entry:'] + \
            ['  ' + '  '.join(str(v) for v in r) for r in entry_rows] + ['', 'By kind of class (% predicted junk / ss / likely):'] + \
            [f'  {r[0]:34s} {r[1]:5d}  ' + '  '.join(f'{v:>4s}' for v in r[2:]) for r in strata_rows]
        open(os.path.join(a.out, 'cv_report.txt'), 'w').write('\n'.join(lines) + '\n')
        cv_html(os.path.join(a.out, 'cv_report.html'), 'Class-average CNN: cross-validation', intro, mt, folds_rows, entry_rows,
                strata_rows, xm, meta, P, lab, a.gallery, snr_line)
        log('\n'.join(lines))
        summary = {'folds': a.folds, 'split': key, 'out_of_fold': mt, 'per_fold': fold_metrics,
                   'ss_auc_by_snr': [{'snr': b, 'n_ss': ns, 'n_likely': nl, 'auc': v} for b, ns, nl, v in by_snr]}
        json.dump(summary, open(os.path.join(a.out, 'cv_metrics.json'), 'w'), indent=1)
    if not a.no_final:
        model, _, _ = train_model(X, lab, pf, np.arange(len(meta)), np.array([], int), a, device, log, 'final')
        cfg['cv'] = {k: summary[k] for k in ('folds', 'split', 'out_of_fold')} if summary else None
        cfg = plain(cfg)                     # plain Python values only, so that torch.load(weights_only=True) reads it
        torch.save({'state': {k: v.cpu() for k, v in model.state_dict().items()}, 'config': cfg}, os.path.join(a.out, 'model.pt'))
        log(f'final model: {os.path.join(a.out, "model.pt")}')


def cmd_predict(a):
    if a.nthr > 0:
        torch.set_num_threads(a.nthr)
    device = pick_device(a.device)
    ck = torch.load(a.model, map_location='cpu', weights_only=True)
    cfg = ck['config']
    nch = 1 if cfg['channels'] == 'avg' else 3
    model = Net(nch, cfg['width'], cfg['dropout'])
    model.load_state_dict(ck['state'])
    model.to(device)
    avg, smpd, _, _ = T.read_mrc(a.stk)
    smpd = a.smpd or smpd
    stem = a.stk[:-4] if a.stk.endswith('.mrc') else a.stk
    ev = T.read_mrc(a.even or stem + '_even.mrc')[0] if nch == 3 else avg
    od = T.read_mrc(a.odd or stem + '_odd.mrc')[0] if nch == 3 else avg
    xs, empty = [], []
    for i in range(len(avg)):
        ch = [np.asarray(im[i], np.float64) for im in (avg, ev, od)]
        empty.append(bool(ch[0].std() == 0))
        if abs(smpd - cfg['smpd']) > 1e-3 or ch[0].shape[-1] != cfg['box']:
            ch = [resample(im, smpd, cfg['smpd'], cfg['box']) for im in ch]
        xs.append(normalise(*ch)[:nch].astype(np.float32))
    X = torch.from_numpy(np.stack(xs))
    P = group_probs(logits_of(model, X, np.arange(len(xs)), device, 32, True))
    pred = P.argmax(1)
    labels = json.load(open(a.labels)) if a.labels else {}
    pop, res = read_class_doc(a.cls, len(P)) if a.cls else (np.full(len(P), np.nan), np.full(len(P), np.nan))
    nf = lambda v, d=0: '' if not np.isfinite(v) else f'{v:.{d}f}'   # noqa: E731
    with open(a.out, 'w') as f:
        f.write('class\tp_junk\tp_ss\tp_likely\tpredicted\tlabel\tpop\tres_A\n')
        for i in range(len(P)):
            g = 'empty' if empty[i] else GROUPS[pred[i]]
            f.write(f'{i + 1}\t' + '\t'.join(f'{v:.4f}' for v in P[i]) + f'\t{g}\t{labels.get(str(i + 1), "")}\t{nf(pop[i])}\t{nf(res[i], 1)}\n')
    emp = np.array(empty)
    counts = {GROUPS[g]: int(((pred == g) & ~emp).sum()) for g in range(3)}
    counts['empty'] = int(emp.sum())
    kept = float(np.nansum(pop[(pred > 0) & ~emp])) if a.cls else None
    total = float(np.nansum(pop)) if a.cls else None
    print(f'{len(P)} classes at {smpd:.3f} A (model {cfg["smpd"]} A, box {cfg["box"]}): ' + ', '.join(f'{v} {k}' for k, v in counts.items())
          + (f'; particles in protein classes {kept:.0f} of {total:.0f}' if a.cls else '') + f'; written to {a.out}')
    json.dump({'stk': os.path.abspath(a.stk), 'model': os.path.abspath(a.model), 'trained_on': cfg['entries'], 'classes': len(P),
               'counts': counts, 'particles_kept': kept, 'particles_total': total,
               'html': os.path.abspath(a.html) if a.html else None, 'tsv': os.path.abspath(a.out),
               'date': time.strftime('%Y-%m-%d %H:%M')},
              open(os.path.splitext(a.out)[0] + '.json', 'w'), indent=1)
    body = f'<p>{html.escape(a.stk)}, {len(P)} classes; model {html.escape(a.model)} (trained on {html.escape(", ".join(cfg["entries"]))}).</p>'
    if a.cls:
        body += (f'<p>Particles in the classes called protein (secondary structure or likely protein): {kept:.0f} of {total:.0f} '
                 f'({100 * kept / max(total, 1):.0f}%). Captions: class, probability of its group, particles, resolution.</p>')
    if labels:
        yes = [i for i in range(len(P)) if labels.get(str(i + 1)) == 'yes' and not empty[i]]
        no = [i for i in range(len(P)) if labels.get(str(i + 1)) == 'no' and not empty[i]]
        a_ss = auc(P[yes, 1], P[no, 1])
        rows = [[lab_name, len(ix)] + [int((pred[ix] == g).sum()) for g in range(3)]
                for lab_name, ix in (('labelled yes', yes), ('labelled no', no))]
        print(f'labels: {len(yes)} yes, {len(no)} no; AUC of P(secondary structure), yes against no: {fmt(a_ss)}')
        for r in rows:
            print(f'  {r[0]:13s} {r[1]:4d}: predicted junk {r[2]}, ss {r[3]}, likely {r[4]}')
        body += f'<p>Your labels: AUC of P(secondary structure), yes against no, {fmt(a_ss)}.</p>'
        body += table(['label', 'classes'] + [f'predicted {g}' for g in GROUPS], rows)
    if a.html:
        for g in range(3):
            ix = [i for i in np.argsort(-P[:, g]) if pred[i] == g and not empty[i]]
            tiles = []
            for i in ix:
                lab_txt = labels.get(str(i + 1), '')
                cap = f'c{i + 1} · p {P[i, g]:.2f}' + (f'<br>{nf(pop[i])} ptcls · {nf(res[i], 1)} Å' if np.isfinite(pop[i]) else '') \
                    + (f'<br>label {html.escape(lab_txt)}' if lab_txt else '')
                tip = f'class {i + 1}: ' + ', '.join(f'{GROUPS[k]} {P[i, k]:.2f}' for k in range(3))
                tiles.append(tile(avg[i], cap, tip))
            body += f'<h2>{GROUPS[g]} ({len(ix)})</h2><div class="g">{"".join(tiles)}</div>'
        open(a.html, 'w').write(page('Class-average CNN: predictions', body))
        print(f'page: {a.html}')


def cmd_summary(a):
    """One page over a run: the set, the training data, the model, the cross-validation, the predictions."""
    md = os.path.abspath(a.model_dir)
    out = os.path.abspath(a.out or os.path.join(md, 'summary.html'))
    pred_dir = os.path.abspath(a.pred_dir or os.path.join(md, 'predictions'))
    link = lambda p, t: f'<a href="{html.escape(os.path.relpath(p, os.path.dirname(out)))}">{html.escape(t)}</a>'   # noqa: E731
    body = (f'<p>{time.strftime("%Y-%m-%d %H:%M")}. Model folder {html.escape(md)}'
            + (f'; simulated set {html.escape(os.path.abspath(a.set))}' if a.set else '') + '.</p>')
    # the simulated set
    if a.set and os.path.exists(os.path.join(a.set, 'entries.tsv')):
        last = {}
        for line in open(os.path.join(a.set, 'entries.tsv')):
            f = line.rstrip('\n').split('\t')
            if f and f[0].startswith('EMD-'):
                last[f[0]] = f                                # the latest round's line
        runs = {}
        rdir = os.path.join(a.set, 'runs')
        for d in os.listdir(rdir) if os.path.isdir(rdir) else []:
            m = NAME.match(d)
            if m and os.path.exists(os.path.join(rdir, d, 'run.json')):
                runs.setdefault(m.group('entry'), []).append('map' if not m.group('kind') else m.group('kind') + m.group('cut'))
        st = [f[4].split(':')[0].split(' in ')[0] for f in last.values()]
        body += (f'<h2>Simulated set</h2><p>{len(last)} entries: ' + ', '.join(f'{st.count(x)} {x}' for x in sorted(set(st)))
                 + f'; {sum(len(v) for v in runs.values())} runs.</p>')
        rows = []
        for e, f in last.items():
            k = runs.get(e, [])
            lp = '/'.join(sorted((x[2:] for x in k if x.startswith('lp')), key=int)) or '-'
            pr = '/'.join(sorted((x[2:] for x in k if x.startswith('pr')), key=int)) or '-'
            g = lambda i: f[i] if len(f) > i else ''   # noqa: E731
            rows.append([e, g(3), g(2), g(9), g(5), g(6), k.count('map'), lp, pr, g(4)[:70]])
        body += table(['entry', 'protein', 'group', 'res (Å)', 'diameter (Å)', 'map-model cc', 'map runs', 'lp', 'pr', 'status'], rows)
    # training data and model
    info_p = os.path.join(md, 'data', 'info.json')
    if os.path.exists(info_p):
        info = json.load(open(info_p))
        body += (f'<h2>Training data</h2><p>{sum(info["counts"].values())} classes from {len(info["entries"])} entries: '
                 + ', '.join(f'{v} {k}' for k, v in info['counts'].items()) + f'. Labels: {html.escape(info["label_rule"])}.</p>')
    mp = os.path.join(md, 'model.pt')
    if os.path.exists(mp):
        cfg = torch.load(mp, map_location='cpu', weights_only=True)['config']
        keys = ('epochs', 'width', 'channels', 'batch', 'lr', 'wd', 'aug_noise', 'aug_noise_min', 'aug_bg', 'box', 'smpd')
        body += '<h2>Model</h2><p>' + link(mp, 'model.pt') + ': ' + ', '.join(f'{k} {cfg[k]}' for k in keys if k in cfg) + '.</p>'
    # cross-validation
    cvp = os.path.join(md, 'cv_metrics.json')
    if os.path.exists(cvp):
        cv = json.load(open(cvp))
        mt = cv['out_of_fold']
        c = np.array(mt['confusion'])
        body += (f'<h2>Cross-validation</h2><p>{cv["folds"]} folds by {cv["split"]}; every class predicted by the model that did not '
                 f'train on its {cv["split"]}. Balanced accuracy {fmt(mt["balanced_accuracy"])}, accuracy {fmt(mt["accuracy"])}; '
                 f'AUC protein vs junk {fmt(mt["auc_protein_vs_junk"])}; AUC secondary structure vs likely protein '
                 f'{fmt(mt["auc_ss_vs_likely"])}. ' + link(os.path.join(md, 'cv_report.html'), 'Full report with galleries') + '.</p>')
        body += table(['label \\ predicted'] + list(GROUPS) + ['recall'],
                      [[GROUPS[i]] + [int(c[i, j]) for j in range(3)] + [fmt(mt['recall'][i], 2)] for i in range(3)]
                      + [['precision'] + [fmt(p, 2) for p in mt['precision']] + ['']])
        if cv.get('ss_auc_by_snr'):
            body += ('<p>Secondary structure vs likely protein within particle-SNR bands: '
                     + '; '.join(f'{b["snr"]}: {fmt(b["auc"])} ({b["n_ss"]}/{b["n_likely"]})' for b in cv['ss_auc_by_snr'])
                     + '. Bands well below the overall AUC would mean the network reads the noise level.</p>')
        body += table(['fold', 'held-out', 'classes', 'balanced accuracy', 'on its training classes', 'AUC protein', 'AUC ss'],
                      [[f['fold'], ' '.join(f['groups']), f['validation']['n'], fmt(f['validation']['balanced_accuracy']),
                        fmt(f['training_balanced_accuracy']), fmt(f['validation']['auc_protein_vs_junk']),
                        fmt(f['validation']['auc_ss_vs_likely'])] for f in cv['per_fold']])
    # predictions
    recs = []
    if os.path.isdir(pred_dir):
        for f in sorted(os.listdir(pred_dir)):
            if f.endswith('.json'):
                r = json.load(open(os.path.join(pred_dir, f)))
                if 'counts' in r:
                    recs.append(r)
    if recs:
        rows = []
        for r in recs:
            k = r['counts']
            kept = (f'{r["particles_kept"]:.0f} of {r["particles_total"]:.0f}' if r.get('particles_total') else '')
            name = os.path.splitext(os.path.basename(r['tsv']))[0]
            rows.append([(link(r['html'], name) if r.get('html') else html.escape(name))
                         + f'<br><span class="cap">{html.escape(r["stk"])}</span>', r['classes'], k.get('junk', 0), k.get('secondary_structure', 0),
                         k.get('likely_protein', 0), k.get('empty', 0), kept, html.escape(os.path.basename(r['model']))])
        body += '<h2>Predictions</h2>' + ('<table><tr>' + ''.join(f'<th>{h}</th>' for h in (
            'class averages', 'classes', 'junk', 'secondary structure', 'likely protein', 'empty', 'particles in protein classes',
            'model')) + '</tr>' + ''.join('<tr>' + ''.join(f'<td>{v}</td>' for v in r) + '</tr>' for r in rows) + '</table>')
    open(out, 'w').write(page('Class-average CNN: summary', body))
    print(f'summary: {out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    t = sub.add_parser('train')
    t.add_argument('--data', required=True)
    t.add_argument('--out', required=True)
    t.add_argument('--folds', type=int, default=5, help='cross-validation folds of whole entries; 0: final model only {5}')
    t.add_argument('--split', choices=('entry', 'run'), default='entry')
    t.add_argument('--no-final', action='store_true', help='cross-validation only')
    t.add_argument('--epochs', type=int, default=40)
    t.add_argument('--batch', type=int, default=32)
    t.add_argument('--lr', type=float, default=1e-3)
    t.add_argument('--wd', type=float, default=0.05, help='AdamW weight decay {0.05}')
    t.add_argument('--width', type=int, default=32, help='channels of the first layer {32}')
    t.add_argument('--dropout', type=float, default=0.3)
    t.add_argument('--channels', choices=('avg+halves', 'avg'), default='avg+halves')
    t.add_argument('--eval-every', type=int, default=1, help='validation loss every N epochs {1}')
    t.add_argument('--gallery', type=int, default=16, help='error tiles per confusion cell in the report {16}')
    t.add_argument('--aug-noise', type=float, default=0.5, help='share of classes given fewer particles\' noise {0.5}')
    t.add_argument('--aug-noise-min', type=float, default=0.05, help='least particle fraction of that augmentation {0.05}')
    t.add_argument('--aug-bg', type=float, default=0.5, help='share of classes given noise outside their mask {0.5}')
    t.add_argument('--seed', type=int, default=1)
    t.add_argument('--device', default='auto', help='auto, cpu, cuda, cuda:1, mps {auto}')
    t.add_argument('--nthr', type=int, default=0, help='CPU threads (0: PyTorch default)')
    p = sub.add_parser('predict')
    p.add_argument('--model', required=True)
    p.add_argument('--stk', required=True, help='class averages (MRC); halves from <stk>_even.mrc, <stk>_odd.mrc')
    p.add_argument('--even')
    p.add_argument('--odd')
    p.add_argument('--smpd', type=float, default=0.0, help='pixel size (A) {from the MRC header}')
    p.add_argument('--labels', help="ss_trainer.py's labels.json")
    p.add_argument('--cls', help="SIMPLE's class document (e.g. classdoc_ranked.txt): populations and resolutions")
    p.add_argument('--out', required=True)
    p.add_argument('--html')
    p.add_argument('--device', default='auto')
    p.add_argument('--nthr', type=int, default=0)
    m = sub.add_parser('summary')
    m.add_argument('--model-dir', required=True)
    m.add_argument('--set', help='the simulated set the model was trained on')
    m.add_argument('--pred-dir', help='folder of predict outputs {MODEL_DIR/predictions}')
    m.add_argument('--out', help='{MODEL_DIR/summary.html}')
    a = ap.parse_args()
    {'train': cmd_train, 'predict': cmd_predict, 'summary': cmd_summary}[a.cmd](a)


if __name__ == '__main__':
    main()
