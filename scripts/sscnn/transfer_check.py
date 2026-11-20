#!/usr/bin/env python3
"""Apply a helix-only model trained on one data set's labels to another data set's class averages.

The model is the trainer's logistic regression on the four helix measurements (no population, no resolution),
standardised with the training set's statistics so that both sets are on one scale. Writes a page of the target's
class averages ranked by predicted p(helices visible), each with its even and odd halves, for judging by eye.

Usage: transfer_check.py --train-stk S --train-cls C --train-mskdiam A --labels labels.json
                         --stk S --cls C --mskdiam A --out page.html
"""
import argparse
import base64
import html
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_trainer as T   # noqa: E402

HELIX = [0, 1, 2, 3]


def features(stk, cls, mskdiam):
    avg, smpd, _, _ = T.read_mrc(stk)
    base = stk[:-4]
    ev, _, _, _ = T.read_mrc(base + '_even.mrc')
    od, _, _, _ = T.read_mrc(base + '_odd.mrc')
    pop, res = T.read_cls(cls, avg.shape[0])
    m = T.Measurer(avg.shape[1], smpd, mskdiam)
    X = np.full((avg.shape[0], 4), np.nan)
    for i in range(avg.shape[0]):
        if avg[i].std() == 0 or (np.isfinite(pop[i]) and pop[i] <= 0):
            continue
        X[i] = m.measure(avg[i], ev[i], od[i])
    return X, pop, res, avg, ev, od


def main():
    ap = argparse.ArgumentParser()
    for k in ('train-stk', 'train-cls', 'labels', 'stk', 'cls', 'out'):
        ap.add_argument('--' + k, required=True)
    ap.add_argument('--train-mskdiam', type=float, required=True)
    ap.add_argument('--mskdiam', type=float, required=True)
    a = ap.parse_args()
    Xt, _, _, _, _, _ = features(a.train_stk, a.train_cls, a.train_mskdiam)
    lab = json.load(open(a.labels))
    pairs = [(int(k) - 1, v) for k, v in lab.items() if v in ('yes', 'no') and np.isfinite(Xt[int(k) - 1]).all()]
    idx = np.array([i for i, _ in pairs])
    y = np.array([1.0 if v == 'yes' else 0.0 for _, v in pairs])
    live_t = np.isfinite(Xt).all(1)
    mu, sd = Xt[live_t].mean(0), Xt[live_t].std(0)
    w = T.fit_logistic((Xt[idx] - mu) / sd, y)
    names = [T.FEATURES[k] for k in HELIX]
    print(f'trained on {int(y.sum())} yes, {int((1 - y).sum())} no; weights: ' +
          ', '.join(f'{n} {c:+.2f}' for n, c in zip(names, w[1:])) + f'; intercept {w[0]:+.2f}')
    X, pop, res, avg, ev, od = features(a.stk, a.cls, a.mskdiam)
    live = np.where(np.isfinite(X).all(1))[0]
    p = T.predict(w, (X[live] - mu) / sd)
    order = live[np.argsort(-p)]
    pmap = dict(zip(live, p))
    print(f'target: {len(live)} classes with content; predicted p >= 0.5 for {int((p >= 0.5).sum())}; '
          f'median p {np.median(p):.2f}')
    print('top 10 (class, p, pop, res):', [(int(i) + 1, round(float(pmap[i]), 2), int(pop[i]), round(float(res[i]), 1)) for i in order[:10]])
    for n, k in zip(names, range(4)):
        print(f'  {n:12s} training median {np.median(Xt[live_t, k]):.3f}   target median {np.median(X[live, k]):.3f}')

    def uri(img, size=192):
        return 'data:image/png;base64,' + base64.b64encode(T.png(img, size)).decode()
    cells = []
    for i in order:
        cells.append(
            f'<figure style="margin:0 0 6px 0;font-size:11px;color:#52514e">'
            f'<img src="{uri(avg[i])}" width="192"><img src="{uri(ev[i], 96)}" width="96" style="vertical-align:top;margin-left:4px">'
            f'<img src="{uri(od[i], 96)}" width="96" style="vertical-align:top;margin-left:2px">'
            f'<figcaption>class {i + 1} · p={pmap[i]:.2f} · n={int(pop[i])} · {res[i]:.1f} Å'
            f' &nbsp;(average · even · odd)</figcaption></figure>')
    page = (f'<!doctype html><meta charset="utf-8"><title>Helix model transfer</title>'
            f'<body style="background:#fcfcfb;color:#0b0b0b;font-family:system-ui,sans-serif;margin:20px">'
            f'<h1 style="font-size:16px;margin:0 0 4px">{html.escape(os.path.basename(a.stk))}: classes ranked by predicted '
            f'p(helices visible), highest first</h1><p style="font-size:12px;color:#52514e;margin:0 0 12px">Model: logistic '
            f'regression on the four helix measurements only, trained on {html.escape(a.labels)} '
            f'({int(y.sum())} yes, {int((1 - y).sum())} no). Population and resolution are not inputs.</p>'
            f'<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(400px,1fr));gap:10px">{"".join(cells)}</div></body>')
    open(a.out, 'w').write(page)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
