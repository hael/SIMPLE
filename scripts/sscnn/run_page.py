#!/usr/bin/env python3
"""A page of a simulated run's classes: each class average beside its helix and strand truth, with its population, junk
share, resolution and visibility numbers (classes.tsv), ordered by population.

Usage: python -I run_page.py RUN_DIR [RUN_DIR ...] --out page.html
"""
import argparse
import base64
import csv
import html
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_trainer as T   # noqa: E402
from sim_truth import load_truth   # noqa: E402


def uri(img, size):
    return 'data:image/png;base64,' + base64.b64encode(T.png(np.asarray(img, np.float64), size)).decode()


def truth_uri(img, ref, size):
    """A truth map on the scale of the full noise-free average (so a weak truth looks weak); grey where empty."""
    s = ref.std() or 1.0
    m = ref.mean()
    g = np.clip((img - m + 3 * s) / (6 * s), 0, 1) * 255
    return uri(g, size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    parts = []
    for rd in a.runs:
        info = json.load(open(os.path.join(rd, 'run.json')))
        rows = list(csv.DictReader(open(os.path.join(rd, 'classes.tsv')), delimiter='\t'))
        cav = T.read_mrc(os.path.join(rd, 'cavgs.mrc'))[0]
        tr = {k: (lambda x: x if x is not None else np.zeros_like(cav))(load_truth(rd, f'{k}_on_cavg'))
              for k in ('full', 'helix', 'strand')}
        rows = [r for r in rows if int(r['pop']) > 0]
        rows.sort(key=lambda r: -int(r['pop']))
        cells = []
        for r in rows:
            c = int(r['class']) - 1
            jf = int(r['njunk']) / max(int(r['pop']), 1)
            border = '#eb6834' if jf > 0.5 else '#e4e3df'
            cells.append(
                f'<figure style="margin:0;font-size:10px;color:#52514e;width:300px">'
                f'<img src="{uri(cav[c], 128)}" width="128" style="border:3px solid {border};border-radius:3px">'
                f'<img src="{truth_uri(tr["helix"][c], tr["full"][c], 72)}" width="72" style="vertical-align:top;margin-left:4px">'
                f'<img src="{truth_uri(tr["strand"][c], tr["full"][c], 72)}" width="72" style="vertical-align:top;margin-left:2px">'
                f'<figcaption>class {c + 1} · n={r["pop"]} · junk {100 * jf:.0f}% · {float(r["res_A"]):.1f} Å<br>'
                f'helix cnr {float(r["cnr_helix_12_6"]):.2f} · strand cnr {float(r["cnr_strand_12_6"]):.2f} '
                f'(12–6 Å), {float(r["cnr_strand_6_4"]):.2f} (6–4 Å)</figcaption></figure>')
        parts.append(f'<h2 style="font-size:14px;margin:16px 0 6px">{html.escape(info["model"])} · particle SNR {info["snr"]:.4g} · '
                     f'{info["nptcls"]} protein + {info["njunk"]} junk particles · {len(rows)} classes with particles</h2>'
                     f'<div style="display:flex;flex-wrap:wrap;gap:10px">{"".join(cells)}</div>')
    page = ('<!doctype html><meta charset="utf-8"><title>Simulated runs</title>'
            '<body style="background:#fcfcfb;color:#0b0b0b;font-family:system-ui,sans-serif;margin:16px">'
            '<p style="font-size:12px;color:#52514e;margin:0">Each class: the class average (orange border: more than half its '
            'particles are junk), then its helix truth and its strand truth on the scale of its noise-free average '
            '(grey: nothing). cnr: the truth\'s contrast against the class\'s noise in a band, from its halves.</p>'
            + ''.join(parts) + '</body>')
    open(a.out, 'w').write(page)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
