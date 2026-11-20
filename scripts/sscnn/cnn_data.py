#!/usr/bin/env python3
"""The CNN's training data: every class average of every simulated run, with its label, as one array.

Per class (classes with no particles are left out):
  protein_frac   1 - (junk particles / particles), the target of the protein head
  label          0 junk (protein_frac < 0.5); otherwise 1 secondary structure present if the run's map is not low-passed
                 (<entry>_r<k>), 2 likely protein if it is (<entry>_lp<L>_r<k>, <entry>_pr<L>_r<k>). The secondary
                 structure label is the run's kind, not a measurement: a full-map class whose helices do not show is
                 still labelled 1 (visibility numbers are kept beside it to find these).
  snr_px         the noise per pixel as the SNR of a 150 A particle (runs of any date compare)
  channels       the class average and its even and odd halves, all three scaled by the average's mean and standard
                 deviation inside the box's inscribed circle (normalise(); predict uses the same function)
Written to --out: x.npy (float16, classes x 3 x box x box), meta.tsv (one row per class), info.json.

Usage: python -I cnn_data.py --set SET_DIR [--set SET_DIR ...] --out DATA_DIR
A SET_DIR is a folder holding runs/ (build_emdb_set.py's --set) or a runs folder itself.
"""
import argparse
import gzip
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ss_trainer as T   # noqa: E402

NAME = re.compile(r'^(?P<entry>.+?)_(?:(?P<kind>lp|pr)(?P<cut>\d+)_)?r(?P<rep>\d+)$')
GROUPS = ('junk', 'secondary_structure', 'likely_protein')
META = ('run', 'entry', 'kind', 'cutoff', 'cls', 'pop', 'njunk', 'protein_frac', 'label', 'junk_kind', 'snr', 'res_A',
        'cnr_helix_12_6', 'cnr_strand_12_6', 'snr_px', 'diam_A')
SNR_REF_DIAM = 150.0      # snr_px: the SNR a 150 A particle would need for the run's noise per pixel (build_sim_set.py)


def normalise(avg, ev, od):
    """The three images scaled together by the average's statistics inside the inscribed circle."""
    n = avg.shape[-1]
    g = np.arange(n) - n // 2
    inner = np.hypot(*np.meshgrid(g, g, indexing='ij')) < n / 2
    m, s = avg[inner].mean(), avg[inner].std()
    s = s if s > 0 else 1.0
    return np.stack([(avg - m) / s, (ev - m) / s, (od - m) / s])


def resample(img, smpd, smpd_out, box_out):
    """An image at smpd_out (Fourier crop or pad to m = floor(n smpd / smpd_out + 1/2) pixels), then cropped or
    padded to box_out. The transforms take their origin at the centre pixel n // 2 (numpy's fftshift convention)
    and the crop or pad keeps it, so the centre stays the centre pixel at every size, odd or even; SIMPLE's
    cavg_cnn_predict computes the same."""
    n = img.shape[-1]
    m = int(np.floor(n * smpd / smpd_out + 0.5))
    if m != n:
        F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))
        G = np.zeros((m, m), complex)
        if m <= n:
            a = n // 2 - m // 2
            G[:, :] = F[a:a + m, a:a + m]
        else:
            b = m // 2 - n // 2
            G[b:b + n, b:b + n] = F
        img = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.ifftshift(G))))
    n = img.shape[-1]
    out = np.zeros((box_out, box_out))
    if n >= box_out:
        a = n // 2 - box_out // 2
        out[:, :] = img[a:a + box_out, a:a + box_out]
    else:
        b = box_out // 2 - n // 2
        out[b:b + n, b:b + n] = img
    return out


def class_junk_kinds(rdir, ncls):
    """Per class, the most common junk kind among its junk particles ('' if none)."""
    j = json.load(open(os.path.join(rdir, 'junk.json')))
    nprot, kinds = j['nprot'], j['kinds']
    cls_re = re.compile(r'(?:^| )class=(\S+)')
    counts = [dict() for _ in range(ncls + 1)]
    with gzip.open(os.path.join(rdir, 'ptcl2D.txt.gz'), 'rt') as f:
        for i, line in enumerate(f):
            if i < nprot:
                continue
            m = cls_re.search(line)
            c = int(float(m.group(1))) if m else 0
            if 1 <= c <= ncls and i - nprot < len(kinds):
                k = kinds[i - nprot]
                counts[c][k] = counts[c].get(k, 0) + 1
    return [max(d, key=d.get) if d else '' for d in counts]


def runs_in(path):
    rd = os.path.join(path, 'runs') if os.path.isdir(os.path.join(path, 'runs')) else path
    return [os.path.join(rd, d) for d in sorted(os.listdir(rd)) if os.path.exists(os.path.join(rd, d, 'run.json'))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--set', action='append', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--box', type=int, default=160, help='box of the CNN input (pixels) {160}')
    ap.add_argument('--smpd', type=float, default=2.0, help='pixel size of the CNN input (A) {2.0}')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = [(os.path.basename(os.path.abspath(s)), r) for s in a.set for r in runs_in(os.path.abspath(s))]
    names = [os.path.basename(r) for _, r in runs]
    xs, rows = [], []
    for sname, rdir in runs:
        base = os.path.basename(rdir)
        name = base if names.count(base) == 1 else f'{sname}:{base}'   # a name in two sets carries its set's
        mt = NAME.match(base)
        if not mt:
            print(f'{name}: name not understood, skipped')
            continue
        info = json.load(open(os.path.join(rdir, 'run.json')))
        ims = [T.read_mrc(os.path.join(rdir, f))[0] for f in ('cavgs.mrc', 'cavgs_even.mrc', 'cavgs_odd.mrc')]
        smpd = info['smpd']
        tab = [l.rstrip('\n').split('\t') for l in open(os.path.join(rdir, 'classes.tsv'))]
        hdr = tab[0]
        jk = class_junk_kinds(rdir, len(ims[0]))
        kind = mt.group('kind') or 'map'
        snr_px = info.get('snr_sim', info['snr']) * (SNR_REF_DIAM / info['diam']) ** 4   # runs before 2026-11-19 included
        nkept = 0
        for c in tab[1:]:
            d = dict(zip(hdr, c))
            i = int(d['class']) - 1
            pop, nj = int(d['pop']), int(d['njunk'])
            if pop == 0:
                continue
            chans = [im[i].astype(np.float64) for im in ims]
            if abs(smpd - a.smpd) > 1e-3 or chans[0].shape[-1] != a.box:
                chans = [resample(im, smpd, a.smpd, a.box) for im in chans]
            xs.append(normalise(*chans).astype(np.float16))
            pf = 1.0 - nj / pop
            label = 0 if pf < 0.5 else (1 if kind == 'map' else 2)
            rows.append((name, mt.group('entry'), kind, mt.group('cut') or '', i + 1, pop, nj, round(pf, 4), label,
                         jk[i + 1] if nj else '', info['snr'], d['res_A'], d['cnr_helix_12_6'], d['cnr_strand_12_6'],
                         round(snr_px, 7), round(info['diam'], 1)))
            nkept += 1
        print(f'{name}: {nkept} classes')
    x = np.stack(xs)
    np.save(os.path.join(a.out, 'x.npy'), x)
    with open(os.path.join(a.out, 'meta.tsv'), 'w') as f:
        f.write('\t'.join(META) + '\n')
        for r in rows:
            f.write('\t'.join(str(v) for v in r) + '\n')
    lab = np.array([r[8] for r in rows])
    counts = {g: int((lab == k).sum()) for k, g in enumerate(GROUPS)}
    json.dump({'box': a.box, 'smpd': a.smpd, 'channels': ['average', 'even', 'odd'], 'groups': GROUPS, 'counts': counts,
               'entries': sorted({r[1] for r in rows}), 'runs': len(runs), 'sets': [os.path.abspath(s) for s in a.set],
               'label_rule': 'junk: protein_frac < 0.5; else secondary structure if the run map is not low-passed, '
                             'likely protein if it is (lp or pr)'},
              open(os.path.join(a.out, 'info.json'), 'w'), indent=1)
    print(f'{len(rows)} classes from {len(runs)} runs of {len({r[1] for r in rows})} entries: '
          + ', '.join(f'{v} {k}' for k, v in counts.items()) + f'; {x.nbytes / 1e6:.0f} MB')


if __name__ == '__main__':
    main()
