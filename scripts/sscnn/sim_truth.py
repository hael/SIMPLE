#!/usr/bin/env python3
"""Ground truth for class averages made by solve2D from simulated particles, one noise-free stack at a time.

  average  each class's average of its particles' noise-free images of one kind (full, helix or strand), each moved by
           the 2D alignment solve2D gave the noisy particle (in-plane rotation e3, shift x, y). Particles from index
           --nprot on are junk: their images are zero, so they dilute the truth as they dilute the class. With
           --label full and no --conv, SIMPLE's convention for applying e3 and the shift (signs, order) is found: of the
           eight combinations, the one whose averages correlate best with solve2D's class averages; it is written to
           conv.json and reused for the other kinds. Writes truth_<label>.mrc.
  summary  per class: population, junk count, the noise-free full average's correlation with the class average, the
           class average's transfer of the truth, shell by shell (H = Re sum C F* / sum |F|^2 of the class average C against
           the full truth F under a soft mask; it carries the units, SIMPLE's restoration filter and the alignment blur; a
           class with too little protein to estimate it, MIN_PROT_H, takes the median of the run's reliable classes),
           applied to all three truths and written as truth_<kind>_on_cavg.mrc (what of each survives in the class average);
           the helix and strand ones also high-passed (BAND_TARGET) as truth_<kind>_band.mrc: only the detail by which
           helices and sheets are recognised, not their blurred silhouette, the training target,
           and the visibility numbers of helices and strands: the contrast of the scaled truth map against the class average's
           noise in a band (that noise taken as at least VIS_FLOOR of the class's total power, so that a band the
           restoration emptied does not inflate the ratio), sqrt(sum |T|^2 / sum |N|^2), N from the halves ((even - odd) / 2 carries the average's noise
           and none of its signal). Bands: 12-6 A (helix rods and spacing; sheets as slabs, paired sheets), and 6-4 A
           (strand spacing, 4.8 A). Writes classes.tsv.

Usage: python -I sim_truth.py average --label full|helix|strand --clean S --cavgs C --ptcl2d P --nprot N --mskdiam A
                                      --out DIR [--conv DIR/conv.json]
       python -I sim_truth.py summary --cavgs C --cls2d C --ptcl2d P --nprot N --mskdiam A --out DIR
"""
import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate_cavgs as sc   # noqa: E402
import ss_trainer as T        # noqa: E402

BANDS = {'12_6': (1 / 12.0, 1 / 6.0), '6_4': (1 / 6.0, 1 / 4.0)}
BAND_TARGET = (15.0, 12.0)   # the band-passed truth: nothing coarser than 15 A, all of what is finer than 12 A (cosine edge)
MIN_PROT_H = 10              # a class's own transfer needs this many protein particles (and at most half junk, cc_full >= 0.5)
VIS_FLOOR = 1e-4             # visibility floor: the band's noise is taken as at least this share of the class's total power
H_CAP = 1.5                  # the transfer at most this times its 40-20 A level; zero where the truth has no power


def load_truth(rdir, key):
    """A run's truth maps truth_<key> (full_on_cavg, helix_on_cavg, strand_on_cavg, helix_band, strand_band, full,
    helix, strand), from truth_<key>.mrc or from truth.npz (--compact); None if neither holds them."""
    p = os.path.join(rdir, f'truth_{key}.mrc')
    if os.path.exists(p):
        return T.read_mrc(p)[0]
    z = os.path.join(rdir, 'truth.npz')
    if os.path.exists(z):
        with np.load(z) as d:
            if key in d:
                return d[key].astype(np.float64)
    return None


def read_ptcl2d(path):
    keys = ('class', 'e3', 'x', 'y', 'state')
    pat = {k: re.compile(r'(?:^| )' + k + r'=(\S+)') for k in keys}
    return np.array([[float(pat[k].search(l).group(1)) if pat[k].search(l) else 0.0 for k in keys] for l in open(path)])


def shift_img(img, sx, sy):
    n = img.shape[0]
    f = np.fft.fftfreq(n)
    FX, FY = np.meshgrid(f, f)
    return np.real(np.fft.ifft2(np.fft.fft2(img) * np.exp(-2j * np.pi * (FX * sx + FY * sy))))


def rotate_img(img, deg):
    n = img.shape[0]
    c = n // 2
    a = np.radians(deg)
    y, x = np.mgrid[0:n, 0:n] - c
    xs = np.cos(a) * x + np.sin(a) * y + c
    ys = -np.sin(a) * x + np.cos(a) * y + c
    x0, y0 = np.floor(xs).astype(int), np.floor(ys).astype(int)
    fx, fy = xs - x0, ys - y0
    out = np.zeros_like(img)
    for dx, dy, w in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)), (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
        xi, yi = x0 + dx, y0 + dy
        ok = (xi >= 0) & (xi < n) & (yi >= 0) & (yi < n)
        out[ok] += img[yi[ok], xi[ok]] * w[ok]
    return out


def transform(img, e3, x, y, conv):
    rs, ss, shift_first = conv
    if shift_first:
        return rotate_img(shift_img(img, ss * x, ss * y), rs * e3)
    return shift_img(rotate_img(img, rs * e3), ss * x, ss * y)


def masked_cc(a, b, mask):
    a, b = a[mask] - a[mask].mean(), b[mask] - b[mask].mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else float('nan')


def setup(a):
    cav, smpd, _, _ = T.read_mrc(a.cavgs)
    P = read_ptcl2d(a.ptcl2d)
    n = cav.shape[1]
    r = np.hypot(*(np.mgrid[0:n, 0:n] - n // 2)) * smpd
    cls = P[:, 0].astype(int)
    live = (P[:, 4] > 0) & (cls >= 1)
    pops = np.bincount(cls[live], minlength=cav.shape[0] + 1)[1:]
    return cav, smpd, P, r <= a.mskdiam / 2, cls, live, pops


def average(a):
    cav, smpd, P, mask, cls, live, pops = setup(a)
    clean, smpd_p, _, _ = T.read_mrc(a.clean)
    if clean.shape[0] != a.nprot or abs(smpd - smpd_p) > 1e-3 or clean.shape[1] != cav.shape[1]:
        raise SystemExit('the clean stack does not match the protein particles or the class averages')
    prot = np.arange(len(P)) < a.nprot
    if a.conv:
        c = json.load(open(a.conv))
        conv = (c['rot_sign'], c['shift_sign'], c['shift_first'])
    elif a.label == 'full':
        convs = [(rs, ss, o) for rs in (1, -1) for ss in (1, -1) for o in (True, False)]
        top = np.argsort(-pops)[:10]
        rng = np.random.default_rng(0)
        scores = []
        for cv in convs:
            cc = []
            for c in top:
                mem = rng.permutation(np.where(live & prot & (cls == c + 1))[0])[:60]
                if len(mem):
                    avg = np.mean([transform(clean[i], P[i, 1], P[i, 2], P[i, 3], cv) for i in mem], axis=0)
                    cc.append(masked_cc(avg, cav[c], mask))
            scores.append(float(np.nanmean(cc)))
        conv = convs[int(np.argmax(scores))]
        json.dump({'rot_sign': conv[0], 'shift_sign': conv[1], 'shift_first': conv[2], 'cc_best': max(scores),
                   'cc_all': scores}, open(os.path.join(a.out, 'conv.json'), 'w'), indent=1)
        print(f'convention: rotation sign {conv[0]:+d}, shift sign {conv[1]:+d}, {"shift then rotate" if conv[2] else "rotate then shift"}; '
              f'mean cc {max(scores):.3f} (others {", ".join(f"{s:.3f}" for i, s in enumerate(scores) if i != int(np.argmax(scores)))})')
    else:
        raise SystemExit('give --conv, or run --label full first')
    acc = np.zeros(cav.shape)
    for i in np.where(live & prot)[0]:
        acc[cls[i] - 1] += transform(clean[i], P[i, 1], P[i, 2], P[i, 3], conv)
    acc /= np.maximum(pops, 1)[:, None, None]
    sc.write_mrc(os.path.join(a.out, f'truth_{a.label}.mrc'), acc, smpd)
    print(f'wrote truth_{a.label}.mrc ({len(cav)} classes)')


def summary(a):
    cav, smpd, P, mask, cls, live, pops = setup(a)
    ev, _, _, _ = T.read_mrc(a.cavgs[:-4] + '_even.mrc')
    od, _, _, _ = T.read_mrc(a.cavgs[:-4] + '_odd.mrc')
    junk = live & (np.arange(len(P)) >= a.nprot)
    njunk = np.bincount(cls[junk], minlength=cav.shape[0] + 1)[1:]
    truth = {}
    for lab in ('full', 'helix', 'strand'):
        t = load_truth(a.out, lab)
        truth[lab] = t if t is not None else np.zeros(cav.shape)
    _, res = T.read_cls(a.cls2d, cav.shape[0])
    n = cav.shape[1]
    f = np.fft.fftfreq(n, d=smpd)
    FR = np.hypot(*np.meshgrid(f, f))
    bands = {k: (FR >= lo) & (FR <= hi) for k, (lo, hi) in BANDS.items()}
    soft = np.clip((a.mskdiam / 2 - np.hypot(*(np.mgrid[0:n, 0:n] - n // 2)) * smpd) / 10.0 + 0.5, 0, 1)
    rows = []
    scaled = {k: np.zeros_like(v) for k, v in truth.items()}
    shell = np.rint(FR * n * smpd).astype(int)
    nsh = shell.max() + 1
    lowband = (np.arange(nsh) >= int(n * smpd / 40.0)) & (np.arange(nsh) <= int(n * smpd / 20.0))
    # each class's transfer; a class with too little protein to estimate it takes the median of the reliable ones
    Hs, reliable = [], []
    for c in range(cav.shape[0]):
        Ff = np.fft.fft2(truth['full'][c] * soft)
        Cf = np.fft.fft2(cav[c] * soft)
        num = np.bincount(shell.ravel(), np.real(Cf * np.conj(Ff)).ravel(), nsh)
        den = np.bincount(shell.ravel(), (np.abs(Ff) ** 2).ravel(), nsh)
        H = np.where(den > 1e-9 * den.max() if den.max() > 0 else False, num / np.where(den > 0, den, 1), 0.0)
        H = np.clip(H, 0, None)
        H[0] = 0.0                                                   # the mean carries no structure
        if lowband.any():
            H = np.minimum(H, H_CAP * np.median(H[lowband]))
        Hs.append(H)
        ok = False
        if pops[c] > 1 and pops[c] - njunk[c] >= MIN_PROT_H and njunk[c] <= 0.5 * pops[c]:
            ok = masked_cc(truth['full'][c], cav[c], mask) >= 0.5
        reliable.append(ok)
    H_run = np.median(np.array(Hs)[np.array(reliable)], axis=0) if any(reliable) else np.zeros(nsh)
    for c in range(cav.shape[0]):
        H = Hs[c] if reliable[c] else H_run
        Nf = np.abs(np.fft.fft2(0.5 * (ev[c] - od[c]) * soft)) ** 2
        ctot = float((np.abs(np.fft.fft2((cav[c] - cav[c].mean()) * soft)) ** 2).sum())
        Hmap = H[shell]
        for k in truth:
            scaled[k][c] = np.real(np.fft.ifft2(np.fft.fft2(truth[k][c]) * Hmap))
        hb = {b: float(np.median(H[np.unique(shell[m])])) for b, m in bands.items()}
        row = {'class': c + 1, 'pop': int(pops[c]), 'njunk': int(njunk[c]), 'res_A': float(res[c]),
               'cc_full': masked_cc(truth['full'][c], cav[c], mask) if pops[c] > 1 else float('nan'),
               'H_own': int(reliable[c]), 'H_12_6': hb['12_6'], 'H_6_4': hb['6_4']}
        for lab in ('helix', 'strand'):
            Tf = np.abs(np.fft.fft2(scaled[lab][c] * soft)) ** 2
            for b, m in bands.items():
                nb = Nf[m].sum() + VIS_FLOOR * ctot
                row[f'cnr_{lab}_{b}'] = float(np.sqrt(Tf[m].sum() / nb)) if nb > 0 and pops[c] > 0 else float('nan')
        rows.append(row)
    derived = {f'{k}_on_cavg': scaled[k] for k in scaled}
    lo, hi = 1.0 / BAND_TARGET[0], 1.0 / BAND_TARGET[1]
    hp = np.where(FR <= lo, 0.0, np.where(FR >= hi, 1.0, 0.5 - 0.5 * np.cos(np.pi * (FR - lo) / (hi - lo))))
    for k in ('helix', 'strand'):
        derived[f'{k}_band'] = np.real(np.fft.ifft2(np.fft.fft2(scaled[k], axes=(1, 2)) * hp[None], axes=(1, 2)))
    z = os.path.join(a.out, 'truth.npz')
    if os.path.exists(z):              # a compact run: the derived maps join the raw ones in truth.npz
        with np.load(z) as d:
            keep = {k: d[k] for k in d.files}
        keep.update({k: v.astype(np.float16) for k, v in derived.items()})
        np.savez_compressed(z, **keep)
    else:
        for k, v in derived.items():
            sc.write_mrc(os.path.join(a.out, f'truth_{k}.mrc'), v, smpd)
    keys = list(rows[0].keys())
    with open(os.path.join(a.out, 'classes.tsv'), 'w') as fh:
        fh.write('\t'.join(keys) + '\n')
        for r in rows:
            fh.write('\t'.join(f'{r[k]:.4g}' if isinstance(r[k], float) else str(r[k]) for k in keys) + '\n')
    live_c = pops > 1
    print(f'{live_c.sum()} classes with particles, {int(((njunk > 0.5 * pops) & live_c).sum())} mostly junk; '
          f'cc_full median {np.nanmedian([r["cc_full"] for r in rows]):.3f}; '
          f'protein classes: cnr helix 12-6 median {np.nanmedian([r["cnr_helix_12_6"] for r in rows if r["njunk"] < 0.5 * r["pop"]]):.2f}, '
          f'cnr strand 12-6 median {np.nanmedian([r["cnr_strand_12_6"] for r in rows if r["njunk"] < 0.5 * r["pop"]]):.2f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='mode', required=True)
    p = sub.add_parser('average')
    for k in ('label', 'clean', 'cavgs', 'ptcl2d', 'out'):
        p.add_argument('--' + k, required=True)
    p.add_argument('--nprot', type=int, required=True)
    p.add_argument('--mskdiam', type=float, required=True)
    p.add_argument('--conv', default='')
    q = sub.add_parser('summary')
    for k in ('cavgs', 'cls2d', 'ptcl2d', 'out'):
        q.add_argument('--' + k, required=True)
    q.add_argument('--nprot', type=int, required=True)
    q.add_argument('--mskdiam', type=float, required=True)
    a = ap.parse_args()
    average(a) if a.mode == 'average' else summary(a)


if __name__ == '__main__':
    main()
