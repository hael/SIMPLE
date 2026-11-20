#!/usr/bin/env python3
"""Negative variants of a prepared map: structure without visible secondary structure.

  lp<L>  the full, helix and strand maps low-passed at L A: the particle's shape, no helices or sheets
  pr<L>  the full map's phases randomised beyond L A, confined to the particle's soft envelope, and rescaled shell by
         shell to the map's own power beyond L: high-resolution texture with the map's power spectrum and no secondary-structure geometry; the helix
         and strand truth low-passed at L (nothing real survives beyond it)
The random phases come from the Fourier transform of a real noise volume, so the map stays real (Friedel symmetric).
Reads DIR/_maps_<name>/ (full.mrc, helix.mrc, strand.mrc, model.json) and writes DIR/_maps_<name>_lp<L>/ and
DIR/_maps_<name>_pr<L>/ in the same layout, for build_sim_set.py --maps-from.

Usage: python -I map_variants.py --maps DIR/_maps_<name> --out DIR --cutoff 15 [--seed 1]
"""
import argparse
import json
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate_cavgs as sc   # noqa: E402


def read_vol(path):
    with open(path, 'rb') as f:
        h = f.read(1024)
        n = int(np.frombuffer(h[:4], '<i4')[0])
        smpd = float(np.frombuffer(h[40:44], '<f4')[0]) / n
        return np.fromfile(f, '<f4', n ** 3).reshape(n, n, n).astype(np.float64), smpd


def freq(n, smpd):
    f = np.fft.fftfreq(n, smpd)
    fr = np.fft.rfftfreq(n, smpd)
    return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + fr[None, None, :] ** 2)


def lowpass(v, smpd, res, edge=0.3):
    """Low-pass with a one-sided cosine edge: all of what is coarser than res * (1 + edge), nothing finer than res."""
    f = freq(v.shape[0], smpd)
    lo, hi = 1.0 / (res * (1 + edge)), 1.0 / res
    w = np.where(f <= lo, 1.0, np.where(f >= hi, 0.0, 0.5 + 0.5 * np.cos(np.pi * (f - lo) / (hi - lo))))
    return np.fft.irfftn(np.fft.rfftn(v) * w, s=v.shape, axes=(0, 1, 2)), w


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--maps', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--cutoff', type=float, default=15.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--kinds', default='lp,pr', help='which variants: lp, pr or both')
    a = ap.parse_args()
    name = os.path.basename(os.path.normpath(a.maps))[len('_maps_'):]
    full, smpd = read_vol(os.path.join(a.maps, 'full.mrc'))
    vols = {k: read_vol(os.path.join(a.maps, f'{k}.mrc'))[0] for k in ('helix', 'strand') if os.path.exists(os.path.join(a.maps, f'{k}.mrc'))}
    model = json.load(open(os.path.join(a.maps, 'model.json')))
    L = a.cutoff
    tag = f'{L:g}'
    kinds = a.kinds.split(',')
    # lp: everything low-passed
    d = os.path.join(a.out, f'_maps_{name}_lp{tag}')
    if 'lp' not in kinds:
        d = None
    if d:
        os.makedirs(d, exist_ok=True)
        sc.write_mrc(os.path.join(d, 'full.mrc'), lowpass(full, smpd, L)[0].astype(np.float32), smpd)
        for k, v in vols.items():
            sc.write_mrc(os.path.join(d, f'{k}.mrc'), lowpass(v, smpd, L)[0].astype(np.float32), smpd)
        json.dump(dict(model, variant=f'low-passed at {L:g} A'), open(os.path.join(d, 'model.json'), 'w'), indent=1)
    if 'pr' not in kinds:
        print(f'{name}: wrote {",".join(kinds)}{tag}')
        return
    # pr: phases randomised beyond L, amplitudes kept, confined to the soft envelope
    rng = np.random.default_rng(a.seed)
    F = np.fft.rfftn(full)
    lo, w = lowpass(full, smpd, L)
    ph = np.fft.rfftn(rng.normal(size=full.shape))
    ph /= np.maximum(np.abs(ph), 1e-30)
    hi = np.fft.irfftn(np.abs(F) * ph * (1 - w), s=full.shape, axes=(0, 1, 2))
    env, _ = lowpass(full, smpd, 20.0)
    env = np.clip(env / (0.1 * env.max()), 0, 1)
    env, _ = lowpass(env, smpd, 15.0)
    env = np.clip(env, 0, 1)
    # the envelope removes most of the texture's power (it filled the box): restore the map's power shell by shell
    G = np.fft.rfftn(hi * env)
    f = freq(full.shape[0], smpd)
    sh = np.rint(f * full.shape[0] * smpd).astype(int)
    nsh = sh.max() + 1
    p_orig = np.bincount(sh.ravel(), (np.abs(F) ** 2 * (1 - w) ** 2).ravel(), nsh)
    p_now = np.bincount(sh.ravel(), (np.abs(G) ** 2).ravel(), nsh)
    gain = np.sqrt(np.where(p_now > 0, p_orig / np.where(p_now > 0, p_now, 1), 0.0))
    G *= gain[sh] * (1 - w)
    pr = lo + np.fft.irfftn(G, s=full.shape, axes=(0, 1, 2))
    d = os.path.join(a.out, f'_maps_{name}_pr{tag}')
    os.makedirs(d, exist_ok=True)
    sc.write_mrc(os.path.join(d, 'full.mrc'), pr.astype(np.float32), smpd)
    for k, v in vols.items():
        sc.write_mrc(os.path.join(d, f'{k}.mrc'), lowpass(v, smpd, L)[0].astype(np.float32), smpd)
    json.dump(dict(model, variant=f'phases randomised beyond {L:g} A'), open(os.path.join(d, 'model.json'), 'w'), indent=1)
    pf = (np.abs(np.fft.rfftn(pr)) ** 2)[w < 1].sum() / (np.abs(F) ** 2)[w < 1].sum()
    print(f'{name}: wrote lp{tag} and pr{tag}; power beyond {L:g} A kept by pr{tag} (after the envelope): {100 * pf:.0f}%')


if __name__ == '__main__':
    main()
