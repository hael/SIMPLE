#!/usr/bin/env python3
"""Compare simulated and real class averages on the trainer's measurements, binned by class resolution.

Usage: python -I realism.py --sim-stk S --sim-cls C --sim-mskdiam A --real-stk S --real-cls C --real-mskdiam A
Each stack needs <stk>_even.mrc and <stk>_odd.mrc beside it; classes with fewer than 20 particles are left out.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_trainer as T   # noqa: E402


def measure(stk, cls, mskdiam):
    avg, smpd, _, _ = T.read_mrc(stk)
    ev, _, _, _ = T.read_mrc(stk[:-4] + '_even.mrc')
    od, _, _, _ = T.read_mrc(stk[:-4] + '_odd.mrc')
    pop, res = T.read_cls(cls, avg.shape[0])
    m = T.Measurer(avg.shape[1], smpd, mskdiam)
    return np.array([(res[i], pop[i], *m.measure(avg[i], ev[i], od[i])) for i in range(avg.shape[0]) if pop[i] >= 20])


def main():
    ap = argparse.ArgumentParser()
    for k in ('sim-stk', 'sim-cls', 'real-stk', 'real-cls'):
        ap.add_argument('--' + k, required=True)
    ap.add_argument('--sim-mskdiam', type=float, required=True)
    ap.add_argument('--real-mskdiam', type=float, required=True)
    a = ap.parse_args()
    s = measure(a.sim_stk, a.sim_cls, a.sim_mskdiam)
    r = measure(a.real_stk, a.real_cls, a.real_mskdiam)
    print('class resolution against population: sim median %.1f A at median %d particles; real median %.1f A at %d'
          % (np.median(s[:, 0]), np.median(s[:, 1]), np.median(r[:, 0]), np.median(r[:, 1])))
    print('                    frc_6_12        ridge_snr       band_ratio')
    for lo, hi in [(5, 9), (9, 12), (12, 16), (16, 25), (25, 60)]:
        x, y = s[(s[:, 0] >= lo) & (s[:, 0] < hi)], r[(r[:, 0] >= lo) & (r[:, 0] < hi)]
        f = lambda z, k: f'{np.median(z[:, k]):7.3f}' if len(z) else '      -'
        print(f'{lo:>2}-{hi:<2} A  sim n={len(x):3d} {f(x, 2)} {f(x, 3)} {f(x, 4)}   | real n={len(y):2d} {f(y, 2)} {f(y, 3)} {f(y, 4)}')


if __name__ == '__main__':
    main()
