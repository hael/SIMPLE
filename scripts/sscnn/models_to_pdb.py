#!/usr/bin/env python3
"""Write an mmCIF assembly as PDB files for SIMPLE's pdb2mrc: all heavy atoms, and the helix and strand atoms alone.

The coordinates are centred (electron-weighted) and moved to the centre of a box of --box voxels at --smpd A, so the
maps pdb2mrc makes from the three files share one frame. Secondary structure as in simulate_cavgs.py (C-alpha geometry).
Prints the model's diameter so that a box can be checked.

Usage: python -I models_to_pdb.py --cif MODEL.cif.gz --box 160 --smpd 2.0 --out DIR
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate_cavgs as sc   # noqa: E402

ELEM = {6.9: 'C', 7.7: 'N', 8.3: 'O', 16.6: 'S', 15.5: 'P', 34.6: 'SE'}


def write_pdb(path, xyz, el, comp, name):
    with open(path, 'w') as f:
        for k, (p, e) in enumerate(zip(xyz, el)):
            nm = name[k][:4]
            nm = (' ' + nm) if len(nm) < 4 and len(e) == 1 else nm
            f.write('ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n'
                    % ((k % 99999) + 1, nm, comp[k][:3], (k // 20) % 9999 + 1, p[0], p[1], p[2], e))
        f.write('END\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cif', required=True)
    ap.add_argument('--box', type=int, required=True)
    ap.add_argument('--smpd', type=float, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    xyz, w, chain, seq, name, comp = sc.read_cif(a.cif)
    xyz = xyz - np.average(xyz, axis=0, weights=w)
    ss, segs = sc.secondary_structure(xyz, chain, seq, name, comp)
    diam = 2 * np.percentile(np.linalg.norm(xyz, axis=1), 99.5)
    el = [ELEM[v] for v in w]
    xyz = xyz + a.box * a.smpd / 2.0
    write_pdb(os.path.join(a.out, 'full.pdb'), xyz, el, comp, name)
    for lab, code in (('helix', 1), ('strand', 2)):
        m = ss == code
        write_pdb(os.path.join(a.out, lab + '.pdb'), xyz[m], [e for e, k in zip(el, m) if k],
                  [c for c, k in zip(comp, m) if k], [n for n, k in zip(name, m) if k])
    print(f'{len(xyz)} heavy atoms, {(ss == 1).sum()} helix, {(ss == 2).sum()} strand; diameter {diam:.1f} A; '
          f'box {a.box * a.smpd:.0f} A')


if __name__ == '__main__':
    main()
