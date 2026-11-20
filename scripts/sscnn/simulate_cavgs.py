#!/usr/bin/env python3
"""Simulated 2D class averages from atomic models, with ground truth for secondary-structure visibility.

For each model (an mmCIF file, the biological assembly) and each simulated class:
  view        a uniformly random 3D orientation
  density     heavy atoms projected with their electron counts (hydrogens folded in), minus the water they displace
              (0.334 e/A^3 over 17.1 A^3 per heavy atom, spread wider), so that low-resolution contrast is protein
              minus solvent, as in ice
  degradation a B-factor envelope exp(-B f^2 / 4), a translational alignment blur, an in-plane rotational blur (an
              average of jittered views), then noise: white, gently low-passed at 4.5 A, scaled so that the halves'
              FRC is 0.5 at the class resolution d_frc (signal and half-noise power equal in that shell), the
              resolution SIMPLE reports for a class
  halves      even and odd carry the same signal and independent noise of twice the average's variance; the average
              is their mean, as a class average and its halves relate
  truth       noise-free maps of the helix atoms and the strand atoms through the same envelope and blur; per class
              the helices' contrast-to-noise over 12-6 A (cnr_helix: their band power over the average's noise power in
              the band), the fraction of helix residues lying within 45 degrees of the image plane, and a label
              visible = cnr_helix >= CNR_VISIBLE (a guessed threshold, kept as metadata so another can be applied)
Secondary structure comes from the C-alpha geometry, the same for every model: helix where d(i,i+3) is 4.6-5.6 A and
d(i,i+4) 5.7-6.8 A (runs of 6 residues or more); strand where the chain is extended (d(i-1,i+1) >= 6.0 A) and an
extended C-alpha of another strand lies 4.2-5.5 A away.

Output in --out: shards <model>.npz (float16 avg, even, odd, helix, strand; metadata arrays), meta.tsv, and a preview
(preview.mrc with preview_even.mrc and preview_odd.mrc, preview_meta.tsv) of --preview classes across all models, which
ss_trainer.py can open.

Usage: python -I simulate_cavgs.py --models DIR --out DIR [--per-model 50] [--box 160] [--smpd 2.0] [--seed 1]
"""
import argparse
import gzip
import os
import struct
import sys
import time

import numpy as np

ZEFF = {'C': 6.9, 'N': 7.7, 'O': 8.3, 'S': 16.6, 'P': 15.5, 'SE': 34.6}
SOLV_E = 0.334 * 17.1        # electrons of water displaced per heavy atom
SIG_ATOM, SIG_SOLV = 0.5, 1.8  # A
DFRC_RANGE = (5.0, 40.0)     # class resolution, the halves' FRC = 0.5 (A), log-uniform
B_RANGE = (50.0, 400.0)      # signal B-factor (A^2), uniform
LP_NOISE = 4.5               # the noise's soft low-pass (A)
TRANS_MAX, ROT_MAX = 2.0, 3.0  # alignment blur: translation sd (A), rotation sd (deg), uniform from 0
NROT_BLUR = 6
CNR_VISIBLE = 1.0
BAND = (1 / 12.0, 1 / 6.0)
AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
      'THR', 'TRP', 'TYR', 'VAL', 'MSE', 'SEC', 'PYL'}


# ------------------------------------------------------------------ model
def read_cif(path):
    """Heavy atoms of ATOM records: coordinates, electron counts, chain keys, residue numbers, atom names, residue names."""
    op = gzip.open if path.endswith('.gz') else open
    cols, rows, in_loop, in_atoms = [], [], False, False
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('loop_'):
                in_loop, in_atoms, cols = True, False, []
                continue
            if in_loop and line.startswith('_atom_site.'):
                cols.append(line.split('.', 1)[1].strip())
                in_atoms = True
                continue
            if in_atoms:
                if line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
                    if rows:
                        break
                    in_atoms = False
                    continue
                t = line.split()
                if len(t) == len(cols):
                    rows.append(t)
    c = {k: i for i, k in enumerate(cols)}
    need = ['group_PDB', 'type_symbol', 'label_atom_id', 'label_comp_id', 'label_asym_id', 'label_seq_id',
            'Cartn_x', 'Cartn_y', 'Cartn_z']
    for k in need:
        if k not in c:
            raise SystemExit(f'{path}: no _atom_site.{k}')
    model = c.get('pdbx_PDB_model_num')
    auth = c.get('auth_asym_id')
    xyz, w, chain, seq, name, comp = [], [], [], [], [], []
    for t in rows:
        if t[c['group_PDB']] != 'ATOM':
            continue
        el = t[c['type_symbol']].upper()
        if el in ('H', 'D') or el not in ZEFF:
            continue
        xyz.append((float(t[c['Cartn_x']]), float(t[c['Cartn_y']]), float(t[c['Cartn_z']])))
        w.append(ZEFF[el])
        chain.append((t[c['label_asym_id']], t[auth] if auth is not None else '', t[model] if model is not None else ''))
        s = t[c['label_seq_id']]
        seq.append(int(s) if s.lstrip('-').isdigit() else -999999)
        name.append(t[c['label_atom_id']].strip('"'))
        comp.append(t[c['label_comp_id']])
    return np.array(xyz), np.array(w), chain, np.array(seq), name, comp


def secondary_structure(xyz, chain, seq, name, comp):
    """Per atom: 1 helix, 2 strand, 0 other; and the helix segments as lists of C-alpha coordinates."""
    ca = [i for i in range(len(name)) if name[i] == 'CA' and comp[i] in AA]
    by_chain = {}
    for i in ca:
        by_chain.setdefault(chain[i], []).append(i)
    res_ss = {}           # (chain, seq) -> 1 / 2
    segments = []
    ext_pts, ext_keys, ext_idx = [], [], []
    for ch, idx in by_chain.items():
        idx.sort(key=lambda i: seq[i])
        s = seq[idx]
        p = xyz[idx]
        n = len(idx)
        hel = np.zeros(n, bool)
        for q in range(n - 4):
            if s[q + 4] - s[q] != 4:
                continue
            d3 = np.linalg.norm(p[q + 3] - p[q])
            d4 = np.linalg.norm(p[q + 4] - p[q])
            if 4.6 <= d3 <= 5.6 and 5.7 <= d4 <= 6.8:
                hel[q:q + 5] = True
        q = 0
        while q < n:                                   # keep helices of six residues or more
            if hel[q]:
                r = q
                while r + 1 < n and hel[r + 1] and s[r + 1] == s[r] + 1:
                    r += 1
                if r - q + 1 >= 6:
                    for k in range(q, r + 1):
                        res_ss[(ch, s[k])] = 1
                    segments.append(p[q:r + 1])
                q = r + 1
            else:
                q += 1
        for q in range(1, n - 1):
            if s[q + 1] - s[q - 1] == 2 and np.linalg.norm(p[q + 1] - p[q - 1]) >= 6.0 and (ch, s[q]) not in res_ss:
                ext_pts.append(p[q])
                ext_keys.append((ch, s[q]))
                ext_idx.append(q)
    if ext_pts:                                        # strands: extended residues with an extended partner nearby
        P = np.array(ext_pts)
        cell = np.floor(P / 6.0).astype(int)
        grid = {}
        for k, cc in enumerate(map(tuple, cell)):
            grid.setdefault(cc, []).append(k)
        for k, cc in enumerate(map(tuple, cell)):
            found = False
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for m in grid.get((cc[0] + dx, cc[1] + dy, cc[2] + dz), ()):
                            if m == k:
                                continue
                            same = ext_keys[m][0] == ext_keys[k][0] and abs(ext_keys[m][1] - ext_keys[k][1]) < 3
                            if not same and 4.2 <= np.linalg.norm(P[m] - P[k]) <= 5.5:
                                found = True
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                res_ss[ext_keys[k]] = 2
    ss = np.array([res_ss.get((chain[i], seq[i]), 0) for i in range(len(name))], dtype=np.int8)
    return ss, segments


# ------------------------------------------------------------------ simulation
def random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def rot_z(deg):
    a = np.radians(deg)
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


class Projector:
    def __init__(self, box, smpd):
        self.box, self.smpd = box, smpd
        fx = np.fft.fftfreq(box, d=smpd)
        FX, FY = np.meshgrid(fx, fx)
        self.f2 = FX ** 2 + FY ** 2
        # the bilinear splat's transfer, divided out
        sx = np.sinc(FX * smpd) ** 2
        sy = np.sinc(FY * smpd) ** 2
        self.unsplat = 1.0 / np.maximum(sx * sy, 0.05)
        self.atom = np.exp(-2 * np.pi ** 2 * SIG_ATOM ** 2 * self.f2)
        self.solv = SOLV_E * np.exp(-2 * np.pi ** 2 * SIG_SOLV ** 2 * self.f2)
        f = np.sqrt(self.f2)
        self.band = (f >= BAND[0]) & (f <= BAND[1])
        c = box // 2
        y, x = np.mgrid[0:box, 0:box] - c
        self.r = np.hypot(x, y) * smpd

    def splat(self, xy, w):
        g = xy / self.smpd + self.box // 2
        x0, y0 = np.floor(g[:, 0]).astype(int), np.floor(g[:, 1]).astype(int)
        fx, fy = g[:, 0] - x0, g[:, 1] - y0
        img = np.zeros((self.box, self.box))
        for dx, dy, ww in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)), (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
            xi, yi = x0 + dx, y0 + dy
            ok = (xi >= 0) & (xi < self.box) & (yi >= 0) & (yi < self.box)
            np.add.at(img, (yi[ok], xi[ok]), (w * ww)[ok])
        return img

    def project_ft(self, xyz, w, Rs):
        """Fourier transform of the mean projection over the rotations Rs: atoms minus displaced water."""
        acc_a = np.zeros((self.box, self.box))
        acc_s = np.zeros((self.box, self.box))
        for R in Rs:
            p = xyz @ R.T
            acc_a += self.splat(p[:, :2], w)
            acc_s += self.splat(p[:, :2], np.ones(len(w)))
        n = len(Rs)
        return (np.fft.fft2(acc_a / n) * self.atom - np.fft.fft2(acc_s / n) * self.solv) * self.unsplat


def write_mrc(path, stack, smpd):
    stack = np.asarray(stack, dtype='<f4')
    nz, ny, nx = stack.shape
    h = bytearray(1024)
    struct.pack_into('<3i', h, 0, nx, ny, nz)
    struct.pack_into('<i', h, 12, 2)
    struct.pack_into('<3i', h, 28, nx, ny, nz)
    struct.pack_into('<3f', h, 40, nx * smpd, ny * smpd, nz * smpd)
    struct.pack_into('<3f', h, 52, 90.0, 90.0, 90.0)
    struct.pack_into('<3i', h, 64, 1, 2, 3)
    struct.pack_into('<3f', h, 76, float(stack.min()), float(stack.max()), float(stack.mean()))
    h[208:212] = b'MAP '
    h[212:216] = bytes([0x44, 0x44, 0, 0])
    struct.pack_into('<f', h, 216, float(stack.std()))
    with open(path, 'wb') as f:
        f.write(bytes(h))
        f.write(stack.tobytes())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--models', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--per-model', type=int, default=50)
    ap.add_argument('--box', type=int, default=160)
    ap.add_argument('--smpd', type=float, default=2.0)
    ap.add_argument('--preview', type=int, default=120)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    pr = Projector(a.box, a.smpd)
    files = sorted(f for f in os.listdir(a.models) if f.endswith('.cif.gz') or f.endswith('.cif'))
    meta_rows, preview = [], []
    hdr = ['model', 'index', 'd_frc_A', 'B', 'trans_sd_A', 'rot_sd_deg', 'snr', 'helix_atom_frac',
           'strand_atom_frac', 'sideon_frac', 'cnr_helix', 'cnr_strand', 'visible', 'diam_A']
    for fn in files:
        t0 = time.time()
        name = fn.split('-')[0].split('.')[0]
        xyz, w, chain, seq, aname, comp = read_cif(os.path.join(a.models, fn))
        if len(xyz) == 0:
            print(f'{name}: no atoms, skipped')
            continue
        xyz = xyz - np.average(xyz, axis=0, weights=w)
        ss, segs = secondary_structure(xyz, chain, seq, aname, comp)
        diam = 2 * np.percentile(np.linalg.norm(xyz, axis=1), 99.5)
        if diam > 0.95 * a.box * a.smpd:
            print(f'{name}: {diam:.0f} A across does not fit the {a.box * a.smpd:.0f} A box, skipped')
            continue
        hfrac, sfrac = (ss == 1).mean(), (ss == 2).mean()
        axes, lens = [], []
        for sg in segs:
            u = np.linalg.svd(sg - sg.mean(0), full_matrices=False)[2][0]
            axes.append(u)
            lens.append(len(sg))
        axes, lens = np.array(axes), np.array(lens)
        rmask = pr.r <= diam / 2 + 10.0
        n = a.per_model
        out = {k: np.zeros((n, a.box, a.box), np.float16) for k in ('avg', 'even', 'odd', 'helix', 'strand')}
        M = {k: np.zeros(n) for k in hdr[2:]}
        for i in range(n):
            R = random_rotation(rng)
            d = np.exp(rng.uniform(*np.log(DFRC_RANGE)))
            B = rng.uniform(*B_RANGE)
            ts = rng.uniform(0, TRANS_MAX)
            rs = rng.uniform(0, ROT_MAX)
            Rs = [rot_z(rng.normal(0, rs)) @ R for _ in range(NROT_BLUR)] if rs > 0.05 else [R]
            env = np.exp(-B * pr.f2 / 4.0) * np.exp(-2 * np.pi ** 2 * ts ** 2 * pr.f2)
            S = pr.project_ft(xyz, w, Rs) * env
            H = pr.project_ft(xyz[ss == 1], w[ss == 1], Rs) * env if (ss == 1).any() else np.zeros_like(S)
            T = pr.project_ft(xyz[ss == 2], w[ss == 2], Rs) * env if (ss == 2).any() else np.zeros_like(S)
            sig = np.real(np.fft.ifft2(S))
            # noise: white, softly low-passed; per half scaled so that its power equals the signal's in the shell 1/d
            nf = 1.0 / (1.0 + (np.sqrt(pr.f2) * LP_NOISE) ** 8)
            Ne = np.fft.fft2(rng.normal(size=sig.shape)) * nf
            No = np.fft.fft2(rng.normal(size=sig.shape)) * nf
            shell = np.abs(np.sqrt(pr.f2) - 1.0 / d) <= 1.0 / (a.box * a.smpd)
            scale = np.sqrt((np.abs(S[shell]) ** 2).mean() / max((np.abs(Ne[shell]) ** 2).mean(), 1e-30))
            ne = np.real(np.fft.ifft2(Ne)) * scale
            no = np.real(np.fft.ifft2(No)) * scale
            navg = 0.5 * (ne + no)
            snr = sig[rmask].var() / navg[rmask].var()
            Np = np.abs(np.fft.fft2(navg)) ** 2
            nb = Np[pr.band].sum()
            out['avg'][i] = sig + navg
            out['even'][i] = sig + ne
            out['odd'][i] = sig + no
            out['helix'][i] = np.real(np.fft.ifft2(H))
            out['strand'][i] = np.real(np.fft.ifft2(T))
            M['d_frc_A'][i], M['B'][i], M['trans_sd_A'][i], M['rot_sd_deg'][i] = d, B, ts, rs
            M['snr'][i] = snr
            M['helix_atom_frac'][i], M['strand_atom_frac'][i] = hfrac, sfrac
            M['sideon_frac'][i] = (lens * (np.abs(axes @ R.T)[:, 2] < np.sin(np.radians(45)))).sum() / lens.sum() if len(lens) else 0.0
            M['cnr_helix'][i] = np.sqrt((np.abs(H[pr.band]) ** 2).sum() / nb) if nb > 0 else 0.0
            M['cnr_strand'][i] = np.sqrt((np.abs(T[pr.band]) ** 2).sum() / nb) if nb > 0 else 0.0
            M['visible'][i] = float(M['cnr_helix'][i] >= CNR_VISIBLE)
            M['diam_A'][i] = diam
            meta_rows.append([name, i] + [M[k][i] for k in hdr[2:]])
            preview.append((name, i))
        np.savez_compressed(os.path.join(a.out, name + '.npz'), smpd=a.smpd, **out, **M)
        print(f'{name}: {len(xyz)} heavy atoms, {diam:.0f} A, helix {100 * hfrac:.0f}% strand {100 * sfrac:.0f}% of atoms, '
              f'{len(segs)} helices; {n} classes, {int(M["visible"].sum())} visible; {time.time() - t0:.0f} s')
        sys.stdout.flush()
    with open(os.path.join(a.out, 'meta.tsv'), 'w') as f:
        f.write('\t'.join(hdr) + '\n')
        for r in meta_rows:
            f.write('\t'.join(str(v) if isinstance(v, str) else (f'{v:.6g}' if isinstance(v, float) else str(v)) for v in r) + '\n')
    # the preview: a random draw across models, as MRC stacks the trainer opens
    pick = rng.choice(len(preview), size=min(a.preview, len(preview)), replace=False)
    shards = {}
    stacks = {k: [] for k in ('avg', 'even', 'odd')}
    with open(os.path.join(a.out, 'preview_meta.tsv'), 'w') as f:
        f.write('class\t' + '\t'.join(hdr) + '\n')
        for c, k in enumerate(pick):
            name, i = preview[k]
            if name not in shards:
                shards[name] = np.load(os.path.join(a.out, name + '.npz'))
            for s in stacks:
                stacks[s].append(shards[name][s][i].astype(np.float32))
            row = meta_rows[k]
            f.write(f'{c + 1}\t' + '\t'.join(str(v) if isinstance(v, str) else (f'{v:.6g}' if isinstance(v, float) else str(v)) for v in row) + '\n')
    write_mrc(os.path.join(a.out, 'preview.mrc'), stacks['avg'], a.smpd)
    write_mrc(os.path.join(a.out, 'preview_even.mrc'), stacks['even'], a.smpd)
    write_mrc(os.path.join(a.out, 'preview_odd.mrc'), stacks['odd'], a.smpd)
    print(f'wrote {len(meta_rows)} classes and a preview of {len(pick)} to {a.out}')


if __name__ == '__main__':
    main()
