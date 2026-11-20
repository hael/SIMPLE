#!/usr/bin/env python3
"""Prepare an EMDB map for build_sim_set.py: the experimental density and its helix part, in the batch's frame.

  density   the mean of the two half maps (unsharpened, unmasked), or the primary map without them; its background
            (the mean outside a sphere of the source box) subtracted, a soft sphere just inside the source box applied
            (so the box's edges do not project), low-passed to the target Nyquist and resampled (trilinear) onto a
            --box cube at --smpd, centred on the fitted model's electron-weighted centre, as models_to_pdb.py centres
  helix     the density within HELIX_RADIUS of the model's helix atoms (soft edge): real helix density, not a model's
  strand    the same for the strand atoms
  check     the correlation, at 8 A and inside the model's sphere, of the density with a density computed from the model:
            a frame error (origin, axis order) shows as a low value; also the share of the density outside the model,
            which is the micelle or nanodisc (and noise)
Writes full.mrc, helix.mrc, strand.mrc and model.json (diameter for the mask, the checks) into --out, the layout build_sim_set.py
reads with --maps-from. Map and model are read as data; run with python -I.

Usage: python -I emdb_prep.py --model MODEL.cif.gz (--half1 H1 --half2 H2 | --map MAP) --out DIR [--box 160] [--smpd 2.0]
"""
import argparse
import gzip
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate_cavgs as sc   # noqa: E402

HELIX_RADIUS = 3.0     # A
SIG_MODEL = 1.5        # A, the model density's atomic blur for the frame check
LP_CHECK = 8.0         # A


def read_ccp4(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rb') as f:
        h = f.read(1024)
        nc, nr, ns, mode = struct.unpack('<4i', h[:16])
        start = struct.unpack('<3i', h[16:28])
        mxyz = struct.unpack('<3i', h[28:40])
        cell = struct.unpack('<3f', h[40:52])
        axes = struct.unpack('<3i', h[64:76])
        nsym = struct.unpack('<i', h[92:96])[0]
        origin = struct.unpack('<3f', h[196:208])
        if mode != 2:
            raise SystemExit(f'{path}: mode {mode}, expected 2')
        f.read(nsym)
        a = np.frombuffer(f.read(nc * nr * ns * 4), dtype='<f4').reshape(ns, nr, nc).astype(np.float64)
    # to [z][y][x]: the array's axes are (section, row, column) = axes[2], axes[1], axes[0]
    order = [3 - axes[2], 3 - axes[1], 3 - axes[0]]          # positions of z, y, x among (sec,row,col) as 0..2
    perm = [order.index(k) for k in (0, 1, 2)]
    a = np.transpose(a, perm)
    nxyz = (nc, nr, ns)
    n_by_axis = {axes[0]: nc, axes[1]: nr, axes[2]: ns}
    s_by_axis = {axes[0]: start[0], axes[1]: start[1], axes[2]: start[2]}
    apix = np.array([cell[0] / mxyz[0], cell[1] / mxyz[1], cell[2] / mxyz[2]])
    st = np.array([s_by_axis[1], s_by_axis[2], s_by_axis[3]], float)
    return a, apix, st, np.array(origin, float), np.array([n_by_axis[1], n_by_axis[2], n_by_axis[3]]), nxyz


def lowpass(v, apix, res):
    F = np.fft.rfftn(v)
    fz = np.fft.fftfreq(v.shape[0], apix[2])
    fy = np.fft.fftfreq(v.shape[1], apix[1])
    fx = np.fft.rfftfreq(v.shape[2], apix[0])
    f = np.sqrt(fz[:, None, None] ** 2 + fy[None, :, None] ** 2 + fx[None, None, :] ** 2)
    edge = 1.0 / res
    filt = np.clip((1.1 * edge - f) / (0.2 * edge), 0, 1)
    return np.fft.irfftn(F * filt, s=v.shape, axes=(0, 1, 2))


def trilinear(v, idx):
    """v[z][y][x] at fractional indices idx (3, N) ordered (x, y, z); zero outside."""
    x, y, z = idx
    x0, y0, z0 = np.floor(x).astype(int), np.floor(y).astype(int), np.floor(z).astype(int)
    out = np.zeros(x.shape)
    nz, ny, nx = v.shape
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                xi, yi, zi = x0 + dx, y0 + dy, z0 + dz
                w = (1 - abs(x - xi)) * (1 - abs(y - yi)) * (1 - abs(z - zi))
                ok = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny) & (zi >= 0) & (zi < nz)
                out[ok] += w[ok] * v[zi[ok], yi[ok], xi[ok]]
    return out


def splat_blur(pts, w, box, smpd, sigma):
    """Atoms (Angstrom, box frame) splatted to the nearest voxel, blurred by a Gaussian of sigma (A)."""
    g = np.rint(pts / smpd).astype(int)
    ok = np.all((g >= 0) & (g < box), axis=1)
    v = np.zeros((box, box, box))
    np.add.at(v, (g[ok, 2], g[ok, 1], g[ok, 0]), w[ok])
    f = np.fft.fftfreq(box, smpd)
    fr = np.fft.rfftfreq(box, smpd)
    f2 = f[:, None, None] ** 2 + f[None, :, None] ** 2 + fr[None, None, :] ** 2
    return np.fft.irfftn(np.fft.rfftn(v) * np.exp(-2 * np.pi ** 2 * sigma ** 2 * f2), s=v.shape, axes=(0, 1, 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--model', required=True)
    ap.add_argument('--half1')
    ap.add_argument('--half2')
    ap.add_argument('--map')
    ap.add_argument('--out', required=True)
    ap.add_argument('--box', type=int, default=160)
    ap.add_argument('--smpd', type=float, default=2.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.half1 and a.half2:
        v1, apix, st, org, nn, _ = read_ccp4(a.half1)
        v2, apix2, st2, org2, _, _ = read_ccp4(a.half2)
        if v1.shape != v2.shape or np.abs(apix - apix2).max() > 1e-4 or np.abs(st - st2).max() > 0 or np.abs(org - org2).max() > 1e-3:
            raise SystemExit('the half maps differ in grid')
        v, source = 0.5 * (v1 + v2), 'mean of the half maps'
    elif a.map:
        v, apix, st, org, nn, _ = read_ccp4(a.map)
        source = 'primary map'
    else:
        raise SystemExit('give --half1 and --half2, or --map')
    if np.abs(apix - apix[0]).max() > 1e-4:
        raise SystemExit('anisotropic voxels are not handled')
    # background and the soft sphere inside the source box
    c_src = (np.array(v.shape[::-1]) - 1) / 2.0                      # x, y, z voxel centre
    zz, yy, xx = np.ogrid[0:v.shape[0], 0:v.shape[1], 0:v.shape[2]]
    r_src = np.sqrt((xx - c_src[0]) ** 2 + (yy - c_src[1]) ** 2 + (zz - c_src[2]) ** 2) * apix[0]
    rmax = (min(v.shape) / 2.0) * apix[0] - 5.0
    v = v - v[r_src > rmax].mean()
    v *= np.clip((rmax - r_src) / 8.0, 0, 1)
    v = lowpass(v, apix, 2.0 * a.smpd)
    # the model: centre and secondary structure as the atomic-model route
    xyz, w, chain, seq, name, comp = sc.read_cif(a.model)
    ss, segs = sc.secondary_structure(xyz, chain, seq, name, comp)
    com = np.average(xyz, axis=0, weights=w)
    diam = 2 * np.percentile(np.linalg.norm(xyz - com, axis=1), 99.5)
    # the target grid: voxel k <-> model coordinate com + (k - box/2) * smpd
    k = (np.arange(a.box) - a.box // 2) * a.smpd
    Z, Y, X = np.meshgrid(k, k, k, indexing='ij')
    P = np.stack([X.ravel() + com[0], Y.ravel() + com[1], Z.ravel() + com[2]])
    idx = (P - org[:, None]) / apix[:, None] - st[:, None]
    full = trilinear(v, idx).reshape(a.box, a.box, a.box)
    # the helix part: density within HELIX_RADIUS of a helix atom
    pts = xyz - com + a.box // 2 * a.smpd
    hel = splat_blur(pts[ss == 1], np.ones((ss == 1).sum()), a.box, a.smpd, 1.0)
    one = splat_blur(np.array([[a.box // 2 * a.smpd] * 3]), np.ones(1), a.box, a.smpd, 1.0)
    c0 = a.box // 2
    prof = one[c0, c0, c0:]                                           # a single atom's blurred profile along x
    thr = np.interp(HELIX_RADIUS / a.smpd, np.arange(len(prof)), prof)
    mask = np.clip((hel - 0.5 * thr) / (0.5 * thr), 0, 1)
    helix = full * mask
    if (ss == 2).any():
        stra = splat_blur(pts[ss == 2], np.ones((ss == 2).sum()), a.box, a.smpd, 1.0)
        strand = full * np.clip((stra - 0.5 * thr) / (0.5 * thr), 0, 1)
    else:
        strand = np.zeros_like(full)
    # the frame check
    model_d = splat_blur(pts, w, a.box, a.smpd, SIG_MODEL)
    fa, fb = lowpass(full, [a.smpd] * 3, LP_CHECK), lowpass(model_d, [a.smpd] * 3, LP_CHECK)
    rr = np.sqrt(sum((g - c0) ** 2 for g in np.indices(full.shape))) * a.smpd
    sph = rr <= diam / 2
    x1, x2 = fa[sph] - fa[sph].mean(), fb[sph] - fb[sph].mean()
    cc = float((x1 * x2).sum() / np.sqrt((x1 * x1).sum() * (x2 * x2).sum()))
    env = splat_blur(pts, np.ones(len(pts)), a.box, a.smpd, 3.0) > thr * 0.2
    pos = np.clip(fa, 0, None)
    outside = float(pos[~env].sum() / pos.sum())
    thr_map = fa.mean() + 3 * fa.std()
    rad = rr[fa > thr_map]
    extent = float(2 * np.percentile(rad, 99.5)) if rad.size else diam
    sc.write_mrc(os.path.join(a.out, 'full.mrc'), full.astype(np.float32), a.smpd)
    sc.write_mrc(os.path.join(a.out, 'helix.mrc'), helix.astype(np.float32), a.smpd)
    sc.write_mrc(os.path.join(a.out, 'strand.mrc'), strand.astype(np.float32), a.smpd)
    info = (f'{source}; voxel {apix[0]:.4f} A, source box {v.shape[2]}; model {len(xyz)} heavy atoms, '
            f'{(ss == 1).sum()} helix, diameter {diam:.1f} A; map extent {extent:.1f} A; '
            f'map-model cc at {LP_CHECK:.0f} A {cc:.3f}; density outside the model {100 * outside:.0f}%; '
            f'helix share of the density {100 * np.clip(helix, 0, None).sum() / np.clip(full, 0, None).sum():.0f}%, '
            f'strand share {100 * np.clip(strand, 0, None).sum() / np.clip(full, 0, None).sum():.0f}%')
    json.dump({'diam': max(diam, extent), 'model_diam': diam, 'map_extent': extent, 'cc_model': cc,
               'outside_frac': outside, 'info': info}, open(os.path.join(a.out, 'model.json'), 'w'), indent=1)
    print(info)


if __name__ == '__main__':
    main()
