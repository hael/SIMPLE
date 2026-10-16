#!/usr/bin/env python3
"""aet_mat2simple.py -- convert atomic electron tomography MATLAB files (GENFIRE/RESIRE data sets) for SINGLE.

Subcommands (numpy only for MATLAB v5 files; h5py for v7.3):
  list   FILE.mat                                    the variables, their types and shapes
  stack  FILE.mat VAR OUT.mrcs --smpd S [--pad N]    3D array (d1 x d2 x nproj) -> MRC stack, MRC x = MATLAB dim 1
  vol    FILE.mat VAR OUT.mrc  --smpd S [--pad N]    3D array -> MRC volume, MRC x = MATLAB dim 1
  angles FILE.mat VAR OUT.txt                        angles -> text for import_tilt_series: a 3 x nproj (or nproj x 3)
                                                     array as 'phi theta psi' lines (angconv=genfire), a vector as one
                                                     tilt per line (angconv=tilt)
  model  FILE.mat VAR OUT.pdb  --smpd S              3 x natoms (or natoms x 3) positions -> PDB in SIMPLE's map
         [--types FILE.mat:VAR --elements Fe,Pt]     convention (A = 0-based voxel * smpd); types 1..k name the
         [--element PT] [--units px|A]               elements in order; px positions are MATLAB 1-based indices of
         [--pad-offset P] [--box-offset B]           a volume padded by P on each side (subtracted), B is added
         [--a-origin-px Z]                           for a padded target box; A positions are converted with Z, the
                                                     voxel index (1-based, padded volume) of their A origin
VAR may name a struct field as STRUCT.FIELD. With image x along the first array dimension of a projection, SINGLE's
angle conversion (angconv=genfire) puts the reconstruction in the frame of their volumes, so their models align
after the origin conversion. Example, the FePt tilt series of Yang et al. (2017), 256 -> 320 box:
  aet_mat2simple.py stack  FePt_Denoised_Projections.mat Mproj projs.mrcs --smpd 0.3725 --pad 320
  aet_mat2simple.py angles FePt_raw_projections_and_angles.mat Angles tilts.txt
  single_exec prg=import_tilt_series projfile=fept.simple stk=projs.mrcs angtab=tilts.txt angconv=tilt smpd=0.3725
"""
import argparse
import struct
import sys
import zlib

import numpy as np

# ---------------------------------------------------------------- MATLAB v5 reader (numeric, char, struct, cell)
_MI = {1: np.int8, 2: np.uint8, 3: np.int16, 4: np.uint16, 5: np.int32, 6: np.uint32,
       7: np.float32, 9: np.float64, 12: np.int64, 13: np.uint64, 16: np.uint8, 17: np.uint16}


def _tag(buf, off):
    t, n = struct.unpack_from('<II', buf, off)
    if t >> 16:                       # small data element: size in the upper half, data in the tag
        return t & 0xFFFF, t >> 16, off + 4, off + 8
    return t, n, off + 8, off + 8 + n + (8 - n % 8) % 8


def _elements(buf):
    off = 0
    while off < len(buf):
        t, n, d, nxt = _tag(buf, off)
        yield t, buf[d:d + n]
        off = nxt


def _numeric(t, raw):
    return np.frombuffer(raw, dtype=_MI[t]).copy()


def _matrix(raw):
    els = list(_elements(raw))
    flags = _numeric(*els[0])
    cls, cplx = flags[0] & 0xFF, (flags[0] >> 11) & 1
    dims = tuple(int(d) for d in _numeric(*els[1]))
    name = bytes(els[2][1]).decode('latin1') if len(els[2][1]) else ''
    if cls == 1:                      # cell
        cells = [_matrix(r)[1] for _, r in els[3:]]
        arr = np.empty(len(cells), dtype=object)
        arr[:] = cells
        return name, arr.reshape(dims, order='F') if int(np.prod(dims)) == len(cells) else arr
    if cls == 2:                      # struct
        flen = int(_numeric(*els[3])[0])
        fraw = bytes(els[4][1])
        fields = [fraw[i:i + flen].split(b'\0')[0].decode() for i in range(0, len(fraw), flen)]
        vals = [_matrix(r)[1] for _, r in els[5:]]
        nel = int(np.prod(dims))
        recs = [{f: vals[k * len(fields) + j] for j, f in enumerate(fields)} for k in range(nel)]
        return name, recs[0] if nel == 1 else recs
    if cls == 4:                      # char
        return name, ''.join(map(chr, _numeric(*els[3])))
    if cls == 5:
        raise SystemExit('sparse MATLAB arrays are not supported')
    re_ = _numeric(*els[3])
    if cplx:
        re_ = re_ + 1j * _numeric(*els[4])
    return name, re_.reshape(dims, order='F')


def _load_v5(path):
    buf = open(path, 'rb').read()
    out, off = {}, 128
    while off < len(buf):
        t, n = struct.unpack_from('<II', buf, off)
        raw = buf[off + 8:off + 8 + n]
        off += 8 + n + ((8 - n % 8) % 8 if t != 15 else 0)
        if t == 15:                   # compressed element: zlib stream holding one tagged element
            raw = zlib.decompress(raw)
            t, n2 = struct.unpack_from('<II', raw, 0)
            raw = raw[8:8 + n2]
        if t == 14:
            name, val = _matrix(raw)
            out[name] = val
    return out


def _load_v73(path):
    try:
        import h5py
    except ImportError:
        raise SystemExit('MATLAB v7.3 files are HDF5: install h5py (pip install h5py)')
    out = {}
    with h5py.File(path, 'r') as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and not name.startswith('#'):
                # HDF5 stores MATLAB's column-major array with the dimensions reversed
                out[name.replace('/', '.')] = np.array(obj).T
        f.visititems(visit)
    return out


def loadmat(path):
    head = open(path, 'rb').read(128)
    if head.startswith(b'\x89HDF') or b'MATLAB 7.3' in head:
        return _load_v73(path)
    return _load_v5(path)


def get_var(path, var):
    d = loadmat(path)
    if var in d:
        return d[var]
    val = d
    for key in var.split('.'):
        if not isinstance(val, dict) or key not in val:
            raise SystemExit(f'{var} is not in {path}; variables: {", ".join(d)}')
        val = val[key]
    return val


# ---------------------------------------------------------------- writers
def write_mrc(fn, a, smpd, ispg):
    """a indexed (x, y, z); x fastest on disk"""
    a = np.asarray(a, np.float32)
    nx, ny, nz = a.shape
    h = bytearray(1024)
    struct.pack_into('<10i', h, 0, nx, ny, nz, 2, 0, 0, 0, nx, ny, nz)
    struct.pack_into('<6f', h, 40, nx * smpd, ny * smpd, nz * smpd, 90., 90., 90.)
    struct.pack_into('<3i', h, 64, 1, 2, 3)
    struct.pack_into('<3f', h, 76, float(a.min()), float(a.max()), float(a.mean()))
    struct.pack_into('<i', h, 88, ispg)
    h[208:212] = b'MAP '
    h[212:216] = bytes([0x44, 0x44, 0, 0])
    struct.pack_into('<f', h, 216, float(a.std()))
    with open(fn, 'wb') as f:
        f.write(h)
        f.write(np.asarray(a, '<f4').ravel(order='F').tobytes())


def pad_xy(a, n, depth_too):
    """zero-pad the first two (and with depth_too the third) dimensions to n, centred as SIMPLE's n/2+1 origin"""
    shape = list(a.shape)
    target = [n, n, n if depth_too else shape[2]]
    if any(t < s for t, s in zip(target, shape)):
        raise SystemExit(f'--pad {n} is smaller than the array {tuple(shape)}')
    out = np.zeros(target, np.float32)
    o = [(t - s) // 2 for t, s in zip(target, shape)]
    out[o[0]:o[0] + shape[0], o[1]:o[1] + shape[1], o[2]:o[2] + shape[2]] = a
    return out


def as_3d(v, what):
    a = np.asarray(v)
    if a.dtype == object or a.ndim != 3:
        raise SystemExit(f'{what} needs a 3D numeric array, got shape {getattr(a, "shape", None)}')
    return np.real(a).astype(np.float32)


# ---------------------------------------------------------------- subcommands
def cmd_list(a):
    for k, v in loadmat(a.mat).items():
        if isinstance(v, np.ndarray) and v.dtype != object:
            rng = f'[{float(np.nanmin(v.real)):.4g}, {float(np.nanmax(v.real)):.4g}]' if v.size else ''
            print(f'{k:32s} {str(v.dtype):8s} {str(v.shape):20s} {rng}')
        elif isinstance(v, dict):
            print(f'{k:32s} struct   fields: {", ".join(v)}')
        else:
            print(f'{k:32s} {type(v).__name__}')


def cmd_stack(a):
    s = as_3d(get_var(a.mat, a.var), 'stack')
    if a.pad:
        s = pad_xy(s, a.pad, False)
    write_mrc(a.out, s, a.smpd, 0)
    print(f'wrote {s.shape[2]} images of {s.shape[0]} x {s.shape[1]} to {a.out} (MRC x = MATLAB dim 1)')


def cmd_vol(a):
    v = as_3d(get_var(a.mat, a.var), 'vol')
    if a.pad:
        v = pad_xy(v, a.pad, True)
    write_mrc(a.out, v, a.smpd, 1)
    print(f'wrote a {v.shape} volume to {a.out} (MRC x = MATLAB dim 1)')


def cmd_angles(a):
    ang = np.real(np.asarray(get_var(a.mat, a.var), float))
    if ang.ndim == 2 and 3 in ang.shape and min(ang.shape) > 1:
        tri = ang if ang.shape[0] == 3 else ang.T
        if ang.shape == (3, 3):
            print('3 x 3 angle array: read as columns = projections (the GENFIRE layout)')
        with open(a.out, 'w') as f:
            f.write('# phi theta psi (degrees, GENFIRE/RESIRE convention), stack order: angconv=genfire\n')
            for p in tri.T:
                f.write(f'{p[0]:12.5f} {p[1]:12.5f} {p[2]:12.5f}\n')
        print(f'wrote {tri.shape[1]} angle triplets to {a.out}; import with angconv=genfire')
    elif ang.size == max(ang.shape):
        with open(a.out, 'w') as f:
            f.write('# single-axis tilt angle (degrees) about image x, stack order: angconv=tilt\n')
            for t in ang.ravel():
                f.write(f'{t:12.5f}\n')
        print(f'wrote {ang.size} tilt angles to {a.out}; import with angconv=tilt')
    else:
        raise SystemExit(f'cannot read angles of shape {ang.shape}: need 3 x nproj, nproj x 3 or a vector')


def cmd_model(a):
    pos = np.real(np.asarray(get_var(a.mat, a.var), float))
    if pos.ndim != 2 or 3 not in pos.shape:
        raise SystemExit(f'positions need a 3 x natoms or natoms x 3 array, got {pos.shape}')
    xyz = pos.T if pos.shape[0] == 3 and pos.shape[1] != 3 else pos
    if a.units == 'A':
        px = xyz / a.smpd + a.a_origin_px
    else:
        px = xyz
    ang = (px - 1. - a.pad_offset + a.box_offset) * a.smpd
    if a.types:
        fn, var = a.types.rsplit(':', 1)
        typ = np.real(np.asarray(get_var(fn, var))).ravel().astype(int)
        if typ.size != ang.shape[0]:
            raise SystemExit(f'{typ.size} types for {ang.shape[0]} atoms')
        els = [e.strip().upper() for e in a.elements.split(',')] if a.elements else []
        if not els or typ.min() < 1 or typ.max() > len(els):
            raise SystemExit(f'types run {typ.min()}..{typ.max()}: give --elements with that many names, in order')
        names = [els[t - 1] for t in typ]
    else:
        names = [a.element.upper()] * ang.shape[0]
    with open(a.out, 'w') as f:
        f.write(f'REMARK converted by aet_mat2simple.py from {a.mat}:{a.var}; A = (0-based voxel) * {a.smpd}\n')
        for i, (p, nm) in enumerate(zip(ang, names)):
            f.write('ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'
                    % ((i + 1) % 100000, ' ' + nm, nm, 'A', 1, p[0], p[1], p[2], 1.0, 0.0, nm))
        f.write('END\n')
    print(f'wrote {ang.shape[0]} atoms to {a.out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('list');   p.add_argument('mat')
    for name in ('stack', 'vol'):
        p = sub.add_parser(name)
        p.add_argument('mat'); p.add_argument('var'); p.add_argument('out')
        p.add_argument('--smpd', type=float, required=True); p.add_argument('--pad', type=int, default=0)
    p = sub.add_parser('angles'); p.add_argument('mat'); p.add_argument('var'); p.add_argument('out')
    p = sub.add_parser('model')
    p.add_argument('mat'); p.add_argument('var'); p.add_argument('out')
    p.add_argument('--smpd', type=float, required=True)
    p.add_argument('--types', default=None, help='FILE.mat:VAR of integer types 1..k')
    p.add_argument('--elements', default=None, help='element names of types 1..k, comma separated')
    p.add_argument('--element', default='X', help='element name of every atom without --types')
    p.add_argument('--units', choices=('px', 'A'), default='px')
    p.add_argument('--pad-offset', dest='pad_offset', type=float, default=0.)
    p.add_argument('--box-offset', dest='box_offset', type=float, default=0.)
    p.add_argument('--a-origin-px', dest='a_origin_px', type=float, default=1.)
    a = ap.parse_args()
    {'list': cmd_list, 'stack': cmd_stack, 'vol': cmd_vol, 'angles': cmd_angles, 'model': cmd_model}[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
