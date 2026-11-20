#!/usr/bin/env python3
"""The class-average CNN for SIMPLE's Fortran inference (simple_cavg_cnn, program cavg_cnn_predict).

  export_fortran.py --model MODEL_DIR/model.pt --out cavg_cnn_weights.bin
  export_fortran.py --fixture DIR
  export_fortran.py --dump-inputs --model MODEL_DIR/model.pt --stk cavgs.mrc --out inputs.bin

Weights file (little-endian, stream): int64 magic and version; int32 width, input channels, box, GroupNorm groups,
blocks; real32 pixel size; then every tensor in forward order (stem convolution, its GroupNorm weight and bias; per
block c1, n1 weight and bias, c2, n2 weight and bias, the skip's 1x1 convolution and its GroupNorm weight and bias;
the head's weight and bias), each as int32 rank, int32 extents in Fortran order (PyTorch's shape reversed) and real32
values in PyTorch's memory order, which is the same array read column-major: a Fortran (w, h, c) image is PyTorch's
(c, h, w).
Inputs file: int32 count, then per input int32 nx, ny, channels and real32 values in the same order.
--fixture writes DIR/weights.raw, DIR/inputs.raw and DIR/logits.txt (.raw: the repository ignores *.bin): a width-8 network with seeded random weights
(GroupNorm weights and biases randomised too), three inputs of 64, 64 and 96 pixels (the network is fully
convolutional), and PyTorch's float32 logits, plain and averaged over the 8 symmetries; it prints T1, ten times the
largest difference between those and the float64 logits.
--dump-inputs writes the input tensors cnn_ss.py predict makes from a stack (empty classes as zeros).
"""
import argparse
import os
import struct
import sys

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ss_trainer as T                      # noqa: E402
from cnn_data import normalise, resample    # noqa: E402
from cnn_ss import Net, d4                  # noqa: E402

MAGIC = 7741227103
VERSION = 1
GN_GROUPS = 8
NBLK = 5
GN_EPS = 1e-5


def tensor_bytes(t):
    a = t.detach().cpu().float().contiguous().numpy().astype('<f4')
    return struct.pack('<i', a.ndim) + struct.pack(f'<{a.ndim}i', *reversed(a.shape)) + a.tobytes()


def write_weights(model, width, cin, box, smpd, path):
    sd = model.state_dict()
    for m in model.modules():
        if isinstance(m, torch.nn.GroupNorm) and (m.num_groups != GN_GROUPS or abs(m.eps - GN_EPS) > 1e-12):
            raise SystemExit(f'a GroupNorm with {m.num_groups} groups and eps {m.eps}: the Fortran reader expects '
                             f'{GN_GROUPS} and {GN_EPS}')
    if len(model.body) != NBLK:
        raise SystemExit(f'{len(model.body)} blocks; the Fortran reader expects {NBLK}')
    keys = ['stem.0.weight', 'stem.1.weight', 'stem.1.bias']
    for b in range(NBLK):
        if model.body[b].skip is None:
            raise SystemExit(f'block {b + 1} has no skip convolution')
        keys += [f'body.{b}.{k}' for k in ('c1.weight', 'n1.weight', 'n1.bias', 'c2.weight', 'n2.weight', 'n2.bias',
                                           'skip.0.weight', 'skip.1.weight', 'skip.1.bias')]
    keys += ['head.weight', 'head.bias']
    if sorted(keys) != sorted(sd.keys()):
        raise SystemExit(f'unexpected tensors: {sorted(set(sd) ^ set(keys))}')
    with open(path, 'wb') as f:
        f.write(struct.pack('<qq', MAGIC, VERSION))
        f.write(struct.pack('<5i', width, cin, box, GN_GROUPS, NBLK))
        f.write(struct.pack('<f', smpd))
        for k in keys:
            f.write(tensor_bytes(sd[k]))


def write_inputs(xs, path):
    """xs: arrays (c, h, w); written as Fortran (w, h, c)."""
    with open(path, 'wb') as f:
        f.write(struct.pack('<i', len(xs)))
        for x in xs:
            c, h, w = x.shape
            f.write(struct.pack('<3i', w, h, c))
            f.write(np.ascontiguousarray(x, dtype='<f4').tobytes())


def logits(model, x, tta):
    with torch.no_grad():
        if tta:
            return torch.stack([model(v) for v in d4(x)]).mean(0)
        return model(x)


def cmd_model(a):
    ck = torch.load(a.model, map_location='cpu', weights_only=True)
    cfg = ck['config']
    cin = 1 if cfg['channels'] == 'avg' else 3
    model = Net(cin, cfg['width'], cfg['dropout'])
    model.load_state_dict(ck['state'])
    model.eval()
    write_weights(model, cfg['width'], cin, cfg['box'], cfg['smpd'], a.out)
    print(f'{a.out}: width {cfg["width"]}, {cin} channels, box {cfg["box"]}, {cfg["smpd"]} A, '
          f'{sum(p.numel() for p in model.parameters())} weights')


def cmd_fixture(a):
    os.makedirs(a.fixture, exist_ok=True)
    torch.manual_seed(1)
    model = Net(3, 8, 0.0)
    with torch.no_grad():
        for m in model.modules():
            if isinstance(m, torch.nn.GroupNorm):
                m.weight.copy_(1 + 0.2 * torch.randn_like(m.weight))
                m.bias.copy_(0.2 * torch.randn_like(m.bias))
    model.eval()
    rng = np.random.default_rng(2)
    xs = []
    for n in (64, 64, 96):
        g = np.arange(n) - n // 2
        y, x = np.meshgrid(g, g, indexing='ij')
        avg = sum(rng.uniform(0.5, 2) * np.exp(-((x - rng.uniform(-n / 4, n / 4)) ** 2 + (y - rng.uniform(-n / 4, n / 4)) ** 2)
                                               / (2 * rng.uniform(2, n / 8) ** 2)) for _ in range(5))
        noise = rng.normal(0, 0.3, (n, n))
        xs.append(normalise(avg + 0.2 * rng.normal(0, 0.3, (n, n)), avg + noise, avg - noise).astype(np.float32))
    write_weights(model, 8, 3, 64, 2.0, os.path.join(a.fixture, 'weights.raw'))
    write_inputs(xs, os.path.join(a.fixture, 'inputs.raw'))
    rows, diff = [], 0.0
    model64 = Net(3, 8, 0.0)
    model64.load_state_dict(model.state_dict())
    model64.double().eval()
    for i, x in enumerate(xs):
        t32 = torch.from_numpy(x)[None]
        p32, s32 = logits(model, t32, False)[0], logits(model, t32, True)[0]
        p64, s64 = logits(model64, t32.double(), False)[0], logits(model64, t32.double(), True)[0]
        diff = max(diff, float((p32.double() - p64).abs().max()), float((s32.double() - s64).abs().max()))
        rows.append(f'{i + 1} ' + ' '.join(f'{v:.9e}' for v in (*p32.tolist(), *s32.tolist())))
    with open(os.path.join(a.fixture, 'logits.txt'), 'w') as f:
        f.write(f'# cnn_ss.Net, width 8, 3 channels, seeds 1 (weights) and 2 (inputs), torch {torch.__version__}; '
                f'export_fortran.py --fixture\n')
        f.write('# input, plain logits 1 and 2, logits 1 and 2 averaged over the 8 symmetries (float32)\n')
        f.write(f'# float32 against float64: largest |difference| {diff:.3e}; T1 = {10 * diff:.3e}\n')
        f.write('\n'.join(rows) + '\n')
    print(f'fixture in {a.fixture}: float32 against float64, largest |difference| {diff:.3e}; T1 = {10 * diff:.3e}')


def cmd_dump(a):
    ck = torch.load(a.model, map_location='cpu', weights_only=True)
    cfg = ck['config']
    avg, smpd, _, _ = T.read_mrc(a.stk)
    stem = a.stk[:-4] if a.stk.endswith('.mrc') else a.stk
    ev, od = T.read_mrc(stem + '_even.mrc')[0], T.read_mrc(stem + '_odd.mrc')[0]
    xs = []
    for i in range(len(avg)):
        ch = [np.asarray(im[i], np.float64) for im in (avg, ev, od)]
        if ch[0].std() == 0:
            xs.append(np.zeros((3, cfg['box'], cfg['box']), np.float32))
            continue
        if abs(smpd - cfg['smpd']) > 1e-3 or ch[0].shape[-1] != cfg['box']:
            ch = [resample(im, smpd, cfg['smpd'], cfg['box']) for im in ch]
        xs.append(normalise(*ch).astype(np.float32))
    write_inputs(xs, a.out)
    print(f'{a.out}: {len(xs)} inputs of {cfg["box"]} px from {a.stk} at {smpd:.4f} A')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--model')
    ap.add_argument('--out')
    ap.add_argument('--fixture')
    ap.add_argument('--dump-inputs', action='store_true')
    ap.add_argument('--stk')
    a = ap.parse_args()
    if a.fixture:
        cmd_fixture(a)
    elif a.dump_inputs:
        cmd_dump(a)
    elif a.model and a.out:
        cmd_model(a)
    else:
        ap.error('give --model and --out, --fixture, or --dump-inputs with --model, --stk and --out')


if __name__ == '__main__':
    main()
