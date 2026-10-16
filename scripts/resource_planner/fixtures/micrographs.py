#!/usr/bin/env python3
"""Generate stationary CTF-filtered Gaussian noise; requires NumPy (Python >=3.6)."""
import argparse
import json
from pathlib import Path
import struct
import numpy as np


def generate(path, size):
    # Frequencies in inverse Angstrom; 300 kV relativistic electron wavelength.
    rng = np.random.RandomState(20261015)
    spectrum = np.fft.rfft2(rng.standard_normal((size, size)).astype(np.float32))
    fy = np.fft.fftfreq(size, d=1.3)[:, None]
    fx = np.fft.rfftfreq(size, d=1.3)[None, :]
    wavelength = 12.2639 / np.sqrt(300000.0 * (1 + 0.97845e-6 * 300000.0))
    # Isotropic 2-micron underfocus avoids an undefined astigmatism-angle gate.
    for first in range(0, size, 128):
        q2 = fx*fx + fy[first:first+128]*fy[first:first+128]
        chi = np.pi*wavelength*q2*(20000.0 - 0.5*wavelength**2*q2*2.7e7)
        spectrum[first:first+128] *= np.sin(chi + np.arcsin(0.1))
    data = np.fft.irfft2(spectrum, s=(size, size)).astype('<f4')
    header = bytearray(1024)
    struct.pack_into('<4i', header, 0, size, size, 1, 2)
    struct.pack_into('<3i', header, 28, size, size, 1)
    struct.pack_into('<3f', header, 40, size*1.3, size*1.3, 1.3)
    struct.pack_into('<3f', header, 52, 90, 90, 90)
    struct.pack_into('<3i', header, 64, 1, 2, 3)
    struct.pack_into('<3f', header, 76, float(data.min()), float(data.max()), float(data.mean()))
    header[208:212] = b'MAP '
    header[212:216] = b'DA\x00\x00'
    with path.open('xb') as stream:
        stream.write(header)
        data.tofile(stream)
    path.with_suffix('.truth.json').write_text(json.dumps({
        'size':size, 'smpd':1.3, 'kv':300, 'cs':2.7, 'fraca':0.1,
        'dfx':2.0, 'dfy':2.0, 'phshift':0, 'seed':20261015,
        'numpy_version':np.__version__, 'model':'CTF-filtered stationary Gaussian noise'}, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--size', type=int, choices=[4096,8192], required=True)
    args = parser.parse_args()
    generate(args.output, args.size)
