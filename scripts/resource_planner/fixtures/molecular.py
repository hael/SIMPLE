#!/usr/bin/env python3
"""External deterministic molecular fixtures; NumPy only, Python >=3.6.

Uses bundled coordinate data, not SIMPLE source files or native simulation code.
Equal-weight atom density with Gaussian smoothing is a calibration model, not
the native atomic scattering-factor implementation.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
if __package__:
    from .validation import write, sampling
else:
    from validation import write, sampling
if __package__:
    from .spatial import read_mrc
else:
    from spatial import read_mrc
if __package__:
    from .particles import header
else:
    from particles import header

DATA = Path(__file__).with_name('data')


def smooth(data, smpd):
    freq = [np.fft.fftfreq(n,d=smpd) for n in data.shape]
    q2 = sum(f.reshape([len(f) if j==i else 1 for j in range(data.ndim)])**2
             for i,f in enumerate(freq))
    return np.fft.ifftn(np.fft.fftn(data)*np.exp(-2*np.pi**2*2.0**2*q2)).real.astype(np.float32)


def density(xyz, box, smpd):
    pixels = np.rint(xyz/smpd+box//2).astype(int)
    if np.any(pixels<0) or np.any(pixels>=box):
        raise ValueError('coordinate model exceeds fixture box')
    shape = (box,)*xyz.shape[1]
    data = np.zeros(shape,dtype=np.float32)
    np.add.at(data,tuple(pixels[:,i] for i in reversed(range(xyz.shape[1]))),1)
    return smooth(data,smpd)


def molecule(root, model, box, smpd):
    metadata = json.loads((DATA/'molecules.json').read_text())['models'][model]
    source = DATA/metadata['file']
    if hashlib.sha256(source.read_bytes()).hexdigest()!=metadata['sha256']:
        raise ValueError('bundled molecular coordinates do not match their manifest')
    with gzip.open(str(source),'rt') as stream:
        xyz = np.loadtxt(stream)
    if xyz.shape != (metadata['atoms'],3) or not np.isfinite(xyz).all():
        raise ValueError('invalid molecular coordinates')
    xyz -= (xyz.min(axis=0)+xyz.max(axis=0))/2
    # A consistent frame/center is shared by the volume and all projections.
    write(root/'truth.mrc',density(xyz,box,smpd),smpd)
    projections = []
    for i in range(100):
        z = 1-2*(i+0.5)/100
        phi = i*np.pi*(3-np.sqrt(5))
        direction = np.array([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
        u = np.cross(direction,[0,0,1]); u /= np.linalg.norm(u)
        v = np.cross(direction,u)
        projections.append(density(xyz.dot(np.column_stack([u,v])),box,smpd))
    write(root/'projections.mrcs',np.asarray(projections),smpd)
    provenance = {'generator':'Python equal-weight atom density, Gaussian sigma 2 A',
                  'model':model,'coordinate_sha256':metadata['sha256'],
                  'numpy_version':np.__version__,'box':box,'smpd':smpd,'projections':100}
    (root/'generator.json').write_text(json.dumps(provenance,indent=2)+'\n')


def movie(root, size, seed):
    source = root/'projections.mrcs'
    particles = read_mrc(source)
    box, smpd = particles.shape[-1], sampling(source)
    rng = np.random.RandomState(seed)
    positions = [(x,y) for y in range(box,size-box+1,2*box) for x in range(box,size-box+1,2*box)]
    if len(particles)>len(positions):
        raise ValueError('too many particles for nonoverlapping fixture placement')
    rng.shuffle(positions)
    base = np.zeros((size,size),dtype=np.float32)
    coords = []
    for particle,(x,y) in zip(particles,positions):
        x0,y0 = x-box//2,y-box//2
        base[y0:y0+box,x0:x0+box] -= particle
        coords.append((x0,y0,box,box))
    base /= max(float(base.std()),1e-12)
    # Stationary specimen background gives the estimator a broad CTF spectrum.
    base += rng.normal(0,1,base.shape).astype(np.float32)
    spectrum = np.fft.rfft2(base)
    fy = np.fft.fftfreq(size,d=smpd)[:,None]
    fx = np.fft.rfftfreq(size,d=smpd)[None,:]
    wavelength = 12.2639/np.sqrt(300000*(1+0.97845e-6*300000))
    for first in range(0,size,128):
        q2 = fx*fx+fy[first:first+128]**2
        chi = np.pi*wavelength*q2*(20000-0.5*wavelength**2*q2*2.7e7)
        spectrum[first:first+128] *= np.sin(chi+np.arcsin(0.1))
    base = np.fft.irfft2(spectrum,s=base.shape).astype(np.float32)
    del spectrum
    noise_sigma = float(base.std())/np.sqrt(0.2)
    shifts = [0,0,0,0,0,1,2,3]
    optimal = np.zeros_like(base)
    with (root/'simulate_movie.mrc').open('xb') as stream:
        stream.write(header(size,len(shifts),smpd))
        for shift in shifts:
            frame = np.roll(base,shift,axis=1)+rng.normal(0,noise_sigma,base.shape).astype(np.float32)
            np.asarray(frame,dtype='<f4').tofile(stream)
            optimal += np.roll(frame,-shift,axis=1)/len(shifts)
    write(root/'optimal_movie_average.mrc',optimal[None,:,:],smpd)
    (root/'simulate_movie.box').write_text(''.join('%d %d %d %d\n'%row for row in coords))
    (root/'simulate_movie_params.txt').write_text('dfx=2 dfy=2 angast=0 '+
        ' '.join('x%d=%d y%d=0'%(i+1,s,i+1) for i,s in enumerate(shifts))+'\n')
    (root/'generator.json').write_text(json.dumps({'seed':seed,'numpy_version':np.__version__,
        'snr':0.2,'stationary_background_sigma':1,'frames':8,'ctf':'isotropic 2 um, 300 kV, Cs 2.7 mm, amplitude contrast 0.1'},indent=2)+'\n')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['molecule','movie'])
    p.add_argument('root',type=Path)
    p.add_argument('--model',choices=['6vxx','1jyx'],default='6vxx')
    p.add_argument('--box',type=int,default=128)
    p.add_argument('--smpd',type=float,default=1.875)
    p.add_argument('--size',type=int,default=4096)
    p.add_argument('--seed',type=int,default=20260923)
    a = p.parse_args()
    if a.action=='molecule':
        molecule(a.root,a.model,a.box,a.smpd)
    else:
        movie(a.root,a.size,a.seed)
