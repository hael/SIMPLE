#!/usr/bin/env python3
"""Seeded eight-shape particle fixture (NumPy, Python >=3.6)."""
import argparse
import json
from pathlib import Path
import struct
import numpy as np


def header(box, count, smpd):
    data=bytearray(1024)
    struct.pack_into('<4i',data,0,box,box,count,2)
    struct.pack_into('<3i',data,28,box,box,count)
    struct.pack_into('<3f',data,40,box*smpd,box*smpd,count*smpd)
    struct.pack_into('<3f',data,52,90,90,90)
    struct.pack_into('<3i',data,64,1,2,3)
    data[208:216]=b'MAP DA\x00\x00'
    return data


def shape(x,y,label):
    # Fixed asymmetric low-frequency shapes exercise alignment without ambiguous
    # rotationally symmetric rings or nearly identical frequency bands.
    rng = np.random.RandomState(8200 + label)
    data = np.zeros_like(x)
    for index in range(6):
        cx, cy = rng.uniform(-0.4, 0.4, 2)
        width = rng.uniform(0.06, 0.12)
        data += rng.uniform(0.5, 1.5) * np.exp(-((x-cx)**2+(y-cy)**2)/(2*width**2))
    return data



def generate(root, count, box):
    root.mkdir(parents=True,exist_ok=False)
    axis=(np.arange(box,dtype=np.float32)-box/2)/(box/2)
    y,x=np.meshgrid(axis,axis,indexing='ij')
    rng=np.random.RandomState(20261015)
    labels=np.arange(count)%8
    rng.shuffle(labels)
    with (root/'particles.mrcs').open('xb') as stream:
        stream.write(header(box,count,1.3))
        for label in labels:
            angle=0.0
            xr=x*np.cos(angle)+y*np.sin(angle)
            yr=-x*np.sin(angle)+y*np.cos(angle)
            particle=shape(xr,yr,int(label))+0.03*rng.standard_normal((box,box))
            particle=(particle-particle.mean())/particle.std()
            particle.astype('<f4').tofile(stream)
    with (root/'truth_classes.mrcs').open('xb') as stream:
        stream.write(header(box,8,1.3))
        for label in range(8):
            data=shape(x,y,label)
            ((data-data.mean())/data.std()).astype('<f4').tofile(stream)
    # Controlled imperfect seed: 20% of labels deliberately belong to the next class.
    # This exercises the supported seeded restart without the expensive random-start phase.
    seed_labels=labels.copy()
    seed_labels[::5]=(seed_labels[::5]+1)%8
    rows=['class=%d state=1 e3=0 eo=%d'%(int(label)+1,i%2) for i,label in enumerate(seed_labels)]
    (root/'seed.txt').write_text('\n'.join(rows)+'\n')
    (root/'truth.json').write_text(json.dumps({'nptcls':count,'box':box,'ncls':8,
        'smpd':1.3,'noise_sigma':0.03,'fixture_version':5,'seed':20261015,'numpy_version':np.__version__,
        'labels':(labels+1).tolist(),'seed_labels':(seed_labels+1).tolist(),'ctf':'no','shifts':'none','rotations':'none'},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('--particles',type=int,default=2000)
    p.add_argument('--box',type=int,default=128)
    args=p.parse_args()
    if args.particles<16 or args.box<32 or args.box%2:
        p.error('requires at least 16 particles and an even box >=32')
    generate(args.output,args.particles,args.box)
