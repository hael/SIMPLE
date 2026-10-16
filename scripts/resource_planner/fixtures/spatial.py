#!/usr/bin/env python3
"""NumPy-only fixture generation and independent spatial validators (Python >=3.6)."""
import argparse
import json
from pathlib import Path
import struct
import numpy as np
if __package__:
    from .particles import header
else:
    from particles import header


def write_mrc(path, data):
    with path.open('xb') as stream:
        stream.write(header(data.shape[-1], data.shape[0], 1.3))
        data.astype('<f4').tofile(stream)


def read_mrc(path):
    with path.open('rb') as stream:
        h=stream.read(1024)
    nx,ny,nz,mode=struct.unpack_from('<4i',h)
    offset=1024+struct.unpack_from('<i',h,92)[0]
    if mode!=2 or min(nx,ny,nz)<1 or offset<1024:
        raise ValueError('expected positive float32 MRC dimensions')
    data=np.memmap(str(path),dtype='<f4',mode='r',offset=offset,shape=(nz,ny,nx))
    if not np.isfinite(data).all(): raise ValueError('nonfinite MRC values')
    return data


def micrograph(root, size, box):
    root.mkdir(parents=True,exist_ok=False)
    rng=np.random.RandomState(20261015)
    data=rng.normal(0,0.03,(size,size)).astype(np.float32)
    yy,xx=np.mgrid[:box,:box]-box//2
    blob=np.exp(-(xx**2+yy**2)/(2*(box*0.10)**2)).astype(np.float32)
    coordinates=[]
    for y in range(box,size-box,2*box):
        for x in range(box,size-box,2*box):
            x0,y0=x-box//2,y-box//2
            data[y0:y0+box,x0:x0+box]-=blob
            coordinates.append([x0,y0,box,box])
    if not coordinates: raise ValueError('micrograph too small for particle box')
    write_mrc(root/'micrograph.mrc',data[None,:,:])
    (root/'micrograph.box').write_text(''.join(' '.join(map(str,row))+'\n' for row in coordinates))
    (root/'truth.json').write_text(json.dumps({'size':size,'box':box,'smpd':1.3,
        'pcontrast':'black','coordinates':coordinates,'particles':len(coordinates),
        'seed':20261015,'numpy_version':np.__version__},indent=2)+'\n')


def volume(root, box):
    root.mkdir(parents=True,exist_ok=False)
    axis=np.linspace(-1,1,box,dtype=np.float32)
    z,y,x=np.meshgrid(axis,axis,axis,indexing='ij')
    data=(np.exp(-((x+.18)**2/.15+(y-.05)**2/.07+(z+.10)**2/.10))
          +.75*np.exp(-((x-.34)**2/.045+(y+.24)**2/.11+(z-.20)**2/.055))
          +.48*np.exp(-((x+.38)**2/.035+(y+.32)**2/.040+(z-.30)**2/.075)))
    data-=data.min();data/=data.max()
    write_mrc(root/'volume.mrc',data)
    (root/'truth.json').write_text(json.dumps({'box':box,'smpd':1.3,'seed':20261015,
        'model':'asymmetric Gaussian volume','numpy_version':np.__version__},indent=2)+'\n')


def correlation(a,b):
    a=np.asarray(a,dtype=np.float64).ravel(); b=np.asarray(b,dtype=np.float64).ravel()
    a-=a.mean();b-=b.mean()
    denom=np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.dot(a,b)/denom) if denom>0 else 0.0


def reconstruction(fixture, output):
    expected=read_mrc(fixture/'volume.mrc'); actual=read_mrc(output)
    if len(set(actual.shape))!=1 or actual.shape[0]<8 or actual.shape[0]>expected.shape[0]:
        raise ValueError('invalid reconstructed volume shape')
    # Compare smooth analytic truth at the output grid if the stage downscaled.
    for axis in range(3):
        if expected.shape[axis]!=actual.shape[axis]:
            old=np.linspace(0,1,expected.shape[axis]);new=np.linspace(0,1,actual.shape[axis])
            expected=np.apply_along_axis(lambda v:np.interp(new,old,v),axis,expected)
    cc=correlation(actual,expected)
    return {'passed':True,'quality_qualified':cc>=0.8,'volume_correlation':cc,
            'output_box':actual.shape[0]}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['micrograph','volume','validate_volume'])
    p.add_argument('root',type=Path)
    p.add_argument('--size',type=int,default=4096)
    p.add_argument('--box',type=int,default=128)
    p.add_argument('--input',type=Path)
    p.add_argument('--result',type=Path)
    args=p.parse_args()
    if args.box<32 or args.box%2: p.error('box must be even and >=32')
    if args.action=='micrograph':micrograph(args.root,args.size,args.box)
    elif args.action=='volume':volume(args.root,args.box)
    else:
        if args.input is None or args.result is None:p.error('validation requires --input and --result')
        result=reconstruction(args.root,args.input)
        args.result.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
