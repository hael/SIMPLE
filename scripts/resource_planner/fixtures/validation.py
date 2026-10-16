#!/usr/bin/env python3
"""Molecular workflow input transforms and external checks (NumPy, Python >=3.6)."""
import argparse
import json
from pathlib import Path
import re
import struct

import numpy as np
if __package__:
    from .particles import header
else:
    from particles import header
if __package__:
    from .spatial import read_mrc, correlation
else:
    from spatial import read_mrc, correlation


def sampling(path):
    with path.open('rb') as stream:
        h = stream.read(1024)
    return struct.unpack_from('<f', h, 40)[0] / struct.unpack_from('<i', h, 28)[0]


def write(path, data, smpd):
    with path.open('xb') as stream:
        stream.write(header(data.shape[-1], data.shape[0], smpd))
        np.asarray(data, dtype='<f4').tofile(stream)


def subset(source, output, count, seed):
    data = read_mrc(source)
    rng = np.random.RandomState(seed)
    indices = np.concatenate([rng.permutation(len(data)) for _ in range((count+len(data)-1)//len(data))])[:count]
    write(output, data[indices], sampling(source))
    output.with_suffix('.indices.json').write_text(json.dumps(indices.tolist())+'\n')


def rows(path, key):
    result = []
    for line in path.read_text().splitlines():
        record = dict(re.findall(r'(\w+)=\s*([^\s]+)', line))
        if key in record:
            result.append(record)
    return result


def boxes(path):
    result = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        values = list(map(float, line.split()[:4]))
        if len(values) != 4 or not np.isfinite(values).all() or min(values[2:]) <= 0:
            raise ValueError('invalid particle boxes')
        result.append(values)
    return result


def unique(work, pattern):
    paths = list(work.rglob(pattern))
    if len(paths) != 1:
        raise ValueError('expected one %s in %s, found %d' % (pattern, work, len(paths)))
    return paths[0]


def matched_fraction(predicted, expected, tolerance):
    candidates = sorted((float(np.hypot(x-a,y-b)),i,j)
        for i,(x,y) in enumerate(predicted) for j,(a,b) in enumerate(expected)
        if np.hypot(x-a,y-b) <= tolerance)
    used_p, used_t = set(), set()
    for _,i,j in candidates:
        if i not in used_p and j not in used_t:
            used_p.add(i); used_t.add(j)
    return len(used_p)/max(1,len(predicted)), len(used_t)/max(1,len(expected))


def validate(workflow, work, truth, preparation):
    result = {'passed': True, 'quality_qualified': False}
    if workflow == 'motion_correct':
        correlations = []
        for movie in truth['movies']:
            actual = read_mrc(unique(work, Path(movie['movie']).stem+'_intg.mrc'))
            expected = read_mrc(Path(movie['optimal']))
            if actual.shape != expected.shape:
                raise ValueError('motion output dimensions differ from simulation truth')
            correlations.append(correlation(actual, expected))
        result.update(minimum_image_correlation=min(correlations), movies=len(correlations),
                      quality_qualified=min(correlations)>=0.98)
    elif workflow == 'ctf_estimate':
        errors = []
        for movie in truth['movies']:
            expected = rows(Path(movie['parameters']), 'dfx')[0]
            fitted = rows(unique(work, Path(movie['movie']).stem+'*_ctf.txt'), 'dfx')[0]
            a = sorted(float(expected[k]) for k in ('dfx','dfy'))
            b = sorted(float(fitted[k]) for k in ('dfx','dfy'))
            if not np.isfinite(b).all():
                raise ValueError('nonfinite CTF estimate')
            errors.append(max(abs(x-y) for x,y in zip(a,b)))
        result.update(max_defocus_error_um=max(errors), quality_qualified=max(errors)<=0.1)
    elif workflow == 'pick':
        scores = []
        for movie in truth['movies']:
            expected = boxes(Path(movie['boxes']))
            predicted = boxes(unique(work, Path(movie['movie']).stem+'*.box'))
            center = lambda rows: [(x+w/2,y+h/2) for x,y,w,h in rows]
            scores.append(matched_fraction(center(predicted), center(expected), truth['box']*0.15))
        result.update(minimum_precision=min(s[0] for s in scores), minimum_recall=min(s[1] for s in scores),
                      quality_qualified=all(min(s)>=0.8 for s in scores))
    elif workflow == 'solve2D':
        candidates = [p for p in work.glob('cavgs_iter*.mrc') if re.fullmatch(r'cavgs_iter\d+\.mrc',p.name)]
        if not candidates:
            raise ValueError('missing class averages')
        data = read_mrc(max(candidates,key=lambda p:int(re.search(r'\d+',p.stem).group())))
        assignments = rows(work/'assignments.log','class')
        expected_count = read_mrc(Path(truth['particle_stack'])).shape[0]
        if not assignments or len(assignments)!=expected_count or data.shape != (8,truth['box'],truth['box']) or not np.any(data):
            raise ValueError('invalid class averages/assignments')
        labels = [float(r['class']) for r in assignments]
        if not all(np.isfinite(v) and v.is_integer() and 1<=v<=8 for v in labels):
            raise ValueError('invalid class labels')
        occupied = len(set(labels))
        # Continuous molecular projections have no eight-class ground truth.
        # Preserve the suite's structural check without presenting it as accuracy.
        result.update(particles_assigned=len(labels), occupied_classes=occupied,
                      quality_qualified=False, structural_passed=occupied>=2,
                      qualification_reason='Molecular class accuracy gate not yet established; structural checks only')
    else:
        expected = read_mrc(Path(truth['volume']))
        options = []
        for tag in ('direct','mirror'):
            cc = float(rows(work/(tag+'_dock.txt'),'cc')[0]['cc'])
            if not np.isfinite(cc):
                raise ValueError('nonfinite docking correlation')
            options.append((cc,tag))
        _,tag = max(options)
        path = work/(tag+'_docked.mrc')
        actual = read_mrc(path)
        if actual.shape != expected.shape or abs(sampling(path)-truth['smpd'])>0.01:
            raise ValueError('reconstruction grid differs from molecular truth')
        cc = correlation(actual, expected)
        # Masked shell FSC against independent molecular truth, after docking both hands.
        n = actual.shape[0]
        z,y,x = np.mgrid[:n,:n,:n]-n//2
        radius = np.sqrt(x*x+y*y+z*z)
        edge = 90/truth['smpd']
        mask = np.where(radius<=edge,1.0,np.where(radius<edge+3,0.5*(1+np.cos(np.pi*(radius-edge)/3)),0))
        a = np.fft.fftn(actual*mask); b = np.fft.fftn(expected*mask)
        f = np.fft.fftfreq(n)*n
        z,y,x = np.meshgrid(f,f,f,indexing='ij')
        shell = np.rint(np.sqrt(x*x+y*y+z*z)).astype(int)
        resolution = None
        for i in range(1,n//2+1):
            keep = shell==i
            denom = np.sqrt(np.sum(abs(a[keep])**2)*np.sum(abs(b[keep])**2))
            fsc = float(np.real(np.sum(a[keep]*np.conj(b[keep])))/denom) if denom>0 else 0
            if fsc<0.143:
                resolution = n*truth['smpd']/i
                break
        if resolution is None:
            resolution = 2*truth['smpd']
        assignments = rows(work/'assignments.log','state')
        expected_count = read_mrc(Path(truth['particle_stack'])).shape[0]
        valid = bool(assignments) and len(assignments)==expected_count and all(float(r['state'])==1 and
            all(np.isfinite(float(r.get(k,'nan'))) for k in ('e1','e2','e3')) for r in assignments)
        result.update(volume_correlation=cc, fsc0143_angstrom=resolution, selected_hand=tag,
                      particles_assigned=len(assignments), quality_qualified=valid and cc>=0.8 and resolution<=40)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['subset','mirror','validate'])
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--count', type=int, default=100)
    p.add_argument('--seed', type=int, default=20260923)
    p.add_argument('--workflow')
    p.add_argument('--truth', type=Path)
    p.add_argument('--preparation', type=Path)
    args = p.parse_args()
    if args.action=='subset':
        subset(args.input,args.output,args.count,args.seed)
    elif args.action=='mirror':
        # SIMPLE's real-space x mirror is a reversal along the first image axis.
        write(args.output,read_mrc(args.input)[:,:,::-1],sampling(args.input))
    else:
        result = validate(args.workflow,args.input,json.loads(args.truth.read_text()),args.preparation)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
