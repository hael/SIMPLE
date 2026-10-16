"""Analytic picking adapter and coordinate quality checks."""
import json
import math
from resource_planner.workflows.common import SpatialAdapter, GENERATOR

def picking_metrics(predicted, expected, tolerance):
    if not expected: raise ValueError('empty picking truth')
    pairs=sorted((math.hypot(x-a,y-b),i,j) for i,(x,y) in enumerate(predicted)
                 for j,(a,b) in enumerate(expected) if math.hypot(x-a,y-b)<=tolerance)
    used_p,used_t=set(),set()
    for distance,i,j in pairs:
        if i not in used_p and j not in used_t:
            used_p.add(i);used_t.add(j)
    precision=len(used_p)/len(predicted) if predicted else 0.0
    recall=len(used_t)/len(expected)
    return {'passed':True,'quality_qualified':precision>=0.8 and recall>=0.8,
            'precision':precision,'recall':recall,'picked':len(predicted),'matched':len(used_p)}


def validate_pick(work, truth):
    paths=list(work.glob('*.box'))
    if len(paths)!=1:raise ValueError('expected one picking box file')
    predicted=[]
    for line in paths[0].read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):continue
        row=list(map(float,line.split()[:4]))
        if len(row)!=4 or not all(math.isfinite(v) for v in row) or min(row[2:])<=0:
            raise ValueError('invalid picking coordinates')
        predicted.append((row[0]+row[2]/2,row[1]+row[3]/2))
    expected=[(x+w/2,y+h/2) for x,y,w,h in truth['coordinates']]
    return picking_metrics(predicted,expected,truth['box']*0.15)


class PickAdapter(SpatialAdapter):
    name='pick'
    qualification='provisional synthetic Gaussian picking; single-micrograph shared-memory mode; precision/recall >=0.8'

    def cases(self,args):return args.sizes

    def science(self,args):
        return {'picker':'new','moldiam':args.box*0.6*1.3,'pcontrast':'black',
                'pick_roi':'no','backgr_subtr':'no','lp':20,'ndev':2.0,'thres':0}

    def metadata(self,args):
        return {'sizes':args.sizes,'threads':args.threads,'workload':{'smpd':1.3,
                'fixture_box':args.box,'particle_spacing_pixels':2*args.box,
                'noise_sigma':0.03,'micrographs_per_run':1,**self.science(args)}}

    def prepare(self,ctx,size):
        folder=ctx.root/str(size);folder.mkdir()
        inputs=folder/'fixture'
        ctx.command([ctx.args.fixture_python,str(GENERATOR),'micrograph',str(inputs),
                     '--size',str(size),'--box',str(ctx.args.box)],folder,'fixture',timeout=300)
        return {'folder':folder,'size':size,'input':inputs,
                'files':[inputs/name for name in ('micrograph.mrc','micrograph.box','truth.json')],
                'truth':json.loads((inputs/'truth.json').read_text())}

    def jobs(self,ctx,fixture,directory,nthr,workers):
        project=ctx.project(directory,'particles',fixture['input']/'micrograph.mrc',particles=False,ctf_flag='no')
        params={**self.science(ctx.args),'fromp':1,'top':1,'nparts':1,'part':1,'nthr':nthr,
                'outfile':f'{self.name}.simple','memreport':'yes','memreport_interval':1}
        argv=[str(ctx.args.simple_exec),f'prg={self.name}',f'projfile={project}',
              'mkdir=no','qsys_name=local',*[f'{k}={v}' for k,v in params.items()]]
        return [(argv,directory)]

    def validate(self,ctx,fixture,work):return validate_pick(work,fixture['truth'])


