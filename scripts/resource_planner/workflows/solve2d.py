"""Particle classification scientific settings, quality gates and profile selection."""
from array import array
from collections import Counter
import json
import math
from pathlib import Path
import re
import statistics
import struct
import sys
from resource_planner.campaign import file_hash

def science(count, box, extr_lim=20, nits_per_stage=5):
    return {'ncls':8,'mskdiam':box*1.3*0.8,'smpd':1.3,'nstages':6,
            'nits_per_stage':nits_per_stage,'extr_lim':extr_lim,'eo_stage':'no','autoscale':'no',
            'refine':'prob','nsample':count,'ctf':'no','rank_cavgs':'no',
            'cache':'no','cls_init':'prev','center':'no','trs':0,'lpstart':6,'lpstop':3}


def assignment_metrics(labels, assigned):
    if len(labels)!=len(assigned) or not labels:
        raise ValueError('particle assignment count does not match external truth')
    joint=Counter(zip(labels,assigned)); truth=Counter(labels); fitted=Counter(assigned)
    choose=lambda n:n*(n-1)/2
    a=sum(choose(n) for n in joint.values())
    b=sum(choose(n) for n in truth.values());c=sum(choose(n) for n in fitted.values())
    expected=b*c/choose(len(labels));denominator=(b+c)/2-expected
    ari=(a-expected)/denominator if denominator else 1.0
    purity=sum(max(joint[t,f] for t in truth) for f in fitted)/len(labels)
    return {'adjusted_rand_index':ari,'purity':purity,'occupied_classes':len(fitted)}


def validate(directory, truth, assignments):
    stacks=sorted(p for p in directory.glob('cavgs_iter*.mrc')
                  if re.fullmatch(r'cavgs_iter\d+\.mrc',p.name))
    if not stacks:
        raise ValueError('missing class-average stack')
    path=max(stacks,key=lambda p:int(re.search(r'\d+',p.stem).group()))
    with path.open('rb') as stream:
        header=stream.read(1024)
        if struct.unpack_from('<4i',header)!=(truth['box'],truth['box'],truth['ncls'],2):
            raise ValueError('incorrect class-average dimensions or MRC mode')
        stream.seek(1024+struct.unpack_from('<i',header,92)[0])
        data=array('f');data.fromfile(stream,truth['box']**2*truth['ncls'])
    if sys.byteorder!='little': data.byteswap()
    if not all(math.isfinite(v) for v in data) or not any(abs(v)>1e-6 for v in data):
        raise ValueError('class averages contain nonfinite or all-zero data')
    assigned=[]
    for line in assignments.read_text().splitlines():
        values=dict(re.findall(r'(\w+)=\s*([^\s]+)',line))
        if 'class' in values:
            cls=float(values['class'])
            if not math.isfinite(cls) or not cls.is_integer() or not 1<=cls<=truth['ncls']:
                raise ValueError('invalid particle class assignment')
            assigned.append(int(cls))
    metrics=assignment_metrics(truth['labels'],assigned)
    if 'seed_labels' in truth:
        metrics['initial_adjusted_rand_index'] = assignment_metrics(truth['labels'], truth['seed_labels'])['adjusted_rand_index']
    # Require external-label recovery before exporting resource recommendations.
    return {'passed':True,'class_averages':str(path),'particles_assigned':len(assigned),
            'quality_qualified':metrics['adjusted_rand_index']>=0.8 and metrics['occupied_classes']==truth['ncls'],**metrics}


def export(root,plan,records):
    summary=[]
    for nthr in plan['threads']:
        measured=[r for r in records if r['nthr']==nthr and not r['warmup']]
        if len(measured)!=plan['repetitions'] or not all(r['validation']['passed'] for r in measured):
            raise ValueError('incomplete or invalid observations')
        summary.append({'nthr':nthr,'nparts':1,'ncunits':1,
            'median_wall_seconds':statistics.median(r['wall_seconds'] for r in measured),
            'timing_cv':(statistics.stdev(r['wall_seconds'] for r in measured) / statistics.mean(r['wall_seconds'] for r in measured)) if len(measured)>1 else None,
            'max_rss_mib':max(r['peak_rss_mib'] for r in measured),
            'minimum_adjusted_rand_index':min(r['validation']['adjusted_rand_index'] for r in measured),
            'median_purity':statistics.median(r['validation']['purity'] for r in measured)})
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not all(r['validation']['quality_qualified'] for r in records): return None
    best=min(r['median_wall_seconds'] for r in summary)
    chosen=min((r for r in summary if r['median_wall_seconds']<=best*1.05),key=lambda r:r['nthr'])
    if file_hash(Path(plan['simple_exec']))!=plan['binary_sha256']:
        raise ValueError('executable changed during benchmark')
    profile={'schema_version':1,'profile_type':'simple_resource_profile','commander':'solve2D',
        'auto_apply':False,'qualification':'provisional synthetic seeded-restart solve2D (20% mislabelled input); external-label ARI >=0.8; de novo starts and real particles unqualified',
        'hardware':plan['hardware'],'cpu_ids':plan['cpus'],
        'binary_sha256':{plan['simple_exec']:plan['binary_sha256']},
        'workload':{'nptcls':plan['particles'],'box':plan['box'],
                    'initialization':plan.get('initialization', 'unspecified'),
                    'fixture':plan.get('fixture', {}),**plan['science']},
        'settings':{'nthr':chosen['nthr'],'nparts':1,'ncunits':1,
                    'job_memory_per_task_mib':max(1024,math.ceil(chosen['max_rss_mib']*1.5/256)*256),
                    'median_wall_seconds':chosen['median_wall_seconds'], 'timing_cv':chosen['timing_cv'],
                    'measured_repetitions':plan['repetitions']},
        'resource_policy':{'execution_mode':'shared','partition_policy_measured':False,
                           'reserve_memory_min_mib':4096,'reserve_memory_fraction':0.1,
                           'explicit_job_parameters_take_precedence':True},
        'evidence':str(root),'seed':20261015}
    path=root/'desktop_profile.json'
    path.write_text(json.dumps(profile,indent=2)+'\n');return path



from resource_planner.adapters import build_command
class Solve2DAdapter:
    name = 'solve2D'
    requires_numpy = True
    shared = True
    rate = None

    def cases(self, args):
        return [args.box]

    def candidates(self, args):
        return [(nthr, 1) for nthr in args.threads]

    def metadata(self, args):
        return {'particles': args.particles, 'box': args.box, 'threads': args.threads,
                'simple_exec': str(args.simple_exec),
                'initialization': 'seeded_restart_20_percent_mislabelled',
                'fixture': {'version': 5, 'noise_sigma': 0.03, 'rotations': 'none',
                            'shifts': 'none', 'shapes': 'asymmetric Gaussian mixtures'},
                'science': science(args.particles, args.box, args.extr_lim, args.nits_per_stage)}

    def prepare(self, ctx, size):
        generator = Path(__file__).resolve().parents[1] / 'fixtures' / 'particles.py'
        ctx.command([ctx.args.fixture_python, str(generator), str(ctx.root / 'fixture'),
                     '--particles', str(ctx.args.particles), '--box', str(size)],
                    ctx.root, 'fixture', timeout=300)
        return {'folder': ctx.root, 'size': size,
                'files': [ctx.root / 'fixture' / name for name in ('particles.mrcs', 'seed.txt', 'truth.json', 'truth_classes.mrcs')],
                'truth': json.loads((ctx.root / 'fixture/truth.json').read_text())}

    def jobs(self, ctx, fixture, directory, nthr, workers):
        ctx.project(directory, 'classes', ctx.root / 'fixture/particles.mrcs', particles=True, ctf_flag='no',
                    oritab=ctx.root / 'fixture/seed.txt')
        case = {'commander': self.name, 'mode': 'shared', 'project': 'classes/classes.simple',
                'work_objects': ctx.args.particles, 'science': self.metadata(ctx.args)['science']}
        return [(build_command(case, {'nthr': nthr}, ctx.args.simple_exec,
                               directory, ctx.root / 'fixture'), directory)]

    def validate(self, ctx, fixture, work):
        ctx.command([str(ctx.args.simple_exec), 'prg=print_project_field',
                     f'projfile={work}/classes/classes.simple', 'oritype=ptcl2D', 'mkdir=no'],
                    work, 'assignments')
        return validate(work, fixture['truth'], work / 'assignments.log')

    def export(self, ctx, plan, records, summary):
        return export(ctx.root, {**plan, 'binary_sha256': plan['binaries'][str(ctx.args.simple_exec)]}, records)

