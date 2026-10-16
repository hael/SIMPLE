"""Shared adapter mechanics and scientific qualification."""
import json
import math
import statistics
from pathlib import Path
from resource_planner.campaign import file_hash
from resource_planner.profiles import write_profile
GENERATOR = Path(__file__).resolve().parents[1] / "fixtures/spatial.py"

class ImageAdapter:
    requires_numpy = False
    shared = False

    def cases(self, args):
        return args.sizes

    def candidates(self, args):
        return args.candidates

    def prepare(self, ctx, size):
        from resource_planner.workflows import motion
        folder = ctx.root / str(size)
        if self.name == 'motion_correct':
            image = motion.generate(folder, size)
            truth = motion.image(folder / 'truth.mrc', size)
        else:
            folder.mkdir()
            image = folder / 'micrograph.mrc'
            generator = Path(__file__).resolve().parents[1] / 'fixtures' / 'micrographs.py'
            ctx.command([ctx.args.fixture_python, str(generator), str(image), '--size', str(size)],
                        folder, 'fixture', timeout=300)
            truth = json.loads(image.with_suffix('.truth.json').read_text())
        project = ctx.project(folder, self.project_name, image, particles=False,
                              ctf_flag='no' if self.name == 'motion_correct' else 'yes')
        files = [image, folder / 'truth.mrc'] if self.name == 'motion_correct' else [image, image.with_suffix('.truth.json')]
        return {'folder': folder, 'project': project, 'truth': truth, 'size': size, 'files': files}

    def jobs(self, ctx, fixture, directory, nthr, workers):
        jobs = []
        for index in range(workers):
            work = directory / f'worker{index}'
            work.mkdir()
            command = [str(ctx.args.private_exec), f'prg={self.name}',
                       f"projfile={fixture['project']}", 'mkdir=no',
                       *[f'{k}={v}' for k, v in self.science.items()],
                       'fromp=1', 'top=1', 'part=1', 'nparts=1', f'nthr={nthr}',
                       f'outfile={self.project_name}.simple', 'memreport=yes', 'memreport_interval=1']
            jobs.append((command, work))
        return jobs

    def validate(self, ctx, fixture, work):
        return self.validator(work, fixture['truth'], fixture['size'])

    def metadata(self, args):
        return {'sizes': args.sizes, 'candidates': args.candidates, 'workload': self.workload}

    def export(self, ctx, plan, records, summary):
        from resource_planner.workflows import ctf
        if self.name == 'ctf_estimate':
            return ctf.export_profile(ctx.root, plan, records, summary)
        if plan['repetitions'] < 3 or not any(w == 1 for _, w in plan['candidates']):
            return None
        return write_profile(ctx.root, simple_exec=ctx.args.simple_exec,
                             simple_path=ctx.env.get('SIMPLE_PATH'))


def export_shared(adapter, ctx, plan, records):
    expected={(size,nthr,repeat) for size in adapter.cases(ctx.args)
              for nthr,_ in adapter.candidates(ctx.args) for repeat in range(plan['repetitions']+1)}
    actual=[(r['size'],r['nthr'],r['repeat']) for r in records]
    if set(actual)!=expected or len(actual)!=len(expected):raise ValueError('incomplete/duplicate calibration observations')
    summary=[]
    for size in adapter.cases(ctx.args):
        for nthr,_ in adapter.candidates(ctx.args):
            rows=[r for r in records if r['size']==size and r['nthr']==nthr and not r['warmup']]
            if len(rows)!=plan['repetitions']:raise ValueError('incomplete calibration observations')
            summary.append({'size':size,'nthr':nthr,'nparts':1,'ncunits':1,
                'median_wall_seconds':statistics.median(r['wall_seconds'] for r in rows),
                'max_rss_mib':max(r['peak_rss_mib'] for r in rows),
                'quality_qualified':all(r['validation']['quality_qualified'] for r in rows)})
    (ctx.root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if plan['repetitions']<3 or not all(r['passed'] and r['validation']['quality_qualified'] for r in records):return None
    for path,digest in plan['binaries'].items():
        if file_hash(Path(path))!=digest:raise ValueError('executable changed during calibration')
    settings={}
    for size in adapter.cases(ctx.args):
        rows=[r for r in summary if r['size']==size]
        fastest=min(r['median_wall_seconds'] for r in rows)
        selected=min((r for r in rows if r['median_wall_seconds']<=fastest*1.05),key=lambda r:r['nthr'])
        settings[str(size)]={'nthr':selected['nthr'],'nparts':1,'ncunits':1,
            'job_memory_per_task_mib':max(1024,math.ceil(selected['max_rss_mib']*1.5/256)*256),
            'median_wall_seconds':selected['median_wall_seconds']}
    profile={'schema_version':1,'profile_type':'simple_resource_profile','commander':adapter.name,
        'auto_apply':False,'qualification':adapter.qualification,'hardware':plan['hardware'],
        'cpu_ids':plan['cpus'],'binary_sha256':plan['binaries'],'workload':plan['workload'],
        'sizes':settings,'evidence':str(ctx.root),'seed':20261015,
        'resource_policy':{'execution_mode':'shared','partition_policy_measured':False,
            'reserve_memory_min_mib':4096,'reserve_memory_fraction':0.1,
            'explicit_job_parameters_take_precedence':True}}
    path=ctx.root/'desktop_profile.json'
    path.write_text(json.dumps(profile,indent=2,allow_nan=False)+'\n');return path


class SpatialAdapter:
    shared=True
    requires_numpy=True
    rate=None

    def candidates(self,args):return [(n,1) for n in args.threads]

    def export(self,ctx,plan,records,summary):return export_shared(self,ctx,plan,records)

    def check_with_numpy(self,ctx,fixture,work,action,path):
        result=work/'spatial_validation.json'
        ctx.command([ctx.args.fixture_python,str(GENERATOR),action,str(fixture['input']),
                     '--input',str(path),'--result',str(result)],work,'validate',timeout=300)
        return json.loads(result.read_text())


