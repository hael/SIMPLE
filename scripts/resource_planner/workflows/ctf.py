"""CTF estimation truth validation and profile selection."""
import json
import math
from pathlib import Path
import re
from resource_planner.campaign import file_hash
from resource_planner.profiles import machine_identity

SCIENCE = {'ctfpatch':'yes','pspecsz':512,'fit_phshift':'no',
           'dfmin':0.5,'dfmax':5.0,'hp':30,'lp':5}


def validate(directory, truth, size):
    docs=list(directory.glob('*_ctf.txt'))
    if len(docs)!=1:
        raise ValueError(f'expected one CTF document, found {docs}')
    values=dict(re.findall(r'(\w+)=\s*([^\s]+)', docs[0].read_text().splitlines()[0]))
    fitted={k:float(values[k]) for k in ('dfx','dfy','ctfres','ctfcc')}
    error=max(abs(fitted[k]-truth[k]) for k in ('dfx','dfy'))
    return {'passed':all(math.isfinite(v) for v in fitted.values()) and error<=0.1,
            'max_defocus_error_um':error, **fitted}


def export_profile(root, plan, records, summary):
    if plan['repetitions']<3 or not all(r['passed'] for r in records):
        return None
    expected={(size,nthr,workers,repeat) for size in plan['sizes']
              for nthr,workers in plan['candidates'] for repeat in range(plan['repetitions']+1)}
    actual=[(r['size'],r['nthr'],r['workers'],r['repeat']) for r in records]
    if set(actual)!=expected or len(actual)!=len(expected):
        raise ValueError('incomplete or duplicate benchmark observations')
    sizes={}
    for size in plan['sizes']:
        rows=[r for r in summary if r['size']==size]
        singles=[r for r in rows if r['workers']==1]
        if not singles:
            return None
        best=max(r['median_micrographs_per_second'] for r in rows)
        batch=min((r for r in rows if r['median_micrographs_per_second']>=best/1.05),
                  key=lambda r:(r['nthr']*r['workers'],r['max_worker_mib']*r['workers']))
        latency=min(singles,key=lambda r:r['median_wave_seconds'])
        sizes[str(size)]={name:{'nthr':row['nthr'],'ncunits':row['workers'],
            'job_memory_per_task_mib':max(1024,math.ceil(row['max_worker_mib']*1.5/256)*256),
            'median_micrographs_per_second':row['median_micrographs_per_second'],
            'median_wave_seconds':row['median_wave_seconds']}
            for name,row in (('batch',batch),('latency',latency))}
    for path,digest in plan['binaries'].items():
        if file_hash(Path(path))!=digest:
            raise ValueError('executable changed during benchmark')
    profile={'schema_version':1,'profile_type':'simple_resource_profile',
             'commander':'ctf_estimate','auto_apply':False,
             'qualification':'provisional synthetic cached-input workers; master throughput unmeasured',
             'hardware':machine_identity(),'cpu_ids':plan['cpus'],
             'binary_sha256':plan['binaries'],'workload':plan['workload'],
             'sizes':sizes,'evidence':str(root),'seed':20261015,
             'resource_policy':{'partition_waves':2,'partition_policy_measured':False,
                                'reserve_memory_min_mib':4096,'reserve_memory_fraction':0.1,
                                'explicit_job_parameters_take_precedence':True}}
    path=root/'desktop_profile.json'
    with path.open('x') as stream:
        json.dump(profile,stream,indent=2);stream.write('\n')
    return path



from resource_planner.workflows.common import ImageAdapter
class CtfAdapter(ImageAdapter):
    name = 'ctf_estimate'
    project_name = 'ctf'
    rate = 'micrographs_per_second'
    requires_numpy = True
    science = SCIENCE
    workload = {'smpd': 1.3, 'kv': 300, 'cs': 2.7, 'fraca': 0.1,
                'synthetic_defocus_um': 2.0, **SCIENCE}
    validator = staticmethod(validate)


