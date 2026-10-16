"""Resource profile identity and legacy motion JSON export."""
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import struct
from resource_planner.campaign import file_hash

def machine_identity():
    models=[line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
            if line.startswith('model name')]
    return {'cpu_model':models[0] if models else 'unknown',
            'logical_cpus':os.cpu_count(),'architecture':platform.machine(),'kernel':platform.release()}


def elf_signature(path):
    """Compare scientific code/data across CMake's build/install RPATH rewrite."""
    data=path.read_bytes()
    if data[:6] != b'\x7fELF\x02\x01':
        raise ValueError('profile qualification requires little-endian ELF64 executables')
    offset=struct.unpack_from('<Q',data,40)[0]
    stride,count,names=struct.unpack_from('<HHH',data,58)
    headers=[struct.unpack_from('<IIQQQQIIQQ',data,offset+i*stride) for i in range(count)]
    entry=headers[names]; table=data[entry[4]:entry[4]+entry[5]]
    result={}
    for entry in headers:
        name=table[entry[0]:].split(b'\0')[0].decode()
        if name in ('.text','.rodata','.data'):
            result[name]=hashlib.sha256(data[entry[4]:entry[4]+entry[5]]).hexdigest()
    if len(result)!=3:
        raise ValueError('missing ELF code/data sections')
    return result


def write_profile(results, installation=None, *, simple_exec=None, simple_path=None):
    plan=json.loads((results/'plan.json').read_text())
    summary=json.loads((results/'summary.json').read_text())
    observations=[json.loads(line) for line in (results/'results.jsonl').read_text().splitlines()]
    if plan['repetitions'] < 3:
        raise ValueError('at least three repetitions are required for a profile')
    for row in summary:
        measurements=[r for r in observations if not r['warmup'] and
                      (r['size'],r['nthr'],r['workers'])==(row['size'],row['nthr'],row['workers'])]
        if len(measurements)!=plan['repetitions'] or any(not r['passed'] for r in measurements):
            raise ValueError('incomplete or failed measurements')
        for measurement in measurements:
            if len(measurement['validation'])!=row['workers']:
                raise ValueError('missing worker validation')
            for check in measurement['validation']:
                # Export additionally gates global recovery against generated
                # translation truth, using the aligner's 0.1-pixel scale.
                if (len(check['shifts_x'])!=8 or len(check['shifts_y'])!=8 or
                        max(abs(x+t) for x,t in zip(check['shifts_x'],plan['known_shifts_x']))>0.1 or
                        max(abs(y) for y in check['shifts_y'])>0.1):
                    raise ValueError('known-motion recovery exceeds 0.1 pixel')
    signatures={}
    for original,digest in plan['binaries'].items():
        installed=installation/'bin'/Path(original).name if installation is not None else Path(original)
        if file_hash(Path(original))!=digest or (installed != Path(original) and
                elf_signature(Path(original))!=elf_signature(installed)):
            raise ValueError('build/install executable mismatch')
        signatures[str(installed)]=file_hash(installed)
    profiles={}
    for size in plan['sizes']:
        rows=[r for r in summary if r['size']==size]
        best=max(r['median_movies_per_second'] for r in rows)
        near=[r for r in rows if r['median_movies_per_second']>=best/1.05]
        batch=min(near,key=lambda r:(r['nthr']*r['workers'],r['max_worker_mib']*r['workers']))
        latency=min((r for r in rows if r['workers']==1),key=lambda r:r['median_wave_seconds'])
        def settings(row):
            return {'nthr':row['nthr'],'ncunits':row['workers'],
                    'job_memory_per_task_mib':max(1024,math.ceil(row['max_worker_mib']*1.5/256)*256),
                    'median_movies_per_second':row['median_movies_per_second'],
                    'median_wave_seconds':row['median_wave_seconds']}
        profiles[str(size)]={'batch':settings(batch),'latency':settings(latency)}
    profile={'schema_version':1,'profile_type':'simple_resource_profile',
             'commander':'motion_correct','auto_apply':False,
             'resource_policy':{'partition_waves':2,'partition_policy_measured':False,
                                'reserve_memory_min_mib':4096,'reserve_memory_fraction':0.1,
                                'explicit_job_parameters_take_precedence':True},
             'qualification':'provisional synthetic cached-input worker measurements; master unmeasured',
             'hardware':machine_identity(),
             'cpu_ids':plan['cpus'],'workload':{'frames':8,'smpd':1.3,'algorithm':'patch','mcpatch':'yes'},
             'sizes':profiles,'simple_exec':str(installation/'bin/simple_exec' if installation is not None else simple_exec),
             'simple_path':str(installation) if installation is not None else simple_path,'binary_sha256':signatures,'seed':20261015,
             'evidence':str(results),'partition_policy':'min(active_movies, 2 * ncunits); unmeasured heuristic',
             'memory_policy':'1.5 times observed maximum worker RSS, rounded up to 256 MiB; reserve >=4 GiB'}
    destination=results/'desktop_profile.json'
    with destination.open('x') as f:
        json.dump(profile,f,indent=2);f.write('\n')
    return destination

