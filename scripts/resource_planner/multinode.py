"""Fixed-workload node comparisons using native SIMPLE partition masters."""
import copy
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

SUPPORTED = {'motion_correct', 'ctf_estimate', 'solve2D', 'solve3D'}


def hosts():
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('node comparison requires a Slurm allocation')
    return subprocess.check_output(['scontrol', 'show', 'hostnames', os.environ['SLURM_JOB_NODELIST']], text=True).split()


def prepare_launch(ctx, directory, jobs, nthr, env):
    args = ctx.args
    selected = args._node_hosts[:args._nodes]
    config = {'hosts': selected, 'threads': nthr, 'cpus_per_node': args.slurm_cpus,
              'job_id': os.environ['SLURM_JOB_ID']}
    (directory / 'launch.json').write_text(json.dumps(config))
    shim = directory / '_launcher'
    shim.mkdir()
    launcher = shim / 'sbatch'
    launcher.write_text('#!/bin/bash\nexec ' + shlex.join([sys.executable, str(Path(__file__).with_name('node_launch.py'))]) + ' "$@"\n')
    launcher.chmod(0o755)
    env.update(PATH=str(shim) + os.pathsep + env['PATH'],
               SIMPLE_CALIBRATION_LAUNCH=str(directory), SIMPLE_CALIBRATION_SUPERVISE='1')
    result = []
    for command, work in jobs:
        project = next(a.split('=', 1)[1] for a in command if a.startswith('projfile='))
        ctx.command([str(args.simple_exec), 'prg=update_project', 'projfile=' + project,
                     'qsys_name=slurm', 'mkdir=no'], work, 'configure_slurm')
        command = [a for a in command if a.split('=', 1)[0] not in ('qsys_name', 'nparts', 'ncunits', 'part')]
        command += ['qsys_name=slurm', 'nparts=' + str(max(args.node_counts)), 'ncunits=' + str(args._nodes)]
        result.append((command, work))
    return result


def placement(directory, expected):
    rows = [json.loads(p.read_text()) for p in directory.glob('*.placement.json')]
    config = json.loads((directory / 'launch.json').read_text())
    actual = {r.get('slurm_node') for r in rows if r.get('exit_code') == 0}
    if not rows or any(r.get('exit_code') != 0 or r.get('job_id') != config['job_id']
                       or not str(r.get('step_id', '')).isdigit() for r in rows) or actual != set(expected):
        raise ValueError('native workers did not complete successfully on every requested node')
    return rows


def compare(adapter, args, root, cpus, runner):
    """Keep each node run and failure independently reviewable; select only qualified runs."""
    root.mkdir()
    results = {}
    for count in args.node_counts:
        trial = copy.copy(args)
        trial._nodes = count
        trial._simulation_cache = {}
        trial.workdir = root / ('nodes' + str(count))
        trial.workdir.mkdir()
        try:
            result = runner(adapter, trial, trial.workdir / adapter.name, cpus)
        except Exception as error:
            result = {'status': 'failed', 'profile': None, 'reason': str(error), 'measurements': []}
        results[str(count)] = result
        (root / 'node_results.json').write_text(json.dumps(results, indent=2) + '\n')
    qualified = {int(n): r for n, r in results.items() if r.get('profile')}
    if not qualified:
        return {'status': 'failed', 'profile': None, 'reason': 'No node count passed placement and scientific gates',
                'node_comparison': results, 'measurements': []}
    chosen = copy.deepcopy(qualified[min(qualified)])
    profile = chosen['profile']
    def settings(p):
        return p.get('sizes', {'default': p.get('settings', {})})
    comparisons = []
    for size in settings(profile):
        options = {n: settings(r['profile'])[size] for n, r in qualified.items() if size in settings(r['profile'])}
        fastest = min(v['median_wall_seconds'] for v in options.values())
        count = min(n for n, v in options.items() if v['median_wall_seconds'] <= fastest * 1.05)
        best = copy.deepcopy(options[count])
        best.update(nodes=count, nparts=max(args.node_counts), ncunits=count)
        baseline = options.get(1, {}).get('median_wall_seconds')
        for n, value in options.items():
            seconds = value['median_wall_seconds']
            comparisons.append({'size': size, 'nodes': n, 'seconds': seconds,
                                'speedup': baseline / seconds if baseline else None,
                                'parallel_efficiency': baseline / (seconds*n) if baseline else None,
                                'node_hours': n*seconds/3600})
        if 'sizes' in profile:
            profile['sizes'][size] = best
        else:
            profile['settings'] = best
    profile['qualification'] += '; verified allocation-local native worker placement; fixed partition count'
    profile['resource_policy']['execution_mode'] = 'distributed'
    chosen['node_comparison'] = results
    chosen['scaling'] = comparisons
    chosen['measurements'] = [{**row, 'nodes': int(n)} for n, r in results.items() for row in r.get('measurements', [])]
    return chosen
