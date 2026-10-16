"""Submit the existing calibration adapters as one Slurm allocation."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


def add_options(parser):
    parser.add_argument('--slurm', metavar='PARTITION', help='submit calibration to this Slurm partition')
    parser.add_argument('--slurm-cpus', type=int, default=16, help='CPUs allocated to one calibration task')
    parser.add_argument('--slurm-memory', type=int, default=65536, metavar='MIB', help='allocation memory in MiB')
    parser.add_argument('--slurm-time', default='01:00:00', help='Slurm wall-time limit')
    parser.add_argument('--slurm-account', help='optional Slurm account')
    parser.add_argument('--slurm-worker', action='store_true', help=__import__('argparse').SUPPRESS)


def allocation():
    if not os.environ.get('SLURM_JOB_ID'):
        return {'scheduler': 'local'}
    return {'scheduler': 'slurm', 'job_id': os.environ['SLURM_JOB_ID'],
            'partition': os.environ.get('SLURM_JOB_PARTITION'),
            'nodes': os.environ.get('SLURM_JOB_NODELIST'),
            'cpus_per_task': os.environ.get('SLURM_CPUS_PER_TASK'),
            'memory_per_node_mib': os.environ.get('SLURM_MEM_PER_NODE'),
            'memory_per_cpu_mib': os.environ.get('SLURM_MEM_PER_CPU'),
            'scope': 'single-node calibration; distributed scaling unmeasured'}


def submit(args):
    if args.slurm_cpus < 1 or args.slurm_memory < 1:
        raise ValueError('Slurm CPU and memory requests must be positive')
    if args.slurm_worker or os.environ.get('SLURM_JOB_ID'):
        raise ValueError('nested Slurm submission is not supported; run calibrate directly in the allocation')
    if args.output.exists() or args.workdir.exists():
        raise ValueError('output profile and work directory must not already exist')
    command = [sys.executable, '-u', str(Path(__file__).resolve().parents[1] / 'hardware_resource_planner.py'),
               'calibrate', '--execute', '--slurm-worker']
    for key, value in vars(args).items():
        if key == 'discard_images' and value is False:
            command.append('--no-discard-images')
            continue
        if key in ('action', 'execute', 'legacy', 'slurm') or key.startswith('slurm_') or value is None or value is False:
            continue
        command.append('--' + key.replace('_', '-'))
        if value is True:
            continue
        values = value if isinstance(value, list) else [value]
        command.extend(str(v.absolute() if isinstance(v, Path) else v) for v in values)
    # srun establishes an allocation-bound affinity before resource discovery.
    script = '#!/bin/bash\nset -euo pipefail\nexec ' + shlex.join([
        'srun', '--ntasks=1', '--cpu-bind=cores', f'--cpus-per-task={args.slurm_cpus}', *command]) + '\n'
    argv = ['sbatch', '--parsable', '--nodes=1', '--ntasks=1', f'--cpus-per-task={args.slurm_cpus}',
            f'--mem={args.slurm_memory}M', f'--time={args.slurm_time}', f'--partition={args.slurm}',
            '--export=ALL', '--job-name=simple-calibration', f'--chdir={Path.cwd()}',
            f'--output={args.output.absolute()}.slurm-%j.log']
    if args.slurm_account:
        argv.append('--account=' + args.slurm_account)
    preview = {'execution_enabled': args.execute, 'submission': argv, 'script': script,
               'output': str(args.output.absolute()), 'resource_discovery': 'inside the Slurm task'}
    if not args.execute:
        print(json.dumps(preview, indent=2))
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(argv, input=script, text=True, capture_output=True, check=True)
    preview['job_id'] = result.stdout.strip()
    receipt = Path(str(args.output) + '.submission.json')
    receipt.write_text(json.dumps(preview, indent=2) + '\n')
    print('Submitted Slurm calibration job: ' + preview['job_id'])
    print('The recommendation table will appear in: ' + str(args.output.absolute()) + '.slurm-' + preview['job_id'].split(';')[0] + '.log')
    print('The job will create JSON for future SIMPLE use at: ' + str(args.output.absolute()))
    print('Submission receipt: ' + str(receipt.absolute()))
    return 0
