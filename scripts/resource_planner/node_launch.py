"""Allocation-local launcher for SIMPLE's generated Slurm worker scripts.

Installed as `sbatch` ONLY in a calibration candidate's child environment.
Never submits a new allocation. Native masters retain partitioning and merging.
"""
import fcntl
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import uuid


def main(argv):
    root = Path(os.environ['SIMPLE_CALIBRATION_LAUNCH']).resolve()
    config = json.loads((root / 'launch.json').read_text())
    if os.environ.get('SLURM_JOB_ID') != config['job_id']:
        raise ValueError('worker launcher requires its original Slurm allocation')
    if argv[0] == '--worker':
        script, token = argv[1:]
        evidence = {'host': socket.gethostname().split('.')[0], 'script': script,
                    'slurm_node': os.environ.get('SLURMD_NODENAME'),
                    'job_id': os.environ['SLURM_JOB_ID'], 'step_id': os.environ.get('SLURM_STEP_ID')}
        (root / (token + '.placement.json')).write_text(json.dumps(evidence))
        code = subprocess.call(['bash', script])
        evidence['exit_code'] = code
        (root / (token + '.placement.json')).write_text(json.dumps(evidence))
        return code
    if len(argv) != 1 or argv[0].startswith('-'):
        raise ValueError('only single generated worker scripts are supported')
    script = Path(argv[0]).resolve()
    if not script.is_relative_to(root):
        raise ValueError('worker script is outside candidate directory')
    content = script.read_text()
    if '#SBATCH --array' in content:
        raise ValueError('array scripts are not supported by calibration launcher')
    match = re.search(r'^#SBATCH --cpus-per-task[= ]+(\d+)', content, re.M)
    threads = int(match.group(1)) if match else config['threads']
    if not 0 < threads <= config['cpus_per_node']:
        raise ValueError('native worker CPU request exceeds allocation')
    with (root / 'launch.counter').open('a+') as counter:
        fcntl.flock(counter, fcntl.LOCK_EX)
        counter.seek(0)
        index = int(counter.read() or '0')
        counter.seek(0); counter.truncate(); counter.write(str(index + 1))
    host = config['hosts'][index % len(config['hosts'])]
    token = uuid.uuid4().hex
    command = ['srun', '--exclusive', '--exact', '--nodes=1', '--ntasks=1',
               '--nodelist=' + host, '--cpu-bind=cores', '--cpus-per-task=' + str(threads),
               sys.executable, str(Path(__file__).resolve()), '--worker', str(script), token]
    with (root / (token + '.srun.log')).open('wb') as log:
        subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
    print('Started allocation-local worker ' + token)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
