"""Process launching, resource limits, evidence writes and cleanup."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from resource_planner.campaign import disk_usage

def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def run_jobs(jobs, env, cpus, timeout, root, disk_limit, label='run'):
    """Run one process group per job; always reap and kill remaining descendants."""
    processes, logs = [], []
    started = time.monotonic()
    next_disk_check = started
    try:
        for command, work in jobs:
            write_json(work / f'{label}_command.json', command)
            log = (work / f'{label}.log').open('wb')
            logs.append(log)
            supervised = env.get('SIMPLE_CALIBRATION_SUPERVISE') == '1'
            launched = [sys.executable, str(Path(__file__).with_name('supervise.py')), *command] if supervised else command
            argv = ['taskset', '-c', ','.join(map(str, cpus)), *launched] if cpus else launched
            processes.append(subprocess.Popen(argv, cwd=work, env=env, stdout=log,
                                              stderr=subprocess.STDOUT, start_new_session=True))
        while True:
            now = time.monotonic()
            if now >= next_disk_check:
                if disk_usage(root) > disk_limit:
                    raise RuntimeError('calibration disk budget exceeded')
                next_disk_check = now + 0.5
            statuses = [p.poll() for p in processes]
            for index, code in enumerate(statuses):
                if code is not None and code != 0:
                    logfile = jobs[index][1] / f'{label}.log'
                    raise RuntimeError(f'calibration process failed ({code}); inspect {logfile}')
            if all(code is not None for code in statuses):
                break
            if now - started >= timeout:
                raise TimeoutError(f'calibration process timeout ({timeout}s)')
            time.sleep(0.01)
        elapsed = time.monotonic() - started
        if disk_usage(root) > disk_limit:
            raise RuntimeError('calibration disk budget exceeded')
        return elapsed
    finally:
        for process in processes:
            try:
                os.killpg(process.pid, signal.SIGTERM if env.get('SIMPLE_CALIBRATION_SUPERVISE') == '1' else signal.SIGKILL)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        for log in logs:
            log.close()

