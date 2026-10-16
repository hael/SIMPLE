"""Linux subreaper for a local SIMPLE command and its detached queue workers."""
import ctypes
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def main(argv):
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), 'cannot enable calibration subreaper')
    children = Path('/proc/self/task/%d/children' % os.getpid())
    if not children.exists():
        raise RuntimeError('Linux /proc child tracking is required for molecular calibration')
    interrupted = []
    def stop(signum, frame):
        interrupted.append(signum)
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    process = subprocess.Popen(argv)
    try:
        while process.poll() is None and not interrupted:
            time.sleep(0.01)
        code = process.returncode if not interrupted else 128+interrupted[0]
    finally:
        # Killing each generation reparents surviving grandchildren here, even
        # if a queue launcher called setsid. Reap before ending the supervisor.
        while True:
            for pid in children.read_text().split():
                try:
                    os.kill(int(pid), signal.SIGKILL)
                except ProcessLookupError:
                    pass
            try:
                while os.waitpid(-1, os.WNOHANG)[0]:
                    pass
            except ChildProcessError:
                break
            time.sleep(0.01)
    return code if code >= 0 else 128-code


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
