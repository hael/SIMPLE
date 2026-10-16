#!/usr/bin/env python3
"""Run standard SIMPLE calibration with no arguments; use subcommands for overrides."""

import os
import shutil
import sys

# Some cluster login environments still call Python 3.6 "python3". Keep this
# bootstrap compatible with it, and use an installed modern interpreter.
def ensure_python():
    if sys.version_info >= (3, 9):
        return
    for minor in range(14, 8, -1):
        executable = shutil.which('python3.' + str(minor))
        if executable:
            os.execv(executable, [executable, *sys.argv])
            return
    raise SystemExit('Calibration requires Python 3.9 or newer on PATH.')


ensure_python()

from resource_planner.cli import main


if __name__ == '__main__':
    sys.exit(main())
