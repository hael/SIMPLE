"""Command-line interface for SIMPLE resource calibration."""
import argparse
import json
from pathlib import Path
import signal
import subprocess
import sys
from resource_planner.campaign import execute, load_plan, render
from resource_planner.calibration import add_options, calibrate
from resource_planner.hardware import discover
from resource_planner.reporting import print_outline

def main(argv=None):
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    commands.add_parser("discover")
    calibration = commands.add_parser("calibrate", help="calibrate supported workflows; dry-run by default")
    add_options(calibration)
    campaign = commands.add_parser("campaign", help="dry-run by default; execution needs --execute")
    campaign.add_argument("manifest", type=Path)
    campaign.add_argument("--output", required=True, type=Path)
    campaign.add_argument("--execute", action="store_true", help="run preparation, SIMPLE and validators")
    argv = list(sys.argv[1:] if argv is None else argv)
    # Invoking the tool itself is the administrator's standard-run request.
    args = parser.parse_args(argv or ['calibrate', '--execute'])
    try:
        if args.action == "calibrate":
            if args.execute:
                print_outline(args)
            return calibrate(args)
        if args.action == "discover":
            print(json.dumps(discover(), indent=2))
            return 0
        plan = load_plan(args.manifest)
        if args.execute:
            return execute(plan, args.output)
        print(json.dumps(render(plan, args.output.resolve()), indent=2))
        return 0
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"resource planner: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("resource planner: interrupted; the current attempt is retained", file=sys.stderr)
        return 130
