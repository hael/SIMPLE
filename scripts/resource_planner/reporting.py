"""Administrator-facing usage outline and qualified recommendation table."""
from pathlib import Path


def print_outline(args):
    print('SIMPLE hardware calibration', flush=True)
    print('  Standard run:    ./hardware_resource_planner.py')
    print('  Full sweep:      ./hardware_resource_planner.py calibrate --full --execute')
    print('  Preview:         ./hardware_resource_planner.py calibrate')
    print('  Slurm run:       ./hardware_resource_planner.py calibrate --slurm csbdevel --execute')
    print('  Select workflows: --workflows motion_correct ctf_estimate pick solve2D solve3D')
    print('  Standard scope:   4K images, 500 particles, 128-pixel boxes; shortened 2D schedule')
    print('  Workload overrides: --sizes 4096 8192 --particles 2000 --box 128')
    print('  Fixtures:         --fixture-source simulated-workflow --molecule 6vxx (or 1jyx)')
    print('                    External Python/NumPy generation; no native simulator or rebuild required.')
    print('  Legacy fixtures:  --fixture-source analytic')
    print('  Molecular runs:   --threads controls candidates; one local job/partition at a time')
    print('  Preparation:      serial prerequisite stages; --fixture-timeout SECONDS (default 900)')
    print('                    Preparation and validation are excluded from reported run times.')
    print('  CPU candidates:   --candidates 4x1 4x2 --threads 4 8 --cpus N')
    print('  Limits:           --repeats N --timeout SECONDS --disk-gib N')
    print('  solve3D limit:    --solve3d-timeout SECONDS (default 1800; explicit --timeout overrides the default)')
    print('  Slurm resources:  --slurm-cpus N --slurm-memory MIB --slurm-time HH:MM:SS --slurm-account NAME')
    print('  Output:           --workdir DIRECTORY --output FILE.json')
    print('  All options:      ./hardware_resource_planner.py calibrate --help')
    print('  Overrides above follow the calibrate subcommand; add --execute to run.')
    print('\nCalibration preset: ' + ('full' if args.full else 'standard (reduced cost)'))
    print('Selected workflows: ' + ', '.join(args.workflows))
    print('Evidence directory: ' + str(args.workdir.absolute()))
    print('JSON destination:   ' + str(args.output.absolute()), flush=True)


def result_table(collection):
    rows = [['Workflow', 'Size/box', 'Mode', 'nthr', 'Workers', 'nparts', 'MiB/worker', 'Status']]
    reasons = []
    for name, result in collection['workflows'].items():
        profile = result.get('profile')
        if not profile or result['status'] != 'provisional':
            rows.append([name, '-', '-', '-', '-', '-', '-', result['status']])
            if result.get('reason'):
                reasons.append(name + ': ' + result['reason'])
            continue
        sizes = profile.get('sizes', {})
        if 'settings' in profile:
            sizes = {str(profile.get('workload', {}).get('box', '-')): profile['settings']}
        for size, settings in sizes.items():
            modes = {key: settings[key] for key in ('batch', 'latency') if key in settings}
            if not modes:
                modes = {'shared': settings}
            for mode, selected in modes.items():
                rows.append([name, str(size), mode, str(selected['nthr']),
                             str(selected['ncunits']), str(selected.get('nparts', 'auto*')),
                             str(selected['job_memory_per_task_mib']), result['status']])
    widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]))]
    lines = ['  '.join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip() for row in rows]
    lines.insert(1, '  '.join('-' * width for width in widths))
    qualified = any(r.get('profile') and r['status'] == 'provisional' for r in collection['workflows'].values())
    notes = ['\nProvisional = passed synthetic quality gates; real workloads remain unqualified.'] if qualified else ['\nNo qualified resource recommendations were produced.']
    if any('auto*' in row for row in rows):
        notes.append('auto* = min(active inputs, 2 x workers); partition heuristic, not a measured optimum.')
    title = 'Calibration recommendations (best tested settings)' if qualified else 'Calibration status'
    return '\n'.join(['\n'+title, *lines, *notes, *reasons])


def print_results(collection, output):
    print(result_table(collection))
    qualified = any(r.get('profile') and r['status'] == 'provisional' for r in collection['workflows'].values())
    if qualified:
        print('\nCreated JSON profile: ' + str(Path(output).absolute()))
        print('This JSON records hardware, measurements and qualified settings for future SIMPLE jobs.')
    else:
        print('\nCreated JSON diagnostic report: ' + str(Path(output).absolute()))
        print('This JSON records workflow status and errors; it contains no qualified settings for future jobs.')
    print('SIMPLE does not read/apply this file automatically yet (auto_apply=false).', flush=True)
