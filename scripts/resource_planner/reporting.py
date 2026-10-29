"""Administrator-facing usage outline and qualified recommendation table."""
from pathlib import Path


def print_outline(args):
    print('SIMPLE hardware calibration', flush=True)
    print('  Standard run:    ./hardware_resource_planner.py')
    print('  Full scope:      4K/8K, 2,000 particles, 5 measured repetitions, broad thread sweep; may take hours')
    print('  Full sweep:      ./hardware_resource_planner.py calibrate --full --execute')
    print('  Preview:         ./hardware_resource_planner.py calibrate')
    print('  Slurm run:       ./hardware_resource_planner.py calibrate --slurm csbdevel --execute')
    print('  Select workflows: --workflows motion_correct ctf_estimate pick solve2D solve3D')
    print('  Standard scope:   4K images, 100 particles, 128-pixel boxes; shortened 2D schedule')
    print('  Workload overrides: --sizes 4096 8192 --particles 2000 --box 128')
    print('  Quick scope:     1 warmup + 1 measured run/thread; known-volume single-stage solve3D')
    print('  Fixtures:         --fixture-source simulated-workflow --molecule 6vxx (or 1jyx)')
    print('                    External Python/NumPy generation; no native simulator or rebuild required.')
    print('  Quality fixtures: Gaussian segmentation picking; asymmetric known-class solve2D (refine=prob)')
    print('  Legacy fixtures:  --fixture-source analytic')
    print('  Molecular runs:   --threads controls candidates; full motion/CTF also sweeps 1/2/4 workers within the CPU/input budget')
    print('  Preparation:      serial prerequisite stages; --fixture-timeout SECONDS (standard 900; full 7200)')
    print('                    Preparation and validation are excluded from reported run times.')
    print('  CPU candidates:   --candidates 4x1 4x2 --threads 4 8 --cpus N')
    print('  Limits:           --repeats N --timeout SECONDS --disk-gib N')
    print('  solve3D limit:    --solve3d-timeout SECONDS (standard 1800; full 21600; explicit --timeout overrides the default)')
    print('  Slurm scaling:    --full compares 1 and 2 nodes; --node-counts 1 2 4 adds four nodes')
    print('                    Picking stays single-node; shared executable storage is required.')
    print('  Slurm resources:  --slurm-cpus N --slurm-memory MIB --slurm-time HH:MM:SS --slurm-account NAME')
    print('  Storage comparison: --io-dirs /usr/local/data /usr/local/cache')
    print('                      --io-size-mib 256 --io-repeats 3 --io-small-files 128')
    print('  Output:           --workdir DIRECTORY --output FILE.json')
    print('  All options:      ./hardware_resource_planner.py calibrate --help')
    print('  Overrides above follow the calibrate subcommand; add --execute to run.')
    print('\nCalibration preset: ' + ('full' if args.full else 'standard (reduced cost)'))
    print('Selected workflows: ' + ', '.join(args.workflows))
    print('Evidence directory: ' + str(args.workdir.absolute()))
    print('JSON destination:   ' + str(args.output.absolute()), flush=True)


def result_table(collection):
    rows = [['Workflow', 'Size/box', 'Mode', 'Nodes', 'nthr', 'Workers', 'nparts', 'MiB/worker', 'Status']]
    reasons = []
    for name, result in collection['workflows'].items():
        profile = result.get('profile')
        if not profile or result['status'] not in ('initial', 'provisional'):
            rows.append([name, '-', '-', '-', '-', '-', '-', '-', result['status']])
            if result.get('reason'):
                reasons.append(name + ': ' + result['reason'])
            continue
        sizes = profile.get('sizes', {})
        if 'settings' in profile:
            sizes = {str(profile.get('workload', {}).get('box', '-')): profile['settings']}
        for size, settings in sizes.items():
            modes = {key: settings[key] for key in ('batch', 'latency') if key in settings}
            if not modes:
                modes = {profile.get('resource_policy', {}).get('execution_mode', 'shared'): settings}
            for mode, selected in modes.items():
                rows.append([name, str(size), mode, str(selected.get('nodes', 1)), str(selected['nthr']),
                             str(selected['ncunits']), str(selected.get('nparts', 'auto*')),
                             str(selected['job_memory_per_task_mib']), result['status']])
    widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]))]
    lines = ['  '.join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip() for row in rows]
    lines.insert(1, '  '.join('-' * width for width in widths))
    qualified = any(r.get('profile') and r['status'] in ('initial', 'provisional') for r in collection['workflows'].values())
    notes = ['\nProvisional = passed synthetic quality gates; real workloads remain unqualified.'] if qualified else ['\nNo qualified resource recommendations were produced.']
    notes.append('Initial = short synthetic screening; repeat with --full for stronger timing evidence.')
    if any('auto*' in row for row in rows):
        notes.append('auto* = min(active inputs, 2 x workers); partition heuristic, not a measured optimum.')
    title = 'Calibration recommendations (best tested settings)' if qualified else 'Calibration status'
    return '\n'.join(['\n'+title, *lines, *notes, *reasons])


def print_results(collection, output):
    print(result_table(collection))
    from resource_planner.storage import report
    report(collection.get('storage_benchmarks', []))
    for name, result in collection['workflows'].items():
        if result.get('scaling'):
            print('\nNode scaling: ' + name)
            print('Size  Nodes  Seconds  Speedup  Efficiency  Node-hours')
            for row in result['scaling']:
                speed = '-' if row['speedup'] is None else f"{row['speedup']:.2f}"
                efficiency = '-' if row['parallel_efficiency'] is None else f"{row['parallel_efficiency']:.1%}"
                print(f"{row['size']}  {row['nodes']}  {row['seconds']:.3f}  {speed}  {efficiency}  {row['node_hours']:.5f}")
        if result.get('distributed_status'):
            print(name + ': ' + result['distributed_status'])
        for nodes, trial in result.get('node_comparison', {}).items():
            if trial.get('reason'):
                print(f"{name}, {nodes} nodes: {trial['reason']}")
    print('\nMeasurements (warmups excluded):')
    print('Workflow  Size/box  Nodes  nthr  Workers  Seconds  Peak MiB  Quality')
    for name, result in collection['workflows'].items():
        for row in result.get('measurements', []):
            if row.get('warmup'):
                continue
            checks = row.get('validation', {})
            checks = checks if isinstance(checks, list) else [checks]
            quality = 'passed' if checks and all(c.get('quality_qualified') is True for c in checks) else 'unqualified'
            peak = row.get('peak_rss_mib', row.get('max_worker_mib', 0))
            print(f"{name}  {row['size']}  {row.get('nodes', 1)}  {row['nthr']}  {row['workers']}  {row['wall_seconds']:.3f}  {peak:.1f}  {quality}")
    qualified = any(r.get('profile') and r['status'] in ('initial', 'provisional') for r in collection['workflows'].values())
    if qualified:
        print('\nCreated JSON profile: ' + str(Path(output).absolute()))
        print('This JSON records hardware, measurements and qualified settings for future SIMPLE jobs.')
    else:
        print('\nCreated JSON diagnostic report: ' + str(Path(output).absolute()))
        print('This JSON records workflow status and errors; it contains no qualified settings for future jobs.')
    print('SIMPLE does not read/apply this file automatically yet (auto_apply=false).', flush=True)
