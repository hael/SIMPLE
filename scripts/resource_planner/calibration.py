"""Shared local calibration execution, evidence collection and CLI options."""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import tempfile

from resource_planner.workflows.registry import ADAPTERS
from resource_planner.campaign import file_hash, telemetry
from resource_planner.hardware import cpu_pool, MIB
from resource_planner.profiles import machine_identity
from resource_planner.slurm import add_options as slurm_options, allocation, submit
from resource_planner.reporting import print_results
from resource_planner.execution import run_jobs, write_json

SEED = 20261015


def adapters_for(args):
    if getattr(args, 'fixture_source', 'analytic') == 'simulated-workflow':
        from resource_planner.workflows.molecular import ADAPTERS as molecular
        return {**ADAPTERS, **{name: value for name, value in molecular.items() if name in ADAPTERS}}
    return ADAPTERS


def fixture_fingerprint(fixture):
    """Hash declared scientific inputs, excluding mutable projects and run logs."""
    return {str(path): file_hash(path) for path in fixture.get('files', [])}




class Context:
    def __init__(self, args, root, cpus):
        self.args, self.root, self.cpus = args, root, cpus
        self.disk_root = root if args.legacy else root.parent
        self.env = os.environ.copy()
        self.env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='FALSE',
                        SIMPLE_SEED=str(SEED))
        if args.simple_prefix is not None:
            self.env["SIMPLE_PATH"] = str(args.simple_prefix.resolve())
        self.env.pop('OMP_PLACES', None)

    def command(self, argv, directory, label, timeout=60):
        env = dict(self.env, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
        return run_jobs([(argv, directory)], env, self.cpus, timeout,
                        self.disk_root, self.args.disk_gib * 1024**3, label)

    def project(self, directory, name, image, particles, ctf_flag, oritab=None):
        self.command([str(self.args.simple_exec), 'prg=new_project', f'projname={name}',
                      'qsys_name=local'], directory, 'new_project')
        project = directory / name / f'{name}.simple'
        if particles:
            inputs = [f'stk={image}']
        else:
            table = directory / 'inputs.txt'
            table.write_text(str(image) + '\n')
            inputs = [f'filetab={table}']
        if oritab is not None:
            inputs.extend([f'oritab={oritab}', 'oritype=ptcl2D'])
        self.command([str(self.args.simple_exec),
                      'prg=import_particles' if particles else 'prg=import_movies',
                      f'projfile={project}', *inputs, 'smpd=1.3', 'kv=300', 'cs=2.7',
                      'fraca=0.1', f'ctf={ctf_flag}', 'mkdir=no'], directory / name, 'import')
        return project


def run_workflow(adapter, args, root, cpus):
    root.mkdir(parents=True, exist_ok=False)
    ctx = Context(args, root, cpus)
    seed = getattr(adapter, 'seed', SEED)
    ctx.env['SIMPLE_SEED'] = str(seed)
    if getattr(adapter, 'supervise', False):
        ctx.env['SIMPLE_CALIBRATION_SUPERVISE'] = '1'
    binaries = [args.simple_exec] if adapter.shared and not getattr(adapter, 'requires_private', False) else [args.simple_exec, args.private_exec]
    timeout = args.solve3d_timeout if adapter.name == 'solve3D' else args.timeout
    plan = {'commander': adapter.name, 'cpus': cpus, 'repetitions': args.repeats,
            'candidate_timeout_seconds': timeout,
            'reproducibility': {'seed': seed, 'native_seed_environment': 'SIMPLE_SEED',
                'fixture_generation_threads': 1, 'fixture_reused_across_candidates': True,
                'scope': 'fixed inputs within the recorded generator/build environment; '
                         'parallel solver outputs and timings are not bitwise deterministic',
                'fixtures': {}},
            'hardware': machine_identity(), 'allocation': allocation(),
            'calibration_preset': 'full' if args.full or args.legacy else 'standard', **adapter.metadata(args),
            'binaries': {str(path): file_hash(path) for path in binaries},
            'source_sha256': {str(p.relative_to(Path(__file__).parent)): file_hash(p) for p in sorted(Path(__file__).parent.rglob('*.py'))}}
    write_json(root / 'plan.json', plan)
    records = []
    for size in adapter.cases(args):
        fixture = adapter.prepare(ctx, size)
        fingerprints = fixture_fingerprint(fixture)
        plan['reproducibility']['fixtures'][str(size)] = fingerprints
        write_json(root / 'plan.json', plan)
        for repeat in range(args.repeats + 1):
            candidates = list(adapter.candidates(args))
            random.Random(seed + repeat).shuffle(candidates)
            for nthr, workers in candidates:
                name = f't{nthr}_r{repeat}' if adapter.shared else f't{nthr}_w{workers}_r{repeat}'
                directory = fixture['folder'] / name
                directory.mkdir()
                jobs = adapter.jobs(ctx, fixture, directory, nthr, workers)
                print(f'{adapter.name}: size={size} threads={nthr} workers={workers} repeat={repeat}', flush=True)
                env = {**ctx.env, 'OMP_NUM_THREADS': str(nthr)}
                elapsed = run_jobs(jobs, env, cpus, timeout, ctx.disk_root, args.disk_gib * 1024**3)
                checks, peaks = [], []
                for _, work in jobs:
                    checks.append(adapter.validate(ctx, fixture, work))
                    memory = telemetry(work)
                    if memory['errors'] or not memory['process_peaks'] or (len(memory['process_peaks']) != 1 and not getattr(adapter, 'aggregate_memory', False)):
                        raise ValueError(f'missing/ambiguous memory telemetry in {work}')
                    peaks.append(sum(p['peak_rss_bytes'] for p in memory['process_peaks']) / MIB)
                result = {'size': size, 'nthr': nthr, 'workers': workers, 'repeat': repeat,
                          'warmup': repeat == 0, 'wall_seconds': elapsed,
                          'passed': all(c['passed'] for c in checks)}
                if adapter.shared:
                    result.update(validation=checks[0], peak_rss_mib=peaks[0])
                else:
                    result.update(validation=checks, max_worker_mib=max(peaks),
                                  sum_worker_peaks_mib=sum(peaks), **{adapter.rate: workers / elapsed})
                write_json(directory / 'result.json', result)
                with (root / 'results.jsonl').open('a') as stream:
                    stream.write(json.dumps(result, allow_nan=False) + '\n')
                records.append(result)
                if not result['passed']:
                    raise ValueError(f'output validation failed in {directory}')
                print(f'  {elapsed:.3f}s, peak per process {max(peaks):.1f} MiB', flush=True)
                if args.discard_images and not adapter.shared:
                    for index, (_, work) in enumerate(jobs):
                        if repeat == 0 and index == 0:
                            continue
                        for path in work.iterdir():
                            if path.suffix in ('.mrc', '.jpg'):
                                path.unlink()
        if fixture_fingerprint(fixture) != fingerprints:
            raise ValueError('calibration inputs changed during the candidate sweep')
    summary = []
    if not adapter.shared:
        for size in adapter.cases(args):
            for nthr, workers in adapter.candidates(args):
                rows = [r for r in records if (r['size'], r['nthr'], r['workers']) == (size, nthr, workers) and not r['warmup']]
                summary.append({'size': size, 'nthr': nthr, 'workers': workers,
                                'median_' + adapter.rate: statistics.median(r[adapter.rate] for r in rows),
                                'median_wave_seconds': statistics.median(r['wall_seconds'] for r in rows),
                                'max_worker_mib': max(r['max_worker_mib'] for r in rows)})
        write_json(root / 'summary.json', summary)
    profile_path = adapter.export(ctx, plan, records, summary)
    if profile_path:
        profile = json.loads(profile_path.read_text())
        profile['allocation'] = allocation()
        profile['calibration_preset'] = plan['calibration_preset']
        profile['reproducibility'] = plan['reproducibility']
        profile['seed'] = seed
        write_json(profile_path, profile)
    reason = None
    if not profile_path:
        if args.repeats < 3:
            reason = 'At least three measured repetitions are required for a profile.'
        else:
            invalid = [r for r in records if isinstance(r.get('validation'), dict)
                       and not r['validation'].get('quality_qualified', True)]
            reason = 'Quality gate not met or a required single-worker candidate is missing.'
            if invalid:
                check = invalid[0]['validation']
                reason = check.get('qualification_reason', 'Scientific quality gate failed.')
                if 'adjusted_rand_index' in check:
                    reason += (f" Class agreement ARI={check['adjusted_rand_index']:.3f} (requires >=0.8);"
                               f" occupied classes={check['occupied_classes']} (requires 8).")
    return {'status': 'provisional' if profile_path else 'unqualified',
            'reason': reason,
            'evidence': str(root), 'profile': json.loads(profile_path.read_text()) if profile_path else None}


def add_options(parser, legacy=None):
    slurm_options(parser)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--full', action='store_true', help='run the extensive 4K/8K calibration sweep')
    parser.add_argument('--output', type=Path, required=bool(legacy),
                        default=None if legacy else Path('output/resource_calibration/standard/profile.json'),
                        help='result directory' if legacy else 'combined JSON profile file')
    parser.add_argument('--simple-prefix', type=Path,
                        help='optional installation override; default: current shell environment')
    parser.add_argument('--simple-exec', type=Path)
    parser.add_argument('--private-exec', type=Path)
    parser.add_argument('--fixture-python', default='python3')
    parser.add_argument('--fixture-source', choices=['simulated-workflow', 'analytic'],
                        default='analytic' if legacy else 'simulated-workflow',
                        help='molecular test-suite recipe (default) or legacy analytic fixtures')
    parser.add_argument('--molecule', choices=['6vxx', '1jyx'], default='6vxx')
    parser.add_argument('--fixture-timeout', type=int, default=900,
                        help='timeout for each serial preparation/validation command, outside measured runs')
    parser.add_argument('--sizes', nargs='+', type=int, help='image sizes (standard: 4096; full: 4096 8192)')
    parser.add_argument('--candidates', nargs='+', help='threads x workers, e.g. 4x1 4x8')
    parser.add_argument('--threads', nargs='+', type=int)
    parser.add_argument('--particles', type=int, help='particle count (standard: 500; full: 2000)')
    parser.add_argument('--box', type=int, default=128)
    parser.add_argument('--extr-lim', type=int, help='2D extremal iterations (standard: 4; full: 20)')
    parser.add_argument('--nits-per-stage', type=int, help='2D iterations per stage (standard: 1; full: 5)')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--timeout', type=int, help='per-process timeout seconds (standard: 180; full: 600)')
    parser.add_argument('--solve3d-timeout', type=int,
                        help='solve3D candidate timeout seconds (default: 1800; explicit --timeout applies unless overridden)')
    parser.add_argument('--disk-gib', type=int, default=20)
    parser.add_argument('--cpus', type=int, help='maximum physical cores from the reserved desktop pool')
    parser.add_argument('--discard-images', action=argparse.BooleanOptionalAction, default=not bool(legacy),
                        help='discard repeated image outputs after validation (standard default: enabled)')
    if legacy:
        parser.set_defaults(workflows=[legacy], legacy=True)
    else:
        parser.add_argument('--workdir', type=Path, default=Path('output/resource_calibration/standard/measurements'),
                            help='new evidence directory (default: output/resource_calibration/standard/measurements)')
        parser.add_argument('--workflows', nargs='+', choices=list(ADAPTERS), default=list(ADAPTERS))
        parser.set_defaults(legacy=False)


def resolve_executable(name, explicit=None, prefix=None, required=False):
    """Keep PATH symlinks intact; never assume a source-tree build directory."""
    if explicit is not None:
        path = explicit.absolute()
    elif prefix is not None:
        path = (prefix / 'bin' / name).absolute()
    else:
        found = shutil.which(name)
        path = Path(found).absolute() if found else None
        if path is None and os.environ.get('SIMPLE_PATH'):
            candidate = Path(os.environ['SIMPLE_PATH']) / 'bin' / name
            if candidate.is_file() and os.access(candidate, os.X_OK):
                path = candidate.absolute()
    if required and (path is None or not path.is_file() or not os.access(path, os.X_OK)):
        raise ValueError(f'{name} unavailable: load the SIMPLE installation environment '
                         '(PATH and SIMPLE_PATH) before calibration')
    # An unresolved command name is useful in a preview, but cannot be executed.
    return path if path is not None else Path(name)


def apply_preset(args):
    """Explicit workload flags always override the chosen preset."""
    if args.solve3d_timeout is None:
        args.solve3d_timeout = args.timeout if args.timeout is not None else 1800
    full = args.full or args.legacy
    defaults = {'sizes': [4096, 8192] if full else [4096],
                'particles': 2000 if full else 500,
                'extr_lim': 20 if full else 4, 'nits_per_stage': 5 if full else 1,
                'timeout': (300 if args.legacy and args.workflows != ['solve2D'] else 600) if full else 180}
    for key, value in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)


def configure(args):
    apply_preset(args)
    cpus = cpu_pool()
    if args.cpus is not None:
        if not 1 <= args.cpus <= len(cpus):
            raise ValueError('--cpus exceeds the usable physical-core pool')
        cpus = cpus[:args.cpus]
    args.simple_exec = resolve_executable('simple_exec', args.simple_exec, args.simple_prefix, args.execute)
    adapters = adapters_for(args)
    args.private_exec = resolve_executable('simple_private_exec', args.private_exec, args.simple_prefix,
                                           args.execute and any(not adapters[n].shared or getattr(adapters[n], 'requires_private', False) for n in args.workflows))
    if args.execute and any(getattr(adapters[n], 'requires_private', False) for n in args.workflows):
        prefix = args.simple_prefix or os.environ.get('SIMPLE_PATH')
        if not prefix:
            raise ValueError('molecular calibration requires SIMPLE_PATH from the installed SIMPLE environment')
        worker = Path(prefix).absolute() / 'bin/simple_private_exec'
        if not worker.is_file() or not worker.samefile(args.private_exec):
            raise ValueError('molecular calibration worker must match SIMPLE_PATH/bin/simple_private_exec; load one consistent installation')
    if args.candidates is None:
        defaults = ['2x1','4x1','8x1','4x4','8x4'] if args.legacy and args.workflows == ['motion_correct'] else ['4x1','8x1','2x8','4x8','8x4']
        if not args.full and not args.legacy:
            defaults = ['4x1', '8x1', '4x2']
        args.candidates = [tuple(map(int, c.split('x'))) for c in defaults]
        args.candidates = [c for c in args.candidates if c[0] * c[1] <= len(cpus)] or [(1, 1)]
    else:
        args.candidates = [tuple(map(int, c.split('x'))) for c in args.candidates]
    args.threads = args.threads if args.threads is not None else [t for t in ([4,8,16] if args.full or args.legacy else [4,8]) if t <= len(cpus)] or [1]
    if (not cpus or args.repeats < 1 or args.timeout < 1 or args.solve3d_timeout < 1 or args.fixture_timeout < 1 or args.disk_gib < 1
            or len(set(args.workflows)) != len(args.workflows)
            or len(set(args.sizes)) != len(args.sizes) or any(s not in (4096,8192) for s in args.sizes)
            or len(set(args.candidates)) != len(args.candidates)
            or any(len(c) != 2 or min(c) < 1 or c[0]*c[1] > len(cpus) for c in args.candidates)
            or len(set(args.threads)) != len(args.threads) or any(t < 1 or t > len(cpus) for t in args.threads)
            or args.particles < 16 or args.box < 32 or args.box % 2 or args.extr_lim < 4 or args.nits_per_stage < 1):
        raise ValueError('invalid workload, repetitions, timeout or resource budget')
    if 'pick' in args.workflows and any(size <= 2*args.box for size in args.sizes):
        raise ValueError('picking micrographs must be larger than twice the particle box')
    return cpus


def calibrate(args):
    if args.slurm:
        if args.legacy:
            raise ValueError('use hardware_resource_planner.py calibrate for Slurm submission')
        return submit(args)
    if args.slurm_worker and (not os.environ.get('SLURM_JOB_ID') or
                              int(os.environ.get('SLURM_JOB_NUM_NODES', '1')) != 1):
        raise ValueError('Slurm worker requires a single-node allocation')
    cpus = configure(args)
    adapters = adapters_for(args)
    args._simulation_cache = {}
    plans = {name: {**adapters[name].metadata(args),
                   'candidate_timeout_seconds': args.solve3d_timeout if name == 'solve3D' else args.timeout}
             for name in args.workflows}
    if not args.execute:
        # Retain the legacy single-workflow plan's top-level fields.
        plan = next(iter(plans.values())) if args.legacy else {'workflows': plans}
        print(json.dumps({**plan, 'execution_enabled': False, 'cpus': cpus,
                          'repetitions': args.repeats, 'allocation': allocation()}, indent=2))
        return 0
    required = [args.simple_exec]
    if any(not adapters[n].shared or getattr(adapters[n], 'requires_private', False) for n in args.workflows):
        required.append(args.private_exec)
    for path in required:
        if not path.is_file() or not os.access(path, os.X_OK):
            raise ValueError(f'executable unavailable: {path}')
    if any(adapters[n].requires_numpy for n in args.workflows):
        subprocess.run([args.fixture_python, '-c', 'import numpy'], check=True, timeout=30)
    if args.legacy:
        if args.fixture_source != 'analytic':
            raise ValueError('use hardware_resource_planner.py calibrate for molecular workflow fixtures')
        result = run_workflow(adapters[args.workflows[0]], args, args.output.resolve(), cpus)
        print(json.dumps(result, indent=2))
        return 0
    output, root = args.output.resolve(), args.workdir.resolve()
    if output == root or output.is_relative_to(root / '_simulation') or any(output.is_relative_to(root / name) for name in args.workflows):
        raise ValueError('combined output must not overlap workflow evidence directories')
    if output.exists() or root.exists():
        raise ValueError('output profile and work directory must not already exist')
    root.mkdir(parents=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    collection = {'schema_version': 1, 'profile_type': 'simple_resource_profile_collection',
                  'auto_apply': False, 'hardware': machine_identity(), 'allocation': allocation(),
                  'calibration_preset': 'full' if args.full else 'standard',
                  'fixture_source': args.fixture_source, 'workflows': {}}
    failed = False
    for name in args.workflows:
        try:
            result = run_workflow(adapters[name], args, root / name, cpus)
        except Exception as error:
            # Isolate workflow failures, including malformed output, while letting
            # KeyboardInterrupt/SystemExit propagate to process-group cleanup.
            failed = True
            result = {'status': 'failed', 'reason': str(error), 'profile': None,
                      'evidence': str(root / name)}
        collection['workflows'][name] = result
        # Publish after each workflow so earlier evidence survives later failure.
        with tempfile.NamedTemporaryFile(mode='w', dir=output.parent,
                                         prefix=output.name + '.', suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            try:
                json.dump(collection, stream, indent=2, allow_nan=False)
                stream.write('\n')
                stream.close()
                temporary.replace(output)
            finally:
                temporary.unlink(missing_ok=True)
    print_results(collection, output)
    return 2 if failed else 0

