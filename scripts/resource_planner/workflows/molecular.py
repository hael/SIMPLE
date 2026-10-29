"""Molecular calibration following the deterministic simulated_workflow suite.

Preparation runs serially once. Candidates copy project snapshots and never
feed their results into the next commander's calibration.
"""
import json
from pathlib import Path
import re
import shutil

from resource_planner.workflows.common import export_shared

SEED = 20260923
STAGES = ('motion_correct', 'ctf_estimate', 'pick', 'solve2D', 'solve3D')
HELPER = Path(__file__).resolve().parents[1] / 'fixtures' / 'validation.py'
GENERATOR = Path(__file__).resolve().parents[1] / 'fixtures' / 'molecular.py'


def sampling(args):
    # Keep the 180 A molecular support inside a box with a generous margin.
    return max(1.3, 240.0 / args.box)


def science(name, args):
    common = {'wfloat16': 'no'}
    if name == 'motion_correct':
        return {**common, 'nparts': 1, 'algorithm': 'patch', 'mcpatch': 'yes', 'downscale': 'no'}
    if name == 'ctf_estimate':
        return {**common, 'nparts': 1, 'ctfpatch': 'yes', 'pspecsz': 512,
                'dfmin': 0.5, 'dfmax': 5, 'hp': 30, 'lp': 5}
    if name == 'pick':
        return {**common, 'nparts': 1, 'picker': 'segdiam', 'pcontrast': 'black', 'moldiam_max': 180}
    if name == 'solve2D':
        return {**common, 'mskdiam': 180, 'ncls': 8, 'extr_lim': args.extr_lim,
                'nits_per_stage': args.nits_per_stage}
    pgrp = 'c3' if args.molecule == '6vxx' else 'd2'
    settings = {**common, 'mskdiam': 180, 'pgrp': pgrp, 'pgrp_start': pgrp}
    if not getattr(args, 'full', False):
        settings.update(nstages=4, balance='none', nsample=args.particles,
                        filt_mode='none', automsk='no', force_lp_range='yes', lpstart=30, lpstop=15)
    return settings


def command(ctx, name, project, nthr):
    return [str(ctx.args.simple_exec), 'prg=' + name, 'projfile=' + str(project),
            'mkdir=no', 'qsys_name=local', 'nthr=' + str(nthr),
            *[f'{key}={value}' for key, value in science(name, ctx.args).items()],
            'memreport=yes', 'memreport_interval=1']


def copy_snapshot(source, destination):
    """Copy mutable project state; native project data paths remain absolute."""
    target = destination / 'calibration.simple'
    shutil.copy2(source, target)
    return target


class Simulation:
    def __init__(self, ctx, size, particles=False):
        self.ctx, self.size = ctx, size
        self.particles = particles
        key = 'particles_box%d' % size if particles else str(size)
        self.root = ctx.root.parent / '_simulation' / key
        self.snapshots, self.outputs = {}, {}
        self.error = None

    def native(self, program, directory, label, **parameters):
        argv = [str(self.ctx.args.simple_exec), 'prg=' + program, 'mkdir=no', 'nthr=1',
                *[f'{k}={v}' for k, v in parameters.items()]]
        self.ctx.command(argv, directory, label, timeout=self.ctx.args.fixture_timeout)

    def generate(self):
        args = self.ctx.args
        self.root.mkdir(parents=True)
        smpd = sampling(args)
        self.ctx.command([args.fixture_python,str(GENERATOR),'molecule',str(self.root),
            '--model',args.molecule,'--box',str(args.box),'--smpd',str(smpd)],
            self.root,'molecule',args.fixture_timeout)
        if self.particles:
            stack = self.root / 'particles.mrcs'
            self.ctx.command([args.fixture_python, str(HELPER), 'subset',
                '--input', str(self.root/'projections.mrcs'), '--output', str(stack),
                '--count', str(args.particles), '--seed', str(SEED)],
                self.root, 'particles', args.fixture_timeout)
            halves = self.root / 'halves.txt'
            halves.write_text(''.join('state=1 eo=%d\n' % (i % 2) for i in range(args.particles)))
            self.truth = {'seed': SEED, 'molecule': args.molecule, 'smpd': smpd,
                          'box': args.box, 'particles_requested': args.particles,
                          'particle_stack': str(stack), 'volume': str(self.root/'truth.mrc'),
                          'particle_model': 'clean molecular projections; ctf=no'}
            (self.root/'truth.json').write_text(json.dumps(self.truth, indent=2)+'\n')
            self.native('new_project', self.root, 'new_project', projname='calibration', qsys_name='local')
            project = self.root / 'calibration/calibration.simple'
            self.native('import_particles', project.parent, 'import', projfile=project,
                        stk=stack, smpd=smpd, ctf='no', kv=300, cs=2.7, fraca=0.1, oritab=halves, oritype='ptcl2D')
            self.snapshots['solve2D'] = project
            if not getattr(args, 'full', False):
                self.snapshots['solve3D'] = project
            return
        # Limit occupied area to 10%, keeping random non-overlapping placement cheap.
        per_movie = max(1, int(self.size**2 * 0.10 / args.box**2))
        counts = [min(per_movie, args.particles-i) for i in range(0, args.particles, per_movie)]
        movies = []
        for index, count in enumerate(counts):
            folder = self.root / ('movie_%03d' % index)
            folder.mkdir()
            self.ctx.command([args.fixture_python, str(HELPER), 'subset',
                '--input', str(self.root/'projections.mrcs'), '--output', str(folder/'projections.mrcs'),
                '--count', str(count), '--seed', str(SEED+index)], folder, 'subset', args.fixture_timeout)
            self.ctx.command([args.fixture_python,str(GENERATOR),'movie',str(folder),
                '--size',str(self.size),'--seed',str(SEED+index)],folder,'simulate',args.fixture_timeout)
            # Distinct basenames are required when a master processes multiple movies.
            movie = folder / ('movie_%03d.mrc' % index)
            (folder/'simulate_movie.mrc').rename(movie)
            movies.append({'movie': str(movie), 'optimal': str(folder/'optimal_movie_average.mrc'),
                           'parameters': str(folder/'simulate_movie_params.txt'),
                           'boxes': str(folder/'simulate_movie.box'), 'particles': count})
        self.truth = {'seed': SEED, 'molecule': args.molecule, 'smpd': smpd, 'box': args.box,
                      'size': self.size, 'frames': 8, 'particles_requested': args.particles,
                      'movies': movies, 'volume': str(self.root/'truth.mrc'),
                      'generator': 'external Python/NumPy; bundled coordinate data; no native simulator',
                      'reference_workflow': 'simple_commanders_test_highlevel_workflow'}
        (self.root/'truth.json').write_text(json.dumps(self.truth, indent=2)+'\n')
        inputs = self.root / 'movies.txt'
        inputs.write_text(''.join(row['movie']+'\n' for row in movies))
        self.native('new_project', self.root, 'new_project', projname='calibration', qsys_name='local')
        project = self.root / 'calibration/calibration.simple'
        self.native('import_movies', project.parent, 'import', projfile=project,
                    filetab=inputs, smpd=smpd, kv=300, cs=2.7, fraca=0.1, ctf='yes')
        self.snapshots['motion_correct'] = project

    def snapshot(self, name):
        if self.error:
            raise RuntimeError(self.error)
        try:
            if not self.snapshots:
                self.generate()
            stages = STAGES[3:] if self.particles else STAGES[:3]
            for previous, following in zip(stages, stages[1:]):
                if name in self.snapshots:
                    break
                if following in self.snapshots:
                    continue
                work = self.root / ('prepare_' + previous)
                work.mkdir()
                project = copy_snapshot(self.snapshots[previous], work)
                print('Preparing deterministic input: ' + previous + ' (one thread)', flush=True)
                self.ctx.command(command(self.ctx, previous, project, 1), work, 'prepare',
                                 timeout=self.ctx.args.fixture_timeout)
                self.outputs[previous] = work
                self.snapshots[following] = project
            return self.snapshots[name]
        except Exception as error:
            self.error = 'Molecular fixture preparation failed: ' + str(error)
            raise

    def files(self, project):
        # Project and all generated data it can reference, but no logs/telemetry.
        suffixes = {'.mrc', '.mrcs', '.simple', '.box', '.txt', '.json', '.star', '.bin'}
        return sorted(p for p in self.root.rglob('*') if p.is_file() and p.suffix in suffixes
                      and not p.name.endswith('_command.json') and p.name != 'metadata.json')


class WorkflowAdapter:
    shared = True
    supervise = True
    requires_private = True
    aggregate_memory = True
    requires_numpy = True
    rate = None
    seed = SEED

    def __init__(self, name):
        self.name = name
        self.qualification = ('provisional molecular simulated_workflow; single-node master; '
                              'real data and multi-node scaling unqualified')
        if name == 'solve2D':
            self.qualification += '; structural class-output checks only, classification accuracy unqualified'

    def cases(self, args):
        return args.sizes if self.name in STAGES[:3] else [args.box]

    def candidates(self, args):
        if getattr(args, '_nodes', None):
            return [(t, 1) for t in args.threads]
        if args.full and self.name in ('motion_correct', 'ctf_estimate'):
            import math
            per_movie = max(1, int(max(args.sizes)**2 * 0.10 / args.box**2))
            movies = math.ceil(args.particles / per_movie)
            return [(t, w) for t in args.threads for w in (1, 2, 4)
                    if w <= movies and t*w <= getattr(args, '_cpu_budget', max(args.threads))]
        return [(t, 1) for t in args.threads]

    def metadata(self, args):
        return {'threads': args.threads, 'tested_candidates': self.candidates(args), 'workload': {'fixture_source': 'simulated_workflow',
            'fixture_generator': 'external_python_equal_weight_atoms', 'fixture_version': 3,
            'molecule': args.molecule, 'seed': SEED, 'box': args.box, 'smpd': sampling(args),
            'particles_requested': args.particles, 'frames': 8, 'sizes': self.cases(args),
            'particle_input': 'direct clean molecular projections, ctf=no' if self.name in STAGES[3:] else None,
            'initialization': ('known_volume_single_stage' if self.name == 'solve3D' and not args.full
                               else 'de_novo' if self.name in STAGES[3:] else 'previous_serial_stage'),
            'execution_scope': 'single-node master; full image workflows sweep local worker/partition counts',
            **science(self.name, args)}}

    def prepare(self, ctx, size):
        cache = ctx.args._simulation_cache
        particles = self.name in STAGES[3:]
        key = ('particles', size) if particles else size
        if key not in cache:
            cache[key] = Simulation(ctx, size, particles=particles)
        simulation = cache[key]
        project = simulation.snapshot(self.name)
        folder = ctx.root / str(size)
        folder.mkdir()
        return {'folder': folder, 'size': size, 'simulation': simulation,
                'project': project, 'files': simulation.files(project), 'truth': simulation.truth}

    def jobs(self, ctx, fixture, directory, nthr, workers):
        project = copy_snapshot(fixture['project'], directory)
        argv = command(ctx, self.name, project, nthr)
        if workers > 1:
            argv = [a for a in argv if not a.startswith('nparts=')]
            argv += [f'nparts={workers}', f'ncunits={workers}']
        if self.name == 'solve3D' and not ctx.args.full:
            argv += ['vol1=' + str(fixture['simulation'].root / 'truth.mrc'), f'nthr_ini3D={nthr}']
        if self.name == 'ctf_estimate':
            fixture['consumables'] = {}
            for path in fixture['simulation'].root.rglob('*_forctf.mrc'):
                backup = directory / (path.name + '.backup')
                shutil.copy2(path, backup)
                fixture['consumables'][path] = backup
        return [(argv, directory)]

    def validate(self, ctx, fixture, work):
        simulation = fixture['simulation']
        for path, backup in fixture.get('consumables', {}).items():
            if not path.exists():
                shutil.copy2(backup, path)
            backup.unlink()
        # Validation and project inspection are deliberately outside the timed command.
        oritype = 'ptcl3D' if self.name == 'solve3D' else 'ptcl2D' if self.name == 'solve2D' else 'mic'
        ctx.command([str(ctx.args.simple_exec), 'prg=print_project_field',
                     'projfile='+str(work/'calibration.simple'), 'oritype='+oritype,
                     'mkdir=no'], work, 'assignments', timeout=ctx.args.fixture_timeout)
        if self.name == 'solve3D':
            volumes = sorted(p for p in work.glob('recvol_state*.mrc') if re.fullmatch(r'recvol_state\d+\.mrc', p.name))
            if len(volumes) != 1:
                raise ValueError('expected one reconstructed molecular state')
            ctx.command([ctx.args.fixture_python, str(HELPER), 'mirror', '--input', str(volumes[0]),
                         '--output', str(work/'mirror.mrc')], work, 'mirror', ctx.args.fixture_timeout)
            for tag, source in [('direct',volumes[0]), ('mirror',work/'mirror.mrc')]:
                simulation.native('dock_volpair', work, 'dock_'+tag, vol1=simulation.root/'truth.mrc',
                    vol2=source, outvol=work/(tag+'_docked.mrc'), outfile=work/(tag+'_dock.txt'),
                    smpd=sampling(ctx.args), mskdiam=180, hp=100, lp=20)
        result = work / 'workflow_validation.json'
        ctx.command([ctx.args.fixture_python, str(HELPER), 'validate', '--workflow', self.name,
                     '--input', str(work), '--truth', str(simulation.root/'truth.json'),
                     '--preparation', str(simulation.root), '--output', str(result)],
                    work, 'validate', ctx.args.fixture_timeout)
        return json.loads(result.read_text())

    def export(self, ctx, plan, records, summary):
        counts = {r['validation'].get('particles_assigned') for r in records}
        counts.discard(None)
        if self.name in ('solve2D', 'solve3D') and len(counts) != 1:
            raise ValueError('particle counts differ across calibration candidates')
        if counts:
            plan = {**plan, 'workload': {**plan['workload'], 'observed_particle_counts': sorted(counts)}}
        return export_shared(self, ctx, plan, records)


ADAPTERS = {name: WorkflowAdapter(name) for name in STAGES}
