"""Molecular preparation contracts; native commands are replaced by recording fakes."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from resource_planner.calibration import add_options, adapters_for, apply_preset, run_jobs
from resource_planner.workflows.molecular import Simulation, WorkflowAdapter, SEED, science


class WorkflowChecks(unittest.TestCase):
    def test_solve3d_timeout_defaults_and_explicit_overrides(self):
        parser = argparse.ArgumentParser()
        add_options(parser)
        for flags, expected in [([], (180, 1800)),
                                (['--full'], (600, 1800)),
                                (['--timeout', '700'], (700, 700)),
                                (['--timeout', '700', '--solve3d-timeout', '1200'], (700, 1200))]:
            args = parser.parse_args(flags)
            apply_preset(args)
            self.assertEqual((args.timeout, args.solve3d_timeout), expected)

    def test_timeout_reaps_detached_local_worker(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            child = "import os,time; from pathlib import Path; Path('worker.pid').write_text(str(os.getpid())); time.sleep(30)"
            parent = 'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c",%r],start_new_session=True); time.sleep(30)' % child
            env = dict(os.environ, SIMPLE_CALIBRATION_SUPERVISE='1')
            with self.assertRaises(TimeoutError):
                run_jobs([([sys.executable,'-c',parent],root)], env, [], 0.5, root, 1024**2)
            pid = int((root/'worker.pid').read_text())
            with self.assertRaises(ProcessLookupError):
                os.kill(pid,0)

    def args(self):
        parser = argparse.ArgumentParser()
        add_options(parser)
        args = parser.parse_args([])
        args.particles, args.box, args.sizes, args.threads = 500, 128, [4096], [2,4]
        args.extr_lim, args.nits_per_stage = 4, 1
        args.simple_exec = Path('/installed/simple_exec')
        args._simulation_cache = {}
        return args

    def test_default_uses_suite_recipe_and_analytic_is_explicit(self):
        args = self.args()
        self.assertEqual(args.fixture_source, 'simulated-workflow')
        adapters = adapters_for(args)
        for name, adapter in adapters.items():
            self.assertIsInstance(adapter, WorkflowAdapter)
            self.assertEqual(adapter.seed, 20260923)
            self.assertEqual(adapter.candidates(args), [(2,1),(4,1)])
        self.assertEqual(science('solve3D', args)['pgrp'], 'c3')
        self.assertNotIn('vol1', science('solve3D', args))
        self.assertNotIn('cls_init', science('solve2D', args))
        args.fixture_source = 'analytic'
        self.assertNotIsInstance(adapters_for(args)['solve2D'], WorkflowAdapter)

    def test_serial_dependencies_once_and_candidate_project_isolation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            class Context:
                args = self.args()
                def __init__(self):
                    self.root = root/'solve2D'
                    self.root.mkdir()
                    self.calls = []
                def command(self, argv, work, label, timeout=60):
                    self.calls.append(argv)
                    project = next(Path(a.split('=',1)[1]) for a in argv if a.startswith('projfile='))
                    project.write_text(project.read_text()+'\n'+argv[1])
            ctx = Context()
            def generate(sim):
                sim.root.mkdir(parents=True)
                project = sim.root/'initial.simple'
                project.write_text('immutable import')
                sim.snapshots['solve2D'] = project
                sim.truth = {'seed': SEED}
            with patch.object(Simulation, 'generate', generate):
                adapter = WorkflowAdapter('solve2D')
                fixture = adapter.prepare(ctx,4096)
                self.assertEqual([a[1] for a in ctx.calls],
                    [])
                self.assertTrue(all('nthr=1' in a for a in ctx.calls))
                original = fixture['project'].read_bytes()
                for i in range(2):
                    work = root/str(i); work.mkdir()
                    argv,_ = adapter.jobs(ctx,fixture,work,4,1)[0]
                    self.assertIn('nthr=4',argv)
                    copy = work/'calibration.simple'
                    self.assertEqual(copy.read_bytes(),original)
                    copy.write_text('candidate changed this project')
                self.assertEqual(fixture['project'].read_bytes(),original)
                simulation = ctx.args._simulation_cache[('particles',4096)]
                self.assertEqual(simulation.snapshot('solve2D'),fixture['project'])
                simulation.snapshot('solve3D')
                simulation.snapshot('solve3D')
                self.assertEqual([a[1] for a in ctx.calls], ['prg=solve2D'])
                self.assertEqual((simulation.root/'initial.simple').read_text(),'immutable import')

    def test_particle_generation_imports_direct_stack_without_image_commands(self):
        with tempfile.TemporaryDirectory() as directory:
            class Context:
                args = self.args()
                root = Path(directory) / 'solve2D'
                calls = []
                def command(self, argv, work, label, timeout=60):
                    self.calls.append(argv)
            ctx = Context()
            simulation = Simulation(ctx, 128, particles=True)
            simulation.generate()
            native = [a for a in ctx.calls if a[1].startswith('prg=')]
            self.assertEqual([a[1] for a in native], ['prg=new_project', 'prg=import_particles'])
            self.assertIn('ctf=no', native[-1])
            self.assertIn('stk='+str(simulation.root/'particles.mrcs'), native[-1])
            self.assertEqual(len((simulation.root/'halves.txt').read_text().splitlines()), 500)
            self.assertEqual(simulation.truth['particle_stack'], str(simulation.root/'particles.mrcs'))
            self.assertFalse(any(a == 'movie' for argv in ctx.calls for a in argv))

    def test_generation_uses_embedded_model_and_unique_movies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            class Context:
                args = self.args()
                def __init__(self):
                    self.root = root/'motion_correct'
                    self.calls = []
                def command(self, argv, work, label, timeout=60):
                    self.calls.append(argv)
                    if label == 'molecule':
                        (work/'truth.mrc').write_bytes(b'fake volume')
                    if label == 'simulate':
                        (work/'simulate_movie.mrc').write_bytes(b'fake movie')
            ctx = Context(); simulation = Simulation(ctx,4096)
            simulation.generate()
            native = [argv for argv in ctx.calls if len(argv)>1 and argv[1].startswith('prg=')]
            self.assertTrue(all('nthr=1' in argv for argv in native))
            self.assertIn('--model',ctx.calls[0])
            self.assertIn('6vxx',ctx.calls[0])
            self.assertFalse(any(any(a.startswith('prg=simulat') or a=='prg=reproject' for a in argv) for argv in ctx.calls))
            self.assertEqual(sum(m['particles'] for m in simulation.truth['movies']),500)
            movies = [m['movie'] for m in simulation.truth['movies']]
            self.assertEqual(len(movies),len(set(Path(p).name for p in movies)))
            self.assertTrue(all(Path(p).exists() for p in movies))

    def test_numpy_subset_replay_and_noise_volume_rejection(self):
        # A separate NumPy interpreter is also supported by the administrator CLI.
        python = os.environ.get('SIMPLE_FIXTURE_PYTHON',sys.executable)
        probe = subprocess.run([python,'-c','import numpy'],capture_output=True)
        if probe.returncode:
            self.skipTest('NumPy fixture interpreter unavailable')
        helper_dir = str(Path(__file__).resolve().parents[1] / "fixtures")
        code = r'''
import sys, tempfile, json
from pathlib import Path
import numpy as np
sys.path.insert(0,sys.argv[1])
from validation import write, subset, read_mrc, validate
from molecular import molecule, movie
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory)
    rng=np.random.RandomState(17)
    source=root/'source.mrcs'
    write(source,rng.normal(size=(8,32,32)),1.875)
    subset(source,root/'a.mrcs',21,20260923)
    subset(source,root/'b.mrcs',21,20260923)
    assert (root/'a.mrcs').read_bytes()==(root/'b.mrcs').read_bytes()
    subset(source,root/'c.mrcs',21,20260924)
    assert (root/'a.mrcs').read_bytes()!=(root/'c.mrcs').read_bytes()
    z,y,x=np.mgrid[:32,:32,:32]-16
    truth=np.exp(-(x*x+y*y+z*z)/25.)
    write(root/'truth.mrc',truth,7.5)
    write(root/'particles.mrcs',truth[:1],7.5)
    for tag in ('direct','mirror'):
        write(root/(tag+'_docked.mrc'),truth,7.5)
        (root/(tag+'_dock.txt')).write_text('cc=0.9\n')
    (root/'assignments.log').write_text('state=1 e1=0 e2=0 e3=0\n')
    result=validate('solve3D',root,{'volume':str(root/'truth.mrc'),'smpd':7.5,'particle_stack':str(root/'particles.mrcs')},root)
    assert result['quality_qualified'], result
    for tag in ('direct','mirror'):
        (root/(tag+'_docked.mrc')).unlink()
        write(root/(tag+'_docked.mrc'),rng.normal(size=truth.shape),7.5)
    result=validate('solve3D',root,{'volume':str(root/'truth.mrc'),'smpd':7.5,'particle_stack':str(root/'particles.mrcs')},root)
    assert not result['quality_qualified'], result
    for model in ('6vxx','1jyx'):
        folders=[root/(model+str(i)) for i in range(2)]
        for folder in folders:
            folder.mkdir()
            molecule(folder,model,32,7.5)
        for filename in ('truth.mrc','projections.mrcs','generator.json'):
            assert (folders[0]/filename).read_bytes()==(folders[1]/filename).read_bytes()
        for folder in folders:
            sample=folder/'projections.mrcs'
            subset(sample,folder/'subset.mrcs',4,20260923)
            sample.unlink()
            (folder/'subset.mrcs').rename(sample)
            movie(folder,256,20260923)
        for filename in ('simulate_movie.mrc','optimal_movie_average.mrc','simulate_movie.box','simulate_movie_params.txt'):
            assert (folders[0]/filename).read_bytes()==(folders[1]/filename).read_bytes()
        frames=read_mrc(folders[0]/'simulate_movie.mrc')
        aligned=np.mean([np.roll(frame,-shift,axis=1) for frame,shift in zip(frames,[0,0,0,0,0,1,2,3])],axis=0)
        np.testing.assert_allclose(aligned,read_mrc(folders[0]/'optimal_movie_average.mrc')[0],atol=1e-6)
'''
        result = subprocess.run([python,'-c',code,helper_dir],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)


if __name__=='__main__':
    unittest.main()
