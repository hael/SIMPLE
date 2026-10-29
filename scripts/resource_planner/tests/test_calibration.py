"""Common calibration lifecycle checks using Python workers, never SIMPLE."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from resource_planner.calibration import add_options, calibrate, run_jobs, write_json, resolve_executable, Context, configure
from resource_planner.workflows.registry import ADAPTERS, MotionAdapter, CtfAdapter, Solve2DAdapter


class FakeAdapter:
    name = 'fake'
    shared = False
    requires_numpy = False
    rate = 'items_per_second'

    def cases(self, args): return [64]
    def candidates(self, args): return [(1, 1), (1, 2)]
    def metadata(self, args): return {'sizes': [64]}

    def prepare(self, ctx, size):
        source = ctx.root / 'fixture.txt'
        source.write_text('fixed scientific input')
        return {'folder': ctx.root, 'size': size, 'files': [source]}

    def jobs(self, ctx, fixture, directory, nthr, workers):
        jobs = []
        for index in range(workers):
            work = directory / str(index)
            work.mkdir()
            code = ("from pathlib import Path; "
                    "Path('product').write_text('correct'); "
                    "Path('memory_usage_fake.csv').write_text('peak_rss_bytes\\n1048576\\n')")
            jobs.append(([sys.executable, '-c', code], work))
        return jobs

    def validate(self, ctx, fixture, work):
        return {'passed': (work / 'product').read_text() == 'correct'}

    def export(self, ctx, plan, records, summary):
        assert len(records) == 8 and len(summary) == 2
        path = ctx.root / 'desktop_profile.json'
        write_json(path, {'commander': self.name, 'auto_apply': False})
        return path


class CalibrationChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        mount = patch('resource_planner.calibration.os.statvfs')
        mount.start().return_value.f_flag = 0
        self.addCleanup(mount.stop)

    def args(self, workflows='fake'):
        parser = argparse.ArgumentParser()
        add_options(parser)
        return parser.parse_args(['--output', str(self.root / 'profile.json'), '--workdir', str(self.root / 'evidence'),
                                  '--workflows', workflows, '--simple-exec', sys.executable,
                                  '--private-exec', sys.executable, '--threads', '1', '--candidates', '1x1', '--repeats', '3', '--execute'])

    def test_standard_full_and_explicit_workloads(self):
        parser = argparse.ArgumentParser()
        add_options(parser)
        with patch('resource_planner.calibration.cpu_pool', return_value=list(range(32))):
            standard = parser.parse_args([])
            configure(standard)
            self.assertEqual(standard.sizes, [4096])
            self.assertEqual(standard.particles, 100)
            self.assertEqual(standard.candidates, [(4,1), (8,1), (4,2)])
            self.assertEqual(standard.threads, [4,8])
            self.assertEqual((standard.extr_lim, standard.nits_per_stage, standard.repeats), (4,1,1))
            full = parser.parse_args(['--full'])
            configure(full)
            self.assertEqual(full.sizes, [4096,8192])
            self.assertEqual(full.particles, 2000)
            self.assertEqual(len(full.candidates), 5)
            self.assertEqual(full.threads, [1,2,4,8,16,32])
            custom = parser.parse_args(['--full','--sizes','4096','--particles','800','--threads','2',
                                        '--candidates','2x1','--timeout','90'])
            configure(custom)
            self.assertEqual((custom.sizes,custom.particles,custom.threads,custom.candidates,custom.timeout),
                             ([4096],800,[2],[(2,1)],90))
        with patch('resource_planner.calibration.cpu_pool', return_value=[0]):
            small = parser.parse_args([])
            configure(small)
            self.assertEqual(small.candidates, [(1,1)])
            self.assertEqual(small.threads, [1])

    def test_environment_executable_discovery(self):
        bindir = self.root / 'bin'
        bindir.mkdir()
        executable = bindir / 'simple_exec'
        executable.symlink_to(sys.executable)
        with patch.dict(os.environ, {'PATH': str(bindir), 'SIMPLE_PATH': '/unused'}):
            self.assertEqual(resolve_executable('simple_exec', required=True), executable)
            self.assertEqual(resolve_executable('simple_exec', Path(sys.executable), required=True),
                             Path(sys.executable).absolute())
        with patch.dict(os.environ, {'PATH': '', 'SIMPLE_PATH': str(self.root)}):
            self.assertEqual(resolve_executable('simple_exec', required=True), executable)
            with self.assertRaisesRegex(ValueError, 'load the SIMPLE installation environment'):
                resolve_executable('simple_private_exec', required=True)
        with patch.dict(os.environ, {'PATH': '', 'SIMPLE_PATH': ''}):
            self.assertEqual(resolve_executable('simple_exec'), Path('simple_exec'))
            with self.assertRaises(ValueError):
                resolve_executable('simple_exec', required=True)

    def test_context_preserves_installation_environment(self):
        with patch.dict(os.environ, {'SIMPLE_PATH': '/installed/simple'}):
            args = self.args('solve2D')
            context = Context(args, self.root, [0])
            self.assertEqual(context.env['SIMPLE_PATH'], '/installed/simple')
            self.assertEqual(context.env['SIMPLE_SEED'], '20261015')
            self.assertEqual(context.env['OMP_NUM_THREADS'], '1')
            args.simple_prefix = self.root
            self.assertEqual(Context(args, self.root, [0]).env['SIMPLE_PATH'], str(self.root))

    def test_unified_execution_and_collection(self):
        with patch.dict(ADAPTERS, {'fake': FakeAdapter()}), patch('resource_planner.calibration.cpu_pool', return_value=sorted(os.sched_getaffinity(0))[:2]):
            args = self.args()
            self.assertEqual(calibrate(args), 0)
        collection = json.loads(args.output.read_text())
        self.assertEqual(collection['profile_type'], 'simple_resource_profile_collection')
        self.assertFalse(collection['auto_apply'])
        self.assertEqual(collection['workflows']['fake']['status'], 'provisional')
        profile = json.loads((args.workdir / 'fake/desktop_profile.json').read_text())
        reproducibility = profile['reproducibility']
        self.assertEqual(reproducibility['seed'], 20261015)
        self.assertEqual(len(reproducibility['fixtures']['64']), 1)
        self.assertEqual(reproducibility, json.loads((args.workdir / 'fake/plan.json').read_text())['reproducibility'])
        rows = [json.loads(s) for s in (args.workdir / 'fake/results.jsonl').read_text().splitlines()]
        self.assertEqual(len(rows), 8)
        self.assertEqual(sum(r['warmup'] for r in rows), 2)
        self.assertTrue(all(r['passed'] for r in rows))

    def test_changed_fixture_cannot_export_recommendation(self):
        class Mutating(FakeAdapter):
            def validate(self, ctx, fixture, work):
                fixture['files'][0].write_text('changed input')
                return super().validate(ctx, fixture, work)
        with patch.dict(ADAPTERS, {'fake': Mutating()}):
            args = self.args()
            self.assertNotEqual(calibrate(args), 0)
        self.assertFalse((args.workdir / 'fake/desktop_profile.json').exists())
        result = json.loads(args.output.read_text())['workflows']['fake']
        self.assertEqual(len(result['measurements']), 8)

    def test_unqualified_keeps_evidence_without_recommendation(self):
        class Unqualified(FakeAdapter):
            def export(self, *args): return None
        with patch.dict(ADAPTERS, {'fake': Unqualified()}):
            args = self.args()
            self.assertEqual(calibrate(args), 0)
        result = json.loads(args.output.read_text())['workflows']['fake']
        self.assertEqual(result['status'], 'unqualified')
        self.assertIsNone(result['profile'])
        self.assertTrue((Path(result['evidence']) / 'summary.json').exists())

    def test_failure_written_to_collection(self):
        class Failed(FakeAdapter):
            def prepare(self, *args): raise ValueError('fixture failed')
        with patch.dict(ADAPTERS, {'fake': Failed()}):
            args = self.args()
            self.assertEqual(calibrate(args), 2)
        result = json.loads(args.output.read_text())['workflows']['fake']
        self.assertEqual(result['status'], 'failed')
        self.assertIn('fixture failed', result['reason'])
        self.assertIsNone(result['profile'])

    def test_failed_workflow_does_not_hide_later_results(self):
        class Failed(FakeAdapter):
            name = 'broken'
            def prepare(self, *args): raise ValueError('fixture failed')
        with patch.dict(ADAPTERS, {'fake': FakeAdapter(), 'broken': Failed()}):
            args = self.args()
            args.workflows = ['broken', 'fake']
            self.assertEqual(calibrate(args), 2)
        workflows = json.loads(args.output.read_text())['workflows']
        self.assertEqual(workflows['broken']['status'], 'failed')
        self.assertEqual(workflows['fake']['status'], 'provisional')

    def test_reported_failure_returns_nonzero_without_exception(self):
        with patch.dict(ADAPTERS, {'fake': FakeAdapter()}):
            args = self.args()
            with patch('resource_planner.calibration.run_workflow', return_value={
                    'status': 'failed', 'reason': 'placement failed', 'profile': None, 'measurements': []}):
                self.assertEqual(calibrate(args), 2)

    def test_partial_node_failure_returns_nonzero(self):
        with patch.dict(ADAPTERS, {'fake': FakeAdapter()}):
            args = self.args()
            with patch('resource_planner.calibration.run_workflow', return_value={
                    'status': 'initial', 'profile': None, 'measurements': [],
                    'node_comparison': {'2': {'status': 'failed', 'reason': 'placement failed'}}}):
                self.assertEqual(calibrate(args), 2)

    def test_solve_adapter_prepares_each_project_separately(self):
        prepared = []
        class Context:
            root = self.root
            args = argparse.Namespace(simple_exec=Path('/simple_exec'), particles=2000, box=128,
                                      extr_lim=20, nits_per_stage=5, threads=[4])
            def project(self, directory, *args, **kwargs):
                prepared.append(directory)
        adapter = Solve2DAdapter()
        for name in ('first', 'second'):
            work = self.root / name
            work.mkdir()
            jobs = adapter.jobs(Context(), {}, work, 4, 1)
            self.assertEqual(jobs[0][0][1], 'prg=solve2D')
            self.assertIn(f'projfile={work}/classes/classes.simple', jobs[0][0])
            self.assertFalse(any(a.startswith('nparts=') for a in jobs[0][0]))
        self.assertEqual(len(set(prepared)), 2)

    def test_timeout_reaps_process(self):
        code = "import os,time; from pathlib import Path; Path('pid').write_text(str(os.getpid())); time.sleep(60)"
        with self.assertRaises(TimeoutError):
            run_jobs([([sys.executable, '-c', code], self.root)], os.environ.copy(), [],
                     0.3, self.root, 1024**2)
        pid = int((self.root / 'pid').read_text())
        with self.assertRaises(ProcessLookupError): os.kill(pid, 0)

    def test_worker_failure_and_disk_limit(self):
        with self.assertRaisesRegex(RuntimeError, 'process failed'):
            run_jobs([([sys.executable, '-c', 'raise SystemExit(7)'], self.root)],
                     os.environ.copy(), [], 5, self.root, 1024**2)
        with self.assertRaisesRegex(RuntimeError, 'disk budget'):
            run_jobs([([sys.executable, '-c', 'import time; time.sleep(60)'], self.root)],
                     os.environ.copy(), [], 5, self.root, 1)

    def test_cli_dry_run_and_unsupported_workflow(self):
        cli = Path(__file__).resolve().parents[2] / 'hardware_resource_planner.py'
        command = [sys.executable, str(cli), 'calibrate', '--output', str(self.root / 'profile.json'),
                   '--workdir', str(self.root / 'evidence')]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(set(json.loads(result.stdout)['workflows']), {'motion_correct', 'ctf_estimate', 'pick', 'solve2D', 'solve3D'})
        self.assertEqual(list(self.root.iterdir()), [])
        bad = subprocess.run(command + ['--workflows', 'not_a_workflow'], capture_output=True, text=True)
        self.assertNotEqual(bad.returncode, 0)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_analytic_workflows_dry_run(self):
        scripts = Path(__file__).resolve().parents[2]
        for name in ('motion_correct', 'ctf_estimate', 'solve2D'):
            result = subprocess.run([sys.executable, str(scripts / 'hardware_resource_planner.py'), 'calibrate',
                                     '--fixture-source', 'analytic', '--workflows', name,
                                     '--output', str(self.root / name)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(json.loads(result.stdout)['execution_enabled'])
            self.assertFalse((self.root / name).exists())

    def test_image_adapter_commands_preserve_science(self):
        class Context:
            args = argparse.Namespace(private_exec=Path('/simple_private_exec'))
        for adapter, key in ((MotionAdapter(), 'mcpatch=yes'), (CtfAdapter(), 'ctfpatch=yes')):
            directory = self.root / adapter.name
            directory.mkdir()
            jobs = adapter.jobs(Context(), {'project': Path('/input.simple')}, directory, 4, 2)
            self.assertEqual(len(jobs), 2)
            for argv, cwd in jobs:
                self.assertEqual(argv[1], 'prg=' + adapter.name)
                self.assertIn(key, argv)
                self.assertIn('nthr=4', argv)
                self.assertIn('nparts=1', argv)
                self.assertTrue(cwd.is_dir())


if __name__ == '__main__': unittest.main()
