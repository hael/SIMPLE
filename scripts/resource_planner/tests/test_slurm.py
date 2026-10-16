"""Slurm submission checks without contacting a scheduler."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, Mock

from resource_planner.calibration import add_options, calibrate
from hardware_resource_planner import main, ensure_python
from resource_planner.hardware import cpu_pool


class SlurmChecks(unittest.TestCase):
    def test_no_arguments_executes_standard_and_explicit_preview_stays_dry(self):
        with patch('resource_planner.cli.calibrate', return_value=0) as run:
            self.assertEqual(main([]), 0)
            args = run.call_args.args[0]
            self.assertTrue(args.execute)
            self.assertEqual(len(args.workflows), 5)
            self.assertEqual(args.output, Path('output/resource_calibration/standard/profile.json'))
            self.assertEqual(main(['calibrate']), 0)
            self.assertFalse(run.call_args.args[0].execute)

    def test_old_python_selects_installed_modern_interpreter(self):
        with patch('hardware_resource_planner.sys.version_info', (3, 6)), \
                patch('hardware_resource_planner.shutil.which', side_effect=lambda name: '/usr/bin/python3.11' if name == 'python3.11' else None), \
                patch('hardware_resource_planner.os.execv') as execute:
            ensure_python()
            self.assertEqual(execute.call_args.args[0], '/usr/bin/python3.11')
        with patch('hardware_resource_planner.sys.version_info', (3, 6)), \
                patch('hardware_resource_planner.shutil.which', return_value=None):
            with self.assertRaisesRegex(SystemExit, 'Python 3.9'):
                ensure_python()

    def test_standard_defaults(self):
        parser = argparse.ArgumentParser()
        add_options(parser)
        args = parser.parse_args([])
        self.assertFalse(args.execute)
        self.assertTrue(args.discard_images)
        self.assertEqual(args.output, Path('output/resource_calibration/standard/profile.json'))
        self.assertEqual(args.workdir, Path('output/resource_calibration/standard/measurements'))
        self.assertEqual(args.workflows, ['motion_correct', 'ctf_estimate', 'pick', 'solve2D', 'solve3D'])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser.parse_args(['--workflows', 'extract'])

    def test_preview_and_submission(self):
        with tempfile.TemporaryDirectory(prefix='calibration space ') as directory:
            parser = argparse.ArgumentParser()
            add_options(parser)
            args = parser.parse_args(['--slurm', 'csbdevel', '--slurm-cpus', '8',
                                      '--solve3d-timeout', '1200',
                                      '--output', directory + '/profile.json',
                                      '--workdir', directory + '/evidence', '--workflows', 'ctf_estimate'])
            with patch.dict(os.environ, {}, clear=True), patch('resource_planner.slurm.subprocess.run') as run:
                with contextlib.redirect_stdout(io.StringIO()) as stream:
                    self.assertEqual(calibrate(args), 0)
                run.assert_not_called()
                preview = json.loads(stream.getvalue())
                self.assertIn('--fixture-source simulated-workflow', preview['script'])
                self.assertIn('--molecule 6vxx', preview['script'])
                self.assertIn('--fixture-timeout 900', preview['script'])
                self.assertIn('--solve3d-timeout 1200', preview['script'])
                self.assertIn('--partition=csbdevel', preview['submission'])
                self.assertIn('--cpu-bind=cores', preview['script'])
                self.assertIn('--slurm-worker', preview['script'])
                self.assertNotIn('--slurm csbdevel', preview['script'])
                self.assertFalse(args.output.exists())
                args.execute = True
                run.return_value = Mock(stdout='12345\n')
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(calibrate(args), 0)
                self.assertEqual(run.call_args.kwargs['input'], preview['script'])
                receipt = json.loads(Path(str(args.output) + '.submission.json').read_text())
                self.assertEqual(receipt['job_id'], '12345')
                args.execute = False
                args.discard_images = False
                with contextlib.redirect_stdout(io.StringIO()) as stream:
                    calibrate(args)
                self.assertIn('--no-discard-images', json.loads(stream.getvalue())['script'])
                args.full = True
                with contextlib.redirect_stdout(io.StringIO()) as stream:
                    calibrate(args)
                self.assertIn('--full', json.loads(stream.getvalue())['script'])

    def test_allocation_caps_cpu_pool_and_requires_grant(self):
        with patch.dict(os.environ, {'SLURM_JOB_ID': '123', 'SLURM_CPUS_PER_TASK': '1'}):
            self.assertEqual(len(cpu_pool()), 1)
        with patch.dict(os.environ, {'SLURM_JOB_ID': '123', 'SLURM_CPUS_PER_TASK': '0'}):
            with self.assertRaises(ValueError):
                cpu_pool()

    def test_worker_rejects_login_node(self):
        parser = argparse.ArgumentParser()
        add_options(parser)
        args = parser.parse_args(['--slurm-worker', '--workdir', '/unused', '--output', '/unused.json'])
        with patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(ValueError, 'single-node'):
            calibrate(args)
