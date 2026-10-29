"""Storage probes preserve user files and distinguish cache-sensitive measurements."""
import argparse
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from resource_planner import storage


class StorageTests(unittest.TestCase):
    def test_measurements_and_cleanup(self):
        with tempfile.TemporaryDirectory() as root:
            sentinel = Path(root)/'keep'
            sentinel.write_text('user data')
            result = storage.benchmark(root, size_mib=1, repeats=2, small_files=3)
            self.assertEqual(len(result['measurements']), 12)
            self.assertEqual(len(result['summary']), 6)
            self.assertTrue(all(r['seconds'] > 0 for r in result['measurements']))
            self.assertEqual({p.name for p in Path(root).iterdir()}, {'keep'})
            self.assertEqual(sentinel.read_text(), 'user data')
            self.assertIsNone(result['shared_across_nodes'])

    def test_failure_cleans_only_probe_files(self):
        with tempfile.TemporaryDirectory() as root:
            sentinel = Path(root)/'keep'
            sentinel.touch()
            with patch.object(storage, 'read_files', side_effect=ValueError('bad content')):
                with self.assertRaisesRegex(ValueError, 'bad content'):
                    storage.benchmark(root, 1, 1, 2)
            self.assertEqual(list(Path(root).iterdir()), [sentinel])

    def test_disk_reserve_before_writing(self):
        with tempfile.TemporaryDirectory() as root:
            with patch.object(storage.shutil, 'disk_usage') as usage:
                usage.return_value.free = 1
                with self.assertRaisesRegex(ValueError, 'free space'):
                    storage.benchmark(root, 1, 1, 1)
            self.assertEqual(list(Path(root).iterdir()), [])

    def test_options_and_failure_isolation(self):
        parser = argparse.ArgumentParser()
        storage.add_options(parser)
        args = parser.parse_args(['--io-dirs', '/missing', '/second'])
        self.assertEqual(args.io_size_mib, 256)
        with patch.object(storage, 'benchmark', side_effect=[OSError('unavailable'), {'status': 'measured'}]):
            with contextlib.redirect_stdout(io.StringIO()):
                results = storage.run(args)
        self.assertEqual([r['status'] for r in results], ['failed', 'measured'])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser.parse_args(['--io-repeats', '0'])

    def test_slurm_forwards_storage_options(self):
        from resource_planner.calibration import add_options, calibrate
        parser = argparse.ArgumentParser()
        add_options(parser)
        with tempfile.TemporaryDirectory() as root:
            args = parser.parse_args(['--slurm', 'csbdevel', '--io-dirs', root,
                                      '--io-size-mib', '8', '--output', root+'/out.json',
                                      '--workdir', root+'/work'])
            with patch.dict('os.environ', {}, clear=True), contextlib.redirect_stdout(io.StringIO()) as out:
                calibrate(args)
            self.assertIn('--io-dirs', out.getvalue())
            self.assertIn('--io-size-mib 8', out.getvalue())
