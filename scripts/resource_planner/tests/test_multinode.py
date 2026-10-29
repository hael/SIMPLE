"""Verify scheduler identities independently of site-specific OS hostnames."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from resource_planner.multinode import placement


class PlacementTests(unittest.TestCase):
    def test_aliases_and_rejected_incomplete_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'launch.json').write_text(json.dumps({'job_id': '123'}))
            rows = [dict(host='fsitgl-hpc127p', slurm_node='cn127', job_id='123', step_id='1', exit_code=0),
                    dict(host='fsitgl-hpc128p', slurm_node='cn128', job_id='123', step_id='2', exit_code=0)]
            for i, row in enumerate(rows):
                (root/f'{i}.placement.json').write_text(json.dumps(row))
            self.assertEqual(len(placement(root, ['cn127', 'cn128'])), 2)
            for field, value in [('slurm_node', 'cn999'), ('job_id', '456'), ('step_id', None), ('exit_code', 1)]:
                changed = {**rows[1], field: value}
                (root/'1.placement.json').write_text(json.dumps(changed))
                with self.assertRaises(ValueError):
                    placement(root, ['cn127', 'cn128'])
            (root/'1.placement.json').unlink()
            with self.assertRaises(ValueError):
                placement(root, ['cn127', 'cn128'])

    def test_worker_records_scheduler_and_os_names_separately(self):
        from resource_planner.node_launch import main
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'launch.json').write_text(json.dumps({'job_id': '123'}))
            env = {'SIMPLE_CALIBRATION_LAUNCH': temp, 'SLURM_JOB_ID': '123',
                   'SLURMD_NODENAME': 'cn127', 'SLURM_STEP_ID': '4'}
            with patch.dict('os.environ', env), patch('socket.gethostname', return_value='fsitgl-hpc127p.site'), \
                    patch('resource_planner.node_launch.subprocess.call', return_value=0):
                self.assertEqual(main(['--worker', 'worker.sh', 'token']), 0)
            row = json.loads((root/'token.placement.json').read_text())
            self.assertEqual(row['host'], 'fsitgl-hpc127p')
            self.assertEqual(row['slurm_node'], 'cn127')
            self.assertEqual(row['exit_code'], 0)
