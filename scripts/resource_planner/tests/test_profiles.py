"""JSON profile export checks; no native commands."""
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch
from resource_planner.profiles import write_profile

class ProfileChecks(unittest.TestCase):
    def test_export_is_json_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            plan={'repetitions':3,'binaries':{'/original/simple_exec':'digest',
                                             '/original/simple_private_exec':'digest'},
                  'sizes':[4096],'cpus':list(range(32)),'known_shifts_x':[0]*8}
            summary=[{'size':4096,'nthr':4,'workers':8,'median_movies_per_second':0.8,
                      'median_wave_seconds':10,'max_worker_mib':1550},
                     {'size':4096,'nthr':8,'workers':1,'median_movies_per_second':0.2,
                      'median_wave_seconds':5,'max_worker_mib':1550}]
            rows=[]
            for s in summary:
                for repeat in range(3):
                    rows.append({**s,'repeat':repeat+1,'warmup':False,'passed':True,
                                 'validation':[{'shifts_x':[0]*8,'shifts_y':[0]*8}]*s['workers']})
            (root/'plan.json').write_text(json.dumps(plan))
            (root/'summary.json').write_text(json.dumps(summary))
            (root/'results.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
            with patch('resource_planner.profiles.file_hash', return_value='digest'), \
                 patch('resource_planner.profiles.elf_signature', return_value={'text':'same'}):
                profile=write_profile(root,root/'install')
            self.assertEqual(json.loads(profile.read_text())['sizes']['4096']['batch']['nthr'],4)
            saved=json.loads(profile.read_text())
            self.assertEqual(saved['commander'],'motion_correct')
            self.assertFalse(saved['auto_apply'])
            self.assertFalse(saved['resource_policy']['partition_policy_measured'])
            self.assertFalse((root/'run_motion_desktop.py').exists())

