"""Independent CTF truth gates and incomplete-export rejection; no SIMPLE runs."""
import json
from pathlib import Path
import tempfile
import unittest

from resource_planner.workflows.ctf import export_profile, validate


class CtfDesktopChecks(unittest.TestCase):
    def test_defocus_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            doc=root/'micrograph_ctf.txt'
            for fitted, expected in ((2.001,True),(2.5,False),(float('nan'),False)):
                doc.write_text(f'dfx={fitted} dfy=2.0 ctfres=3.0 ctfcc=0.5\n')
                self.assertEqual(validate(root,{'dfx':2.,'dfy':2.},4096)['passed'],expected)

    def test_missing_document_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                validate(Path(temporary),{'dfx':2.,'dfy':2.},4096)

    def test_profile_requires_complete_measurements(self):
        plan={'repetitions':3,'sizes':[4096],'candidates':[(4,1)]}
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                export_profile(Path(temporary),plan,[],[])
            self.assertFalse((Path(temporary)/'desktop_profile.json').exists())

    def test_json_profile_and_memory_margin(self):
        plan={'repetitions':3,'sizes':[4096],'candidates':[(4,1)],
              'binaries':{},'cpus':[0,1,2,3],'workload':{'smpd':1.3}}
        records=[{'size':4096,'nthr':4,'workers':1,'repeat':i,'passed':True} for i in range(4)]
        summary=[{'size':4096,'nthr':4,'workers':1,'max_worker_mib':1000,
                  'median_micrographs_per_second':0.25,'median_wave_seconds':4}]
        with tempfile.TemporaryDirectory() as temporary:
            path=export_profile(Path(temporary),plan,records,summary)
            profile=json.loads(path.read_text())
            self.assertEqual(profile['commander'],'ctf_estimate')
            self.assertFalse(profile['auto_apply'])
            self.assertEqual(profile['sizes']['4096']['batch']['job_memory_per_task_mib'],1536)
            self.assertFalse(profile['resource_policy']['partition_policy_measured'])


if __name__=='__main__':
    unittest.main()
