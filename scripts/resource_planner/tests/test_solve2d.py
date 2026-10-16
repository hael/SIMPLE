"""Solve2D routing, label-invariant quality metrics and class-stack validation."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from resource_planner.adapters import build_command, validate_case
from resource_planner.workflows.solve2d import assignment_metrics, export, science, validate


class Solve2DChecks(unittest.TestCase):
    def test_shared_routing(self):
        case={'commander':'solve2D','mode':'shared','work_objects':2000,
              'project':'classes.simple','science':science(2000,128),'candidates':[{'nthr':8}]}
        validate_case(case,32)
        command=build_command(case,{'nthr':8},Path('/simple'),Path('/work'),Path('/fixture'))
        self.assertEqual(command[1],'prg=solve2D')
        for key in ('outstk=','part=','nparts=','fromp=','top='):
            self.assertFalse(any(a.startswith(key) for a in command))
        case['mode']='master'
        with self.assertRaises(ValueError):validate_case(case,32)

    def test_label_permutation(self):
        metrics=assignment_metrics([1,1,2,2],[8,8,3,3])
        self.assertEqual(metrics['adjusted_rand_index'],1)
        self.assertEqual(metrics['purity'],1)
        mixed=assignment_metrics([1,1,2,2],[1,2,1,2])
        self.assertAlmostEqual(mixed['adjusted_rand_index'],-0.5)
        self.assertEqual(mixed['purity'],0.5)
        with self.assertRaises(ValueError):assignment_metrics([1,2],[1])

    def test_default_off(self):
        with tempfile.TemporaryDirectory() as temporary:
            output=Path(temporary)/'untouched'
            result=subprocess.run([sys.executable,str(Path(__file__).resolve().parents[2]/'hardware_resource_planner.py'),
                                   'calibrate','--fixture-source','analytic','--workflows','solve2D','--full','--output',str(output)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(json.loads(result.stdout)['workflows']['solve2D']['particles'],2000)
            self.assertFalse(output.exists())

    def test_qualified_profile_selection(self):
        plan={'threads':[4,8],'repetitions':3,'simple_exec':'/simple',
              'binary_sha256':'digest','hardware':{},'cpus':list(range(8)),
              'particles':2000,'box':128,'science':science(2000,128)}
        records=[{'nthr':nthr,'warmup':False,'wall_seconds':seconds,'peak_rss_mib':800,
                  'validation':{'passed':True,'quality_qualified':True,
                                'adjusted_rand_index':0.95,'purity':0.97}}
                 for nthr,seconds in ((4,10.4),(8,10)) for _ in range(3)]
        with tempfile.TemporaryDirectory() as temporary, \
                patch('resource_planner.workflows.solve2d.file_hash',return_value='digest'):
            profile=json.loads(export(Path(temporary),plan,records).read_text())
            self.assertEqual(profile['settings']['nthr'],4)  # Within 5%, uses fewer threads.
            self.assertEqual(profile['settings']['nparts'],1)
            self.assertEqual(profile['settings']['job_memory_per_task_mib'],1280)
            self.assertFalse(profile['auto_apply'])

    def test_quality_failure_blocks_recommendations(self):
        records=[{'nthr':4,'warmup':False,'wall_seconds':10,'peak_rss_mib':100,
                  'validation':{'passed':True,'quality_qualified':False,
                                'adjusted_rand_index':0.3,'purity':0.5}} for _ in range(3)]
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            self.assertIsNone(export(root,{'threads':[4],'repetitions':3},records))
            self.assertTrue((root/'summary.json').exists())
            self.assertFalse((root/'desktop_profile.json').exists())

    def test_missing_class_averages(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError,'missing class-average'):
                validate(Path(temporary),{},Path(temporary)/'assignments')


if __name__=='__main__':unittest.main()
