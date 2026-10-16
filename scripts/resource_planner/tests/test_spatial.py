"""Spatial calibration adapter contracts; no SIMPLE executions."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest

from resource_planner.workflows.picking import PickAdapter, picking_metrics, validate_pick
from resource_planner.workflows.solve3d import Solve3DAdapter
from resource_planner.workflows.common import export_shared
from resource_planner.adapters import build_command, validate_case


class FakeContext:
    def __init__(self, root):
        self.root=root
        self.args=argparse.Namespace(sizes=[4096],threads=[4],box=128,particles=2000,
                                    simple_exec=Path('/simple_exec'),fixture_python='python3',timeout=600)
        self.commands=[]
    def command(self,argv,directory,label,timeout=60):
        self.commands.append((argv,directory,label))
    def project(self,directory,name,image,particles,ctf_flag,oritab=None):
        path=directory/name
        path.mkdir()
        return path/(name+'.simple')


class SpatialChecks(unittest.TestCase):
    def test_native_simulation_uses_one_thread_and_balanced_import_halves(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);ctx=FakeContext(root)
            (root/'fixture').mkdir()
            fixture=Solve3DAdapter().prepare(ctx,128)
            simulation=next(argv for argv,_,label in ctx.commands if label=='simulate')
            self.assertIn('prg=simulate_particles',simulation)
            self.assertIn('nthr=1',simulation)
            halves=(root/'fixture/halves.txt').read_text().splitlines()
            self.assertEqual(len(halves),2000)
            self.assertEqual(halves.count('state=1 eo=0'),1000)
            self.assertEqual(halves.count('state=1 eo=1'),1000)
            self.assertIn(root/'fixture/particles.mrcs',fixture['files'])

    def test_picking_duplicates_are_false_positives(self):
        expected=[(10,10),(100,100)]
        self.assertTrue(picking_metrics(expected,expected,5)['quality_qualified'])
        result=picking_metrics([(10,10)]*4,expected,5)
        self.assertEqual(result['matched'],1)
        self.assertEqual(result['precision'],0.25)
        self.assertFalse(result['quality_qualified'])
        self.assertFalse(picking_metrics([],expected,5)['quality_qualified'])

    def test_pick_box_origin_conversion_and_nonfinite_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);box=root/'picked.box'
            box.write_text('0 0 20 20\n90 90 20 20\n')
            truth={'box':20,'coordinates':[[0,0,20,20],[90,90,20,20]]}
            self.assertTrue(validate_pick(root,truth)['quality_qualified'])
            box.write_text('nan 0 20 20\n')
            with self.assertRaises(ValueError):validate_pick(root,truth)

    def test_pick_and_3d_commands_are_isolated_shared_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);ctx=FakeContext(root)
            fixture={'input':root/'fixture'}
            for adapter in (PickAdapter(),Solve3DAdapter()):
                work=root/adapter.name;work.mkdir()
                argv,cwd=adapter.jobs(ctx,fixture,work,4,1)[0]
                self.assertEqual(argv[1],'prg='+adapter.name)
                self.assertIn('nthr=4',argv)
                if adapter.name == 'solve3D':
                    self.assertFalse(any(a.startswith('nparts=') for a in argv))
                else:
                    self.assertIn('nparts=1', argv)
                    self.assertIn('part=1', argv)
                self.assertEqual(cwd,work)
                if adapter.name == 'pick':
                    self.assertIn('fromp=1',argv);self.assertIn('top=1',argv)
                if adapter.name=='solve3D':
                    self.assertIn('nstages=4',argv)
                    self.assertIn('vol1='+str(fixture['input']/'volume.mrc'),argv)
                    self.assertIn('balance=none',argv)
                    self.assertEqual(ctx.commands,[])

    def test_generic_campaign_capabilities(self):
        for name,science in [('pick',{'picker':'new','moldiam':100,'pcontrast':'black','lp':20}),
                             ('solve3D',{'pgrp':'c1','nstates':1,'nstages':1,'nsample':100,'mskdiam':100,'vol1':'/truth.mrc'})]:
            case={'commander':name,'mode':'shared','science':science,'work_objects':10,
                  'project':'project.simple','candidates':[{'nthr':4}]}
            validate_case(case,4)
            argv=build_command(case,{'nthr':4},Path('/simple'),Path('/work'),Path('/fixture'))
            self.assertEqual(argv[1],'prg='+name)
            if name!='solve3D':self.assertIn('top=10',argv)

    def test_export_preserves_quality_and_scope(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);ctx=FakeContext(root);adapter=PickAdapter()
            plan={'repetitions':3,'binaries':{},'hardware':{},'cpus':[0,1,2,3],
                  'workload':adapter.metadata(ctx.args)['workload']}
            records=[{'size':4096,'nthr':4,'repeat':i,'warmup':i==0,'wall_seconds':2,
                      'peak_rss_mib':800,'passed':True,'validation':{'quality_qualified':False}} for i in range(4)]
            self.assertIsNone(export_shared(adapter,ctx,plan,records))
            self.assertTrue((root/'summary.json').exists())
            self.assertFalse((root/'desktop_profile.json').exists())
            for row in records:row['validation']['quality_qualified']=True
            path=export_shared(adapter,ctx,plan,records)
            saved=json.loads(path.read_text())
            self.assertEqual(saved['commander'],'pick')
            self.assertEqual(saved['sizes']['4096']['job_memory_per_task_mib'],1280)
            self.assertFalse(saved['auto_apply'])
            with self.assertRaises(ValueError):export_shared(adapter,ctx,plan,records[:-1])


if __name__=='__main__':unittest.main()
