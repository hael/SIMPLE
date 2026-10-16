"""Hermetic Python checks; no SIMPLE executable or compiler is invoked."""

import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

from resource_planner.adapters import build_command, validate_case
from resource_planner.campaign import execute, load_plan, run_process


class CampaignChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture = self.root / "fixture"
        fixture.mkdir()
        (fixture / "input").write_text("immutable\n")
        # Invoke Python directly: temporary filesystems may be mounted noexec.
        self.fake_command = [sys.executable, '-c', "from pathlib import Path; Path('product').write_text('correct')"]
        self.case = {"id": "correction", "commander": "ctf_correct", "mode": "shared",
                     "work_objects": 8, "project": "input.simple",
                     "science": {"ctf_correct_mode": "phaseflip"},
                     "candidates": [{"nthr": 1}],
                     "prepare": [[sys.executable, "-c", "from pathlib import Path; Path('input.simple').write_text('fixture')"]],
                     "validate": [sys.executable, "-c",
                         "import json; from pathlib import Path; print(json.dumps({'passed':Path('product').read_text()=='correct','validator_version':'fixture-1','metrics':{'pixel_error':0}}))"]}
        self.manifest = {"schema_version": 1, "backend": "local", "fixture": "fixture",
                         "executables": {"simple": sys.executable, "private": sys.executable},
                         "seed": 17, "repetitions": 3,
                         "budget": {"cpus": 1, "wall_seconds": 100, "cpu_seconds": 100,
                                    "attempt_seconds": 10, "disk_bytes": 1000000},
                         "cases": [self.case]}
        self.path = self.root / "campaign.json"

    def plan(self):
        self.path.write_text(json.dumps(self.manifest))
        return load_plan(self.path)

    def test_dry_run_has_no_side_effects(self):
        self.plan()
        output = self.root / "output"
        result = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[2] / "hardware_resource_planner.py"),
                                 "campaign", str(self.path), "--output", str(output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(output.exists())
        self.assertFalse(json.loads(result.stdout)["execution_enabled"])

    def test_success_resume_and_identity(self):
        mock = patch('resource_planner.campaign.build_command', return_value=self.fake_command)
        mock.start()
        self.addCleanup(mock.stop)
        plan = self.plan()
        output = self.root / "output"
        self.assertEqual(execute(plan, output), 0)
        original = (output / "observations.jsonl").read_bytes()
        self.assertEqual(execute(plan, output), 0)
        self.assertEqual((output / "observations.jsonl").read_bytes(), original)
        records = [json.loads(line) for line in original.splitlines()]
        finished = [r for r in records if r['event'] == 'finished']
        self.assertEqual(len(finished), 4)  # One warmup plus three independent repeats.
        self.assertTrue(all(r['qualification'] == 'unqualified' for r in finished))
        self.assertEqual((self.root / "fixture/input").read_text(), "immutable\n")
        (self.root / "fixture/input").write_text("changed")
        with self.assertRaisesRegex(ValueError, "resume"):
            execute(self.plan(), output)

    def test_wrong_science_output_fails(self):
        self.case['validate'] = [sys.executable, '-c',
            "import json; print(json.dumps({'passed':False,'validator_version':'v1','metrics':{'error':1}}))"]
        with patch('resource_planner.campaign.build_command', return_value=self.fake_command):
            self.assertEqual(execute(self.plan(), self.root / 'output'), 2)

    def test_capabilities_and_budgets(self):
        self.case['candidates'][0]['nparts'] = 1
        with self.assertRaisesRegex(ValueError, 'exactly'):
            self.plan()
        del self.case['candidates'][0]['nparts']
        self.manifest['budget']['wall_seconds'] = 39
        with self.assertRaisesRegex(ValueError, 'budget'):
            self.plan()

    def test_worker_and_master_routing(self):
        for name, science in [('motion_correct', {'algorithm':'patch','mcpatch':'yes','downscale':'yes'}),
                              ('ctf_estimate', {'ctfpatch':'yes','pspecsz':512,'fit_phshift':'no',
                                                'dfmin':0.5,'dfmax':5,'hp':30,'lp':5})]:
            case = {**copy.deepcopy(self.case), 'commander':name, 'mode':'worker', 'science':science}
            validate_case(case, 1)
            worker = build_command(case, {'nthr':1}, '/simple', self.root, self.root)
            self.assertEqual(worker[1], 'prg=' + name)
            self.assertIn('fromp=1', worker)
            self.assertIn('top=8', worker)
            case['mode'] = 'master'
            candidate = {'nthr':1,'nparts':4,'ncunits':1,'job_memory_per_task':1024}
            case['candidates'] = [candidate]
            validate_case(case, 1)
            master = build_command(case, candidate, '/simple', self.root, self.root)
            self.assertEqual(master[1], 'prg=' + name)
            self.assertIn('nparts=4', master)
            self.assertFalse(any(arg.startswith(('part=', 'fromp=', 'top=')) for arg in master))

    def test_timeout(self):
        outcome = run_process([sys.executable, '-c', 'import time; time.sleep(30)'], self.root,
                              os.environ.copy(), time.monotonic() + 0.2, self.root/'out', self.root/'err',
                              self.root, 1000000)
        self.assertEqual(outcome['reason'], 'timeout')
        self.assertNotEqual(outcome['returncode'], 0)

    def test_timeout_stops_descendant(self):
        marker = self.root / 'escaped'
        child = f"import time; from pathlib import Path; time.sleep(0.6); Path({str(marker)!r}).touch()"
        parent = f"import subprocess,sys,time; subprocess.Popen([sys.executable,'-c',{child!r}]); time.sleep(30)"
        run_process([sys.executable, '-c', parent], self.root, os.environ.copy(),
                    time.monotonic() + 0.2, self.root/'out', self.root/'err', self.root, 1000000)
        time.sleep(0.7)
        self.assertFalse(marker.exists())

    def test_fixture_symlinks_rejected(self):
        (self.root / 'fixture/link').symlink_to(self.root / 'fixture/input')
        with self.assertRaisesRegex(ValueError, 'symlinks'):
            self.plan()

    def test_master_execution_rejected_before_output(self):
        self.case.update(commander='motion_correct', mode='master',
                         science={'algorithm':'iso','mcpatch':'no','downscale':'yes'},
                         candidates=[{'nthr':1,'nparts':2,'ncunits':1,'job_memory_per_task':1024}])
        output = self.root / 'output'
        with self.assertRaisesRegex(ValueError, 'render-only'):
            execute(self.plan(), output)
        self.assertFalse(output.exists())

    def test_failed_retries_consume_budget(self):
        self.manifest['budget'].update(wall_seconds=40, cpu_seconds=40)
        self.case['prepare'] = [[sys.executable, '-c', 'raise SystemExit(1)']]
        plan = self.plan()
        output = self.root / 'output'
        self.assertEqual(execute(plan, output), 2)
        with self.assertRaisesRegex(ValueError, 'budget exhausted'):
            execute(plan, output)


if __name__ == '__main__':
    unittest.main()
