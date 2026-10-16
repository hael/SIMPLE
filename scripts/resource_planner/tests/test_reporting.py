"""Console reports must not turn failed quality gates into recommendations."""
import contextlib
import io
import unittest
from resource_planner.reporting import result_table, print_results


class ReportingChecks(unittest.TestCase):
    def test_all_failed_report_does_not_claim_calibrated_settings(self):
        collection = {'workflows': {'motion_correct': {'status': 'failed', 'profile': None,
                                                       'reason': 'missing simulate_molecule'}}}
        with contextlib.redirect_stdout(io.StringIO()) as stream:
            print_results(collection, 'failure.json')
        report = stream.getvalue()
        self.assertIn('JSON diagnostic report', report)
        self.assertIn('contains no qualified settings', report)
        self.assertNotIn('Provisional =', report)
        self.assertNotIn('auto*', report)
        self.assertNotIn('Created JSON profile:', report)

    def test_table_preserves_modes_and_qualification(self):
        settings = {'nthr':4,'ncunits':2,'job_memory_per_task_mib':2560}
        collection = {'workflows': {
            'motion_correct': {'status':'provisional','profile':{'sizes':{'4096':{'batch':settings,'latency':{**settings,'ncunits':1}}}}},
            'ctf_estimate': {'status':'provisional','profile':{'sizes':{'4096':{**settings,'ncunits':1,'nparts':1}}}},
            'pick': {'status':'unqualified','profile':None,'reason':'quality gate failed'},
            'solve3D': {'status':'failed','profile':None,'reason':'missing output'}}}
        report = result_table(collection)
        self.assertIn('batch', report)
        self.assertIn('latency', report)
        self.assertIn('auto*', report)
        self.assertIn('MiB/worker', report)
        self.assertIn('quality gate failed', report)
        row = next(row for row in report.splitlines() if row.startswith('pick '))
        self.assertIn('unqualified', row)
        self.assertNotIn('2560', row)
        with contextlib.redirect_stdout(io.StringIO()) as stream:
            print_results(collection, 'result.json')
        self.assertIn('Created JSON profile:', stream.getvalue())
        self.assertIn('does not read/apply', stream.getvalue())
