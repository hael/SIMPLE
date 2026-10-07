"""Stream snapshots: 2D and 3D particle sets, their requests to the stream, the stages' reports,
and the stage folders they are linked from."""

import os
import tempfile
from importlib import import_module

from django.apps import apps
from django.template.loader import render_to_string
from django.test import SimpleTestCase, TestCase
from django.urls import reverse
from django.utils import timezone

from ..data_structures.streamjob import StreamJob, snapshot_stage_dir
from ..models import JobModel, ProjectModel, WorkspaceModel
from ..views import stream_views

snapshot_migration = import_module("..migrations.0008_snapshot_types", __package__)


class StreamSnapshotTests(TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.project = ProjectModel.objects.create(
            name="project",
            dirc=self.tempdir.name,
            date=timezone.now(),
        )
        os.mkdir(os.path.join(self.tempdir.name, ".workspace_1"))
        self.workspace = WorkspaceModel.objects.create(
            proj=self.project,
            dirc=".workspace_1",
            user="tester",
            jcnt=1,
        )
        self.jobmodel = JobModel.objects.create(
            dset=self.workspace,
            cdat=timezone.now(),
            disp=1,
            dirc="1_simple_stream",
            pckg="simple_stream",
            prog="stream",
        )

    def _streamjob(self):
        return StreamJob(id=self.jobmodel.id)

    def test_snapshot_stage_dirs(self):
        self.assertEqual(snapshot_stage_dir("snapshot2D"), "classification_2D")
        self.assertEqual(snapshot_stage_dir("snapshot3D"), "solve3D_multistate")
        self.assertIsNone(snapshot_stage_dir("snapshot"))
        self.assertIsNone(snapshot_stage_dir("final"))

    def test_2D_snapshot_is_typed_snapshot2D(self):
        streamjob = self._streamjob()
        self.assertTrue(streamjob.snapshot_classification_2D([3, 5], 12))
        particle_set = streamjob.get_jobmodel().particle_sets_stats["particle_sets"][0]
        self.assertEqual(particle_set["type"], "snapshot2D")

    def test_3D_snapshot_request_is_recorded_and_sent_once(self):
        streamjob = self._streamjob()
        self.assertTrue(streamjob.snapshot_solve3D([3, 1]))
        particle_set = streamjob.get_jobmodel().particle_sets_stats["particle_sets"][0]
        self.assertEqual(particle_set["type"], "snapshot3D")
        self.assertEqual(particle_set["states"], [1, 3])
        self.assertEqual(particle_set["filename"], "snapshot_1.simple")
        first = streamjob.get_master_update()
        self.assertEqual(first["snapshot3D"], {"id": 1, "selection": [1, 3], "filename": "snapshot_1.simple"})
        second = streamjob.get_master_update()
        self.assertNotIn("snapshot3D", second)

    def test_3D_snapshot_needs_a_state(self):
        streamjob = self._streamjob()
        self.assertFalse(streamjob.snapshot_solve3D([]))
        self.assertNotIn("snapshot3D", streamjob.get_master_update())

    def test_3D_report_fills_the_set_with_a_bare_file_name(self):
        streamjob = self._streamjob()
        streamjob.snapshot_solve3D([1, 3])
        written = "/data/1_simple_stream/solve3D_multistate/snapshots/snapshot_1/snapshot_1.simple"
        self.assertTrue(streamjob.update_stats({"solve3D_multistate": {
            "stage": "running solve3D_addon",
            "snapshot": {"id": 1, "snapshot_filename": written, "snapshot_nptcls": 9000,
                         "snapshot_time": 123, "states": [1, 3]},
        }}))
        jobmodel = streamjob.get_jobmodel()
        particle_set = jobmodel.particle_sets_stats["particle_sets"][0]
        self.assertEqual(particle_set["filename"], "snapshot_1.simple")
        self.assertEqual(particle_set["path"], written)
        self.assertEqual(particle_set["nptcls"], 9000)
        self.assertEqual(particle_set["time"], 123)
        self.assertNotIn("snapshot", jobmodel.solve3D_multistate_stats)

    def test_unwritten_report_leaves_nothing_to_link(self):
        streamjob = self._streamjob()
        streamjob.snapshot_solve3D([2])
        streamjob.update_stats({"solve3D_multistate": {
            "snapshot": {"id": 1, "snapshot_filename": "", "snapshot_nptcls": 0,
                         "snapshot_time": 123, "states": [2]},
        }})
        particle_set = streamjob.get_jobmodel().particle_sets_stats["particle_sets"][0]
        self.assertEqual(particle_set["nptcls"], 0)
        self.assertNotIn("filename", particle_set)

    def test_selection_after_finish_needs_the_multistate_project(self):
        streamjob = self._streamjob()
        self.assertFalse(streamjob.selection_solve3D([1]))

    def test_migration_retypes_2D_snapshots(self):
        self.jobmodel.particle_sets_stats = {"particle_sets": [
            {"id": 2, "type": "snapshot", "filename": "snapshot_2.simple"},
            {"id": 1, "type": "final", "filename": "deselected.txt"},
        ]}
        self.jobmodel.save()
        snapshot_migration.snapshot_to_snapshot2D(apps, None)
        self.jobmodel.refresh_from_db()
        types = [x["type"] for x in self.jobmodel.particle_sets_stats["particle_sets"]]
        self.assertEqual(types, ["snapshot2D", "final"])


class StateListParsingTests(TestCase):
    def test_parse_state_list(self):
        self.assertEqual(stream_views._parse_state_list("[3, 1, 3]", "test"), [1, 3])
        for raw in ("", "not json", "[]", "[0]", "[true]", "[1.5]", "{\"a\": 1}"):
            with self.subTest(raw=raw):
                self.assertIsNone(stream_views._parse_state_list(raw, "test"))


class StateSelectorTemplateTests(SimpleTestCase):
    """The 3D state tiles: a click views a state, a checkbox selects it, and the selection survives
    the streaming zoom page's reload and is sent once."""

    def _stream_context(self, status):
        return {
            "jobid": 7,
            "status": status,
            "jobstats": {"state_volumes": [{"state": 1}, {"state": 2}, {"state": 3}]},
        }

    def test_stream_tiles_select_by_checkbox(self):
        rendered = render_to_string("includes/_cls3D_state_selector.html", self._stream_context("running"))
        for state in (1, 2, 3):
            self.assertIn(f'data-state-include="{state}"', rendered)
        self.assertEqual(rendered.count('type="checkbox" checked'), 3)
        self.assertIn('onchange="toggleStateInclude(this)"', rendered)
        self.assertIn(f'action="{reverse("nice_lite:snapshot_stream_solve3D")}"', rendered)
        self.assertIn('data-submit-label="create snapshot of"', rendered)
        # viewing a state leaves the selection alone
        select_state = rendered.split("const selectState", 1)[1].split("const selectCls3DStage", 1)[0]
        self.assertNotIn("disabledbutton", select_state)
        self.assertNotIn("checked", select_state)
        # the states left out survive a reload, and a selection is sent once
        self.assertIn("sessionStorage.setItem", rendered)
        self.assertIn("cls3D_excluded:7:", rendered)
        self.assertIn("element.disabled = true", rendered)

    def test_finished_stream_posts_a_final_selection(self):
        rendered = render_to_string("includes/_cls3D_state_selector.html", self._stream_context("finished"))
        self.assertIn(f'action="{reverse("nice_lite:select_stream_solve3D")}"', rendered)
        self.assertIn('data-submit-label="create final particle set of"', rendered)
        self.assertNotIn(reverse("nice_lite:snapshot_stream_solve3D"), rendered)

    def test_batch_stage_groups_select_by_checkbox(self):
        rendered = render_to_string("includes/_cls3D_state_selector.html", {
            "jobid": 9,
            "jobstats": {"cls3D": {"stage1": [{"state": 1}, {"state": 2}]}},
            "cls3d_stages": [{"key": "stage1", "states": [{"state": 1}, {"state": 2}]}],
        })
        self.assertIn('data-stage="stage1"', rendered)
        self.assertEqual(rendered.count('type="checkbox" checked'), 2)
        self.assertIn(f'action="{reverse("nice_lite:batch_cls3D_selection", args=[9])}"', rendered)
        self.assertIn('data-submit-label="save selection of"', rendered)
