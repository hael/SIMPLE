import os
import signal
import struct
import subprocess
import tempfile
from unittest.mock import ANY, patch

from django.test import TestCase, override_settings
from django.utils import timezone

from ..data_structures import batchjob as batchjob_module
from ..data_structures import simple as simple_module
from ..data_structures.batchjob import BatchJob
from ..data_structures.simple import SIMPLEBatch
from ..data_structures.workspace import Workspace
from ..models import JobModel, ProjectModel, WorkspaceModel
from ..views import job_builder_views


class BatchJobLifecycleTests(TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.project = ProjectModel.objects.create(
            name="project",
            dirc=self.tempdir.name,
            date=timezone.now(),
        )
        self.workspace_dir = os.path.join(self.tempdir.name, ".workspace_1")
        os.mkdir(self.workspace_dir)
        with open(os.path.join(self.workspace_dir, "workspace.simple"), "w", encoding="utf-8"):
            pass
        self.workspace_model = WorkspaceModel.objects.create(
            proj=self.project,
            dirc=".workspace_1",
            user="tester",
        )
        self.workspace = Workspace(self.workspace_model.id)

    def test_new_persists_batch_identity_and_launches_selected_commander(self):
        launcher = batchjob_module.SIMPLEBatch
        with patch.object(launcher, "loadUIJSON", return_value=True), patch.object(launcher, "start", return_value=True) as start:
            job = BatchJob()
            created = job.new(
                self.workspace,
                "single",
                "pick",
                {"mode": "fast"},
                display_name="Pick Particles",
            )

        self.assertTrue(created)
        jobmodel = JobModel.objects.get(id=job.id)
        self.assertEqual(jobmodel.args, {"mode": "fast"})
        self.assertEqual(jobmodel.name, "Pick Particles")
        self.assertEqual(jobmodel.master_stats["job_type"], "batch")
        self.assertEqual(jobmodel.master_stats["package"], "single")
        self.assertEqual(jobmodel.master_stats["program"], "pick")
        self.assertEqual(jobmodel.status, "queued")
        start.assert_called_once_with(
            {"mode": "fast"},
            os.path.join(self.workspace_dir, "1_pick"),
            self.workspace_dir,
            "pick",
            job.id,
        )
        loaded_job = BatchJob(id=job.id)
        self.assertEqual(loaded_job.prog, "pick")
        self.assertEqual(loaded_job.name, "Pick Particles")

    def test_batch_detail_resolves_named_new_project_in_job_directory(self):
        job_dir = os.path.join(self.workspace_dir, "1_new_project")
        os.mkdir(job_dir)
        project_path = os.path.join(job_dir, "test.simple")
        with open(project_path, "w", encoding="utf-8"):
            pass
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_new_project",
            args={"projname": "test"},
            status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "new_project"},
        )

        job = BatchJob(id=jobmodel.id)

        self.assertEqual(job.get_safe_job_dir(), job_dir)
        self.assertEqual(job.get_result_project_path(), project_path)

    def test_batch_detail_rejects_ambiguous_or_external_result_projects(self):
        job_dir = os.path.join(self.workspace_dir, "1_demo")
        os.mkdir(job_dir)
        for filename in ("first.simple", "second.simple"):
            with open(os.path.join(job_dir, filename), "w", encoding="utf-8"):
                pass
        outside_project = os.path.join(self.tempdir.name, "outside.simple")
        with open(outside_project, "w", encoding="utf-8"):
            pass
        os.symlink(outside_project, os.path.join(job_dir, "external.simple"))
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_demo",
            status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "demo"},
        )
        job = BatchJob(id=jobmodel.id)

        self.assertIsNone(job.get_result_project_path())

        workspace_project = os.path.join(job_dir, "workspace.simple")
        with open(workspace_project, "w", encoding="utf-8"):
            pass
        self.assertEqual(job.get_result_project_path(), workspace_project)

    def test_batch_detail_returns_bounded_logs_and_artifact_previews(self):
        job_dir = os.path.join(self.workspace_dir, "1_motion_correct")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "stdout.log"), "w", encoding="utf-8") as stdout_file:
            stdout_file.write("0123456789")
        with open(os.path.join(job_dir, "movie_thumb.jpg"), "wb") as image_file:
            image_file.write(b"image")
        with open(os.path.join(job_dir, "movie_intg.mrc"), "wb") as mrc_file:
            mrc_file.write(b"mrc")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_motion_correct",
            status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "motion_correct"},
        )

        job = BatchJob(id=jobmodel.id)
        logs = job.get_log_tails(max_bytes=5)
        artifacts = job.get_artifact_summary()

        self.assertEqual(logs[0]["text"], "56789")
        self.assertTrue(logs[0]["truncated"])
        self.assertFalse(logs[1]["exists"])
        self.assertEqual(
            {entry["extension"]: entry["count"] for entry in artifacts["counts"]},
            {"JPG": 1, "MRC": 1},
        )
        self.assertEqual(artifacts["images"][0]["name"], "movie_thumb.jpg")
        previews = artifacts["images"][0]["previews"]
        self.assertEqual(
            [(preview["kind"], preview["crop"]) for preview in previews],
            [
                ("power_spectrum", "left"),
                ("micrograph", "right"),
            ],
        )
        self.assertEqual(
            [preview["visibility_group"] for preview in previews],
            ["motion", "motion"],
        )
        self.assertNotIn("hidden_by_default", previews[0])
        self.assertTrue(previews[1]["hidden_by_default"])

    def test_ctf_artifacts_pair_diagnostics_with_source_motion_thumbnails(self):
        motion_dir = os.path.join(self.workspace_dir, "1_motion_correct")
        os.mkdir(motion_dir)
        motion_thumbnail = os.path.join(motion_dir, "movie_thumb.jpg")
        with open(motion_thumbnail, "wb") as image_file:
            image_file.write(b"motion")
        motion_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_motion_correct",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "motion_correct",
            },
        )

        ctf_dir = os.path.join(self.workspace_dir, "2_ctf_estimate")
        os.mkdir(ctf_dir)
        diagnostic = os.path.join(ctf_dir, "movie_ctf_estimate_diag.jpg")
        with open(diagnostic, "wb") as image_file:
            image_file.write(b"ctf")
        ctf_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=2,
            dirc="2_ctf_estimate",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "ctf_estimate",
                "source": {"type": "batch_job", "batch_job_id": motion_job.id},
            },
        )

        job = BatchJob(id=ctf_job.id)
        artifacts = job.get_artifact_summary()

        self.assertEqual(artifacts["images"][0]["path"], diagnostic)
        previews = artifacts["images"][0]["previews"]
        self.assertEqual(
            [(preview["kind"], preview["crop"]) for preview in previews],
            [
                ("image", "full"),
                ("micrograph", "right"),
            ],
        )
        self.assertEqual(previews[1]["name"], "movie_thumb.jpg")
        self.assertEqual(previews[1]["path"], motion_thumbnail)
        self.assertTrue(previews[1]["hidden_by_default"])
        self.assertEqual(previews[1]["visibility_group"], "ctf_micrograph")

    def test_ctf_artifacts_do_not_pair_with_a_non_motion_source(self):
        source_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(source_dir)
        with open(os.path.join(source_dir, "movie_thumb.jpg"), "wb") as image_file:
            image_file.write(b"not motion")
        source_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "import_movies",
            },
        )

        ctf_dir = os.path.join(self.workspace_dir, "2_ctf_estimate")
        os.mkdir(ctf_dir)
        with open(os.path.join(ctf_dir, "movie_ctf_estimate_diag.jpg"), "wb") as image_file:
            image_file.write(b"ctf")
        ctf_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=2,
            dirc="2_ctf_estimate",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "ctf_estimate",
                "source": {"type": "batch_job", "batch_job_id": source_job.id},
            },
        )

        artifacts = BatchJob(id=ctf_job.id).get_artifact_summary()

        self.assertEqual(
            [preview["kind"] for preview in artifacts["images"][0]["previews"]],
            ["image"],
        )

    def test_pick_previews_follow_source_chain_and_read_owned_artifacts(self):
        motion_dir = os.path.join(self.workspace_dir, "1_motion_correct")
        os.mkdir(motion_dir)
        thumbnail_path = os.path.join(motion_dir, "movie_thumb.jpg")
        with open(thumbnail_path, "wb") as thumbnail_file:
            thumbnail_file.write(b"image")
        with open(os.path.join(motion_dir, "movie_intg.mrc"), "wb") as mrc_file:
            mrc_file.write(struct.pack("<3i", 4096, 3072, 1))
        motion_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_motion_correct",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "motion_correct",
            },
        )

        ctf_dir = os.path.join(self.workspace_dir, "2_ctf_estimate")
        os.mkdir(ctf_dir)
        ctf_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=2,
            dirc="2_ctf_estimate",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "ctf_estimate",
                "source": {"type": "batch_job", "batch_job_id": motion_job.id},
            },
        )

        pick_dir = os.path.join(self.workspace_dir, "3_pick")
        os.mkdir(pick_dir)
        with open(os.path.join(pick_dir, "movie_intg.box"), "w", encoding="utf-8") as box_file:
            box_file.write("-10 -20 40 60 99.0\ninvalid record\n100 200 20 40\n")
        pick_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=3,
            dirc="3_pick",
            status="finished",
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "pick",
                "source": {"type": "batch_job", "batch_job_id": ctf_job.id},
            },
        )

        previews = BatchJob(id=pick_job.id).get_pick_micrograph_previews()

        self.assertEqual(previews, [{
            "path": thumbnail_path,
            "number": 1,
            "xdim": 4096,
            "ydim": 3072,
            "boxes": [
                {"x": 10.0, "y": 10.0, "width": 40.0, "height": 60.0},
                {"x": 110.0, "y": 220.0, "width": 20.0, "height": 40.0},
            ],
        }])

    def test_batch_detail_rejects_job_directory_outside_workspace(self):
        outside_dir = os.path.join(self.tempdir.name, "outside_job")
        os.mkdir(outside_dir)
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="../outside_job",
            status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "demo"},
        )
        job = BatchJob(id=jobmodel.id)

        self.assertIsNone(job.get_safe_job_dir())
        self.assertIsNone(job.get_result_project_path())
        self.assertEqual(job.get_log_tails(), [])
        self.assertEqual(job.get_artifact_summary(), {"counts": [], "images": []})

    def test_new_marks_persisted_job_failed_when_dispatch_fails(self):
        launcher = batchjob_module.SIMPLEBatch
        with patch.object(launcher, "loadUIJSON", return_value=True), patch.object(launcher, "start", return_value=False):
            job = BatchJob()
            created = job.new(self.workspace, "simple", "import_movies", {})

        self.assertFalse(created)
        jobmodel = JobModel.objects.get(id=job.id)
        self.assertEqual(jobmodel.status, "failed")
        self.assertEqual(jobmodel.master_status, "failed")

    def test_stop_terminates_verified_process_group_and_marks_job_stopped(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "nice.pid"), "w", encoding="utf-8") as pid_file:
            pid_file.write("4321\n")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="running",
            master_status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )
        job = BatchJob(id=jobmodel.id)

        with patch.object(job, "_get_local_process_group", return_value=4321), patch.object(batchjob_module.os, "killpg") as killpg:
            stopped = job.stop()

        self.assertTrue(stopped)
        killpg.assert_called_once_with(4321, signal.SIGTERM)
        jobmodel.refresh_from_db()
        self.assertEqual(jobmodel.status, "stopped")
        self.assertEqual(jobmodel.master_status, "stopped")

    def test_stop_rejects_non_running_job_without_signalling(self):
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="finished",
            master_status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        with patch.object(batchjob_module.os, "killpg") as killpg:
            stopped = BatchJob(id=jobmodel.id).stop()

        self.assertFalse(stopped)
        killpg.assert_not_called()

    def test_stop_rejects_pid_from_outside_job_directory(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "nice.pid"), "w", encoding="utf-8") as pid_file:
            pid_file.write("4321\n")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="running",
            master_status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        with patch.object(batchjob_module.os, "readlink", return_value=self.workspace_dir), patch.object(batchjob_module.os, "getpgid", return_value=4321), patch.object(batchjob_module.os, "killpg") as killpg:
            stopped = BatchJob(id=jobmodel.id).stop()

        self.assertFalse(stopped)
        killpg.assert_not_called()

    def test_late_heartbeat_does_not_resurrect_stopped_job(self):
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="stopped",
            master_status="stopped",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        response = BatchJob(id=jobmodel.id).updateStats({"job_heartbeat": {}}, None, None)

        self.assertEqual(response, {})
        jobmodel.refresh_from_db()
        self.assertEqual(jobmodel.status, "stopped")
        self.assertEqual(jobmodel.master_status, "stopped")

    def test_terminal_callback_uses_reported_finished_status(self):
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="running",
            master_status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        BatchJob(id=jobmodel.id).updateStats(
            {"job": {"status": "finished", "terminate": True}},
            None,
            None,
        )

        jobmodel.refresh_from_db()
        self.assertEqual(jobmodel.status, "finished")
        self.assertEqual(jobmodel.master_status, "finished")

    def test_terminal_callback_preserves_reported_failure(self):
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="running",
            master_status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        BatchJob(id=jobmodel.id).updateStats(
            {"job": {"status": "failed", "terminate": True}},
            None,
            None,
        )

        jobmodel.refresh_from_db()
        self.assertEqual(jobmodel.status, "failed")
        self.assertEqual(jobmodel.master_status, "failed")

    def test_late_callback_does_not_replace_terminal_status(self):
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="failed",
            master_status="failed",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        response = BatchJob(id=jobmodel.id).updateStats(
            {"job": {"status": "finished", "terminate": True}},
            None,
            None,
        )

        jobmodel.refresh_from_db()
        self.assertEqual(response, {})
        self.assertEqual(jobmodel.status, "failed")
        self.assertEqual(jobmodel.master_status, "failed")

    def test_reconcile_local_completion_marks_normal_exited_job_finished(self):
        job_dir = os.path.join(self.workspace_dir, "1_new_project")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "job.script"), "w", encoding="utf-8") as script_file:
            script_file.write("#!/bin/sh\n# localtemplate\n")
        with open(os.path.join(job_dir, "nice.pid"), "w", encoding="utf-8") as pid_file:
            pid_file.write("4321\n")
        with open(os.path.join(job_dir, "stdout.log"), "w", encoding="utf-8") as stdout_file:
            stdout_file.write("**** NEW_PROJECT NORMAL STOP ****\n")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_new_project",
            status="queued",
            master_status="queued",
            master_stats={"job_type": "batch", "package": "simple", "program": "new_project"},
        )
        job = BatchJob(id=jobmodel.id)

        with patch.object(job, "_local_process_is_running", return_value=False):
            changed = job.reconcile_local_completion()

        self.assertTrue(changed)
        jobmodel.refresh_from_db()
        self.assertEqual(jobmodel.status, "finished")
        self.assertEqual(jobmodel.master_status, "finished")

    def test_reconcile_local_completion_keeps_active_job_running(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "job.script"), "w", encoding="utf-8") as script_file:
            script_file.write("#!/bin/sh\n# localtemplate\n")
        with open(os.path.join(job_dir, "nice.pid"), "w", encoding="utf-8") as pid_file:
            pid_file.write("4321\n")
        with open(os.path.join(job_dir, "stdout.log"), "w", encoding="utf-8") as stdout_file:
            stdout_file.write("**** IMPORT_MOVIES NORMAL STOP ****\n")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="running",
            master_status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )
        job = BatchJob(id=jobmodel.id)

        with patch.object(job, "_local_process_is_running", return_value=True):
            changed = job.reconcile_local_completion()

        self.assertFalse(changed)
        jobmodel.refresh_from_db()
        self.assertEqual(jobmodel.status, "running")

    def test_queued_job_can_delete_after_local_process_exits(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "job.script"), "w", encoding="utf-8") as script_file:
            script_file.write("#!/bin/sh\n# localtemplate\n")
        with open(os.path.join(job_dir, "nice.pid"), "w", encoding="utf-8") as pid_file:
            pid_file.write("4321\n")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="queued",
            master_status="queued",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )
        job = BatchJob(id=jobmodel.id)

        with patch.object(job, "_local_process_is_running", return_value=False):
            can_delete = job.queued_job_can_delete()

        self.assertTrue(can_delete)

    def test_delete_rejects_active_local_queued_job(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "job.script"), "w", encoding="utf-8") as script_file:
            script_file.write("#!/bin/sh\n# localtemplate\n")
        with open(os.path.join(job_dir, "nice.pid"), "w", encoding="utf-8") as pid_file:
            pid_file.write("4321\n")
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="queued",
            master_status="queued",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )
        job = BatchJob(id=jobmodel.id)

        with patch.object(job, "_local_process_is_running", return_value=True):
            deleted = job.delete(None, self.workspace)

        self.assertFalse(deleted)
        self.assertTrue(os.path.isdir(job_dir))
        self.assertTrue(JobModel.objects.filter(id=jobmodel.id).exists())

    def test_delete_rejects_running_batch_job(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="running",
            master_status="running",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        deleted = BatchJob(id=jobmodel.id).delete(None, self.workspace)

        self.assertFalse(deleted)
        self.assertTrue(os.path.isdir(job_dir))
        self.assertTrue(JobModel.objects.filter(id=jobmodel.id).exists())

    def test_delete_permanently_removes_batch_directory_without_trash(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        with open(os.path.join(job_dir, "result.txt"), "w", encoding="utf-8"):
            pass
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="1_import_movies",
            status="finished",
            master_status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )
        self.workspace_model.jcnt = 1
        self.workspace_model.save(update_fields=["jcnt"])

        deleted = BatchJob(id=jobmodel.id).delete(None, self.workspace)

        self.assertTrue(deleted)
        self.assertFalse(os.path.exists(job_dir))
        self.assertFalse(os.path.exists(os.path.join(self.workspace_dir, "TRASH")))
        self.assertFalse(JobModel.objects.filter(id=jobmodel.id).exists())
        self.workspace_model.refresh_from_db()
        self.assertEqual(self.workspace_model.jcnt, 0)

    def test_delete_rejects_batch_directory_outside_workspace(self):
        outside_dir = os.path.join(self.tempdir.name, "outside")
        os.mkdir(outside_dir)
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            dirc="../outside",
            status="finished",
            master_status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        deleted = BatchJob(id=jobmodel.id).delete(None, self.workspace)

        self.assertFalse(deleted)
        self.assertTrue(os.path.isdir(outside_dir))
        self.assertTrue(JobModel.objects.filter(id=jobmodel.id).exists())

    def test_new_initializes_missing_workspace_project_before_dispatch(self):
        os.remove(os.path.join(self.workspace_dir, "workspace.simple"))
        launcher = batchjob_module.SIMPLEBatch
        project_initializer = batchjob_module.SIMPLEProject

        with patch.object(project_initializer, "create", return_value=True) as create, patch.object(launcher, "loadUIJSON", return_value=True), patch.object(launcher, "start", return_value=True):
            created = BatchJob().new(self.workspace, "simple", "import_movies", {})

        self.assertTrue(created)
        create.assert_called_once_with()

    def test_new_keeps_failed_record_when_job_directory_cannot_be_created(self):
        launcher = batchjob_module.SIMPLEBatch
        with patch.object(batchjob_module, "ensure_directory", return_value=False), patch.object(launcher, "loadUIJSON", return_value=True), patch.object(launcher, "start") as start:
            job = BatchJob()
            created = job.new(self.workspace, "simple", "import_movies", {})

        self.assertFalse(created)
        jobmodel = JobModel.objects.get(id=job.id)
        self.assertEqual(jobmodel.status, "failed")
        self.assertEqual(jobmodel.master_status, "failed")
        self.workspace_model.refresh_from_db()
        self.assertEqual(self.workspace_model.jcnt, 1)
        start.assert_not_called()

    def test_new_launches_from_explicit_stream_snapshot_and_records_source(self):
        snapshot_dir = os.path.join(
            self.workspace_dir,
            "2_simple_stream",
            "classification_2D",
            "snapshots",
            "snapshot_1",
        )
        os.makedirs(snapshot_dir)
        snapshot_path = os.path.join(snapshot_dir, "snapshot_1.simple")
        with open(snapshot_path, "w", encoding="utf-8"):
            pass
        source = {
            "type": "stream_snapshot",
            "stream_job_id": 12,
            "particle_set_id": 1,
            "filename": "snapshot_1.simple",
        }

        launcher = batchjob_module.SIMPLEBatch
        with (
            patch.object(launcher, "loadUIJSON", return_value=True),
            patch.object(launcher, "start", return_value=True) as start,
        ):
            job = BatchJob()
            created = job.new(
                self.workspace,
                "simple",
                "cluster2D",
                {"nthr": "8"},
                parent_proj=snapshot_path,
                source=source,
            )

        self.assertTrue(created)
        start.assert_called_once_with(
            {"nthr": "8"},
            os.path.join(self.workspace_dir, "1_cluster2D"),
            self.workspace_dir,
            "cluster2D",
            job.id,
            parent_proj=snapshot_path,
        )
        jobmodel = JobModel.objects.get(id=job.id)
        self.assertEqual(jobmodel.master_stats["source"], source)
        self.assertEqual(BatchJob(id=job.id).source, source)

    @override_settings(NICE_LITE_BATCH_PROJECT_INHERITANCE=True)
    def test_batch_project_sources_list_latest_completed_job_first(self):
        expected = []
        for disp, name in ((1, "Import Movie Data"), (2, "Correct Movie Motion")):
            job_dir = os.path.join(self.workspace_dir, f"{disp}_job")
            os.mkdir(job_dir)
            project_path = os.path.join(job_dir, "workspace.simple")
            with open(project_path, "w", encoding="utf-8"):
                pass
            jobmodel = JobModel.objects.create(
                dset=self.workspace_model,
                cdat=timezone.now(),
                disp=disp,
                name=name,
                dirc=f"{disp}_job",
                status="finished",
                master_status="finished",
                master_stats={"job_type": "batch", "package": "simple", "program": "demo"},
            )
            expected.insert(0, {
                "key": f"job:{jobmodel.id}",
                "label": f"job {disp} - {name}",
            })

        failed_dir = os.path.join(self.workspace_dir, "3_failed")
        os.mkdir(failed_dir)
        with open(os.path.join(failed_dir, "workspace.simple"), "w", encoding="utf-8"):
            pass
        JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=3,
            name="Failed Job",
            dirc="3_failed",
            status="failed",
            master_status="failed",
            master_stats={"job_type": "batch", "package": "simple", "program": "demo"},
        )

        sources = job_builder_views._collect_batch_job_sources(self.workspace)

        self.assertEqual(sources, expected)
        self.assertEqual(
            job_builder_views._default_batch_project_file(self.workspace),
            os.path.join(self.workspace_dir, "2_job", "workspace.simple"),
        )

    @override_settings(NICE_LITE_BATCH_PROJECT_INHERITANCE=True)
    def test_batch_job_project_source_resolves_owned_completed_project(self):
        job_dir = os.path.join(self.workspace_dir, "1_import_movies")
        os.mkdir(job_dir)
        project_path = os.path.join(job_dir, "workspace.simple")
        with open(project_path, "w", encoding="utf-8"):
            pass
        jobmodel = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=1,
            name="Import Movie Data",
            dirc="1_import_movies",
            status="finished",
            master_status="finished",
            master_stats={"job_type": "batch", "package": "simple", "program": "import_movies"},
        )

        resolved_path, metadata, error = job_builder_views._resolve_batch_project_source(
            self.workspace,
            f"job:{jobmodel.id}",
        )

        self.assertIsNone(error)
        self.assertEqual(resolved_path, project_path)
        self.assertEqual(metadata, {
            "type": "batch_job",
            "batch_job_id": jobmodel.id,
        })

    def test_batch_project_file_resolver_accepts_only_simple_files_in_workspace(self):
        project_path = os.path.join(self.workspace_dir, "selected.simple")
        with open(project_path, "w", encoding="utf-8"):
            pass
        text_path = os.path.join(self.workspace_dir, "selected.txt")
        with open(text_path, "w", encoding="utf-8"):
            pass
        outside_path = os.path.join(self.tempdir.name, "outside.simple")
        with open(outside_path, "w", encoding="utf-8"):
            pass

        resolved_path, metadata, error = job_builder_views._resolve_batch_project_file(
            self.workspace,
            project_path,
        )

        self.assertIsNone(error)
        self.assertEqual(resolved_path, project_path)
        self.assertEqual(metadata, {
            "type": "project_file",
            "filename": "selected.simple",
        })
        self.assertEqual(
            job_builder_views._resolve_batch_project_file(self.workspace, text_path),
            (None, None, "input project must be a .simple file"),
        )
        self.assertEqual(
            job_builder_views._resolve_batch_project_file(self.workspace, outside_path),
            (None, None, "batch project file is outside the selected workspace"),
        )

    def test_new_rejects_explicit_project_outside_workspace(self):
        outside_project = os.path.join(self.tempdir.name, "outside.simple")
        with open(outside_project, "w", encoding="utf-8"):
            pass

        launcher = batchjob_module.SIMPLEBatch
        with patch.object(launcher, "start") as start:
            created = BatchJob().new(
                self.workspace,
                "simple",
                "cluster2D",
                {},
                parent_proj=outside_project,
            )

        self.assertFalse(created)
        self.assertEqual(JobModel.objects.count(), 0)
        start.assert_not_called()

    def test_snapshot_sources_are_resolved_from_workspace_metadata(self):
        snapshot_dir = os.path.join(
            self.workspace_dir,
            "2_simple_stream",
            "classification_2D",
            "snapshots",
            "snapshot_3",
        )
        os.makedirs(snapshot_dir)
        snapshot_path = os.path.join(snapshot_dir, "snapshot_3.simple")
        with open(snapshot_path, "w", encoding="utf-8"):
            pass
        stream_job = JobModel.objects.create(
            dset=self.workspace_model,
            cdat=timezone.now(),
            disp=2,
            dirc="2_simple_stream",
            particle_sets_stats={
                "particle_sets": [{
                    "id": 3,
                    "name": "particle set 3",
                    "type": "snapshot",
                    "filename": "snapshot_3.simple",
                    "time": 123,
                }],
            },
        )

        sources = job_builder_views._collect_batch_snapshot_sources(self.workspace)
        project_path, metadata, error = job_builder_views._resolve_batch_project_source(
            self.workspace,
            f"snapshot:{stream_job.id}:3",
        )

        self.assertEqual(sources, [{
            "key": f"snapshot:{stream_job.id}:3",
            "label": "stream 2 - particle set 3",
        }])
        self.assertIsNone(error)
        self.assertEqual(project_path, snapshot_path)
        self.assertEqual(metadata, {
            "type": "stream_snapshot",
            "stream_job_id": stream_job.id,
            "particle_set_id": 3,
            "filename": "snapshot_3.simple",
        })


class SimpleBatchDispatchTests(TestCase):
    def test_lsf_submission_waits_for_scheduler_acceptance(self):
        dispatch = type("Dispatch", (), {"scmd": "bsub"})()

        with tempfile.TemporaryDirectory() as job_dir:
            script_path = os.path.join(job_dir, "job.script")
            with open(script_path, "w", encoding="utf-8"):
                pass
            with patch.object(simple_module.subprocess, "run") as run:
                simple_module._submit(dispatch, script_path, job_dir)

        run.assert_called_once_with(
            ["bsub", "-L", "/bin/sh"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=job_dir,
            stdin=ANY,
            text=True,
            check=True,
        )
        self.assertTrue(run.call_args.kwargs["stdin"].closed)

    def test_start_fails_when_lsf_rejects_submission(self):
        with tempfile.TemporaryDirectory() as parent_dir:
            parent_proj = os.path.join(parent_dir, "workspace.simple")
            with open(parent_proj, "w", encoding="utf-8"):
                pass
            base_dir = os.path.join(parent_dir, "import_movies")
            os.mkdir(base_dir)
            dispatch = type("Dispatch", (), {
                "tplt": "#!/bin/sh\nXXXSIMPLEXXX",
                "scmd": "bsub",
                "url": "localhost:8000",
            })()

            with (
                patch.object(SIMPLEBatch, "loadUIJSON", return_value=True),
                patch.object(simple_module.DispatchModel.objects, "filter") as dispatch_filter,
                patch.object(simple_module.shutil, "which", return_value="/usr/bin/bsub"),
                patch.object(
                    simple_module.subprocess,
                    "run",
                    side_effect=subprocess.CalledProcessError(1, ["bsub"]),
                ),
            ):
                dispatch_filter.return_value.last.return_value = dispatch
                started = SIMPLEBatch(pckg="simple").start(
                    {},
                    base_dir,
                    parent_dir,
                    "import_movies",
                    9,
                )

        self.assertFalse(started)

    @override_settings(NICE_LITE_BATCH_STATUS_CALLBACKS=False)
    def test_start_dispatches_corresponding_executable_and_quotes_values(self):
        with tempfile.TemporaryDirectory() as parent_dir:
            parent_proj = os.path.join(parent_dir, "parent project.simple")
            with open(parent_proj, "w", encoding="utf-8"):
                pass

            dispatch = type("Dispatch", (), {
                "tplt": "#!/bin/sh\n# CPU XXXNCPUXXX\nXXXSIMPLEXXX",
                "scmd": "sh",
                "url": "http://localhost:8000",
            })()

            for package, executable in (("simple", "simple_exec"), ("single", "single_exec")):
                with self.subTest(package=package):
                    base_dir = os.path.join(parent_dir, package)
                    os.mkdir(base_dir)
                    with patch.object(SIMPLEBatch, "loadUIJSON", return_value=True), patch.object(simple_module.DispatchModel.objects, "filter") as dispatch_filter, patch.object(simple_module.shutil, "which", return_value="/bin/sh"), patch.object(simple_module, "_submit") as submit:
                        dispatch_filter.return_value.last.return_value = dispatch
                        launcher = SIMPLEBatch(pckg=package)
                        started = launcher.start(
                            {"input": "path with spaces", "nthr": "8"},
                            base_dir,
                            parent_dir,
                            "demo_commander",
                            9,
                            parent_proj=parent_proj,
                        )

                    self.assertTrue(started)
                    with open(os.path.join(base_dir, "job.script"), encoding="utf-8") as script:
                        content = script.read()
                    self.assertIn(f"{executable} prg=demo_commander input='path with spaces' nthr=8", content)
                    self.assertIn(f"cp -v '{parent_proj}' workspace.simple", content)
                    self.assertIn("# CPU 8", content)
                    self.assertNotIn("nice_status_callback", content)
                    submit.assert_called_once()

    @override_settings(NICE_LITE_BATCH_STATUS_CALLBACKS=True)
    def test_status_callbacks_wrap_every_batch_package(self):
        with tempfile.TemporaryDirectory() as parent_dir:
            parent_proj = os.path.join(parent_dir, "workspace.simple")
            with open(parent_proj, "w", encoding="utf-8"):
                pass

            dispatch = type("Dispatch", (), {
                "tplt": "#!/bin/sh\nXXXNCPUXXX\nXXXSIMPLEXXX",
                "scmd": "sh",
                "url": "http://localhost:8000",
            })()

            for package, executable in (("simple", "simple_exec"), ("single", "single_exec")):
                with self.subTest(package=package):
                    base_dir = os.path.join(parent_dir, package)
                    os.mkdir(base_dir)
                    with (
                        patch.object(SIMPLEBatch, "loadUIJSON", return_value=True),
                        patch.object(simple_module.DispatchModel.objects, "filter") as dispatch_filter,
                        patch.object(simple_module.shutil, "which", return_value="/bin/sh"),
                        patch.object(simple_module, "_submit") as submit,
                    ):
                        dispatch_filter.return_value.last.return_value = dispatch
                        started = SIMPLEBatch(pckg=package).start(
                            {},
                            base_dir,
                            parent_dir,
                            "demo_commander",
                            9,
                            parent_proj=parent_proj,
                        )

                    self.assertTrue(started)
                    with open(os.path.join(base_dir, "job.script"), encoding="utf-8") as script:
                        content = script.read()

                    running = '{"jobid":9,"job_heartbeat":{}}'
                    finished = '{"jobid":9,"job":{"status":"finished","terminate":true}}'
                    failed = '{"jobid":9,"job":{"status":"failed","terminate":true}}'
                    command = f"{executable} prg=demo_commander"

                    self.assertIn("nice_status_callback()", content)
                    self.assertIn("X-Worker-Token: ${NICE_LITE_WORKER_TOKEN}", content)
                    self.assertIn("2>> nice_status.log", content)
                    self.assertIn('exit "$nice_job_exit"', content)
                    self.assertLess(content.index(running), content.index(command))
                    self.assertLess(content.index(command), content.index("nice_job_exit=$?"))
                    self.assertLess(content.index("nice_job_exit=$?"), content.index(finished))
                    self.assertLess(content.index("nice_job_exit=$?"), content.index(failed))
                    submit.assert_called_once()

    @override_settings(NICE_LITE_BATCH_STATUS_CALLBACKS=True)
    def test_new_project_dispatch_creates_project_in_job_directory(self):
        with tempfile.TemporaryDirectory() as parent_dir:
            parent_proj = os.path.join(parent_dir, "workspace.simple")
            with open(parent_proj, "w", encoding="utf-8"):
                pass
            base_dir = os.path.join(parent_dir, "new_project")
            os.mkdir(base_dir)
            dispatch = type("Dispatch", (), {
                "tplt": "#!/bin/sh\nXXXSIMPLEXXX",
                "scmd": "sh",
                "url": "http://localhost:8000",
            })()

            with (
                patch.object(SIMPLEBatch, "loadUIJSON", return_value=True),
                patch.object(simple_module.DispatchModel.objects, "filter") as dispatch_filter,
                patch.object(simple_module.shutil, "which", return_value="/bin/sh"),
                patch.object(simple_module, "_submit") as submit,
            ):
                dispatch_filter.return_value.last.return_value = dispatch
                started = SIMPLEBatch(pckg="simple").start(
                    {"projname": "new"},
                    base_dir,
                    parent_dir,
                    "new_project",
                    3,
                )

            self.assertTrue(started)
            with open(os.path.join(base_dir, "job.script"), encoding="utf-8") as script:
                content = script.read()
            self.assertIn("simple_exec prg=new_project projname=new dir=. mkdir=no", content)
            self.assertNotIn("cp -v", content)
            self.assertNotIn("prg=update_project", content)
            self.assertNotIn("projfile=workspace.simple", content)
            self.assertIn("nice_status_callback()", content)
            self.assertIn('{"jobid":3,"job":{"status":"finished","terminate":true}}', content)
            submit.assert_called_once()
