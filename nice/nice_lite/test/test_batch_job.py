import os
import signal
import subprocess
import tempfile
from unittest.mock import ANY, patch

from django.test import TestCase
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
                    submit.assert_called_once()

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
            submit.assert_called_once()
