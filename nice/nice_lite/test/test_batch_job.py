import os
import tempfile
from unittest.mock import patch

from django.test import TestCase
from django.utils import timezone

from ..data_structures import batchjob as batchjob_module
from ..data_structures import simple as simple_module
from ..data_structures.batchjob import BatchJob
from ..data_structures.simple import SIMPLEBatch
from ..data_structures.workspace import Workspace
from ..models import JobModel, ProjectModel, WorkspaceModel


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
            created = job.new(self.workspace, "single", "pick", {"mode": "fast"})

        self.assertTrue(created)
        jobmodel = JobModel.objects.get(id=job.id)
        self.assertEqual(jobmodel.args, {"mode": "fast"})
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
        self.assertEqual(BatchJob(id=job.id).prog, "pick")

    def test_new_marks_persisted_job_failed_when_dispatch_fails(self):
        launcher = batchjob_module.SIMPLEBatch
        with patch.object(launcher, "loadUIJSON", return_value=True), patch.object(launcher, "start", return_value=False):
            job = BatchJob()
            created = job.new(self.workspace, "simple", "import_movies", {})

        self.assertFalse(created)
        jobmodel = JobModel.objects.get(id=job.id)
        self.assertEqual(jobmodel.status, "failed")
        self.assertEqual(jobmodel.master_status, "failed")

    def test_new_initializes_missing_workspace_project_before_dispatch(self):
        os.remove(os.path.join(self.workspace_dir, "workspace.simple"))
        launcher = batchjob_module.SIMPLEBatch
        project_initializer = batchjob_module.SIMPLEProject

        with patch.object(project_initializer, "create", return_value=True) as create, patch.object(launcher, "loadUIJSON", return_value=True), patch.object(launcher, "start", return_value=True):
            created = BatchJob().new(self.workspace, "simple", "import_movies", {})

        self.assertTrue(created)
        create.assert_called_once_with()


class SimpleBatchDispatchTests(TestCase):
    def test_start_dispatches_corresponding_executable_and_quotes_values(self):
        with tempfile.TemporaryDirectory() as parent_dir:
            parent_proj = os.path.join(parent_dir, "parent project.simple")
            with open(parent_proj, "w", encoding="utf-8"):
                pass

            dispatch = type("Dispatch", (), {
                "tplt": "#!/bin/sh\nXXXSIMPLEXXX",
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
                            {"input": "path with spaces"},
                            base_dir,
                            parent_dir,
                            "demo_commander",
                            9,
                            parent_proj=parent_proj,
                        )

                    self.assertTrue(started)
                    with open(os.path.join(base_dir, "job.script"), encoding="utf-8") as script:
                        content = script.read()
                    self.assertIn(f"{executable} prg=demo_commander input='path with spaces'", content)
                    self.assertIn(f"cp -v '{parent_proj}' workspace.simple", content)
                    submit.assert_called_once()
