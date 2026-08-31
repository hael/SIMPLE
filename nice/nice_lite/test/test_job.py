import os
import tempfile

from django.test import TestCase
from django.utils import timezone

from ..data_structures.job import Job
from ..models import JobModel, ProjectModel, WorkspaceModel


class JobLifecycleTests(TestCase):
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
        self.workspace = WorkspaceModel.objects.create(
            proj=self.project,
            dirc=".workspace_1",
            user="tester",
            jcnt=1,
        )

    def test_delete_resets_counter_and_uses_unique_trash_name_when_number_is_reused(self):
        job_dir = os.path.join(self.workspace_dir, "1_simple_stream")
        os.mkdir(job_dir)
        trash_dir = os.path.join(self.workspace_dir, "TRASH")
        os.mkdir(trash_dir)
        os.mkdir(os.path.join(trash_dir, "1_simple_stream"))
        jobmodel = JobModel.objects.create(
            dset=self.workspace,
            cdat=timezone.now(),
            disp=1,
            dirc="1_simple_stream",
        )
        job_id = jobmodel.id

        deleted = Job(id=job_id).delete()

        self.assertTrue(deleted)
        self.assertFalse(os.path.exists(job_dir))
        self.assertTrue(os.path.isdir(os.path.join(trash_dir, f"1_simple_stream_{job_id}")))
        self.assertFalse(JobModel.objects.filter(id=job_id).exists())
        self.workspace.refresh_from_db()
        self.assertEqual(self.workspace.jcnt, 0)
