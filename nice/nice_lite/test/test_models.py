"""Model smoke tests for core NICE Lite ORM relationships."""

from django.test import TestCase
from django.utils import timezone

from ..models import DispatchModel
from ..models import JobModel
from ..models import ProjectModel
from ..models import WorkspaceModel


class _ProjectFixtureMixin:
    test_project_name = "testproject"
    test_project_desc = "test project description"
    test_project_dirc = "/tmp"

    def _create_project(self):
        return ProjectModel.objects.create(
            name=self.test_project_name,
            desc=self.test_project_desc,
            dirc=self.test_project_dirc,
            date=timezone.now(),
        )

    def _assert_project_fields(self, project):
        self.assertEqual(project.name, self.test_project_name)
        self.assertEqual(project.desc, self.test_project_desc)
        self.assertEqual(project.dirc, self.test_project_dirc)


class ProjectModelTest(_ProjectFixtureMixin, TestCase):
    def setUp(self):
        self._create_project()

    def test_project_lookup_by_id(self):
        project = ProjectModel.objects.get(id=1)
        self._assert_project_fields(project)


class WorkspaceModelTest(_ProjectFixtureMixin, TestCase):
    def setUp(self):
        project = self._create_project()
        WorkspaceModel.objects.create(proj=project)

    def test_workspace_lookup_by_id(self):
        workspace = WorkspaceModel.objects.get(id=1)
        self.assertEqual(workspace.id, 1)
        self._assert_project_fields(workspace.proj)


class JobModelTest(_ProjectFixtureMixin, TestCase):
    def setUp(self):
        project = self._create_project()
        workspace = WorkspaceModel.objects.create(proj=project)
        JobModel.objects.create(dset=workspace, cdat=timezone.now())

    def test_job_lookup_by_id(self):
        job = JobModel.objects.get(id=1)
        self.assertEqual(job.dset.id, 1)
        self._assert_project_fields(job.dset.proj)


class DispatchModelTest(TestCase):
    def setUp(self):
        DispatchModel.objects.create()

    def test_dispatch_lookup_by_id(self):
        dispatch = DispatchModel.objects.get(id=1)
        self.assertEqual(dispatch.id, 1)
