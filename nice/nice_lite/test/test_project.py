import tempfile

from django.http  import HttpRequest
from django.test  import TestCase
from django.utils import timezone

from ..data_structures.project import Project
from ..models                  import ProjectModel
from ..models                  import WorkspaceModel


class ProjectTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def _assert_project_fields(self, project):
    self.assertEqual(project.name, self.test_project_name)
    self.assertEqual(project.desc, self.test_project_desc)
    self.assertEqual(project.dirc, self.test_project_dirc)


  def setUp(self):
    project = ProjectModel.objects.create(
      name=self.test_project_name,
      desc=self.test_project_desc,
      dirc=self.test_project_dirc,
      date=timezone.now(),
    )
    WorkspaceModel.objects.create(proj=project)
    WorkspaceModel.objects.create(proj=project)


  def test_project_init(self):
    project = Project()
    self.assertEqual(project.id, 0)


  def test_project_init_by_request(self):
    request = HttpRequest()
    request.POST["selected_project_id"] = "1"
    project = Project(request=request)
    self._assert_project_fields(project)
  

  def test_project_init_by_id(self):
    project = Project(project_id=1)
    self._assert_project_fields(project)


  def test_project_new(self):
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      self.assertEqual(project.id, 2)
      self.assertEqual(project.name, self.test_project_name)


  def test_project_trash(self):
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      self.assertEqual(project.id, 2)
      self.assertTrue(project.ensureTrashfolder(), "failed to generate trash folder")


  def test_project_workspaces_list_contains_existing_workspace(self):
    project = Project(project_id=1)
    workspace_ids = [workspace.id for workspace in project.workspaces_list]
    self.assertIn(1, workspace_ids, "project workspaces_list missing workspace 1")
    self._assert_project_fields(project)
  

  def test_project_workspaces_list_excludes_missing_workspace(self):
    project = Project(project_id=1)
    workspace_ids = [workspace.id for workspace in project.workspaces_list]
    self.assertNotIn(3, workspace_ids, "project workspaces_list unexpectedly contains workspace 3")
    self._assert_project_fields(project)
    