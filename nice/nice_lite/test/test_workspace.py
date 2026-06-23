import tempfile

from django.test  import TestCase
from django.http  import HttpRequest
from django.utils import timezone

from ..data_structures.project    import Project
from ..data_structures.workspace    import Workspace
from ..models                     import ProjectModel
from ..models                     import WorkspaceModel
from .test_helpers                import assertProject
from .test_helpers                import assertWorkspaceAlias

class WorkspaceAliasTest(TestCase):

  test_project_name   = "testproject"
  test_project_desc   = "test project description"
  test_project_dirc   = "/tmp"
  test_workspace_user = "testuser"
  test_workspace_name = "testworkspace"
  test_workspace_desc = "test workspace description"

  def setUp(self):
    # test projectmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    WorkspaceModel.objects.create(proj=project)

  def test_workspace_init(self):
    workspace = Workspace()
    assertWorkspaceAlias(workspace)

  def test_workspace_init_by_request(self):
    request = HttpRequest()
    request.POST["selected_workspace_id"] = "1"
    workspace = Workspace(1)
    assertWorkspaceAlias(workspace, id=1)
  
  def test_workspace_init_by_id(self):
    workspace = Workspace(1)
    assertWorkspaceAlias(workspace, id=1)

  def test_workspace_new(self):
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspaceAlias(workspace, id=2, user=self.test_workspace_user)

  def test_workspace_trash(self):
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspaceAlias(workspace, id=2, user=self.test_workspace_user)
      self.assertEqual(workspace.get_trashdir() is not None, True, "failed to generate trash folder")

  def test_workspace_delete(self):
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspaceAlias(workspace, id=2, user=self.test_workspace_user)
      workspace.delete()
      workspace = Workspace(2)
      assertWorkspaceAlias(workspace, id=0)

  def test_workspace_rename(self):
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_workspace_name"] = self.test_workspace_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspaceAlias(workspace, id=2, user=self.test_workspace_user)
      self.assertEqual(workspace.rename(self.test_workspace_name), True,  "failed to rename workspace")
      workspace = Workspace(2)
      assertWorkspaceAlias(workspace, id=2)

  def test_workspace_update_description(self):
    workspace = Workspace(1)
    assertWorkspaceAlias(workspace, id=1)
    self.assertEqual(workspace.updateDescription(self.test_workspace_desc), True,  "failed to update workspace description")
    workspace = Workspace(1)
    assertWorkspaceAlias(workspace, id=1)
