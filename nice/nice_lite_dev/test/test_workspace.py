import tempfile

from django.test  import TestCase
from django.http  import HttpRequest
from django.utils import timezone

from ..data_structures.project    import Project
from ..data_structures.workspace  import Workspace
from ..models                     import ProjectModel
from ..models                     import DatasetModel
from ..models                     import WorkspaceModel
from .test_helpers                import assertProject
from .test_helpers                import assertWorkspace

class WorkspaceTest(TestCase):

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
    # test init empty workspace
    workspace = Workspace()
    assertWorkspace(workspace)

  def test_workspace_init_by_request(self):
    # test projectmodel object retrieval using request
    request = HttpRequest()
    request.POST["selected_workspace_id"] = "1"
    workspace = Workspace(request=request)
    assertWorkspace(workspace, id=1)
  
  def test_workspace_init_by_id(self):
    # test workspacemodel object retrieval using id
    workspace = Workspace(workspace_id=1)
    assertWorkspace(workspace, id=1)

  def test_workspace_new(self):
    # test init empty workspace
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspace(workspace, id=2, user=self.test_workspace_user)

  def test_workspace_trash(self):
    # test create workspace trash folder
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspace(workspace, id=2, user=self.test_workspace_user)
      self.assertEqual(workspace.ensureTrashfolder(project), True,  "failed to generate trash folder")

  def test_workspace_delete(self):
    # test delete workspace
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      workspace = Workspace()
      workspace.new(project, user=self.test_workspace_user)
      assertWorkspace(workspace, id=2, user=self.test_workspace_user)
      workspace.delete(project)
      workspace = Workspace(workspace_id=2)
      assertWorkspace(workspace, id=2)
      self.assertEqual(workspace.disp, 0, "deleted workspace os still present")

  def test_workspace_rename(self):
    # test renaming workspace
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
      assertWorkspace(workspace, id=2, user=self.test_workspace_user)
      self.assertEqual(workspace.rename(request, project), True,  "failed to rename workspace")
      workspace = Workspace(workspace_id=2)
      assertWorkspace(workspace, id=2, name=self.test_workspace_name)

  def test_workspace_update_description(self):
    # test workspacemodel object retrieval using id
    request = HttpRequest()
    request.POST["new_workspace_description"] = self.test_workspace_desc
    workspace = Workspace(workspace_id=1)
    assertWorkspace(workspace, id=1)
    self.assertEqual(workspace.updateDescription(request), True,  "failed to update workspace description")
    workspace = Workspace(workspace_id=1)
    assertWorkspace(workspace, id=1, desc=self.test_workspace_desc)
