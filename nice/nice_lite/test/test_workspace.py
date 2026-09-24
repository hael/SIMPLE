import os
import tempfile
from unittest.mock import patch

from django.http  import HttpRequest
from django.test  import TestCase
from django.utils import timezone

from ..data_structures import simple as simple_module
from ..data_structures.project   import Project
from ..data_structures.workspace import Workspace
from ..models                    import JobModel
from ..models                    import ProjectModel
from ..models                    import WorkspaceModel
from .test_helpers               import assertProject
from .test_helpers               import assertWorkspace

class WorkspaceTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"
  test_workspace_user = "testuser"
  test_workspace_name = "testworkspace"
  test_workspace_desc = "test workspace description"

  def setUp(self):
    simple_exec_patch = patch.object(
      simple_module.subprocess,
      "run",
      side_effect=self._create_workspace_project,
    )
    self.simple_exec = simple_exec_patch.start()
    self.addCleanup(simple_exec_patch.stop)

    project = ProjectModel.objects.create(
      name=self.test_project_name,
      desc=self.test_project_desc,
      dirc=self.test_project_dirc,
      date=timezone.now(),
    )
    WorkspaceModel.objects.create(proj=project)

  @staticmethod
  def _create_workspace_project(command, check):
    workspace_dir = next(arg.removeprefix("dir=") for arg in command if arg.startswith("dir="))
    with open(os.path.join(workspace_dir, "workspace.simple"), "w", encoding="utf-8"):
      pass

  def test_workspace_init(self):
    # test init empty workspace
    workspace = Workspace()
    assertWorkspace(workspace)

  def test_workspace_init_by_id(self):
    # test workspacemodel object retrieval using id
    workspace = Workspace(id=1)
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
      self.assertTrue(workspace.new(project, user=self.test_workspace_user))
      assertWorkspace(workspace, id=2, user=self.test_workspace_user)
      workspace_path = os.path.join(project.get_absdir(), ".workspace_2")
      self.simple_exec.assert_called_once_with(
        ["simple_exec", "prg=new_project", "projname=workspace", "dir=" + workspace_path],
        check=True,
      )
      self.assertTrue(os.path.isfile(os.path.join(workspace_path, "workspace.simple")))

  def test_workspace_new_fails_when_project_file_creation_fails(self):
    self.simple_exec.side_effect = None
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      workspace = Workspace()

      self.assertFalse(workspace.new(project, user=self.test_workspace_user))
      self.assertFalse(WorkspaceModel.objects.filter(id=2).exists())
      self.assertFalse(os.path.lexists(os.path.join(project.get_absdir(), "ds_new_workspace_1")))

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
      self.assertTrue(workspace.delete(), "failed to delete workspace")
      workspace = Workspace(id=2)
      self.assertEqual(workspace.id, 0, "deleted workspace is still present")

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
      self.assertEqual(workspace.rename(self.test_workspace_name), True, "failed to rename workspace")
      workspace = Workspace(id=2)
      assertWorkspace(workspace, id=2, name=self.test_workspace_name)

  def test_workspace_update_description(self):
    # test workspacemodel object retrieval using id
    request = HttpRequest()
    request.POST["new_workspace_description"] = self.test_workspace_desc
    workspace = Workspace(id=1)
    assertWorkspace(workspace, id=1)
    self.assertEqual(
      workspace.updateDescription(request.POST["new_workspace_description"]),
      True,
      "failed to update workspace description",
    )
    workspace = Workspace(id=1)
    assertWorkspace(workspace, id=1, desc=self.test_workspace_desc)

  def test_workspace_reconciles_job_counter_to_highest_surviving_job(self):
    workspace_model = WorkspaceModel.objects.get(id=1)
    workspace_model.jcnt = 7
    workspace_model.save(update_fields=["jcnt"])
    JobModel.objects.create(
      dset=workspace_model,
      cdat=timezone.now(),
      disp=2,
    )
    JobModel.objects.create(
      dset=workspace_model,
      cdat=timezone.now(),
      disp=5,
    )

    workspace = Workspace(id=workspace_model.id)
    self.assertTrue(workspace.reconcile_job_counter())
    workspace_model.refresh_from_db()
    self.assertEqual(workspace_model.jcnt, 5)

    JobModel.objects.filter(dset=workspace_model, disp=5).delete()
    self.assertTrue(workspace.reconcile_job_counter())
    workspace_model.refresh_from_db()
    self.assertEqual(workspace_model.jcnt, 2)

    JobModel.objects.filter(dset=workspace_model).delete()
    self.assertTrue(workspace.reconcile_job_counter())
    workspace_model.refresh_from_db()
    self.assertEqual(workspace_model.jcnt, 0)
