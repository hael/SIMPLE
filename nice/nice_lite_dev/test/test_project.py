import tempfile

from django.test  import TestCase
from django.http  import HttpRequest
from django.utils import timezone

from ..data_structures.project    import Project
from ..models                     import ProjectModel
from ..models                     import DatasetModel
from ..models                     import WorkspaceModel
from .test_helpers                import assertProject

class ProjectTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def setUp(self):
    # test projectmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    DatasetModel.objects.create(proj=project)
    WorkspaceModel.objects.create(proj=project)

  def test_project_init(self):
    # test init empty project
    project = Project()
    assertProject(project)

  def test_project_init_by_request(self):
    # test projectmodel object retrieval using request
    request = HttpRequest()
    request.POST["selected_project_id"] = "1"
    project = Project(request=request)
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)
  
  def test_project_init_by_id(self):
    # test projectmodel object retrieval using id
    project = Project(project_id=1)
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)

  def test_project_new(self):
    # test init empty project
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)

  def test_project_trash(self):
    # test create project trash folder
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      self.assertEqual(project.ensureTrashfolder(), True,  "failed to generate trash folder")

  def test_project_contains_dataset(self):
    # test project contains dataset returns true when exists
    project = Project(project_id=1)
    self.assertEqual(project.containsDataset(1), True, "project doesn't contain dataset")
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)
  
  def test_project_contains_workspace(self):
    # test project contains workspace returns true when exists
    project = Project(project_id=1)
    self.assertEqual(project.containsWorkspace(1), True, "project doesn't contain workspace") 
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)

  def test_project_doesnt_contain_dataset(self):
    # test project contains dataset returns false when not exists
    project = Project(project_id=1)
    self.assertEqual(project.containsDataset(2), False, "project contains dataset")
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)
  
  def test_project_doesnt_contain_workspace(self):
    # test project contains workspace returns false when not exists
    project = Project(project_id=1)
    self.assertEqual(project.containsWorkspace(2), False, "project contains workspace") 
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)
    