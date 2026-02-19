import tempfile

from django.test  import TestCase
from django.http  import HttpRequest
from django.utils import timezone

from ..data_structures.project    import Project
from ..data_structures.dataset    import Dataset
from ..models                     import ProjectModel
from ..models                     import DatasetModel
from ..models                     import WorkspaceModel
from .test_helpers                import assertProject
from .test_helpers                import assertDataset

class DatasetTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"
  test_dataset_user = "testuser"
  test_dataset_name = "testdataset"
  test_dataset_desc = "test dataset description"

  def setUp(self):
    # test projectmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    DatasetModel.objects.create(proj=project)

  def test_dataset_init(self):
    # test init empty dataset
    dataset = Dataset()
    assertDataset(dataset)

  def test_dataset_init_by_request(self):
    # test projectmodel object retrieval using request
    request = HttpRequest()
    request.POST["selected_dataset_id"] = "1"
    dataset = Dataset(request=request)
    assertDataset(dataset, id=1)
  
  def test_dataset_init_by_id(self):
    # test datasetmodel object retrieval using id
    dataset = Dataset(dataset_id=1)
    assertDataset(dataset, id=1)

  def test_dataset_new(self):
    # test init empty dataset
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      dataset = Dataset()
      dataset.new(project, user=self.test_dataset_user)
      assertDataset(dataset, id=2, user=self.test_dataset_user)

  def test_dataset_trash(self):
    # test create dataset trash folder
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      dataset = Dataset()
      dataset.new(project, user=self.test_dataset_user)
      assertDataset(dataset, id=2, user=self.test_dataset_user)
      self.assertEqual(dataset.test_dataset_trash(project), True,  "failed to generate trash folder")

  def test_dataset_delete(self):
    # test delete dataset
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      dataset = Dataset()
      dataset.new(project, user=self.test_dataset_user)
      assertDataset(dataset, id=2, user=self.test_dataset_user)
      dataset.delete(project)
      dataset = Dataset(dataset_id=2)
      self.assertEqual(dataset.proj, 0, "deleted dataset os still present")

  def test_dataset_rename(self):
    # test renaming dataset
    with tempfile.TemporaryDirectory() as tmpdirc:
      request = HttpRequest()
      request.POST["new_project_name"] = self.test_project_name
      request.POST["new_dataset_name"] = self.test_dataset_name
      request.POST["new_project_dirc"] = tmpdirc
      project = Project()
      project.new(request)
      assertProject(project, name=self.test_project_name, id=2)
      dataset = Dataset()
      dataset.new(project, user=self.test_dataset_user)
      assertDataset(dataset, id=2, user=self.test_dataset_user)
      self.assertEqual(dataset.rename(request, project), True,  "failed to rename dataset")
      dataset = Dataset(dataset_id=2)
      assertDataset(dataset, id=2, name=self.test_dataset_name)

  def test_dataset_update_description(self):
    # test datasetmodel object retrieval using id
    request = HttpRequest()
    request.POST["new_dataset_description"] = self.test_dataset_desc
    dataset = Dataset(dataset_id=1)
    assertDataset(dataset, id=1)
    self.assertEqual(dataset.updateDescription(request), True,  "failed to update dataset description")
    dataset = Dataset(dataset_id=1)
    assertDataset(dataset, id=1, desc=self.test_dataset_desc)
