from django.test  import TestCase
from django.utils import timezone

from ..models      import ProjectModel
from ..models      import DatasetModel
from ..models      import WorkspaceModel
from ..models      import JobModel
from ..models      import JobClassicModel
from ..models      import DispatchModel
from .test_helpers import assertProject
from .test_helpers import assertWorkspace
from .test_helpers import assertDataset
from .test_helpers import assertDispatch

class ProjectModelTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def setUp(self):
    # test projectmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())

  def test_project_lookup_by_id(self):
    # test projectmodel object retrieval using id
    project = ProjectModel.objects.get(id=1)
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)


class DatasetModelTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def setUp(self):
    # test datasetmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    DatasetModel.objects.create(proj=project)

  def test_dataset_lookup_by_id(self):
    # test datasetmodel object retrieval using id
    dataset = DatasetModel.objects.get(id=1)
    assertDataset(dataset, id=1)
    project = dataset.proj
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)


class WorkspaceModelTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def setUp(self):
    # test workspacemodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    WorkspaceModel.objects.create(proj=project)

  def test_workspace_lookup_by_id(self):
    # test workspacemodel object retrieval using id
    workspace = WorkspaceModel.objects.get(id=1)
    assertWorkspace(workspace, id=1)
    project   = workspace.proj
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)


class JobModelTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def setUp(self):
    # test jobmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    DatasetModel.objects.create(proj=project)
    dataset = DatasetModel.objects.get(id=1)
    JobModel.objects.create(dset=dataset, cdat=timezone.now())

  def test_job_lookup_by_id(self):
    # test jobmodel object retrieval using id
    job     = JobModel.objects.get(id=1)
    dataset = job.dset
    assertDataset(dataset, id=1)
    project = dataset.proj
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)


class JobClassicModelTest(TestCase):

  test_project_name = "testproject"
  test_project_desc = "test project description"
  test_project_dirc = "/tmp"

  def setUp(self):
    # test jobmodel object creation with non-optional fields
    ProjectModel.objects.create(name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc, date=timezone.now())
    project = ProjectModel.objects.get(id=1)
    WorkspaceModel.objects.create(proj=project)
    workspace = WorkspaceModel.objects.get(id=1)
    JobClassicModel.objects.create(wspc=workspace, cdat=timezone.now())

  def test_job_classic_lookup_by_id(self):
    # test jobmodelclassic object retrieval using id
    jobclassic = JobClassicModel.objects.get(id=1)
    workspace  = jobclassic.wspc
    project    = workspace.proj
    assertWorkspace(workspace, id=1)
    assertProject(project, name=self.test_project_name, desc=self.test_project_desc, dirc=self.test_project_dirc)
    

class DispatchModelTest(TestCase):

  def setUp(self):
    # test dispatchmodel object creation with non-optional fields
    DispatchModel.objects.create()

  def test_dispatch_lookup_by_id(self):
    # test jobmodelclassic object retrieval using id
    dispatch = DispatchModel.objects.get(id=1)
    assertDispatch(dispatch, id=1)
