# global imports
import os
from django.utils import timezone

# local imports
from ..models import ProjectModel, DatasetModel, WorkspaceModel

class Project:

    id       = 0
    name     = ""
    desc     = ""
    dirc     = ""
    datasets_list   = []
    workspaces_list = []
    trashfolder = ""

    def __init__(self, project_id=None, request=None):
        if project_id is not None:
            self.id = project_id
        elif request is not None:
            self.setIDFromRequest(request)
        if self.id > 0:
            self.load()

    def setIDFromRequest(self, request):
        if "selected_project_id" in request.POST:
            test_id_str = request.POST["selected_project_id"]
        else:
            test_id_str = request.COOKIES.get('selected_project_id', 'none')
        if test_id_str.isnumeric():
            self.id = int(test_id_str)

    def load(self):
        projectmodel = ProjectModel.objects.filter(id=self.id).first()
        if projectmodel is not None:
            self.name = projectmodel.name
            self.desc = projectmodel.desc
            self.dirc = projectmodel.dirc
            datasetmodels = DatasetModel.objects.filter(proj=self.id)
            self.datasets_list = []
            for datasetmodel in datasetmodels:
                self.datasets_list.append({"id":datasetmodel.id, "name":datasetmodel.name})
            workspacemodels = WorkspaceModel.objects.filter(proj=self.id)
            self.workspaces_list = []
            for workspacemodel in workspacemodels:
                self.workspaces_list.append({"id":workspacemodel.id, "name":workspacemodel.name})
        else:
            self.id = 0

    def new(self, request):

        if "new_project_name" in request.POST:
            new_project_name = request.POST["new_project_name"]
        else:
            return False

        if "new_project_dirc" in request.POST:
            new_project_dirc = request.POST["new_project_dirc"]
        else:
            return False

        if not os.path.exists(new_project_dirc):
            return False
        
        if not os.path.isdir(new_project_dirc):
            return False
        
        new_project_path = os.path.join(new_project_dirc, new_project_name.replace(" ", "_"))

        try:
            os.makedirs(new_project_path, exist_ok=True)
        except OSError as error:
            print("Directory '%s' can not be created")
            return False
        
        projectmodel = ProjectModel(name=new_project_name, dirc=new_project_path, date=timezone.now())
        projectmodel.save()
        self.id = projectmodel.id
        if(self.id == 0):
            return False
        self.load()
        return True

    def containsDataset(self, dataset_id):
        for dataset in self.datasets_list:
            if dataset["id"] == dataset_id:
                return True
        return False
    
    def containsWorkspace(self, workspace_id):
        for workspace in self.workspaces_list:
            if workspace["id"] == workspace_id:
                return True
        return False
    
    def ensureTrashfolder(self):
        if not os.path.exists(self.dirc):
            return False
        if not os.path.isdir(self.dirc):
            return False
        trashfolder = os.path.join(self.dirc, "TRASH")
        if os.path.isdir(trashfolder):
            self.trashfolder = trashfolder
            return True
        try:
            os.makedirs(trashfolder, exist_ok=True)
        except OSError as error:
            print("Directory '%s' can not be created")
            return False
        self.trashfolder = trashfolder
        return True
        