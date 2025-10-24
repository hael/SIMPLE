# global imports
import os
import re
from django.utils import timezone

# local imports
from ..models import ProjectModel, DatasetModel, WorkspaceModel

class Dataset:

    id       = 0
    disp     = 0
    name     = ""
    desc     = ""
    dirc     = ""
    link     = ""
    cdat     = ""
    mdat     = ""
    user     = ""
    proj     = 0
    trashfolder = ""

    def __init__(self, dataset_id=None, request=None):
        if dataset_id is not None:
            self.id = dataset_id
        elif request is not None:
            self.setIDFromRequest(request)
        if self.id > 0:
            self.load()

    def setIDFromRequest(self, request):
        if "selected_dataset_id" in request.POST:
            test_id_str = request.POST["selected_dataset_id"]
        else:
            test_id_str = request.COOKIES.get('selected_dataset_id', 'none')
        if test_id_str.isnumeric():
            self.id = int(test_id_str)

    def load(self):
        datasetmodel = DatasetModel.objects.filter(id=self.id).first()
        if datasetmodel is not None:
            self.name = datasetmodel.name
            self.dirc = datasetmodel.dirc
            self.link = datasetmodel.link
            self.cdat = datasetmodel.cdat
            self.mdat = datasetmodel.mdat
            self.proj = datasetmodel.proj
            self.desc = datasetmodel.desc
            self.disp = datasetmodel.disp
            self.user = datasetmodel.user

    def new(self, project, user):
        projectmodel = ProjectModel.objects.filter(id=project.id).first()
        if projectmodel is None:
            return False
      
        if not os.path.exists(project.dirc):
            return False
        
        if not os.path.isdir(project.dirc):
            return False
        
        datasetmodels = DatasetModel.objects.filter(proj=projectmodel)
        self.disp = datasetmodels.count() + 1
        new_dataset_name = "new dataset " + str(self.disp)
        datasetmodel = DatasetModel(proj=projectmodel, name=new_dataset_name, disp=self.disp, cdat=timezone.now(), mdat=timezone.now(), user=user)
        datasetmodel.save()
        self.id = datasetmodel.id
        if self.id == 0:
            return False
        
        new_dataset_dirc = ".dataset_" + str(self.id)
        new_dataset_link = "ds_" + new_dataset_name.replace(" ", "_")
        new_dataset_path = os.path.join(project.dirc, new_dataset_dirc)

        try:
            os.makedirs(new_dataset_path)
        except OSError as error:
            print("Directory '%s' can not be created")
            return False
        
        try:
            os.symlink(new_dataset_path, os.path.join(project.dirc, new_dataset_link))
        except FileExistsError:
            print("Symlink already exists.")
        except PermissionError:
            print("Permission denied: You might need admin rights.")
        except OSError as e:
            print("OS error occurred:", e)
                
        datasetmodel.dirc = new_dataset_dirc
        datasetmodel.link = new_dataset_link
        datasetmodel.save()
        self.load()
        project.load()
        return True#

    def delete(self, project):
        datasetmodel = DatasetModel.objects.filter(id=self.id).first()
        if project.ensureTrashfolder() and datasetmodel is not None:
            dataset_path = os.path.join(project.dirc, self.dirc)
            dataset_link = os.path.join(project.dirc, self.link)
            trash_path   = os.path.join(project.trashfolder, self.dirc)
            trash_link   = os.path.join(project.trashfolder, self.link)
            if not os.path.exists(dataset_path):
                return 
            if not os.path.isdir(dataset_path):
                return 
            if not os.path.exists(dataset_link):
                return 
            if not os.path.islink(dataset_link):
                return 
            if os.path.exists(trash_path):
                return 
            if os.path.isdir(trash_path):
                return
            if os.path.exists(trash_link):
                return 
            if os.path.islink(trash_link):
                return
            try:
                os.rename(dataset_path, trash_path)
            except OSError as error:
                print("Directory '%s' can not be renamed")
                return
            try:
                os.remove(dataset_link)
                os.symlink(trash_path, trash_link)
            except OSError as error:
                print("Link '%s' can not be renamed")
                return
            datasetmodel.delete()
        workspacemodels = WorkspaceModel.objects.filter(proj=project.id)
        datasetmodels   = DatasetModel.objects.filter(proj=project.id)
        if workspacemodels.count() + datasetmodels.count() == 0:
            projectmodel = ProjectModel.objects.filter(id=project.id).first()
            projectmodel.delete()

    def rename(self, request, project):
        if "new_dataset_name" in request.POST:
            new_dataset_name = request.POST["new_dataset_name"]
            datasetmodel = DatasetModel.objects.filter(id=self.id).first()
            datasetmodel.name = new_dataset_name
            # strip any non-alphanumeric characters(except _)
            new_dataset_name = new_dataset_name.replace(" ", "_")
            new_dataset_name = "ds_" + re.sub(r'\W+', '', new_dataset_name)
            current_dataset_link = os.path.join(project.dirc, self.link)
            new_dataset_link     = os.path.join(project.dirc, new_dataset_name)
            try :
                os.rename(current_dataset_link, new_dataset_link)
                print("Source path renamed to destination path successfully.")
            except IsADirectoryError:
                print("Source is a file but destination is a directory.")
                return False
            except NotADirectoryError:
                print("Source is a directory but destination is a file.")

            except PermissionError:
                print("Operation not permitted")
                return False
            except OSError as error:
                print(error)
                return False
            self.link = new_dataset_name
            datasetmodel.link = new_dataset_name
            datasetmodel.save()
            return True
        return False

    def updateDescription(self, request):
        if "new_dataset_description" in request.POST:
            new_dataset_description = request.POST["new_dataset_description"]
            datasetmodel = DatasetModel.objects.filter(id=self.id).first()
            datasetmodel.desc = new_dataset_description
            datasetmodel.save()
            return True
        return False
    
    def ensureTrashfolder(self, project):
        if not os.path.exists(os.path.join(project.dirc, self.dirc)):
            return False
        if not os.path.isdir(os.path.join(project.dirc, self.dirc)):
            return False
        trashfolder = os.path.join(project.dirc, self.dirc, "TRASH")
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


