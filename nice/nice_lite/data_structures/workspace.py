# global imports
import os
import re
import copy
from django.utils import timezone
from django.template.loader import render_to_string

# local imports
from ..models import ProjectModel, WorkspaceModel, JobClassicModel, DatasetModel
from .jobclassic import JobClassic
from .simple import SIMPLEProject

class Workspace:

    id       = 0
    disp     = 0
    name     = ""
    desc     = ""
    dirc     = ""
    link     = ""
    cdat     = ""
    mdat     = ""
    user     = ""
    nstr     = {}
    nstrhtml = {}
    request  = None
    trashfolder = ""

    def __init__(self, workspace_id=None, request=None):
        if workspace_id is not None:
            self.id = workspace_id
        elif request is not None:
            self.setIDFromRequest(request)
            self.request = request
        if self.id > 0:
            self.load()

    def setIDFromRequest(self, request):
        if "selected_workspace_id" in request.POST:
            test_id_str = request.POST["selected_workspace_id"]
        else:
            test_id_str = request.COOKIES.get('selected_workspace_id', 'none')
        if test_id_str.isnumeric():
            self.id = int(test_id_str)

    def load(self):
        workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
        if workspacemodel is not None:
            self.name     = workspacemodel.name
            self.dirc     = workspacemodel.dirc
            self.link     = workspacemodel.link
            self.cdat     = workspacemodel.cdat
            self.mdat     = workspacemodel.mdat
            self.nstr     = workspacemodel.nstr
            self.desc     = workspacemodel.desc
            self.disp     = workspacemodel.disp
            self.user     = workspacemodel.user
            self.nstrhtml = copy.deepcopy(self.nstr)
            self.addNodesHTML()

    def new(self, project, user):
        projectmodel = ProjectModel.objects.filter(id=project.id).first()
        if projectmodel is None:
            return False
      
        if not os.path.exists(project.dirc):
            return False
        
        if not os.path.isdir(project.dirc):
            return False
        
        workspacemodels = WorkspaceModel.objects.filter(proj=projectmodel)
        self.disp = workspacemodels.count() + 1
        new_workspace_name = "new workspace " + str(self.disp)
        workspacemodel = WorkspaceModel(proj=projectmodel, disp=self.disp, name=new_workspace_name, cdat=timezone.now(), mdat=timezone.now(), user=user)
        workspacemodel.save()
        self.id = workspacemodel.id
        if self.id == 0:
            return False
        
        new_workspace_dirc = ".workspace_" + str(self.id)
        new_workspace_link = "ws_" + new_workspace_name.replace(" ", "_")
        new_workspace_path = os.path.join(project.dirc, new_workspace_dirc)

        try:
            os.makedirs(new_workspace_path)
        except OSError as error:
            print("Directory '%s' can not be created")
            return False
        
        try:
            os.symlink(new_workspace_path, os.path.join(project.dirc, new_workspace_link))
        except FileExistsError:
            print("Symlink already exists.")
        except PermissionError:
            print("Permission denied: You might need admin rights.")
        except OSError as e:
            print("OS error occurred:", e)

        workspacemodel.dirc = new_workspace_dirc
        workspacemodel.link = new_workspace_link
        workspacemodel.nstr = {
            "children":[
                {
                    "type"      : "new",
                    "jobid"     : 0,
                    "innerHTML" : ""
                }
            ]
        }
        workspacemodel.save()
        self.load()
        project.load()
        simpleproject = SIMPLEProject(new_workspace_path)
        simpleproject.create()
        return True
    
    def addNodesHTML(self):

        def addNodeHTML(obj):
            if "children" in obj:
                for child in obj["children"]:
                    item = addNodeHTML(child)
            if "innerHTML" in obj and "type" in obj and "jobid" in obj:
                if "job" in obj["type"]:
                    jobclassic = JobClassic(id=obj["jobid"])
                    context = { 
                        'id'     : obj["jobid"],
                        'status' : jobclassic.status,
                        'disp'   : jobclassic.disp,
                        'name'   : jobclassic.name,
                        'desc'   : jobclassic.desc,
                        'pckg'   : jobclassic.pckg,
                    }
                    obj["innerHTML"] = render_to_string('nice_classic/jobnode.html', context, request=self.request).replace("\n", "").replace('"','\\"')
                elif "new" in obj["type"]:
                    context = { 
                        'id': obj["jobid"], 
                    }
                    obj["innerHTML"] = render_to_string('nice_classic/newjobnode.html', context, request=self.request).replace("\n", "").replace('"','\\"')

        if self.request != None:
            addNodeHTML(self.nstrhtml)
        
    def addChild(self, parentid, jobid):

        def findParent(obj):
            if "jobid" in obj:
                if parentid == obj["jobid"]:
                    return obj
            if "children" in obj:
                for child in obj["children"]:
                    item = findParent(child)
                    if item is not None:
                        return item
            return None

        workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
        if workspacemodel is not None:
            updated = workspacemodel.nstr
            newchild = {
                    "jobid"     : jobid,
                    "type"      : "job",
                    "innerHTML" : "",
                    "children"  : [
                        {
                            "type"      : "new",
                            "jobid"     : jobid,
                            "innerHTML" : ""
                        }
                    ]
            }
            if parentid == 0:
                updated["children"].insert(0,newchild)
            else:
                parent = findParent(updated)
                if parent is not None:
                    parent["children"].insert(0,newchild)
            workspacemodel.nstr = updated
            workspacemodel.save()
            self.nstr = updated

    def removeChild(self, jobid):

        def findChild(obj):
            if "jobid" in obj:
                if jobid == obj["jobid"]:
                    return obj
            if "children" in obj:
                for child in obj["children"]:
                    item = findChild(child)
                    if item is not None:
                        return item
            return None
        
        def pruneDeleted(obj):
            if "children" in obj:
                for i, child in enumerate(obj["children"]):
                    if "type" in child and "children" in child and child["type"] == "deleted" and len(child["children"]) == 0:
                        del obj["children"][i]
                        return True
                    else:
                        rtn = pruneDeleted(child)
                        if rtn is not None:
                            return rtn
            return None

        workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
        if workspacemodel is not None:
            updated = workspacemodel.nstr
            jobobj  = findChild(updated)
            if jobobj is not None:
                jobobj["type"] = "deleted"
                if "children" in jobobj:
                    for i, child in enumerate(jobobj["children"]):
                        if child["type"] == "new":
                            del jobobj["children"][i]
            pruned = None
            while pruned is not True:
                pruned = pruneDeleted(updated) 
            workspacemodel.nstr = updated
            workspacemodel.save()
            self.nstr = updated

    def delete(self, project):
        workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
        if project.ensureTrashfolder() and workspacemodel is not None:
            workspace_path = os.path.join(project.dirc, self.dirc)
            workspace_link = os.path.join(project.dirc, self.link)
            trash_path   = os.path.join(project.trashfolder, self.dirc)
            trash_link   = os.path.join(project.trashfolder, self.link)
            if not os.path.exists(workspace_path):
                return 
            if not os.path.isdir(workspace_path):
                return 
            if not os.path.exists(workspace_link):
                return 
            if not os.path.islink(workspace_link):
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
                os.rename(workspace_path, trash_path)
            except OSError as error:
                print("Directory '%s' can not be renamed")
                return
            try:
                os.remove(workspace_link)
                os.symlink(trash_path, trash_link)
            except OSError as error:
                print("Link '%s' can not be renamed")
                return
            workspacemodel.delete()
        workspacemodels = WorkspaceModel.objects.filter(proj=project.id)
        datasetmodels   = DatasetModel.objects.filter(proj=project.id)
        if workspacemodels.count() + datasetmodels.count() > 0:
            project.delete()

    def rename(self, request, project):
        if "new_workspace_name" in request.POST:
            new_workspace_name = request.POST["new_workspace_name"]
            workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
            workspacemodel.name = new_workspace_name
            # strip any non-alphanumeric characters(except _)
            new_workspace_name = new_workspace_name.replace(" ", "_")
            new_workspace_name = "ws_" + re.sub(r'\W+', '', new_workspace_name)
            current_workspace_link = os.path.join(project.dirc, self.link)
            new_workspace_link     = os.path.join(project.dirc, new_workspace_name)
            try :
                os.rename(current_workspace_link, new_workspace_link)
                print("Source path renamed to destination path successfully.")
            except IsADirectoryError:
                print("Source is a file but destination is a directory.")
                return False
            except NotADirectoryError:
                print("Source is a directory but destination is a file.")
                return False
            except PermissionError:
                print("Operation not permitted")
                return False
            except OSError as error:
                print(error)
                return False
            self.link = new_workspace_name
            workspacemodel.link = new_workspace_name
            workspacemodel.save()
            return True
        return False

    def updateDescription(self, request):
        if "new_workspace_description" in request.POST:
            new_workspace_description = request.POST["new_workspace_description"]
            workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
            workspacemodel.desc = new_workspace_description
            workspacemodel.save()
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
    
    
