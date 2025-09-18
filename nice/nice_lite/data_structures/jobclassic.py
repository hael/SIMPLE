# global imports
import os
import datetime
from django.utils import timezone

# local imports
from ..models import JobClassicModel, WorkspaceModel, ProjectModel
from .simple  import SIMPLE

class JobClassic:

    id       = 0
    disp     = 0
    prnt     = 0
    name     = ""
    desc     = ""
    dirc     = ""
    cdat     = ""
    prog     = ""
    args     = {}
    wspc     = 0
    pckg     = None
    status   = "unknown"

    def __init__(self, pckg=None, id=None, request=None):
        if pckg is not None:
            self.pckg = pckg
        if id is not None:
            self.id = id
            self.load()
        elif request is not None:
            self.setIDFromRequest(request)
            self.load()

    def setIDFromRequest(self, request):
        if "jobid" in request.POST:
            test_id_str = request.POST["jobid"]
        elif "jobid" in request.GET:
            test_id_str = request.GET["jobid"]
        if test_id_str.isnumeric():
            self.id = int(test_id_str)
    
    def getAbsDir(self):
        return os.path.join(self.wspc.proj.dirc, self.wspc.dirc, self.dirc)

    def new(self, request, project, workspace, parentid, package, jobtype):
        self.args = {} # ensure empty
        self.prnt = parentid
        for key, value in request.POST.items():
            if "csrfmiddlewaretoken" not in key and value is not "":
                self.args[key] = value
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        jobmodels = JobClassicModel.objects.filter(wspc=workspacemodel)
        self.disp = jobmodels.count() + 1
        self.pckg = package
        self.prog = jobtype
        self.name = jobtype.replace("_", " ")
        jobmodel = JobClassicModel(wspc=workspacemodel, cdat=timezone.now(), disp=self.disp, args=self.args, pckg=self.pckg, prog=self.prog, name=self.name)
        jobmodel.save()
        self.id = jobmodel.id
        self.dirc = str(self.disp) + "_" + self.prog
        if not self.createDir(os.path.join(project.dirc, workspace.dirc)):
            return False
        jobmodel.dirc = self.dirc
        jobmodel.status = "queued"
        jobmodel.save()
        workspace.addChild(self.prnt, self.id)
        parentjob = JobClassicModel.objects.filter(id=parentid).first()
        simple = SIMPLE(pckg=self.pckg)
        if not simple.start(self.args, os.path.join(project.dirc, workspace.dirc, self.dirc), os.path.join(project.dirc, workspace.dirc, parentjob.dirc), self.prog, self.id):
            return False
        return True
    
    def linkParticleSet(self, project, workspace, set_proj):
        self.args = {} # ensure empty
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        jobmodels = JobClassicModel.objects.filter(wspc=workspacemodel)
        self.disp = jobmodels.count() + 1
        self.name = "particle set"
        jobmodel = JobClassicModel(wspc=workspacemodel, cdat=timezone.now(), disp=self.disp, args={}, pckg="simple_stream", name=self.name)
        jobmodel.save()
        self.id = jobmodel.id
        self.dirc = str(self.disp) + "_"  + "link_particle_set" 
        if not self.createDir(os.path.join(project.dirc, workspace.dirc)):
            return False
        if not self.createLink(set_proj, os.path.join(project.dirc, workspace.dirc, self.dirc, "workspace.simple")):
            return False
        jobmodel.dirc = self.dirc
        jobmodel.status = "finished"
        jobmodel.save()
        workspace.addChild(self.prnt, self.id)
        return True
    
    def load(self):
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            self.name = jobmodel.name
            self.desc = jobmodel.desc
            self.dirc = jobmodel.dirc
            self.cdat = jobmodel.cdat
            self.args = jobmodel.args
            self.wspc = jobmodel.wspc
            self.disp = jobmodel.disp
            self.prog = jobmodel.prog
            self.pckg = jobmodel.pckg
            self.status = jobmodel.status
            
    def createDir(self, parent_dir):
        if self.dirc is "":
            return False
        
        if not os.path.exists(parent_dir):
            return False
        
        if not os.path.isdir(parent_dir):
            return False
    
        new_dir_path = os.path.join(parent_dir, self.dirc)

        try:
            os.makedirs(new_dir_path, exist_ok=False)
        except OSError as error:
            print("Directory '%s' can not be created", new_dir_path)
            return False
        return True
    
    def createLink(self, source, destination):
        if not os.path.exists(source):
            return False
        
        if not os.path.isfile(source):
            return False
        try:
            os.symlink(source, destination)
        except OSError as error:
            print("Symlink '%s' can not be created", source, destination)
            return False
        return True
    
    def updateDescription(self, request):
        if "new_job_description" in request.POST:
            new_job_description = request.POST["new_job_description"]
            self.desc = new_job_description
            jobmodel = JobClassicModel.objects.filter(id=self.id).first()
            jobmodel.desc = self.desc
            jobmodel.save()

    def markComplete(self):
        self.status = "finished"
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        jobmodel.status = self.status
        jobmodel.save()

    