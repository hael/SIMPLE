# global imports
import os
import datetime
from django.utils import timezone

# local imports
from ..models   import JobClassicModel, WorkspaceModel, ProjectModel
from .simple    import SIMPLE, SIMPLEProjFile

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
            if "csrfmiddlewaretoken" not in key and value != "":
                self.args[key] = value
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        #jobmodels = JobClassicModel.objects.filter(wspc=workspacemodel)
        #self.disp = jobmodels.count() + 1
        self.disp = workspacemodel.jcnt + 1
        workspacemodel.jcnt = self.disp
        workspacemodel.save()
        self.pckg = package
        self.prog = jobtype
        self.name = jobtype.replace("_", " ")
        jobmodel = JobClassicModel(wspc=workspacemodel, cdat=timezone.now(), disp=self.disp, args=self.args, pckg=self.pckg, prog=self.prog, name=self.name, prnt=self.prnt)
        jobmodel.save()
        self.id = jobmodel.id
        self.dirc = str(self.disp) + "_" + self.prog
        if not self.createDir(os.path.join(project.dirc, workspace.dirc)):
            return False
        jobmodel.dirc = self.dirc
        jobmodel.status = "queued"
        jobmodel.save()
        workspace.addChild(self.prnt, self.id)
        parent_dir = os.path.join(project.dirc, workspace.dirc)
        if parentid > 0:
            parentjob = JobClassicModel.objects.filter(id=parentid).first()
            if parentjob.status != "finished":
                return False
            parent_dir = os.path.join(project.dirc, workspace.dirc, parentjob.dirc)
        simple = SIMPLE(pckg=self.pckg)
        if not simple.start(self.args, os.path.join(project.dirc, workspace.dirc, self.dirc), parent_dir, self.prog, self.id):
            return False
        return True
    
    def newSelection(self, request, project, workspace, parentid):
        self.args = {} # ensure empty
        self.prnt = parentid
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        #jobmodels = JobClassicModel.objects.filter(wspc=workspacemodel)
        #self.disp = jobmodels.count() + 1
        self.disp = workspacemodel.jcnt + 1
        workspacemodel.jcnt = self.disp
        workspacemodel.save()
        self.pckg = 'simple'
        self.prog = 'selection'
        self.name = 'user selection'
        jobmodel = JobClassicModel(wspc=workspacemodel, cdat=timezone.now(), disp=self.disp, pckg=self.pckg, prog=self.prog, name=self.name, prnt=self.prnt)
        jobmodel.save()
        self.id = jobmodel.id
        self.dirc = str(self.disp) + "_" + self.prog
        if not self.createDir(os.path.join(project.dirc, workspace.dirc)):
            return False
        if "deselected_cls2D" in request.POST:
            self.args["deselfile"] = "deselected.txt"
            self.args["oritype"]   = "cls2D"
            with open(os.path.join(project.dirc, workspace.dirc, self.dirc, self.args["deselfile"]), "w") as f:
                for deselected in request.POST["deselected_cls2D"].split(','):
                    f.write(deselected + '\n')
        elif "deselected_mic" in request.POST:
            self.args["deselfile"] = "deselected.txt"
            self.args["oritype"]   = "mic"
            with open(os.path.join(project.dirc, workspace.dirc, self.dirc, self.args["deselfile"]), "w") as f:
                for deselected in request.POST["deselected_mic"].split(','):
                    f.write(deselected + '\n')
        jobmodel.args = self.args
        jobmodel.dirc = self.dirc
        jobmodel.status = "queued"
        jobmodel.save()
        workspace.addChild(self.prnt, self.id)
        parentjob = JobClassicModel.objects.filter(id=parentid).first()
        if parentjob.status != "finished":
            return False
        simple = SIMPLE(pckg=self.pckg)
        if not simple.start(self.args, os.path.join(project.dirc, workspace.dirc, self.dirc), os.path.join(project.dirc, workspace.dirc, parentjob.dirc), self.prog, self.id):
            return False
        return True
    
    def linkParticleSet(self, project, workspace, set_proj):
        self.args = {} # ensure empty
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        #jobmodels = JobClassicModel.objects.filter(wspc=workspacemodel)
        #self.disp = jobmodels.count() + 1
        self.disp = workspacemodel.jcnt + 1
        workspacemodel.jcnt = self.disp
        workspacemodel.save()
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
    
    def linkParticleSetFinal(self, project, workspace, set_proj, set_desel):
        self.args = {} # ensure empty
        self.prnt = 0
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        #jobmodels = JobClassicModel.objects.filter(wspc=workspacemodel)
        #self.disp = jobmodels.count() + 1
        self.disp = workspacemodel.jcnt + 1
        workspacemodel.jcnt = self.disp
        workspacemodel.save()
        self.pckg = 'simple'
        self.prog = 'selection'
        self.name = "particle set"
        jobmodel = JobClassicModel(wspc=workspacemodel, cdat=timezone.now(), disp=self.disp, pckg=self.pckg, prog=self.prog, name=self.name, prnt=self.prnt)
        jobmodel.save()
        self.id = jobmodel.id
        self.dirc = str(self.disp) + "_" + self.prog
        if not self.createDir(os.path.join(project.dirc, workspace.dirc)):
            return False
        self.args["deselfile"] = set_desel
        self.args["oritype"]   = "cls2D"
        jobmodel.args = self.args
        jobmodel.dirc = self.dirc
        jobmodel.status = "queued"
        jobmodel.save()
        workspace.addChild(self.prnt, self.id)
        simple = SIMPLE(pckg=self.pckg)
        if not simple.start(self.args, os.path.join(project.dirc, workspace.dirc, self.dirc), os.path.join(project.dirc, workspace.dirc), self.prog, self.id, parent_proj=set_proj):
            return False
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
            self.args = jobmodel.args
            self.prnt = jobmodel.prnt
            self.status = jobmodel.status
            
    def createDir(self, parent_dir):
        if self.dirc == "":
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

    def markComplete(self, project, workspace):
        self.status = "finished"
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        jobmodel.status = self.status
        jobmodel.save()
        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        workspacemodel.mdat = timezone.now()
        workspacemodel.save()
        for childjob in JobClassicModel.objects.filter(prnt=self.id):
            if childjob.status == "queued":
                simple = SIMPLE(pckg=childjob.pckg)
                simple.start(childjob.args, os.path.join(project.dirc, workspace.dirc, childjob.dirc), os.path.join(project.dirc, workspace.dirc, self.dirc), childjob.prog, childjob.id)

    def getProjectStats(self):
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        projfile = os.path.join(jobmodel.wspc.proj.dirc, jobmodel.wspc.dirc, self.dirc, "workspace.simple")
        simpleprojfile = SIMPLEProjFile(projfile)
        return simpleprojfile.getGlobalStats()

    def getProjectFieldStats(self, oritype, fromp=None, top=None, sortkey=None, sortasc=None, hist=False, boxes=False, plotkey=None):
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        projfile = os.path.join(jobmodel.wspc.proj.dirc, jobmodel.wspc.dirc, self.dirc, "workspace.simple")
        simpleprojfile = SIMPLEProjFile(projfile)
        return simpleprojfile.getFieldStats(oritype, fromp, top, sortkey, sortasc, hist, boxes, plotkey)

    def updateStats(self, stats_json, project, workspace):
        jobmodel  = JobClassicModel.objects.filter(id=self.id).first()
        response  = {}
        if "job_heartbeat" in stats_json:
            jobmodel.heartbeat = timezone.now()
        if "job" in stats_json and "terminate" in stats_json["job"]:
            jobmodel.status = "finished"
            for childjob in JobClassicModel.objects.filter(prnt=self.id):
                if childjob.status == "queued":
                    simple = SIMPLE(pckg=childjob.pckg)
                    simple.start(childjob.args, os.path.join(project.dirc, workspace.dirc, childjob.dirc), os.path.join(project.dirc, workspace.dirc, self.dirc), childjob.prog, childjob.id)
        else:
            jobmodel.status = "running"
            response = jobmodel.update
            jobmodel.update = {}
        jobmodel.save()
        return response

    def delete(self, project, workspace):
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if workspace.ensureTrashfolder(project) and jobmodel is not None:
            workspace.removeChild(self.id)
            job_path   = os.path.join(project.dirc, workspace.dirc, self.dirc)
            trash_path = os.path.join(workspace.trashfolder, self.dirc)
            if not os.path.exists(job_path):
                return 
            if not os.path.isdir(job_path):
                return 
            if os.path.exists(trash_path):
                return 
            if os.path.isdir(trash_path):
                return
            try:
                os.rename(job_path, trash_path)
            except OSError as error:
                print("Directory '%s' can not be renamed")
                return
            jobmodel.delete()
    
    def update_description(self, description):
        self.desc = description
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        jobmodel.desc = self.desc
        jobmodel.save()