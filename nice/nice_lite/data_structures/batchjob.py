# global imports
import os

# django imports
from django.utils import timezone

# local imports
from ..helpers import directory_exists, ensure_directory, print_error
from ..models import WorkspaceModel
from .simple import SIMPLEBatch, SIMPLEProjFile
from .job import Job


class BatchJob(Job):
    """Classic (non-stream) SIMPLE job attached to a workspace."""

    def __init__(self, pckg=None, id=None, request=None):
        super().__init__(id=None)
        self.disp = 0
        self.prnt = 0
        self.name = ""
        self.dirc = ""
        self.cdat = ""
        self.prog = ""
        self.args = {}
        self.wspc = None
        self.pckg = pckg
        self.status = "unknown"

        if id is not None:
            self.id = id
            self.load()
        elif request is not None:
            self.setIDFromRequest(request)
            if self.id > 0:
                self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def setIDFromRequest(self, request):
        self.set_id_from_request(request, key="jobid")

    def load(self):
        """Populate fields from DB. Resets to empty state if not found."""
        self.jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if self.jobmodel is None:
            self.id = 0
            self.absdir = None
            return

        self.name = self.jobmodel.name
        self.desc = self.jobmodel.desc
        self.dirc = self.jobmodel.dirc
        self.cdat = self.jobmodel.cdat
        self.args = self.jobmodel.args
        self.wspc = self.jobmodel.wspc
        self.disp = self.jobmodel.disp
        self.prog = self.jobmodel.prog
        self.pckg = self.jobmodel.pckg
        self.prnt = self.jobmodel.prnt
        self.status = self.jobmodel.status
        self.absdir = self.get_absdir()

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_absdir(self):
        if self.wspc is None:
            return None
        return os.path.join(self.wspc.proj.dirc, self.wspc.dirc, self.dirc)

    def getAbsDir(self):
        """Backward-compatible alias."""
        return self.get_absdir()

    def get_jobmodel(self):
        return self.jobmodel

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _create_dir(self, parent_dir):
        if self.dirc == "":
            print_error("_create_dir: empty dir name")
            return False
        if not directory_exists(parent_dir):
            print_error("_create_dir: parent directory missing")
            return False

        new_dir_path = os.path.join(parent_dir, self.dirc)
        if directory_exists(new_dir_path):
            print_error("_create_dir: destination already exists")
            return False

        return ensure_directory(new_dir_path)

    def createDir(self, parent_dir):
        """Backward-compatible alias."""
        return self._create_dir(parent_dir)

    def createLink(self, source, destination):
        if not os.path.exists(source):
            print_error("createLink: source missing")
            return False
        if not os.path.isfile(source):
            print_error("createLink: source is not a file")
            return False
        try:
            os.symlink(source, destination)
        except OSError:
            print_error("createLink: symlink creation failed")
            return False
        return True

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def linkParticleSet(self, project, workspace, set_proj):
        self.args = {}

        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        if workspacemodel is None:
            print_error("linkParticleSet: workspace not found")
            return False

        self.disp = workspacemodel.jcnt + 1
        workspacemodel.jcnt = self.disp
        workspacemodel.save()

        self.name = "particle set"
        self.dirc = str(self.disp) + "_link_particle_set"
        workspace_dir = os.path.join(project.dirc, workspacemodel.dirc)
        if not self._create_dir(workspace_dir):
            return False

        link_dst = os.path.join(workspace_dir, self.dirc, "workspace.simple")
        if not self.createLink(set_proj, link_dst):
            return False

        jobmodel = JobClassicModel(
            wspc=workspacemodel,
            cdat=timezone.now(),
            disp=self.disp,
            args={},
            pckg="simple_stream",
            name=self.name,
            dirc=self.dirc,
            status="finished",
        )
        jobmodel.save()

        self.id = jobmodel.id
        self.jobmodel = jobmodel
        self.wspc = workspacemodel
        self.status = "finished"
        self.absdir = self.get_absdir()
        return True

    def linkParticleSetFinal(self, project, workspace, set_proj, set_desel):
        self.args = {}
        self.prnt = 0

        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        if workspacemodel is None:
            print_error("linkParticleSetFinal: workspace not found")
            return False

        self.disp = workspacemodel.jcnt + 1
        workspacemodel.jcnt = self.disp
        workspacemodel.save()

        self.pckg = "simple"
        self.prog = "selection"
        self.name = "particle set"
        self.dirc = str(self.disp) + "_" + self.prog
        workspace_dir = os.path.join(project.dirc, workspacemodel.dirc)
        if not self._create_dir(workspace_dir):
            return False

        self.args["deselfile"] = set_desel
        self.args["oritype"] = "cls2D"

        jobmodel = JobClassicModel(
            wspc=workspacemodel,
            cdat=timezone.now(),
            disp=self.disp,
            pckg=self.pckg,
            prog=self.prog,
            name=self.name,
            prnt=self.prnt,
            args=self.args,
            dirc=self.dirc,
            status="queued",
        )
        jobmodel.save()

        self.id = jobmodel.id
        self.jobmodel = jobmodel
        self.wspc = workspacemodel
        self.status = "queued"
        self.absdir = self.get_absdir()

        simple = SIMPLEBatch(pckg=self.pckg)
        return simple.start(
            self.args,
            os.path.join(workspace_dir, self.dirc),
            workspace_dir,
            self.prog,
            self.id,
            parent_proj=set_proj,
        )

    # ------------------------------------------------------------------
    # Updates / completion
    # ------------------------------------------------------------------

    def updateDescription(self, request):
        if "new_job_description" in request.POST:
            return self.set_description(request.POST["new_job_description"])
        return False

    def markComplete(self, project, workspace):
        del project, workspace
        self.status = "finished"
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return False
        jobmodel.status = self.status
        jobmodel.save()
        return True

    def updateStats(self, stats_json, project, workspace):
        del project, workspace
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return {}

        response = {}
        if "job_heartbeat" in stats_json:
            jobmodel.heartbeat = timezone.now()

        if "job" in stats_json and "terminate" in stats_json["job"]:
            jobmodel.status = "finished"
        else:
            jobmodel.status = "running"
            response = jobmodel.update
            jobmodel.update = {}

        jobmodel.save()
        return response

    # ------------------------------------------------------------------
    # Projfile helpers
    # ------------------------------------------------------------------

    def getProjectStats(self):
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return None
        projfile = os.path.join(jobmodel.wspc.proj.dirc, jobmodel.wspc.dirc, self.dirc, "workspace.simple")
        return SIMPLEProjFile(projfile).getGlobalStats()

    def getProjectFieldStats(self, oritype, fromp=None, top=None, sortkey=None, sortasc=None, hist=False, boxes=False, plotkey=None):
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return None
        projfile = os.path.join(jobmodel.wspc.proj.dirc, jobmodel.wspc.dirc, self.dirc, "workspace.simple")
        return SIMPLEProjFile(projfile).getFieldStats(oritype, fromp, top, sortkey, sortasc, hist, boxes, plotkey)

    # ------------------------------------------------------------------
    # Deletion
    # ------------------------------------------------------------------

    def delete(self, project, workspace):
        del project
        jobmodel = JobClassicModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return False

        workspace_dir = os.path.join(jobmodel.wspc.proj.dirc, jobmodel.wspc.dirc)
        trash_dir = os.path.join(workspace_dir, "TRASH")
        if not ensure_directory(trash_dir):
            return False

        job_path = os.path.join(workspace_dir, self.dirc)
        trash_path = os.path.join(trash_dir, self.dirc)
        if not directory_exists(job_path):
            return False
        if os.path.exists(trash_path):
            return False

        try:
            os.rename(job_path, trash_path)
        except OSError:
            print_error("delete: failed to move job to trash")
            return False

        jobmodel.delete()
        return True
