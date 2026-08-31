# global imports
import os
import shutil
import signal
import time

# django imports
from django.db import transaction
from django.utils import timezone

# local imports
from ..helpers import directory_exists, ensure_directory, print_error
from ..models import JobModel, WorkspaceModel
from .simple import SIMPLEBatch, SIMPLEProjFile, SIMPLEProject
from .job import Job
from .workspace import Workspace


class BatchJob(Job):
    """Classic (non-stream) SIMPLE job attached to a workspace."""

    TERMINAL_STATUSES = frozenset(("finished", "failed", "stopped"))
    DELETABLE_STATUSES = TERMINAL_STATUSES | frozenset(("queued",))

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
        self.source = None

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
        self.jobmodel = JobModel.objects.filter(id=self.id).first()
        metadata = self.jobmodel.master_stats if self.jobmodel is not None else None
        if not isinstance(metadata, dict) or metadata.get("job_type") != "batch":
            self.jobmodel = None
            self.id = 0
            self.absdir = None
            return

        self.name = self.jobmodel.name
        self.desc = self.jobmodel.desc
        self.dirc = self.jobmodel.dirc
        self.cdat = self.jobmodel.cdat
        self.args = self.jobmodel.args
        self.wspc = self.jobmodel.dset
        self.disp = self.jobmodel.disp
        self.prog = metadata.get("program", "")
        self.pckg = metadata.get("package", "")
        self.prnt = metadata.get("parent", 0)
        self.status = self.jobmodel.status
        self.source = metadata.get("source")
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

    @staticmethod
    def _metadata(pckg, prog, parent=0, source=None):
        metadata = {
            "job_type": "batch",
            "package": pckg,
            "program": prog,
            "parent": parent,
        }
        if isinstance(source, dict):
            metadata["source"] = dict(source)
        return metadata

    def new(self, workspace, pckg, prog, args, parent_proj=None, source=None, display_name=None):
        """Create and launch a SIMPLE or SINGLE batch job."""
        if workspace is None or pckg not in ("simple", "single") or not prog:
            print_error("new: invalid batch job configuration")
            return False
        if not isinstance(args, dict):
            print_error("new: args is not a dictionary")
            return False

        workspacemodel = workspace.get_workspacemodel()
        workspace_dir = workspace.get_absdir()
        if workspacemodel is None or not directory_exists(workspace_dir):
            print_error("new: workspace is unavailable")
            return False

        explicit_parent_proj = parent_proj is not None
        if explicit_parent_proj:
            try:
                workspace_root = os.path.realpath(workspace_dir)
                parent_proj = os.path.realpath(parent_proj)
                source_in_workspace = os.path.commonpath((workspace_root, parent_proj)) == workspace_root
            except (TypeError, ValueError):
                source_in_workspace = False
            if not source_in_workspace or not os.path.isfile(parent_proj):
                print_error("new: batch project source is unavailable")
                return False
        else:
            parent_proj = os.path.join(workspace_dir, "workspace.simple")
            if not os.path.isfile(parent_proj):
                if not SIMPLEProject(workspace_dir).create():
                    print_error("new: failed to initialize workspace.simple")
                    return False

        self.pckg = pckg
        self.prog = prog
        if isinstance(display_name, str) and display_name.strip():
            self.name = display_name.strip()
        else:
            self.name = prog.replace("_", " ")
        self.args = dict(args)

        # Reserve the display counter under a row lock so two near-simultaneous
        # Start clicks cannot choose the same job directory.
        with transaction.atomic():
            locked_workspace = WorkspaceModel.objects.select_for_update().get(pk=workspacemodel.pk)
            self.disp = locked_workspace.jcnt + 1
            self.dirc = f"{self.disp}_{prog}"
            jobmodel = JobModel(
                dset=locked_workspace,
                cdat=timezone.now(),
                disp=self.disp,
                args=self.args,
                name=self.name,
                dirc=self.dirc,
                status="queued",
                master_status="queued",
                master_stats=self._metadata(pckg, prog, source=source),
                master_heartbeat=0,
            )
            jobmodel.save()
            locked_workspace.jcnt = self.disp
            locked_workspace.save()

        self.id = jobmodel.id
        self.jobmodel = jobmodel
        self.wspc = locked_workspace
        self.status = "queued"
        self.absdir = self.get_absdir()

        if not self._create_dir(workspace_dir):
            self.status = "failed"
            jobmodel.status = "failed"
            jobmodel.master_status = "failed"
            jobmodel.save()
            return False

        simple = SIMPLEBatch(pckg=pckg)
        launch_options = {"parent_proj": parent_proj} if explicit_parent_proj else {}
        if simple.start(self.args, self.absdir, workspace_dir, prog, self.id, **launch_options):
            return True

        self.status = "failed"
        jobmodel.status = "failed"
        jobmodel.master_status = "failed"
        jobmodel.save()
        return False

    def linkParticleSet(self, project, workspace, set_proj):
        self.args = {}

        workspacemodel = WorkspaceModel.objects.filter(id=workspace.id).first()
        if workspacemodel is None:
            print_error("linkParticleSet: workspace not found")
            return False

        self.disp = workspacemodel.jcnt + 1
        self.name = "particle set"
        self.dirc = str(self.disp) + "_link_particle_set"
        workspace_dir = os.path.join(project.dirc, workspacemodel.dirc)
        if not self._create_dir(workspace_dir):
            return False

        link_dst = os.path.join(workspace_dir, self.dirc, "workspace.simple")
        if not self.createLink(set_proj, link_dst):
            return False

        jobmodel = JobModel(
            dset=workspacemodel,
            cdat=timezone.now(),
            disp=self.disp,
            args={},
            name=self.name,
            dirc=self.dirc,
            status="finished",
            master_status="finished",
            master_stats=self._metadata("simple_stream", "link_particle_set"),
        )
        jobmodel.save()
        workspacemodel.jcnt = self.disp
        workspacemodel.save()

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
        self.pckg = "simple"
        self.prog = "selection"
        self.name = "particle set"
        self.dirc = str(self.disp) + "_" + self.prog
        workspace_dir = os.path.join(project.dirc, workspacemodel.dirc)
        if not self._create_dir(workspace_dir):
            return False

        self.args["deselfile"] = set_desel
        self.args["oritype"] = "cls2D"

        jobmodel = JobModel(
            dset=workspacemodel,
            cdat=timezone.now(),
            disp=self.disp,
            name=self.name,
            args=self.args,
            dirc=self.dirc,
            status="queued",
            master_status="queued",
            master_stats=self._metadata(self.pckg, self.prog, self.prnt),
        )
        jobmodel.save()
        workspacemodel.jcnt = self.disp
        workspacemodel.save()

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

    def _local_process_is_running(self, pid):
        """Return True while the generated local job process still exists."""
        try:
            process_dir = os.path.realpath(os.readlink(f"/proc/{pid}/cwd"))
        except OSError:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return False
            except PermissionError:
                return True
            return True
        return process_dir == os.path.realpath(self.get_absdir() or "")

    def queued_job_can_delete(self):
        """Return True for a stale local queued job whose process has exited."""
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None or jobmodel.status != "queued":
            return False

        job_dir = self.get_absdir()
        if job_dir is None:
            return False

        try:
            with open(os.path.join(job_dir, "job.script"), encoding="utf-8") as script_file:
                if "# localtemplate" not in script_file.read(4096):
                    return False
            with open(os.path.join(job_dir, "nice.pid"), encoding="utf-8") as pid_file:
                pid = int(pid_file.read(32).strip())
        except (OSError, ValueError):
            return False

        return pid > 1 and not self._local_process_is_running(pid)

    def reconcile_local_completion(self):
        """Mark a locally dispatched job finished after a verified normal exit.

        This is a conservative fallback for local jobs whose final NICE callback
        did not reach Django. Scheduler jobs and abnormal local exits are left
        unchanged for their authoritative status path to handle.
        """
        if self.jobmodel is None or self.status not in ("queued", "running"):
            return False

        job_dir = self.get_absdir()
        if job_dir is None:
            return False

        try:
            with open(os.path.join(job_dir, "job.script"), encoding="utf-8") as script_file:
                if "# localtemplate" not in script_file.read(4096):
                    return False
            with open(os.path.join(job_dir, "nice.pid"), encoding="utf-8") as pid_file:
                pid = int(pid_file.read(32).strip())
        except (OSError, ValueError):
            return False

        if pid <= 1 or self._local_process_is_running(pid):
            return False

        try:
            stdout_path = os.path.join(job_dir, "stdout.log")
            stdout_size = os.path.getsize(stdout_path)
            with open(stdout_path, "rb") as stdout_file:
                stdout_file.seek(max(0, stdout_size - 65536))
                stdout_tail = stdout_file.read().decode("utf-8", errors="replace")
        except OSError:
            return False

        normal_stop = any(
            line.startswith("**** ") and line.endswith(" NORMAL STOP ****")
            for line in stdout_tail.splitlines()
        )
        if not normal_stop:
            return False

        with transaction.atomic():
            jobmodel = JobModel.objects.select_for_update().filter(id=self.id).first()
            if jobmodel is None or jobmodel.status not in ("queued", "running"):
                return False
            jobmodel.status = "finished"
            jobmodel.master_status = "finished"
            jobmodel.save(update_fields=("status", "master_status"))

        self.jobmodel = jobmodel
        self.status = "finished"
        return True

    def _get_local_process_group(self, pid):
        """Return an isolated process group owned by this job, or ``None``."""
        job_dir = self.get_absdir()
        if job_dir is None:
            print_error("stop: job directory is unavailable")
            return None

        try:
            process_dir = os.path.realpath(os.readlink(f"/proc/{pid}/cwd"))
            job_dir = os.path.realpath(job_dir)
            process_group = os.getpgid(pid)
        except OSError:
            print_error("stop: job process is not available on this host")
            return None

        if process_dir != job_dir:
            print_error("stop: pid does not belong to the job directory")
            return None
        if process_group != pid:
            print_error("stop: job process group is not isolated")
            return None
        return process_group

    def stop(self):
        """Stop a running local batch job by terminating its process group."""
        with transaction.atomic():
            jobmodel = JobModel.objects.select_for_update().filter(id=self.id).first()
            if jobmodel is None or jobmodel.status != "running":
                print_error("stop: batch job is not running")
                return False

            pid_path = os.path.join(self.get_absdir() or "", "nice.pid")
            try:
                with open(pid_path, encoding="utf-8") as pid_file:
                    pid = int(pid_file.read(32).strip())
            except (OSError, ValueError):
                print_error("stop: nice.pid is missing or invalid")
                return False
            if pid <= 1:
                print_error("stop: nice.pid contains an unsafe process id")
                return False

            process_group = self._get_local_process_group(pid)
            if process_group is None:
                return False
            try:
                os.killpg(process_group, signal.SIGTERM)
            except OSError:
                print_error("stop: failed to terminate the job process group")
                return False

            self.status = "stopped"
            jobmodel.status = "stopped"
            jobmodel.master_status = "stopped"
            jobmodel.master_update = {}
            jobmodel.save(update_fields=("status", "master_status", "master_update"))
            return True

    def updateDescription(self, request):
        if "new_job_description" in request.POST:
            return self.set_description(request.POST["new_job_description"])
        return False

    def markComplete(self, project, workspace):
        del project, workspace
        self.status = "finished"
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return False
        jobmodel.status = self.status
        jobmodel.master_status = self.status
        jobmodel.save()
        return True

    def updateStats(self, stats_json, project, workspace):
        del project, workspace
        with transaction.atomic():
            jobmodel = JobModel.objects.select_for_update().filter(id=self.id).first()
            if jobmodel is None:
                return {}
            # Script-level and in-process callbacks can race at shutdown. Once
            # one of them records a terminal state, do not let a later callback
            # resurrect the job or replace the original outcome.
            if jobmodel.status in self.TERMINAL_STATUSES:
                return {}

            response = {}
            if "job_heartbeat" in stats_json:
                jobmodel.master_heartbeat = int(time.time())

            job_stats = stats_json.get("job")
            if isinstance(job_stats, dict) and "terminate" in job_stats:
                terminal_status = job_stats.get("status")
                if not isinstance(terminal_status, str) or terminal_status not in self.TERMINAL_STATUSES:
                    terminal_status = "finished"
                jobmodel.status = terminal_status
                jobmodel.master_status = terminal_status
            else:
                jobmodel.status = "running"
                jobmodel.master_status = "running"
                response = jobmodel.master_update or {}
                jobmodel.master_update = {}

            jobmodel.save()
            return response

    # ------------------------------------------------------------------
    # Projfile helpers
    # ------------------------------------------------------------------

    def getProjectStats(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return None
        projfile = os.path.join(jobmodel.dset.proj.dirc, jobmodel.dset.dirc, self.dirc, "workspace.simple")
        return SIMPLEProjFile(projfile).getGlobalStats()

    def getProjectFieldStats(self, oritype, fromp=None, top=None, sortkey=None, sortasc=None, hist=False, boxes=False, plotkey=None):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return None
        projfile = os.path.join(jobmodel.dset.proj.dirc, jobmodel.dset.dirc, self.dirc, "workspace.simple")
        return SIMPLEProjFile(projfile).getFieldStats(oritype, fromp, top, sortkey, sortasc, hist, boxes, plotkey)

    # ------------------------------------------------------------------
    # Deletion
    # ------------------------------------------------------------------

    def delete(self, project, workspace):
        """Permanently remove a batch job and reset its workspace job counter."""
        del project
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is None:
            return False
        if jobmodel.status not in self.DELETABLE_STATUSES:
            print_error("delete: batch job status is not deletable")
            return False
        if jobmodel.status == "queued" and not self.queued_job_can_delete():
            print_error("delete: queued batch job is still active or unverifiable")
            return False

        workspace_dir = os.path.realpath(
            os.path.join(jobmodel.dset.proj.dirc, jobmodel.dset.dirc)
        )
        if not directory_exists(workspace_dir):
            return False

        job_dirc = jobmodel.dirc
        if (
            not isinstance(job_dirc, str)
            or not job_dirc.strip()
            or os.path.isabs(job_dirc)
        ):
            return False

        raw_job_path = os.path.abspath(os.path.join(workspace_dir, job_dirc))
        job_path = os.path.realpath(raw_job_path)
        try:
            path_is_safe = (
                job_path != workspace_dir
                and os.path.commonpath((workspace_dir, job_path)) == workspace_dir
            )
        except ValueError:
            path_is_safe = False
        if not path_is_safe or os.path.islink(raw_job_path):
            print_error("delete: unsafe batch job directory")
            return False
        if not directory_exists(job_path):
            return False

        try:
            shutil.rmtree(job_path)
        except OSError:
            print_error("delete: failed to permanently remove batch job directory")
            return False

        workspace_id = jobmodel.dset_id
        jobmodel.delete()
        if (
            getattr(workspace, "id", None) != workspace_id
            or not callable(getattr(workspace, "reconcile_job_counter", None))
        ):
            workspace = Workspace(workspace_id)
        workspace.reconcile_job_counter()
        return True
