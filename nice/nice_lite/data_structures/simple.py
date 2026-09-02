"""Launchers and project-file helpers for SIMPLE stream and classic workflows."""

# global imports
import json
import os
import shlex
import shutil
import stat
import subprocess

# local imports
from ..helpers import *
from ..features import batch_status_callbacks_enabled
from ..models import DispatchModel


# ------------------------------------------------------------------
# Shared scheduler submission helper
# ------------------------------------------------------------------

def _submit(dispatchmodel, script_path, cwd):
    """
    Submit a pre-written shell script via the configured scheduler.

    Supports two modes selected by dispatchmodel.scmd:
      - 'bsub'  : pipes the script into bsub -L /bin/sh (LSF)
      - anything else : launches the script directly as a detached process

    LSF submission waits for bsub to report scheduler acceptance. The file
    handle for bsub stdin is managed via a context manager to prevent resource
    leaks.
    """
    if dispatchmodel.scmd == 'bsub':
        with open(script_path, 'r') as f:
            subprocess.run(
                ['bsub', '-L', '/bin/sh'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=cwd,
                stdin=f,
                text=True,
                check=True
            )
    else:
        subprocess.Popen(
            [dispatchmodel.scmd, script_path, '&'],
            cwd=cwd,
            start_new_session=True
        )


def _classic_status_callback_wrapper(jobid, endpoint):
    """Return shell fragments that report one classic job's full lifecycle.

    The wrapper is shared by SIMPLE and SINGLE jobs. Callback failures are
    logged without changing the scientific program's exit status. When the
    worker token is exported to the compute environment it is forwarded in
    the same header accepted by the NICE API.
    """
    running_payload = json.dumps(
        {"jobid": jobid, "job_heartbeat": {}},
        separators=(",", ":"),
    )
    finished_payload = json.dumps(
        {"jobid": jobid, "job": {"status": "finished", "terminate": True}},
        separators=(",", ":"),
    )
    failed_payload = json.dumps(
        {"jobid": jobid, "job": {"status": "failed", "terminate": True}},
        separators=(",", ":"),
    )
    quoted_endpoint = shlex.quote(endpoint)
    curl_options = (
        "--silent --show-error --fail --connect-timeout 3 --max-time 10 "
        "--output /dev/null --header 'Content-Type: application/json'"
    )
    token_request = (
        f"        curl {curl_options}"
        ' --header "X-Worker-Token: ${NICE_LITE_WORKER_TOKEN}"'
        f' --data "$1" {quoted_endpoint} 2>> nice_status.log\n'
    )
    anonymous_request = (
        f"        curl {curl_options}"
        f' --data "$1" {quoted_endpoint} 2>> nice_status.log\n'
    )
    prefix = (
        "(\n"
        "nice_status_callback() {\n"
        "    if [ -n \"${NICE_LITE_WORKER_TOKEN:-}\" ]; then\n"
        f"{token_request}"
        "    else\n"
        f"{anonymous_request}"
        "    fi\n"
        "}\n"
        f"nice_status_callback {shlex.quote(running_payload)} || true\n"
        "nice_run_job() {\n"
    )
    suffix = (
        "}\n"
        "if nice_run_job; then\n"
        "    nice_job_exit=0\n"
        "else\n"
        "    nice_job_exit=$?\n"
        "fi\n"
        "if [ \"$nice_job_exit\" -eq 0 ]; then\n"
        f"    nice_status_callback {shlex.quote(finished_payload)} || true\n"
        "else\n"
        f"    nice_status_callback {shlex.quote(failed_payload)} || true\n"
        "fi\n"
        "exit \"$nice_job_exit\"\n"
        ")\n"
    )
    return prefix, suffix


# ------------------------------------------------------------------
# Stream job launcher
# ------------------------------------------------------------------

class SIMPLEStream:
    """
    Manages the launch and restart of a SIMPLE stream processing job.

    A stream job runs the multi-stage cryo-EM pipeline (preprocessing,
    optics assignment, picking, classification) under a single master process
    (simple_stream prg=master). Communication with the master uses HTTP via
    the niceprocid / niceserver arguments written into the dispatch script.

    Typical usage:
        stream = SIMPLEStream(args)
        stream.loadUIJSON()          # parse the stream UI definition
        stream.start(absdir, jobid)  # write script and submit
    """

    # class-level constants — shared across all instances
    ui_cmd            = ["simple_private_exec", "prg=print_ui_json"]
    executable        = "simple_stream prg=master"
    tplt_simple_motif = "XXXSIMPLEXXX"   # placeholder in dispatch template for the command
    tplt_nthr_motif   = "XXXNCPUXXX"     # placeholder in dispatch template for thread count
    nthr_master       = 4                # number of threads for the master process

    def __init__(self, args=None):
        self.ui          = None    # parsed UI definition dict, populated by loadUIJSON()
        self.absdir      = None    # absolute path to the job directory
        self.args        = args    # user-supplied argument dict
        self.skip_refgen = False   # True when pickrefs are supplied externally
        self.jobid       = 0
        self.dispid      = 0
        if args is not None and "pickrefs" in args:
            self.skip_refgen = True

    # ------------------------------------------------------------------
    # UI loading
    # ------------------------------------------------------------------

    def loadUIJSON(self):
        """
        Populate self.ui by running simple_private_exec and parsing its JSON output.

        The JSON is expected to contain a "master" key whose sections (excluding
        "program") are flattened into a single user_inputs list.
        Returns True on success, False on any failure.
        """
        try:
            result  = subprocess.run(
                self.ui_cmd,
                capture_output=True,
                check=True,
                text=True
            )
            ui_json = json.loads(result.stdout)
            if "master" not in ui_json:
                print_error("'master' key missing from UI JSON")
                self.ui = None
                return False
            # flatten all non-program sections into a single user_inputs list
            ui_stream = {"user_inputs": []}
            for section, inputs in ui_json["master"].items():
                if section != "program":
                    ui_stream["user_inputs"].extend(inputs)
            self.ui = ui_stream
        except subprocess.CalledProcessError as cpe:
            print_error(cpe.stderr)
            print_error("Failed to call ui json")
            self.ui = None
            return False
        except Exception:
            print_error("Failed to load ui json")
            self.ui = None
            return False
        return True

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_ui(self):
        return self.ui

    def get_skip_refgen(self):
        return self.skip_refgen

    # ------------------------------------------------------------------
    # Launch / restart
    # ------------------------------------------------------------------

    def start(self, absdir, jobid, dispid):
        """
        Validate state, record the job directory and id, then dispatch.
        loadUIJSON() must have been called successfully before start().
        Returns True on success, False on any failure.
        """
        if self.ui is None:
            print_error("ui is none")
            return False
        if self.args is None:
            print_error("args is none")
            return False
        if absdir is None:
            print_error("absdir is none")
            return False
        if jobid <= 0:
            print_error("invalid job id")
            return False
        if dispid <= 0:
            print_error("invalid dispatch id")
            return False
        self.absdir = absdir
        self.jobid  = jobid
        self.dispid = dispid
        if not directory_exists(self.absdir):
            return False
        return self.dispatch()

    def restart(self, base_dir, processname):
        """
        Re-submit an existing sub-process script without rebuilding it.

        Looks up processname in self.ui["processes"] and submits its
        pre-written .script file via the active dispatch configuration.
        Returns True on success, False on any failure.
        """
        if self.ui is None or "processes" not in self.ui:
            print_error("ui has no processes")
            return False
        if not base_dir or not os.path.isdir(base_dir):
            return False
        dispatchmodel = DispatchModel.objects.filter(active=True).last()
        if dispatchmodel is None:
            return False
        for process in self.ui["processes"]:
            if process["name"] == processname:
                script_path = os.path.join(base_dir, process["name"] + ".script")
                if shutil.which(dispatchmodel.scmd) is None:
                    return False
                try:
                    _submit(dispatchmodel, script_path, base_dir)
                except Exception as e:
                    print_error(str(e))
                    return False
        return True

    def dispatch(self):
        """
        Build the stream master dispatch script from the active template and submit it.

        Steps:
          1. Load the active DispatchModel from the DB.
          2. Build the simple_stream command string from user_inputs and static args.
          3. Substitute the command and thread-count placeholders into the template.
          4. Write the filled script to stream_master.script in the job directory.
          5. Submit via _submit().

        Returns True on success, False on any failure.
        """
        dispatchmodel = DispatchModel.objects.filter(active=True).last()
        if dispatchmodel is None:
            print_error("dispatchmodel is none")
            return False
        if self.tplt_simple_motif not in dispatchmodel.tplt:
            print_error("simple motif missing from dispatch template")
            return False
 
        # build the command string: executable + user args + static args + API args
        command_string = self.executable
        if "user_inputs" in self.ui:
            nicedispid_found = False
            for user_input in self.ui["user_inputs"]:
                if user_input["key"] in self.args and user_input["keytype"] != "hidden":
                    arg_key = user_input["key"]
                    arg_val = self.args[arg_key]
                    if "nicedispid" == arg_key:
                        if str(arg_val).strip() != "0":
                            nicedispid_found = True
                            command_string += " " + arg_key + "=" + str(arg_val)
                    else:
                        command_string += " " + arg_key + "=" + str(arg_val)
            # if nicedispid is not supplied by the user or is 0, add it to the command string
            if not nicedispid_found:
                command_string += " nicedispid=" + str(self.dispid)
        command_string += " mkdir=no nparts=1 nthr=" + str(self.nthr_master)
        command_string += " niceprocid=" + str(self.jobid) + " niceserver=" + dispatchmodel.url + "/api"
        command_string += " >> stream_master.log 2>> stream_master.error\n"

        # fill template placeholders and normalise line endings
        dispatch_script = dispatchmodel.tplt
        dispatch_script = dispatch_script.replace(self.tplt_nthr_motif,   str(self.nthr_master))
        dispatch_script = dispatch_script.replace(self.tplt_simple_motif, command_string)
        dispatch_script = dispatch_script.replace("\r\n", "\n")

        dispatch_script_path = os.path.join(self.absdir, "stream_master.script")
        try:
            with open(dispatch_script_path, "w") as scriptfile:
                scriptfile.write(dispatch_script)
        except IOError:
            print_error("File '%s' can not be written")
            return False
        
        try:
            os.chmod(dispatch_script_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        except OSError as e:
            print_error(str(e))
            return False

        if shutil.which(dispatchmodel.scmd) is None:
            print_error("dispatch command doesnt exist")
            return False
        try:
            _submit(dispatchmodel, dispatch_script_path, self.absdir)
        except Exception as e:
            print_error(str(e))
            return False
        return True


# ------------------------------------------------------------------
# Classic (non-stream) job launcher
# ------------------------------------------------------------------

class SIMPLEBatch:
    """
    Manages the launch of a single classic SIMPLE job (simple_exec or single_exec).

    Unlike SIMPLEStream, each job corresponds to one program call rather than a
    pipeline. Except for new_project, the dispatch script copies the workspace
    project file into the job directory, runs update_project, then executes the
    program. A new_project job instead creates its project directly in the job
    directory without inheriting workspace.simple.

    Typical usage:
        simple = SIMPLEBatch(pckg="simple")
        simple.start(args, base_dir, parent_dir, jobtype, jobid)
    """

    # class-level constants
    ui_cmd            = ["simple_private_exec", "prg=print_ui_json"]
    tplt_simple_motif = "XXXSIMPLEXXX"
    tplt_nthr_motif   = "XXXNCPUXXX"

    def __init__(self, pckg=None):
        self.ui          = {}   # parsed UI definition dict
        self.base_dir    = ""   # absolute path to the job directory
        self.args        = {}   # user-supplied argument dict
        self.jobtype     = ""   # program name passed as prg=
        self.jobid       = 0
        self.parent_proj = ""   # path to the parent workspace.simple project file
        self.executable  = ""   # "simple_exec" or "single_exec"
        if pckg == "simple":
            self.executable = "simple_exec"
        elif pckg == "single":
            self.executable = "single_exec"
        self.loadUIJSON()

    # ------------------------------------------------------------------
    # UI loading
    # ------------------------------------------------------------------

    def loadUIJSON(self):
        """
        Populate self.ui by running simple_private_exec and parsing its JSON output.

        The JSON is expected to contain a "master" key whose sections (excluding
        "program") are flattened into a single user_inputs list.
        Returns True on success, False on any failure.
        """
        try:
            result  = subprocess.run(
                self.ui_cmd,
                capture_output=True,
                check=True,
                text=True
            )
            ui_json = json.loads(result.stdout)

            self.ui = ui_json
        except subprocess.CalledProcessError as cpe:
            print_error(cpe.stderr)
            print_error("Failed to call ui json")
            self.ui = None
            return False
        except Exception:
            print_error("Failed to load ui json")
            self.ui = None
            return False
        return True
    
    # ------------------------------------------------------------------
    # Launch
    # ------------------------------------------------------------------

    def start(self, args, base_dir, parent_dir, jobtype, jobid, parent_proj=None):
        """
        Set job parameters, validate the filesystem state, and dispatch.

        parent_proj defaults to <parent_dir>/workspace.simple when not supplied.
        Returns True on success, False on any failure.
        """
        self.base_dir    = base_dir
        self.args        = args
        self.jobtype     = jobtype
        self.jobid       = jobid
        self.parent_proj = parent_proj if parent_proj is not None else os.path.join(parent_dir, "workspace.simple")
        if not self.base_dir:
            return False
        if not os.path.isdir(self.base_dir):
            return False
        if not os.path.isdir(parent_dir):
            return False
        if not os.path.isfile(self.parent_proj):
            return False
        if not self.executable or not self.jobtype:
            return False
        return self.dispatch()

    def dispatch(self):
        """
        Build the job dispatch script from the active template and submit it.

        Steps:
          1. Load the active DispatchModel from the DB.
          2. Build the command, inheriting the parent project when required.
          3. Substitute placeholders into the template.
          4. Write the filled script to job.script in the job directory.
          5. Submit via _submit().

        Returns True on success, False on any failure.
        """
        dispatchmodel = DispatchModel.objects.filter(active=True).last()
        if dispatchmodel is None:
            return False
        if self.tplt_simple_motif not in dispatchmodel.tplt:
            return False

        creates_project = self.executable == "simple_exec" and self.jobtype == "new_project"
        status_endpoint = dispatchmodel.url + "/api_classic"
        status_prefix = ""
        status_suffix = ""
        if batch_status_callbacks_enabled():
            status_prefix, status_suffix = _classic_status_callback_wrapper(
                self.jobid,
                status_endpoint,
            )

        command_string = status_prefix
        setup_failure_guard = " || return $?" if status_suffix else ""
        if not creates_project:
            command_string += (
                "cp -v " + shlex.quote(self.parent_proj)
                + " workspace.simple" + setup_failure_guard + "\n"
            )
            command_string += (
                "simple_exec prg=update_project projfile=workspace.simple"
                + setup_failure_guard + "\n"
            )
        command_string += self.executable + " prg=" + shlex.quote(self.jobtype)
        for key, val in self.args.items():
            command_string += " " + key + "=" + shlex.quote(str(val))
        if creates_project and "dir" not in self.args:
            command_string += " dir=."
        command_string += " mkdir=no"
        if not creates_project:
            command_string += " projfile=workspace.simple"
        command_string += " niceprocid=" + str(self.jobid) + " niceserver=" + shlex.quote(status_endpoint)
        command_string += " >> stdout.log 2>> stderr.log\n"
        command_string += status_suffix

        # fill template placeholders and normalise line endings
        try:
            scheduler_nthr = max(1, int(float(self.args.get("nthr", 1))))
        except (TypeError, ValueError, OverflowError):
            scheduler_nthr = 1

        dispatch_script = dispatchmodel.tplt
        dispatch_script = dispatch_script.replace(self.tplt_nthr_motif,   str(scheduler_nthr))
        dispatch_script = dispatch_script.replace(self.tplt_simple_motif, command_string)
        dispatch_script = dispatch_script.replace("\r\n", "\n")

        dispatch_script_path = os.path.join(self.base_dir, "job.script")
        try:
            with open(dispatch_script_path, "w") as scriptfile:
                scriptfile.write(dispatch_script)
        except IOError:
            print_error("File '%s' can not be written" % dispatch_script_path)
            return False
        try:
            os.chmod(dispatch_script_path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        except OSError as e:
            print_error(str(e))
            return False

        if shutil.which(dispatchmodel.scmd) is None:
            return False
        try:
            _submit(dispatchmodel, dispatch_script_path, self.base_dir)
        except Exception as e:
            print_error(str(e))
            return False
        return True
    
    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_ui(self):
        return self.ui


# ------------------------------------------------------------------
# Project file utilities
# ------------------------------------------------------------------

class SIMPLEProject:
    """
    Initialises a new SIMPLE project file (workspace.simple) in a directory.

    Used when a workspace first needs to run a batch job, so downstream jobs
    can copy a valid project file.
    """

    cmd = ["simple_exec", "prg=new_project", "projname=workspace"]

    def __init__(self, absdir):
        self.absdir = absdir  # directory in which to create the project file

    def create(self):
        """Run simple_exec to create workspace.simple in self.absdir."""
        if self.absdir is None or not os.path.isdir(self.absdir):
            return False
        cmd = self.cmd + ["dir=" + self.absdir]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as cpe:
            print_error(str(cpe))
            return False
        except OSError as error:
            print_error(str(error))
            return False
        return os.path.isfile(os.path.join(self.absdir, "workspace.simple"))


class SIMPLEProjFile:
    """
    Read-only interface to a SIMPLE project file (workspace.simple).

    Wraps simple_exec prg=print_project_info and prg=print_project_field to
    expose project-level and per-field statistics as Python dicts.
    """

    cmd_header = ["simple_exec", "prg=print_project_info",  "json=yes"]
    cmd_field  = ["simple_exec", "prg=print_project_field", "json=yes"]
    command_timeout = 30

    def __init__(self, projfile):
        self.projfile = projfile  # absolute path to the .simple project file
        self.ui       = {}        # most recently returned stats dict

    def _run(self, cmd):
        """Run a simple_exec command and return the parsed JSON, or {} on failure."""
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                check=True,
                text=True,
                timeout=self.command_timeout,
            )
            return json.loads(result.stdout)
        except subprocess.CalledProcessError as cpe:
            print_error(cpe.stderr)
            return {}
        except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError) as error:
            print_error(str(error))
            return {}

    def _check_projfile(self):
        """Return True if self.projfile exists and is a regular file or symlink."""
        if self.projfile is None:
            return False
        return os.path.exists(self.projfile) and (
            os.path.isfile(self.projfile) or os.path.islink(self.projfile)
        )

    def getGlobalStats(self):
        """
        Return project-level statistics (micrograph/particle counts, CTF summary, etc.)
        as a dict. Returns {} if the project file is missing or the command fails.
        """
        if not self._check_projfile():
            return {}
        cmd     = self.cmd_header + ["projfile=" + self.projfile]
        self.ui = self._run(cmd)
        return self.ui

    def getFieldStats(self, oritype, fromp=None, top=None, sortkey=None,
                      sortasc=True, hist=False, boxes=False, plotkey=None):
        """
        Return per-field statistics for the given orientation type.

        oritype   — orientation type (e.g. "ptcl2D", "mic")
        fromp/top — pagination range (both must be supplied together)
        sortkey   — field name to sort by
        sortasc   — sort ascending when True (default)
        hist      — include histogram data
        boxes     — include box coordinates
        plotkey   — field to use as plot x-axis (requires sortkey)

        Returns a dict, or {} if the project file is missing or the command fails.
        """
        if oritype is None or not self._check_projfile():
            return {}
        cmd = self.cmd_field + ["oritype=" + oritype, "projfile=" + self.projfile]
        if fromp is not None and top is not None:
            cmd += ["fromp=" + str(fromp), "top=" + str(top)]
        if sortkey is not None:
            cmd.append("sort=" + sortkey)
        if not sortasc:
            cmd.append("sort_asc=no")
        if hist:
            cmd.append("hist=yes")
        if boxes:
            cmd.append("boxes=yes")
        if plotkey is not None and sortkey is not None:
            cmd.append("plot_key=" + plotkey)
        self.ui = self._run(cmd)
        return self.ui
