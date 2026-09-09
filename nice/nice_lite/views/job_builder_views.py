"""Job-builder view for Stream and batch modes.

This module renders ``jobbuilder.html`` and prepares:
- stream-specific user inputs (optionally prefilled from a selected stream job)
- SIMPLE and SINGLE program catalogs, optionally prefilled from a batch job
- completed batch projects and stream snapshots that can seed a batch job
"""

# global imports
import copy
import logging
import os

# django imports
from django import forms
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect, render
from django.views.decorators.http import require_POST

# local imports
from ..data_structures.batchjob import BatchJob
from ..data_structures.simple import SIMPLEBatch, SIMPLEStream
from ..data_structures.streamjob import StreamJob
from ..data_structures.workspace import Workspace
from ..models import JobModel
from ..helpers import (
    clear_checksum_cookies,
    get_integer,
    get_job_id,
    get_project_id,
    get_workspace_id,
)


_BATCH_LAUNCHER_KEYS = {"prg", "projfile", "mkdir", "niceprocid", "niceserver"}
_BATCH_WORKSPACE_SOURCE = "workspace"
_BATCH_JOB_SOURCE_PREFIX = "job"
_BATCH_SNAPSHOT_SOURCE_PREFIX = "snapshot"
_BATCH_CLASS_SELECTION_FILENAME = "class_selection.txt"


logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Internal Helpers
# ------------------------------------------------------------------


class _BatchProgramArgumentsForm(forms.Form):
    """Build typed batch fields from the authoritative SIMPLE UI metadata."""

    def __init__(self, data, program_cfg):
        super().__init__(data=data)
        self.input_metadata = []
        for section_name, section_inputs in program_cfg.items():
            if section_name == "program" or not isinstance(section_inputs, list):
                continue
            for user_input in section_inputs:
                if not isinstance(user_input, dict):
                    continue
                key = user_input.get("key")
                if (
                    not isinstance(key, str)
                    or key == ""
                    or key in _BATCH_LAUNCHER_KEYS
                ):
                    continue
                self.fields[key] = self._field_for_input(key, user_input)
                self.input_metadata.append(user_input)

    @staticmethod
    def _field_for_input(key, user_input):
        required = bool(user_input.get("required"))
        required_error = f"missing required input: {key}"
        keytype = user_input.get("keytype")
        if keytype == "int":
            invalid_error = f"invalid integer value for input: {key}"
            return forms.IntegerField(
                required=required,
                min_value=0,
                error_messages={
                    "required": required_error,
                    "invalid": invalid_error,
                    "min_value": invalid_error,
                },
            )
        if keytype in ("float", "num"):
            invalid_error = f"invalid numeric value for input: {key}"
            field_options = {
                "required": required,
                "error_messages": {
                    "required": required_error,
                    "invalid": invalid_error,
                    "min_value": invalid_error,
                },
            }
            if keytype == "float":
                field_options["min_value"] = 0
            return forms.FloatField(**field_options)
        return forms.CharField(
            required=required,
            strip=True,
            error_messages={"required": required_error},
        )


def _is_job_accessible(jobmodel, username=None):
    """Return True when job model resolves and belongs to the authenticated user."""
    if jobmodel is None:
        return False
    dset = getattr(jobmodel, "dset", None)
    if dset is None:
        return False
    owner = (getattr(dset, "user", "") or "").strip()
    if username is None:
        return True
    if owner == "":
        return False
    return owner == username


def _collect_programs(
    batchui,
    executable_name,
    prefill_program=None,
    prefill_args=None,
    rerun_of=None,
    class_selection_of=None,
):
    """Collect program metadata and optionally prefill one saved command."""
    programs = []
    program_inputs = []
    prefill_args = prefill_args if isinstance(prefill_args, dict) else {}

    if not isinstance(batchui, dict):
        return programs, program_inputs

    for prg, prg_cfg in batchui.items():
        if not isinstance(prg_cfg, dict):
            continue

        program_meta = prg_cfg.get("program")
        if not isinstance(program_meta, dict):
            continue

        executable = program_meta.get("executable")
        if executable not in (executable_name, "all"):
            continue

        sections = []
        for section_name, section_inputs in prg_cfg.items():
            if section_name == "program":
                continue
            visible_inputs = []
            if isinstance(section_inputs, list):
                for entry in section_inputs:
                    if (
                        not isinstance(entry, dict)
                        or entry.get("key") in _BATCH_LAUNCHER_KEYS
                    ):
                        continue
                    user_input = copy.deepcopy(entry)
                    key = user_input.get("key")
                    if prg == prefill_program and key in prefill_args:
                        value = prefill_args[key]
                        options = user_input.get("options")
                        if isinstance(options, list):
                            value = next(
                                (
                                    option for option in options
                                    if str(option) == str(value)
                                ),
                                value,
                            )
                        user_input["value"] = value
                    visible_inputs.append(user_input)
            if visible_inputs:
                sections.append({
                    "name": section_name,
                    "inputs": visible_inputs,
                })

        display_name = program_meta.get("display_name") or prg.replace("_", " ")
        programs.append({
            "prg": prg,
            "disp": display_name,
            "desc": program_meta.get("summary", ""),
        })
        program_input = {
            "prg": prg,
            "disp": display_name,
            "requirements": program_meta.get("requirements", []),
            "sections": sections,
        }
        if (
            prg == prefill_program
            and isinstance(rerun_of, int)
            and not isinstance(rerun_of, bool)
            and rerun_of > 0
        ):
            program_input["rerun_of"] = rerun_of
        if (
            prg == prefill_program
            and isinstance(class_selection_of, int)
            and not isinstance(class_selection_of, bool)
            and class_selection_of > 0
        ):
            program_input["class_selection_of"] = class_selection_of
        program_inputs.append(program_input)

    return programs, program_inputs


def _get_batch_program(batchui, package, program):
    """Return a registered program definition for the selected executable."""
    executable = {"simple": "simple_exec", "single": "single_exec"}.get(package)
    if executable is None or not isinstance(batchui, dict):
        return None
    program_cfg = batchui.get(program)
    if not isinstance(program_cfg, dict):
        return None
    program_meta = program_cfg.get("program")
    if not isinstance(program_meta, dict):
        return None
    if program_meta.get("executable") not in (executable, "all"):
        return None
    return program_cfg


def _collect_batch_args(post, program_cfg):
    """Allow-list and normalize values against authoritative UI JSON metadata."""
    form = _BatchProgramArgumentsForm(post, program_cfg)
    if not form.is_valid():
        for user_input in form.input_metadata:
            key = user_input.get("key")
            if key in form.errors:
                return None, str(form.errors[key][0])
        return None, "invalid batch input"

    args = {}
    for user_input in form.input_metadata:
        key = user_input.get("key")
        raw_value = post.get(key, "")
        value = raw_value.strip() if isinstance(raw_value, str) else str(raw_value).strip()
        cleaned_value = form.cleaned_data.get(key)
        if value == "":
            if key == "nthr" and user_input.get("has_default") and user_input.get("default") is not None:
                # The scheduler needs the effective thread count even when the
                # user accepts the UI default.
                try:
                    cleaned_value = form.fields[key].clean(user_input["default"])
                except forms.ValidationError as error:
                    return None, str(error.messages[0])
            else:
                # Leave other untouched optional inputs out of the command.
                # SIMPLE will apply its own defaults, and mutually exclusive
                # inputs such as import_movies filetab/dir_movies must not both
                # be materialized.
                continue

        keytype = user_input.get("keytype")
        if keytype == "int":
            value = str(cleaned_value)

        options = user_input.get("options")
        if isinstance(options, list) and options and value not in [str(option) for option in options]:
            return None, f"invalid value for input: {key}"
        args[key] = value
    return args, None


def _validate_batch_requirements(program_cfg, args):
    """Validate cross-field requirements published in the batch UI metadata."""
    program_meta = program_cfg.get("program") if isinstance(program_cfg, dict) else None
    requirements = program_meta.get("requirements") if isinstance(program_meta, dict) else None
    if not isinstance(requirements, list):
        return None

    # The batch launcher always supplies the selected project source.
    supplied_keys = {"projfile"}
    supplied_keys.update(
        key for key, value in args.items()
        if isinstance(key, str) and str(value).strip() != ""
    )

    for requirement in requirements:
        if not isinstance(requirement, dict):
            continue
        keys = [key for key in requirement.get("keys", []) if isinstance(key, str)]
        if not keys:
            continue

        min_selected = requirement.get("min_selected", 0)
        max_selected = requirement.get("max_selected", len(keys))
        if not isinstance(min_selected, int) or not isinstance(max_selected, int):
            continue

        selected = sum(key in supplied_keys for key in keys)
        if selected < min_selected or selected > max_selected:
            help_text = requirement.get("help")
            if isinstance(help_text, str) and help_text.strip():
                return help_text.strip()
            label = requirement.get("label") or "Input requirement"
            return f"{label}: select between {min_selected} and {max_selected} inputs"
    return None


def _validate_batch_program_args(program, args):
    """Validate program-specific contracts not expressible in UI field metadata."""
    if program != "import_movies":
        return None

    has_filetab = bool(str(args.get("filetab", "")).strip())
    has_movies_dir = bool(str(args.get("dir_movies", "")).strip())
    if has_filetab and has_movies_dir:
        return "choose either a movie-file list or an input movies directory, not both"
    if not has_filetab and not has_movies_dir:
        return "choose a movie-file list or an input movies directory"
    return None


def _is_batch_job(jobmodel):
    """Return True when a shared JobModel record represents a batch job."""
    metadata = getattr(jobmodel, "master_stats", None)
    return isinstance(metadata, dict) and metadata.get("job_type") == "batch"


def _resolve_class_selection_prefill(jobmodel):
    """Resolve the fixed infile and output project for one ab initio 2D job."""
    batch_job = BatchJob(id=jobmodel.id)
    loaded_jobmodel = batch_job.get_jobmodel()
    if loaded_jobmodel is None or loaded_jobmodel.id != jobmodel.id:
        return None, None, "selected ab initio 2D job is unavailable"

    job_dir = batch_job.get_safe_job_dir()
    result_project = batch_job.get_result_project_path()
    if job_dir is None or result_project is None:
        return None, None, "selected ab initio 2D output is unavailable"

    infile_path = os.path.abspath(
        os.path.join(job_dir, _BATCH_CLASS_SELECTION_FILENAME)
    )
    resolved_infile = os.path.realpath(infile_path)
    try:
        infile_is_safe = os.path.commonpath((job_dir, resolved_infile)) == job_dir
    except ValueError:
        infile_is_safe = False
    if (
        not infile_is_safe
        or os.path.islink(infile_path)
        or not os.path.isfile(infile_path)
    ):
        return None, None, "class selection infile is unavailable"
    return result_project, resolved_infile, None


def _resolve_class_selection_source(workspace_obj, source_id, username):
    """Revalidate one ab initio 2D source before launching selection."""
    if (
        not isinstance(source_id, int)
        or isinstance(source_id, bool)
        or source_id <= 0
    ):
        return None, None, "invalid 2D class selection source"

    jobmodel = JobModel.objects.filter(
        id=source_id,
        dset_id=workspace_obj.get_id(),
    ).first()
    metadata = (
        jobmodel.master_stats
        if jobmodel is not None and isinstance(jobmodel.master_stats, dict)
        else {}
    )
    if (
        not _is_job_accessible(jobmodel, username)
        or not _is_batch_job(jobmodel)
        or jobmodel.status != "finished"
        or metadata.get("package") != "simple"
        or metadata.get("program") != "abinitio2D"
    ):
        return None, None, "invalid 2D class selection source"

    project_path, infile_path, error = _resolve_class_selection_prefill(jobmodel)
    if error is not None:
        return None, None, error
    return project_path, infile_path, None


def _batch_job_source(jobmodel, workspace_dir):
    """Build a safe source descriptor for one completed batch job project."""
    if not _is_batch_job(jobmodel) or jobmodel.status != "finished":
        return None
    if not isinstance(workspace_dir, str) or not isinstance(jobmodel.dirc, str):
        return None
    if jobmodel.dirc in ("", ".", "..") or jobmodel.dirc != os.path.basename(jobmodel.dirc):
        return None

    workspace_root = os.path.realpath(workspace_dir)
    project_path = os.path.realpath(os.path.join(workspace_root, jobmodel.dirc, "workspace.simple"))
    try:
        if os.path.commonpath((workspace_root, project_path)) != workspace_root:
            return None
    except ValueError:
        return None
    if not os.path.isfile(project_path):
        return None

    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    job_name = (jobmodel.name or "").strip() or metadata.get("program") or jobmodel.dirc
    return {
        "key": f"{_BATCH_JOB_SOURCE_PREFIX}:{jobmodel.id}",
        "label": f"job {jobmodel.disp} - {job_name}",
        "path": project_path,
        "metadata": {
            "type": "batch_job",
            "batch_job_id": jobmodel.id,
        },
    }


def _snapshot_source(jobmodel, particle_set, workspace_dir):
    """Build a safe, completed stream-snapshot source descriptor."""
    if _is_batch_job(jobmodel) or not isinstance(particle_set, dict):
        return None
    if not isinstance(workspace_dir, str) or not isinstance(jobmodel.dirc, str):
        return None
    if particle_set.get("type") != "snapshot":
        return None
    if "time" not in particle_set and "ctime" not in particle_set:
        return None

    set_id = particle_set.get("id")
    filename = particle_set.get("filename")
    if not isinstance(set_id, int) or set_id <= 0:
        return None
    if not isinstance(filename, str) or filename != os.path.basename(filename):
        return None
    snapshot_dir = os.path.splitext(filename)[0]
    if not filename.endswith(".simple") or snapshot_dir in ("", ".", ".."):
        return None

    workspace_root = os.path.realpath(workspace_dir)
    project_path = os.path.realpath(os.path.join(
        workspace_root,
        jobmodel.dirc,
        "classification_2D",
        "snapshots",
        snapshot_dir,
        filename,
    ))
    try:
        if os.path.commonpath((workspace_root, project_path)) != workspace_root:
            return None
    except ValueError:
        return None
    if not os.path.isfile(project_path):
        return None

    snapshot_name = particle_set.get("name") or filename
    return {
        "key": f"{_BATCH_SNAPSHOT_SOURCE_PREFIX}:{jobmodel.id}:{set_id}",
        "label": f"stream {jobmodel.disp} - {snapshot_name}",
        "path": project_path,
        "metadata": {
            "type": "stream_snapshot",
            "stream_job_id": jobmodel.id,
            "particle_set_id": set_id,
            "filename": filename,
        },
    }


def _collect_batch_job_sources(workspace_obj):
    """List completed batch projects, newest first, for the selected workspace."""
    workspace_dir = workspace_obj.get_absdir()
    if not isinstance(workspace_dir, str):
        return []

    sources = []
    jobs = JobModel.objects.filter(dset=workspace_obj.get_id()).order_by("-id")
    for jobmodel in jobs:
        source = _batch_job_source(jobmodel, workspace_dir)
        if source is not None:
            sources.append({"key": source["key"], "label": source["label"]})
    return sources


def _collect_batch_snapshot_sources(workspace_obj):
    """List completed snapshot projects available in the selected workspace."""
    workspace_dir = workspace_obj.get_absdir()
    if not isinstance(workspace_dir, str):
        return []

    sources = []
    jobs = JobModel.objects.filter(dset=workspace_obj.get_id()).order_by("-id")
    for jobmodel in jobs:
        stats = jobmodel.particle_sets_stats
        particle_sets = stats.get("particle_sets") if isinstance(stats, dict) else None
        if not isinstance(particle_sets, list):
            continue
        for particle_set in particle_sets:
            source = _snapshot_source(jobmodel, particle_set, workspace_dir)
            if source is not None:
                sources.append({"key": source["key"], "label": source["label"]})
    return sources


def _resolve_batch_project_source(workspace_obj, source_key):
    """Resolve an allow-listed batch project source without trusting a path from POST."""
    if source_key in (None, "", _BATCH_WORKSPACE_SOURCE):
        return None, None, None
    if not isinstance(source_key, str):
        return None, None, "invalid batch project source"

    parts = source_key.split(":")
    if len(parts) == 2 and parts[0] == _BATCH_JOB_SOURCE_PREFIX:
        if not parts[1].isdecimal():
            return None, None, "invalid batch project source"
        jobmodel = JobModel.objects.filter(
            id=int(parts[1]),
            dset=workspace_obj.get_id(),
        ).first()
        source = _batch_job_source(jobmodel, workspace_obj.get_absdir())
        if source is None:
            return None, None, "batch job project is unavailable"
        return source["path"], source["metadata"], None

    if len(parts) != 3 or parts[0] != _BATCH_SNAPSHOT_SOURCE_PREFIX:
        return None, None, "invalid batch project source"
    if not parts[1].isdecimal() or not parts[2].isdecimal():
        return None, None, "invalid batch project source"

    job_id = int(parts[1])
    set_id = int(parts[2])
    jobmodel = JobModel.objects.filter(id=job_id, dset=workspace_obj.get_id()).first()
    if jobmodel is None or _is_batch_job(jobmodel):
        return None, None, "invalid batch project source"

    stats = jobmodel.particle_sets_stats
    particle_sets = stats.get("particle_sets") if isinstance(stats, dict) else None
    if not isinstance(particle_sets, list):
        return None, None, "invalid batch project source"
    particle_set = next(
        (
            entry for entry in particle_sets
            if isinstance(entry, dict) and entry.get("id") == set_id
        ),
        None,
    )
    source = _snapshot_source(jobmodel, particle_set, workspace_obj.get_absdir())
    if source is None:
        return None, None, "batch snapshot is unavailable"
    return source["path"], source["metadata"], None


def _default_batch_source_key(workspace_obj):
    """Return the newest eligible batch project key or the workspace seed."""
    sources = _collect_batch_job_sources(workspace_obj)
    if sources:
        return sources[0]["key"]
    return _BATCH_WORKSPACE_SOURCE


def _default_batch_project_file(workspace_obj):
    """Return the inherited batch project path, falling back to workspace.simple."""
    workspace_dir = workspace_obj.get_absdir()
    if not isinstance(workspace_dir, str):
        return ""

    jobs = JobModel.objects.filter(dset=workspace_obj.get_id()).order_by("-id")
    for jobmodel in jobs:
        source = _batch_job_source(jobmodel, workspace_dir)
        if source is not None:
            return source["path"]

    workspace_project = os.path.realpath(os.path.join(workspace_dir, "workspace.simple"))
    if os.path.isfile(workspace_project):
        return workspace_project
    return ""


def _batch_project_metadata_for_path(workspace_obj, project_path):
    """Describe a selected project file, retaining known job/snapshot provenance."""
    workspace_dir = workspace_obj.get_absdir()
    jobs = JobModel.objects.filter(dset=workspace_obj.get_id()).order_by("-id")
    for jobmodel in jobs:
        source = _batch_job_source(jobmodel, workspace_dir)
        if source is not None and source["path"] == project_path:
            return source["metadata"]

        stats = getattr(jobmodel, "particle_sets_stats", None)
        particle_sets = stats.get("particle_sets") if isinstance(stats, dict) else None
        if not isinstance(particle_sets, list):
            continue
        for particle_set in particle_sets:
            source = _snapshot_source(jobmodel, particle_set, workspace_dir)
            if source is not None and source["path"] == project_path:
                return source["metadata"]

    return {
        "type": "project_file",
        "filename": os.path.relpath(project_path, os.path.realpath(workspace_dir)),
    }


def _resolve_batch_project_file(workspace_obj, project_file):
    """Validate and resolve a posted .simple project inside the selected workspace."""
    if project_file in (None, ""):
        return None, None, None
    if not isinstance(project_file, str):
        return None, None, "invalid batch project file"

    project_file = project_file.strip()
    if project_file == "":
        return None, None, None
    if os.path.splitext(project_file)[1].lower() != ".simple":
        return None, None, "input project must be a .simple file"

    workspace_dir = workspace_obj.get_absdir()
    if not isinstance(workspace_dir, str):
        return None, None, "batch project file is unavailable"
    workspace_root = os.path.realpath(workspace_dir)
    if not os.path.isabs(project_file):
        project_file = os.path.join(workspace_root, project_file)
    project_path = os.path.realpath(project_file)
    try:
        if os.path.commonpath((workspace_root, project_path)) != workspace_root:
            return None, None, "batch project file is outside the selected workspace"
    except ValueError:
        return None, None, "invalid batch project file"
    if not os.path.isfile(project_path):
        return None, None, "batch project file is unavailable"

    return (
        project_path,
        _batch_project_metadata_for_path(workspace_obj, project_path),
        None,
    )


def resolve_recorded_batch_project(workspace_obj, metadata):
    """Resolve persisted batch input provenance for a safe rerun form."""
    if not isinstance(metadata, dict):
        return None, None, "batch job metadata is unavailable"
    if metadata.get("program") == "new_project":
        return None, None, None

    source = metadata.get("source")
    if source is None:
        parent_id = metadata.get("parent")
        if (
            isinstance(parent_id, int)
            and not isinstance(parent_id, bool)
            and parent_id > 0
        ):
            return _resolve_batch_project_source(
                workspace_obj,
                f"{_BATCH_JOB_SOURCE_PREFIX}:{parent_id}",
            )
        workspace_dir = workspace_obj.get_absdir()
        if not isinstance(workspace_dir, str):
            return None, None, "batch project file is unavailable"
        project_path = os.path.join(workspace_dir, "workspace.simple")
        resolved_path, _resolved_source, error = _resolve_batch_project_file(
            workspace_obj,
            project_path,
        )
        return resolved_path, None, error
    if not isinstance(source, dict):
        return None, None, "batch project source is invalid"

    source_type = source.get("type")
    if source_type == "batch_job":
        source_id = source.get("batch_job_id")
        if not isinstance(source_id, int) or isinstance(source_id, bool) or source_id <= 0:
            return None, None, "batch project source is invalid"
        return _resolve_batch_project_source(
            workspace_obj,
            f"{_BATCH_JOB_SOURCE_PREFIX}:{source_id}",
        )
    if source_type == "stream_snapshot":
        stream_id = source.get("stream_job_id")
        set_id = source.get("particle_set_id")
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value <= 0
            for value in (stream_id, set_id)
        ):
            return None, None, "batch project source is invalid"
        return _resolve_batch_project_source(
            workspace_obj,
            f"{_BATCH_SNAPSHOT_SOURCE_PREFIX}:{stream_id}:{set_id}",
        )
    if source_type == "project_file":
        filename = source.get("filename")
        if not isinstance(filename, str) or not filename.strip():
            return None, None, "batch project source is invalid"
        return _resolve_batch_project_file(workspace_obj, filename)
    if source_type == "workspace":
        workspace_dir = workspace_obj.get_absdir()
        if not isinstance(workspace_dir, str):
            return None, None, "batch project file is unavailable"
        return _resolve_batch_project_file(
            workspace_obj,
            os.path.join(workspace_dir, "workspace.simple"),
        )
    return None, None, "batch project source is invalid"


def _is_workspace_accessible(workspace_obj, project_id, username):
    """Return True when the selected workspace belongs to project and user."""
    if workspace_obj is None or not project_id:
        return False
    workspacemodel = workspace_obj.get_workspacemodel()
    if workspacemodel is None or not workspace_obj.in_project(project_id):
        return False
    return (workspacemodel.user or "").strip() == username


# ------------------------------------------------------------------
# Views
# ------------------------------------------------------------------

@login_required(login_url="/login")
def view_job_builder(request):
    """Render the job builder, optionally prefilled from an existing job."""
    template = "jobbuilder.html"
    jobid = get_job_id(request)
    streamui = None
    batchui = None
    stream_args = None
    selected_batch_jobmodel = None
    selected_class_selection_jobmodel = None
    class_selection_prefill_requested = request.GET.get("class_selection") == "1"
    clear_selected_job_cookie = False

    if jobid is not None:
        streamjob = StreamJob(jobid)
        streamjobmodel = streamjob.get_jobmodel()
        if not _is_job_accessible(streamjobmodel, request.user.username):
            messages.add_message(request, messages.ERROR, "selected job is not accessible")
            # Drop stale invalid selection state to avoid repeated access errors.
            clear_selected_job_cookie = True
        elif _is_batch_job(streamjobmodel):
            metadata = streamjobmodel.master_stats
            if class_selection_prefill_requested:
                if (
                    streamjobmodel.status != "finished"
                    or metadata.get("package") != "simple"
                    or metadata.get("program") != "abinitio2D"
                ):
                    messages.add_message(
                        request,
                        messages.ERROR,
                        "selected batch job cannot create a 2D class selection",
                    )
                    clear_selected_job_cookie = True
                else:
                    selected_class_selection_jobmodel = streamjobmodel
            else:
                if (
                    streamjobmodel.status not in BatchJob.RERUNNABLE_STATUSES
                    or metadata.get("package") not in ("simple", "single")
                    or not isinstance(metadata.get("program"), str)
                    or not metadata["program"]
                    or not isinstance(streamjobmodel.args, dict)
                ):
                    messages.add_message(
                        request,
                        messages.ERROR,
                        "selected batch job cannot be rerun",
                    )
                    clear_selected_job_cookie = True
                else:
                    selected_batch_jobmodel = streamjobmodel
        elif isinstance(streamjobmodel.args, dict):
            stream_args = streamjobmodel.args

    simplestream = SIMPLEStream()
    if simplestream.loadUIJSON():
        # Copy UI metadata before injecting request-specific values.
        streamui = copy.deepcopy(simplestream.get_ui())
    else:
        messages.add_message(request, messages.ERROR, "failed to read stream ui JSON")
    simplebatch = SIMPLEBatch()
    if simplebatch.loadUIJSON():
        batchui = simplebatch.get_ui()
    else:
        messages.add_message(request, messages.ERROR, "failed to read batch ui JSON")

    context = {
        "default_batch_project_file": "",
    }
    if isinstance(streamui, dict):
        user_inputs = streamui.get("user_inputs")
        if isinstance(user_inputs, list):
            if isinstance(stream_args, dict):
                for user_input in user_inputs:
                    if not isinstance(user_input, dict):
                        continue
                    key = user_input.get("key")
                    if key in stream_args:
                        user_input["value"] = stream_args[key]
            context["stream_user_inputs"] = user_inputs

    workspace_id = get_workspace_id(request)
    project_id = get_project_id(request)
    workspace_obj = None
    if workspace_id is not None and project_id is not None:
        selected_workspace = Workspace(workspace_id)
        if _is_workspace_accessible(
            selected_workspace,
            project_id,
            request.user.username,
        ):
            workspace_obj = selected_workspace
            context["default_batch_project_file"] = _default_batch_project_file(workspace_obj)

    batch_prefill = None
    if selected_class_selection_jobmodel is not None:
        if (
            workspace_obj is None
            or selected_class_selection_jobmodel.dset_id != workspace_id
        ):
            messages.add_message(
                request,
                messages.ERROR,
                "selected batch job is outside the current workspace",
            )
            clear_selected_job_cookie = True
        elif _get_batch_program(batchui, "simple", "selection") is None:
            messages.add_message(
                request,
                messages.ERROR,
                "2D class selection is no longer available",
            )
            clear_selected_job_cookie = True
        else:
            parent_proj, infile_path, error = _resolve_class_selection_prefill(
                selected_class_selection_jobmodel
            )
            if error is not None:
                logger.error("job_builder_class_selection: %s", error)
                messages.add_message(request, messages.ERROR, error)
                clear_selected_job_cookie = True
            else:
                context["default_batch_project_file"] = parent_proj
                batch_prefill = {
                    "class_selection_of": selected_class_selection_jobmodel.id,
                    "package": "simple",
                    "program": "selection",
                    "args": {
                        "infile": infile_path,
                        "oritype": "cls2D",
                    },
                }
                context["batch_prefill"] = batch_prefill
    elif selected_batch_jobmodel is not None:
        metadata = selected_batch_jobmodel.master_stats
        package = metadata["package"]
        program = metadata["program"]
        if workspace_obj is None or selected_batch_jobmodel.dset_id != workspace_id:
            messages.add_message(
                request,
                messages.ERROR,
                "selected batch job is outside the current workspace",
            )
            clear_selected_job_cookie = True
        elif _get_batch_program(batchui, package, program) is None:
            messages.add_message(
                request,
                messages.ERROR,
                "selected batch program is no longer available",
            )
            clear_selected_job_cookie = True
        else:
            parent_proj, _source, error = resolve_recorded_batch_project(
                workspace_obj,
                metadata,
            )
            if error is not None:
                logger.error("job_builder_rerun: %s", error)
                messages.add_message(request, messages.ERROR, error)
                parent_proj = ""
            context["default_batch_project_file"] = parent_proj or ""
            batch_prefill = {
                "job_id": selected_batch_jobmodel.id,
                "package": package,
                "program": program,
                "args": selected_batch_jobmodel.args,
            }
            context["batch_prefill"] = batch_prefill

    if isinstance(batchui, dict):
        prefill_package = batch_prefill.get("package") if batch_prefill else None
        prefill_program = batch_prefill.get("program") if batch_prefill else None
        prefill_args = batch_prefill.get("args") if batch_prefill else None
        rerun_of = batch_prefill.get("job_id") if batch_prefill else None
        class_selection_of = (
            batch_prefill.get("class_selection_of") if batch_prefill else None
        )
        simple_programs, simple_program_inputs = _collect_programs(
            batchui,
            "simple_exec",
            prefill_program=prefill_program if prefill_package == "simple" else None,
            prefill_args=prefill_args,
            rerun_of=rerun_of,
            class_selection_of=class_selection_of,
        )
        single_programs, single_program_inputs = _collect_programs(
            batchui,
            "single_exec",
            prefill_program=prefill_program if prefill_package == "single" else None,
            prefill_args=prefill_args,
            rerun_of=rerun_of,
            class_selection_of=class_selection_of,
        )
        context["simple_programs"] = simple_programs
        context["simple_program_inputs"] = simple_program_inputs
        context["single_programs"] = single_programs
        context["single_program_inputs"] = single_program_inputs

    response = render(request, template, context)
    if clear_selected_job_cookie:
        response.delete_cookie(key="selected_job_id")

    # Ensure this page starts with fresh checksums after prior iframe/navigation updates.
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
@require_POST
def view_create_batch(request):
    """Validate the selected UI program and launch its batch commander job."""
    workspace_id = get_workspace_id(request)
    project_id = get_project_id(request)
    workspace_obj = Workspace(workspace_id)
    if not _is_workspace_accessible(workspace_obj, project_id, request.user.username):
        logger.error("create_batch: invalid workspace access for workspace %s", workspace_id)
        messages.add_message(request, messages.ERROR, "invalid workspace selection")
        return redirect("nice_lite:workspace")

    package = request.POST.get("package", "")
    program = request.POST.get("program", "")
    simplebatch = SIMPLEBatch(pckg=package)
    program_cfg = _get_batch_program(simplebatch.get_ui(), package, program)
    if program_cfg is None:
        logger.error("create_batch: unknown %s program %s", package, program)
        messages.add_message(request, messages.ERROR, "invalid batch program selection")
        return redirect("nice_lite:workspace")

    rerun_jobmodel = None
    if "rerun_of" in request.POST:
        rerun_of = get_integer(request.POST, "rerun_of", silent=True)
        if rerun_of is not None:
            rerun_jobmodel = JobModel.objects.filter(
                id=rerun_of,
                dset_id=workspace_id,
            ).first()
        rerun_metadata = (
            rerun_jobmodel.master_stats
            if rerun_jobmodel is not None
            and isinstance(rerun_jobmodel.master_stats, dict)
            else {}
        )
        if (
            rerun_jobmodel is None
            or not _is_job_accessible(rerun_jobmodel, request.user.username)
            or not _is_batch_job(rerun_jobmodel)
            or rerun_jobmodel.status not in BatchJob.RERUNNABLE_STATUSES
            or rerun_metadata.get("package") != package
            or rerun_metadata.get("program") != program
        ):
            logger.error("create_batch: invalid rerun source %s", rerun_of)
            messages.add_message(request, messages.ERROR, "invalid batch rerun selection")
            return redirect("nice_lite:workspace")

    class_selection_project = None
    class_selection_infile = None
    class_selection_of = None
    if "class_selection_of" in request.POST:
        class_selection_of = get_integer(
            request.POST,
            "class_selection_of",
            silent=True,
        )
        if (
            package != "simple"
            or program != "selection"
            or rerun_jobmodel is not None
        ):
            class_selection_of = None
        if class_selection_of is None:
            logger.error("create_batch: invalid 2D class selection source")
            messages.add_message(
                request,
                messages.ERROR,
                "invalid 2D class selection source",
            )
            return redirect("nice_lite:workspace")
        (
            class_selection_project,
            class_selection_infile,
            error,
        ) = _resolve_class_selection_source(
            workspace_obj,
            class_selection_of,
            request.user.username,
        )
        if error is not None:
            logger.error("create_batch: %s", error)
            messages.add_message(request, messages.ERROR, error)
            return redirect("nice_lite:workspace")

    submitted_args = request.POST
    if class_selection_of is not None:
        # These values describe the source project's class segment. Resolve
        # them again on POST instead of trusting editable or stale form state.
        submitted_args = request.POST.copy()
        submitted_args["infile"] = class_selection_infile
        submitted_args["oritype"] = "cls2D"

    args, error = _collect_batch_args(submitted_args, program_cfg)
    if error is not None:
        logger.error("create_batch: %s", error)
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    error = _validate_batch_requirements(program_cfg, args)
    if error is not None:
        logger.error("create_batch: %s", error)
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    error = _validate_batch_program_args(program, args)
    if error is not None:
        logger.error("create_batch: %s", error)
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    project_file = class_selection_project
    if project_file is None:
        project_file = request.POST.get("batch_project_file")
    if project_file in (None, ""):
        project_file = _default_batch_project_file(workspace_obj)
    parent_proj, source_metadata, error = _resolve_batch_project_file(
        workspace_obj,
        project_file,
    )
    if error is not None:
        logger.error("create_batch: %s", error)
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    batchjob = BatchJob()
    launch_options = {}
    program_meta = program_cfg.get("program", {})
    if rerun_jobmodel is not None:
        if isinstance(rerun_jobmodel.name, str) and rerun_jobmodel.name.strip():
            launch_options["display_name"] = rerun_jobmodel.name.strip()
        launch_options["description"] = rerun_jobmodel.desc
        launch_options["rerun_of"] = rerun_jobmodel.id
    else:
        display_name = program_meta.get("display_name")
        if isinstance(display_name, str) and display_name.strip():
            launch_options["display_name"] = display_name.strip()
    if parent_proj is not None:
        launch_options.update({"parent_proj": parent_proj, "source": source_metadata})
    if not batchjob.new(workspace_obj, package, program, args, **launch_options):
        logger.error("failed to create new batch job")
        messages.add_message(request, messages.ERROR, "failed to create batch job")
    else:
        messages.add_message(request, messages.SUCCESS, "batch job queued")
    return redirect("nice_lite:workspace")
