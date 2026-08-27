"""Job-builder view for stream mode.

This module renders ``jobbuilder.html`` and prepares:
- stream-specific user inputs (optionally prefilled from a selected stream job)
- SIMPLE and SINGLE program catalogs derived from batch UI JSON metadata
- completed stream snapshots that can explicitly seed a batch job
"""

# global imports
import copy
import math
import os

# django imports
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
    get_job_id,
    get_project_id,
    get_workspace_id,
    print_error,
)


_BATCH_LAUNCHER_KEYS = {"prg", "projfile", "mkdir", "niceprocid", "niceserver"}
_BATCH_WORKSPACE_SOURCE = "workspace"
_BATCH_SNAPSHOT_SOURCE_PREFIX = "snapshot"


# ------------------------------------------------------------------
# Internal Helpers
# ------------------------------------------------------------------


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


def _collect_programs(batchui, executable_name):
    """Collect program metadata and section inputs for a target executable."""
    programs = []
    program_inputs = []

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
            visible_inputs = [
                entry for entry in section_inputs
                if isinstance(entry, dict) and entry.get("key") not in _BATCH_LAUNCHER_KEYS
            ] if isinstance(section_inputs, list) else []
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
        program_inputs.append({
            "prg": prg,
            "disp": display_name,
            "sections": sections,
        })

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
    inputs = []
    for section_name, section_inputs in program_cfg.items():
        if section_name != "program" and isinstance(section_inputs, list):
            inputs.extend(entry for entry in section_inputs if isinstance(entry, dict))

    args = {}
    for user_input in inputs:
        key = user_input.get("key")
        if not isinstance(key, str) or key == "" or key in _BATCH_LAUNCHER_KEYS:
            continue
        value = post.get(key, "").strip()
        if value == "":
            if user_input.get("required"):
                return None, f"missing required input: {key}"
            if key == "nthr" and user_input.get("has_default") and user_input.get("default") is not None:
                # The scheduler needs the effective thread count even when the
                # user accepts the UI default.
                value = str(user_input["default"])
            else:
                # Leave other untouched optional inputs out of the command.
                # SIMPLE will apply its own defaults, and mutually exclusive
                # inputs such as import_movies filetab/dir_movies must not both
                # be materialized.
                continue

        keytype = user_input.get("keytype")
        if keytype == "int":
            try:
                numeric_value = float(value)
                if not math.isfinite(numeric_value) or not numeric_value.is_integer() or numeric_value < 0:
                    raise ValueError
                value = str(int(numeric_value))
            except (TypeError, ValueError, OverflowError):
                return None, f"invalid integer value for input: {key}"
        elif keytype in ("float", "num"):
            try:
                numeric_value = float(value)
                if not math.isfinite(numeric_value) or (keytype == "float" and numeric_value < 0):
                    raise ValueError
            except (TypeError, ValueError, OverflowError):
                return None, f"invalid numeric value for input: {key}"

        options = user_input.get("options")
        if isinstance(options, list) and value not in [str(option) for option in options]:
            return None, f"invalid value for input: {key}"
        args[key] = value
    return args, None


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
    """Render stream job-builder page for a new job or from an existing stream job."""
    template = "jobbuilder.html"
    jobid = get_job_id(request)
    streamui = None
    batchui = None
    args = None
    clear_selected_job_cookie = False

    if jobid is not None:
        streamjob = StreamJob(jobid)
        streamjobmodel = streamjob.get_jobmodel()
        if not _is_job_accessible(streamjobmodel, request.user.username):
            messages.add_message(request, messages.ERROR, "selected job is not accessible")
            # Drop stale invalid selection state to avoid repeated access errors.
            clear_selected_job_cookie = True
        elif isinstance(streamjobmodel.args, dict):
            args = streamjobmodel.args

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

    context = {"batch_snapshot_sources": []}
    if isinstance(streamui, dict):
        user_inputs = streamui.get("user_inputs")
        if isinstance(user_inputs, list):
            if isinstance(args, dict):
                for user_input in user_inputs:
                    if not isinstance(user_input, dict):
                        continue
                    key = user_input.get("key")
                    if key in args:
                        user_input["value"] = args[key]
            context["stream_user_inputs"] = user_inputs
    if isinstance(batchui, dict):
        simple_programs, simple_program_inputs = _collect_programs(batchui, "simple_exec")
        single_programs, single_program_inputs = _collect_programs(batchui, "single_exec")
        context["simple_programs"] = simple_programs
        context["simple_program_inputs"] = simple_program_inputs
        context["single_programs"] = single_programs
        context["single_program_inputs"] = single_program_inputs

    workspace_id = get_workspace_id(request)
    project_id = get_project_id(request)
    if workspace_id is not None and project_id is not None:
        workspace_obj = Workspace(workspace_id)
        if _is_workspace_accessible(workspace_obj, project_id, request.user.username):
            context["batch_snapshot_sources"] = _collect_batch_snapshot_sources(workspace_obj)

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
        print_error(f"create_batch: invalid workspace access for workspace {workspace_id}")
        messages.add_message(request, messages.ERROR, "invalid workspace selection")
        return redirect("nice_lite:workspace")

    package = request.POST.get("package", "")
    program = request.POST.get("program", "")
    simplebatch = SIMPLEBatch(pckg=package)
    program_cfg = _get_batch_program(simplebatch.get_ui(), package, program)
    if program_cfg is None:
        print_error(f"create_batch: unknown {package} program {program}")
        messages.add_message(request, messages.ERROR, "invalid batch program selection")
        return redirect("nice_lite:workspace")

    args, error = _collect_batch_args(request.POST, program_cfg)
    if error is not None:
        print_error(f"create_batch: {error}")
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    error = _validate_batch_program_args(program, args)
    if error is not None:
        print_error(f"create_batch: {error}")
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    parent_proj, source_metadata, error = _resolve_batch_project_source(
        workspace_obj,
        request.POST.get("batch_source", _BATCH_WORKSPACE_SOURCE),
    )
    if error is not None:
        print_error(f"create_batch: {error}")
        messages.add_message(request, messages.ERROR, error)
        return redirect("nice_lite:workspace")

    batchjob = BatchJob()
    launch_options = {}
    if parent_proj is not None:
        launch_options = {"parent_proj": parent_proj, "source": source_metadata}
    if not batchjob.new(workspace_obj, package, program, args, **launch_options):
        print_error("failed to create new batch job")
        messages.add_message(request, messages.ERROR, "failed to create batch job")
    else:
        messages.add_message(request, messages.SUCCESS, "batch job queued")
    return redirect("nice_lite:workspace")
