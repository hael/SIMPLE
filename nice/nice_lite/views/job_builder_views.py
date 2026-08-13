"""Job-builder view for stream mode.

This module renders ``jobbuilder.html`` and prepares:
- stream-specific user inputs (optionally prefilled from a selected stream job)
- SIMPLE and SINGLE program catalogs derived from batch UI JSON metadata
"""

# global imports
import copy

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
from ..helpers import (
    clear_checksum_cookies,
    get_job_id,
    get_project_id,
    get_workspace_id,
    print_error,
)


_BATCH_LAUNCHER_KEYS = {"prg", "projfile", "mkdir", "niceprocid", "niceserver"}


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
    """Allow-list submitted values against authoritative UI JSON metadata."""
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
            continue
        options = user_input.get("options")
        if isinstance(options, list) and value not in options:
            return None, f"invalid value for input: {key}"
        args[key] = value
    return args, None


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

    context = {}
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

    batchjob = BatchJob()
    if not batchjob.new(workspace_obj, package, program, args):
        print_error("failed to create new batch job")
        messages.add_message(request, messages.ERROR, "failed to create batch job")
    return redirect("nice_lite:workspace")
