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
from django.shortcuts import render

# local imports
from ..data_structures.simple import SIMPLEBatch, SIMPLEStream
from ..data_structures.streamjob import StreamJob
from ..helpers import clear_checksum_cookies, get_job_id


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
            if isinstance(section_inputs, list) and len(section_inputs) > 0:
                sections.append({
                    "name": section_name,
                    "inputs": section_inputs,
                })

        display_name = prg.replace("_", " ")
        programs.append({
            "prg": prg,
            "disp": display_name,
            "desc": program_meta.get("descr_short", ""),
        })
        program_inputs.append({
            "prg": prg,
            "disp": display_name,
            "sections": sections,
        })

    return programs, program_inputs


# ------------------------------------------------------------------
# Views
# ------------------------------------------------------------------

@login_required(login_url="/login/")
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