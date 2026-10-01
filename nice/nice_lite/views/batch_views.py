"""Batch detail views and lifecycle actions for SIMPLE and SINGLE jobs."""

import json
import logging
import math
import os
import struct

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import FileResponse, HttpResponse, HttpResponseBadRequest, JsonResponse
from django.shortcuts import redirect, render
from django.template.loader import render_to_string
from django.urls import reverse
from django.views.decorators.cache import cache_control
from django.views.decorators.http import require_GET, require_POST

from ..data_structures.batchjob import BatchJob
from ..data_structures.class_selection import (
    ClassSelectionError,
    load_batch_class_selection,
    validate_deselected_class_ids,
)
from ..data_structures.project import Project
from ..data_structures.simple import SIMPLEBatch
from ..data_structures.workspace import Workspace
from ..helpers import clear_checksum_cookies, get_job_id
from ..models import JobModel


logger = logging.getLogger(__name__)


_BATCH_LAUNCHER_KEYS = {"prg", "projfile", "mkdir", "niceprocid", "niceserver"}
_BATCH_MICROGRAPH_PAGE_SIZE = 50
# manual picking has no established particle box size yet; use a fixed placeholder for written box files.
_MANUALPICK_BOX_SIZE_PX = 100


def _get_accessible_batch_job(request, log_context, job_id=None):
    """Return an owned batch job and model selected by the request."""
    resolved_job_id = job_id if job_id is not None else get_job_id(request)
    if (
        not isinstance(resolved_job_id, int)
        or isinstance(resolved_job_id, bool)
        or resolved_job_id <= 0
    ):
        logger.error("%s: missing job id", log_context)
        return None, None

    batch_job = BatchJob(id=resolved_job_id)
    jobmodel = batch_job.get_jobmodel()
    workspace = getattr(jobmodel, "dset", None)
    owner = (getattr(workspace, "user", "") or "").strip()
    if jobmodel is None or owner != request.user.username:
        logger.error("%s: access denied for job %s", log_context, resolved_job_id)
        return None, None
    return batch_job, jobmodel


def _commander_arguments(package, program):
    """Return ordered input metadata from the selected commander UI."""
    if package == "simple":
        executable = "simple_exec"
    elif package == "single":
        executable = "single_exec"
    else:
        return []
    if not isinstance(program, str) or not program:
        return []

    batch_ui = SIMPLEBatch(pckg=package).get_ui()
    if not isinstance(batch_ui, dict):
        return []
    program_cfg = batch_ui.get(program)
    if not isinstance(program_cfg, dict):
        return []
    program_meta = program_cfg.get("program")
    if (
        not isinstance(program_meta, dict)
        or program_meta.get("executable") not in (executable, "all")
    ):
        return []

    arguments = []
    seen_keys = set()
    for section_name, section_inputs in program_cfg.items():
        if section_name == "program" or not isinstance(section_inputs, list):
            continue
        for user_input in section_inputs:
            if not isinstance(user_input, dict):
                continue
            key = user_input.get("key")
            label = user_input.get("label")
            if (
                not isinstance(key, str)
                or not key
                or key in _BATCH_LAUNCHER_KEYS
                or key in seen_keys
            ):
                continue
            seen_keys.add(key)
            arguments.append({
                "key": key,
                "label": label.strip() if isinstance(label, str) and label.strip() else key,
                "has_default": bool(user_input.get("has_default")) and user_input.get("default") is not None,
                "default": user_input.get("default"),
                "visibility": user_input.get("visibility") or "standard",
            })
    return arguments


def _recorded_project_file(jobmodel):
    """Return the workspace-relative project source recorded for one job."""
    if (
        getattr(jobmodel, "pckg", None) == "simple"
        and getattr(jobmodel, "prog", None) in ("new_project", "reproject")
    ):
        # SIMPLE dispatch deliberately omits projfile for these programs.
        return None

    metadata = getattr(jobmodel, "master_stats", None)
    if not isinstance(metadata, dict):
        return None
    source = metadata.get("source")
    if not isinstance(source, dict):
        return None

    source_type = source.get("type")
    if source_type == "project_file":
        filename = source.get("filename")
        if isinstance(filename, str) and filename.strip():
            return filename
        return None
    if source_type == "workspace":
        return "workspace.simple"
    if source_type not in ("batch_job", "stream_snapshot"):
        return None

    source_id_key = (
        "batch_job_id" if source_type == "batch_job" else "stream_job_id"
    )
    source_id = source.get(source_id_key)
    workspace_id = getattr(jobmodel, "dset_id", None)
    if (
        not isinstance(source_id, int)
        or isinstance(source_id, bool)
        or source_id <= 0
        or not isinstance(workspace_id, int)
        or isinstance(workspace_id, bool)
        or workspace_id <= 0
    ):
        return None

    source_job = JobModel.objects.filter(
        id=source_id,
        dset_id=workspace_id,
    ).first()
    source_dir = getattr(source_job, "dirc", None)
    if (
        not isinstance(source_dir, str)
        or source_dir in ("", ".", "..")
        or source_dir != os.path.basename(source_dir)
    ):
        return None

    if source_type == "batch_job":
        if getattr(source_job, "pckg", None) not in ("simple", "single"):
            return None
        return os.path.join(source_dir, "workspace.simple")

    filename = source.get("filename")
    particle_set_id = source.get("particle_set_id")
    if (
        getattr(source_job, "pckg", None) in ("simple", "single")
        or not isinstance(filename, str)
        or filename != os.path.basename(filename)
        or not filename.endswith(".simple")
        or not isinstance(particle_set_id, int)
        or isinstance(particle_set_id, bool)
        or particle_set_id <= 0
    ):
        return None
    snapshot_dir = os.path.splitext(filename)[0]
    return os.path.join(
        source_dir,
        "classification_2D",
        "snapshots",
        snapshot_dir,
        filename,
    )


def _argument_rows(jobmodel):
    """Build submitted/default/unset rows from saved args and commander UI."""
    raw_saved_args = jobmodel.args if isinstance(jobmodel.args, dict) else {}
    saved_args = {str(key): value for key, value in raw_saved_args.items()}
    definitions = _commander_arguments(
        jobmodel.pckg,
        jobmodel.prog,
    )
    arguments = []
    known_keys = set()

    project_file = _recorded_project_file(jobmodel)
    if project_file is not None:
        arguments.append({
            "key": "projfile",
            "label": "Project file",
            "value": project_file,
            "origin": "submitted",
            "submitted": True,
            "visibility": "standard",
        })
        known_keys.add("projfile")

    for definition in definitions:
        key = definition["key"]
        known_keys.add(key)
        if key in saved_args:
            value = saved_args[key]
            origin = "submitted"
        elif definition["has_default"]:
            value = definition["default"]
            origin = "default"
        else:
            value = None
            origin = "unset"
        arguments.append({
            "key": key,
            "label": definition["label"],
            "value": value,
            "origin": origin,
            "submitted": origin == "submitted",
            "visibility": definition["visibility"],
        })

    # Preserve legacy or newly removed commander arguments recorded on the job.
    for key in sorted(saved_args):
        if key not in known_keys:
            arguments.append({
                "key": key,
                "label": key,
                "value": saved_args[key],
                "origin": "submitted",
                "submitted": True,
                "visibility": "standard",
            })
    # Submitted arguments are surfaced first; stable sort preserves definition order otherwise.
    arguments.sort(key=lambda argument: not argument["submitted"])
    return arguments


def _present_volume_kinds(volume_outputs):
    """Return the {key, label} kinds actually present in volume_outputs, in fixed order."""
    present = {volume["kind"] for volume in volume_outputs if "kind" in volume}
    return [
        {"key": key, "label": label}
        for key, label in BatchJob.VOLUME_KINDS
        if key in present
    ]


def _cls3d_stages(jobstats, stage_volume_outputs):
    """Group cls3D per-stage states with that stage's own browser-safe volume outputs.

    jobstats["cls3D"] is a dict of stage key -> list of state entries (reprojtiles/
    FSC/oridist/volpath/lppath/pprocpath/pprocmirrpath); stage_volume_outputs entries
    carry a matching "stage" key (see BatchJob.get_stage_volume_outputs), so each
    stage's Mol* volume <select> only ever offers that stage's own volumes. Entries
    keep their "path" (like streaming) so the viewer serves volumes via the
    project-scoped nice_lite:volume route instead of a job-scoped one.
    """
    cls3d = jobstats.get("cls3D") if isinstance(jobstats, dict) else None
    if not isinstance(cls3d, dict):
        return []
    stages = []
    for stage_key, states in cls3d.items():
        if not isinstance(states, list):
            continue
        outputs = [
            output for output in stage_volume_outputs
            if output.get("stage") == stage_key
        ]
        stages.append({
            "key": stage_key,
            "states": states,
            "volume_outputs": outputs,
            "volume_kinds": _present_volume_kinds(outputs),
        })
    return stages


def _positive_page_number(request, query_name):
    """Return a positive page query value, defaulting safely to one."""
    try:
        page = int(request.GET.get(query_name, "1"))
    except (TypeError, ValueError):
        return 1
    return max(1, page)


def _positive_tile_width(request):
    """Return a positive tile_width query value, or None to use the template default."""
    try:
        width = int(request.GET.get("tile_width", ""))
    except (TypeError, ValueError):
        return None
    return width if width > 0 else None


def _volume_viewer_requested(request):
    """Return True only for the explicit, default-off volume-viewer key."""
    return request.GET.get("volume_viewer") == "1"


def _deselected_class_ids(request):
    try:
        deselected_ids = json.loads(request.POST.get("deselected_class_ids", ""))
    except (TypeError, json.JSONDecodeError) as error:
        raise ClassSelectionError("Selection data is missing or invalid.") from error
    if not isinstance(deselected_ids, list):
        raise ClassSelectionError("Selection data is missing or invalid.")
    return deselected_ids

def _deselected_mic_ids(request):
    try:
        deselected_ids = json.loads(request.POST.get("deselected_mic_ids", ""))
    except (TypeError, json.JSONDecodeError) as error:
        raise ClassSelectionError("Selection data is missing or invalid.") from error
    if not isinstance(deselected_ids, list):
        raise ClassSelectionError("Selection data is missing or invalid.")
    return deselected_ids


def _selected_states(request):
    """Parse and validate the posted 'selected_states' JSON array of state numbers to keep."""
    try:
        selected_states = json.loads(request.POST.get("selected_states", ""))
    except (TypeError, json.JSONDecodeError) as error:
        raise ClassSelectionError("Selection data is missing or invalid.") from error
    if not isinstance(selected_states, list) or not selected_states:
        raise ClassSelectionError("At least one state must remain selected.")
    states = sorted({
        state for state in selected_states
        if isinstance(state, int) and not isinstance(state, bool) and state >= 1
    })
    if not states:
        raise ClassSelectionError("Selection data is missing or invalid.")
    return states

def _manual_pick_box_coordinates(request):
    """Parse and validate manually-picked box centers posted as a JSON array of {x, y}."""
    try:
        boxes = json.loads(request.POST.get("boxes", ""))
    except (TypeError, json.JSONDecodeError):
        return None
    if not isinstance(boxes, list):
        return None
    coordinates = []
    for box in boxes:
        if not isinstance(box, dict):
            continue
        x, y = box.get("x"), box.get("y")
        if isinstance(x, bool) or isinstance(y, bool):
            continue
        if not isinstance(x, (int, float)) or not isinstance(y, (int, float)):
            continue
        coordinates.append((x, y))
    return coordinates

def _resolve_manualpick_boxfile(batch_job, boxfile_path):
    """Validate a mic's stored boxfile path resolves inside the job directory, without requiring it to already exist."""
    job_dir = batch_job.get_safe_job_dir()
    if job_dir is None or not isinstance(boxfile_path, str) or not boxfile_path.strip():
        return None
    resolved = os.path.realpath(boxfile_path)
    try:
        is_safe = os.path.commonpath((job_dir, resolved)) == job_dir
    except ValueError:
        is_safe = False
    if not is_safe or os.path.islink(boxfile_path):
        return None
    return resolved

def _translate_mic_field_records(raw_micrographs):
    """Translate raw project-file mic field names to the names _mic_selection_element.html expects."""
    micrographs = []
    for record in raw_micrographs:
        if not isinstance(record, dict):
            continue
        micrographs.append({
            "path": record.get("thumb"),
            "ctfimg": record.get("ctfjpg"),
            "dfx": record.get("dfx"),
            "dfy": record.get("dfy"),
            "ctfres": record.get("ctfres"),
            "i": record.get("n"),
            "xdim": record.get("xdim"),
            "ydim": record.get("ydim"),
            "boxes": record.get("boxes"),
        })
    return micrographs

def _log_parts(logtext):
    """Split log text into text/image parts on ">>> JPEG" markers."""
    parts = []
    logpart_str = ""
    for line in logtext.splitlines():
        if ">>> JPEG " in line:
            split_line = line.split()
            if len(split_line) >= 3:
                parts.append({"text": logpart_str})
                parts.append({"image": split_line[2]})
                logpart_str = ""
            else:
                # Keep malformed marker lines in text output instead of crashing.
                logpart_str += line + "\n"
        else:
            logpart_str += line + "\n"
    if logpart_str != "":
        parts.append({"text": logpart_str})
    return parts


def _batch_overview_context(
    batchjob,
    jobmodel,
    volume_viewer_requested=False,
):
    """Shared overview/logs/arguments context for the batch and manual-picker detail views."""
    log_by_name = {entry["name"]: entry for entry in batchjob.get_log_tails()}
    stdout_entry = log_by_name.get("stdout.log", {})
    stderr_entry = log_by_name.get("stderr.log", {})
    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    jobstats = metadata.get("project_metadata", {})
    arguments = _argument_rows(jobmodel)
    show_volume_viewer = (
        volume_viewer_requested
        and jobmodel.prog in BatchJob.VOLUME_VIEWER_PROGRAMS
        and jobmodel.status != "queued"
    )
    stage_volume_outputs = (
        batchjob.get_stage_volume_outputs(jobstats) if show_volume_viewer else []
    )

    return {
        "jobid"  : jobmodel.id,
        "disp"   : jobmodel.disp,
        "name"   : jobmodel.name,
        "desc"   : jobmodel.desc,
        "prog"   : jobmodel.prog,
        "proj"   : jobmodel.dset.proj.name,
        "dset"   : jobmodel.dset.name,
        "args"   : jobmodel.args,
        "created": jobmodel.cdat,
        "folder" : batchjob.get_absdir(),
        "jobstats": jobstats,
        "cls3d_stages": _cls3d_stages(jobstats, stage_volume_outputs),
        "log"    : _log_parts(stdout_entry["text"]) if stdout_entry.get("exists") else [],
        "error"  : stderr_entry.get("text") if stderr_entry.get("exists") else None,
        "arguments": arguments,
        "submitted_argument_count": sum(argument["submitted"] for argument in arguments),
        "volume_viewer_requested": show_volume_viewer,
        "volume_outputs": [
            {key: value for key, value in output.items() if key != "path"}
            for output in stage_volume_outputs
        ],
    }

def _manualpick_overview_context(batchjob, jobmodel):
    """Overview context for the manual-picker view, with a live first page of mic tiles.

    jobstats.micrographs is a snapshot frozen at manualpick job creation (the job
    finishes immediately and never streams further updates), so it never reflects
    picks saved afterwards. Refresh it here with a direct, boxes=True project-file
    read so already-picked coordinates show up on load.
    """
    context = _batch_overview_context(batchjob, jobmodel)
    jobstats = dict(context["jobstats"])

    project_stats = batchjob.getProjectStats()
    mic_stats = project_stats.get("mic") if isinstance(project_stats, dict) else None
    total = mic_stats.get("n") if isinstance(mic_stats, dict) else None
    if isinstance(total, int) and not isinstance(total, bool) and total > 0:
        top = min(_BATCH_MICROGRAPH_PAGE_SIZE, total)
        mic_field_stats = batchjob.getProjectFieldStats("mic", fromp=1, top=top, boxes=True)
        raw_micrographs = mic_field_stats.get("data") if isinstance(mic_field_stats, dict) else None
        if isinstance(raw_micrographs, list):
            jobstats["micrographs"] = _translate_mic_field_records(raw_micrographs)

    context["jobstats"] = jobstats
    return context


@login_required(login_url="/login")
@require_GET
def view_batch(request, jobid):
    """Returns batch view."""
    template = "nice_batch/batchview.html"
    batchjob, jobmodel = _get_accessible_batch_job(request, job_id=jobid, log_context="view_batch")
    if batchjob is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    response = render(
        request,
        template,
        _batch_overview_context(
            batchjob,
            jobmodel,
            volume_viewer_requested=_volume_viewer_requested(request),
        ),
    )

    response.set_cookie(key="selected_project_id", value=jobmodel.dset.proj_id)
    response.set_cookie(key="selected_workspace_id", value=jobmodel.dset_id)
    # Ensure Back renders the checksum-gated workspace instead of returning 204.
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
@require_GET
def view_batch_manual_picker(request, jobid):
    """Returns the manual-picker detail view for a batch job."""
    template = "nice_batch/manual_picker.html"
    batchjob, jobmodel = _get_accessible_batch_job(request, job_id=jobid, log_context="view_batch_manual_picker")
    if batchjob is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")
    
    response = render(request, template, _manualpick_overview_context(batchjob, jobmodel))

    response.set_cookie(key="selected_project_id", value=jobmodel.dset.proj_id)
    response.set_cookie(key="selected_workspace_id", value=jobmodel.dset_id)
    # Ensure Back renders the checksum-gated workspace instead of returning 204.
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
@require_GET
@cache_control(private=True, max_age=300, no_transform=True)
def view_batch_volume_data(request, jobid, volume_name):
    """Stream one declared, owned 3D MRC output for Mol*."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_volume_data",
        job_id=jobid,
    )
    if (
        batch_job is None
        or jobmodel.status == "queued"
        or jobmodel.prog not in BatchJob.VOLUME_VIEWER_PROGRAMS
    ):
        return HttpResponse(status=404)

    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    jobstats = metadata.get("project_metadata", {})
    volume = next(
        (
            output
            for output in (
                batch_job.get_volume_outputs()
                + batch_job.get_stage_volume_outputs(jobstats)
            )
            if output.get("name") == volume_name
        ),
        None,
    )
    if volume is None:
        return HttpResponse(status=404)

    try:
        volume_file = open(volume["path"], "rb")
    except OSError:
        return HttpResponse(status=404)

    response = FileResponse(
        volume_file,
        as_attachment=False,
        filename=volume["name"],
        content_type="application/octet-stream",
    )
    response["X-Content-Type-Options"] = "nosniff"
    return response


@login_required(login_url="/login")
@require_POST
def view_batch_class_2D_selection(request, jobid):
    """Create and launch a new cls2D-deselection batch job from the current selection."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_class_2D_selection",
        job_id=jobid,
    )
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    try:
        deselected_class_ids = _deselected_class_ids(request)
        selection = load_batch_class_selection(
            batch_job.get_result_project_path(),
            jobmodel.dset.proj.dirc,
            jobmodel.id,
        )
        deselected_ids = validate_deselected_class_ids(selection, deselected_class_ids)
    except (ClassSelectionError, OSError, OverflowError, struct.error) as error:
        logger.warning(
            "batch classification 2D selection failed for job %s: %s",
            jobmodel.id,
            error,
        )
        messages.add_message(request, messages.ERROR, f"selection failed: {error}")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    selectionjob = BatchJob()
    project = Project(id=jobmodel.dset.proj.id)
    workspace = Workspace(jobmodel.dset.id)
    if not selectionjob.createClassDeselection(
        project,
        workspace,
        batch_job.get_result_project_path(),
        deselected_ids,
    ):
        logger.warning(
            "batch classification 2D selection job creation failed for job %s",
            jobmodel.id,
        )
        messages.add_message(request, messages.ERROR, "failed to create classification selection job")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    return redirect("nice_lite:workspace")


@login_required(login_url="/login")
@require_POST
def view_batch_class_3D_selection(request, jobid):
    """Create and launch a new cls3D state-selection batch job from the current selection."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_class_3D_selection",
        job_id=jobid,
    )
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    try:
        selected_states = _selected_states(request)
    except ClassSelectionError as error:
        logger.warning(
            "batch classification 3D selection failed for job %s: %s",
            jobmodel.id,
            error,
        )
        messages.add_message(request, messages.ERROR, f"selection failed: {error}")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    selectionjob = BatchJob()
    project = Project(id=jobmodel.dset.proj.id)
    workspace = Workspace(jobmodel.dset.id)
    if not selectionjob.createStateSelection(
        project,
        workspace,
        batch_job.get_result_project_path(),
        selected_states,
    ):
        logger.warning(
            "batch classification 3D selection job creation failed for job %s",
            jobmodel.id,
        )
        messages.add_message(request, messages.ERROR, "failed to create state selection job")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    return redirect("nice_lite:workspace")


@login_required(login_url="/login")
@require_POST
def view_batch_micrograph_selection(request, jobid):
    """Create and launch a new mic-deselection batch job from the current selection."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_micrograph_selection",
        job_id=jobid,
    )
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    try:
        deselected_micrograph_ids = _deselected_mic_ids(request)
    except ClassSelectionError as error:
        logger.warning(
            "batch micrograph selection failed for job %s: %s",
            jobmodel.id,
            error,
        )
        messages.add_message(request, messages.ERROR, f"selection failed: {error}")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    jobstats = metadata.get("project_metadata", {})
    micrographs = jobstats.get("micrographs") or jobstats.get("latest_micrographs") or []
    # jobstats["micrographs"] is capped at 50 entries; nmics is the true project-wide count.
    nmics = jobstats.get("nmics")
    if isinstance(nmics, int) and not isinstance(nmics, bool) and nmics > 0:
        total_micrographs = nmics
    else:
        total_micrographs = len(micrographs)

    result_project = batch_job.get_result_project_path()
    if result_project is None or total_micrographs == 0:
        logger.warning(
            "batch micrograph selection has no source project for job %s",
            jobmodel.id,
        )
        messages.add_message(request, messages.ERROR, "no micrographs available for selection")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    selectionjob = BatchJob()
    project = Project(id=jobmodel.dset.proj.id)
    workspace = Workspace(jobmodel.dset.id)
    if not selectionjob.createMicrographDeselection(
        project,
        workspace,
        result_project,
        deselected_micrograph_ids,
        total_micrographs,
    ):
        logger.warning(
            "batch micrograph selection job creation failed for job %s",
            jobmodel.id,
        )
        messages.add_message(request, messages.ERROR, "failed to create micrograph selection job")
        return redirect("nice_lite:view_batch", jobid=jobmodel.id)

    return redirect("nice_lite:workspace")


@login_required(login_url="/login")
@require_GET
def view_batch_micrographs_page(request, jobid):
    """Return one page of micrograph entries as JSON, read directly from the project file."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_micrographs_page",
        job_id=jobid,
    )
    if batch_job is None:
        return HttpResponse(status=404)

    project_stats = batch_job.getProjectStats()
    mic_stats = project_stats.get("mic") if isinstance(project_stats, dict) else None
    total = mic_stats.get("n") if isinstance(mic_stats, dict) else None
    if not isinstance(total, int) or isinstance(total, bool) or total <= 0:
        return JsonResponse({
            "html": "",
            "total": 0,
            "page": 1,
            "pages": 0,
            "first_micrograph": 0,
            "last_micrograph": 0,
        })

    pages = math.ceil(total / _BATCH_MICROGRAPH_PAGE_SIZE)
    page = min(_positive_page_number(request, "page"), pages)
    fromp = ((page - 1) * _BATCH_MICROGRAPH_PAGE_SIZE) + 1
    top = min(page * _BATCH_MICROGRAPH_PAGE_SIZE, total)

    sort_key = request.GET.get("sort") or None
    sort_asc = request.GET.get("asc", "1") != "0"

    mic_field_stats = batch_job.getProjectFieldStats(
        "mic", fromp=fromp, top=top, sortkey=sort_key, sortasc=sort_asc, boxes=True
    )
    raw_micrographs = mic_field_stats.get("data") if isinstance(mic_field_stats, dict) else None
    if not isinstance(raw_micrographs, list):
        raw_micrographs = []

    micrographs = _translate_mic_field_records(raw_micrographs)

    tiles_html = render_to_string(
        "includes/_mic_page_tiles.html",
        {
            "micrographs": micrographs,
            "jobstats": jobmodel.master_stats.get("project_metadata", {}) if isinstance(jobmodel.master_stats, dict) else {},
            "selectable": request.GET.get("selectable") == "1",
            "hide_pspec": request.GET.get("hide_pspec") == "1",
            "tile_width": _positive_tile_width(request),
        },
        request=request,
    )

    return JsonResponse({
        "html": tiles_html,
        "total": total,
        "page": page,
        "pages": pages,
        "first_micrograph": fromp,
        "last_micrograph": top,
    })


@login_required(login_url="/login")
@require_POST
def view_batch_save_manual_pick_boxes(request, jobid):
    """Overwrite one micrograph's box file with its manually picked coordinates."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_save_manual_pick_boxes",
        job_id=jobid,
    )
    if batch_job is None:
        return HttpResponse(status=404)

    try:
        mic_index = int(request.POST.get("mic", ""))
    except (TypeError, ValueError):
        return HttpResponseBadRequest("invalid micrograph index")
    if mic_index <= 0:
        return HttpResponseBadRequest("invalid micrograph index")

    coordinates = _manual_pick_box_coordinates(request)
    if coordinates is None:
        return HttpResponseBadRequest("invalid box data")

    mic_field_stats = batch_job.getProjectFieldStats("mic", fromp=mic_index, top=mic_index)
    records = mic_field_stats.get("data") if isinstance(mic_field_stats, dict) else None
    # print_project_field's fromto range is only applied when fromp < top, so a single-index
    # query (fromp == top) silently returns every mic; match on "n" instead of assuming records[0].
    record = next(
        (r for r in records if isinstance(r, dict) and r.get("n") == mic_index),
        None,
    ) if isinstance(records, list) else None
    boxfile_path = _resolve_manualpick_boxfile(
        batch_job, record.get("boxfile") if isinstance(record, dict) else None
    )
    if boxfile_path is None:
        return HttpResponse(status=404)

    half_box = _MANUALPICK_BOX_SIZE_PX / 2
    try:
        with open(boxfile_path, "w", encoding="utf-8") as box_file:
            box_file.writelines(
                f"{round(x - half_box)} {round(y - half_box)} "
                f"{_MANUALPICK_BOX_SIZE_PX} {_MANUALPICK_BOX_SIZE_PX} -3\n"
                for x, y in coordinates
            )
    except OSError:
        logger.warning(
            "failed to write manual-pick box file for job %s mic %s",
            jobmodel.id,
            mic_index,
        )
        return HttpResponse(status=500)

    return JsonResponse({"saved": True, "count": len(coordinates)})


@login_required(login_url="/login")
@require_POST
def view_batch_terminate(request):
    """Terminates batch job and refreshes page by redirecting to workspace view."""
    batch_job, _jobmodel = _get_accessible_batch_job(request, log_context="terminate_batch")
    if batch_job is None:
        return redirect("nice_lite:workspace")
    batch_job.terminate()
    return redirect("nice_lite:workspace")

@login_required(login_url="/login")
@require_POST
def view_batch_stop(request):
    """Stop an owned running batch job and return to its workspace."""
    batch_job, jobmodel = _get_accessible_batch_job(request, "stop_batch")
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
    elif jobmodel.status != "running":
        messages.add_message(request, messages.ERROR, "batch job is not running")
    elif batch_job.stop():
        messages.add_message(request, messages.INFO, "batch job stopped successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to stop batch job")
    return redirect("nice_lite:workspace")


@login_required(login_url="/login")
@require_POST
def view_batch_mark_finished(request):
    """Mark an owned queued or failed batch job as finished."""
    batch_job, jobmodel = _get_accessible_batch_job(request, "mark_batch_finished")
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
    elif request.POST.get("mark_finished") != "1":
        messages.add_message(request, messages.ERROR, "batch finish confirmation is missing")
    elif jobmodel.status not in BatchJob.MANUALLY_FINISHABLE_STATUSES:
        messages.add_message(request, messages.ERROR, "batch job cannot be marked finished")
    elif batch_job.markComplete(None, None):
        messages.add_message(request, messages.INFO, "batch job marked finished successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to mark batch job finished")
    return redirect("nice_lite:workspace")


@login_required(login_url="/login")
@require_POST
def view_batch_rerun(request):
    """Open the job builder with an owned batch job selected, regardless of status."""
    batch_job, jobmodel = _get_accessible_batch_job(request, "rerun_batch")
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    if jobmodel.pckg not in ("simple", "single"):
        messages.add_message(request, messages.ERROR, "batch job cannot be rerun")
        return redirect("nice_lite:workspace")

    return redirect(reverse(
        "nice_lite:new_stream",
        query={"selected_job_id": jobmodel.id},
    ))


@login_required(login_url="/login")
@require_POST
def view_batch_delete(request):
    """Permanently delete an owned deletable batch job and return to its workspace."""
    batch_job, jobmodel = _get_accessible_batch_job(request, "delete_batch")
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
    elif jobmodel.status not in BatchJob.DELETABLE_STATUSES:
        messages.add_message(request, messages.ERROR, "batch job is not complete")
    elif jobmodel.status == "queued" and not batch_job.queued_job_can_delete():
        messages.add_message(request, messages.ERROR, "queued batch job is still active or cannot be verified")
    else:
        project = Project(id=jobmodel.dset.proj_id)
        workspace = Workspace(jobmodel.dset_id)
        if batch_job.delete(project, workspace):
            messages.add_message(request, messages.INFO, "batch job permanently deleted successfully")
        else:
            messages.add_message(request, messages.ERROR, "failed to delete batch job")
    return redirect("nice_lite:workspace")
