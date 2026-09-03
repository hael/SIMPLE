"""Lifecycle actions for classic SIMPLE and SINGLE batch jobs."""

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect, render
from django.views.decorators.http import require_GET, require_POST

from ..data_structures.batchjob import BatchJob
from ..data_structures.project import Project
from ..data_structures.simple import SIMPLEBatch, SIMPLEProjFile
from ..data_structures.workspace import Workspace
from ..features import batch_job_controls_enabled
from ..helpers import clear_checksum_cookies, get_job_id, print_error


_PROJECT_SECTION_LABELS = {
    "mic": "micrographs",
    "stk": "stacks",
    "ptcl2D": "2D particles",
    "cls2D": "2D classes",
    "ptcl3D": "3D particles",
    "cls3D": "3D classes",
}
_BATCH_LAUNCHER_KEYS = {"prg", "projfile", "mkdir", "niceprocid", "niceserver"}
_BATCH_PICK_PREVIEW_LIMIT = 20
_BATCH_PICK_COORDINATE_LIMIT = 1500


def _get_accessible_batch_job(request, log_context, job_id=None):
    """Return an owned batch job and model selected by the request."""
    resolved_job_id = job_id if job_id is not None else get_job_id(request)
    if (
        not isinstance(resolved_job_id, int)
        or isinstance(resolved_job_id, bool)
        or resolved_job_id <= 0
    ):
        print_error(f"{log_context}: missing job id")
        return None, None

    batch_job = BatchJob(id=resolved_job_id)
    jobmodel = batch_job.get_jobmodel()
    workspace = getattr(jobmodel, "dset", None)
    owner = (getattr(workspace, "user", "") or "").strip()
    if jobmodel is None or owner != request.user.username:
        print_error(f"{log_context}: access denied for job {resolved_job_id}")
        return None, None
    return batch_job, jobmodel


def _source_label(jobmodel):
    """Return a concise display label for the persisted batch input source."""
    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    source = metadata.get("source")
    if not isinstance(source, dict):
        if metadata.get("program") == "new_project":
            return "none (creates a project)"
        return "workspace project"

    source_type = source.get("type")
    if source_type == "batch_job":
        return f"batch job {source.get('batch_job_id', '')}".strip()
    if source_type == "stream_snapshot":
        stream_id = source.get("stream_job_id", "")
        set_id = source.get("particle_set_id", "")
        return f"stream {stream_id}, particle set {set_id}".strip(", ")
    if source_type == "project_file":
        return source.get("filename") or "selected project file"
    return str(source_type or "workspace project").replace("_", " ")


def _project_sections(project_stats):
    """Normalize project orientation counts for the batch-detail template."""
    if not isinstance(project_stats, dict):
        return []

    sections = []
    for key, value in project_stats.items():
        if not isinstance(value, dict) or "n" not in value:
            continue
        sections.append({
            "key": key,
            "label": _PROJECT_SECTION_LABELS.get(key, key.replace("_", " ")),
            "count": value.get("n", 0),
            "info": value.get("info", ""),
        })
    return sections


def _pick_micrograph_previews(project_reader, project_stats):
    """Return recent picked micrographs in the stream preview data shape."""
    mic_stats = project_stats.get("mic") if isinstance(project_stats, dict) else None
    micrograph_count = mic_stats.get("n") if isinstance(mic_stats, dict) else None
    if (
        not isinstance(micrograph_count, int)
        or isinstance(micrograph_count, bool)
        or micrograph_count <= 0
    ):
        return []

    from_micrograph = max(1, micrograph_count - _BATCH_PICK_PREVIEW_LIMIT + 1)
    field_stats = project_reader.getFieldStats(
        "mic",
        fromp=from_micrograph,
        top=micrograph_count,
        boxes=True,
    )
    raw_micrographs = field_stats.get("data") if isinstance(field_stats, dict) else None
    if not isinstance(raw_micrographs, list):
        return []

    previews = []
    for micrograph in raw_micrographs:
        if not isinstance(micrograph, dict):
            continue
        path = micrograph.get("thumb")
        xdim = micrograph.get("xdim")
        ydim = micrograph.get("ydim")
        if (
            not isinstance(path, str)
            or not path
            or not isinstance(xdim, (int, float))
            or isinstance(xdim, bool)
            or xdim <= 0
            or not isinstance(ydim, (int, float))
            or isinstance(ydim, bool)
            or ydim <= 0
        ):
            continue

        boxes = []
        raw_boxes = micrograph.get("boxes")
        if isinstance(raw_boxes, list):
            for box in raw_boxes[:_BATCH_PICK_COORDINATE_LIMIT]:
                if not isinstance(box, dict):
                    continue
                xcoord = box.get("x")
                ycoord = box.get("y")
                if (
                    isinstance(xcoord, (int, float))
                    and not isinstance(xcoord, bool)
                    and isinstance(ycoord, (int, float))
                    and not isinstance(ycoord, bool)
                ):
                    normalized_box = {"x": xcoord, "y": ycoord}
                    width = box.get("width")
                    height = box.get("height")
                    if (
                        isinstance(width, (int, float))
                        and not isinstance(width, bool)
                        and width > 0
                        and isinstance(height, (int, float))
                        and not isinstance(height, bool)
                        and height > 0
                    ):
                        normalized_box.update({"width": width, "height": height})
                    boxes.append(normalized_box)

        previews.append({
            "path": path,
            "number": micrograph.get("n"),
            "xdim": xdim,
            "ydim": ydim,
            "boxes": boxes,
        })
    return previews


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
            })
    return arguments


def _argument_rows(jobmodel, metadata):
    """Build submitted/default/unset rows from saved args and commander UI."""
    raw_saved_args = jobmodel.args if isinstance(jobmodel.args, dict) else {}
    saved_args = {str(key): value for key, value in raw_saved_args.items()}
    definitions = _commander_arguments(
        metadata.get("package"),
        metadata.get("program"),
    )
    arguments = []
    known_keys = set()
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
            })
    return arguments


def _batch_detail_context(batch_job, jobmodel):
    """Assemble validated batch metadata, logs, artifacts, and project summary."""
    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    result_project = batch_job.get_result_project_path()
    project_stats = {}
    pick_micrographs = []
    if jobmodel.status == "finished" and result_project is not None:
        project_reader = SIMPLEProjFile(result_project)
        project_stats = project_reader.getGlobalStats()
    if jobmodel.status == "finished" and metadata.get("program") == "pick":
        pick_micrographs = batch_job.get_pick_micrograph_previews(
            max_previews=_BATCH_PICK_PREVIEW_LIMIT,
            max_coordinates=_BATCH_PICK_COORDINATE_LIMIT,
        )
        if not pick_micrographs and result_project is not None:
            pick_micrographs = _pick_micrograph_previews(project_reader, project_stats)

    artifact_summary = batch_job.get_artifact_summary()
    artifact_images = artifact_summary.get("images", [])
    arguments = _argument_rows(jobmodel, metadata)
    return {
        "jobid": jobmodel.id,
        "disp": jobmodel.disp,
        "name": jobmodel.name,
        "desc": jobmodel.desc,
        "status": jobmodel.status,
        "created": jobmodel.cdat,
        "package": metadata.get("package", ""),
        "program": metadata.get("program", ""),
        "project": jobmodel.dset.proj.name,
        "workspace": jobmodel.dset.name,
        "job_dir": batch_job.get_safe_job_dir() or "unavailable",
        "source_label": _source_label(jobmodel),
        "arguments": arguments,
        "submitted_argument_count": sum(argument["submitted"] for argument in arguments),
        "result_project": result_project,
        "project_sections": _project_sections(project_stats),
        "project_summary_available": bool(project_stats),
        "pick_micrographs": pick_micrographs,
        "pick_box_overlay_available": any(
            isinstance(box.get("width"), (int, float))
            and not isinstance(box.get("width"), bool)
            and box.get("width") > 0
            and isinstance(box.get("height"), (int, float))
            and not isinstance(box.get("height"), bool)
            and box.get("height") > 0
            for micrograph in pick_micrographs
            if isinstance(micrograph, dict)
            for box in micrograph.get("boxes", [])
            if isinstance(box, dict)
        ),
        "logs": batch_job.get_log_tails(),
        "artifact_counts": artifact_summary.get("counts", []),
        "artifact_images": artifact_images,
        "motion_artifact_toggle_available": (
            metadata.get("program") == "motion_correct"
            and any(
                preview.get("visibility_group") == "motion"
                for image in artifact_images
                for preview in image.get("previews", [])
            )
        ),
        "ctf_artifact_micrograph_toggle_available": (
            metadata.get("program") == "ctf_estimate"
            and any(
                preview.get("hidden_by_default")
                for image in artifact_images
                for preview in image.get("previews", [])
            )
        ),
        "batch_job_controls_enabled": batch_job_controls_enabled(),
        "auto_refresh": jobmodel.status in ("queued", "running"),
    }


def _controls_available(request):
    """Reject lifecycle writes while the opt-in feature is disabled."""
    if batch_job_controls_enabled():
        return True
    messages.add_message(request, messages.ERROR, "batch job controls are disabled")
    return False


@login_required(login_url="/login")
@require_GET
def view_batch(request, jobid):
    """Render the owned SIMPLE/SINGLE batch detail page."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch",
        job_id=jobid,
    )
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    response = render(
        request,
        "nice_classic/batchview.html",
        _batch_detail_context(batch_job, jobmodel),
    )
    response.set_cookie(key="selected_project_id", value=jobmodel.dset.proj_id)
    response.set_cookie(key="selected_workspace_id", value=jobmodel.dset_id)
    # Ensure Back renders the checksum-gated workspace instead of returning 204.
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
@require_POST
def view_batch_stop(request):
    """Stop an owned running batch job and return to its workspace."""
    if not _controls_available(request):
        return redirect("nice_lite:workspace")

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
def view_batch_delete(request):
    """Permanently delete an owned deletable batch job and return to its workspace."""
    if not _controls_available(request):
        return redirect("nice_lite:workspace")

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
