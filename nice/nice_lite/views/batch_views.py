"""Lifecycle actions for classic SIMPLE and SINGLE batch jobs."""

import json
import logging
import os
import struct

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core import signing
from django.core.paginator import Paginator
from django.http import HttpResponse
from django.shortcuts import redirect, render
from django.urls import reverse
from django.views.decorators.cache import cache_control
from django.views.decorators.http import require_GET, require_POST

from ..data_structures.batchjob import BatchJob
from ..data_structures.class_selection import (
    ClassSelectionError,
    batch_class_selection_available,
    deselected_class_ids,
    load_batch_class_selection,
)
from ..data_structures.mrc import render_mrc_particle_png
from ..data_structures.movie import movie_preview_supported, read_movie_dimensions
from ..data_structures.project import Project
from ..data_structures.simple import SIMPLEBatch, SIMPLEProjFile
from ..data_structures.workspace import Workspace
from ..helpers import clear_checksum_cookies, get_job_id


logger = logging.getLogger(__name__)


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
_BATCH_PARTICLE_PAGE_SIZE = 40
_BATCH_MOVIE_PAGE_SIZE = 40
_BATCH_MOVIE_THUMBNAIL_SALT = "nice-lite.batch-movie-thumbnail"


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


def _empty_movie_page():
    """Return the stable empty shape used by the import movie gallery."""
    return {
        "movies": [],
        "total": 0,
        "page": 1,
        "pages": 0,
        "page_numbers": [],
        "ellipsis": Paginator.ELLIPSIS,
        "has_previous": False,
        "has_next": False,
        "first_movie": 0,
        "last_movie": 0,
    }


def _has_positive_dimensions(item, width_key="width", height_key="height"):
    """Return whether a rendered output item has usable pixel dimensions."""
    if not isinstance(item, dict):
        return False
    return all(
        isinstance(item.get(key), (int, float))
        and not isinstance(item.get(key), bool)
        and item.get(key) > 0
        for key in (width_key, height_key)
    )


def _movie_thumbnail_token(job_id, movie_path):
    """Sign an imported movie path so it cannot be replaced in the URL."""
    return signing.Signer(salt=_BATCH_MOVIE_THUMBNAIL_SALT).sign_object(
        {
            "job_id": job_id,
            "path": movie_path,
            "version": BatchJob.MOVIE_THUMBNAIL_CACHE_VERSION,
        },
        compress=True,
    )


def _import_movie_page(job_id, batch_job, page=1, page_size=40):
    """Return one metadata-only page from an import job's submitted source."""
    empty_page = _empty_movie_page()
    movie_paths = batch_job.get_import_movie_paths()
    if not isinstance(movie_paths, list):
        return empty_page
    total = len(movie_paths)
    if (
        total < 1
        or not isinstance(page_size, int)
        or isinstance(page_size, bool)
        or page_size < 1
    ):
        return empty_page

    page_size = min(page_size, 100)
    paginator = Paginator(movie_paths, page_size)
    page_obj = paginator.get_page(page)
    first_movie = page_obj.start_index()
    last_movie = page_obj.end_index()

    movies = []
    for offset, movie_path in enumerate(page_obj.object_list):
        preview_available = movie_preview_supported(movie_path)
        dimensions = read_movie_dimensions(movie_path)
        movies.append({
            "number": first_movie + offset,
            "name": os.path.basename(movie_path) or movie_path,
            "width": dimensions[0] if dimensions is not None else None,
            "height": dimensions[1] if dimensions is not None else None,
            "preview_available": preview_available,
            "token": (
                _movie_thumbnail_token(job_id, movie_path)
                if preview_available else ""
            ),
        })

    return {
        "movies": movies,
        "total": total,
        "page": page_obj.number,
        "pages": paginator.num_pages,
        "page_numbers": list(paginator.get_elided_page_range(
            page_obj.number,
            on_each_side=2,
            on_ends=1,
        )),
        "ellipsis": Paginator.ELLIPSIS,
        "has_previous": page_obj.has_previous(),
        "previous_page": (
            page_obj.previous_page_number() if page_obj.has_previous() else None
        ),
        "has_next": page_obj.has_next(),
        "next_page": page_obj.next_page_number() if page_obj.has_next() else None,
        "first_movie": first_movie,
        "last_movie": last_movie,
    }


def _batch_detail_context(
    batch_job,
    jobmodel,
    particle_page=1,
    movie_page=1,
    class_selector_requested=False,
):
    """Assemble validated batch metadata, logs, artifacts, and project summary."""
    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    result_project = batch_job.get_result_project_path()
    project_stats = {}
    project_reader = None
    pick_micrographs = []
    if jobmodel.status == "finished" and result_project is not None:
        project_reader = SIMPLEProjFile(result_project)
        project_stats = project_reader.getGlobalStats()
    cls2d_stats = project_stats.get("cls2D") if isinstance(project_stats, dict) else None
    reported_cls2d_available = (
        isinstance(cls2d_stats, dict)
        and isinstance(cls2d_stats.get("n"), int)
        and not isinstance(cls2d_stats.get("n"), bool)
        and cls2d_stats["n"] > 0
    )
    class_selector_available = (
        jobmodel.status == "finished"
        and result_project is not None
        and (
            reported_cls2d_available
            or batch_class_selection_available(
                result_project,
                jobmodel.dset.proj.dirc,
            )
        )
    )
    batch_class_selector = None
    class_selector_error = ""
    if class_selector_requested and class_selector_available:
        try:
            selection = load_batch_class_selection(
                result_project,
                jobmodel.dset.proj.dirc,
                jobmodel.id,
            )
            batch_class_selector = {
                "classes": selection.classes,
                "class_count": len(selection.classes),
                "stack_name": selection.stack_name,
                "width": selection.width,
                "height": selection.height,
                "initial_selected_class_ids": selection.initial_selected_class_ids,
                "browser_data": selection.browser_data(),
            }
        except (ClassSelectionError, OSError, OverflowError, struct.error) as error:
            logger.warning(
                "batch class selector unavailable for job %s: %s",
                jobmodel.id,
                error,
            )
            class_selector_error = str(error)
    if jobmodel.status == "finished" and metadata.get("program") == "pick":
        pick_micrographs = batch_job.get_pick_micrograph_previews(
            max_previews=_BATCH_PICK_PREVIEW_LIMIT,
            max_coordinates=_BATCH_PICK_COORDINATE_LIMIT,
        )
        if not pick_micrographs and result_project is not None:
            pick_micrographs = _pick_micrograph_previews(project_reader, project_stats)

    artifact_summary = batch_job.get_artifact_summary()
    artifact_counts = list(artifact_summary.get("counts", []))
    artifact_images = artifact_summary.get("images", [])
    class_selector_replaces_artifact_previews = (
        metadata.get("program") == "abinitio2D"
        and batch_class_selector is not None
    )
    particle_stack_page = {}
    if metadata.get("program") in BatchJob.PARTICLE_STACK_PROGRAMS:
        particle_stack_page = batch_job.get_particle_stack_page(
            page=particle_page,
            page_size=_BATCH_PARTICLE_PAGE_SIZE,
        )
        if particle_stack_page.get("stacks"):
            artifact_counts.append({
                "extension": "MRCS",
                "count": len(particle_stack_page["stacks"]),
            })
    import_movie_page = {}
    if (
        metadata.get("program") == "import_movies"
        and jobmodel.status == "finished"
    ):
        import_movie_page = _import_movie_page(
            jobmodel.id,
            batch_job,
            page=movie_page,
            page_size=_BATCH_MOVIE_PAGE_SIZE,
        )
    output_dimensions_available = (
        _has_positive_dimensions(batch_class_selector)
        or any(
            _has_positive_dimensions(micrograph, "xdim", "ydim")
            for micrograph in pick_micrographs
        )
        or any(
            _has_positive_dimensions(preview)
            for image in artifact_images
            for preview in image.get("previews", [])
        )
        or any(
            _has_positive_dimensions(particle)
            for particle in particle_stack_page.get("particles", [])
        )
        or any(
            _has_positive_dimensions(movie)
            for movie in import_movie_page.get("movies", [])
        )
    )
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
        "class_selector_available": class_selector_available,
        "class_selector_requested": class_selector_requested,
        "batch_class_selector": batch_class_selector,
        "class_selector_error": class_selector_error,
        "class_selector_replaces_artifact_previews": (
            class_selector_replaces_artifact_previews
        ),
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
        "artifact_counts": artifact_counts,
        "artifact_images": artifact_images,
        "particle_stack_page": particle_stack_page,
        "import_movie_page": import_movie_page,
        "output_dimensions_available": output_dimensions_available,
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
        "auto_refresh": jobmodel.status in ("queued", "running"),
    }


def _positive_page_number(request, query_name):
    """Return a positive page query value, defaulting safely to one."""
    try:
        page = int(request.GET.get(query_name, "1"))
    except (TypeError, ValueError):
        return 1
    return max(1, page)


def _class_selector_requested(request):
    """Return True only for the explicit, default-off batch selector key."""
    return request.GET.get("class_selector") == "1"


def _class_selector_redirect(job_id):
    batch_url = reverse("nice_lite:view_batch", args=(job_id,))
    return f"{batch_url}?class_selector=1#batch_class_selector"


def _selected_class_ids(request):
    try:
        selected_ids = json.loads(request.POST.get("selected_class_ids", ""))
    except (TypeError, json.JSONDecodeError) as error:
        raise ClassSelectionError("Selection data is missing or invalid.") from error
    if not isinstance(selected_ids, list):
        raise ClassSelectionError("Selection data is missing or invalid.")
    return selected_ids


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
        _batch_detail_context(
            batch_job,
            jobmodel,
            particle_page=_positive_page_number(request, "particle_page"),
            movie_page=_positive_page_number(request, "movie_page"),
            class_selector_requested=_class_selector_requested(request),
        ),
    )
    response.set_cookie(key="selected_project_id", value=jobmodel.dset.proj_id)
    response.set_cookie(key="selected_workspace_id", value=jobmodel.dset_id)
    # Ensure Back renders the checksum-gated workspace instead of returning 204.
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
@require_GET
@cache_control(private=True, max_age=300, no_transform=True)
def view_batch_class_thumbnail(request, jobid, stack_index):
    """Render one class average from an owned finished batch result."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_class_thumbnail",
        job_id=jobid,
    )
    if batch_job is None or jobmodel.status != "finished":
        return HttpResponse(status=404)

    result_project = batch_job.get_result_project_path()
    try:
        selection = load_batch_class_selection(
            result_project,
            jobmodel.dset.proj.dirc,
            jobmodel.id,
        )
    except (ClassSelectionError, OSError, OverflowError, struct.error):
        return HttpResponse(status=404)
    if stack_index not in {
        entry["stack_index"] for entry in selection.classes
    }:
        return HttpResponse(status=404)

    thumbnail = render_mrc_particle_png(
        selection.stack_path,
        stack_index,
        max_size=512,
    )
    if thumbnail is None:
        return HttpResponse(status=404)
    response = HttpResponse(thumbnail, content_type="image/png")
    response["X-Content-Type-Options"] = "nosniff"
    return response


@login_required(login_url="/login")
@require_POST
def view_batch_class_deselection_export(request, jobid):
    """Download a validated, one-based deselection list for a batch result."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_class_deselection_export",
        job_id=jobid,
    )
    if batch_job is None or jobmodel.status != "finished":
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    try:
        selection = load_batch_class_selection(
            batch_job.get_result_project_path(),
            jobmodel.dset.proj.dirc,
            jobmodel.id,
        )
        deselected_ids = deselected_class_ids(
            selection,
            _selected_class_ids(request),
        )
    except (ClassSelectionError, OSError, OverflowError, struct.error) as error:
        logger.warning(
            "batch class deselection export failed for job %s: %s",
            jobmodel.id,
            error,
        )
        messages.add_message(request, messages.ERROR, f"selection export failed: {error}")
        return redirect(_class_selector_redirect(jobmodel.id))

    response = HttpResponse(
        "".join(f"{class_id}\n" for class_id in deselected_ids),
        content_type="text/plain; charset=utf-8",
    )
    response["Content-Disposition"] = (
        f'attachment; filename="batch_{jobmodel.id}_deselected_classes.txt"'
    )
    response["X-Content-Type-Options"] = "nosniff"
    return response


@login_required(login_url="/login")
@require_GET
@cache_control(private=True, max_age=300, no_transform=True)
def view_batch_particle_thumbnail(request, jobid, stack_name, particle_index):
    """Render one owned extract particle on demand without writing a thumbnail."""
    batch_job, _ = _get_accessible_batch_job(
        request,
        "view_batch_particle_thumbnail",
        job_id=jobid,
    )
    if batch_job is None:
        return HttpResponse(status=404)

    thumbnail = batch_job.get_particle_thumbnail(stack_name, particle_index)
    if thumbnail is None:
        return HttpResponse(status=404)
    return HttpResponse(thumbnail, content_type="image/png")


@login_required(login_url="/login")
@require_GET
@cache_control(private=True, max_age=300, no_transform=True)
def view_batch_movie_thumbnail(request, jobid, token):
    """Return one signed import movie from the on-demand WebP cache."""
    batch_job, _ = _get_accessible_batch_job(
        request,
        "view_batch_movie_thumbnail",
        job_id=jobid,
    )
    if batch_job is None:
        return HttpResponse(status=404)

    try:
        payload = signing.Signer(
            salt=_BATCH_MOVIE_THUMBNAIL_SALT,
        ).unsign_object(token)
    except signing.BadSignature:
        return HttpResponse(status=404)
    if (
        not isinstance(payload, dict)
        or payload.get("job_id") != jobid
        or not isinstance(payload.get("path"), str)
        or payload.get("version") != BatchJob.MOVIE_THUMBNAIL_CACHE_VERSION
    ):
        return HttpResponse(status=404)

    thumbnail = batch_job.get_import_movie_thumbnail(payload["path"])
    if thumbnail is None:
        return HttpResponse(status=404)
    return HttpResponse(thumbnail, content_type="image/webp")


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
