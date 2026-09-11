"""Lifecycle actions for classic SIMPLE and SINGLE batch jobs."""

import json
import logging
import math
import os
import struct
import tempfile

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core import signing
from django.core.paginator import Paginator
from django.http import FileResponse, HttpResponse
from django.shortcuts import redirect, render
from django.urls import reverse
from django.views.decorators.cache import cache_control
from django.views.decorators.http import require_GET, require_POST

from ..data_structures.batchjob import BatchJob
from ..data_structures.class_selection import (
    ClassSelectionError,
    SIMPLEProjectFileReader,
    batch_class_selection_available,
    class_selection_flags,
    validate_deselected_class_ids,
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
_BATCH_CLASS_SELECTION_FILENAME = "class_selection.txt"


def _positive_finite_number(value):
    """Return a normalized positive finite numeric value, otherwise None."""
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or number <= 0:
        return None
    rounded = round(number, 6)
    return int(rounded) if rounded.is_integer() else rounded


def _project_sampling_distance(project_path):
    """Read the output sampling distance directly from a SIMPLE project."""
    if not isinstance(project_path, (str, os.PathLike)):
        return None
    try:
        reader = SIMPLEProjectFileReader(project_path)
        for oritype in ("out", "stk", "mic"):
            for record in reader.read_records(oritype):
                sampling_distance = _positive_finite_number(
                    record.get("smpd") if isinstance(record, dict) else None
                )
                if sampling_distance is not None:
                    return sampling_distance
    except (ClassSelectionError, OSError, OverflowError, struct.error):
        return None
    return None


def _project_picked_particle_count(project_path):
    """Sum per-micrograph picked-particle counts from a SIMPLE project."""
    if not isinstance(project_path, (str, os.PathLike)):
        return None
    try:
        records = SIMPLEProjectFileReader(project_path).read_records("mic")
    except (ClassSelectionError, OSError, OverflowError, struct.error):
        return None
    if not records:
        return None

    total = 0
    for record in records:
        value = record.get("nptcls") if isinstance(record, dict) else None
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return None
        if isinstance(value, int):
            if value < 0:
                return None
            total += value
            continue
        if not math.isfinite(value) or value < 0 or not value.is_integer():
            return None
        total += int(value)
    return total


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
        if jobmodel.prog == "new_project":
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


def _submitted_mask_diameter(jobmodel):
    """Return a positive finite submitted mask diameter in angstroms."""
    raw_args = jobmodel.args if isinstance(jobmodel.args, dict) else {}
    return _positive_finite_number(raw_args.get("mskdiam"))


def _class_overlay_settings(
    jobmodel,
    selection,
    project_sampling_distance=None,
):
    """Return class-overlay sizes derived from SIMPLE project sampling."""
    sampling_distance = getattr(selection, "sampling_distance", None)
    if (
        isinstance(sampling_distance, bool)
        or not isinstance(sampling_distance, (int, float))
        or not math.isfinite(sampling_distance)
        or sampling_distance <= 0
    ):
        sampling_distance = _positive_finite_number(project_sampling_distance)

    mask_diameter_angstroms = _submitted_mask_diameter(jobmodel)
    settings = {
        "sampling_distance": sampling_distance,
        "overlay_size": None,
        "overlay_display_size": None,
        "overlay_unit": "pixels",
    }
    if sampling_distance is None or mask_diameter_angstroms is None:
        return settings

    overlay_size = round(mask_diameter_angstroms / sampling_distance, 6)
    settings.update({
        "overlay_size": overlay_size,
        "overlay_display_size": mask_diameter_angstroms,
        "overlay_unit": "angstroms",
    })
    return settings


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
    volume_viewer_requested=False,
):
    """Assemble validated batch metadata, logs, artifacts, and project summary."""
    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    result_project = batch_job.get_result_project_path()
    project_sampling_distance = _project_sampling_distance(result_project)
    project_stats = {}
    project_reader = None
    pick_micrographs = []
    pick_particle_count = None
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
                **_class_overlay_settings(
                    jobmodel,
                    selection,
                    project_sampling_distance=project_sampling_distance,
                ),
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
    if jobmodel.status == "finished" and jobmodel.prog == "pick":
        pick_particle_count = _project_picked_particle_count(result_project)
        pick_micrographs = batch_job.get_pick_micrograph_previews(
            max_previews=_BATCH_PICK_PREVIEW_LIMIT,
            max_coordinates=_BATCH_PICK_COORDINATE_LIMIT,
        )

    artifact_summary = batch_job.get_artifact_summary()
    artifact_counts = list(artifact_summary.get("counts", []))
    artifact_images = (
        []
        if jobmodel.prog == "pick"
        else artifact_summary.get("images", [])
    )
    class_selector_replaces_artifact_previews = (
        jobmodel.prog == "abinitio2D"
        and batch_class_selector is not None
    )
    volume_outputs = []
    if jobmodel.status == "finished" and jobmodel.prog == "abinitio3D":
        volume_outputs = batch_job.get_volume_outputs()
    public_volume_outputs = [
        {
            key: value
            for key, value in volume.items()
            if key != "path"
        }
        for volume in volume_outputs
    ]
    particle_stack_page = {}
    if jobmodel.prog in BatchJob.MRC_STACK_PREVIEW_PROGRAMS:
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
        jobmodel.prog == "import_movies"
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
            _has_positive_dimensions(micrograph)
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
    arguments = _argument_rows(jobmodel)
    return {
        "jobid": jobmodel.id,
        "disp": jobmodel.disp,
        "name": jobmodel.name,
        "desc": jobmodel.desc,
        "status": jobmodel.status,
        "created": jobmodel.cdat,
        "package": jobmodel.pckg,
        "program": jobmodel.prog,
        "project": jobmodel.dset.proj.name,
        "workspace": jobmodel.dset.name,
        "job_dir": batch_job.get_safe_job_dir() or "unavailable",
        "source_label": _source_label(jobmodel),
        "arguments": arguments,
        "submitted_argument_count": sum(argument["submitted"] for argument in arguments),
        "result_project": result_project,
        "project_sampling_distance": project_sampling_distance,
        "project_sections": _project_sections(project_stats),
        "project_summary_available": bool(project_stats),
        "class_selector_available": class_selector_available,
        "class_selector_requested": class_selector_requested,
        "batch_class_selector": batch_class_selector,
        "class_selector_error": class_selector_error,
        "class_selector_replaces_artifact_previews": (
            class_selector_replaces_artifact_previews
        ),
        "volume_viewer_available": bool(public_volume_outputs),
        "volume_viewer_requested": volume_viewer_requested,
        "batch_volume_viewer": (
            public_volume_outputs if volume_viewer_requested else []
        ),
        "pick_micrographs": pick_micrographs,
        "pick_particle_count": pick_particle_count,
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
            jobmodel.prog == "motion_correct"
            and any(
                preview.get("visibility_group") == "motion"
                for image in artifact_images
                for preview in image.get("previews", [])
            )
        ),
        "ctf_artifact_micrograph_toggle_available": (
            jobmodel.prog == "ctf_estimate"
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


def _volume_viewer_requested(request):
    """Return True only for the explicit, default-off volume-viewer key."""
    return request.GET.get("volume_viewer") == "1"


def _class_selector_redirect(job_id):
    batch_url = reverse("nice_lite:view_batch", args=(job_id,))
    return f"{batch_url}?class_selector=1#batch_class_selector"


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


@login_required(login_url="/login")
@require_GET
def view_batch_dj(request, jobid):
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
def view_batch(request, jobid):
    """Returns batch view."""
    template = "nice_batch/batchview.html"
    batchjob, jobmodel = _get_accessible_batch_job(request, job_id=jobid, log_context="view_batch")
    if batchjob is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    log_by_name = {entry["name"]: entry for entry in batchjob.get_log_tails()}
    stdout_entry = log_by_name.get("stdout.log", {})
    stderr_entry = log_by_name.get("stderr.log", {})
    metadata = jobmodel.master_stats if isinstance(jobmodel.master_stats, dict) else {}
    arguments = _argument_rows(jobmodel)

    context = {
        "jobid"  : jobmodel.id,
        "disp"   : jobmodel.disp,
        "desc"   : jobmodel.desc,
        "proj"   : jobmodel.dset.proj.name,
        "dset"   : jobmodel.dset.name,
        "args"   : jobmodel.args,
        "created": jobmodel.cdat,
        "folder" : batchjob.get_absdir(),
        "jobstats": metadata.get("project_metadata", {}),
        "log"    : _log_parts(stdout_entry["text"]) if stdout_entry.get("exists") else [],
        "error"  : stderr_entry.get("text") if stderr_entry.get("exists") else None,
        "arguments": arguments,
        "submitted_argument_count": sum(argument["submitted"] for argument in arguments),
        "jobstats": metadata.get("project_metadata", {}),
    }
    
    response = render(request, template, context)

    response.set_cookie(key="selected_project_id", value=jobmodel.dset.proj_id)
    response.set_cookie(key="selected_workspace_id", value=jobmodel.dset_id)
    # Ensure Back renders the checksum-gated workspace instead of returning 204.
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
@require_GET
@cache_control(private=True, max_age=300, no_transform=True)
def view_batch_volume_data(request, jobid, volume_name):
    """Stream one owned ab initio 3D MRC output for Mol*."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_volume_data",
        job_id=jobid,
    )
    if (
        batch_job is None
        or jobmodel.status != "finished"
        or jobmodel.prog != "abinitio3D"
    ):
        return HttpResponse(status=404)

    volume = next(
        (
            output
            for output in batch_job.get_volume_outputs()
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
def view_batch_class_selection_export(request, jobid):
    """Download validated, project-ordered 1/0 class-selection state."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "view_batch_class_selection_export",
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
        selection_flags = class_selection_flags(
            selection,
            _deselected_class_ids(request),
        )
    except (ClassSelectionError, OSError, OverflowError, struct.error) as error:
        logger.warning(
            "batch class selection export failed for job %s: %s",
            jobmodel.id,
            error,
        )
        messages.add_message(request, messages.ERROR, f"selection export failed: {error}")
        return redirect(_class_selector_redirect(jobmodel.id))

    response = HttpResponse(
        "".join(f"{state}\n" for state in selection_flags),
        content_type="text/plain; charset=utf-8",
    )
    response["Content-Disposition"] = (
        f'attachment; filename="batch_{jobmodel.id}_class_selection.txt"'
    )
    response["X-Content-Type-Options"] = "nosniff"
    return response


def _save_batch_class_selection_infile(batch_job, selection_flags):
    """Atomically replace the fixed class-selection infile in its source job."""
    job_dir = batch_job.get_safe_job_dir()
    if job_dir is None:
        raise ClassSelectionError("The ab initio 2D job directory is unavailable.")

    infile_path = os.path.abspath(
        os.path.join(job_dir, _BATCH_CLASS_SELECTION_FILENAME)
    )
    resolved_infile = os.path.realpath(infile_path)
    try:
        infile_is_safe = os.path.commonpath((job_dir, resolved_infile)) == job_dir
    except ValueError:
        infile_is_safe = False
    if not infile_is_safe:
        raise ClassSelectionError("The class-selection infile path is unsafe.")

    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=".class_selection.",
            dir=job_dir,
            delete=False,
        ) as temporary_file:
            temporary_path = temporary_file.name
            temporary_file.write(
                "".join(f"{state}\n" for state in selection_flags)
            )
        os.replace(temporary_path, infile_path)
    except OSError:
        if temporary_path is not None:
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
        raise
    return infile_path


@login_required(login_url="/login")
@require_POST
def view_batch_class_selection_run(request, jobid):
    """Save class state beside ab initio 2D output and prefill selection."""
    batch_job, jobmodel = _get_accessible_batch_job(
        request,
        "run_batch_class_selection",
        job_id=jobid,
    )
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
        return redirect("nice_lite:workspace")

    if (
        jobmodel.status != "finished"
        or jobmodel.pckg != "simple"
        or jobmodel.prog != "abinitio2D"
    ):
        messages.add_message(
            request,
            messages.ERROR,
            "2D class selection requires a finished ab initio 2D job",
        )
        return redirect(_class_selector_redirect(jobmodel.id))

    result_project = batch_job.get_result_project_path()
    try:
        selection = load_batch_class_selection(
            result_project,
            jobmodel.dset.proj.dirc,
            jobmodel.id,
        )
        selected_ids = _deselected_class_ids(request)
        selection_flags = class_selection_flags(selection, selected_ids)
        if not selected_ids:
            raise ClassSelectionError(
                "Select at least one class before running a selection job."
            )
        _save_batch_class_selection_infile(batch_job, selection_flags)
    except (ClassSelectionError, OSError, OverflowError, struct.error) as error:
        logger.warning(
            "batch class selection job validation failed for job %s: %s",
            jobmodel.id,
            error,
        )
        messages.add_message(request, messages.ERROR, f"selection job failed: {error}")
        return redirect(_class_selector_redirect(jobmodel.id))

    messages.add_message(
        request,
        messages.SUCCESS,
        "class selection infile saved; review and start the selection job",
    )
    return redirect(reverse(
        "nice_lite:workspace",
        query={
            "selected_job_id": jobmodel.id,
            "class_selection": "1",
            "selected_project_id": jobmodel.dset.proj_id,
            "selected_workspace_id": jobmodel.dset_id,
        },
    ))


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
@cache_control(private=True, max_age=300, no_transform=True)
def view_batch_particle_thumbnail(request, jobid, stack_name, particle_index):
    """Render one owned output-stack image on demand without writing a thumbnail."""
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
