"""JSON and media API endpoints for NICE Lite.

This module serves:
- worker heartbeat ingestion endpoints for stream and classic jobs
- authenticated project-scoped image reads for UI payloads

Heartbeat endpoints are intentionally CSRF-exempt because they are called by
non-browser worker processes.
"""

# global imports
import hmac
import json
import os

# django imports
from django.http                    import HttpResponse
from django.http                    import JsonResponse
from django.views.decorators.http   import require_GET
from django.views.decorators.http   import require_POST
from django.views.decorators.cache  import cache_control
from django.views.decorators.csrf   import csrf_exempt
from django.contrib.auth.decorators import login_required

# local imports
from .helpers                       import dict_present
from .helpers                       import get_project_id
from .helpers                       import print_error
from .data_structures.batchjob      import BatchJob
from .data_structures.project       import Project
from .data_structures.streamjob     import StreamJob
from .data_structures.workspace     import Workspace
from .models                        import ProjectModel


# ------------------------------------------------------------------
# Internal Helpers
# ------------------------------------------------------------------

def _parse_request_json(request):
    """Return decoded JSON object from request body, or ``None`` on failure."""
    try:
        return json.loads(request.body.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return None


def _is_known_job(jobmodel):
    """Return True when the resolved job model exists in the database."""
    return jobmodel is not None


def _get_valid_job_id(payload):
    """Return strictly valid job id integer from payload, otherwise ``None``."""
    if "jobid" not in payload:
        return None

    jobid = payload["jobid"]
    # Reject booleans explicitly since bool is a subclass of int in Python.
    if isinstance(jobid, bool) or not isinstance(jobid, int):
        return None
    if jobid <= 0:
        return None
    return jobid


def _has_numeric_version(payload):
    """Return True when payload has numeric version (int/float, not bool)."""
    if "version" not in payload:
        return False
    version = payload["version"]
    if isinstance(version, bool):
        return False
    return isinstance(version, (int, float))


def _is_worker_authorized(request):
    """Return True when worker request provides a valid shared secret.

    Auth is opt-in: when ``NICE_LITE_WORKER_TOKEN`` is unset, checks are bypassed.
    Supported sources:
    - ``X-Worker-Token`` header
    - ``Authorization: Bearer <token>`` header
    """
    expected_token = (os.environ.get("NICE_LITE_WORKER_TOKEN") or "").strip()
    if expected_token == "":
        return True

    provided_token = (request.headers.get("X-Worker-Token") or "").strip()
    if provided_token == "":
        authorization = (request.headers.get("Authorization") or "").strip()
        if authorization.lower().startswith("bearer "):
            provided_token = authorization[7:].strip()

    if provided_token == "":
        return False
    return hmac.compare_digest(provided_token, expected_token)


def _resolve_safe_image_path(request, src):
    """Return project-scoped absolute file path, or ``None`` when invalid."""
    selected_project_id = get_project_id(request)
    if selected_project_id is None:
        return None

    username = request.user.username
    accessible_project_ids = set(
        ProjectModel.objects.filter(workspacemodel__user=username).distinct().values_list("id", flat=True)
    )
    if selected_project_id not in accessible_project_ids:
        return None

    project = Project(id=selected_project_id)
    if not project.absdir:
        return None

    base_dir = os.path.realpath(project.absdir)
    normalized = (src or "").strip()
    if normalized == "":
        return None

    if not os.path.isabs(normalized):
        normalized = os.path.join(base_dir, normalized)
    normalized = os.path.realpath(normalized)

    try:
        # Keep image access constrained to the selected project root.
        if os.path.commonpath([normalized, base_dir]) != base_dir:
            return None
    except ValueError:
        return None
    return normalized


def _is_allowed_image_path(path):
    """Return True when path suffix is one of supported image extensions."""
    allowed_extensions = {".jpg", ".jpeg", ".png", ".gif", ".webp"}
    extension = os.path.splitext(path)[1].lower()
    return extension in allowed_extensions


def _image_content_type_from_path(path):
    """Return image MIME type inferred from filename extension."""
    extension = os.path.splitext(path)[1].lower()
    extension_to_content_type = {
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".png": "image/png",
        ".gif": "image/gif",
        ".webp": "image/webp",
    }
    return extension_to_content_type.get(extension, "application/octet-stream")


# ------------------------------------------------------------------
# API Endpoints
# ------------------------------------------------------------------

@csrf_exempt
@require_POST
def index(request):
    """Handle stream heartbeat updates from worker clients.

    Expected payload keys include:
    - ``version`` (number)
    - ``jobid`` (positive integer)
    - ``stream_heartbeat`` (dict)
    """
    response = {}
    if not _is_worker_authorized(request):
        print_error("unauthorized worker request")
        return JsonResponse(response, status=403)

    request_json = _parse_request_json(request)
    if request_json is None:
        print_error("failed to load json from request")
        return JsonResponse(response, status=400)
    if not isinstance(request_json, dict):
        print_error("request payload is not a JSON object")
        return JsonResponse(response, status=400)

    if not _has_numeric_version(request_json):
        print_error("version missing from request")
        return JsonResponse(response, status=400)
    jobid = _get_valid_job_id(request_json)
    if jobid is None:
        print_error("jobid missing from request")
        return JsonResponse(response, status=400)

    streamjob = StreamJob(id=jobid)
    jobmodel = streamjob.get_jobmodel()
    if not _is_known_job(jobmodel):
        print_error(f"unknown stream job {jobid}")
        return JsonResponse(response, status=404)

    if dict_present(request_json, "stream_heartbeat"):
        if streamjob.update_stats(request_json):
            response = streamjob.get_master_update()
            return JsonResponse(response)
        return JsonResponse(response, status=400)

    if dict_present(request_json, "heartbeat"):
        print_error("classic heartbeat is not supported on this endpoint")
        return JsonResponse(response, status=400)

    print_error("unknown job heartbeat type")
    return JsonResponse(response, status=400)


@csrf_exempt
@require_POST
def index_classic(request):
    """Handle classic heartbeat updates from worker clients.

    Expected payload keys include:
    - ``jobid`` (positive integer)
    - classic stats/update fields consumed by ``BatchJob.updateStats``
    """
    response = {}
    if not _is_worker_authorized(request):
        print_error("unauthorized worker request")
        return JsonResponse(response, status=403)

    request_json = _parse_request_json(request)
    if request_json is None:
        print_error("failed to load json from classic request")
        return JsonResponse(response, status=400)
    if not isinstance(request_json, dict):
        print_error("classic payload is not a JSON object")
        return JsonResponse(response, status=400)

    jobid = _get_valid_job_id(request_json)
    if jobid is None:
        print_error("jobid missing from classic request")
        return JsonResponse(response, status=400)

    job = BatchJob(id=jobid)
    if not _is_known_job(job.jobmodel):
        print_error(f"unknown classic job {jobid}")
        return JsonResponse(response, status=404)

    workspace = Workspace(job.jobmodel.wspc_id)
    project = Project(job.jobmodel.wspc.proj_id)
    response = job.updateStats(request_json, project, workspace)
    return JsonResponse(response)


@login_required(login_url="/login")
@require_GET
@cache_control(max_age=300, must_revalidate=True, no_transform=True)
def image(request, src):
    """Serve supported image content constrained to the selected project root."""
    safe_path = _resolve_safe_image_path(request, src)
    if safe_path is None:
        print_error("invalid image path request")
        return HttpResponse(status=404)
    if not _is_allowed_image_path(safe_path):
        print_error("unsupported image extension")
        return HttpResponse(status=404)

    try:
        content_type = _image_content_type_from_path(safe_path)
        with open(safe_path, "rb") as image_file:
            return HttpResponse(image_file.read(), content_type=content_type)
    except OSError:
        print_error("image IO error")
        return HttpResponse(status=404)

