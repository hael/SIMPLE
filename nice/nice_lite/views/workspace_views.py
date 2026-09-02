"""Workspace and workspace-jobs views for stream mode.

This module serves two coupled HTML payloads:
- ``workspace.html``: parent shell containing metadata, controls, and the jobs iframe.
- ``jobs.html``: iframe payload containing stream cards.

It also exposes write endpoints for workspace delete/rename/description updates.
The refresh endpoint reconciles externally removed job directories, cards, and
the workspace job counter.
"""

# global imports
import json
import hashlib
import os

# django imports
from django.contrib                 import messages
from django.shortcuts               import redirect, render
from django.contrib.auth.decorators import login_required
from django.views.decorators.http   import require_POST

# local imports
from ..models                    import JobModel
from ..data_structures.batchjob  import BatchJob
from ..data_structures.streamjob import StreamJob
from ..data_structures.workspace import Workspace
from ..features                  import (
    batch_job_controls_enabled,
    workspace_job_refresh_enabled,
)
from ..helpers                   import (
    HttpResponseNoContent,
    clear_checksum_cookies,
    get_integer,
    get_project_id,
    get_string,
    get_workspace_id,
    print_error,
)

# ------------------------------------------------------------------
# Internal Helpers
# ------------------------------------------------------------------

def _is_workspace_accessible(workspace_obj, project_id, username=None):
    """Return True when workspace is valid for project and requester identity.

    Rules:
    - workspace must resolve to a model
    - selected project id must be non-zero and match the workspace project
    - when ``username`` is provided, workspace owner must match it
    """
    if workspace_obj is None or project_id == 0:
        return False
    workspacemodel = workspace_obj.get_workspacemodel()
    if workspacemodel is None or not workspace_obj.in_project(project_id):
        return False

    owner = (workspacemodel.user or "").strip()
    if username is None:
        return True
    if owner == "":
        return False
    return owner == username


def _is_batch_job(jobmodel):
    """Return True for classic jobs stored in the shared JobModel table."""
    return isinstance(jobmodel.master_stats, dict) and jobmodel.master_stats.get("job_type") == "batch"


def _reconcile_local_batch_completions(jobs):
    """Refresh verified local completions before rendering batch controls."""
    changed = False
    for jobmodel in jobs:
        if (
            _is_batch_job(jobmodel)
            and jobmodel.status in ("queued", "running")
            and BatchJob(id=jobmodel.id).reconcile_local_completion()
        ):
            changed = True
    return changed


def _remove_missing_job_records(workspace_obj):
    """Remove records for missing job directories and reset the job counter."""
    workspace_dir = workspace_obj.get_absdir()
    if workspace_dir is None:
        return 0
    workspace_dir = os.path.realpath(workspace_dir)
    if not os.path.isdir(workspace_dir):
        print_error("refresh_workspace_jobs: workspace directory is unavailable")
        return 0

    missing_job_ids = []
    jobs = JobModel.objects.filter(dset=workspace_obj.id)
    for jobmodel in jobs:
        if not isinstance(jobmodel.dirc, str) or not jobmodel.dirc.strip():
            continue
        try:
            job_dir = os.path.realpath(os.path.join(workspace_dir, jobmodel.dirc))
            if os.path.commonpath((workspace_dir, job_dir)) != workspace_dir:
                print_error(f"refresh_workspace_jobs: unsafe directory for job {jobmodel.id}")
                continue
        except (OSError, ValueError):
            print_error(f"refresh_workspace_jobs: invalid directory for job {jobmodel.id}")
            continue
        if not os.path.isdir(job_dir):
            missing_job_ids.append(jobmodel.id)

    if missing_job_ids:
        JobModel.objects.filter(dset=workspace_obj.id, id__in=missing_job_ids).delete()
    workspace_obj.reconcile_job_counter()
    return len(missing_job_ids)


def _filter_by_selection(stats, latest):
    """Reduce ``latest`` to entries whose ``idx`` is in ``stats["selection"]``, when present."""
    selection = stats.get("selection")
    if isinstance(selection, list):
        selected_idxs = set(selection)
        latest = [entry for entry in latest if isinstance(entry, dict) and entry.get("idx") in selected_idxs]
    return latest


def _normalize_latest_cls2d(jobs):
    """Sort ``latest_cls2D`` by population when classification stats are well-formed.

    When ``stats`` carries a ``selection`` array of indices, ``latest_cls2D`` is
    first reduced to only the entries whose ``idx`` is in that selection.
    """
    for jobmodel in jobs:
        # sort latest_cls2D in classification_2D_stats
        stats = jobmodel.classification_2D_stats
        if isinstance(stats, dict):
            latest = stats.get("latest_cls2D")
            if isinstance(latest, list):
                latest = _filter_by_selection(stats, latest)
                stats["latest_cls2D"] = sorted(latest, key=lambda entry: entry.get("res", 0) if isinstance(entry, dict) else 0, reverse=False)
        # sort latest_cls2D in particle_sieving_stats
        stats = jobmodel.particle_sieving_stats
        if isinstance(stats, dict):
            latest = stats.get("latest_cls2D")
            if isinstance(latest, list):
                latest = _filter_by_selection(stats, latest)
                stats["latest_cls2D"] = sorted(latest, key=lambda entry: entry.get("res", 0) if isinstance(entry, dict) else 0, reverse=False)

# ------------------------------------------------------------------
# Views
# ------------------------------------------------------------------

@login_required(login_url="/login")
def view_workspace(request):
    """Render parent workspace payload, or return 204 when unchanged."""
    workspace_id = get_workspace_id(request)
    project_id = get_project_id(request)
    workspace_obj = Workspace(workspace_id)
    response = HttpResponseNoContent()

    if not _is_workspace_accessible(workspace_obj, project_id, request.user.username):
        print_error("workspace is not accessible for selected project/user")
        messages.add_message(request, messages.ERROR, "workspace is not accessible")
        return response

    # Include stream statuses in checksum seed so parent iframe updates on state changes.
    jobs = list(JobModel.objects.filter(dset=workspace_obj.id).order_by("id"))
    if batch_job_controls_enabled() and _reconcile_local_batch_completions(jobs):
        jobs = list(JobModel.objects.filter(dset=workspace_obj.id).order_by("id"))
    jobstats = "|".join(
        job.status if _is_batch_job(job) else StreamJob(id=job.id).get_status()
        for job in jobs
    )

    workspacemodel = workspace_obj.get_workspacemodel()
    projectmodel = workspacemodel.proj
    context = {
        "current_project_id": project_id,
        "current_workspace_id": workspace_obj.get_id(),
        "current_project_name": projectmodel.name,
        "current_workspace_name": workspacemodel.name,
        "created": workspacemodel.cdat,
        "modified": workspacemodel.mdat,
        "user": workspacemodel.user,
        "folder": workspace_obj.get_linkpath(),
        "description": workspacemodel.desc,
        "jobstats": jobstats,
        "workspace_job_refresh_enabled": workspace_job_refresh_enabled(),
    }

    # Render only when payload changed to avoid unnecessary parent iframe redraws.
    checksum = hashlib.md5(json.dumps(context, sort_keys=True, default=str).encode()).hexdigest()
    old_checksum = request.COOKIES.get("workspace_checksum", "none")
    if old_checksum == "none" or old_checksum != checksum:
        _normalize_latest_cls2d(jobs)
        context["jobs"] = jobs
        response = render(request, "workspace.html", context)
        response.set_cookie(key="workspace_checksum", value=checksum)
        response.delete_cookie(key="workspace_jobs_checksum")

    response.set_cookie(key="selected_project_id", value=project_id)
    response.set_cookie(key="selected_workspace_id", value=workspace_obj.id)
    return response


@login_required(login_url="/login")
def view_workspace_jobs(request):
    """Render jobs iframe payload, or return 204 when unchanged."""
    workspace_id  = get_workspace_id(request)
    project_id    = get_project_id(request)
    workspace_obj = Workspace(workspace_id)
    response      = HttpResponseNoContent()

    if not _is_workspace_accessible(workspace_obj, project_id, request.user.username):
        return render(request, "jobs.html", {"jobs": []})

    jobs = JobModel.objects.filter(dset=workspace_obj.id).order_by("id")
    controls_enabled = batch_job_controls_enabled()
    if controls_enabled and _reconcile_local_batch_completions(jobs):
        jobs = JobModel.objects.filter(dset=workspace_obj.id).order_by("id")

    # Checksum-gate iframe redraws using current DB state for all jobs in workspace.
    checksum_payload = {
        "jobs": list(jobs.values()),
        "batch_job_controls_enabled": controls_enabled,
    }
    checksum = hashlib.md5(json.dumps(checksum_payload, sort_keys=True, default=str).encode()).hexdigest()
    old_checksum = request.COOKIES.get("workspace_jobs_checksum", "none")
    if old_checksum == "none" or old_checksum != checksum:
        _normalize_latest_cls2d(jobs)
        response = render(request, "jobs.html", {
            "jobs": jobs,
            "batch_job_controls_enabled": controls_enabled,
        })
        response.set_cookie(key="workspace_jobs_checksum", value=checksum)

    return response


@login_required(login_url="/login")
@require_POST
def view_refresh_workspace_jobs(request):
    """Reconcile externally deleted job directories and the workspace job counter."""
    response = redirect("nice_lite:workspace")
    if not workspace_job_refresh_enabled():
        messages.add_message(request, messages.ERROR, "workspace job refresh is disabled")
        return response

    workspace_id = get_workspace_id(request)
    project_id = get_project_id(request)
    workspace_obj = Workspace(workspace_id)
    if not _is_workspace_accessible(workspace_obj, project_id, request.user.username):
        messages.add_message(request, messages.ERROR, "invalid workspace selection")
        return response

    removed_count = _remove_missing_job_records(workspace_obj)
    if removed_count == 0:
        messages.add_message(request, messages.INFO, "job cards are up to date")
    else:
        messages.add_message(
            request,
            messages.INFO,
            f"removed {removed_count} job card{'s' if removed_count != 1 else ''} with missing directories",
        )
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
def view_delete_workspace(request):
    """Delete a workspace and redirect back to workspace landing page."""
    deleteworkspaceid = get_integer(request.POST, "delete_workspace_id")
    workspace = Workspace(deleteworkspaceid)
    workspacemodel = workspace.get_workspacemodel()
    # Derive project from target workspace, not selected-project cookie.
    project_id = workspacemodel.proj.id if workspacemodel is not None else 0

    if not _is_workspace_accessible(workspace, project_id, request.user.username):
        messages.add_message(request, messages.ERROR, "invalid workspace selection")
    elif workspace.delete():
        messages.add_message(request, messages.INFO, "workspace deleted successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to delete workspace")
    response = redirect("nice_lite:index")
    return response


@login_required(login_url="/login")
def view_update_workspace_name(request):
    """Rename selected workspace and redirect back to workspace landing page."""
    workspaceid = get_workspace_id(request)
    workspace = Workspace(workspaceid)
    workspacemodel = workspace.get_workspacemodel()
    # Derive project from target workspace, not selected-project cookie.
    project_id = workspacemodel.proj.id if workspacemodel is not None else 0
    newworkspacename = get_string(request.POST, "new_workspace_name")

    if not _is_workspace_accessible(workspace, project_id, request.user.username):
        messages.add_message(request, messages.ERROR, "invalid workspace selection")
    elif workspace.rename(newworkspacename):
        messages.add_message(request, messages.INFO, "workspace renamed successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to rename workspace")
    response = redirect("nice_lite:index")
    return response


@login_required(login_url="/login")
def view_update_workspace_description(request):
    """Update selected workspace description and redirect back to workspace page."""
    workspaceid = get_workspace_id(request)
    workspace = Workspace(workspaceid)
    workspacemodel = workspace.get_workspacemodel()
    # Derive project from target workspace, not selected-project cookie.
    project_id = workspacemodel.proj.id if workspacemodel is not None else 0
    newworkspacedescription = get_string(request.POST, "new_workspace_description")

    if not _is_workspace_accessible(workspace, project_id, request.user.username):
        messages.add_message(request, messages.ERROR, "invalid workspace selection")
    elif workspace.updateDescription(newworkspacedescription):
        messages.add_message(request, messages.INFO, "description updated successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to update description")
    response = redirect("nice_lite:index")
    return response
