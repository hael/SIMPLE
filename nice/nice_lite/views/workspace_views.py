"""Workspace and workspace-jobs views for stream mode.

This module serves two coupled HTML payloads:
- ``workspace.html``: parent shell containing metadata, controls, and the jobs iframe.
- ``jobs.html``: iframe payload containing stream cards.

It also exposes write endpoints for workspace delete/rename/description updates.
"""

# global imports
import json
import hashlib

# django imports
from django.contrib                 import messages
from django.shortcuts               import redirect, render
from django.contrib.auth.decorators import login_required

# local imports
from ..models                    import JobModel
from ..data_structures.streamjob import StreamJob
from ..data_structures.workspace import Workspace
from ..helpers                   import (
    HttpResponseNoContent,
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


def _normalize_latest_cls2d(jobs):
    """Sort ``latest_cls2D`` by population when classification stats are well-formed."""
    for jobmodel in jobs:
        stats = jobmodel.classification_2D_stats
        if not isinstance(stats, dict):
            continue
        latest = stats.get("latest_cls2D")
        if not isinstance(latest, list):
            continue
        stats["latest_cls2D"] = sorted(latest, key=lambda entry: entry.get("pop", 0) if isinstance(entry, dict) else 0, reverse=True)

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
    job_ids = JobModel.objects.filter(dset=workspace_obj.id).order_by("id").values_list("id", flat=True)
    jobstats = "|".join(StreamJob(id=jobid).get_status() for jobid in job_ids)

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
    }

    # Render only when payload changed to avoid unnecessary parent iframe redraws.
    checksum = hashlib.md5(json.dumps(context, sort_keys=True, default=str).encode()).hexdigest()
    old_checksum = request.COOKIES.get("workspace_checksum", "none")
    if old_checksum == "none" or old_checksum != checksum:
        jobs = JobModel.objects.filter(dset=workspace_obj.id).order_by("id")
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
    workspace_id = get_workspace_id(request)
    project_id = get_project_id(request)
    workspace_obj = Workspace(workspace_id)
    response = HttpResponseNoContent()

    if not _is_workspace_accessible(workspace_obj, project_id, request.user.username):
        return render(request, "jobs.html", {"jobs": []})

    jobs = JobModel.objects.filter(dset=workspace_obj.id).order_by("id")

    # Checksum-gate iframe redraws using current DB state for all jobs in workspace.
    checksum_payload = list(jobs.values())
    checksum = hashlib.md5(json.dumps(checksum_payload, sort_keys=True, default=str).encode()).hexdigest()
    old_checksum = request.COOKIES.get("workspace_jobs_checksum", "none")
    if old_checksum == "none" or old_checksum != checksum:
        _normalize_latest_cls2d(jobs)
        response = render(request, "jobs.html", {"jobs": jobs})
        response.set_cookie(key="workspace_jobs_checksum", value=checksum)

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
    response = redirect("nice_lite:workspace")
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
    response = redirect("nice_lite:workspace")
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
    response = redirect("nice_lite:workspace")
    return response
