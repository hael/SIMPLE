"""Workspace landing view for stream mode."""

import hashlib
import json

from django.contrib                 import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts               import redirect, render

from ..data_structures.streamjob import StreamJob
from ..data_structures.workspace import Workspace
from ..helpers                   import HttpResponseNoContent, get_integer, get_project_id, get_string, get_workspace_id, print_error
from ..models                    import JobModel

@login_required(login_url="/login/")
def view_workspace(request):
    """Return workspace panel content for the selected workspace id."""
    workspace_id  = get_workspace_id(request)
    project_id    = get_project_id(request)
    workspace_obj = Workspace(workspace_id)
    response      = HttpResponseNoContent()

    if workspace_obj is None:
        print_error("workspace is None")
        messages.add_message(request, messages.ERROR, "workspace is None")
        return response
    if project_id == 0:
        print_error("projectid is zero")
        messages.add_message(request, messages.ERROR, "projectid is zero")
        return response
    if not workspace_obj.in_project(project_id):
        print_error("workspace is not in project")
        messages.add_message(request, messages.ERROR, "workspace is not in project")
        return response

    # Fold process statuses into a checksum seed so the iframe refreshes on state change.
    jobs = JobModel.objects.filter(dset=workspace_obj.id)
    jobstats = ""
    for jobmodel in jobs:
        streamjob = StreamJob(id=jobmodel.id)
        jobstats += streamjob.get_status()

    workspacemodel = workspace_obj.get_workspacemodel()
    projectmodel = workspacemodel.proj
    context = {
        "current_project_id"     : project_id,
        "current_workspace_id"   : workspace_obj.get_id(),
        "current_project_name"   : projectmodel.name,
        "current_workspace_name" : workspacemodel.name,
        "created"                : workspacemodel.cdat,
        "modified"               : workspacemodel.mdat,
        "user"                   : workspacemodel.user,
        "folder"                 : workspace_obj.get_linkpath(),
        "description"            : workspacemodel.desc,
        "jobstats"               : jobstats,
    }

    # Render only when the payload changed; otherwise return 204 to avoid redundant updates.
    checksum = hashlib.md5(json.dumps(context, sort_keys=True, default=str).encode()).hexdigest()
    old_checksum = request.COOKIES.get("workspace_checksum", "none")
    if old_checksum == "none" or old_checksum != checksum:
        jobs = JobModel.objects.filter(dset=workspace_obj.id)
        for jobmodel in jobs:
            if "latest_cls2D" in jobmodel.classification_2D_stats:
                jobmodel.classification_2D_stats["latest_cls2D"] = sorted(
                    jobmodel.classification_2D_stats["latest_cls2D"],
                    key=lambda d: d["pop"],
                    reverse=True,
                )
        context["jobs"] = jobs
        response = render(request, "nice_stream/workspace.html", context)
        response.set_cookie(key="workspace_checksum", value=checksum)

    response.set_cookie(key="selected_project_id", value=project_id)
    response.set_cookie(key="selected_workspace_id", value=workspace_obj.id)
    return response


@login_required(login_url="/login/")
def view_delete_workspace(request):
    """Delete a workspace and return to stream landing page."""
    deleteworkspaceid = get_integer(request.POST, "delete_workspace_id")
    workspace = Workspace(deleteworkspaceid)
    if workspace.delete():
        messages.add_message(request, messages.INFO, "workspace deleted successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to delete workspace")
    response = redirect("nice_lite:stream")
    return response


@login_required(login_url="/login/")
def view_update_workspace_name(request):
    """Rename selected workspace and return to stream landing page."""
    workspaceid = get_workspace_id(request)
    workspace = Workspace(workspaceid)
    newworkspacename = get_string(request.POST, "new_workspace_name")
    if workspace.rename(newworkspacename):
        messages.add_message(request, messages.INFO, "workspace renamed successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to rename workspace")
    response = redirect("nice_lite:stream")
    return response


@login_required(login_url="/login/")
def view_update_workspace_description(request):
    """Update selected workspace description and return to stream landing page."""
    workspaceid = get_workspace_id(request)
    workspace = Workspace(workspaceid)
    newworkspacedescription = get_string(request.POST, "new_workspace_description")
    if workspace.updateDescription(newworkspacedescription):
        messages.add_message(request, messages.INFO, "description updated successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to update description")
    response = redirect("nice_lite:stream")
    return response


# Backward-compatible alias for existing routes/imports.
workspace = view_workspace
delete_workspace = view_delete_workspace
update_workspace_name = view_update_workspace_name
update_workspace_description = view_update_workspace_description