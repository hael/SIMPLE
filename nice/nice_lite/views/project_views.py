"""Project creation and project-setup page views.

This module serves:
- ``newproject.html`` rendering for project creation UI
- write endpoint that creates a project and an initial workspace
"""

# global imports
import os
import shutil

# django imports
from django.contrib                 import messages
from django.shortcuts               import redirect, render
from django.views.decorators.http   import require_POST
from django.contrib.auth.decorators import login_required

# local imports
from ..data_structures.project   import Project
from ..data_structures.workspace import Workspace
from ..helpers                   import clear_checksum_cookies, get_string, print_error


# ------------------------------------------------------------------
# Views
# ------------------------------------------------------------------


@login_required(login_url="/login")
@require_POST
def view_create_project(request):
    """
    Create a new project with one empty workspace, set session cookies,
    and redirect to the stream shell.
    """
    response = redirect("nice_lite:index")
    username = request.user.username
    projname = get_string(request.POST, "new_project_name")
    projdirc = get_string(request.POST, "new_project_dirc")
    if projname is None:
        messages.add_message(request, messages.ERROR, "Project name is missing")
        return response
    if projdirc is None:
        messages.add_message(request, messages.ERROR, "Project directory is missing")
        return response

    project = Project()
    project_path = os.path.join(projdirc, str(projname).replace(" ", "_"))
    project_path_existed = os.path.isdir(project_path)
    if not project.new(projname, projdirc):
        print_error("failed to create new project")
        messages.add_message(request, messages.ERROR, "Failed to create new project")
        return response

    workspace = Workspace()
    if not workspace.new(project, username):
        print_error("failed to create new workspace")
        # Roll back project DB entry and freshly created project directory.
        projectmodel = project.get_projectmodel()
        if projectmodel is not None:
            projectmodel.delete()
        if not project_path_existed and os.path.isdir(project_path):
            try:
                shutil.rmtree(project_path)
            except OSError:
                print_error(f"failed to remove rolled-back project directory {project_path}")
        messages.add_message(request, messages.ERROR, "Failed to create new workspace")
        return response

    response.set_cookie(key="selected_project_id", value=project.id)
    response.set_cookie(key="selected_workspace_id", value=workspace.id)
    clear_checksum_cookies(request, response)
    return response


@login_required(login_url="/login")
def view_new_project(request, mode):
    """Render the create-new-project page."""
    template = "newproject.html"
    context = {"mode": mode}
    response = render(request, template, context)
    # New-project page should always clear stale checksums from prior views.
    clear_checksum_cookies(request, response)
    return response
