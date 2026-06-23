from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect, render

from ..data_structures.project import Project
from ..data_structures.workspace import Workspace
from ..helpers import clear_checksum_cookies, get_string, print_error


@login_required(login_url="/login")
def view_create_project(request):
    """
    Create a new project with one empty workspace, set session cookies,
    and redirect to the stream view.
    """
    response = redirect('nice_lite:stream')
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
    if not project.new(projname, projdirc):
        print_error("failed to create new project")
        messages.add_message(request, messages.ERROR, "Failed to create new project")
        return response

    workspace = Workspace()
    if not workspace.new(project, username):
        print_error("failed to create new workspace")
        messages.add_message(request, messages.ERROR, "Failed to create new workspace")
        return response

    response.set_cookie(key='selected_project_id', value=project.id)
    response.set_cookie(key='selected_workspace_id', value=workspace.id)
    clear_checksum_cookies(request, response)
    return response


# Backward-compatible alias for existing imports/routes.
create_project = view_create_project

@login_required(login_url="/login")
def view_new_project(request, mode):
    """Render the create-new-project page."""
    template = "newproject.html"
    context = {
        "mode" : mode
    }
    response = render(request, template, context)
    clear_checksum_cookies(request,response)
    return response

