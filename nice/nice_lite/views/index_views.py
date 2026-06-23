"""Workspace landing view for stream mode."""

from django.contrib                 import messages
from django.contrib.auth            import authenticate, login
from django.contrib.auth            import logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms      import AuthenticationForm
from django.shortcuts               import redirect, render
from django.urls                    import reverse

from ..data_structures.project   import Project
from ..data_structures.workspace import Workspace
from ..helpers                   import clear_checksum_cookies, get_project_id, get_workspace_id
from ..models                    import ProjectModel, WorkspaceModel

def view_user_login(request):
    """Authenticate user and render login page on failure."""
    template = "login.html"
    if request.method == "POST":
        username = request.POST.get("username", "")
        password = request.POST.get("password", "")
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect("nice_lite:index")
        messages.add_message(request, messages.ERROR, "Invalid credentials")
    form     = AuthenticationForm()
    response = render(request, template, {'form':form, 'title':'log in'})
    return response

def view_user_logout(request):
    """Log out the authenticated user and redirect to login."""
    logout(request)
    return redirect("nice_lite:login")

@login_required(login_url="/login")
def view_index(request):
    """Render the top-level stream workspace page."""
    template = "workspace.html"

    # Resolve current UI selection state from cookies/params.
    projectid = get_project_id(request)
    workspaceid = get_workspace_id(request)
    username = request.user.username
    iframeurl = None

    # Populate top-level selectors for projects and workspaces.
    projects = ProjectModel.objects.all()
    workspaces = []

    if username is None:
        messages.add_message(request, messages.ERROR, "username is none")

    if projectid is None:
        if len(projects) == 0:
            messages.add_message(request, messages.INFO, "please create a project")
        else:
            messages.add_message(request, messages.INFO, "please select a project")
    else:
        workspaces = WorkspaceModel.objects.filter(proj=projectid)

    if workspaceid is None:
        if projectid is not None:
            if len(workspaces) == 0:
                messages.add_message(request, messages.INFO, "please create a workspace")
            else:
                messages.add_message(request, messages.INFO, "please select a workspace")
    else:
        workspace = Workspace(workspaceid)
        # Guard against stale workspace ids that do not belong to the selected project.
        if not workspace.in_project(projectid):
            workspaceid = None
            messages.add_message(request, messages.INFO, "please select a workspace")

    # Sentinel ids from the UI drive special navigation/creation paths.
    if projectid == -1:
        iframeurl = reverse("nice_lite:new_project", args=["stream"])
    elif workspaceid == -1 and username is not None:
        messages.add_message(request, messages.INFO, "created new workspace")
        project = Project(projectid)
        new_workspace = Workspace()
        new_workspace.new(project, username)
        workspaceid = new_workspace.get_id()
        workspaces = WorkspaceModel.objects.filter(proj=projectid)
        iframeurl = reverse("nice_lite:workspace", query={"selected_workspace_id": workspaceid})
    elif workspaceid is not None and workspaceid > 0:
        iframeurl = reverse("nice_lite:workspace", query={"selected_workspace_id": workspaceid})

    context = {
        "current_project_id"   : projectid,
        "current_workspace_id" : workspaceid,
        "projects"             : projects,
        "workspaces"           : workspaces,
        "iframeurl"            : iframeurl,
    }

    # Persist current selection to keep the stream UI stateful across requests.
    response = render(request, template, context)
    response.set_cookie(key="selected_project_id", value=projectid)
    response.set_cookie(key="selected_workspace_id", value=workspaceid)
    response.set_cookie(key="mode", value="stream")
    clear_checksum_cookies(request, response)
    return response


