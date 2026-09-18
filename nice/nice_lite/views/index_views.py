"""Top-level index and auth views for stream mode.

This module serves:
- ``login.html`` authentication entry/exit endpoints
- ``index.html`` shell selection logic for project/workspace navigation
"""

# django imports
from django.contrib                 import messages
from django.contrib.auth            import login, logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms      import AuthenticationForm
from django.shortcuts               import redirect, render
from django.urls                    import reverse
from django.views.decorators.http   import require_POST

# local imports
from ..helpers                   import clear_checksum_cookies, get_project_id, get_workspace_id
from ..models                    import ProjectModel, WorkspaceModel


# ------------------------------------------------------------------
# Authentication Views
# ------------------------------------------------------------------

def view_user_login(request):
    """Authenticate user and render login page on failure."""
    template = "login.html"
    form = AuthenticationForm()

    if request.method == "POST":
        form = AuthenticationForm(request, data=request.POST)
        if form.is_valid():
            login(request, form.get_user())
            return redirect("nice_lite:index")
        messages.add_message(request, messages.ERROR, "Invalid credentials")
    
    response = render(request, template, {"form": form, "title": "log in"})
    return response


@login_required(login_url="/login")
@require_POST
def view_user_logout(request):
    """Log out the authenticated user and redirect to login."""
    logout(request)
    return redirect("nice_lite:login")


# ------------------------------------------------------------------
# Index View
# ------------------------------------------------------------------

@login_required(login_url="/login")
def view_index(request):
    """Render the top-level stream workspace page."""
    template = "index.html"

    # Resolve current UI selection state from cookies/params.
    projectid = get_project_id(request)
    project_only_request = (
        request.method == "GET"
        and "selected_project_id" in request.GET
        and "selected_workspace_id" not in request.GET
    )
    workspaceid = None if project_only_request else get_workspace_id(request)
    username = request.user.username
    iframeurl = None
    projectmodel = None

    # Show only projects that have at least one workspace owned by this user.
    projects = ProjectModel.objects.filter(workspacemodel__user=username).distinct()
    project_ids = set(projects.values_list("id", flat=True))
    workspaces = []

    if username is None:
        messages.add_message(request, messages.ERROR, "username is none")

    if projectid is not None and projectid != -1 and projectid not in project_ids:
        projectid = None
        workspaceid = None
        messages.add_message(request, messages.INFO, "please select a project")

    if projectid is None:
        if len(projects) == 0:
            messages.add_message(request, messages.INFO, "please create a project")
        else:
            messages.add_message(request, messages.INFO, "please select a project")
    elif projectid > 0:
        projectmodel = projects.filter(id=projectid).first()
        workspaces = WorkspaceModel.objects.filter(proj=projectid, user=username)

    if workspaceid is None:
        if projectid is not None:
            if len(workspaces) == 0:
                messages.add_message(request, messages.INFO, "please create a workspace")
            else:
                messages.add_message(request, messages.INFO, "please select a workspace")
    else:
        workspace_model = WorkspaceModel.objects.filter(id=workspaceid, proj=projectid, user=username).first()
        # Guard against stale/foreign workspace ids, including the legacy GET
        # creation sentinel. Workspace creation is a POST-only action.
        if workspace_model is None:
            workspaceid = None
            messages.add_message(request, messages.INFO, "please select a workspace")

    # Sentinel ids from the UI drive special navigation/creation paths.
    if projectid == -1:
        # Project sentinel routes user to the new-project page.
        iframeurl = reverse("nice_lite:new_project", args=["stream"])
    elif workspaceid is not None and workspaceid > 0:
        iframeurl = reverse("nice_lite:workspace", query={"selected_project_id": projectid, "selected_workspace_id": workspaceid})

    context = {
        "current_project_id": projectid,
        "current_workspace_id": workspaceid,
        "project": projectmodel,
        "projects": projects,
        "workspaces": workspaces,
        "iframeurl": iframeurl,
    }

    # Persist current selection to keep the stream UI stateful across requests.
    # Opening the new-project form is only a temporary navigation state, so keep
    # the last valid project/workspace cookies for the Back action to restore.
    response = render(request, template, context)
    if projectid != -1:
        response.set_cookie(key="selected_project_id", value=projectid)
        response.set_cookie(key="selected_workspace_id", value=workspaceid)
    response.set_cookie(key="mode", value="stream")
    clear_checksum_cookies(request, response)
    return response
