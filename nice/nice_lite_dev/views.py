"""
views.py — shared views used by both stream and classic modes.

Covers:
  - Authentication : login, logout
  - Entry point    : index (mode-cookie dispatch)
  - Project setup  : new_project, create_project
  - Utilities      : file_browser
"""

# global imports
from django.shortcuts               import redirect, render
from django.contrib                 import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth            import login, logout, authenticate

# local imports
from .helpers                    import *
from .app_views.coreview         import LoginView
from .app_views.newprojectview   import NewProjectView
from .app_views.filebrowserview  import FileBrowserView
from .data_structures.project    import Project
from .data_structures.dataset    import Dataset


# ------------------------------------------------------------------
# Authentication
# ------------------------------------------------------------------

def user_login(request):
    """
    Authenticate the user on POST. Redirects to stream on success;
    re-renders the login page with an error message on failure.
    """
    if request.method == 'POST':
        username = request.POST.get('username', '')
        password = request.POST.get('password', '')
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('nice_lite:stream')
        messages.add_message(request, messages.ERROR, "Invalid credentials")
    loginview = LoginView(request)
    return loginview.render()


def user_logout(request):
    """Log out the authenticated user and redirect to the login page."""
    logout(request)
    return redirect('nice_lite:login')


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

@login_required(login_url="/login")
def index(request):
    """
    Dispatch to stream or classic based on the mode cookie.
    Falls back to the mode-selection page when the cookie is absent.
    """
    mode = request.COOKIES.get('mode', 'none')
    if mode == "classic":
        return redirect('nice_lite:classic')
    elif mode == "stream":
        return redirect('nice_lite:stream')
    return render(request, "index.html", {})


# ------------------------------------------------------------------
# Project setup
# ------------------------------------------------------------------

@login_required(login_url="/login")
def new_project(request, mode):
    """Render the create-new-project page."""
    newprojectview = NewProjectView(request, mode)
    return newprojectview.render()


@login_required(login_url="/login")
def create_project(request):
    """
    Create a new project with one empty dataset, set session cookies,
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
    dataset = Dataset()
    if not dataset.new(project, username):
        print_error("failed to create new dataset")
        messages.add_message(request, messages.ERROR, "Failed to create new dataset")
        return response
    response.set_cookie(key='selected_project_id', value=project.id)
    response.set_cookie(key='selected_dataset_id', value=dataset.id)
    clear_checksum_cookies(request, response)
    return response


# ------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------

@login_required(login_url="/login")
def file_browser(request, browser_type, path=None):
    """
    Render the file browser page.

    Path resolution order:
      1. URL path segment (passed directly as an argument)
      2. 'selectedpath' GET parameter
      3. Root directory of the currently selected project (from cookie)
      4. Filesystem root '/'
    """
    if path is None:
        if "selectedpath" in request.GET:
            path = request.GET['selectedpath']
        elif request.COOKIES.get('selected_project_id', 'none') != 'none':
            project = Project(id=int(request.COOKIES['selected_project_id']))
            path = project.absdir or "/"
        else:
            path = "/"
    filebrowserview = FileBrowserView(request, browser_type, path)
    return filebrowserview.render()
