"""File browser view.

This module renders ``filebrowser.html`` for authenticated users. Normal file
and directory browsing is scoped to the selected project root; explicit external
directory pickers can browse before or outside a project.
"""

# global imports
import os

# django imports
from django.contrib.auth.decorators import login_required
from django.shortcuts               import render

# local imports
from ..data_structures.project   import Project
from ..helpers                   import get_project_id
from ..models                    import ProjectModel


# ------------------------------------------------------------------
# Views
# ------------------------------------------------------------------


@login_required(login_url="/login")
def view_file_browser(request, type, path=None):
    """
    Render the file browser page.

    Path resolution order (project-scoped):
      1. URL path segment (passed directly as an argument)
      2. 'selectedpath' GET parameter
      3. Root directory of the currently selected project

    Browsing is restricted to the selected project root for the current user,
    except for explicit external directory pickers.
    """
    template = "filebrowser.html"
    known_file_extensions = {
        ".spi", ".mrc", ".mrcs", ".dm3", ".bin", ".gain", ".eer", ".tiff", ".tif", ".jpeg",
        ".jpg", ".txt", ".simple", ".dat", ".img", ".map", ".head", ".ctf", ".raw", ".sbin",
        ".dbin", ".asc", ".box", ".pdb", ".star", ".hdf", ".pdf", ".ps",
    }
    path_isdir = None
    parentdir = None
    error = False
    errortext = ""
    files = []
    dirs = []

    username = request.user.username
    purpose = request.GET.get("purpose", "")
    is_unrestricted_picker = purpose == "external_input" or (type == "dir" and purpose == "project_root")
    selected_path = request.GET.get("selectedpath")
    path_was_requested = path is not None or selected_path is not None
    using_remembered_path = request.GET.get("remembered") == "1"
    selected_project_id = get_project_id(request)
    accessible_project_ids = set(
        ProjectModel.objects.filter(workspacemodel__user=username).distinct().values_list("id", flat=True)
    )

    if is_unrestricted_picker:
        if path is None:
            path = selected_path if selected_path is not None else os.path.expanduser("~")

        path = (path or "").strip()
        if path == "":
            path = os.path.expanduser("~")
        path = os.path.realpath(path)
        if using_remembered_path and not os.path.isdir(path):
            path = os.path.realpath(os.path.expanduser("~"))
    elif selected_project_id is None or selected_project_id not in accessible_project_ids:
        error = True
        errortext = "invalid project selection"
        path = ""
    else:
        project = Project(id=selected_project_id)
        if not project.absdir:
            error = True
            errortext = "project root is unavailable"
            path = ""
        else:
            base_dir = os.path.realpath(project.absdir)

            if path is None:
                path = selected_path if selected_path is not None else base_dir

            path = (path or "").strip()
            if path == "":
                path = base_dir
            elif not os.path.isabs(path):
                path = os.path.join(base_dir, path)

            path = os.path.realpath(path)
            # Keep browsing constrained to the selected project root.
            if os.path.commonpath([path, base_dir]) != base_dir:
                if using_remembered_path:
                    path = base_dir
                else:
                    error = True
                    errortext = "path outside project"
                    path = base_dir
            elif using_remembered_path and not os.path.isdir(path):
                path = base_dir

    try:
        if not error:
            if not os.path.exists(path):
                error = True
                errortext = "path does not exist"
            if os.path.isdir(path):
                path_isdir = True
                parentdir = os.path.dirname(path)
                if (
                    not is_unrestricted_picker
                    and selected_project_id is not None
                    and selected_project_id in accessible_project_ids
                    and not error
                ):
                    project = Project(id=selected_project_id)
                    base_dir = os.path.realpath(project.absdir) if project.absdir else ""
                    if base_dir and os.path.commonpath([parentdir, base_dir]) != base_dir:
                        parentdir = base_dir
            else:
                path_isdir = False
    except OSError:
        if not error:
            error = True
            errortext = "error"

    if not error and path_isdir:
        try:
            contents = os.listdir(path)
            for entry in contents:
                # ignore hidden files/folders
                if not entry or entry.startswith('.'):
                    continue
                ext = os.path.splitext(entry)[1].lower()
                if ext in known_file_extensions:
                    files.append(entry)
                elif os.path.isdir(os.path.join(path, entry)):
                    dirs.append(entry)
                else:
                    files.append(entry)
        except OSError:
            error = True
            errortext = "permission denied"

    files.sort()
    dirs.sort()
    context = {
        "type": type,
        "purpose": purpose,
        "path": path,
        "parentdir": parentdir,
        "error": error,
        "errortext": errortext,
        "files": files,
        "dirs": dirs,
        "remember_directory": path if not error and path_isdir and path_was_requested else "",
    }
    response = render(request, template, context)
    return response
