import os

from django.shortcuts               import render
from django.contrib.auth.decorators import login_required

from ..data_structures.project   import Project


@login_required(login_url="/login")
def view_file_browser(request, type, path=None):
    """
    Render the file browser page.

    Path resolution order:
      1. URL path segment (passed directly as an argument)
      2. 'selectedpath' GET parameter
      3. Root directory of the currently selected project (from cookie)
      4. Filesystem root '/'
    """
    template              = "filebrowser.html"
    known_file_extensions = [ ".spi", ".mrc", ".mrcs", ".dm3", ".bin", ".gain", ".eer", ".tiff", ".tif", ".jpeg", 
                              ".jpg", ".txt", ".simple", ".dat", ".img", ".map", ".head", ".ctf", ".raw", ".sbin",
                              ".dbin", ".asc", ".box", ".dat", ".pdb", ".star", ".hdf", ".pdf", ".ps"]
    path_isdir  = None
    parentdir   = None
    error       = False
    errortext   = ""
    files       = []
    dirs        = []

    if path is None:
        if "selectedpath" in request.GET:
            path = request.GET['selectedpath']
        elif request.COOKIES.get('selected_project_id', 'none') != 'none':
            project = Project(id=int(request.COOKIES['selected_project_id']))
            path = project.absdir or "/"
        else:
            path = "/"
    
    path = path.replace('//', '/') # deal with multiple / in path
    if path[0] != '/':
        path = '/' + path # fix missing leading / when using proxy
    try:
        if not os.path.exists(path):
            error     = True
            errortext = "path does not exist"
        if os.path.isdir(path):
            path_isdir = True
            parentdir  = os.path.dirname(path)
        else:
            path_isdir = False
    except OSError:
        error     = True
        errortext = "error"

    if not error and path_isdir:
        try:
            contents = os.listdir(path)
            for entry in contents:
                # ignore hidden files/folders
                if entry[0] == '.':
                    continue
                ext = os.path.splitext(entry)[1].lower()
                if ext in known_file_extensions:
                    files.append(entry)
                elif os.path.isdir(os.path.join(path, entry)):
                    dirs.append(entry)
                else:
                    files.append(entry)
        except OSError:
            error     = True
            errortext = "permission denied"
        
    files.sort()
    dirs.sort()
    context = {
        "type"      : type,
        "path"      : path,
        "parentdir" : parentdir,
        "error"     : error,
        "errortext" : errortext,
        "files"     : files,
        "dirs"      : dirs
    }
    response = render(request, template, context)
    return response
