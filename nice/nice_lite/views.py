'''
contains common views to both stream and classic views
'''

# global imports
from django.shortcuts               import redirect, render
from django.http                    import HttpResponse
from django.contrib.auth.decorators import login_required
from django.contrib.auth            import login, logout, authenticate

# local imports
from .app_views.coreview         import LoginView
from .app_views.newprojectview   import NewProjectView
from .app_views.filebrowserview  import FileBrowserView
from .data_structures.project    import Project
from .data_structures.dataset    import Dataset
from .data_structures.workspace  import Workspace

@login_required(login_url="/login")
def index(request):
    ''' 
    return post-authentication select stream or classic page. Only displayed
    when mode cookie is unset
    '''
    mode = request.COOKIES.get('mode', 'none')
    if mode == "classic":
        response = redirect('nice_lite:classic')
    elif mode == "stream":
        response = redirect('nice_lite:stream')   
    else:
        response = render(request, "index.html", {})
    return response

def user_login(request):
    '''
    if user authentication fails return login page else
    redirect to stream
    '''
    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username = username, password = password)
        if user is not None:
            login(request, user)
            response = redirect('nice_lite:stream')
            return response
    loginview = LoginView(request)
    return loginview.render()

def user_logout(request):
    '''
    logout authenticated user
    '''
    logout(request)
    return HttpResponse("SUCCESSFUL LOGOUT.")

@login_required(login_url="/login/")
def new_project(request, caller):
    '''
    returns create new project page
    '''
    newprojectview = NewProjectView(request, caller)
    return newprojectview.render()

@login_required(login_url="/login/")
def file_browser(request, type, path=None):
    '''
    returns file browser page
    '''
    if path is None:
        if "selectedpath" in request.GET:
            path = request.GET['selectedpath']
        elif request.COOKIES.get('selected_project_id', 'none') != 'none' :
            project = Project(request=request)
            path = project.dirc
        else:
            path = "/"
    filebrowserview = FileBrowserView(request, type, path)
    return filebrowserview.render()

@login_required(login_url="/login/")
def create_project(request):
    '''
    creates new project containing a new empty dataset and workspace. 
    sets cookies to new project, workspace and dataset and redirects to 
    stream view to load
    '''
    project = Project()
    project.new(request)
    dataset = Dataset()
    dataset.new(project, request.user.username)
    workspace = Workspace()
    workspace.new(project, request.user.username)
    response = redirect('nice_lite:stream')
    response.set_cookie(key='selected_project_id',   value=project.id)
    response.set_cookie(key='selected_dataset_id',   value=dataset.id)
    response.set_cookie(key='selected_workspace_id', value=workspace.id)
    return response
