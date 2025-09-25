# global imports
from django.shortcuts               import redirect, render
from django.http                    import HttpResponse
from django.contrib.auth.decorators import login_required
from django.contrib.auth            import login, logout, authenticate

# local imports
from .app_views.coreview       import CoreViewClassic
from .app_views.newprojectview import NewProjectView
from .app_views.workspaceview  import WorkspaceView
from .app_views.newjobview     import NewJobView, NewJobTypeView

from .data_structures.project  import Project
from .data_structures.workspace  import Workspace
from .data_structures.jobclassic      import JobClassic

@login_required(login_url="/login")
def classic(request):
    coreviewclassic = CoreViewClassic(request)
    return coreviewclassic.render()

@login_required(login_url="/login/")
def workspace(request):
    workspaceview = WorkspaceView(request)
    return workspaceview.render()

@login_required(login_url="/login/")
def new_job_type(request, parentid):
    newjobtypeview = NewJobTypeView(request, parentid)
    return newjobtypeview.render()

@login_required(login_url="/login/")
def new_job(request, parentid, package, jobtype):
    newjobview = NewJobView(request, parentid, package, jobtype)
    return newjobview.render()

@login_required(login_url="/login/")
def rerun_job(request, parentid):
    job = JobClassic(id=parentid)
    newjobview = NewJobView(request, parentid, job.pckg, job.prog, args=job.args)
    return newjobview.render()

@login_required(login_url="/login/")
def create_workspace(request, projectid):
    project = Project(project_id=projectid)
    workspace = Workspace()
    workspace.new(project)
    response = redirect('nice_lite:classic')
    response.set_cookie(key='selected_project_id',   value=project.id)
    response.set_cookie(key='selected_workspace_id', value=workspace.id)
    return response
    
@login_required(login_url="/login/")
def delete_workspace(request, workspaceid):
    project = Project(request=request)
    workspace = Workspace(workspace_id=workspaceid)
    workspace.delete(project)
    response = redirect('nice_lite:classic')
    return response

@login_required(login_url="/login/")
def create_job(request, parentid, package, jobtype):
    project   = Project(request=request)
    workspace = Workspace(request=request)
    job = JobClassic()
    job.new(request, project, workspace, parentid, package, jobtype)
    response = redirect('nice_lite:workspace')
    return response

@login_required(login_url="/login/")
def update_workspace_name(request):
    project   = Project(request=request)
    workspace = Workspace(request=request)
    workspace.rename(request, project)
    response = redirect('nice_lite:classic')
    return response

@login_required(login_url="/login/")
def update_workspace_description(request):
    workspace = Workspace(request=request)
    workspace.updateDescription(request)
    response = redirect('nice_lite:classic')
    return response

@login_required(login_url="/login/")
def update_job_description(request, jobid):
    job = JobClassic(id=jobid)
    job.updateDescription(request)
    response = redirect('nice_lite:workspace')
    return response

@login_required(login_url="/login/")
def mark_job_complete(request, jobid):
    job = JobClassic(id=jobid)
    job.markComplete()
    response = redirect('nice_lite:workspace')
    return response