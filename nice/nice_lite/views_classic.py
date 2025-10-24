import time
# global imports
from django.shortcuts               import redirect
from django.http                    import HttpResponse
from django.contrib.auth.decorators import login_required

# local imports
from .app_views.coreview       import CoreViewClassic
from .app_views.workspaceview  import WorkspaceView
from .app_views.jobview        import JobView, JobViewMicrographs, JobViewMicrographsHistogram, JobViewCls2D, JobViewCls2DHistogram
from .app_views.jobview        import JobViewLogs, JobViewMicrographsPlot, JobViewCls2DPlot
from .app_views.newjobview     import NewJobView, NewJobTypeView

from .data_structures.project    import Project
from .data_structures.workspace  import Workspace
from .data_structures.jobclassic import JobClassic

@login_required(login_url="/login")
def classic(request):
    coreviewclassic = CoreViewClassic(request)
    return coreviewclassic.render()

@login_required(login_url="/login/")
def workspace(request):
    workspaceview = WorkspaceView(request)
    return workspaceview.render()

@login_required(login_url="/login/")
def view_job(request, jobid):
    jobview = JobView(request, jobid)
    return jobview.render()

@login_required(login_url="/login/")
def view_job_micrographs(request, jobid):
    sort_micrographs_key = None
    if "sort_micrographs_key" in request.POST:
        sort_micrographs_key = request.POST["sort_micrographs_key"]
    sort_micrographs_asc = True
    if "sort_micrographs_asc" in request.POST and request.POST["sort_micrographs_asc"] == "false":
        sort_micrographs_asc = False 
    fromp = None
    top   = None
    if "fromp" in request.POST and "top" in request.POST:
        fromp = request.POST["fromp"]
        top   = request.POST["top"]
    page = 1
    if "page" in request.POST:
        page = request.POST["page"]
    jobviewmicrographs = JobViewMicrographs(request, jobid, sort_micrographs_key, sort_micrographs_asc)
    return jobviewmicrographs.render(fromp, top, page)

@login_required(login_url="/login/")
def view_job_micrographs_histogram(request, jobid):
    sort_micrographs_key = None
    if "sort_micrographs_key" in request.POST:
        sort_micrographs_key = request.POST["sort_micrographs_key"]
    jobviewmicrographs = JobViewMicrographsHistogram(request, jobid, sort_micrographs_key)
    return jobviewmicrographs.render()

@login_required(login_url="/login/")
def view_job_micrographs_plot(request, jobid):
    sort_micrographs_key = None
    if "sort_micrographs_key" in request.POST:
        sort_micrographs_key = request.POST["sort_micrographs_key"]
    plot_micrographs_key = None
    if "plot_micrographs_key" in request.POST:
        plot_micrographs_key = request.POST["plot_micrographs_key"]  
    jobviewmicrographs = JobViewMicrographsPlot(request, jobid, sort_micrographs_key, plot_micrographs_key)
    return jobviewmicrographs.render()

@login_required(login_url="/login/")
def view_job_cls2D(request, jobid):
    sort_cls2d_key = None
    if "sort_cls2d_key" in request.POST:
        sort_cls2d_key = request.POST["sort_cls2d_key"]
    sort_cls2d_asc = True
    if "sort_cls2d_asc" in request.POST and request.POST["sort_cls2d_asc"] == "false":
        sort_cls2d_asc = False 
    jobviewcls2d = JobViewCls2D(request, jobid, sort_cls2d_key, sort_cls2d_asc)
    return jobviewcls2d.render()

@login_required(login_url="/login/")
def view_job_cls2D_histogram(request, jobid):
    sort_cls2d_key = None
    if "sort_cls2d_key" in request.POST:
        sort_cls2d_key = request.POST["sort_cls2d_key"]
    jobviewcls2d = JobViewCls2DHistogram(request, jobid, sort_cls2d_key)
    return jobviewcls2d.render()

@login_required(login_url="/login/")
def view_job_cls2D_plot(request, jobid):
    sort_cls2d_key = None
    if "sort_cls2d_key" in request.POST:
        sort_cls2d_key = request.POST["sort_cls2d_key"]
    plot_cls2d_key = None
    if "plot_cls2d_key" in request.POST:
        plot_cls2d_key = request.POST["plot_cls2d_key"]  
    jobviewcls2d = JobViewCls2DPlot(request, jobid, sort_cls2d_key, plot_cls2d_key)
    return jobviewcls2d.render()

@login_required(login_url="/login/")
def view_job_logs(request, jobid):
    jobviewlogs = JobViewLogs(request, jobid)
    return jobviewlogs.render()

@login_required(login_url="/login/")
def select_job_micrographs(request, jobid):
    time.sleep(1) # gives user time to see message on button
    project      = Project(request=request)
    workspace    = Workspace(request=request)
    selectionjob = JobClassic()
    selectionjob.newSelection(request, project, workspace, jobid)
    response = redirect('nice_lite:view_job_micrographs', jobid)
    return response

@login_required(login_url="/login/")
def select_job_cls2D(request, jobid):
    time.sleep(1) # gives user time to see message on button
    project      = Project(request=request)
    workspace    = Workspace(request=request)
    selectionjob = JobClassic()
    selectionjob.newSelection(request, project, workspace, jobid)
    response = redirect('nice_lite:view_job_cls2D', jobid)
    return response

@login_required(login_url="/login/")
def new_job_type(request, parentid):
    newjobtypeview = NewJobTypeView(request, parentid)
    return newjobtypeview.render()

@login_required(login_url="/login/")
def new_job(request, parentid, package, jobtype):
    newjobview = NewJobView(request, parentid, package, jobtype)
    return newjobview.render()

@login_required(login_url="/login/")
def rerun_job(request, jobid):
    job = JobClassic(id=jobid)
    newjobview = NewJobView(request, job.prnt, job.pckg, job.prog, args=job.args)
    return newjobview.render()

@login_required(login_url="/login/")
def create_workspace(request, projectid):
    project = Project(project_id=projectid)
    workspace = Workspace()
    workspace.new(project, request.user.username)
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
    project   = Project(request=request)
    workspace = Workspace(request=request)
    job = JobClassic(id=jobid)
    job.markComplete(project, workspace)
    response = redirect('nice_lite:workspace')
    return response

@login_required(login_url="/login/")
def delete_job(request, jobid):
    project   = Project(request=request)
    workspace = Workspace(request=request)
    job = JobClassic(id=jobid)
    job.delete(project, workspace)
    response = redirect('nice_lite:workspace')
    return response