# global imports
import os
import time
import pathlib
from django.shortcuts               import redirect, render
from django.http                    import HttpResponse
from django.contrib.auth.decorators import login_required
from django.contrib.auth            import login, logout, authenticate

# local imports
from .app_views.coreview         import CoreViewStream, LoginView
from .app_views.newprojectview   import NewProjectView
from .app_views.datasetview      import DatasetView
from .app_views.newstreamview    import NewStreamView
from .app_views.streamview       import StreamView, StreamViewLogs, StreamViewMovies, StreamViewPreprocess, StreamViewOptics, StreamViewInitialPick, StreamViewGeneratePickrefs
from .app_views.streamview       import StreamViewReferencePicking, StreamViewSieveParticles, StreamViewClassification2D, StreamViewParticleSets
from .data_structures.project    import Project
from .data_structures.dataset    import Dataset
from .data_structures.workspace  import Workspace
from .data_structures.job        import Job
from .data_structures.jobclassic import JobClassic

@login_required(login_url="/login")
def index(request):
    mode = request.COOKIES.get('mode', 'none')
    if mode == "classic":
        response = redirect('nice_lite:classic')
    elif mode == "stream":
        response = redirect('nice_lite:stream')   
    else:
        response = render(request, "index.html", {})
    return response

@login_required(login_url="/login")
def stream(request):
    coreviewstream = CoreViewStream(request)
    return coreviewstream.render()

@login_required(login_url="/login/")
def new_project(request, caller):
    newprojectview = NewProjectView(request, caller)
    return newprojectview.render()

@login_required(login_url="/login/")
def create_project(request):
    project = Project()
    project.new(request)
    dataset = Dataset()
    dataset.new(project)
    response = redirect('nice_lite:stream')
    response.set_cookie(key='selected_project_id', value=project.id)
    response.set_cookie(key='selected_dataset_id', value=dataset.id)
    return response

@login_required(login_url="/login/")
def create_dataset(request):
    project = Project(request=request)
    dataset = Dataset()
    dataset.new(project)
    response = redirect('nice_lite:stream')
    response.set_cookie(key='selected_project_id', value=project.id)
    response.set_cookie(key='selected_dataset_id', value=dataset.id)
    return response

@login_required(login_url="/login/")
def delete_dataset(request, datasetid):
    project = Project(request=request)
    dataset = Dataset(dataset_id=datasetid)
    dataset.delete(project)
    response = redirect('nice_lite:stream')
    return response

@login_required(login_url="/login/")
def dataset(request):
    datasetview = DatasetView(request)
    return datasetview.render()

@login_required(login_url="/login/")
def new_stream(request):
    newstreamview = NewStreamView(request)
    return newstreamview.render()

@login_required(login_url="/login/")
def rerun_stream(request, parentid):
    job = Job(id=parentid)
    newstreamview = NewStreamView(request, args=job.args)
    return newstreamview.render()

@login_required(login_url="/login/")
def terminate_stream(request, jobid):
    job = Job(id=jobid)
    job.terminate()
    response = redirect('nice_lite:dataset')
    return response

@login_required(login_url="/login/")
def delete_stream(request, jobid):
    job = Job(id=jobid)
    job.delete()
    response = redirect('nice_lite:dataset')
    return response

@login_required(login_url="/login/")
def create_stream(request):
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job()
    job.new(request, project, dataset)
    response = redirect('nice_lite:view_stream', jobid=job.id)
    return response

@login_required(login_url="/login/")
def view_stream(request, jobid):
    streamview = StreamView(request, jobid)
    return streamview.render()

@login_required(login_url="/login/")
def view_stream_movies(request, jobid=None, jobidzoom=None):
    streamviewmovies = StreamViewMovies(request, jobid, jobidzoom)
    return streamviewmovies.render()

@login_required(login_url="/login/")
def view_stream_preprocess(request, jobid=None, jobidzoom=None):
    streamviewpreprocess = StreamViewPreprocess(request, jobid, jobidzoom)
    return streamviewpreprocess.render()

@login_required(login_url="/login/")
def term_stream_preprocess(request, jobid):
    job = Job(id=jobid)
    job.terminate_preprocess()
    response = redirect('nice_lite:view_stream_preprocess', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_optics(request, jobid):
    streamviewoptics = StreamViewOptics(request, jobid)
    return streamviewoptics.render()

@login_required(login_url="/login/")
def term_stream_optics(request, jobid):
    job = Job(id=jobid)
    job.terminate_optics()
    response = redirect('nice_lite:view_stream_optics', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_initial_pick(request, jobid):
    streamviewinitialpick = StreamViewInitialPick(request, jobid)
    return streamviewinitialpick.render()

@login_required(login_url="/login/")
def term_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    job.terminate_intial_pick()
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    diameter = request.POST["diameter"]
    job.select_moldiam_initial_pick(diameter)
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def refine_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    refine_diameter = request.POST["refine_diameter"]
    job.update_moldiam_refine_initial_pick(refine_diameter)
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def increase_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    job.update_moldiam_refine_initial_pick(-1)
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def decrease_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    job.update_moldiam_refine_initial_pick(-2)
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_generate_pickrefs(request, jobid):
    streamviewgeneratepickrefs = StreamViewGeneratePickrefs(request, jobid)
    return streamviewgeneratepickrefs.render()

@login_required(login_url="/login/")
def term_stream_generate_pickrefs(request, jobid):
    job = Job(id=jobid)
    job.terminate_generate_pickrefs()
    response = redirect('nice_lite:view_stream_generate_pickrefs', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_refs_stream_generate_pickrefs(request, jobid):
    job = Job(id=jobid)
    final_selection         = [int(numeric_string) for numeric_string in request.POST["final_selection"].split(',')]
    final_selection_source  = request.POST["final_selection_source"]
    final_selection_boxsize = request.POST["final_selection_boxsize"]
    job.select_refs_generate_pickrefs(final_selection, final_selection_source, final_selection_boxsize)
    response = redirect('nice_lite:view_stream_generate_pickrefs', jobid=jobid)
    return response

@login_required(login_url="/login/")
def snapshot_stream_classification_2D(request, jobid):
    job = Job(id=jobid)
    snapshot_selection  = [int(numeric_string) for numeric_string in request.POST["snapshot_selection"].split(',')]
    snapshot_iteration = request.POST["snapshot_iteration"]
    job.snapshot_classification_2D(snapshot_selection, snapshot_iteration)
    response = redirect('nice_lite:view_stream_classification_2D', jobid=jobid)
    time.sleep(2) #sleep for user to see message
    return response

@login_required(login_url="/login/")
def view_stream_reference_picking(request, jobid):
    streamviewreferencepicking = StreamViewReferencePicking(request, jobid)
    return streamviewreferencepicking.render()

@login_required(login_url="/login/")
def term_stream_reference_picking(request, jobid):
    job = Job(id=jobid)
    job.terminate_reference_picking()
    response = redirect('nice_lite:view_stream_reference_picking', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_sieve_particles(request, jobid):
    streamviewsieveparticles = StreamViewSieveParticles(request, jobid)
    return streamviewsieveparticles.render()

@login_required(login_url="/login/")
def term_stream_sieve_particles(request, jobid):
    job = Job(id=jobid)
    job.terminate_sieve_particles()
    response = redirect('nice_lite:view_stream_sieve_particles', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_stream_sieve_particles(request, jobid):
    job = Job(id=jobid)
    accepted_cls2D = [int(numeric_string) for numeric_string in request.POST["accepted_cls2D"].split(',')]
    rejected_cls2D = [int(numeric_string) for numeric_string in request.POST["rejected_cls2D"].split(',')]
    job.select_sieve_particles(accepted_cls2D, rejected_cls2D)
    response = redirect('nice_lite:view_stream_sieve_particles', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_classification_2D(request, jobid):
    streamviewclassification2D = StreamViewClassification2D(request, jobid)
    return streamviewclassification2D.render()

@login_required(login_url="/login/")
def term_stream_classification_2D(request, jobid):
    job = Job(id=jobid)
    job.terminate_classification_2D()
    response = redirect('nice_lite:view_stream_classification_2D', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_logs(request, jobid, log, error):
    streamviewlogs = StreamViewLogs(request, jobid, log, error)
    return streamviewlogs.render()

@login_required(login_url="/login/")
def view_stream_particle_sets(request, jobid):
    streamviewparticlesets = StreamViewParticleSets(request, jobid)
    return streamviewparticlesets.render()

def user_login(request):
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
    logout(request)
    return HttpResponse("SUCCESSFUL LOGOUT.")
    
@login_required(login_url="/login/")
def update_stream_description(request, jobid):
    job = Job(id=jobid)
    new_stream_description = request.POST["new_stream_description"]
    job.update_description(new_stream_description)
    return HttpResponse(status=204)
    
@login_required(login_url="/login/")
def update_preprocess_ctfres(request, jobid):
    job = Job(id=jobid)
    ctfres = request.POST["ctfres"]
    job.update_ctfres(ctfres)
    return HttpResponse(status=204)
    
@login_required(login_url="/login/")
def update_preprocess_astig(request, jobid):
    job = Job(id=jobid)
    astigmatism = request.POST["astigmatism"]
    job.update_astigmatism(astigmatism)
    return HttpResponse(status=204)

@login_required(login_url="/login/")
def update_preprocess_icescore(request, jobid):
    job = Job(id=jobid)
    icescore = request.POST["icescore"]
    job.update_icescore(icescore)
    return HttpResponse(status=204)
    
@login_required(login_url="/login/")
def update_classification_2D_mskdiam(request, jobid):
    time.sleep(1) #sleep for user to see message
    job = Job(id=jobid)
    mskdiam = request.POST["mskdiam"]
    job.update_mskdiam(mskdiam)
    response = redirect('nice_lite:view_stream_classification_2D', jobid=jobid)
    return response

@login_required(login_url="/login/")
def link_stream_particle_set(request, jobid, setid, filename):
    # this ties stream and classic together !!!
    link_workspace_id = int(request.POST["link_workspace_id"])
    streamjob  = Job(id=jobid)
    classicjob = JobClassic()
    project    = Project(request=request)
    workspace  = Workspace(workspace_id=link_workspace_id)
    dataset    = Dataset(request=request)
    set_proj   = os.path.join(project.dirc, dataset.dirc, streamjob.dirc, "classification_2D", "snapshots", pathlib.Path(filename).stem, filename)
    classicjob.linkParticleSet(project, workspace, set_proj)
    classicjob.update_description("from " + dataset.name + "->" + str(streamjob.id) + " stream->particle set " + str(setid))
    response = redirect('nice_lite:view_stream_particle_sets', jobid=jobid)
    return response

@login_required(login_url="/login/")
def update_dataset_name(request):
    project = Project(request=request)
    dataset = Dataset(request=request)
    dataset.rename(request, project)
    response = redirect('nice_lite:stream')
    return response

@login_required(login_url="/login/")
def update_dataset_description(request):
    dataset = Dataset(request=request)
    dataset.updateDescription(request)
    response = redirect('nice_lite:stream')
    return response