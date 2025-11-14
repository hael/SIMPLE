'''
contains views associated with stream view
'''

# global imports
import os
import time
import pathlib
from django.shortcuts               import redirect
from django.http                    import HttpResponse
from django.contrib.auth.decorators import login_required

# local imports
from .app_views.coreview         import CoreViewStream
from .app_views.datasetview      import DatasetView
from .app_views.newstreamview    import NewStreamView
from .app_views.streamview       import StreamView, StreamViewMovies, StreamViewPreprocess, StreamViewOptics, StreamViewInitialPick, StreamViewGeneratePickrefs
from .app_views.streamview       import StreamViewReferencePicking, StreamViewSieveParticles, StreamViewClassification2D, StreamViewParticleSets
from .data_structures.project    import Project
from .data_structures.dataset    import Dataset
from .data_structures.workspace  import Workspace
from .data_structures.job        import Job
from .data_structures.jobclassic import JobClassic

@login_required(login_url="/login")
def stream(request):
    '''
    returns stream core view
    '''
    coreviewstream = CoreViewStream(request)
    return coreviewstream.render()

@login_required(login_url="/login/")
def create_dataset(request):
    '''
    creates a new dataset, sets cookies and returns redirect to load
    '''
    project = Project(request=request)
    dataset = Dataset()
    dataset.new(project, request.user.username)
    response = redirect('nice_lite:stream')
    response.set_cookie(key='selected_project_id', value=project.id)
    response.set_cookie(key='selected_dataset_id', value=dataset.id)
    return response

@login_required(login_url="/login/")
def delete_dataset(request, datasetid):
    '''
    deletes a dataset from project and returns redirect to reload page
    '''
    project = Project(request=request)
    dataset = Dataset(dataset_id=datasetid)
    dataset.delete(project)
    response = redirect('nice_lite:stream')
    return response

@login_required(login_url="/login/")
def dataset(request):
    '''
    returns dataset view for dataset set in selected_dataset_id cookie
    '''
    datasetview = DatasetView(request)
    return datasetview.render()

@login_required(login_url="/login/")
def new_stream(request):
    '''
    returns create stream page
    '''
    newstreamview = NewStreamView(request)
    return newstreamview.render()

@login_required(login_url="/login/")
def rerun_stream(request, parentid):
    '''
    returns create stream page using arguments from previous stream job
    with id=parentid
    '''
    job = Job(id=parentid)
    newstreamview = NewStreamView(request, args=job.args)
    return newstreamview.render()

@login_required(login_url="/login/")
def terminate_stream(request, jobid):
    '''
    terminates stream and refreshes page by redirecting to dataset view
    '''
    job = Job(id=jobid)
    job.terminate()
    response = redirect('nice_lite:dataset')
    return response

@login_required(login_url="/login/")
def delete_stream(request, jobid):
    '''
    deletes stream and refreshes page by redirecting to dataset view
    '''
    job = Job(id=jobid)
    project = Project(request=request)
    dataset = Dataset(request=request)
    job.delete(project, dataset)
    response = redirect('nice_lite:dataset')
    return response

@login_required(login_url="/login/")
def create_stream(request):
    '''
    starts a new stream and refreshes page by redirecting to dataset view
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job()
    job.new(request, project, dataset)
    response = redirect('nice_lite:view_stream', jobid=job.id)
    return response

@login_required(login_url="/login/")
def view_stream(request, jobid):
    '''
    returns stream view 
    '''
    streamview = StreamView(request, jobid)
    return streamview.render()

@login_required(login_url="/login/")
def view_stream_movies(request, jobid=None, jobidzoom=None):
    '''
    returns movies panel in stream view
    '''
    streamviewmovies = StreamViewMovies(request, jobid, jobidzoom)
    return streamviewmovies.render()

@login_required(login_url="/login/")
def view_stream_preprocess(request, jobid=None, jobidzoom=None):
    '''
    returns preprocess panel in stream view
    '''
    streamviewpreprocess = StreamViewPreprocess(request, jobid, jobidzoom)
    return streamviewpreprocess.render()

@login_required(login_url="/login/")
def term_stream_preprocess(request, jobid):
    '''
    sends termination signal to preprocess process and reloads preprocess panel
    by redirection
    '''
    job = Job(id=jobid)
    job.terminate_preprocess()
    response = redirect('nice_lite:view_stream_preprocess', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_preprocess(request, jobid):
    '''
    starts terminated preprocess process and reloads preprocess panel
    by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_preprocess(project, dataset)
    response = redirect('nice_lite:view_stream_preprocess', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_optics(request, jobid=None, jobidzoom=None):
    '''
    returns optics panel in stream view
    '''
    streamviewoptics = StreamViewOptics(request, jobid, jobidzoom)
    return streamviewoptics.render()

@login_required(login_url="/login/")
def term_stream_optics(request, jobid):
    '''
    sends termination signal to optics process and reloads optics panel
    by redirection
    '''
    job = Job(id=jobid)
    job.terminate_optics()
    response = redirect('nice_lite:view_stream_optics', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_optics(request, jobid):
    '''
    starts terminated optics process and reloads optics panel
    by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_optics(project, dataset)
    response = redirect('nice_lite:view_stream_optics', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_initial_pick(request, jobid=None, jobidzoom=None):
    '''
    returns initial picking panel in stream view
    '''
    streamviewinitialpick = StreamViewInitialPick(request, jobid, jobidzoom)
    return streamviewinitialpick.render()

@login_required(login_url="/login/")
def term_stream_initial_pick(request, jobid):
    '''
    sends termination signal to initial picking process and reloads intitial picking
    panel by redirection
    '''
    job = Job(id=jobid)
    job.terminate_initial_pick()
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_initial_pick(request, jobid):
    '''
    starts terminated initial picking process and reloads intitial picking
    panel by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_initial_pick(project, dataset)
    response = redirect('nice_lite:view_stream_initial_pick', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    diameter = request.POST["diameter"]
    job.select_moldiam_initial_pick(diameter)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    return response

@login_required(login_url="/login/")
def refine_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    refine_diameter = request.POST["refine_diameter"]
    job.update_moldiam_refine_initial_pick(refine_diameter)
    response = redirect('nice_lite:view_stream_initial_pick_zoom', jobidzoom=jobid)
    return response

@login_required(login_url="/login/")
def increase_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    job.update_moldiam_refine_initial_pick(-1)
    response = redirect('nice_lite:view_stream_initial_pick_zoom', jobidzoom=jobid)
    return response

@login_required(login_url="/login/")
def decrease_moldiam_stream_initial_pick(request, jobid):
    job = Job(id=jobid)
    job.update_moldiam_refine_initial_pick(-2)
    response = redirect('nice_lite:view_stream_initial_pick_zoom', jobidzoom=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_generate_pickrefs(request, jobid=None, jobidzoom=None):
    '''
    returns reference generation panel in stream view
    '''
    streamviewgeneratepickrefs = StreamViewGeneratePickrefs(request, jobid, jobidzoom)
    return streamviewgeneratepickrefs.render()

@login_required(login_url="/login/")
def term_stream_generate_pickrefs(request, jobid):
    '''
    sends termination signal to pickrefs generation process and reloads pickrefs
    generation panel by redirection
    '''
    job = Job(id=jobid)
    job.terminate_generate_pickrefs()
    response = redirect('nice_lite:view_stream_generate_pickrefs', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_generate_pickrefs(request, jobid):
    '''
    restarts terminated pickrefs generation process and reloads pickrefs
    generation panel by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_generate_pickrefs(project, dataset)
    response = redirect('nice_lite:view_stream_generate_pickrefs', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_refs_stream_generate_pickrefs(request, jobid):
    job = Job(id=jobid)
    final_selection         = [int(numeric_string) for numeric_string in request.POST["final_selection"].split(',')]
    final_selection_source  = request.POST["final_selection_source"]
    job.select_refs_generate_pickrefs(final_selection, final_selection_source)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    return response

@login_required(login_url="/login/")
def stream_regenerate_pickrefs(request, jobid):
    job = Job(id=jobid)
    job.regenerate_pickrefs()
    response = redirect('nice_lite:view_stream', jobid=jobid)
    return response

@login_required(login_url="/login/")
def snapshot_stream_classification_2D(request, jobid):
    job = Job(id=jobid)
    snapshot_selection  = [int(numeric_string) for numeric_string in request.POST["snapshot_selection"].split(',')]
    snapshot_iteration = request.POST["snapshot_iteration"]
    job.snapshot_classification_2D(snapshot_selection, snapshot_iteration)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    time.sleep(2) #sleep for user to see message
    return response

@login_required(login_url="/login/")
def select_stream_classification_2D(request, jobid):
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    final_deselection     = [int(numeric_string) for numeric_string in request.POST["final_deselection"].split(',')]
    final_selection_ptcls = request.POST["final_selection_ptcls"]
    job.selection_classification_2D(final_deselection, project, dataset, final_selection_ptcls)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    time.sleep(2) #sleep for user to see message
    return response

@login_required(login_url="/login/")
def view_stream_reference_picking(request, jobid=None, jobidzoom=None):
    '''
    returns reference picking panel in stream view
    '''
    streamviewreferencepicking = StreamViewReferencePicking(request, jobid, jobidzoom)
    return streamviewreferencepicking.render()

@login_required(login_url="/login/")
def term_stream_reference_picking(request, jobid):
    '''
    sends termination signal to reference picking process and reloads 
    reference picking panel by redirection
    '''
    job = Job(id=jobid)
    job.terminate_reference_picking()
    response = redirect('nice_lite:view_stream_reference_picking', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_reference_picking(request, jobid):
    '''
    starts terminated reference picking process and reloads 
    reference picking panel by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_reference_picking(project, dataset)
    response = redirect('nice_lite:view_stream_reference_picking', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_sieve_particles(request, jobid=None, jobidzoom=None):
    '''
    returns particle sieving panel in stream view
    '''
    streamviewsieveparticles = StreamViewSieveParticles(request, jobid, jobidzoom)
    return streamviewsieveparticles.render()

@login_required(login_url="/login/")
def term_stream_sieve_particles(request, jobid):
    '''
    sends termination signal to particle sieving process and reloads particle sieving
    panel by redirection
    '''
    job = Job(id=jobid)
    job.terminate_sieve_particles()
    response = redirect('nice_lite:view_stream_sieve_particles', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_sieve_particles(request, jobid):
    '''
    starts terminated particle sieving process and reloads particle sieving
    panel by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_sieve_particles(project, dataset)
    response = redirect('nice_lite:view_stream_sieve_particles', jobid=jobid)
    return response

# @login_required(login_url="/login/")
def select_stream_sieve_particles(request, jobid):
    job = Job(id=jobid)
    accepted_cls2D = [int(numeric_string) for numeric_string in request.POST["accepted_cls2D"].split(',')]
    if len(request.POST["rejected_cls2D"]) > 0:
        rejected_cls2D = [int(numeric_string) for numeric_string in request.POST["rejected_cls2D"].split(',')]
    else:
        rejected_cls2D = ['0']
    job.select_sieve_particles(accepted_cls2D, rejected_cls2D)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_stream_classification_2D(request, jobid=None, jobidzoom=None):
    '''
    returns 2D classification panel in stream view
    '''
    streamviewclassification2D = StreamViewClassification2D(request, jobid, jobidzoom)
    return streamviewclassification2D.render()

@login_required(login_url="/login/")
def term_stream_classification_2D(request, jobid):
    '''
    sends termination signal to 2D classification process and reloads 2D classification panel
    by redirection
    '''
    job = Job(id=jobid)
    job.terminate_classification_2D()
    response = redirect('nice_lite:view_stream_classification_2D', jobid=jobid)
    return response

@login_required(login_url="/login/")
def restart_stream_classification_2D(request, jobid):
    '''
    starts terminated 2D classification process and reloads 2D classification panel
    by redirection
    '''
    project = Project(request=request)
    dataset = Dataset(request=request)
    job = Job(id=jobid)
    job.restart_classification_2D(project, dataset)
    response = redirect('nice_lite:view_stream_classification_2D', jobid=jobid)
    return response

@login_required(login_url="/login/")
def view_logs(request, jobid, log, error):
    streamviewlogs = StreamViewLogs(request, jobid, log, error)
    return streamviewlogs.render()

@login_required(login_url="/login/")
def view_stream_particle_sets(request, jobid):
    '''
    returns particle sets panel in stream view
    '''
    streamviewparticlesets = StreamViewParticleSets(request, jobid)
    return streamviewparticlesets.render()

@login_required(login_url="/login/")
def update_stream_description(request, jobid):
    '''
    update description for give stream jobid
    '''
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
    response = redirect('nice_lite:view_stream_classification_2D_zoom', jobidzoom=jobid)
    return response

@login_required(login_url="/login/")
def link_stream_particle_set(request, jobid, setid, filename, type):
    # this ties stream and classic together !!!
    link_workspace_id = int(request.POST["link_workspace_id"])
    streamjob  = Job(id=jobid)
    classicjob = JobClassic()
    project    = Project(request=request)
    workspace  = Workspace(workspace_id=link_workspace_id)
    dataset    = Dataset(request=request)
    if type == "snapshot":
        set_proj   = os.path.join(project.dirc, dataset.dirc, streamjob.dirc, "classification_2D", "snapshots", pathlib.Path(filename).stem, filename)
        classicjob.linkParticleSet(project, workspace, set_proj)
    elif type == "final":
        set_proj  = os.path.join(project.dirc, dataset.dirc, streamjob.dirc, "classification_2D", "stream_abinitio2D.simple")
        set_desel = os.path.join(project.dirc, dataset.dirc, streamjob.dirc, filename)
        classicjob.linkParticleSetFinal(project, workspace, set_proj, set_desel)
    classicjob.update_description("from " + dataset.name + "->" + str(streamjob.id) + " stream->particle set " + str(setid))    
    #response = redirect('nice_lite:view_stream_particle_sets', jobid=jobid)
    response = redirect('nice_lite:classic')
    response.set_cookie(key='selected_project_id',   value=project.id)
    response.set_cookie(key='selected_workspace_id', value=workspace.id)
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