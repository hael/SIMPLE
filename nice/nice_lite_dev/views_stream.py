'''
contains views associated with stream view
'''

# global imports
import os
import time
import pathlib
from django.shortcuts               import redirect
from django.contrib                 import messages
from django.http                    import HttpResponse, HttpResponseRedirect
from django.urls                    import reverse
from django.contrib.auth.decorators import login_required

# local imports
from .helpers                    import *
from .app_views.streamindex      import StreamIndex
from .app_views.datasetview      import DatasetView
from .app_views.newstreamview    import NewStreamView
from .app_views.streamview       import StreamView, StreamViewMovies, StreamViewPreprocess, StreamViewOptics, StreamViewInitialPick, StreamViewGeneratePickrefs
from .app_views.streamview       import StreamViewReferencePicking, StreamViewSieveParticles, StreamViewClassification2D, StreamViewParticleSets
from .data_structures.project    import Project
from .data_structures.dataset    import Dataset
from .data_structures.workspace  import Workspace
from .data_structures.streamjob  import StreamJob
from .data_structures.jobclassic import JobClassic

@login_required(login_url="/login")
def stream_index(request):
    '''
    returns stream index view
    '''
    projectid   = get_project_id(request)
    datasetid   = get_dataset_id(request)
    username    = request.user.username
    streamindex = StreamIndex(request, projectid, datasetid, username)
    return streamindex.render()

@login_required(login_url="/login/")
def delete_dataset(request):
    '''
    deletes a dataset from project and returns redirect to reload page
    '''
    deletedatasetid = get_integer(request.POST, "delete_dataset_id")
    dataset         = Dataset(deletedatasetid)
    if dataset.delete():
        messages.add_message(request, messages.INFO, "dataset deleted successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to delete dataset")
    response = redirect('nice_lite:stream')
    return response

@login_required(login_url="/login/")
def dataset(request):
    '''
    returns dataset view for dataset set in selected_dataset_id cookie
    '''
    datasetid   = get_dataset_id(request)
    projectid   = get_project_id(request)
    dataset     = Dataset(datasetid)
    datasetview = DatasetView(request, dataset, projectid)
    return datasetview.render()

@login_required(login_url="/login/")
def new_stream(request):
    '''
    returns create stream page
    '''
    jobid = get_job_id(request)
    if jobid is None:
        newstreamview = NewStreamView(request)
    else:
        streamjob      = StreamJob(jobid)
        streamjobmodel = streamjob.get_jobmodel()
        newstreamview  = NewStreamView(request, streamjobmodel.args)
    return newstreamview.render()

@login_required(login_url="/login/")
def create_stream(request):
    '''
    starts a new stream and refreshes page by redirecting to dataset view
    '''
    args = {}
    for key, value in request.POST.items():
        if "csrfmiddlewaretoken" not in key and value != "":
            args[key] = value
    datasetid = get_dataset_id(request)
    dataset   = Dataset(datasetid)
    streamjob = StreamJob()
    if not streamjob.new(dataset, args):
        print_error("failed to create new stream job")
    response = redirect('nice_lite:view_stream', jobid=streamjob.id)
    return response

@login_required(login_url="/login/")
def terminate_stream(request):
    '''
    terminates stream and refreshes page by redirecting to dataset view
    '''
    jobid = get_job_id(request)
    streamjob = StreamJob(jobid)
    streamjob.terminate_master()
    response = redirect('nice_lite:dataset')
    return response

@login_required(login_url="/login/")
def terminate_stream_process(request):
    '''
    terminate stream process and refreshes panel by redirect
    '''
    jobid                  = get_job_id(request)
    term_preprocess        = string_present(request.POST,        "terminate_preprocess", silent=True)
    term_optics_assignment = string_present(request.POST, "terminate_optics_assignment", silent=True)
    term_generate_pickrefs = string_present(request.POST, "terminate_generate_pickrefs", silent=True)
    streamjob              = StreamJob(jobid)
    streamjob.terminate_process(term_preprocess, term_optics_assignment, term_generate_pickrefs)
    if term_preprocess:
        response = redirect('nice_lite:view_stream_preprocess', jobid=jobid)
    elif term_optics_assignment:
        response = HttpResponseRedirect(reverse('nice_lite:view_stream_optics', query={"selected_job_id" : jobid}))    #reverse("nice_lite:dataset", query={"selected_dataset_id" : self.datasetid})
    elif term_generate_pickrefs:
        response = HttpResponseRedirect(reverse('nice_lite:view_stream_generate_pickrefs', query={"selected_job_id" : jobid}))
    return response

@login_required(login_url="/login/")
def restart_stream_process(request):
    '''
    restart stream process and refresh panel by redirect
    '''
    jobid                     = get_job_id(request)
    restart_preprocess        = string_present(request.POST,        "restart_preprocess", silent=True)
    restart_optics_assignment = string_present(request.POST, "restart_optics_assignment", silent=True)
    restart_generate_pickrefs = string_present(request.POST, "restart_generate_pickrefs", silent=True)
    streamjob                 = StreamJob(jobid)
    streamjob.restart_process(restart_preprocess, restart_optics_assignment, restart_generate_pickrefs)
    if restart_preprocess:
        response = redirect('nice_lite:view_stream_preprocess', jobid=jobid)
    elif restart_optics_assignment:
        response = HttpResponseRedirect(reverse('nice_lite:view_stream_optics', query={"selected_job_id" : jobid}))
    elif restart_generate_pickrefs:
        response = HttpResponseRedirect(reverse('nice_lite:view_stream_generate_pickrefs', query={"selected_job_id" : jobid}))
    return response

@login_required(login_url="/login/")
def delete_stream(request):
    '''
    deletes stream and refreshes page by redirecting to dataset view
    '''
    jobid     = get_job_id(request)
    streamjob = StreamJob(jobid)
    streamjob.delete()
    response = redirect('nice_lite:dataset')
    return response

@login_required(login_url="/login/")
def view_stream(request, jobid):
    '''
    returns stream view 
    '''
    streamview = StreamView(request, jobid)
    return streamview.render()

@login_required(login_url="/login/")
def view_stream_movies(request):
    '''
    returns movies panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewmovies = StreamViewMovies(request, jobid)
    return streamviewmovies.render()

@login_required(login_url="/login/")
def view_stream_preprocess(request):
    '''
    returns preprocess panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewpreprocess = StreamViewPreprocess(request, jobid, None)
    return streamviewpreprocess.render()

@login_required(login_url="/login/")
def view_stream_preprocess_zoom(request):
    '''
    returns preprocess zoom panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewpreprocess = StreamViewPreprocess(request, None, jobid)
    return streamviewpreprocess.render()

@login_required(login_url="/login/")
def view_stream_optics(request):
    '''
    returns optics panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewoptics = StreamViewOptics(request, jobid, None)
    return streamviewoptics.render()

@login_required(login_url="/login/")
def view_stream_optics_zoom(request):
    '''
    returns optics zoom panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewoptics = StreamViewOptics(request, None, jobid)
    return streamviewoptics.render()

@login_required(login_url="/login/")
def view_stream_initial_pick(request):
    '''
    returns initial picking panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewinitialpick = StreamViewInitialPick(request, jobid, None)
    return streamviewinitialpick.render()

@login_required(login_url="/login/")
def view_stream_initial_pick_zoom(request):
    '''
    returns initial picking panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewinitialpick = StreamViewInitialPick(request, None, jobid)
    return streamviewinitialpick.render()

@login_required(login_url="/login/")
def view_stream_generate_pickrefs(request):
    '''
    returns reference generation panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewgeneratepickrefs = StreamViewGeneratePickrefs(request, jobid, None)
    return streamviewgeneratepickrefs.render()

@login_required(login_url="/login/")
def view_stream_generate_pickrefs_zoom(request):
    '''
    returns reference generation panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewgeneratepickrefs = StreamViewGeneratePickrefs(request, None, jobid)
    return streamviewgeneratepickrefs.render()

@login_required(login_url="/login/")
def view_stream_reference_picking(request):
    '''
    returns reference picking panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewreferencepicking = StreamViewReferencePicking(request, jobid, None)
    return streamviewreferencepicking.render()

@login_required(login_url="/login/")
def view_stream_reference_picking_zoom(request):
    '''
    returns reference picking panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewreferencepicking = StreamViewReferencePicking(request, None, jobid)
    return streamviewreferencepicking.render()

@login_required(login_url="/login/")
def view_stream_sieve_particles(request):
    '''
    returns particle sieving panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewsieveparticles = StreamViewSieveParticles(request, jobid, None)
    return streamviewsieveparticles.render()

@login_required(login_url="/login/")
def view_stream_sieve_particles_zoom(request):
    '''
    returns particle sieving panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewsieveparticles = StreamViewSieveParticles(request, None, jobid)
    return streamviewsieveparticles.render()

@login_required(login_url="/login/")
def view_stream_classification_2D(request):
    '''
    returns 2D classification panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewclassification2D = StreamViewClassification2D(request, jobid, None)
    return streamviewclassification2D.render()

@login_required(login_url="/login/")
def view_stream_classification_2D_zoom(request):
    '''
    returns 2D classification panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewclassification2D = StreamViewClassification2D(request, None, jobid)
    return streamviewclassification2D.render()

@login_required(login_url="/login/")
def view_stream_particle_sets(request):
    '''
    returns particle sets panel in stream view
    '''
    jobid = get_job_id(request)
    streamviewparticlesets = StreamViewParticleSets(request, jobid)
    return streamviewparticlesets.render()

@login_required(login_url="/login/")
def regenerate_pickrefs(request):
    jobid = get_job_id(request)
    streamjob = StreamJob(jobid)
    streamjob.update(increase_nmics=True)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_pickrefs(request):
    jobid = get_job_id(request)
    streamjob = StreamJob(jobid)
    final_selection         = [int(numeric_string) for numeric_string in request.POST["final_selection"].split(',')]
    streamjob.select_pickrefs(final_selection)
    response = redirect('nice_lite:view_stream', jobid=jobid)
    return response

@login_required(login_url="/login/")
def select_stream_sieve_particles(request):
    jobid = get_job_id(request)
    streamjob = StreamJob(jobid)
    accepted_cls2D = [int(numeric_string) for numeric_string in request.POST["accepted_cls2D"].split(',')]
    streamjob.select_sieve_particles(accepted_cls2D)
    streamviewsieveparticles = StreamViewSieveParticles(request, None, jobid)
    return streamviewsieveparticles.render()

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
def update_stream_description(request):
    '''
    update description for give stream jobid
    '''
    jobid           = get_job_id(request)
    new_description = get_string(request.POST, "new_stream_description")
    if jobid is not None:
        streamjob = StreamJob(jobid)
        streamjob.update_description(new_description)
    return HttpResponseNoContent()
    
@login_required(login_url="/login/")
def update_stream_parameters(request):
    jobid = get_job_id(request)
    if jobid is not None:
        streamjob    = StreamJob(jobid)
        ctfres       = get_float(request.POST, "ctfres",      silent=True)
        astigmatism  = get_float(request.POST, "astigmatism", silent=True)
        icescore     = get_float(request.POST, "icescore",    silent=True)
        streamjob.update(ctfres, astigmatism, icescore)
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
    datasetid      = get_dataset_id(request)
    dataset        = Dataset(datasetid)
    newdatasetname = get_string(request.POST, "new_dataset_name")
    if dataset.rename(newdatasetname):
        messages.add_message(request, messages.INFO, "dataset renamed successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to rename dataset")
    response = redirect('nice_lite:stream')
    return response

@login_required(login_url="/login/")
def update_dataset_description(request):
    datasetid             = get_dataset_id(request)
    dataset               = Dataset(datasetid)
    newdatasetdescription = get_string(request.POST, "new_dataset_description")
    if dataset.updateDescription(newdatasetdescription):
        messages.add_message(request, messages.INFO, "description updated successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to update description")
    response = redirect('nice_lite:stream')
    return response