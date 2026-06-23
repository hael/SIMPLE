import os
import json
import pathlib
import hashlib

from django.contrib.auth.decorators import login_required
from django.http                    import HttpResponse, HttpResponseRedirect
from django.shortcuts               import redirect, render
from django.urls                    import reverse

from .job_builder_views     import view_job_builder_stream
from ..helpers              import *
from ..models               import WorkspaceModel

from ..data_structures.batchjob  import BatchJob
from ..data_structures.project   import Project
from ..data_structures.streamjob import StreamJob
from ..data_structures.workspace import Workspace


def _render_if_changed(request, template, context, checksum_cookie):
    """Render template only when checksum changed; otherwise return HTTP 204."""
    checksum = hashlib.md5(json.dumps(context, sort_keys=True).encode()).hexdigest()
    old_checksum = request.COOKIES.get(checksum_cookie, "none")
    if old_checksum == checksum:
        return HttpResponseNoContent()
    response = render(request, template, context)
    response.set_cookie(key=checksum_cookie, value=checksum)
    return response


def _get_jobmodel_from_request(request):
    """Return selected stream job model from request cookies/params, or None."""
    jobid = get_job_id(request)
    if jobid is None:
        return None
    streamjob = StreamJob(id=jobid)
    return streamjob.get_jobmodel()


def _get_jobmodel_and_dir_from_request(request):
    """Return selected stream job model and absolute directory, or (None, None)."""
    jobid = get_job_id(request)
    if jobid is None:
        return None, None
    streamjob = StreamJob(id=jobid)
    jobmodel = streamjob.get_jobmodel()
    if jobmodel is None:
        return None, None
    return jobmodel, streamjob.get_absdir()


@login_required(login_url="/login/")
def view_stream_create_stream(request):
    """Starts a new stream and refreshes page by redirecting to stream view."""
    args = {}
    for key, value in request.POST.items():
        if "csrfmiddlewaretoken" not in key and value != "":
            args[key] = value
    workspaceid = get_workspace_id(request)
    workspace_obj = Workspace(workspaceid)
    streamjob = StreamJob()
    if not streamjob.new(workspace_obj, args):
        print_error("failed to create new stream job")
    response = redirect("nice_lite:view_stream", jobid=streamjob.id)
    return response


@login_required(login_url="/login/")
def view_stream_terminate_stream(request):
    """Terminates stream and refreshes page by redirecting to workspace view."""
    jobid = get_job_id(request)
    streamjob = StreamJob(jobid)
    streamjob.terminate_master()
    response = redirect("nice_lite:workspace")
    return response


@login_required(login_url="/login/")
def view_stream_terminate_stream_process(request):
    """Terminate stream process and refresh panel by redirect."""
    jobid = get_job_id(request)
    term_preprocess = string_present(request.POST, "terminate_preprocess", silent=True)
    term_optics_assignment = string_present(request.POST, "terminate_optics_assignment", silent=True)
    term_generate_pickrefs = string_present(request.POST, "terminate_generate_pickrefs", silent=True)
    streamjob = StreamJob(jobid)
    streamjob.terminate_process(term_preprocess, term_optics_assignment, term_generate_pickrefs)
    if term_preprocess:
        response = redirect("nice_lite:view_stream_preprocess", jobid=jobid)
    elif term_optics_assignment:
        response = HttpResponseRedirect(reverse("nice_lite:view_stream_optics", query={"selected_job_id": jobid}))
    elif term_generate_pickrefs:
        response = HttpResponseRedirect(reverse("nice_lite:view_stream_generate_pickrefs", query={"selected_job_id": jobid}))
    return response


@login_required(login_url="/login/")
def view_stream_restart_stream_process(request):
    """Restart stream process and refresh panel by redirect."""
    jobid = get_job_id(request)
    restart_preprocess = string_present(request.POST, "restart_preprocess", silent=True)
    restart_optics_assignment = string_present(request.POST, "restart_optics_assignment", silent=True)
    restart_generate_pickrefs = string_present(request.POST, "restart_generate_pickrefs", silent=True)
    streamjob = StreamJob(jobid)
    streamjob.restart_process(restart_preprocess, restart_optics_assignment, restart_generate_pickrefs)
    if restart_preprocess:
        response = redirect("nice_lite:view_stream_preprocess", jobid=jobid)
    elif restart_optics_assignment:
        response = HttpResponseRedirect(reverse("nice_lite:view_stream_optics", query={"selected_job_id": jobid}))
    elif restart_generate_pickrefs:
        response = HttpResponseRedirect(reverse("nice_lite:view_stream_generate_pickrefs", query={"selected_job_id": jobid}))
    return response


@login_required(login_url="/login/")
def view_stream_delete_stream(request):
    """Deletes stream and refreshes page by redirecting to workspace view."""
    jobid = get_job_id(request)
    streamjob = StreamJob(jobid)
    streamjob.delete()
    response = redirect("nice_lite:workspace")
    return response


@login_required(login_url="/login/")
def view_stream(request, jobid):
    """Returns stream view."""
    template = "nice_stream/streamview.html"
    jobmodel = None

    if jobid is not None:
        streamjob = StreamJob(id=jobid)
        jobmodel  = streamjob.get_jobmodel()
        
    context = {
        "jobid" : jobmodel.id,
        "disp"  : jobmodel.disp,
        "proj"  : jobmodel.dset.proj.name,
        "dset"  : jobmodel.dset.name,
        "args"  : jobmodel.args
    }
    response = render(request, template, context)
    for cookie in  request.COOKIES:
        if "checksum" in cookie:
            response.delete_cookie(key=cookie)
    return response


@login_required(login_url="/login/")
def view_stream_movies(request):
    """Returns movies panel in stream view."""
    template        = "nice_stream/panelmovies.html"
    checksum_cookie = "panel_movies_checksum"
    jobmodel = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobstats" : jobmodel.preprocessing_stats,
        "args"     : jobmodel.args,
    }
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_preprocess(request):
    """Returns preprocess panel in stream view."""
    template        = "nice_stream/panelpreprocess.html"
    checksum_cookie = "panel_preprocessing_checksum"
    jobmodel = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : jobmodel.preprocessing_stats,
        "status"    : jobmodel.preprocessing_status,
    }

    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_preprocess_zoom(request):
    """Returns preprocess zoom panel in stream view."""
    template       = "nice_stream/zoompreprocess.html"
    checksum_cookie = "panel_preprocessing_checksum"
    logfile         = "preprocessing.log"
    errfile         = "preprocessing.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : jobmodel.preprocessing_stats,
        "status"    : jobmodel.preprocessing_status,
        "log"       : [],
        "error"     : ""
    }
        
    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            context["log"].append({"text":str(logtext, errors='replace')})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')

    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_optics(request):
    """Returns optics panel in stream view."""
    template        = "nice_stream/paneloptics.html"
    checksum_cookie = "panel_optics_checksum"
    jobmodel = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : jobmodel.optics_assignment_stats,
        "status"    : jobmodel.optics_assignment_status,     
    }

    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_optics_zoom(request):
    """Returns optics zoom panel in stream view."""
    template         = "nice_stream/zoomoptics.html"
    checksum_cookie = "panel_optics_checksum"
    logfile         = "optics_assignment.log"
    errfile         = "optics_assignment.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)
    
    if jobmodel is None:
        return HttpResponseNoContent()
    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : jobmodel.optics_assignment_stats,
        "status"    : jobmodel.optics_assignment_status,     
        "log"       : [],
        "error"     : ""
    }

    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            context["log"].append({"text":str(logtext, errors='replace')})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')

    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_initial_pick(request):
    """Returns initial picking panel in stream view."""
    template        = "nice_stream/panelinitialpick.html"
    checksum_cookie = "panel_initialpick_checksum"
    jobmodel = _get_jobmodel_from_request(request)
   
    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobid"        : jobmodel.id,
        "displayid"    : jobmodel.disp,
        "jobstats"     : jobmodel.initial_picking_stats,
        "status"       : jobmodel.initial_picking_status
    }

    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_initial_pick_zoom(request):
    """Returns initial picking panel in stream view."""
    template        = "nice_stream/zoominitialpick.html"
    checksum_cookie = "panel_initialpick_checksum"
    logfile         = "opening_2D.log"
    errfile         = "opening_2D.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()

    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : jobmodel.initial_picking_stats,
        "status"    : jobmodel.initial_picking_status,
        "log"       : [],
        "error"     : ""
    }

    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            context["log"].append({"text":str(logtext, errors='replace')})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')

    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_generate_pickrefs(request):
    """Returns reference generation panel in stream view."""
    template        = "nice_stream/panelgeneratepickrefs.html"
    checksum_cookie = "panel_pickrefs_checksum"
    jobmodel = _get_jobmodel_from_request(request)
    
    if jobmodel is None:
        return HttpResponseNoContent()
    
    stats = jobmodel.generate_pickrefs_stats

    if "latest_cls2D" in stats:
        stats["latest_cls2D"].sort(key=lambda d: d['res'])

    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : stats,
        "status"    : jobmodel.generate_pickrefs_status,    
    }
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_generate_pickrefs_zoom(request):
    """Returns reference generation panel in stream view."""
    template        = "nice_stream/zoomgeneratepickrefs.html"
    checksum_cookie = "panel_pickrefs_checksum"
    logfile         = "opening_2D.log"
    errfile         = "opening_2D.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)
    
    if jobmodel is None:
        return HttpResponseNoContent()
    
    stats = jobmodel.generate_pickrefs_stats

    if "latest_cls2D" in stats:
        stats["latest_cls2D"].sort(key=lambda d: d['res'])

    context = {
        "jobid"     : jobmodel.id,
        "displayid" : jobmodel.disp,
        "jobstats"  : stats,
        "status"    : jobmodel.generate_pickrefs_status,    
        "log"       : [],
        "error"     : ""
    }

    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            logtext_str = str(logtext, errors='replace')
            logpart_str = ""
            for line in logtext_str.splitlines():
                if ">>> JPEG " in line:
                    context["log"].append({"text":logpart_str})
                    split_line = line.split()
                    context["log"].append({"image":split_line[2]})
                    logpart_str = ""
                else:
                    logpart_str += line + '\n'
            if logpart_str != "":
                context["log"].append({"text":logpart_str})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')

    return _render_if_changed(request, template, context, checksum_cookie)

@login_required(login_url="/login/")
def view_stream_reference_picking(request):
    """Returns reference picking panel in stream view."""
    template        = "nice_stream/panelreferencepicking.html"
    checksum_cookie = "panel_refpick_checksum"
    jobmodel        = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobid"    : jobmodel.id,
        "jobstats" : jobmodel.reference_picking_stats,
        "status"   : jobmodel.reference_picking_status,         
    }
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_reference_picking_zoom(request):
    """Returns reference picking panel in stream view."""
    template         = "nice_stream/zoomreferencepicking.html"
    checksum_cookie  = "panel_refpick_checksum"
    logfile          = "reference_based_picking.log"
    errfile          = "reference_based_picking.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    context = {
        "jobid"    : jobmodel.id,
        "jobstats" : jobmodel.reference_picking_stats,
        "status"   : jobmodel.reference_picking_status, 
        "log"       : [],
        "error"     : ""      
    }
 
    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            context["log"].append({"text":str(logtext, errors='replace')})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_sieve_particles(request):
    """Returns particle sieving panel in stream view."""
    template        = "nice_stream/panelsieveparticles.html"
    checksum_cookie = "panel_sieve_checksum"
    jobmodel        = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    # sort cls2D on res, and annotate ref entries with matching latest entry
    stats = jobmodel.particle_sieving_stats
    if "ref_cls2D" in stats:
        stats["ref_cls2D"].sort(key=lambda d: d['res'])
        if "latest_cls2D" in stats:
            latest_by_idx = {d["idx"]: d for d in stats["latest_cls2D"]}
            for ref in stats["ref_cls2D"]:
                if ref["idx"] in latest_by_idx:
                    ref["latest"] = latest_by_idx[ref["idx"]]
        if "ref_selection" in jobmodel.master_update:
            ref_selection = set(jobmodel.master_update["ref_selection"])
            for ref in stats["ref_cls2D"]:
                ref["selected"] = ref["idx"] in ref_selection

    context = {
        "jobid"    : jobmodel.id,
        "jobstats" : jobmodel.particle_sieving_stats,
        "status"   : jobmodel.particle_sieving_status, 
    }
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_sieve_particles_zoom(request):
    """Returns particle sieving panel in stream view."""
    template         = "nice_stream/zoomsieveparticles.html"
    checksum_cookie  = "panel_sieve_checksum"
    logfile          = "particle_sieving.log"
    errfile          = "particle_sieving.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()

    # sort cls2D on res, and annotate ref entries with matching latest entry
    stats = jobmodel.particle_sieving_stats

    if "ref_cls2D" in stats:
        stats["ref_cls2D"].sort(key=lambda d: d['res'])
        if "latest_cls2D" in stats:
            latest_by_idx = {d["idx"]: d for d in stats["latest_cls2D"]}
            for ref in stats["ref_cls2D"]:
                if ref["idx"] in latest_by_idx:
                    ref["latest"] = latest_by_idx[ref["idx"]]
        if "ref_selection" in jobmodel.master_update:
            ref_selection = set(jobmodel.master_update["ref_selection"])
            for ref in stats["ref_cls2D"]:
                ref["selected"] = ref["idx"] in ref_selection

    context = {
        "jobid"    : jobmodel.id,
        "jobstats" : jobmodel.particle_sieving_stats,
        "status"   : jobmodel.particle_sieving_status, 
        "log"      : [],
        "error"    : "" 
    }

    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            logtext_str = str(logtext, errors='replace')
            logpart_str = ""
            for line in logtext_str.splitlines():
                if ">>> JPEG " in line:
                    context["log"].append({"text":logpart_str})
                    split_line = line.split()
                    context["log"].append({"image":split_line[2]})
                    logpart_str = ""
                else:
                    logpart_str += line + '\n'
            if logpart_str != "":
                context["log"].append({"text":logpart_str})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_classification_2D(request):
    """Returns 2D classification panel in stream view."""
    template        = "nice_stream/panelclassification2D.html"
    checksum_cookie = "panel_cls2D_checksum"  
    jobmodel        = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    # Remove empty classes (pop == 0), then sort by resolution.
    if "latest_cls2D" in jobmodel.classification_2D_stats:
        latest_cls2d = []
        for cls2d in jobmodel.classification_2D_stats["latest_cls2D"]:
            try:
                if float(cls2d.get("pop", 0)) == 0.0:
                    continue
            except (TypeError, ValueError, AttributeError):
                pass
            latest_cls2d.append(cls2d)
        jobmodel.classification_2D_stats["latest_cls2D"] = sorted(latest_cls2d, key=lambda d: d['res'], reverse=False)
        
    context = {
        "jobid"    : jobmodel.id,
        "jobstats" : jobmodel.classification_2D_stats,
        "status"   : jobmodel.classification_2D_status,  
    }
    return _render_if_changed(request, template, context, checksum_cookie)
    


@login_required(login_url="/login/")
def view_stream_classification_2D_zoom(request):
    """Returns 2D classification zoom panel in stream view."""
    template        = "nice_stream/zoomclassification2D.html"
    checksum_cookie = "panel_cls2D_checksum"  
    logfile         = "classification_2D.log"
    errfile         = "classification_2D.error"
    jobmodel, jobdir = _get_jobmodel_and_dir_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    # Remove empty classes (pop == 0), then sort by resolution.
    if "latest_cls2D" in jobmodel.classification_2D_stats:
        latest_cls2d = []
        for cls2d in jobmodel.classification_2D_stats["latest_cls2D"]:
            try:
                if float(cls2d.get("pop", 0)) == 0.0:
                    continue
            except (TypeError, ValueError, AttributeError):
                pass
            latest_cls2d.append(cls2d)
        jobmodel.classification_2D_stats["latest_cls2D"] = sorted(latest_cls2d, key=lambda d: d['res'], reverse=False)
        
    context = {
        "jobid"    : jobmodel.id,
        "jobstats" : jobmodel.classification_2D_stats,
        "status"   : jobmodel.classification_2D_status,  
        "log"      : [],
        "error"    : "" 
    }

    logfile = os.path.join(jobdir, logfile)
    errfile = os.path.join(jobdir, errfile)
    if os.path.exists(logfile) and os.path.isfile(logfile):
        with open(logfile, 'rb') as f:
            logtext = f.read()
            logtext_str = str(logtext, errors='replace')
            logpart_str = ""
            for line in logtext_str.splitlines():
                if ">>> JPEG " in line:
                    context["log"].append({"text":logpart_str})
                    split_line = line.split()
                    context["log"].append({"image":split_line[2]})
                    logpart_str = ""
                else:
                    logpart_str += line + '\n'
            if logpart_str != "":
                context["log"].append({"text":logpart_str})
    if os.path.exists(errfile) and os.path.isfile(errfile):
        with open(errfile, 'rb') as f:
            errortext = f.read()
            context["error"] = str(errortext, errors='replace')
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_particle_sets(request):
    """Returns particle sets panel in stream view."""
    template        = "nice_stream/panelparticlesets.html"
    checksum_cookie = "panel_particlesets_checksum" 
    jobmodel        = _get_jobmodel_from_request(request)

    if jobmodel is None:
        return HttpResponseNoContent()
    
    projectid       = jobmodel.dset.proj.id
    workspacemodels = WorkspaceModel.objects.filter(proj=projectid)
    workspaces = []
    for workspacemodel in workspacemodels:
        workspaces.append({
            "id"   : workspacemodel.id, 
            "name" : workspacemodel.name
        })

    context = {
        "jobid"      : jobmodel.id,
        "jobstats"   : jobmodel.particle_sets_stats,
        "workspaces" : workspaces
    }
    return _render_if_changed(request, template, context, checksum_cookie)


@login_required(login_url="/login/")
def view_stream_logs(request, jobid, log, error):
    """Fallback logs endpoint that routes back to the stream view."""
    return redirect("nice_lite:view_stream", jobid=jobid)


@login_required(login_url="/login/")
def view_stream_update_description(request):
    """Update description for given stream job id."""
    jobid           = get_job_id(request)
    new_description = get_string(request.POST, "new_stream_description")
    if jobid is not None:
        streamjob = StreamJob(jobid)
        streamjob.set_description(new_description)
    return HttpResponseNoContent()


@login_required(login_url="/login/")
def view_stream_update_parameters(request):
    """Update stream process parameters."""
    jobid = get_job_id(request)
    if jobid is not None:
        streamjob   = StreamJob(jobid)
        ctfres      = get_float(request.POST, "ctfres", silent=True)
        astigmatism = get_float(request.POST, "astigmatism", silent=True)
        icescore    = get_float(request.POST, "icescore", silent=True)
        streamjob.update(ctfres, astigmatism, icescore)
    return HttpResponse(status=204)


@login_required(login_url="/login/")
def view_stream_link_particle_set(request, jobid, setid, filename, type):
    """Tie stream output particle sets to classic jobs."""
    link_workspace_id = get_integer(request.POST, "link_workspace_id")
    streamjob         = StreamJob(id=jobid)
    classicjob        = BatchJob()
    project           = Project(request=request)
    workspace         = Workspace(request=request)
    link_workspace    = Workspace(link_workspace_id) if link_workspace_id is not None else workspace

    if type == "snapshot":
        set_proj = os.path.join(
            project.dirc,
            workspace.dirc,
            streamjob.dirc,
            "classification_2D",
            "snapshots",
            pathlib.Path(filename).stem,
            filename,
        )
        classicjob.linkParticleSet(project, link_workspace, set_proj)
    elif type == "final":
        set_proj  = os.path.join(project.dirc, workspace.dirc, streamjob.dirc, "classification_2D", "stream_abinitio2D.simple")
        set_desel = os.path.join(project.dirc, workspace.dirc, streamjob.dirc, filename)
        classicjob.linkParticleSetFinal(project, link_workspace, set_proj, set_desel)

    classicjob.set_description("from " + workspace.name + "->" + str(streamjob.id) + " stream->particle set " + str(setid))
    response = redirect("nice_lite:stream")
    response.set_cookie(key="selected_project_id", value=project.id)
    response.set_cookie(key="selected_workspace_id", value=link_workspace.id)
    return response


@login_required(login_url="/login/")
def view_stream_select_pickrefs(request):
    """Store the picking-reference selection and redirect to the stream view."""
    jobid         = get_job_id(request)
    streamjob     = StreamJob(jobid)
    raw_selection = request.POST.get("final_selection", "")
    if not raw_selection:
        print_error("select_pickrefs: final_selection missing")
        return redirect("nice_lite:view_stream", jobid=jobid)
    final_selection = raw_selection
    if not streamjob.select_pickrefs(final_selection):
        print_error(f"select_pickrefs: failed for job {jobid}")
    return redirect("nice_lite:view_stream", jobid=jobid)


@login_required(login_url="/login/")
def view_stream_update_classification_2D_mskdiam(request):
    """Validate and store mask diameter for 2D classification."""
    jobid     = get_job_id(request)
    streamjob = StreamJob(jobid)
    mskdiam   = get_float(request.POST, "mskdiam")
    if mskdiam is None:
        print_error("update_classification_2D_mskdiam: mskdiam missing or non-numeric")
    elif mskdiam <= 0:
        print_error(f"update_classification_2D_mskdiam: invalid mskdiam {mskdiam} (must be > 0)")
    else:
        if not streamjob.update_mskdiam(mskdiam):
            print_error(f"update_classification_2D_mskdiam: update_mskdiam failed for job {jobid}")
    streamviewclassification2D = StreamViewClassification2D(request, None, jobid)
    return streamviewclassification2D.render()


@login_required(login_url="/login/")
def view_stream_snapshot_classification_2D(request):
    """Record a 2D-classification snapshot particle set and redirect to stream view."""
    jobid              = get_job_id(request)
    streamjob          = StreamJob(jobid)
    raw_selection      = request.POST.get("snapshot_selection", "")
    snapshot_iteration = get_integer(request.POST, "snapshot_iteration")
    if not raw_selection:
        print_error("snapshot_stream_classification_2D: snapshot_selection missing")
        return redirect("nice_lite:view_stream", jobid=jobid)
    if snapshot_iteration is None:
        print_error("snapshot_stream_classification_2D: snapshot_iteration missing or non-integer")
        return redirect("nice_lite:view_stream", jobid=jobid)
    snapshot_selection = [int(s) for s in raw_selection.split(",")]
    if not streamjob.snapshot_classification_2D(snapshot_selection, snapshot_iteration):
        print_error(f"snapshot_stream_classification_2D: failed for job {jobid}")
    return redirect("nice_lite:view_stream", jobid=jobid)


@login_required(login_url="/login/")
def view_stream_select_classification_2D(request):
    """Record final 2D-classification particle selection and redirect to stream view."""
    jobid                 = get_job_id(request)
    streamjob             = StreamJob(jobid)
    raw_deselection       = request.POST.get("final_deselection", "")
    final_selection_ptcls = get_integer(request.POST, "final_selection_ptcls")
    if not raw_deselection:
        print_error("select_stream_classification_2D: final_deselection missing")
        return redirect("nice_lite:view_stream", jobid=jobid)
    if final_selection_ptcls is None:
        print_error("select_stream_classification_2D: final_selection_ptcls missing or non-integer")
        return redirect("nice_lite:view_stream", jobid=jobid)
    final_deselection = [int(s) for s in raw_deselection.split(",")]
    if not streamjob.selection_classification_2D(final_deselection, final_selection_ptcls):
        print_error(f"select_stream_classification_2D: failed for job {jobid}")
    return redirect("nice_lite:view_stream", jobid=jobid)



