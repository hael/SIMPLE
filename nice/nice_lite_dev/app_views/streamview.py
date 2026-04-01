import os
import json
import hashlib
from django.shortcuts import render
from django.http import HttpResponse

from ..helpers                   import *
from ..data_structures.streamjob import StreamJob
from ..data_structures.project   import Project
from ..data_structures.simple    import SIMPLEStream

class StreamView:

    template = "nice_stream/streamview.html"
    jobmodel = None

    def __init__(self, request, jobid):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        
    def render(self):
        context = {
            "jobid" : self.jobmodel.id,
            "disp"  : self.jobmodel.disp,
            "proj"  : self.jobmodel.dset.proj.name,
            "dset"  : self.jobmodel.dset.name,
            "args"  : self.jobmodel.args
        }
        response = render(self.request, self.template, context)
        for cookie in  self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response
    
class StreamViewMovies:

    template        = "nice_stream/panelmovies.html"
    checksum_cookie = "panel_movies_checksum"
    jobmodel        = None

    def __init__(self, request, jobid):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        
    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        context = {
            "jobstats" : self.jobmodel.preprocessing_stats,
            "args"     : self.jobmodel.args,
        }
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponseNoContent()
    
class StreamViewPreprocess:

    template        = "nice_stream/panelpreprocess.html"
    templatezoom    = "nice_stream/zoompreprocess.html"
    checksum_cookie = "panel_preprocessing_checksum"
    logfile         = "preprocessing.log"
    errfile         = "preprocessing.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True

    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        context = {
            "jobid"     : self.jobmodel.id,
            "displayid" : self.jobmodel.disp,
            "jobstats"  : self.jobmodel.preprocessing_stats,
            "status"    : self.jobmodel.preprocessing_status,
        }
        if self.zoom:
            context["log"]   = ""
            context["error"] = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
            if os.path.exists(logfile) and os.path.isfile(logfile):
                with open(logfile, 'rb') as f:
                    logtext = f.read()
                    context["log"] = str(logtext, errors='replace')
            if os.path.exists(errfile) and os.path.isfile(errfile):
                with open(errfile, 'rb') as f:
                    errortext = f.read()
                    context["error"] = str(errortext, errors='replace')
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponseNoContent()
    
class StreamViewOptics:

    template        = "nice_stream/paneloptics.html"
    templatezoom    = "nice_stream/zoomoptics.html"
    checksum_cookie = "panel_optics_checksum"
    logfile         = "optics_assignment.log"
    errfile         = "optics_assignment.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True
        
    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        context = {
            "jobid"     : self.jobmodel.id,
            "displayid" : self.jobmodel.disp,
            "jobstats"  : self.jobmodel.optics_assignment_stats,
            "status"    : self.jobmodel.optics_assignment_status,     
        }
        if self.zoom:
            context["log"]   = ""
            context["error"] = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
            if os.path.exists(logfile) and os.path.isfile(logfile):
                with open(logfile, 'rb') as f:
                    logtext = f.read()
                    context["log"] = str(logtext, errors='replace')
            if os.path.exists(errfile) and os.path.isfile(errfile):
                with open(errfile, 'rb') as f:
                    errortext = f.read()
                    context["error"] = str(errortext, errors='replace')
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponseNoContent()

class StreamViewInitialPick:

    template        = "nice_stream/panelinitialpick.html"
    templatezoom    = "nice_stream/zoominitialpick.html"
    checksum_cookie = "panel_initialpick_checksum"
    logfile         = "opening_2D.log"
    errfile         = "opening_2D.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True

    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        context = {
            "jobid"        : self.jobmodel.id,
            "displayid"    : self.jobmodel.disp,
            "jobstats"     : self.jobmodel.initial_picking_stats,
            "status"       : self.jobmodel.initial_picking_status
        }
        if self.zoom:
            context["log"]   = ""
            context["error"] = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
            if os.path.exists(logfile) and os.path.isfile(logfile):
                with open(logfile, 'rb') as f:
                    logtext = f.read()
                    context["log"] = str(logtext, errors='replace')
            if os.path.exists(errfile) and os.path.isfile(errfile):
                with open(errfile, 'rb') as f:
                    errortext = f.read()
                    context["error"] = str(errortext, errors='replace')
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponseNoContent()

class StreamViewGeneratePickrefs:

    template        = "nice_stream/panelgeneratepickrefs.html"
    templatezoom    = "nice_stream/zoomgeneratepickrefs.html"
    checksum_cookie = "panel_pickrefs_checksum"
    logfile         = "opening_2D.log"
    errfile         = "opening_2D.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True
        
    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        context = {
            "jobid"     : self.jobmodel.id,
            "displayid" : self.jobmodel.disp,
            "jobstats"  : self.jobmodel.generate_pickrefs_stats,
            "status"    : self.jobmodel.generate_pickrefs_status,    
        }
        if self.zoom:
            context["log"]       = []
            context["error"]     = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
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
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponseNoContent()

class StreamViewReferencePicking:

    template        = "nice_stream/panelreferencepicking.html"
    templatezoom    = "nice_stream/zoomreferencepicking.html"
    checksum_cookie = "panel_refpick_checksum"
    logfile         = "reference_based_picking.log"
    errfile         = "reference_based_picking.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True
        
    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        context = {
            "jobid"    : self.jobmodel.id,
            "jobstats" : self.jobmodel.reference_picking_stats,
            "status"   : self.jobmodel.reference_picking_status,         
        }
        if self.zoom:
            context["log"]   = ""
            context["error"] = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
            if os.path.exists(logfile) and os.path.isfile(logfile):
                with open(logfile, 'rb') as f:
                    logtext = f.read()
                    context["log"] = str(logtext, errors='replace')
            if os.path.exists(errfile) and os.path.isfile(errfile):
                with open(errfile, 'rb') as f:
                    errortext = f.read()
                    context["error"] = str(errortext, errors='replace')
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponseNoContent()

class StreamViewSieveParticles:

    template        = "nice_stream/panelsieveparticles.html"
    templatezoom    = "nice_stream/zoomsieveparticles.html"
    checksum_cookie = "panel_sieve_checksum"
    logfile         = "particle_sieving.log"
    errfile         = "particle_sieving.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True
        
    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        # sort cls2D on res, and annotate ref entries with matching latest entry
        stats = self.jobmodel.particle_sieving_stats
        if "ref_cls2D" in stats:
            stats["ref_cls2D"].sort(key=lambda d: d['res'])
            if "latest_cls2D" in stats:
                latest_by_idx = {d["idx"]: d for d in stats["latest_cls2D"]}
                for ref in stats["ref_cls2D"]:
                    if ref["idx"] in latest_by_idx:
                        ref["latest"] = latest_by_idx[ref["idx"]]
            if "ref_selection" in self.jobmodel.master_update:
                ref_selection = set(self.jobmodel.master_update["ref_selection"])
                for ref in stats["ref_cls2D"]:
                    ref["selected"] = ref["idx"] in ref_selection
        context = {
            "jobid"    : self.jobmodel.id,
            "jobstats" : self.jobmodel.particle_sieving_stats,
            "status"   : self.jobmodel.particle_sieving_status, 
        }
        if self.zoom:
            context["log"]   = []
            context["error"] = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
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
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponse(status=204)
    
class StreamViewClassification2D:

    template        = "nice_stream/panelclassification2D.html"
    templatezoom    = "nice_stream/zoomclassification2D.html"
    checksum_cookie = "panel_cls2D_checksum"  
    logfile         = "classification_2D.log"
    errfile         = "classification_2D.error"
    zoom            = False
    jobdir          = None
    jobmodel        = None

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            streamjob     = StreamJob(id=jobid)
            self.jobmodel = streamjob.get_jobmodel()
        elif jobidzoom is not None:
            streamjob     = StreamJob(id=jobidzoom)
            self.jobmodel = streamjob.get_jobmodel()
            self.template = self.templatezoom
            self.jobdir   = streamjob.get_absdir()
            self.zoom     = True
        
    def render(self):
        if self.jobmodel is None:
            return HttpResponseNoContent()
        # sort cls2D on pop
        if "latest_cls2D" in self.jobmodel.classification_2D_stats:
            self.jobmodel.classification_2D_stats["latest_cls2D"] = sorted(self.jobmodel.classification_2D_stats["latest_cls2D"], key=lambda d: d['pop'], reverse=True)
        context = {
            "jobid"    : self.jobmodel.id,
            "jobstats" : self.jobmodel.classification_2D_stats,
            "status"   : self.jobmodel.classification_2D_status,  
        }
        if self.zoom:
            context["log"]   = []
            context["error"] = ""
            logfile = os.path.join(self.jobdir, self.logfile)
            errfile = os.path.join(self.jobdir, self.errfile)
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
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponse(status=204)
    
class StreamViewParticleSets:

    template        = "nice_stream/panelparticlesets.html"
    checksum_cookie = "panel_particlesets_checksum" 

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        project = Project(project_id=self.job.dset.proj.id)
        context = {
            "jobid"      : self.job.id,
            "jobstats"   : self.job.particle_sets_stats,
            "workspaces" : project.workspaces_list
        }
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponse(status=204)