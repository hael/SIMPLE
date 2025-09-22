import os
from django.shortcuts import render

from ..data_structures.job import Job
from ..data_structures.project import Project
from ..data_structures.simple  import SIMPLEStream

class StreamView:

    template = "nice_lite/streamview.html"
   # simplestream = None

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        context = {
            "jobid" : self.job.id,
            "proj"  : self.job.dset.proj.name,
            "dset"  : self.job.dset.name,
            "args"  : self.job.args
        }
        response = render(self.request, self.template, context)
        return response
    
class StreamViewLogs:

    template = "nice_lite/streamviewlogs.html"
   # simplestream = None

    def __init__(self, request, jobid, log, error):
        self.request = request
        self.job     = Job(id=jobid)
        jobdir = self.job.getAbsDir()
        self.log   = os.path.join(jobdir, log)
        self.error = os.path.join(jobdir, error)

    def render(self):
        context = {
            "jobid" : self.job.id
        }
        if os.path.exists(self.log) and os.path.isfile(self.log):
            with open(self.log, 'rb') as f:
                logtext = f.read()
                context["log"] = str(logtext, errors='replace')
        if os.path.exists(self.error) and os.path.isfile(self.error):
            with open(self.error, 'rb') as f:
                errortext = f.read()
                context["error"] = str(errortext, errors='replace')
        response = render(self.request, self.template, context)
        return response
    
class StreamViewMovies:

    template     = "nice_stream/panelmovies.html"

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
        elif jobidzoom is not None:
            self.job = Job(id=jobidzoom)
            self.template = self.templatezoom
        
    def render(self):
        context = {
            "jobstats" : self.job.preprocessing_stats,
            "args"     : self.job.args,
        }
        response = render(self.request, self.template, context)
        return response   
    
class StreamViewPreprocess:

    template     = "nice_stream/panelpreprocess.html"
    templatezoom = "nice_stream/zoompreprocess.html"
    logfile  = "preprocessing.log"
    errfile  = "preprocessing.error"

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
        elif jobidzoom is not None:
            self.job = Job(id=jobidzoom)
            self.template = self.templatezoom
        
    def render(self):
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.preprocessing_stats,
            "status"   : self.job.preprocessing_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile    
        }
        response = render(self.request, self.template, context)
        return response
    
class StreamViewOptics:

    template = "nice_stream/paneloptics.html"
    logfile  = "optics_assignment.log"
    errfile  = "optics_assignment.error"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.optics_assignment_stats,
            "status"   : self.job.optics_assignment_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile      
        }
        response = render(self.request, self.template, context)
        return response

class StreamViewInitialPick:

    template = "nice_stream/panelinitialpick.html"
    logfile  = "initial_picking.log"
    errfile  = "initial_picking.error"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.initial_picking_stats,
            "status"   : self.job.initial_picking_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile      
        }
        response = render(self.request, self.template, context)
        return response

class StreamViewGeneratePickrefs:

    template = "nice_stream/panelgeneratepickrefs.html"
    logfile  = "generate_picking_refs.log"
    errfile  = "generate_picking_refs.error"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.generate_pickrefs_stats,
            "status"   : self.job.generate_pickrefs_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile       
        }
        response = render(self.request, self.template, context)
        return response

class StreamViewReferencePicking:

    template = "nice_stream/panelreferencepicking.html"
    logfile  = "reference_based_picking.log"
    errfile  = "reference_based_picking.error"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.reference_picking_stats,
            "status"   : self.job.reference_picking_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile            
        }
        response = render(self.request, self.template, context)
        return response

class StreamViewSieveParticles:

    template = "nice_stream/panelsieveparticles.html"    
    logfile  = "particle_sieving.log"
    errfile  = "particle_sieving.error"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        # sort cls2D on pop
        if "accepted_cls2D" in self.job.particle_sieving_stats:
            self.job.particle_sieving_stats["accepted_cls2D"] = sorted(self.job.particle_sieving_stats["accepted_cls2D"], key=lambda d: d['pop'], reverse=True)
        if "rejected_cls2D" in self.job.particle_sieving_stats:
            self.job.particle_sieving_stats["rejected_cls2D"] = sorted(self.job.particle_sieving_stats["rejected_cls2D"], key=lambda d: d['pop'], reverse=True)
        if "latest_accepted_cls2D" in self.job.particle_sieving_stats:
            self.job.particle_sieving_stats["latest_accepted_cls2D"] = sorted(self.job.particle_sieving_stats["latest_accepted_cls2D"], key=lambda d: d['pop'], reverse=True)
        if "latest_rejected_cls2D" in self.job.particle_sieving_stats:
            self.job.particle_sieving_stats["latest_rejected_cls2D"] = sorted(self.job.particle_sieving_stats["latest_rejected_cls2D"], key=lambda d: d['pop'], reverse=True)   
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.particle_sieving_stats,
            "status"   : self.job.particle_sieving_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile   
        }
        response = render(self.request, self.template, context)
        return response
    
class StreamViewClassification2D:

    template = "nice_stream/panelclassification2D.html"    
    logfile  = "classification_2D.log"
    errfile  = "classification_2D.error"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        # sort cls2D on pop
        if "latest_cls2D" in self.job.classification_2D_stats:
            self.job.classification_2D_stats["latest_cls2D"] = sorted(self.job.classification_2D_stats["latest_cls2D"], key=lambda d: d['pop'], reverse=True)
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.classification_2D_stats,
            "status"   : self.job.classification_2D_status,
            "logfile"  : self.logfile,
            "errfile"  : self.errfile   
        }
        response = render(self.request, self.template, context)
        return response
    
class StreamViewParticleSets:

    template = "nice_stream/panelparticlesets.html"    

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        project = Project(project_id=self.job.dset.proj.id)
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.particle_sets_stats,
            "workspaces" : project.workspaces_list
        }
        response = render(self.request, self.template, context)
        return response