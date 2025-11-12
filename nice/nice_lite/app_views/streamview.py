import os
import json
import hashlib
from django.shortcuts import render
from django.http import HttpResponse

from ..data_structures.job import Job
from ..data_structures.project import Project
from ..data_structures.simple  import SIMPLEStream

class StreamView:

    template = "nice_stream/streamview.html"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(id=jobid)
        
    def render(self):
        context = {
            "jobid" : self.job.id,
            "disp"  : self.job.disp,
            "proj"  : self.job.dset.proj.name,
            "dset"  : self.job.dset.name,
            "args"  : self.job.args
        }
        response = render(self.request, self.template, context)
        for cookie in  self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response
    
class StreamViewMovies:

    template        = "nice_stream/panelmovies.html"
    checksum_cookie = "panel_movies_checksum"

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
        hash = hashlib.md5(json.dumps(context, sort_keys=True).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
            return response 
        else:
            return HttpResponse(status=204)
    
class StreamViewPreprocess:

    template        = "nice_stream/panelpreprocess.html"
    templatezoom    = "nice_stream/zoompreprocess.html"
    checksum_cookie = "panel_preprocessing_checksum"
    logfile         = "preprocessing.log"
    errfile         = "preprocessing.error"
    zoom            = False
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True

    def render(self):
        context = {
            "jobid"     : self.job.id,
            "displayid" : self.job.disp,
            "jobstats"  : self.job.preprocessing_stats,
            "status"    : self.job.preprocessing_status,
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
            return HttpResponse(status=204)
    
class StreamViewOptics:

    template        = "nice_stream/paneloptics.html"
    templatezoom    = "nice_stream/zoomoptics.html"
    checksum_cookie = "panel_optics_checksum"
    logfile         = "optics_assignment.log"
    errfile         = "optics_assignment.error"
    zoom            = False
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True
        
    def render(self):
        context = {
            "jobid"     : self.job.id,
            "displayid" : self.job.disp,
            "jobstats"  : self.job.optics_assignment_stats,
            "status"    : self.job.optics_assignment_status,     
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
            return HttpResponse(status=204)

class StreamViewInitialPick:

    template        = "nice_stream/panelinitialpick.html"
    templatezoom    = "nice_stream/zoominitialpick.html"
    checksum_cookie = "panel_initialpick_checksum"
    logfile         = "opening_2D.log"
    errfile         = "opening_2D.error"
    zoom            = False
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True
        
    def render(self):
        start_pc = 0
        end_pc   = 0
        if "picking_diameters" in self.job.initial_picking_stats:
            self.job.initial_picking_stats["picking_diameters"] = sorted(list(set(self.job.initial_picking_stats["picking_diameters"])))
        if "refinement_diameters" in self.job.initial_picking_stats:
            self.job.initial_picking_stats["refinement_diameters"] = sorted(list(set(self.job.initial_picking_stats["refinement_diameters"])))   
        if "picking_diameters" in self.job.initial_picking_stats and "refinement_diameters" in self.job.initial_picking_stats:
            if len(self.job.initial_picking_stats["picking_diameters"]) > 1 and len(self.job.initial_picking_stats["refinement_diameters"]) > 1:
                start_pc = self.job.initial_picking_stats["picking_diameters"].index(self.job.initial_picking_stats["refinement_diameters"][0])  * (100 / (len(self.job.initial_picking_stats["picking_diameters"]) - 1))
                end_pc   = self.job.initial_picking_stats["picking_diameters"].index(self.job.initial_picking_stats["refinement_diameters"][-1]) * (100 / (len(self.job.initial_picking_stats["picking_diameters"]) - 1))
        context = {
            "jobid"        : self.job.id,
            "displayid"    : self.job.disp,
            "jobstats"     : self.job.initial_picking_stats,
            "status"       : self.job.initial_picking_status,
            "refine_start" : start_pc,
            "refine_width" : end_pc - start_pc,  
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
            return HttpResponse(status=204)

class StreamViewGeneratePickrefs:

    template        = "nice_stream/panelgeneratepickrefs.html"
    templatezoom    = "nice_stream/zoomgeneratepickrefs.html"
    checksum_cookie = "panel_pickrefs_checksum"
    logfile         = "opening_2D.log"
    errfile         = "opening_2D.error"
    zoom            = False
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True
        
    def render(self):
        context = {
            "jobid"     : self.job.id,
            "displayid" : self.job.disp,
            "jobstats"  : self.job.generate_pickrefs_stats,
            "status"    : self.job.generate_pickrefs_status,    
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
            return HttpResponse(status=204)

class StreamViewReferencePicking:

    template        = "nice_stream/panelreferencepicking.html"
    templatezoom    = "nice_stream/zoomreferencepicking.html"
    checksum_cookie = "panel_refpick_checksum"
    logfile         = "reference_based_picking.log"
    errfile         = "reference_based_picking.error"
    zoom            = False
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True
        
    def render(self):
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.reference_picking_stats,
            "status"   : self.job.reference_picking_status,         
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
            return HttpResponse(status=204)

class StreamViewSieveParticles:

    template        = "nice_stream/panelsieveparticles.html"
    templatezoom    = "nice_stream/zoomsieveparticles.html"
    checksum_cookie = "panel_sieve_checksum"
    logfile         = "particle_sieving.log"
    errfile         = "particle_sieving.error"
    zoom            = False
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True
        
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
    jobdir          = ""

    def __init__(self, request, jobid, jobidzoom):
        self.request = request
        if jobid is not None:
            self.job = Job(id=jobid)
            
        elif jobidzoom is not None:
            self.job             = Job(id=jobidzoom)
            self.template        = self.templatezoom
            self.jobdir          = self.job.getAbsDir()
            self.zoom            = True
        
    def render(self):
        # sort cls2D on pop
        if "latest_cls2D" in self.job.classification_2D_stats:
            self.job.classification_2D_stats["latest_cls2D"] = sorted(self.job.classification_2D_stats["latest_cls2D"], key=lambda d: d['pop'], reverse=True)
        context = {
            "jobid"    : self.job.id,
            "jobstats" : self.job.classification_2D_stats,
            "status"   : self.job.classification_2D_status,  
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