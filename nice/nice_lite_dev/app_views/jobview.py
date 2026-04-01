import os
import math

from django.shortcuts import render

from ..data_structures.jobclassic import JobClassic
#from ..data_structures.simple  import SIMPLEStream

class JobView:

    template = "nice_classic/jobview.html"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = JobClassic(id=jobid)
        
    def render(self):
        projstats = self.job.getProjectStats()
        context = {
            "jobid"     : self.job.id,
            "prog"      : self.job.prog,
            "disp"      : self.job.disp,
            "name"      : self.job.name,
            "status"    : self.job.status,
            "projstats" : projstats,
        }
        response = render(self.request, self.template, context)
        for cookie in self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response

class JobViewMicrographs:

    template = "nice_classic/panelmicrographs.html"
    sortkey  = "n"
    sortasc  = True

    def __init__(self, request, jobid, sort_micrographs_key, sort_micrographs_asc):
        self.request = request
        self.job     = JobClassic(id=jobid)
        if sort_micrographs_key is not None:
            self.sortkey = sort_micrographs_key
        if sort_micrographs_asc is not None:
            self.sortasc = sort_micrographs_asc

    def render(self, fromp=None, top=None, page=None):
        pagen     = 50
        projstats = self.job.getProjectStats()
        pages     = []
        if "mic" in projstats and "n" in projstats["mic"]:
            npages = math.ceil(projstats["mic"]["n"] / pagen)
            for n in range(npages):
                page_fromp = (n * pagen) + 1
                page_top   = (n + 1) * pagen
                if page_top > projstats["mic"]["n"]: page_top = projstats["mic"]["n"]
                pages.append({
                    "n"     : n + 1,
                    "fromp" : page_fromp, 
                    "top"   : page_top
                })
        if fromp is None and top is None:
            fromp = pages[page - 1]["fromp"]
            top   = pages[page - 1]["top"]
        micsstats = self.job.getProjectFieldStats('mic', fromp=fromp, top=top, sortkey=self.sortkey, sortasc=self.sortasc, boxes=True)
        if self.sortkey is not None and "data" in micsstats:
            for micstat in micsstats["data"]:
               micstat["sortval"] = round(micstat[self.sortkey], -int(math.floor(math.log10(abs(micstat[self.sortkey])))) + (2)) # 3 sig fig
        context = {
            "jobid"       : self.job.id,
            "sortkey"     : self.sortkey,
            "pages"       : pages,
            "sortasc"     : self.sortasc,
            "currentpage" : int(page),
            "projstats"   : projstats,
            "mics"        : micsstats,
        }
        response = render(self.request, self.template, context)
        return response

class JobViewMicrographsHistogram:

    template = "nice_classic/panelmicrographshistogram.html"
    sortkey  = "ctfres"

    def __init__(self, request, jobid, sort_micrographs_key):
        self.request = request
        self.job     = JobClassic(id=jobid)
        if sort_micrographs_key is not None:
            self.sortkey = sort_micrographs_key

    def render(self, fromp=None, top=None, page=None):
        projstats = self.job.getProjectStats()
        micsstats = self.job.getProjectFieldStats('mic', sortkey=self.sortkey, hist=True)
        context = {
            "jobid"       : self.job.id,
            "sortkey"     : self.sortkey,
            "projstats"   : projstats,
            "mics"        : micsstats,
        }
        response = render(self.request, self.template, context)
        return response
    
class JobViewMicrographsPlot:

    template = "nice_classic/panelmicrographsplot.html"
    sortkey  = "n"
    plotkey  = "ctfres"

    def __init__(self, request, jobid, sort_micrographs_key, plot_micrographs_key):
        self.request = request
        self.job     = JobClassic(id=jobid)
        if sort_micrographs_key is not None:
            self.sortkey = sort_micrographs_key
        if plot_micrographs_key is not None:
            self.plotkey = plot_micrographs_key

    def render(self, fromp=None, top=None, page=None):
        projstats = self.job.getProjectStats()
        micsstats = self.job.getProjectFieldStats('mic', sortkey=self.sortkey, plotkey=self.plotkey)
        context = {
            "jobid"       : self.job.id,
            "sortkey"     : self.sortkey,
            "plotkey"     : self.plotkey,
            "projstats"   : projstats,
            "mics"        : micsstats,
        }
        response = render(self.request, self.template, context)
        return response
    
class JobViewCls2D:

    template = "nice_classic/panelcls2d.html"
    sortkey  = "n"
    sortasc  = True

    def __init__(self, request, jobid, sort_micrographs_key, sort_micrographs_asc):
        self.request = request
        self.job     = JobClassic(id=jobid)
        if sort_micrographs_key is not None:
            self.sortkey = sort_micrographs_key
        if sort_micrographs_asc is not None:
            self.sortasc = sort_micrographs_asc

    def render(self, fromp=None, top=None, page=None):
        projstats = self.job.getProjectStats()
        cls2dstats = self.job.getProjectFieldStats('cls2D', fromp=None, top=None, sortkey=self.sortkey, sortasc=self.sortasc)
        if self.sortkey is not None and "data" in cls2dstats:
            for cls2dstat in cls2dstats["data"]:
                cls2dstat["sortval"] = round(cls2dstat[self.sortkey], -int(math.floor(math.log10(abs(cls2dstat[self.sortkey])))) + (2)) # 3 sig fig
                if "thumbnx" in cls2dstat:
                    spritew = cls2dstat["thumbnx"]
                else:
                    spritew = math.floor(math.sqrt(cls2dstat["thumbn"]))
                if "thumbny" in cls2dstat:
                    spriteh = cls2dstat["thumbny"]
                else:
                    spriteh = math.ceil(cls2dstat["thumbn"] / spritew) 
                cls2dstat["spritew"] = spritew * 100
                cls2dstat["spriteh"] = spriteh * 100
                x = (cls2dstat["thumbidx"] - 1) % spritew
                y = math.floor((cls2dstat["thumbidx"] - 1) / spritew)
                cls2dstat["spritex"] = x * (100 / (spritew - 1)) 
                cls2dstat["spritey"] = y * (100 / (spriteh - 1))
        context = {
            "jobid"       : self.job.id,
            "sortkey"     : self.sortkey,
            "sortasc"     : self.sortasc,
            "projstats"   : projstats,
            "cls2ds"      : cls2dstats,
        }
        response = render(self.request, self.template, context)
        return response
    
class JobViewCls2DHistogram:

    template = "nice_classic/panelcls2dhistogram.html"
    sortkey  = "pop"

    def __init__(self, request, jobid, sort_cls2d_key):
        self.request = request
        self.job     = JobClassic(id=jobid)
        if sort_cls2d_key is not None:
            self.sortkey = sort_cls2d_key

    def render(self, fromp=None, top=None, page=None):
        projstats  = self.job.getProjectStats()
        cls2dstats = self.job.getProjectFieldStats('cls2D', sortkey=self.sortkey, hist=True)
        context = {
            "jobid"       : self.job.id,
            "sortkey"     : self.sortkey,
            "projstats"   : projstats,
            "cls2ds"      : cls2dstats,
        }
        response = render(self.request, self.template, context)
        return response
    
class JobViewCls2DPlot:

    template = "nice_classic/panelcls2dplot.html"
    sortkey  = "n"
    plotkey  = "res"

    def __init__(self, request, jobid, sort_cls2d_key, plot_cls2d_key):
        self.request = request
        self.job     = JobClassic(id=jobid)
        if sort_cls2d_key is not None:
            self.sortkey = sort_cls2d_key
        if plot_cls2d_key is not None:
            self.plotkey = plot_cls2d_key

    def render(self, fromp=None, top=None, page=None):
        projstats = self.job.getProjectStats()
        cls2dstats = self.job.getProjectFieldStats('cls2D', sortkey=self.sortkey, plotkey=self.plotkey)
        context = {
            "jobid"       : self.job.id,
            "sortkey"     : self.sortkey,
            "plotkey"     : self.plotkey,
            "projstats"   : projstats,
            "cls2ds"      : cls2dstats,
        }
        response = render(self.request, self.template, context)
        return response   
    
class JobViewLogs:

    template = "nice_classic/panellogs.html"

    def __init__(self, request, jobid):
        self.request = request
        self.job     = JobClassic(id=jobid)

    def render(self):
        jobdir = self.job.getAbsDir()
        stdout = ""
        stderr = ""
        stdout_file = os.path.join(jobdir, "stdout.log")
        stderr_file = os.path.join(jobdir, "stderr.log")
        if os.path.exists(stdout_file) and os.path.isfile(stdout_file):
             with open(stdout_file, 'rb') as f:
                logtext = f.read()
                stdout  = str(logtext, errors='replace')
        if os.path.exists(stderr_file) and os.path.isfile(stderr_file):
             with open(stderr_file, 'rb') as f:
                logtext = f.read()
                stderr  = str(logtext, errors='replace')
        context = {
            "jobid"  : self.job.id,
            "stdout" : stdout,
            "stderr" : stderr
        }
        response = render(self.request, self.template, context)
        return response