import os
from django.shortcuts import render

from ..data_structures.job import Job
#from ..data_structures.simple  import SIMPLEStream

class JobView:

    template = "nice_classic/jobview.html"
   # simplestream = None

    def __init__(self, request, jobid):
        self.request = request
        self.job     = Job(type='stream', id=jobid)
        
    def render(self):
        context = {
            "jobid" : self.job.id,
            "proj"  : self.job.dset.proj.name,
            "dset"  : self.job.dset.name
        }
        response = render(self.request, self.template, context)
        return response
    
