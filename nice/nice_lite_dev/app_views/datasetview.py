# global imports
import os
import json
import hashlib

# django imports
from django.http      import HttpResponse
from django.contrib   import messages
from django.shortcuts import render

# local imports
from ..helpers                   import *
from ..models                    import JobModel
from ..data_structures.project   import Project
from ..data_structures.dataset   import Dataset
from ..data_structures.streamjob import StreamJob

class DatasetView:

    template        = "nice_stream/dataset.html"
    checksum_cookie = "dataset_checksum"
    dataset         = None
    request         = None
    projectid       = 0

    def __init__(self, request, dataset=None, projectid=None):
        self.request = request
        if dataset is not None:
            self.dataset = dataset
        if projectid is not None:
            self.projectid = projectid
    
    def render(self):
        response = HttpResponseNoContent()
        if self.dataset is None:
            print_error("dataset is None")
            messages.add_message(self.request, messages.ERROR, "dataset is None")
            return response
        if self.projectid == 0:
            print_error("projectid is zero")
            messages.add_message(self.request, messages.ERROR, "projectid is zero")
            return response
        if not self.dataset.in_project(self.projectid):
            print_error("dataset is not in project")
            messages.add_message(self.request, messages.ERROR, "dataset is not in project")
            return response
        # get jobs
        jobs = JobModel.objects.filter(dset=self.dataset.id)
        jobstats = "" # forces hash change on status change
        for jobmodel in jobs:
            # update status
            job = StreamJob(id=jobmodel.id)
            jobstats = jobstats + job.get_status()
        datasetmodel = self.dataset.get_datasetmodel()
        projectmodel = datasetmodel.proj
        context = {"current_project_id"   : self.projectid,
                   "current_dataset_id"   : self.dataset.get_id(),
                   "current_project_name" : projectmodel.name,
                   "current_dataset_name" : datasetmodel.name,
                   "created"              : datasetmodel.cdat, 
                   "modified"             : datasetmodel.mdat,
                   "user"                 : datasetmodel.user,
                   "folder"               : self.dataset.get_linkpath(),
                   "description"          : datasetmodel.desc,
                   "jobstats"             : jobstats
                  }
        hash = hashlib.md5(json.dumps(context, sort_keys=True, default=str).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        if(old_checksum == "none" or old_checksum != checksum):
            jobs = JobModel.objects.filter(dset=self.dataset.id)
            for jobmodel in jobs:
                # sort cls2D
                if "latest_cls2D" in jobmodel.classification_2D_stats:
                    jobmodel.classification_2D_stats["latest_cls2D"] = sorted(jobmodel.classification_2D_stats["latest_cls2D"], key=lambda d: d['pop'], reverse=True)
            context["jobs"] = jobs
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
        response.set_cookie(key='selected_project_id', value=self.projectid)
        response.set_cookie(key='selected_dataset_id', value=self.dataset.id)
        return response