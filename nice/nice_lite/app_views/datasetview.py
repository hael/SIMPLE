import os
import json
import hashlib
from django.http import HttpResponse
from django.shortcuts import render

from ..models import JobModel
from ..data_structures.project import Project
from ..data_structures.dataset import Dataset
from ..data_structures.job     import Job

class DatasetView:

    template        = "nice_stream/dataset.html"
    checksum_cookie = "dataset_checksum"

    def __init__(self, request, project=None, dataset=None):
        if project is not None:
            self.project = project
        else:
            self.project = Project(request=request)
        if dataset is not None:
            self.dataset = dataset
        else:
            self.dataset = Dataset(request=request)
        self.request = request
    
    def render(self):
        # check dataset is in project, zero if not
        if not self.project.containsDataset(self.dataset.id):
            self.dataset.id = 0

        # get jobs
        jobs = JobModel.objects.filter(dset=self.dataset.id)
        jobstats = "" # forces hash change on status change
        for jobmodel in jobs:
            # update status
            job = Job(id=jobmodel.id)
            jobstats = jobstats + job.status

        context = {"current_project_id"   : self.project.id,
                   "current_dataset_id"   : self.dataset.id,
                   "current_project_name" : self.project.name,
                   "current_dataset_name" : self.dataset.name,
                   "created"              : self.dataset.cdat, 
                   "modified"             : self.dataset.mdat,
                   "user"                 : self.dataset.user,
                   "folder"               : os.path.join(self.project.dirc, self.dataset.link),
                   "description"          : self.dataset.desc,
                   "jobstats"             : jobstats
                  }
        hash = hashlib.md5(json.dumps(context, sort_keys=True, default=str).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        response = HttpResponse(status=204)
        if(old_checksum == "none" or old_checksum != checksum):
            jobs = JobModel.objects.filter(dset=self.dataset.id)
            for jobmodel in jobs:
                # sort cls2D
                if "latest_cls2D" in jobmodel.classification_2D_stats:
                    jobmodel.classification_2D_stats["latest_cls2D"] = sorted(jobmodel.classification_2D_stats["latest_cls2D"], key=lambda d: d['pop'], reverse=True)
            context["jobs"] = jobs
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
        response.set_cookie(key='selected_project_id', value=self.project.id)
        response.set_cookie(key='selected_dataset_id', value=self.dataset.id)
        return response