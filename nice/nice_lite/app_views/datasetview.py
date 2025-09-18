import os
from django.shortcuts import render

from ..models import JobModel
from ..data_structures.project import Project
from ..data_structures.dataset import Dataset
from ..data_structures.job     import Job

class DatasetView:

    template = "nice_lite/dataset.html"

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
        for jobmodel in jobs:
            job = Job(id=jobmodel.id)

        context = {"current_project_id"   : self.project.id,
                   "current_dataset_id"   : self.dataset.id,
                   "current_project_name" : self.project.name,
                   "current_dataset_name" : self.dataset.name,
                   "created"              : self.dataset.cdat, 
                   "modified"             : self.dataset.mdat,
                   "folder"               : os.path.join(self.project.dirc, self.dataset.link),
                   "description"          : self.dataset.desc,
                   "jobs"                 : jobs,
                  }
        
        response = render(self.request, self.template, context)
        response.set_cookie(key='selected_project_id', value=self.project.id)
        response.set_cookie(key='selected_dataset_id', value=self.dataset.id)
        return response