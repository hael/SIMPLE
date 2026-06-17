# django imports
from django.urls               import reverse
from django.shortcuts          import render
from django.contrib            import messages
from django.contrib.auth.forms import AuthenticationForm

# local imports
from ..helpers                   import clear_checksum_cookies
from ..models                    import ProjectModel, DatasetModel
from ..data_structures.dataset   import Dataset
from ..data_structures.project   import Project

class StreamIndex:

    template   = "nice_stream/stream.html"
    request    = None
    iframeurl  = None
    username   = None
    projects   = []
    datasets   = []
    projectid  = 0
    datasetid  = 0
    
    def __init__(self, request, projectid=None, datasetid=None, username=None):
        self.request  = request
        self.username = username
        self.projects = ProjectModel.objects.all()
        if username is None:
            messages.add_message(self.request, messages.ERROR, "username is none")
        if projectid is None:
            if len(self.projects) == 0:
                messages.add_message(self.request, messages.INFO, "please create a project")
            else:
                messages.add_message(self.request, messages.INFO, "please select a project")
        else:
            self.projectid = projectid
            self.datasets  = DatasetModel.objects.filter(proj=projectid)
        if datasetid is None:
            if projectid is not None:
                if len(self.datasets) == 0:
                    messages.add_message(self.request, messages.INFO, "please create a dataset")
                else:
                    messages.add_message(self.request, messages.INFO, "please select a dataset")
        else:
            dataset = Dataset(datasetid)
            if dataset.in_project(projectid):
                self.datasetid = datasetid
            else:
                messages.add_message(self.request, messages.INFO, "please select a dataset")
        if projectid == -1:
            self.iframeurl = reverse("nice_lite:new_project", args=["stream"])
        elif datasetid == -1 and self.username is not None:
            messages.add_message(self.request, messages.INFO, "created new dataset")
            project     = Project(self.projectid)
            new_dataset = Dataset()
            new_dataset.new(project, self.username)
            self.datasetid = new_dataset.get_id()
            self.datasets  = DatasetModel.objects.filter(proj=projectid)
            self.iframeurl = reverse("nice_lite:dataset", query={"selected_dataset_id" : self.datasetid})
        elif self.datasetid > 0:
            self.iframeurl = reverse("nice_lite:dataset", query={"selected_dataset_id" : self.datasetid})
    
    def render(self):
        context = {"current_project_id" : self.projectid,
                   "current_dataset_id" : self.datasetid,
                   "projects"           : self.projects,
                   "datasets"           : self.datasets,
                   "iframeurl"          : self.iframeurl
                  }
        response = render(self.request, self.template, context)
        response.set_cookie(key='selected_project_id', value=self.projectid)
        response.set_cookie(key='selected_dataset_id', value=self.datasetid)
        response.set_cookie(key='mode',                value="stream"      )    
        clear_checksum_cookies(self.request, response)
        return response