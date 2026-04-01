from django.shortcuts import render
from django.contrib.auth.forms import AuthenticationForm

from ..helpers import *
from ..models import ProjectModel, DatasetModel
from ..data_structures.project import Project
from ..data_structures.dataset import Dataset
from ..data_structures.workspace import Workspace

class CoreViewClassic:

    template = "nice_classic/classic.html"

    def __init__(self, request, project=None, workspace=None):
        if project is not None:
            self.project = project
        else:
            self.project = Project(request=request)
        if workspace is not None:
            self.workspace = workspace
        else:
            self.workspace = Workspace(request=request)
        self.request = request
    
    def render(self):
        # check workspace is in project, zero if not
        if not self.project.containsWorkspace(self.workspace.id):
            self.workspace.id = 0
        
        context = {"current_project_id"   : self.project.id,
                   "current_workspace_id" : self.workspace.id,
                   "projects"             : ProjectModel.objects.all(),
                   "workspaces"           : self.project.workspaces_list,
                  }
        
        response = render(self.request, self.template, context)
        response.set_cookie(key='selected_project_id', value=self.project.id)
        response.set_cookie(key='selected_workspace_id', value=self.workspace.id)
        response.set_cookie(key='mode',                value="classic")
        for cookie in  self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response

class LoginView:

    template = "login.html"
    
    def __init__(self, request):
        self.request = request

    def render(self):
        form     = AuthenticationForm()
        response = render(self.request, self.template, {'form':form, 'title':'log in'})
        return response