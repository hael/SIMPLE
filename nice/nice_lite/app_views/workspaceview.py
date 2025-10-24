import json
import os
import hashlib

from django.http import HttpResponse
from django.shortcuts import render

from ..models import JobModel
from ..data_structures.project   import Project
from ..data_structures.workspace import Workspace

class WorkspaceView:

    template = "nice_classic/workspace.html"
    checksum_cookie = "workspace_checksum"

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
        
        # get jobs

        context = {"current_project_id"     : self.project.id,
                   "current_workspace_id"   : self.workspace.id,
                   "current_project_name"   : self.project.name,
                   "current_workspace_name" : self.workspace.name,
                   "created"                : self.workspace.cdat, 
                   "modified"               : self.workspace.mdat,
                   "user"                   : self.workspace.user,
                   "folder"                 : os.path.join(self.project.dirc, self.workspace.link),
                   "description"            : self.workspace.desc,
                   "nodes"                  : json.dumps(self.workspace.nstr)
                  }
        
        hash = hashlib.md5(json.dumps(context, sort_keys=True, default=str).encode())
        checksum = hash.hexdigest()
        old_checksum = self.request.COOKIES.get(self.checksum_cookie, 'none')
        response = HttpResponse(status=204)
        if(old_checksum == "none" or old_checksum != checksum):
            context["nodes"] =  json.dumps(self.workspace.nstrhtml)
            response = render(self.request, self.template, context)
            response.set_cookie(key=self.checksum_cookie, value=checksum)
        response.set_cookie(key='selected_project_id', value=self.project.id)
        response.set_cookie(key='selected_workspace_id', value=self.workspace.id)
        return response