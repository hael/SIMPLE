from django.shortcuts import render

from ..data_structures.project import Project

class NewProjectView:

    template = "nice_lite/newproject.html"

    def __init__(self, request, caller):
        self.request = request
        self.caller  = caller
        self.project = Project(request=request)
        
    def render(self):
        context = {
            "caller" : self.caller
        }
        response = render(self.request, self.template, context)
        return response