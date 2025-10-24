from django.shortcuts import render

from ..data_structures.project import Project

class NewProjectView:

    template = "newproject.html"

    def __init__(self, request, caller):
        self.request = request
        self.caller  = caller
        self.project = Project(request=request)
        
    def render(self):
        context = {
            "caller" : self.caller
        }
        response = render(self.request, self.template, context)
        for cookie in  self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response