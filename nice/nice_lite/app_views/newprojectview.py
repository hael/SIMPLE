# django imports
from django.shortcuts import render

# local imports
from ..helpers                 import clear_checksum_cookies

class NewProjectView:

    template = "newproject.html"
    request  = None
    mode     = ""
    
    def __init__(self, request, mode):
        self.request = request
        self.mode    = mode
        
    def render(self):
        context = {
            "mode" : self.mode
        }
        response = render(self.request, self.template, context)
        clear_checksum_cookies(self.request,response)
        return response