from django.shortcuts import render


from ..data_structures.project import Project
from ..data_structures.simple  import SIMPLEStream

class NewStreamView:

    template = "nice_lite/newstream.html"
    simplestream = None

    def __init__(self, request):
        self.request = request
        self.simplestream = SIMPLEStream()
        
    def render(self):
        context = {}
        if self.simplestream.ui is not None:
            if "user_inputs" in self.simplestream.ui:
                context["user_inputs"] = self.simplestream.ui["user_inputs"]
        response = render(self.request, self.template, context)
        return response