from django.shortcuts import render


from ..data_structures.project import Project
from ..data_structures.simple  import SIMPLEStream

class NewStreamView:

    template = "nice_lite/newstream.html"
    simplestream = None

    def __init__(self, request, args=None):
        self.request = request
        self.simplestream = SIMPLEStream()
        self.args = args

    def render(self):
        context = {}
        if self.simplestream.ui is not None:
            if "user_inputs" in self.simplestream.ui:
                if self.args is not None:
                    for input in self.simplestream.ui["user_inputs"]:
                        if input["key"] in self.args:
                            input["value"] = self.args[input["key"]]
                context["user_inputs"] = self.simplestream.ui["user_inputs"]
        response = render(self.request, self.template, context)
        return response