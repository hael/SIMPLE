# Django imports
from django.shortcuts import render
from django.contrib   import messages

# Local imports
from ..helpers                 import *
from ..data_structures.project import Project
from ..data_structures.simple  import SIMPLEStream

class NewStreamView:

    template = "nice_stream/newstream.html"
    ui       = None
    args     = None

    def __init__(self, request, args=None):
        self.request = request
        simplestream = SIMPLEStream()
        if simplestream.loadUIJSON():
            self.ui = simplestream.get_ui()
        else:
            messages.add_message(self.request, messages.ERROR, "failed to read ui JSON")
        if args is not None:
            self.args = args

    def render(self):
        context = {}
        if self.ui is not None:
            if "user_inputs" in self.ui:
                if self.args is not None:
                    for input in self.ui ["user_inputs"]:
                        if input["key"] in self.args:
                            input["value"] = self.args[input["key"]]
            context["user_inputs"] = self.ui["user_inputs"]
        print(context)
        response = render(self.request, self.template, context)
        clear_checksum_cookies(self.request, response)
        return response