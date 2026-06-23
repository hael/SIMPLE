"""Workspace landing view for stream mode."""

from django.contrib                 import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts               import render

from ..data_structures.simple    import SIMPLEStream
from ..data_structures.streamjob import StreamJob
from ..helpers                   import clear_checksum_cookies, get_job_id

@login_required(login_url="/login/")
def view_job_builder_stream(request):
    """Render stream job-builder page for a new job or from an existing stream job."""
    template     = "nice_stream/newstream.html"
    jobid        = get_job_id(request)
    ui           = None
    args         = None
    if jobid is not None:
        streamjob      = StreamJob(jobid)
        streamjobmodel = streamjob.get_jobmodel()
        args           = streamjobmodel.args
    simplestream = SIMPLEStream()
    if simplestream.loadUIJSON():
        ui = simplestream.get_ui()
    else:
        messages.add_message(request, messages.ERROR, "failed to read ui JSON")
    context = {}
    if ui is not None:
        if "user_inputs" in ui:
            if args is not None:
                for input in ui["user_inputs"]:
                    if input["key"] in args:
                        input["value"] = args[input["key"]]
        context["user_inputs"] = ui["user_inputs"]
    response = render(request, template, context)
    clear_checksum_cookies(request, response)
    return response