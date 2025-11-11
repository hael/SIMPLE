# global imports
import json
from django.http import JsonResponse, HttpResponse
from django.contrib.auth.decorators import login_required
from django.views.decorators.cache import cache_control

# local imports
from .data_structures.workspace  import Workspace
from .data_structures.project    import Project
from .data_structures.job        import Job
from .data_structures.jobclassic import JobClassic

def index(request):
    response = {}
    try:
        request_json = json.loads(request.body.decode('utf-8'))
        print("REQUEST JSON", request_json)
    except:
        return JsonResponse(response)
    if "jobid" in request_json:
        job = Job(id=request_json["jobid"])
        response = job.updateStats(request_json)
    return JsonResponse(response)

def index_classic(request):
    response = {}
    try:
        request_json = json.loads(request.body.decode('utf-8'))
    except:
        return JsonResponse(response)
    if "jobid" in request_json:
        job = JobClassic(id=request_json["jobid"])
        workspace = Workspace(workspace_id=job.wspc.id)
        project   = Project(project_id=job.wspc.proj.id)
        response = job.updateStats(request_json, project, workspace)
    return JsonResponse(response)

@login_required(login_url="/login/")
@cache_control(max_age=300, must_revalidate=True, no_transform=True)
def image(request, src):
    response = HttpResponse()
    try:
        with open(src, "rb") as f:
            response = HttpResponse(f.read(), content_type="image/jpeg")
            return response
    except IOError:
        response = HttpResponse(content_type="image/jpeg")
        print("Image IO error :", src)
        return response
    else:
        print("Image error :", src)
        return response

