# global imports
import json
from django.http                    import JsonResponse, HttpResponse
from django.contrib.auth.decorators import login_required
from django.views.decorators.cache  import cache_control
from django.views.decorators.csrf   import csrf_exempt

# local imports
from .helpers                    import *
from .data_structures.workspace    import Workspace
from .data_structures.project    import Project
from .data_structures.batchjob import BatchJob
from .data_structures.streamjob  import StreamJob

@csrf_exempt
def index(request):
    response = {}
    try:
        request_json = json.loads(request.body.decode('utf-8'))
        print("REQUEST JSON", request_json)
    except:
        print_error("failed to load json from request")
        return JsonResponse(response, status=400)
    if not float_present(request_json, "version"):
        print_error("version missing from request")
        return JsonResponse(response, status=400)  
    if not int_present(request_json, "jobid"):
        print_error("jobid missing from request")
        return JsonResponse(response, status=400)
    version = request_json["version"]
    jobid   = request_json["jobid"]
    if dict_present(request_json, "stream_heartbeat"):
        streamjob = StreamJob(id=jobid)
        if streamjob.update_stats(request_json):
           response = streamjob.get_master_update()
           print("API RESPONSE", response)
        else:
            return JsonResponse(response, status=400)
    elif dict_present(request_json, "heartbeat"):
        print("deal with classic job")
    else:
        print_error("unknown job heartbeat type")
        return JsonResponse(response, status=400)
    print("RESPONSE", response)
    return JsonResponse(response)

def index_classic(request):
    response = {}
    try:
        request_json = json.loads(request.body.decode('utf-8'))
    except:
        return JsonResponse(response)
    if "jobid" in request_json:
        job = BatchJob(id=request_json["jobid"])
        if job.jobmodel is None:
            return JsonResponse(response)
        workspace = Workspace(job.jobmodel.wspc_id)
        project = Project(job.jobmodel.wspc.proj_id)
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

