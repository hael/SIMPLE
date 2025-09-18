# global imports
import json
from django.http import JsonResponse, HttpResponse
from django.contrib.auth.decorators import login_required
from django.views.decorators.cache import cache_control

# local imports
from .data_structures.project  import Project
from .data_structures.dataset  import Dataset
from .data_structures.job      import Job

def index(request):
    response = {}
    try:
        request_json = json.loads(request.body.decode('utf-8'))
    except:
        return JsonResponse(response)
    if "jobid" in request_json:
        job = Job(id=request_json["jobid"])
        response = job.updateStats(request_json)
    return JsonResponse(response)

@login_required(login_url="/login/")
@cache_control(max_age=3600, must_revalidate=True, no_transform=True)
def image(request, src):
    try:
        with open(src, "rb") as f:
            response = HttpResponse(f.read(), content_type="image/jpeg")

            return response
    except IOError:
        red = Image.new('RGBA', (1, 1), (255,0,0,0))
        response = HttpResponse(content_type="image/jpeg")
        red.save(response, "JPEG")
        return response
    else:
        red = Image.new('RGBA', (1, 1), (255,0,0,0))
        response = HttpResponse(content_type="image/jpeg")
        red.save(response, "JPEG")
        return response 
