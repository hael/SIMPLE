# global imports
from django.contrib import admin

# local imports
from .models import ProjectModel, WorkspaceModel, JobModel, DispatchModel

admin.site.register(ProjectModel)
admin.site.register(WorkspaceModel)
admin.site.register(JobModel)
admin.site.register(DispatchModel)