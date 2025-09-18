# global imports
from django.contrib import admin

# local imports
from .models import ProjectModel, DatasetModel, JobModel, DispatchModel, WorkspaceModel, JobClassicModel

admin.site.register(ProjectModel)
admin.site.register(DatasetModel)
admin.site.register(JobModel)
admin.site.register(JobClassicModel)
admin.site.register(DispatchModel)
admin.site.register(WorkspaceModel)