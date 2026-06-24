"""Django admin registrations for NICE Lite models.

This module exposes core project, workspace, job, and dispatch records in the
admin site for operational inspection and maintenance tasks.
"""

# django imports
from django.contrib import admin

# local imports
from .models import DispatchModel
from .models import JobModel
from .models import ProjectModel
from .models import WorkspaceModel


# ------------------------------------------------------------------
# Admin Registrations
# ------------------------------------------------------------------

admin.site.register(ProjectModel)
admin.site.register(WorkspaceModel)
admin.site.register(JobModel)
admin.site.register(DispatchModel)