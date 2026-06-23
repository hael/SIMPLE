"""WorkspaceModel definition."""

from django.db     import models
from .projectmodel import ProjectModel

class WorkspaceModel(models.Model):
    """Workspace metadata stored under a project."""

    disp = models.IntegerField(default=0)
    jcnt = models.IntegerField(default=0)
    proj = models.ForeignKey(ProjectModel, on_delete=models.CASCADE)
    name = models.CharField(max_length=200, default='')
    desc = models.CharField(max_length=200, default='')
    dirc = models.CharField(max_length=200, default='')
    link = models.CharField(max_length=200, default='')
    user = models.CharField(max_length=50, default='')
    cdat = models.DateTimeField(auto_now_add=True, blank=True)
    mdat = models.DateTimeField(auto_now_add=True, blank=True)
    # Backward-compatibility payload used by legacy stream state.
    nstr = models.JSONField(default=dict, null=True, blank=True)
