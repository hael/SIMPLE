"""ProjectModel definition."""

from django.db import models

class ProjectModel(models.Model):
    """Top-level project metadata and root directory mapping."""

    name = models.CharField(max_length=200)
    date = models.DateTimeField()
    desc = models.CharField(max_length=1024, default='')
    dirc = models.CharField(max_length=1024, default='')
