"""DispatchModel definition."""

from django.db import models

class DispatchModel(models.Model):
    """Stores command-dispatch templates and activation state."""

    scmd = models.CharField(max_length=200, default='')
    url  = models.CharField(max_length=200, default='')
    tplt = models.TextField(default='')
    active = models.BooleanField(default=False)
