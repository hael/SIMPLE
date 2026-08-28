"""Validated NICE Lite feature flags."""

from django.conf import settings
from django.core.exceptions import ImproperlyConfigured


def batch_job_controls_enabled():
    """Return the explicitly configured batch-control feature state."""
    enabled = getattr(settings, "NICE_LITE_BATCH_JOB_CONTROLS", False)
    if not isinstance(enabled, bool):
        raise ImproperlyConfigured("NICE_LITE_BATCH_JOB_CONTROLS must be a boolean")
    return enabled


def workspace_job_refresh_enabled():
    """Return the explicitly configured workspace-card refresh state."""
    enabled = getattr(settings, "NICE_LITE_WORKSPACE_JOB_REFRESH", False)
    if not isinstance(enabled, bool):
        raise ImproperlyConfigured("NICE_LITE_WORKSPACE_JOB_REFRESH must be a boolean")
    return enabled
