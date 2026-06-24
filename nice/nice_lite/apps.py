"""Django application configuration for NICE Lite.

This module defines app metadata used by Django app discovery and default
model primary-key field configuration.
"""

# django imports
from django.apps import AppConfig


class NiceLiteConfig(AppConfig):
    """Application config for the ``nice_lite`` Django app."""

    default_auto_field = "django.db.models.BigAutoField"
    name = "nice_lite"
