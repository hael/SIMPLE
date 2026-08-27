"""WSGI entry point for source-tree NICE development."""

import os

from django.core.wsgi import get_wsgi_application


os.environ.setdefault("DJANGO_SETTINGS_MODULE", "nice_project.settings_local")

application = get_wsgi_application()
