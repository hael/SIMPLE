"""
WSGI config for nice_project project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/5.2/howto/deployment/wsgi/
"""

import os

from django.core.wsgi import get_wsgi_application

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'nice_project.settings')
os.environ["PATH"] = '/home/jcaesar/Code/SIMPLE/build/bin' + os.pathsep + os.environ["PATH"]
# add any desired environment variables here
#os.environ["LD_LIBRARY_PATH"] = '/path/to/lib' + os.pathsep + os.environ["LD_LIBRARY_PATH"]
#


application = get_wsgi_application()
