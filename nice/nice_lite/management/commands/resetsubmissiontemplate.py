import os
from django.core.management.base import BaseCommand, CommandError

from ...models import DispatchModel

class Command(BaseCommand):
    help = "Resets the submission template to run locally"

    def handle(self, *args, **options):
        dispatchmodels = DispatchModel.objects.all()
        for dispatchmodel in dispatchmodels:
            dispatchmodel.active = False
            dispatchmodel.save()
        template="#"   # needs to be populated
        newdispatchmodel = DispatchModel(scmd="bash", tplt=template, url="localhost:8000", active=True)
        newdispatchmodel.save()