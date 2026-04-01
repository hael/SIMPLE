import os
from django.core.management.base import BaseCommand, CommandError

from ...models import DispatchModel

class Command(BaseCommand):
    help = "Resets the submission template to run locally"

    def add_arguments(self, parser):
        parser.add_argument("simple_path", nargs="?", type=str)

    def handle(self, *args, **options):
        if options["simple_path"] is None:
            raise CommandError("simple_path missing from arguments")
        localtemplate = f"""#!/usr/bin/env bash
# localtemplate
export SIMPLE_PATH={options["simple_path"]}
export SIMPLE_EMAIL=''
export SIMPLE_QSYS=local

echo $$ > nice.pid
echo "Running on $HOSTNAME"
XXXSIMPLEXXX"""
        dispatchmodels = DispatchModel.objects.all()
        for dispatchmodel in dispatchmodels:
            dispatchmodel.delete()
        newdispatchmodel = DispatchModel(scmd="nohup", tplt=localtemplate, url="localhost:8000", active=True)
        newdispatchmodel.save()
