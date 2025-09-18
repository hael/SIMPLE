import os
from django.core.management.base import BaseCommand, CommandError

from ...models import DispatchModel

class Command(BaseCommand):
    help = "Updates the job submission template"

    def add_arguments(self, parser):
        parser.add_argument("command", nargs="?", type=str)
        parser.add_argument("url",     nargs="?", type=str)
        parser.add_argument("script",  nargs="?", type=str)

    def handle(self, *args, **options):
        if options["command"] is None:
            raise CommandError("command missing from arguments")
        if options["url"] is None:
            raise CommandError("url missing from arguments")
        if options["script"] is None:
            raise CommandError("script missing from arguments")
        if not os.path.exists(options["script"]):
            raise CommandError(options["script"], "does not exist")
        if not os.path.isfile(options["script"]):
            raise CommandError(options["script"], "is not a file")
        try:
            f = open(options["script"], 'r')
        except OSError:
            raise CommandError("Could not open/read file:", fname)
        except:
            raise CommandError("Failed to open file:", fname)
        template = f.read()
        f.close()
        dispatchmodels = DispatchModel.objects.all()
        for dispatchmodel in dispatchmodels:
            dispatchmodel.active = False
            dispatchmodel.save()     
        newdispatchmodel = DispatchModel(scmd=options["command"], tplt=template, url=options["url"], active=True)
        newdispatchmodel.save()