import os
from django.core.management.base import BaseCommand, CommandError
from django.contrib.auth.models import User

from ...models import DispatchModel

class Command(BaseCommand):
    help = "Ensures local user is present - do not use on non-local installs"

    def handle(self, *args, **options):
        if not User.objects.filter(username=os.getlogin()).exists():
            User.objects.create_superuser(os.getlogin())
        localuser = User.objects.get(username=os.getlogin())
        localuser.set_password("simple")
        localuser.save()
        
       