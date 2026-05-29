import os
from django.core.management.base import BaseCommand, CommandError
from django.contrib.auth.models import User

from ...models import DispatchModel

class Command(BaseCommand):
    help = "Ensures local user is present - do not use on non-local installs"

    def handle(self, *args, **options):
        user_name = os.environ["USER"] or os.environ["LOGNAME"] or os.environ["USERNAME"] or os.getlogin()
        if not User.objects.filter(username=user_name).exists():
            User.objects.create_superuser(user_name)
        localuser = User.objects.get(username=user_name)
        localuser.set_password("simple")
        localuser.save()
        
       