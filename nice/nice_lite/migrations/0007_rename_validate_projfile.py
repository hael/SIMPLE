# Release 4: the program validate_projfile became fix_projfile, which converts a project
# from an earlier release to the current stack-index contract. Rewrites stored program
# names; 0006 may already have been applied, so this is a separate migration.

from django.db import migrations

# old program name -> current program name
RENAMED_PROGS = {
    'validate_projfile': 'fix_projfile',
}


def rename_progs(apps, schema_editor):
    JobModel = apps.get_model('nice_lite', 'JobModel')
    for old, new in RENAMED_PROGS.items():
        JobModel.objects.filter(prog=old).update(prog=new)


def restore_progs(apps, schema_editor):
    JobModel = apps.get_model('nice_lite', 'JobModel')
    for old, new in RENAMED_PROGS.items():
        JobModel.objects.filter(prog=new).update(prog=old)


class Migration(migrations.Migration):

    dependencies = [
        ('nice_lite', '0006_rename_programs'),
    ]

    operations = [
        migrations.RunPython(rename_progs, restore_progs),
    ]
